"""
GPU → CPU (DtoH) mesh field access for SCHISM (gschism).

Each CUDA device holds a partition of the mesh (local Node / Side / Element arrays).
This module reads named fields from device memory and assembles them into full-mesh
NumPy arrays in host RAM, indexed in global SCHISM order (0 .. nNodes-1, etc.).

Used by writeOutput and other I/O; not for plotting or debugging (see debug.monitor).

Variable names and types come from schism_consts.varInfoDict.
Mapped fields use get*VarAddr (SoA view into local*Data); some fields (e.g. outputCoords)
use get*VarFunc kernels that write by global idx.

Single-GPU and multi-GPU use the same logic: always call init and initGlobalIndices,
then get*Var loops over every device in configCuda.devices, reads the local SoA slice,
and scatters into the full host array via nodeGlobalIndices / elemGlobalIndices /
sideGlobalIndices. With one device, those lists have length 1 (typically arange(n)).

main.py fills the index lists from partPopulate the same way regardless of partition count.
"""

from __future__ import annotations

import numpy as np
import cupy as cp
import python.cuda_functions as cufunc
from python.schism_consts import varInfoDict




# Per-partition global ids: nodeGlobalIndices[m][i] = global node id of local node i on GPU m.
# Set by initGlobalIndices (required for both single- and multi-GPU).



class DeviceToHost:
    def __init__(self, nNodes, nElements, nSides, nSideBnds, numVrt, 
                 globalNodeIndices, globalElemIndices, globalSideIndices, configCuda):
        """Store full-mesh dimensions and vertical level count."""
        self.numNodes = nNodes
        self.numElements = nElements
        self.numSides = nSides
        self.numSideBnds = nSideBnds
        self.nvrt = numVrt
        self.configCuda = configCuda

        # Allocate buffers for each device. These will store a variable in packed index order.
        self.bufferGPU = []   
        self.bufferCPU = []
        for m, device_id in enumerate(self.configCuda.devices):
            nMax = max(len(globalNodeIndices[m]),
                       len(globalElemIndices[m]),
                       len(globalSideIndices[m]))
            maxFieldBytes = max(32, numVrt * 16)  # 3D Vector2 is the largest packed field
            with cp.cuda.Device(device_id):                
                self.bufferGPU.append(cp.empty(nMax*maxFieldBytes, dtype=cp.uint8))
                self.bufferCPU.append(np.empty(nMax*maxFieldBytes, dtype=np.uint8))

        # Allocate a single result buffer for the final assembly.
        totScratchBytes = sum(len(buffer) for buffer in self.bufferCPU)
        self.result = np.empty(totScratchBytes, dtype=np.uint8)

        # Per-partition global SCHISM indices for scatter after GPU read.
        self.globalNodeIndices = [cp.asnumpy(a).astype(np.int32, copy=False) for a in globalNodeIndices]
        self.globalElemIndices = [cp.asnumpy(a).astype(np.int32, copy=False) for a in globalElemIndices]
        self.globalSideIndices = [cp.asnumpy(a).astype(np.int32, copy=False) for a in globalSideIndices]




    def getViewLayout(self, varType, is3D, numItems, structSize, kind):
        """
        Shape, strides, and dtype for a CuPy ndarray view into a struct array on one device.

        kind is "node", "side", or "element".
        structSize is sizeof(Node), sizeof(Side), or sizeof(Element) in bytes.
        numItems is the local count on that device (np, ns, or ne), not the global count.

        Returns (shape, strides, dtype).
        """
        raise NotImplementedError


    def varLayout(self, varName, k=-1):
        """
        Packed field layout from varInfoDict.

        Returns (varCode, varType, is3D, nbytes, dtype, numComponents).
        nbytes is bytes written per local entity and must match the kernel.
        """
        varInfo = varInfoDict[varName]
        varCode = list(varInfoDict.keys()).index(varName)
        varType = varInfo['type']
        is3D = varInfo['is3D']

        if not is3D and k != -1:
            raise ValueError("k cannot be provided for 2D variables")

        levels = self.nvrt if (is3D and k == -1) else 1

        if varType == 'double':
            nbytes = 8 * levels
            dtype = np.float64
            numComponents = levels
        elif varType == 'Vector2':
            nbytes = 16 * levels
            dtype = np.float64
            numComponents = 2 * levels
        elif varType in ('int', 'bool'):
            nbytes = 4
            dtype = np.int32
            numComponents = 1
        elif varType == 'FaceIndices':
            nbytes = 32
            dtype = np.uint64
            numComponents = 4
        elif varType == 'pointer':
            nbytes = 8
            dtype = np.uint64
            numComponents = 1
        else:
            raise ValueError(varType)

        return varCode, varType, is3D, nbytes, dtype, numComponents


    def getVar(self, varName, k, globalIndices, numItems, kernel):
        """
        Read a field from all GPU partitions into one host array in global SCHISM order.

        globalIndices[m] maps local packed index on device m to the global id.
        numItems is the full-mesh count (numNodes / numSides / numElements).
        kernel writes the packed local field into bufferGPU[m].
        """
        varCode, varType, is3D, nbytes, dtype, numComponents = self.varLayout(varName, k)
        raw = self.result[:numItems * nbytes].view(dtype)
        if is3D and varType == 'Vector2':
            result = raw.reshape(numItems, self.nvrt, 2)
        elif is3D and varType == 'double':
            result = raw.reshape(numItems, self.nvrt)
        else:
            result = raw.reshape(numItems, numComponents)

        for m, device_id in enumerate(self.configCuda.devices):
            gidx = globalIndices[m]
            sizePart = len(gidx)

            with cp.cuda.Device(device_id):
                kernel(varCode, self.bufferGPU[m], k)
                cp.asnumpy(self.bufferGPU[m][:sizePart*nbytes], out=self.bufferCPU[m][:sizePart*nbytes])

            partRaw = self.bufferCPU[m][:sizePart * nbytes].view(dtype)
            if is3D and varType == 'Vector2':
                arrayPart = partRaw.reshape(sizePart, self.nvrt, 2)
            elif is3D and varType == 'double':
                arrayPart = partRaw.reshape(sizePart, self.nvrt)
            else:
                arrayPart = partRaw.reshape(sizePart, numComponents)
            result[gidx] = arrayPart

        if varType == 'Vector2' and not is3D:
            return result.T
        if is3D and varType == 'Vector2':
            return result.transpose(2, 1, 0)
        if is3D and varType == 'double':
            return result.T
        return result.squeeze()


    def getNodeVar(self, varName, k=-1):
        """
        Read a nodal field from all GPU partitions into one host array.

        Returns a NumPy array in global node order (length numNodes along the last axis).
        """
        return self.getVar(varName, k, self.globalNodeIndices, self.numNodes, cufunc.getNodeVar)


    def getSideVar(self, varName, k=-1):
        """
        Read a side field from all GPU partitions into one host array.

        Returns a NumPy array in global side order (length numSides along the last axis).
        """
        return self.getVar(varName, k, self.globalSideIndices, self.numSides, cufunc.getSideVar)


    def getElementVar(self, varName, k=-1):
        """
        Read an element field from all GPU partitions into one host array.

        Returns a NumPy array in global element order (length numElements along the last axis).
        """
        return self.getVar(varName, k, self.globalElemIndices, self.numElements, cufunc.getElementVar)
