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


# Global mesh sizes (full grid, all partitions).
numNodes = 0
numElements = 0
numSides = 0
numSideBnds = 0
nvrt = 0

# Per-partition global ids: nodeGlobalIndices[m][i] = global node id of local node i on GPU m.
# Set by initGlobalIndices (required for both single- and multi-GPU).
nodeGlobalIndices = None
elemGlobalIndices = None
sideGlobalIndices = None


def init(nNodes, nElements, nSides, nSideBnds, numVrt):
    """Store full-mesh dimensions and vertical level count."""
    global numNodes, numElements, numSides, numSideBnds, nvrt
    numNodes = nNodes
    numElements = nElements
    numSides = nSides
    numSideBnds = nSideBnds
    nvrt = numVrt


def initGlobalIndices(nodeIndices, elemIndices, sideIndices):
    """
    Per-partition global SCHISM indices for scatter after GPU read.

    nodeIndices[m] is a 1D int array of length np_m (global node id of each local node on GPU m).
    Same pattern for elemIndices and sideIndices.

    Called once from main after partPopulate. With one GPU, nodeIndices is a one-element list
    (e.g. [np.arange(numNodes)]). With several GPUs, one entry per device from globalNodeIndices.
    """
    global nodeGlobalIndices, elemGlobalIndices, sideGlobalIndices
    nodeGlobalIndices = [np.asarray(a, dtype=np.int32) for a in nodeIndices]
    elemGlobalIndices = [np.asarray(a, dtype=np.int32) for a in elemIndices]
    sideGlobalIndices = [np.asarray(a, dtype=np.int32) for a in sideIndices]


def getViewLayout(varType, is3D, numItems, structSize, kind):
    """
    Shape, strides, and dtype for a CuPy ndarray view into a struct array on one device.

    kind is "node", "side", or "element".
    structSize is sizeof(Node), sizeof(Side), or sizeof(Element) in bytes.
    numItems is the local count on that device (np, ns, or ne), not the global count.

    Returns (shape, strides, dtype).
    """
    raise NotImplementedError


def getNodeVar(varName, k=-1):
    """
    Read a nodal field from all GPU partitions into one host array.

    varName must appear in varInfoDict (e.g. "eta", "depth", "outputCoords").
    k selects a vertical level for 3D fields; must be -1 for 2D fields.

    Returns a NumPy array in global node order (length numNodes along the last axis).
    """

    varInfo = varInfoDict[varName]
    varCode = list(varInfoDict.keys()).index(varName)
    varType = varInfo['type']
    is3D = varInfo['is3D']


    if not is3D and k != -1:
        raise ValueError("k cannot be provided for 2D variables")

    # Gets information about the location and size of the variable in the GPU memory
    # [This should be lists for each device]
    offset      = cp.zeros(1, dtype=cp.uint64)   
    addr        = cp.zeros(1, dtype=cp.uint64)
    sizeNode    = cp.zeros(1, dtype=cp.uint64)
    sizeType    = cp.zeros(1, dtype=cp.uint64)
    cufunc.getNodeVarAddr(varCode, addr, offset, sizeNode, sizeType)
    sizeNode = sizeNode.item()

    if offset.get() != cp.uint32(cp.int32(-1)):
        # If offset is not -1, then the variable is directly mapped in memory

        if varType == 'double':
            if is3D:
                shape = (nvrt, numNodes)
                strides = (8, sizeNode)  # 8 bytes for double, then stride for next element
            else:
                shape = (numNodes,)
                strides = (sizeNode,)
            dtype = cp.float64
        elif varType == 'Vector2':
            if is3D:
                shape = (2, nvrt, numNodes)  # XXXXX  check 2*nvrt
                strides = (8, 16, sizeNode)  # 8 bytes for double, then 8*nvrt for next level, then stride for next element
            else:
                shape = (2, numNodes)
                strides = (8, sizeNode)
            dtype = cp.float64
        elif varType == 'int' or varType == 'bool':
            if is3D:
                shape = (nvrt, numNodes)
                strides = (4, sizeNode)  # 4 bytes for int, then stride for next element
            else:
                shape = (numNodes,)
                strides = (sizeNode,)
            dtype = cp.int32
        else:
            raise ValueError(f"Invalid variable type: {varType}")
            
        view = cp.ndarray(shape=shape, dtype=dtype, strides=strides, order='F',
                          memptr=cp.cuda.MemoryPointer(cp.cuda.UnownedMemory(addr.item() + offset.item(), 0, owner=None), 0)) 

        return view

    elif varName == "outputCoords":
            res = cp.zeros((2, numNodes), dtype=cp.float64, order='F')
            cufunc.getNodeVarFunc(varCode, res)
            return res
    else:
        raise ValueError(f"Variable {varName} is not mapped in memory")    


def getSideVar(varName, k=-1):
    """
    Read a side field from all GPU partitions into one host array.

    Returns a NumPy array in global side order (length numSides along the last axis).
    """
    raise NotImplementedError


def getElementVar(varName, k=-1):
    """
    Read an element field from all GPU partitions into one host array.

    Returns a NumPy array in global element order (length numElements along the last axis).
    """
    raise NotImplementedError
