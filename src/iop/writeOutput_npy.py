import os
import cupy as cp
import numpy as np

import netCDF4 as nc
import debug.monitor as monitor 

import threading
import queue

# Add these global variables
save_thread = None
save_queue = queue.Queue()
save_lock = threading.Lock()


path = None
nNodes = None
nElements = None
nSides = None
nLayers = None
sizeNode = None
sizeElem = None
sizeSide = None



def initializeOutput(basePath, numNodes, numElements, numSides, numLayers, nodeSize, elemSize, sideSize):
    """Initialize some module-level global variables to simplify calls to the writeOutput functions."""
    
    global path, nNodes, nElements, nSides, nLayers, sizeNode, sizeElem, sizeSide
    path = basePath
    nNodes = numNodes
    nElements = numElements
    nSides = numSides
    nLayers = numLayers
    sizeNode = nodeSize
    sizeElem = elemSize
    sizeSide = sideSize


def _saveDataThread(data, numFile):
    """Background thread for saving data to file."""
    try:
        filepath = os.path.join(path, f'res_{numFile}.npy')
        np.save(filepath, data)
        print(f"[Save] Data saved to {filepath}")
        
        # Clean up memory
        del data
        
    except Exception as e:
        print(f"[Save] Error saving data: {e}")
    finally:
        print("[Save] Save thread finished")



def writeOutput(it, core, opt):
    # Write the output files
    global sizeNode, save_thread

    ihfskip = core.ihfskip
    nspool  = core.nspool
    dt      = core.dt
    timesPerFile = ihfskip // nspool


    if (it % nspool) != 0:
        return

    print(f"it = {it}, nspool = {nspool}, ihfskip = {ihfskip}")
    isNewFile = (it - nspool) % ihfskip == 0
    numFile = it//nspool
    idxTime = ((it // nspool) - 1) % timesPerFile

    data = {'numFile': numFile, 'idxTime': idxTime, 'isNewFile': isNewFile}

    # Horizontal geometric variables are only written once per file
    if idxTime == 0:
        view = monitor.getNodeVar("x")  # x is a 2D vector
        xy = view.get()
        data['SCHISM_hgrid_node_x'] = xy[0, :].copy()
        data['SCHISM_hgrid_node_y'] = xy[1, :].copy()
        del xy

        view = monitor.getSideVar("x")  # x is a 2D vector
        xy = view.get()
        data['SCHISM_hgrid_edge_x'] = xy[0, :].copy()
        data['SCHISM_hgrid_edge_y'] = xy[1, :].copy()
        del xy

        view = monitor.getElementVar("x")  # x is a 2D vector
        xy = view.get()
        data['SCHISM_hgrid_face_x'] = xy[0, :].copy()
        data['SCHISM_hgrid_face_y'] = xy[1, :].copy()
        del xy

        view = monitor.getElementVar("FaceIndices")
        mask = (view[:,:] != 0)
        faceIndices = view.astype(cp.int64).get()
        firstPointer = view[mask].min().get()
        mask = mask.get()
        # To convert pointers to indices, we subtract the first pointer and then divide by the size of the node and add 1 (for 1-indexing)
        faceIndices[mask] = (faceIndices[mask] - firstPointer)//sizeNode + 1
        faceIndices[~mask] = -1
        data['SCHISM_hgrid_face_nodes'] = faceIndices.T.copy()
        del faceIndices, mask

        view = monitor.getSideVar("n1")
        n1 = view.get()
        firstPointer = view.min().get()
        # To convert pointers to indices, we subtract the first pointer and then divide by the size of the node and add 1 (for 1-indexing)
        n1 = (n1 - firstPointer)//sizeNode + 1

        view = monitor.getSideVar("n2")
        n2_cpu = view.get()
        firstPointer = view.min().get()
        # To convert pointers to indices, we subtract the first pointer and then divide by the size of the node and add 1 (for 1-indexing)
        n2 = (n2_cpu - firstPointer)//sizeNode + 1        

        data['SCHISM_hgrid_edge_nodes'] = np.zeros((nSides, 2), dtype=np.int32)
        data['SCHISM_hgrid_edge_nodes'][:, 0] = n1
        data['SCHISM_hgrid_edge_nodes'][:, 1] = n2
        del n1, n2


        view = monitor.getNodeVar("depth")
        data['depth'] = view.get()
 

        data['minimum_depth'] = 0.0001  #XXXXXXXXXXX Should be h0


    # For each variable, get the data from GPU (using monitor) and store it in the NetCDF file.
    # Explicitly delete the data from CPU after storing it in the NetCDF file.
    data['elevation']      = monitor.getNodeVar("eta").get()

    data['dryFlagNode']    = monitor.getNodeVar("isdry").astype(cp.int32).get()

    data['dryFlagSide']    = monitor.getSideVar("isdry").astype(cp.int32).get()

    data['dryFlagElement'] = monitor.getElementVar("isdry").astype(cp.int32).get()


    view = monitor.getNodeVar("u")
    u = view[0,:,:].get()
    v = view[1,:,:].get()
    data['horizontalVelX'] = u.T.copy()
    data['horizontalVelY'] = v.T.copy()
    del u, v


    data['zCoordinates'] = monitor.getNodeVar("z").get().T

    data['time'] = it*dt

    # Save velocity components and time to res.npy
    # np.save(os.path.join(path, f'res_{numFile}.npy'), data)

    with save_lock:
        if save_thread is not None and save_thread.is_alive():
            print("[Save] Waiting for previous save thread to finish")
            save_thread.join()  # Wait until it finishes
            print("[Save] Previous save thread finished")

    # Launch new thread for saving
    with save_lock:
        save_thread = threading.Thread(target=_saveDataThread, args=(data, numFile))
        save_thread.start()







    




