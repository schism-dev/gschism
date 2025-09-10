from python.schism_consts import *
import python.cuda_functions as cufunc
import cupy as cp
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
import matplotlib
from python.schism_consts import varInfoDict

matplotlib.use('TkAgg') 



debugPlot = cp.zeros([10000,5], dtype=cp.float64, order='F') 
debugPlotIdx = cp.zeros(1, dtype=cp.int32, order='F') 
debugPlotIdx[0] = 0


fig = None
ax  = None


numNodes = 0
numElements = 0
numSides = 0
numSideBnds = 0
nvrt = 0

def init(nNodes, nElements, nSides, nSideBnds, numVrt):
    global numNodes, numElements, numSides, numSideBnds, nvrt
    numNodes = nNodes
    numElements = nElements
    numSides = nSides
    numSideBnds = nSideBnds
    nvrt = numVrt



def getElementVar(varName, k = -1):
    varInfo = varInfoDict[varName]
    varCode = list(varInfoDict.keys()).index(varName)
    varType = varInfo['type']
    is3D = varInfo['is3D']

    # if is3D and k == -1:
    #     raise ValueError("k must be provided for 3D variables")

    if not is3D and k != -1:
        raise ValueError("k cannot be provided for 2D variables")


    offset      = cp.zeros(1, dtype=cp.uint64)   
    addr        = cp.zeros(1, dtype=cp.uint64)
    sizeElement = cp.zeros(1, dtype=cp.uint64)
    sizeType    = cp.zeros(1, dtype=cp.uint64)
    cufunc.getElementVarAddr(varCode, addr, offset, sizeElement, sizeType)
    sizeElement = sizeElement.item()

    if offset.get() != -1:
        # If offset is not -1, then the variable is directly mapped in memory

        if varType == 'double':
            if is3D:
                shape = (nvrt, numElements)
                strides = (8, sizeElement)  # 8 bytes for double, then stride for next element
            else:
                shape = (numElements,)
                strides = (sizeElement,)
            dtype = cp.float64
        elif varType == 'Vector2':
            if is3D:
                shape = (2, nvrt, numElements)
                strides = (8, 8*nvrt, sizeElement)  # 8 bytes for double, then 8*nvrt for next level, then stride for next element
            else:
                shape = (2, numElements)
                strides = (8, sizeElement)
            dtype = cp.float64
        elif varType == 'int' or varType == 'bool':
            if is3D:
                shape = (nvrt, numElements)
                strides = (4, sizeElement)  # 4 bytes for int, then stride for next element
            else:
                shape = (numElements,)
                strides = (sizeElement,)
            dtype = cp.int32
        elif varType == 'FaceIndices':
            shape = (4, numElements)
            strides = (8, sizeElement,) # 8 bytes for the pointer
            dtype = cp.uint64
        elif varType == 'pointer':
            shape = (numElements,)
            strides = (sizeElement,)
            dtype = cp.uint64
        else:
            raise ValueError(f"Invalid variable type: {varType}")
            
        view = cp.ndarray(shape=shape, dtype=dtype, strides=strides, order='F',
                          memptr=cp.cuda.MemoryPointer(cp.cuda.UnownedMemory(addr.item() + offset.item(), 0, owner=None), 0)) 

        return view

    

    # var = cp.zeros(numElements, dtype=cp.float64)
    # cufunc.getElementVar(varCode, var, k, j)
    # return var.get()    
    raise ValueError(f"Variable {varName} is not mapped in memory")


def getSideVar(varName, k = -1):
    varInfo = varInfoDict[varName]
    varCode = list(varInfoDict.keys()).index(varName)
    varType = varInfo['type']
    is3D = varInfo['is3D']

    # if is3D and k == -1:
    #     raise ValueError("k must be provided for 3D variables")

    if not is3D and k != -1:
        raise ValueError("k cannot be provided for 2D variables")


    offset      = cp.zeros(1, dtype=cp.uint64)   
    addr        = cp.zeros(1, dtype=cp.uint64)
    sizeSide    = cp.zeros(1, dtype=cp.uint64)
    sizeType    = cp.zeros(1, dtype=cp.uint64)
    cufunc.getSideVarAddr(varCode, addr, offset, sizeSide, sizeType)
    sizeSide = sizeSide.item()

    if offset.get() != -1:
        # If offset is not -1, then the variable is directly mapped in memory

        if varType == 'double':
            if is3D:
                shape = (nvrt, numSides)
                strides = (8, sizeSide)  # 8 bytes for double, then stride for next element
            else:
                shape = (numSides,)
                strides = (sizeSide,)
            dtype = cp.float64
        elif varType == 'Vector2':
            if is3D:
                shape = (2, nvrt, numSides)
                strides = (8, 8*nvrt, sizeSide)  # 8 bytes for double, then 8*nvrt for next level, then stride for next element
            else:
                shape = (2, numSides)
                strides = (8, sizeSide)
            dtype = cp.float64
        elif varType == 'int' or varType == 'bool':
            if is3D:
                shape = (nvrt, numSides)
                strides = (4, sizeSide)  # 4 bytes for int, then stride for next element
            else:
                shape = (numSides,)
                strides = (sizeSide,)
            dtype = cp.int32
        elif varType == 'pointer':
            if is3D:
                shape = (nvrt, numSides)
                strides = (4, sizeSide)  # 4 bytes for int, then stride for next element
            else:
                shape = (numSides,)
                strides = (sizeSide,)
            dtype = cp.uint64            
        else:
            raise ValueError(f"Invalid variable type: {varType}")
            
        view = cp.ndarray(shape=shape, dtype=dtype, strides=strides, order='F',
                          memptr=cp.cuda.MemoryPointer(cp.cuda.UnownedMemory(addr.item() + offset.item(), 0, owner=None), 0)) 

        return view

    

    # var = cp.zeros(numElements, dtype=cp.float64)
    # cufunc.getElementVar(varCode, var, k, j)
    # return var.get()    
    raise ValueError(f"Variable {varName} is not mapped in memory")

def getNodeVar(varName, k = -1):
    varInfo = varInfoDict[varName]
    varCode = list(varInfoDict.keys()).index(varName)
    varType = varInfo['type']
    is3D = varInfo['is3D']

    # if is3D and k == -1:
    #     raise ValueError("k must be provided for 3D variables")

    if not is3D and k != -1:
        raise ValueError("k cannot be provided for 2D variables")


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
                strides = (8, 8*nvrt, sizeNode)  # 8 bytes for double, then 8*nvrt for next level, then stride for next element
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

  
    raise ValueError(f"Variable {varName} is not mapped in memory")


# XXXX REMOVE THESE FUNCTIONS and the corresponding ones in the kernels
def getElementVarX(varname, k = -1, j = -1):
    varInfo = varInfoDict[varname]
    var = cp.zeros(numElements, dtype=cp.float64)
    cufunc.getElementVar(varCode, var, k, j)
    return var.get()

def getNodeVarX(varname, k = -1):
    varInfo = varInfoDict[varname]
    var = cp.zeros(numNodes, dtype=cp.float64)
    cufunc.getNodeVar(varCode, var, k)
    return var.get()

def getSideVarX(varname, k = -1):
    varInfo = varInfoDict[varname]
    var = cp.zeros(numSides, dtype=cp.float64)
    cufunc.getSideVar(varCode, var, k)
    return var.get()

def getSideBndVar(varname, k = -1):
    varInfo = varInfoDict[varname]
    var = cp.zeros(numSideBnds, dtype=cp.float64)
    cufunc.getSideBndVar(varCode, var, k)
    return var.get()    


def intToRGB(color_array):
    """Convert 0xRRGGBB integer to matplotlib (r, g, b) tuple in [0, 1]."""
    """Convert an array of 0xRRGGBB ints to an array of (r,g,b) tuples in [0,1]."""
    color_array = np.asarray(color_array).astype(np.int32)
    r = (color_array >> 16) & 0xFF
    g = (color_array >> 8) & 0xFF
    b = color_array & 0xFF
    return np.stack([r, g, b], axis=1) / 255.0  # shape: (N, 3)

def plotMesh():
    global fig, ax

    res = True
    
    if fig is None: 
        fig, ax = plt.subplots(figsize=(12, 8))
    
    # Remove all existing plots
    ax.clear()

    x = getNodeVar("x")
    y = getNodeVar("y")
        
        # Plot all nodes
    plt.plot(x, y, 'ko', alpha=0.93, markersize=2.3)

    n1 = getSideVar("n1").astype(int)
    n2 = getSideVar("n2").astype(int)
    sx = getSideVar("x")
    sy = getSideVar("y")
    ux = getSideVar("u", 0)
    uy = getSideVar("v", 0)
    eta = getNodeVar("eta")
    x  = getNodeVar("x")
    y  = getNodeVar("y")

    # Build the segments array for LineCollection
    segments = np.stack([np.column_stack([x[n1], y[n1]]), np.column_stack([x[n2], y[n2]])], axis=1)

    # Plots eta at the element centers. Eta is shown as the radius of the circle.
    plt.scatter(x, y, color='b', s=(eta*200 + 20))

    # Create the LineCollection
    lc = LineCollection(segments, colors='r', alpha=0.53, linewidths=0.5)
    ax.add_collection(lc)

    # Use a fixed scale factor for consistent arrow length representation
    # Adjust this value to make arrows longer or shorter while maintaining 1 m/s = constant length
    fixed_scale = 10.001  # Adjust this value as needed for your visualization
    
    plt.quiver(sx, sy, fixed_scale*ux, fixed_scale*uy, angles='xy', scale_units='inches', scale=1, color='b', width=0.003)

    # Plot the backtracked curves.
    segments_np = cp.asnumpy(debugPlot[:debugPlotIdx, :])  # shape: (N, 4)

    # Filter out rows with any inf or nan values
    mask = np.all(np.isfinite(segments_np), axis=1)
    segments_np = segments_np[mask]

    # Reshape to (N, 2, 2): N segments, each with 2 endpoints (x1,y1), (x2,y2)
    segments = segments_np[:, :4].reshape(-1, 2, 2)

    # segments_np is (N, 4): [x1, y1, x2, y2] for each segment
    x1 = segments_np[:, 0]
    y1 = segments_np[:, 1]
    x2 = segments_np[:, 2]
    y2 = segments_np[:, 3]      
    color = segments_np[:, 4]

    # Find the segments of length larger than 100.0     
    length = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    mask = length > 100.0
    color[mask] = 0x0000FF
    if mask.any():
        res = False


    # Plot start points (x1, y1) and end points (x2, y2)
    plt.scatter(x1, y1, c=intToRGB(color), s=16)
    plt.scatter(x2, y2, c=intToRGB(color), s=16)

    # Create the LineCollection
    lc = LineCollection(segments, colors='g', linewidths=1.0, alpha=0.8)

    # Add to your axes (assuming ax is your matplotlib axes)
    ax.add_collection(lc)


    ax.autoscale()  # Adjust axes to fit the lines



    plt.show(block=False)
    plt.pause(0.2)
    # fig = None

    return res


def plotSideVar(varname, level=-1):
    global fig, ax
    
    # if fig is None: 
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Remove all existing plots
    ax.clear()

    # Get node coordinates for plotting
    x = getNodeVar("x")
    y = getNodeVar("y")
    
    # Plot all nodes
    plt.plot(x, y, 'ko', alpha=0.93, markersize=2.3)

    # Get side data
    n1 = getSideVar("n1").astype(int)
    n2 = getSideVar("n2").astype(int)
    sx = getSideVar("x")
    sy = getSideVar("y")
    
    # Get the variable to plot
    var = getSideVar(varname, level)
    
    # Build the segments array for LineCollection
    segments = np.stack([np.column_stack([x[n1], y[n1]]), np.column_stack([x[n2], y[n2]])], axis=1)

    # Create the LineCollection with colors based on the variable
    lc = LineCollection(segments, cmap='viridis', alpha=1.0, linewidths=0.5)
    lc.set_array(var)
    ax.add_collection(lc)
    
    # Add colorbar
    plt.colorbar(lc, label=varname)

    # Plot side centers with variable values
    scatter = plt.scatter(sx, sy, c=var, cmap='viridis', s=5, edgecolors='none')

    ax.autoscale()
    plt.title(f'Side Variable: {varname}')
    plt.show(block=True)
    # plt.pause(0.0)
    
    return True

def plotNodeVar(varname, level=-1):
    global fig, ax
    
    # Create or clear the figure
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.clear()
    
    # Get node coordinates for plotting
    x = getNodeVar("x")
    y = getNodeVar("y")
    
    # Get the variable to plot
    var = getNodeVar(varname, level)
    
    # Create scatter plot of nodes with colors based on the variable
    scatter = plt.scatter(x, y, c=var, cmap='viridis', s=5, 
                         edgecolors='none',  # Remove black edges
                         vmin=-0.1, vmax=0.1) 
    
    # Add colorbar
    plt.colorbar(scatter, label=varname)
    
    # Optional: Plot node connections if needed
    # If you want to show the mesh structure:
    n1 = getSideVar("n1").astype(int)
    n2 = getSideVar("n2").astype(int)
    segments = np.stack([np.column_stack([x[n1], y[n1]]), 
                        np.column_stack([x[n2], y[n2]])], axis=1)
    lc = LineCollection(segments, colors='gray', alpha=0.2, linewidths=0.5)
    ax.add_collection(lc)
    
    ax.autoscale()
    plt.title(f'Node Variable: {varname}')
    plt.show(block=True)
    
    return True