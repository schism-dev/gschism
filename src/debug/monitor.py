from python.schism_consts import *
import python.cuda_functions as cufunc
import cupy as cp
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
import matplotlib

matplotlib.use('TkAgg') 



debugPlot = cp.zeros([10000,5], dtype=cp.float64, order='F') 
debugPlotIdx = cp.zeros(1, dtype=cp.int32, order='F') 
debugPlotIdx[0] = 0


fig = None
ax  = None

# Create dictionary mapping variable names to indices
varCodeDict = {var: idx for idx, var in enumerate(varCodes.split(','))}

numNodes = 0
numElements = 0
numSides = 0

def init(nNodes, nElements, nSides):
    global numNodes, numElements, numSides
    numNodes = nNodes
    numElements = nElements
    numSides = nSides

def getElementVar(varname, k = -1, j = -1):
    varCode = varCodeDict[varname]
    var = cp.zeros(numElements, dtype=cp.float64)
    cufunc.getElementVar(varCode, var, k, j)
    return var.get()

def getNodeVar(varname, k = -1):
    varCode = varCodeDict[varname]
    var = cp.zeros(numNodes, dtype=cp.float64)
    cufunc.getNodeVar(varCode, var, k)
    return var.get()

def getSideVar(varname, k = -1):
    varCode = varCodeDict[varname]
    var = cp.zeros(numSides, dtype=cp.float64)
    cufunc.getSideVar(varCode, var, k)
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
