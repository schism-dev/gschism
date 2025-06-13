from python.schism_consts import *
import python.cuda_functions as cufunc
import cupy as cp
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection


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

def getElementVar(varname, k = -1):
    varCode = varCodeDict[varname]
    var = cp.zeros(numElements, dtype=cp.float64)
    cufunc.getElementVar(varCode, var, k)
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




def plotMesh():
    fig, ax = plt.subplots(figsize=(12, 8))

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

    # Build the segments array for LineCollection
    segments = np.stack([np.column_stack([x[n1], y[n1]]), np.column_stack([x[n2], y[n2]])], axis=1)

    # Create the LineCollection
    lc = LineCollection(segments, colors='r', alpha=0.53, linewidths=0.5)
    ax.add_collection(lc)

    max_u = max(np.max(np.abs(ux)), np.max(np.abs(uy)))
    plt.quiver(sx, sy, 1*ux/max_u, 1*uy/max_u, angles='xy', scale_units='inches', scale=1, color='b', width=0.003)

    ax.autoscale()  # Adjust axes to fit the lines



    plt.show()
    plt.pause(0.0)
