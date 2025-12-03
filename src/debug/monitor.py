from python.schism_consts import *
import python.cuda_functions as cufunc
import cupy as cp
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
import matplotlib

if not VISTA:
    import paramiko
    import re

# matplotlib.use('TkAgg') 
x, y, elem_x, elem_y,  sx, sy = None, None, None, None, None, None
nodesInfo = [40798, 41152, 41486, 40796, 41150, 41484, 40797, 41151, 41485]
sidesInfo = [120390, 121448, 121447, 121450, 122462, 121446, 121452]
elemsInfo = [69889, 69888, 70516, 69891, 70514, 70518]
fortranDebuggerInitialized = False
chan = None
continue_execution = False
plotCreated = False

debugPlot = cp.zeros([10000,6], dtype=cp.float64, order='F') 
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
                strides = (8, 16, sizeSide)  # 8 bytes for double, then 8*nvrt for next level, then stride for next element
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

        if k != -1:
            view = view[:, k, :]

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

  
    raise ValueError(f"Variable {varName} is not mapped in memory")


# XXXX REMOVE THESE FUNCTIONS and the corresponding ones in the kernels


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


def setupFortranDebugger():
    global fortranDebuggerInitialized, chan

    if fortranDebuggerInitialized:
        return
    fortranDebuggerInitialized = True
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect("femto.sciclone.wm.edu", username="jgonzalezondin", password="Zam3loide2025!" )
    chan = ssh.invoke_shell()  # Create chan FIRST
    
    # Send srun command and wait for it to complete
    chan.send("srun --pty -p batch -n 1 -t 02:00:00 bash -l\n")
    # Wait for srun to allocate resources and start the interactive session
    import time
    time.sleep(2)  # Give it a moment to start
    # Wait for prompt (you may need to adjust this based on actual output)
    output = ""
    while "bash" not in output.lower() and len(output) < 1000:
        if chan.recv_ready():
            output += chan.recv(1024).decode('utf-8', errors='replace')
        else:
            time.sleep(0.1)
    
    # Now send the unset command
    chan.send("unset PMI_RANK PMI_SIZE PMI_FD PMIX_NAMESPACE PMIX_RANK PMIX_SERVER_URI2\n")
    time.sleep(0.5)  # Brief wait for command to execute
    _ = chan.recv(1024).decode('utf-8', errors='replace')  # Clear any output
    
    chan.send("cd /sciclone/home/jgonzalezondin/data10/RUN02d\n")
    time.sleep(0.5)
    _ = chan.recv(1024).decode('utf-8', errors='replace')
    
    chan.send("gdb -q -i=mi ./pschism_FEMTO_OLDIO_BLD_STANDALONE_SH_MEM_COMM_TVD-VL\n")

    # Wait for gdb prompt and show output
    output = chan.recv(2048).decode()
    print("=== Initial gdb startup ===")
    print(output)

    # Attach or run as needed
    # Set breakpoint first
    # chan.send("-break-insert solver_subs.F90:316\n") # end solver
    chan.send("-break-insert schism_step.F90:4505\n") # after btrack
    # chan.send("-break-insert schism_step.F90:6821\n") # Before shapiro
    # chan.send("-break-insert schism_step.F90:6887\n") # After shapiro
    # chan.send("-break-insert schism_step.F90:8271\n") # Before levels0
    # chan.send("-break-insert schism_step.F90:8276\n") # After levels0
    # Receive confirmation
    _ = chan.recv(4096).decode('utf-8', errors='replace')
    # Add condition
    chan.send(f"-break-condition 1 it .eq. {1}\n")
    _ = chan.recv(4096).decode('utf-8', errors='replace')

    chan.send("-exec-run\n")

    chan.send("-gdb-set print elements 0\n")  # 0 means unlimited
    chan.send("-gdb-set print repeats 0\n")     # Don't show repeated values

    # After the breakpoint hits, increase the max value size
    chan.send("-gdb-set max-value-size unlimited\n")
    # Or set it to a large number
    chan.send("-gdb-set max-value-size 10000000\n")
    # Wait for confirmation
    _ = chan.recv(4096).decode('utf-8', errors='replace')


def fortranDebugPlot(it):

    global fig, ax, plotCreated,x, y, elem_x, elem_y,  sx, sy, continue_execution

    setupFortranDebugger()

    # Set the break condition
    chan.send(f"-break-condition 1 it .eq. {it}\n")
    _ = chan.recv(4096).decode('utf-8', errors='replace')
    # Now continue execution
    chan.send("-exec-continue\n")
    _ = chan.recv(4096).decode('utf-8', errors='replace')

    res = True
    
    if fig is None: 
        fig, ax = plt.subplots(figsize=(12, 8))
    
    # # Remove all existing plots
    # if plotCreated:
    #     ax.clear()

    # Connect the click event handler
    if not hasattr(fig, '_click_handler_connected'):
        fig.canvas.mpl_connect('button_press_event', on_click)
        fig._click_handler_connected = True

    # Connect the key press event handler
    if not hasattr(fig, '_key_handler_connected'):
        fig.canvas.mpl_connect('key_press_event', on_key_press)
        fig._key_handler_connected = True        

   # Get node data
    isNodeDry = getNodeVar("isdry").get() == 1.0
    isSegDry = getSideVar("isdry").get()
    
    node_x_var = getNodeVar("x")
    x = node_x_var[0, :].get()
    y = node_x_var[1, :].get()
    
    nu = getNodeVar("u").get()
    nux = nu[0, 1, :]
    nuy = nu[1, 1, :]

    # Get element center coordinates and dry status
    elem_x_var = getElementVar("x")
    elem_x = elem_x_var[0, :].get()
    elem_y = elem_x_var[1, :].get()
    elem_isDry = getElementVar("isdry").get() == 1.0

    # Plot element centers with different colors based on dry/wet status
    elem_dry_mask = elem_isDry.astype(bool)
    elem_wet_mask = ~elem_dry_mask

    # Get side data
    n1 = getSideVar("n1").astype(int).get()
    n2 = getSideVar("n2").astype(int).get()

    n1 = (n1 - n1.min()) 
    tmp = np.sort(np.unique(n1))
    tmplen = tmp[1] - tmp[0]   # Size of the node structure in bytes
    n1 = n1 // tmplen
    n2 = (n2 - n2.min()) // tmplen
    
    side_x_var = getSideVar("x")
    sx = side_x_var[0, :].get()
    sy = side_x_var[1, :].get()
    
    side_u_var = getSideVar("U")
    ux = side_u_var[0, :].get()
    uy = side_u_var[1, :].get()

    side_sdbt_var = getSideVar("btu")
    sbux = side_sdbt_var[0, 0, :].get()
    sbuy = side_sdbt_var[1, 0,:].get()    

    side_su_var = getSideVar("u")
    sux = side_su_var[0, 0,:].get()
    suy = side_su_var[1, 0,:].get()       
    
    eta = getNodeVar("eta").get()

    if not plotCreated:
        # Build the segments array for LineCollection
        # Each segment needs to be shape (2, 2): [[x1, y1], [x2, y2]]
        p1 = np.column_stack([x[n1], y[n1]])  # Start points: shape (N, 2)
        p2 = np.column_stack([x[n2], y[n2]])  # End points: shape (N, 2)
        segments = np.stack([p1, p2], axis=1)  # Shape (N, 2, 2)

        # Plots eta at the element centers. Eta is shown as the radius of the circle.
        # plt.scatter(x, y, color='b', s=(eta*200 + 20))
        # Plot dry nodes in gray
        if np.any(isNodeDry):
            plt.plot(x[isNodeDry], y[isNodeDry], 'ro', alpha=0.93, markersize=2.3)
            # plt.scatter(x[isNodeDry], y[isNodeDry], c='gray', s=3, alpha=0.8, label='Dry nodes')
        
        # Plot wet nodes in blue
        if np.any(~isNodeDry):
            plt.plot(x[~isNodeDry], y[~isNodeDry], 'bo', alpha=0.93, markersize=2.3)
            # plt.scatter(x[~isNodeDry], y[~isNodeDry], c='blue', s=3, alpha=0.8, label='Wet nodes')    


        # Plot dry element centers as gray 'x'
        if np.any(elem_dry_mask):
            plt.scatter(elem_x[elem_dry_mask], elem_y[elem_dry_mask], 
                    c='darkgray', marker='x', s=20, alpha=0.8, label='Dry elements')
        
        # Plot wet element centers as red 'x'
        if np.any(elem_wet_mask):
            plt.scatter(elem_x[elem_wet_mask], elem_y[elem_wet_mask], 
                    c='red', marker='x', s=20, alpha=0.8, label='Wet elements')        

        # Create colors for segments based on dry/wet status
        segment_colors = np.where(isSegDry, 'gray', 'r')  # gray for dry, red for wet

        # Create the LineCollection
        lc = LineCollection(segments, colors=segment_colors, alpha=0.53, linewidths=0.5)
        ax.add_collection(lc)

        # # Use a fixed scale factor for consistent arrow length representation
        # # Adjust this value to make arrows longer or shorter while maintaining 1 m/s = constant length
        # fixed_scale = 3.001  # Adjust this value as needed for your visualization
        
        # plt.quiver(sx, sy, fixed_scale*ux, fixed_scale*uy, angles='xy', scale_units='inches', scale=1, color='b', width=0.003)
        # plt.quiver(x, y, fixed_scale*nux, fixed_scale*nuy, angles='xy', scale_units='inches', scale=1, color='r', width=0.003)

        showLabels = False
        
        # Plot labels if requested
        if showLabels:
            # Get indices
            for i in range(len(x)):
                ax.text(x[i], y[i], f'N{i}', fontsize=6, color='black', 
                    ha='center', va='bottom', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
            
            # Plot element indices
            for i in range(len(elem_x)):
                ax.text(elem_x[i], elem_y[i], f'E{i}', fontsize=6, color='darkred', 
                    ha='center', va='top', bbox=dict(boxstyle='round,pad=0.2', facecolor='yellow', alpha=0.5))
            
            # Plot side indices (at side centers)
            for i in range(len(sx)):
                ax.text(sx[i], sy[i], f'S{i}', fontsize=5, color='darkblue', 
                    ha='center', va='center', bbox=dict(boxstyle='round,pad=0.15', facecolor='cyan', alpha=0.4))

    
        ax.autoscale()  # Adjust axes to fit the lines



        plt.show(block=False)
        plotCreated = True


    import time
    while not continue_execution and plt.get_fignums():
        fig.canvas.flush_events()  # Process pending events
        time.sleep(0.01) 

    # Print the debugging information
    print(f"Debugging information for iteration {it}:")
    for nodeIdx in nodesInfo:

        u = evaluate_gdb_expression(f"uu2(1, {nodeIdx+1})")
        v = evaluate_gdb_expression(f"vv2(1, {nodeIdx+1})")
        eta2 = evaluate_gdb_expression(f"eta2({nodeIdx+1})")
        isdry = evaluate_gdb_expression(f"idry({nodeIdx+1})")
        print(f"Node {nodeIdx}: x = {x[nodeIdx]}, y = {y[nodeIdx]}, isdry = {isdry}-{int(isNodeDry[nodeIdx])}")        
        print(f"    u2 = ({u}, {v}) - ({nux[nodeIdx]}, {nuy[nodeIdx]})")
        print(f"    eta2 = {eta2}, {eta[nodeIdx]}")
    print("--------------------------------")
    for elemIdx in elemsInfo:
        isdry = evaluate_gdb_expression(f"idry_e({elemIdx+1})")
        print(f"Element {elemIdx}: x = {elem_x[elemIdx]}, y = {elem_y[elemIdx]}, isdry = {isdry}-{int(elem_isDry[elemIdx])}")
    print("--------------------------------")
    for sideIdx in sidesInfo:
        u = evaluate_gdb_expression(f"su2(1, {sideIdx+1})")
        v = evaluate_gdb_expression(f"sv2(1, {sideIdx+1})")
        sdbtu = evaluate_gdb_expression(f"sdbt(1,1, {sideIdx+1})")
        sdbtv = evaluate_gdb_expression(f"sdbt(2,1, {sideIdx+1})")        
        isdry = evaluate_gdb_expression(f"idry_s({sideIdx+1})")
        print(f"Side {sideIdx}: x = {sx[sideIdx]}, y = {sy[sideIdx]}, isdry = {isdry}-{int(isSegDry[sideIdx])}")
        print(f"    u    = ({u}, {v}) - ({sux[sideIdx]}, {suy[sideIdx]}")
        print(f"    sdbt = ({sdbtu}, {sdbtv}) - ({sbux[sideIdx]}, {sbuy[sideIdx]})")
    print("--------------------------------")
    continue_execution = False
    # fig = None

    return res    


def on_key_press(event):
    """Clear all selected points when delete key is pressed."""
    global nodesInfo, elemsInfo, sidesInfo, continue_execution 
    
    # Check if delete key was pressed
    if event.key == 'delete' or event.key == 'del':
        # Clear all selected items
        nodesInfo.clear()
        elemsInfo.clear()
        sidesInfo.clear()
        
        # Remove all highlight markers
        if hasattr(ax, '_highlight_markers'):
            for marker in ax._highlight_markers:
                marker.remove()
            ax._highlight_markers.clear()
        
        # Redraw the plot
        if fig is not None:
            fig.canvas.draw()

    elif event.key == 'c':
        continue_execution = True


def on_click(event):
    """Display information about the closest node, side, or element when clicking."""
    if event.inaxes != ax:
        return
    
    click_x, click_y = event.xdata, event.ydata
    
    # Calculate distances to nodes
    node_distances = np.sqrt((x - click_x)**2 + (y - click_y)**2)
    closest_node_idx = np.argmin(node_distances)
    node_dist = node_distances[closest_node_idx]
    
    # Calculate distances to element centers
    elem_distances = np.sqrt((elem_x - click_x)**2 + (elem_y - click_y)**2)
    closest_elem_idx = np.argmin(elem_distances)
    elem_dist = elem_distances[closest_elem_idx]
    
    # Calculate distances to side centers
    side_distances = np.sqrt((sx - click_x)**2 + (sy - click_y)**2)
    closest_side_idx = np.argmin(side_distances)
    side_dist = side_distances[closest_side_idx]
    
    # Find the overall closest
    min_dist = min(node_dist, elem_dist, side_dist)
    
    if min_dist == node_dist:
        closest_type = "Node"
        closest_idx = closest_node_idx
        if closest_idx not in nodesInfo:
            nodesInfo.append(closest_idx)
            
    elif min_dist == elem_dist:
        closest_type = "Element"
        closest_idx = closest_elem_idx
        if closest_idx not in elemsInfo:
            elemsInfo.append(closest_idx)
    else:
        closest_type = "Side"
        closest_idx = closest_side_idx
        if closest_idx not in sidesInfo:
            sidesInfo.append(closest_idx)
    

    
    # Plot all selected nodes, elements, and sides, removing the previously plotted ones
    # Remove previous highlight markers if they exist
    if hasattr(ax, '_highlight_markers'):
        for marker in ax._highlight_markers:
            marker.remove()
        ax._highlight_markers.clear()
    else:
        ax._highlight_markers = []
    
    # Plot new markers and track them
    for node in nodesInfo:
        marker = ax.plot(x[node], y[node], 'y*', markersize=15, markeredgecolor='k', markeredgewidth=1)[0]
        ax._highlight_markers.append(marker)
    for elem in elemsInfo:
        marker = ax.plot(elem_x[elem], elem_y[elem], 'r*', markersize=15, markeredgecolor='k', markeredgewidth=1)[0]
        ax._highlight_markers.append(marker)
    for side in sidesInfo:
        marker = ax.plot(sx[side], sy[side], 'b*', markersize=15, markeredgecolor='k', markeredgewidth=1)[0]
        ax._highlight_markers.append(marker)


    fig.canvas.draw()


def plotMesh():
    global fig, ax

    res = True
    
    if fig is None: 
        fig, ax = plt.subplots(figsize=(12, 8))
    
    # Remove all existing plots
    ax.clear()

    isNodeDry = getNodeVar("isdry").get() == 1.0
    isSegDry = getSideVar("isdry").get()
    x = getNodeVar("x")
    y = x[1, :].get()
    x = x[0, :].get()
    nu = getNodeVar("u").get()
    nux = nu[0, 1, :]
    nuy = nu[1, 1, :]



    # Get element center coordinates and dry status
    elem_x = getElementVar("x")
    elem_y = elem_x[1, :].get()
    elem_x = elem_x[0, :].get()
    elem_isDry = getElementVar("isdry").get() == 1.0

    # Plot element centers with different colors based on dry/wet status
    elem_dry_mask = elem_isDry.astype(bool)
    elem_wet_mask = ~elem_dry_mask
        
        
        # Plot all nodes
    # plt.plot(x, y, 'ko', alpha=0.93, markersize=2.3)

    n1 = getSideVar("n1").astype(int).get()
    n2 = getSideVar("n2").astype(int).get()
    n1 = (n1 - n1.min())//1584
    n2 = (n2 - n2.min())//1584
    sx = getSideVar("x")
    sy = sx[1, :].get()
    sx = sx[0, :].get()
    u = getSideVar("U")
    ux = u[0, :].get()
    uy = u[1, :].get()
    eta = getNodeVar("eta").get()



    # Build the segments array for LineCollection
    segments = np.stack([np.column_stack([x[n1], y[n1]]), np.column_stack([x[n2], y[n2]])], axis=1)

    # Plots eta at the element centers. Eta is shown as the radius of the circle.
    # plt.scatter(x, y, color='b', s=(eta*200 + 20))
    # Plot dry nodes in gray
    if np.any(isNodeDry):
        plt.plot(x[isNodeDry], y[isNodeDry], 'ro', alpha=0.93, markersize=2.3)
        # plt.scatter(x[isNodeDry], y[isNodeDry], c='gray', s=3, alpha=0.8, label='Dry nodes')
    
    # Plot wet nodes in blue
    if np.any(~isNodeDry):
        plt.plot(x[~isNodeDry], y[~isNodeDry], 'bo', alpha=0.93, markersize=2.3)
        # plt.scatter(x[~isNodeDry], y[~isNodeDry], c='blue', s=3, alpha=0.8, label='Wet nodes')    


    # Plot dry element centers as gray 'x'
    if np.any(elem_dry_mask):
        plt.scatter(elem_x[elem_dry_mask], elem_y[elem_dry_mask], 
                   c='darkgray', marker='x', s=20, alpha=0.8, label='Dry elements')
    
    # Plot wet element centers as red 'x'
    if np.any(elem_wet_mask):
        plt.scatter(elem_x[elem_wet_mask], elem_y[elem_wet_mask], 
                   c='red', marker='x', s=20, alpha=0.8, label='Wet elements')        

    # Create colors for segments based on dry/wet status
    segment_colors = np.where(isSegDry, 'gray', 'r')  # gray for dry, red for wet

    # Create the LineCollection
    lc = LineCollection(segments, colors=segment_colors, alpha=0.53, linewidths=0.5)
    ax.add_collection(lc)

    # Use a fixed scale factor for consistent arrow length representation
    # Adjust this value to make arrows longer or shorter while maintaining 1 m/s = constant length
    fixed_scale = 3.001  # Adjust this value as needed for your visualization
    
    plt.quiver(sx, sy, fixed_scale*ux, fixed_scale*uy, angles='xy', scale_units='inches', scale=1, color='b', width=0.003)
    plt.quiver(x, y, fixed_scale*nux, fixed_scale*nuy, angles='xy', scale_units='inches', scale=1, color='r', width=0.003)

    showLabels = True
    
    # Plot labels if requested
    if showLabels:
        # Get indices
        for i in range(len(x)):
            ax.text(x[i], y[i], f'N{i}', fontsize=6, color='black', 
                   ha='center', va='bottom', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
        
        # Plot element indices
        for i in range(len(elem_x)):
            ax.text(elem_x[i], elem_y[i], f'E{i}', fontsize=6, color='darkred', 
                   ha='center', va='top', bbox=dict(boxstyle='round,pad=0.2', facecolor='yellow', alpha=0.5))
        
        # Plot side indices (at side centers)
        for i in range(len(sx)):
            ax.text(sx[i], sy[i], f'S{i}', fontsize=5, color='darkblue', 
                   ha='center', va='center', bbox=dict(boxstyle='round,pad=0.15', facecolor='cyan', alpha=0.4))

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
    marker = segments_np[:, 5]
    # Find the segments of length larger than 100.0     
    length = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    mask = length > 100.0
    # color[mask] = 0x0000FF
    if mask.any():
        res = False


    # Plot start points (x1, y1) and end points (x2, y2)
    unique_markers = np.unique(marker)
    for m in unique_markers:
        mask = marker == m
        plt.scatter(x1[mask], y1[mask], c=intToRGB(color[mask]), 
                    s=16, marker=chr(m.astype(int)))
    for m in unique_markers:
        mask = marker == m
        plt.scatter(x2[mask], y2[mask], c=intToRGB(color[mask]), 
                    s=16, marker=chr(m.astype(int)))

    # Create the LineCollection
    lc = LineCollection(segments, colors='g', linewidths=1.0, alpha=0.8)

    # Add to your axes (assuming ax is your matplotlib axes)
    ax.add_collection(lc)


    ax.autoscale()  # Adjust axes to fit the lines



    plt.show(block=False)
    plt.pause(0.0)
    fig = None

    return res


def plotSideVar(varname, level=-1):
    global fig, ax
    
    # if fig is None: 
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Remove all existing plots
    ax.clear()

    # Get node coordinates for plotting
    x = getNodeVar("x")[0,:].get()
    y = getNodeVar("x")[1,:].get()
    
    # Plot all nodes
    ax.plot(x, y, 'ko', alpha=0.93, markersize=2.3)

    # Get side data
    n1 = getSideVar("n1").astype(int).get()
    n2 = getSideVar("n2").astype(int).get()
    sx = getSideVar("x")[0,:].get()
    sy = getSideVar("x")[1,:].get()

    n1 = (n1 - n1.min())//1584 # XXXXX fix
    n2 = (n2 - n2.min())//1584
    
    # Get the variable to plot
    var = getSideVar(varname, level).get()
    if var.ndim == 2:
        var = np.sqrt(var[0,:]**2 + var[1,:]**2)
    
    # Build the segments array for LineCollection
    segments = np.stack([np.column_stack([x[n1], y[n1]]), np.column_stack([x[n2], y[n2]])], axis=1)

    # Create the LineCollection with colors based on the variable
    lc = LineCollection(segments, cmap='viridis', alpha=1.0, linewidths=0.5)
    lc.set_array(var)
    ax.add_collection(lc)
    
    # Add colorbar
    plt.colorbar(lc, label=varname)

    # Plot side centers with variable values
    scatter = ax.scatter(sx, sy, c=var, cmap='viridis', s=5, edgecolors='none')

    ax.autoscale()
    plt.title(f'Side Variable: {varname}')
    plt.show(block=True)
    plt.pause(0.0)
    
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


def evaluate_gdb_expression(expression, timeout=2.0):
    """
    Evaluate an expression in GDB and return the result.
    
    Args:
        expression: The expression to evaluate (e.g., "uu2(123)" or "it")
        timeout: Maximum time to wait for response in seconds
    
    Returns:
        str: The value returned by GDB, or None if error/timeout
    """
    global chan
    
    if chan is None:
        print("Error: GDB channel not initialized")
        return None

    # Wait for and receive the response
    import time
    import time
    while chan.recv_ready():
        _ = chan.recv(4096).decode('utf-8', errors='replace')
    time.sleep(0.01)  # Brief pause to ensure buffer is clear
    
    # Send the evaluation command
    chan.send(f'-data-evaluate-expression "{expression}"\n')
    

    response = ""
    start_time = time.time()
    
    while time.time() - start_time < timeout:
        if chan.recv_ready():
            chunk = chan.recv(4096).decode('utf-8', errors='replace')
            response += chunk
            # GDB MI responses end with ^done or ^error
            if '^done' in response or '^error' in response:
                break
        else:
            time.sleep(0.02)
    
    # Parse GDB MI format: ^done,value="<value>"
    if '^done' in response:
        match = re.search(r'value="([^"]+)"', response)
        if match:
            return match.group(1)
        else:
            print(f"Warning: Could not parse GDB response: {response[:200]}")
            return None
    elif '^error' in response:
        match = re.search(r'msg="([^"]+)"', response)
        error_msg = match.group(1) if match else "Unknown error"
        print(f"GDB Error evaluating '{expression}': {error_msg}")
        return None
    else:
        print(f"Warning: No valid response from GDB for '{expression}': {response[:200]}")
        return None