import os
import cupy as cp
import threading
import queue

import netCDF4 as nc
import debug.monitor as monitor 


path = None
nNodes = None
nElements = None
nSides = None
nLayers = None
sizeNode = None
sizeElem = None
sizeSide = None
dataC = None
dataU = None
dataV = None
dataZ = None
fileNum = 0


def initializeOutput(basePath, numNodes, numElements, numSides, numLayers, nodeSize, elemSize, sideSize):
    
    global path, nNodes, nElements, nSides, nLayers, sizeNode, sizeElem, sizeSide
    path = basePath
    nNodes = numNodes
    nElements = numElements
    nSides = numSides
    nLayers = numLayers
    sizeNode = nodeSize
    sizeElem = elemSize
    sizeSide = sideSize


output_thread = None
output_queue = queue.Queue()
output_lock = threading.Lock()

def writeOutput(it, core, opt):
    """Modified writeOutput with threading support."""
    global dataC, dataU, dataV, dataZ, fileNum, sizeNode, output_thread, output_queue, output_lock
    
    # Check if previous call is still working
    with output_lock:
        if output_thread is not None and output_thread.is_alive():
            output_thread.join()  # Wait until it finishes
    
    ihfskip = core.ihfskip
    nspool  = core.nspool
    dt      = core.dt
    timesPerFile = ihfskip // nspool

    if (it % nspool) != 0:
        return

    if (it-nspool) % ihfskip == 0:
        if dataC is not None:
            dataC.close()
            dataU.close()
            dataV.close()
            dataZ.close()

        fileNum += 1
        ifile=(it-1)/ihfskip+1
        dataC, dataU, dataV = createNetCDFFileHeaders(fileNum, nNodes, nElements, nSides, nLayers, opt)

    if dataC is None:
        return

    idxTime = ((it // nspool) - 1) % timesPerFile

    # Get all memory into CPU (this is the expensive part)
    output_data = {}
    
    # Horizontal geometric variables (only written once per file)
    if idxTime == 0:
        view = monitor.getNodeVar("x")
        output_data['xy_cpu'] = view.get()
        
        view = monitor.getSideVar("x")
        output_data['xy_side_cpu'] = view.get()
        
        view = monitor.getElementVar("x")
        output_data['xy_elem_cpu'] = view.get()
        
        view = monitor.getElementVar("FaceIndices")
        mask = (view[:,:] != 0)
        output_data['faceIndices_cpu'] = view.astype(cp.int64).get()
        output_data['firstPointer_face'] = view[mask].min().get()
        output_data['mask'] = mask.get()
        
        view = monitor.getSideVar("n1")
        output_data['n1_cpu'] = view.get()
        output_data['firstPointer_n1'] = view.min().get()
        
        view = monitor.getSideVar("n2")
        output_data['n2_cpu'] = view.get()
        output_data['firstPointer_n2'] = view.min().get()
        
        view = monitor.getNodeVar("depth")
        output_data['depth_cpu'] = view.get()

    # Time-dependent variables
    view = monitor.getNodeVar("eta")
    output_data['eta_cpu'] = view.get()
    
    view = monitor.getNodeVar("isdry")
    output_data['isdry_node_cpu'] = view.astype(cp.int32).get()
    
    view = monitor.getSideVar("isdry")
    output_data['isdry_side_cpu'] = view.astype(cp.int32).get()
    
    view = monitor.getElementVar("isdry")
    output_data['isdry_elem_cpu'] = view.astype(cp.int32).get()
    
    view = monitor.getNodeVar("u")
    output_data['u_cpu'] = view[0,:,:].get()
    output_data['v_cpu'] = view[1,:,:].get()
    
    view = monitor.getNodeVar("z")
    output_data['z_cpu'] = view.get()
    
    # Add metadata
    output_data['it'] = it
    output_data['core'] = core
    output_data['opt'] = opt
    output_data['idxTime'] = idxTime
    output_data['dt'] = dt
    output_data['timesPerFile'] = timesPerFile

    # Launch thread for writing (GPU can continue now)
    with output_lock:
        output_thread = threading.Thread(target=_writeOutputThread, args=(output_data,))
        output_thread.start()

def _writeOutputThread(output_data):
    """Background thread for actual file writing."""
    global dataC, dataU, dataV, dataZ, sizeNode
    
    try:
        it = output_data['it']
        core = output_data['core']
        opt = output_data['opt']
        idxTime = output_data['idxTime']
        dt = output_data['dt']
        timesPerFile = output_data['timesPerFile']
        
        # Horizontal geometric variables (only written once per file)
        if idxTime == 0:
            xy_cpu = output_data['xy_cpu']
            SCHISM_hgrid_node_x_var = dataC.variables['SCHISM_hgrid_node_x']
            SCHISM_hgrid_node_x_var[:] = xy_cpu[0, :]
            SCHISM_hgrid_node_y_var = dataC.variables['SCHISM_hgrid_node_y']
            SCHISM_hgrid_node_y_var[:] = xy_cpu[1, :]
            del xy_cpu

            xy_side_cpu = output_data['xy_side_cpu']
            SCHISM_hgrid_edge_x_var = dataC.variables['SCHISM_hgrid_edge_x']
            SCHISM_hgrid_edge_x_var[:] = xy_side_cpu[0, :]
            SCHISM_hgrid_edge_y_var = dataC.variables['SCHISM_hgrid_edge_y']
            SCHISM_hgrid_edge_y_var[:] = xy_side_cpu[1, :]
            del xy_side_cpu

            xy_elem_cpu = output_data['xy_elem_cpu']
            SCHISM_hgrid_face_x_var = dataC.variables['SCHISM_hgrid_face_x']
            SCHISM_hgrid_face_x_var[:] = xy_elem_cpu[0, :]
            SCHISM_hgrid_face_y_var = dataC.variables['SCHISM_hgrid_face_y']
            SCHISM_hgrid_face_y_var[:] = xy_elem_cpu[1, :]
            del xy_elem_cpu

            faceIndices_cpu = output_data['faceIndices_cpu']
            firstPointer = output_data['firstPointer_face']
            mask = output_data['mask']
            faceIndices_cpu[mask] = (faceIndices_cpu[mask] - firstPointer)//sizeNode + 1
            faceIndices_cpu[~mask] = -1
            nSCHISM_hgrid_face_nodes_var = dataC.variables['SCHISM_hgrid_face_nodes']
            nSCHISM_hgrid_face_nodes_var[:] = faceIndices_cpu.T
            del faceIndices_cpu, mask

            n1_cpu = output_data['n1_cpu']
            firstPointer = output_data['firstPointer_n1']
            n1_cpu[:] = (n1_cpu - firstPointer)//sizeNode + 1

            n2_cpu = output_data['n2_cpu']
            firstPointer = output_data['firstPointer_n2']
            n2_cpu[:] = (n2_cpu - firstPointer)//sizeNode + 1

            nSCHISM_hgrid_edge_nodes_var = dataC.variables['SCHISM_hgrid_edge_nodes']
            nSCHISM_hgrid_edge_nodes_var[:, 0] = n1_cpu
            nSCHISM_hgrid_edge_nodes_var[:, 1] = n2_cpu
            del n1_cpu, n2_cpu

            depth_cpu = output_data['depth_cpu']
            depth_var = dataC.variables['depth']
            depth_var[:] = depth_cpu
            del depth_cpu

            min_depth_var = dataC.variables['minimum_depth']
            min_depth_var[:] = 0.0001

        # Time-dependent variables
        eta_cpu = output_data['eta_cpu']
        elevation_var = dataC.variables['elevation']
        elevation_var[idxTime, :] = eta_cpu
        del eta_cpu

        isdry_node_cpu = output_data['isdry_node_cpu']
        dryFlagNode_var = dataC.variables['dryFlagNode']
        dryFlagNode_var[idxTime, :] = isdry_node_cpu
        del isdry_node_cpu

        isdry_side_cpu = output_data['isdry_side_cpu']
        dryFlagSide_var = dataC.variables['dryFlagSide']
        dryFlagSide_var[idxTime, :] = isdry_side_cpu
        del isdry_side_cpu

        isdry_elem_cpu = output_data['isdry_elem_cpu']
        dryFlagElement_var = dataC.variables['dryFlagElement']
        dryFlagElement_var[idxTime, :] = isdry_elem_cpu
        del isdry_elem_cpu

        u_cpu = output_data['u_cpu']
        v_cpu = output_data['v_cpu']
        horizontalVelX_var = dataU.variables['horizontalVelX']
        horizontalVelX_var[idxTime, :, :] = u_cpu.T
        horizontalVelY_var = dataV.variables['horizontalVelY']
        horizontalVelY_var[idxTime, :, :] = v_cpu.T
        del u_cpu, v_cpu

        z_cpu = output_data['z_cpu']
        zCoordinates_var = dataZ.variables['zCoordinates']
        zCoordinates_var[idxTime, :, :] = z_cpu.T
        del z_cpu

        # Time variables
        time_var = dataC.variables['time']
        time_var[idxTime] = it*dt
        time_var = dataU.variables['time']
        time_var[idxTime] = it*dt
        time_var = dataV.variables['time']
        time_var[idxTime] = it*dt
        time_var = dataZ.variables['time']
        time_var[idxTime] = it*dt

    finally:
        # Clean up all CPU memory at the end
        for key, value in output_data.items():
            if hasattr(value, 'shape'):  # Check if it's an array
                del value
        output_data.clear()



def closeOutputFiles():
    # Close the output files
    pass


def createNetCDFFileHeaders(fileNum, nNodes, nElements, nSides, nLayers, opt):
    """
    Create the NetCDF file "header", meaning the appropriate dimensions, variables, and attributes.
    The data variable is stored in the global variable data.
    """

    global dataC, dataU, dataV, dataZ
    
    filename = 'out2d_%i.nc' % fileNum
    filenameU = 'horizontalVelX_%i.nc' % fileNum
    filenameV = 'horizontalVelY_%i.nc' % fileNum
    filenameZ = 'zCoordinates_%i.nc' % fileNum

    # Remove the file if it exists before creating the NetCDF file
    filepath = os.path.join(path, filename)
    if os.path.exists(filepath):
        os.remove(filepath)

    filepathU = os.path.join(path, filenameU)
    if os.path.exists(filepathU):
        os.remove(filepathU)

    filepathV = os.path.join(path, filenameV)
    if os.path.exists(filepathV):
        os.remove(filepathV)

    filepathZ = os.path.join(path, filenameZ)
    if os.path.exists(filepathZ):
        os.remove(filepathZ)

    # Create the NetCDF files
    dataC = nc.Dataset(filepath,  "w", format="NETCDF4")
    dataU = nc.Dataset(filepathU, "w", format="NETCDF4")
    dataV = nc.Dataset(filepathV, "w", format="NETCDF4")
    dataZ = nc.Dataset(filepathZ, "w", format="NETCDF4")
    
    
    # Add dimensions to the file
    for data in [dataC, dataU, dataV, dataZ]:
        data.createDimension('nSCHISM_hgrid_node',          nNodes)
        data.createDimension('nSCHISM_hgrid_face',          nElements)
        data.createDimension('nSCHISM_hgrid_edge',          nSides)
        data.createDimension('nMaxSCHISM_hgrid_face_nodes', 4)
        data.createDimension('nSCHISM_vgrid_layers',        nLayers)
        data.createDimension('one',                         1)
        data.createDimension('two',                         2)
        data.createDimension('time',                        None)

        # Add global attributes
        data.title = "GSCHISM Output File"
        data.institution = "VIMS"
        data.source = "GSCHISM Model"
        data.history = "Created by GSCHISM"
        data.references = "SCHISM documentation"
        data.comment = "Output from GSCHISM simulation"

        # Add coordinate system attributes
        data.Conventions = "CF-1.6"
        data.coordinate_system = "Cartesian"
        
        # # Add time attributes
        # data.time_coverage_start = "2023-01-01T00:00:00Z"
        # data.time_coverage_end = "2023-01-02T00:00:00Z"

        # Create variables
        # Time variable
        time_var = data.createVariable('time', 'f8', ('time',))
        time_var.i23d = cp.int32(0)

        
        # Format the base_date string
        time_var.base_date = f" {opt.start_year:4d} {opt.start_month:2d} {opt.start_day:2d} {opt.start_hour:8.2f} {opt.utc_start:8.2f}"
        time_var.units = f"seconds since {opt.start_year:04d}-{opt.start_month:02d}-{opt.start_day:02d}T{opt.start_hour:02.0f}:00:00+00"
        time_var.standard_name = "time"
        time_var.axis = "T"
    

    # Minimum depth
    min_depth_var = dataC.createVariable('minimum_depth', 'f8', ('one',))
    min_depth_var.units = "m"
    
    # SCHISM grid topology
    schism_grid_var = dataC.createVariable('SCHISM_hgrid', 'S1', ('one',))
    schism_grid_var.long_name = "Topology data of 2d unstructured mesh"
    schism_grid_var.topology_dimension = 2
    schism_grid_var.cf_role = "mesh_topology"
    schism_grid_var.node_coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    schism_grid_var.edge_coordinates = "SCHISM_hgrid_edge_x SCHISM_hgrid_edge_y"
    schism_grid_var.face_coordinates = "SCHISM_hgrid_face_x SCHISM_hgrid_face_y"
    schism_grid_var.edge_node_connectivity = "SCHISM_hgrid_edge_nodes"
    schism_grid_var.face_node_connectivity = "SCHISM_hgrid_face_nodes"
    
    # Node coordinates
    node_x_var = dataC.createVariable('SCHISM_hgrid_node_x', 'f8', ('nSCHISM_hgrid_node',))
    node_x_var.axis = "X"
    node_x_var.location = "node"
    node_x_var.mesh = "SCHISM_hgrid"
    node_x_var.units = "m"
    node_x_var.standard_name = "projection_x_coordinate"
    
    node_y_var = dataC.createVariable('SCHISM_hgrid_node_y', 'f8', ('nSCHISM_hgrid_node',))
    node_y_var.axis = "Y"
    node_y_var.location = "node"
    node_y_var.mesh = "SCHISM_hgrid"
    node_y_var.units = "m"
    node_y_var.standard_name = "projection_y_coordinate"
    
    # Depth and bottom index
    depth_var = dataC.createVariable('depth', 'f4', ('nSCHISM_hgrid_node',))
    depth_var.units = "m"
    depth_var.axis = "Z"
    depth_var.positive = "down"
    depth_var.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    depth_var.location = "node"
    depth_var.grid_mapping = "crs"
    depth_var.mesh = "SCHISM_hgrid"
    
    bottom_idx_var = dataC.createVariable('bottom_index_node', 'i4', ('nSCHISM_hgrid_node',))
    bottom_idx_var.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    bottom_idx_var.location = "node"
    bottom_idx_var.grid_mapping = "crs"
    bottom_idx_var.mesh = "SCHISM_hgrid"
    
    # Face coordinates
    face_x_var = dataC.createVariable('SCHISM_hgrid_face_x', 'f8', ('nSCHISM_hgrid_face',))
    face_x_var.axis = "X"
    face_x_var.location = "face"
    face_x_var.mesh = "SCHISM_hgrid"
    face_x_var.units = "m"
    face_x_var.standard_name = "projection_x_coordinate"
    
    face_y_var = dataC.createVariable('SCHISM_hgrid_face_y', 'f8', ('nSCHISM_hgrid_face',))
    face_y_var.axis = "Y"
    face_y_var.location = "face"
    face_y_var.mesh = "SCHISM_hgrid"
    face_y_var.units = "m"
    face_y_var.standard_name = "projection_y_coordinate"
    
    # Edge coordinates
    edge_x_var = dataC.createVariable('SCHISM_hgrid_edge_x', 'f8', ('nSCHISM_hgrid_edge',))
    edge_x_var.axis = "X"
    edge_x_var.location = "edge"
    edge_x_var.mesh = "SCHISM_hgrid"
    edge_x_var.units = "m"
    edge_x_var.standard_name = "projection_x_coordinate"
    
    edge_y_var = dataC.createVariable('SCHISM_hgrid_edge_y', 'f8', ('nSCHISM_hgrid_edge',))
    edge_y_var.axis = "Y"
    edge_y_var.location = "edge"
    edge_y_var.mesh = "SCHISM_hgrid"
    edge_y_var.units = "m"
    edge_y_var.standard_name = "projection_y_coordinate"
    
    # Connectivity arrays
    face_nodes_var = dataC.createVariable('SCHISM_hgrid_face_nodes', 'i4', ('nSCHISM_hgrid_face', 'nMaxSCHISM_hgrid_face_nodes'), fill_value=-1)
    face_nodes_var.start_index = 1
    face_nodes_var.cf_role = "face_node_connectivity"
    
    edge_nodes_var = dataC.createVariable('SCHISM_hgrid_edge_nodes', 'i4', ('nSCHISM_hgrid_edge', 'two'), fill_value=-1)
    edge_nodes_var.start_index = 1
    edge_nodes_var.cf_role = "edge_node_connectivity"
    
    # Time-dependent variables
    dry_flag_node_var = dataC.createVariable('dryFlagNode', 'f4', ('time', 'nSCHISM_hgrid_node'))
    dry_flag_node_var.i23d = cp.int32(1)
    dry_flag_node_var.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    dry_flag_node_var.location = "node"
    dry_flag_node_var.grid_mapping = "crs"
    dry_flag_node_var.mesh = "SCHISM_hgrid"
    
    elevation_var = dataC.createVariable('elevation', 'f4', ('time', 'nSCHISM_hgrid_node'))
    elevation_var.i23d = cp.int32(1)
    elevation_var.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    elevation_var.location = "node"
    elevation_var.grid_mapping = "crs"
    elevation_var.mesh = "SCHISM_hgrid"
    
    wind_speed_x_var = dataC.createVariable('windSpeedX', 'f4', ('time', 'nSCHISM_hgrid_node'))
    wind_speed_x_var.i23d = cp.int32(1)
    wind_speed_x_var.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    wind_speed_x_var.location = "node"
    wind_speed_x_var.grid_mapping = "crs"
    wind_speed_x_var.mesh = "SCHISM_hgrid"
    
    wind_speed_y_var = dataC.createVariable('windSpeedY', 'f4', ('time', 'nSCHISM_hgrid_node'))
    wind_speed_y_var.i23d = cp.int32(1)
    wind_speed_y_var.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    wind_speed_y_var.location = "node"
    wind_speed_y_var.grid_mapping = "crs"
    wind_speed_y_var.mesh = "SCHISM_hgrid"
    
    dry_flag_element_var = dataC.createVariable('dryFlagElement', 'f4', ('time', 'nSCHISM_hgrid_face'))
    dry_flag_element_var.i23d = cp.int32(4)
    dry_flag_element_var.coordinates = "SCHISM_hgrid_face_x SCHISM_hgrid_face_y"
    dry_flag_element_var.location = "face"
    dry_flag_element_var.grid_mapping = "crs"
    dry_flag_element_var.mesh = "SCHISM_hgrid"
    
    dry_flag_side_var = dataC.createVariable('dryFlagSide', 'f4', ('time', 'nSCHISM_hgrid_edge'))
    dry_flag_side_var.i23d = cp.int32(7)
    dry_flag_side_var.coordinates = "SCHISM_hgrid_edge_x SCHISM_hgrid_edge_y"
    dry_flag_side_var.location = "edge"
    dry_flag_side_var.grid_mapping = "crs"
    dry_flag_side_var.mesh = "SCHISM_hgrid"


    horizontalVelX = dataU.createVariable('horizontalVelX', 'f4', ('time', 'nSCHISM_hgrid_node', 'nSCHISM_vgrid_layers'))
    horizontalVelX.i23d = cp.int32(2)
    horizontalVelX.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    horizontalVelX.location = "node"
    horizontalVelX.grid_mapping = "crs"
    horizontalVelX.mesh = "SCHISM_hgrid" 

    horizontalVelY = dataV.createVariable('horizontalVelY', 'f4', ('time', 'nSCHISM_hgrid_node', 'nSCHISM_vgrid_layers'))
    horizontalVelY.i23d = cp.int32(2)
    horizontalVelY.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    horizontalVelY.location = "node"
    horizontalVelY.grid_mapping = "crs"
    horizontalVelY.mesh = "SCHISM_hgrid"    


    # Add zCoordinates variable to dataZ file
    zCoordinates = dataZ.createVariable('zCoordinates', 'f4', ('time', 'nSCHISM_hgrid_node', 'nSCHISM_vgrid_layers'))
    zCoordinates.i23d = cp.int32(2)
    zCoordinates.missing_value = 9.96921e+36
    zCoordinates.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    zCoordinates.location = "node"
    zCoordinates.grid_mapping = "crs"
    zCoordinates.mesh = "SCHISM_hgrid"

    
    return dataC, dataU, dataV

    




