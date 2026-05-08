import os
import cupy as cp
import numpy as np
import netCDF4 as nc
import debug.monitor as monitor 
import time


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
dataW = None
dataZ = None
fileNum = 0
outputTime = 0.0

paramsCompression = {'zlib': False}


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





def writeOutput(it, core, opt):
    # Write the output files
    global dataC, dataU, dataV, dataW, dataZ, fileNum, sizeNode, outputTime
    

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
            dataW.close()   
            pass

        fileNum += 1

        ifile=(it-1)/ihfskip+1
        dataC, dataU, dataV, dataW, dataZ = createNetCDFFileHeaders(fileNum, nNodes, nElements, nSides, nLayers, opt)


    if dataC is None:
        return

    outputTime -= time.perf_counter()


    idxTime = ((it // nspool) - 1) % timesPerFile

    # Horizontal geometric variables are only written once per file
    if idxTime == 0:

        xy_cpu = monitor.getNodeVar("outputCoords").get()
        dataC.variables['SCHISM_hgrid_node_x'][:] = xy_cpu[0, :]
        dataC.variables['SCHISM_hgrid_node_y'][:] = xy_cpu[1, :]
        del xy_cpu

        xy_cpu = monitor.getSideVar("outputCoords").get()
        dataC.variables['SCHISM_hgrid_edge_x'][:] = xy_cpu[0, :]
        dataC.variables['SCHISM_hgrid_edge_y'][:] = xy_cpu[1, :]
        del xy_cpu

        xy_cpu = monitor.getElementVar("outputCoords").get()
        dataC.variables['SCHISM_hgrid_face_x'][:] = xy_cpu[0, :]
        dataC.variables['SCHISM_hgrid_face_y'][:] = xy_cpu[1, :]
        del xy_cpu

        view = monitor.getElementVar("FaceIndices")
        mask = (view[:,:] != 0)
        faceIndices_cpu = view.astype(cp.int64).get()
        firstPointer = view[mask].min().get()
        mask = mask.get()
        # To convert pointers to indices, we subtract the first pointer and then divide by the size of the node and add 1 (for 1-indexing)
        faceIndices_cpu[mask] = (faceIndices_cpu[mask] - firstPointer)//sizeNode + 1
        faceIndices_cpu[~mask] = -1
        nSCHISM_hgrid_face_nodes_var = dataC.variables['SCHISM_hgrid_face_nodes']
        nSCHISM_hgrid_face_nodes_var[:] = faceIndices_cpu.T
        del faceIndices_cpu, mask

        view = monitor.getSideVar("n1")
        n1_cpu = view.get()
        firstPointer = view.min().get()
        # To convert pointers to indices, we subtract the first pointer and then divide by the size of the node and add 1 (for 1-indexing)
        n1_cpu[:] = (n1_cpu - firstPointer)//sizeNode + 1

        view = monitor.getSideVar("n2")
        n2_cpu = view.get()
        firstPointer = view.min().get()
        # To convert pointers to indices, we subtract the first pointer and then divide by the size of the node and add 1 (for 1-indexing)
        n2_cpu[:] = (n2_cpu - firstPointer)//sizeNode + 1        

        nSCHISM_hgrid_edge_nodes_var = dataC.variables['SCHISM_hgrid_edge_nodes']
        nSCHISM_hgrid_edge_nodes_var[:, 0] = n1_cpu
        nSCHISM_hgrid_edge_nodes_var[:, 1] = n2_cpu
        del n1_cpu, n2_cpu


        view = monitor.getNodeVar("depth")
        depth_cpu = view.get()
        depth_var = dataC.variables['depth']
        depth_var[:] = depth_cpu
        del depth_cpu

        view = monitor.getNodeVar("kb")
        kb_cpu = view.get()
        kb_var = dataC.variables['bottom_index_node']
        kb_var[:] = kb_cpu
        del kb_cpu

        min_depth_var = dataC.variables['minimum_depth']
        min_depth_var[:] = 0.0001

        


    # For each variable, get the data from GPU (using monitor) and store it in the NetCDF file.
    # Explicitly delete the data from CPU after storing it in the NetCDF file.
    view = monitor.getNodeVar("eta")
    eta_cpu = view.get()
    elevation_var = dataC.variables['elevation']
    elevation_var[idxTime, :] = eta_cpu
    del eta_cpu


    view = monitor.getNodeVar("u_air")
    u_air_cpu = view.get()
    windSpeedX_var = dataC.variables['windSpeedX']
    windSpeedX_var[idxTime, :] = u_air_cpu[0, :]
    windSpeedY_var = dataC.variables['windSpeedY']
    windSpeedY_var[idxTime, :] = u_air_cpu[1, :]
    del u_air_cpu

    view = monitor.getNodeVar("isdry")
    isdry_cpu = view.astype(cp.int32).get()
    dryFlagNode_var = dataC.variables['dryFlagNode']
    dryFlagNode_var[idxTime, :] = isdry_cpu
    del isdry_cpu

    view = monitor.getSideVar("isdry")
    isdry_cpu = view.astype(cp.int32).get()
    dryFlagSide_var = dataC.variables['dryFlagSide']
    dryFlagSide_var[idxTime, :] = isdry_cpu
    del isdry_cpu    

    view = monitor.getElementVar("isdry")
    isdry_cpu = view.astype(cp.int32).get()
    dryFlagElement_var = dataC.variables['dryFlagElement']
    dryFlagElement_var[idxTime, :] = isdry_cpu
    del isdry_cpu    

    view = monitor.getNodeVar("u")
    u_cpu = view[0,:,:].get()
    v_cpu = view[1,:,:].get()

    # Mask depth layers above kb (dry/above bed) as NaN
    kb_cpu = monitor.getNodeVar("kb").get()
    nvrt_u = u_cpu.shape[0]
    k_2d = np.arange(nvrt_u, dtype=kb_cpu.dtype)[:, np.newaxis]
    mask_above_kb = k_2d < kb_cpu[np.newaxis, :]
    u_cpu[mask_above_kb] = np.nan
    v_cpu[mask_above_kb] = np.nan

    horizontalVelX_var = dataU.variables['horizontalVelX']
    horizontalVelX_var[idxTime, :, :] = u_cpu.T

    horizontalVelY_var = dataV.variables['horizontalVelY']
    horizontalVelY_var[idxTime, :, :] = v_cpu.T
    del u_cpu, v_cpu

    view = monitor.getNodeVar("z")
    z_cpu = view.get()
    zCoordinates_var = dataZ.variables['zCoordinates']
    zCoordinates_var[idxTime, :, :] = z_cpu.T
    del z_cpu

    view = monitor.getNodeVar("w")
    w_cpu = view.get()
    verticalVel_var = dataW.variables['verticalVelocity']
    verticalVel_var[idxTime, :, :] = w_cpu.T
    del w_cpu
   

    time_var = dataC.variables['time']
    time_var[idxTime] = it*dt

    time_var = dataU.variables['time']
    time_var[idxTime] = it*dt

    time_var = dataV.variables['time']
    time_var[idxTime] = it*dt

    time_var = dataW.variables['time']
    time_var[idxTime] = it*dt

    time_var = dataZ.variables['time']
    time_var[idxTime] = it*dt

    outputTime += time.perf_counter()


def getOutputTime():
    global outputTime
    return outputTime

def createNetCDFFileHeaders(fileNum, nNodes, nElements, nSides, nLayers, opt):
    """
    Create the NetCDF file "header", meaning the appropriate dimensions, variables, and attributes.
    The data variable is stored in the global variable data.
    """

    global dataC, dataU, dataV, dataW, dataZ
    
    filename = 'out2d_%i.nc' % fileNum
    filenameU = 'horizontalVelX_%i.nc' % fileNum
    filenameV = 'horizontalVelY_%i.nc' % fileNum
    filenameW = 'verticalVelocity_%i.nc' % fileNum
    filenameZ = 'zCoordinates_%i.nc' % fileNum

    isGeographic = getattr(opt, "ics", 1) == 2

    # Create the directory if it doesn't exist
    os.makedirs(path, exist_ok=True)

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

    filepathW = os.path.join(path, filenameW)
    if os.path.exists(filepathW):
        os.remove(filepathW)

    filepathZ = os.path.join(path, filenameZ)
    if os.path.exists(filepathZ):
        os.remove(filepathZ)

    # Create the NetCDF files
    dataC = nc.Dataset(filepath,  "w", format="NETCDF4")
    dataU = nc.Dataset(filepathU, "w", format="NETCDF4")
    dataV = nc.Dataset(filepathV, "w", format="NETCDF4")
    dataW = nc.Dataset(filepathW, "w", format="NETCDF4")
    dataZ = nc.Dataset(filepathZ, "w", format="NETCDF4")
    
    
    # Add dimensions to the file
    for data in [dataC, dataU, dataV, dataW, dataZ]:
        data.createDimension('nSCHISM_hgrid_node',          nNodes)
        data.createDimension('nSCHISM_hgrid_face',          nElements)
        data.createDimension('nSCHISM_hgrid_edge',          nSides)
        data.createDimension('nMaxSCHISM_hgrid_face_nodes', 4)
        data.createDimension('nSCHISM_vgrid_layers',        nLayers)
        data.createDimension('one',                         1)
        data.createDimension('two',                         2)
        data.createDimension('time',                        None)

                # CF grid mapping variable for ics==2
        if getattr(opt, "ics", 1) > 1:
            crs_var = data.createVariable("crs", "i4", ("one",))
            crs_var.long_name = "Coordinate reference system (CRS) definition"
            crs_var.grid_mapping_name = "latitude_longitude"
            crs_var.longitude_of_prime_meridian = 0.0
            crs_var.semi_major_axis = 6378137.0
            crs_var.inverse_flattening = 298.257223563
            crs_var[:] = 0  # value is irrelevant; attributes carry the meaning

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
    node_x_var.mesh = "SCHISM_hgrid"
    node_x_var.axis = "X"
    node_x_var.location = "node"
    if isGeographic:
        # Node coordinates are lon/lat
        node_x_var.units = "degree_E"
        node_x_var.standard_name = "longitude"

    else:
        # Node coordinates are projection meters
        node_x_var.units = "m"
        node_x_var.standard_name = "projection_x_coordinate"


    node_y_var = dataC.createVariable('SCHISM_hgrid_node_y', 'f8', ('nSCHISM_hgrid_node',))
    node_y_var.axis = "Y"
    node_y_var.location = "node"
    node_y_var.mesh = "SCHISM_hgrid"
    if isGeographic:
        node_y_var.units = "degree_N"
        node_y_var.standard_name = "latitude"
    else:
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

    if isGeographic:
        face_x_var.units = "degree_E"
        face_x_var.standard_name = "longitude"
    else:
        face_x_var.units = "m"
        face_x_var.standard_name = "projection_x_coordinate"

    face_y_var = dataC.createVariable('SCHISM_hgrid_face_y', 'f8', ('nSCHISM_hgrid_face',))
    face_y_var.axis = "Y"
    face_y_var.location = "face"
    face_y_var.mesh = "SCHISM_hgrid"

    if isGeographic:
        face_y_var.units = "degree_N"
        face_y_var.standard_name = "latitude"
    else:
        face_y_var.units = "m"
        face_y_var.standard_name = "projection_y_coordinate"
    
    # Edge coordinates
    edge_x_var = dataC.createVariable('SCHISM_hgrid_edge_x', 'f8', ('nSCHISM_hgrid_edge',))
    edge_x_var.axis = "X"
    edge_x_var.location = "edge"
    edge_x_var.mesh = "SCHISM_hgrid"

    if isGeographic:
        edge_x_var.units = "degree_E"
        edge_x_var.standard_name = "longitude"
    else:
        edge_x_var.units = "m"
        edge_x_var.standard_name = "projection_x_coordinate"

    edge_y_var = dataC.createVariable('SCHISM_hgrid_edge_y', 'f8', ('nSCHISM_hgrid_edge',))
    edge_y_var.axis = "Y"
    edge_y_var.location = "edge"
    edge_y_var.mesh = "SCHISM_hgrid"

    if isGeographic:
        edge_y_var.units = "degree_N"
        edge_y_var.standard_name = "latitude"
    else:
        edge_y_var.units = "m"
        edge_y_var.standard_name = "projection_y_coordinate"
    
    # Connectivity arrays
    face_nodes_var = dataC.createVariable('SCHISM_hgrid_face_nodes', 'i4', ('nSCHISM_hgrid_face', 'nMaxSCHISM_hgrid_face_nodes'), fill_value=-1, **paramsCompression)
    face_nodes_var.start_index = 1
    face_nodes_var.cf_role = "face_node_connectivity"
    
    edge_nodes_var = dataC.createVariable('SCHISM_hgrid_edge_nodes', 'i4', ('nSCHISM_hgrid_edge', 'two'), fill_value=-1, **paramsCompression)
    edge_nodes_var.start_index = 1
    edge_nodes_var.cf_role = "edge_node_connectivity"
    
    # Time-dependent variables
    dry_flag_node_var = dataC.createVariable('dryFlagNode', 'f4', ('time', 'nSCHISM_hgrid_node'), **paramsCompression)
    dry_flag_node_var.i23d = cp.int32(1)
    dry_flag_node_var.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    dry_flag_node_var.location = "node"
    dry_flag_node_var.grid_mapping = "crs"
    dry_flag_node_var.mesh = "SCHISM_hgrid"
    
    elevation_var = dataC.createVariable('elevation', 'f4', ('time', 'nSCHISM_hgrid_node'), chunksizes=(1,nNodes), **paramsCompression)
    elevation_var.i23d = cp.int32(1)
    elevation_var.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    elevation_var.location = "node"
    elevation_var.grid_mapping = "crs"
    elevation_var.mesh = "SCHISM_hgrid"
    
    wind_speed_x_var = dataC.createVariable('windSpeedX', 'f4', ('time', 'nSCHISM_hgrid_node'), chunksizes=(1,nNodes), **paramsCompression)
    wind_speed_x_var.i23d = cp.int32(1)
    wind_speed_x_var.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    wind_speed_x_var.location = "node"
    wind_speed_x_var.grid_mapping = "crs"
    wind_speed_x_var.mesh = "SCHISM_hgrid"
    
    wind_speed_y_var = dataC.createVariable('windSpeedY', 'f4', ('time', 'nSCHISM_hgrid_node'), chunksizes=(1,nNodes), **paramsCompression)
    wind_speed_y_var.i23d = cp.int32(1)
    wind_speed_y_var.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    wind_speed_y_var.location = "node"
    wind_speed_y_var.grid_mapping = "crs"
    wind_speed_y_var.mesh = "SCHISM_hgrid"
    
    dry_flag_element_var = dataC.createVariable('dryFlagElement', 'f4', ('time', 'nSCHISM_hgrid_face'), chunksizes=(1,nElements), **paramsCompression)
    dry_flag_element_var.i23d = cp.int32(4)
    dry_flag_element_var.coordinates = "SCHISM_hgrid_face_x SCHISM_hgrid_face_y"
    dry_flag_element_var.location = "face"
    dry_flag_element_var.grid_mapping = "crs"
    dry_flag_element_var.mesh = "SCHISM_hgrid"
    
    dry_flag_side_var = dataC.createVariable('dryFlagSide', 'f4', ('time', 'nSCHISM_hgrid_edge'), chunksizes=(1,nSides), **paramsCompression)
    dry_flag_side_var.i23d = cp.int32(7)
    dry_flag_side_var.coordinates = "SCHISM_hgrid_edge_x SCHISM_hgrid_edge_y"
    dry_flag_side_var.location = "edge"
    dry_flag_side_var.grid_mapping = "crs"
    dry_flag_side_var.mesh = "SCHISM_hgrid"


    horizontalVelX = dataU.createVariable('horizontalVelX', 'f4', ('time', 'nSCHISM_hgrid_node', 'nSCHISM_vgrid_layers'), chunksizes=(1,nNodes,nLayers), **paramsCompression)
    horizontalVelX.i23d = cp.int32(2)
    horizontalVelX.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    horizontalVelX.location = "node"
    horizontalVelX.grid_mapping = "crs"
    horizontalVelX.mesh = "SCHISM_hgrid" 

    horizontalVelY = dataV.createVariable('horizontalVelY', 'f4', ('time', 'nSCHISM_hgrid_node', 'nSCHISM_vgrid_layers'), chunksizes=(1,nNodes,nLayers), **paramsCompression)
    horizontalVelY.i23d = cp.int32(2)
    horizontalVelY.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    horizontalVelY.location = "node"
    horizontalVelY.grid_mapping = "crs"
    horizontalVelY.mesh = "SCHISM_hgrid"    

    verticalVel = dataW.createVariable('verticalVelocity', 'f4', ('time', 'nSCHISM_hgrid_node', 'nSCHISM_vgrid_layers'), chunksizes=(1,nNodes,nLayers), **paramsCompression)
    verticalVel.i23d = cp.int32(2)
    verticalVel.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    verticalVel.location = "node"
    verticalVel.grid_mapping = "crs"
    verticalVel.mesh = "SCHISM_hgrid" 

    # Add zCoordinates variable to dataZ file
    zCoordinates = dataZ.createVariable('zCoordinates', 'f4', ('time', 'nSCHISM_hgrid_node', 'nSCHISM_vgrid_layers'))
    zCoordinates.i23d = cp.int32(2)
    zCoordinates.missing_value = 9.96921e+36
    zCoordinates.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    zCoordinates.location = "node"
    zCoordinates.grid_mapping = "crs"
    zCoordinates.mesh = "SCHISM_hgrid"

    
    return dataC, dataU, dataV, dataW, dataZ

    




