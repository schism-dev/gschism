import os
import cupy as cp

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
data = None
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


def writeOutput(it, core):
    # Write the output files
    global data, fileNum, sizeNode

    ihfskip = core.ihfskip
    nspool  = core.nspool
    dt      = core.dt


    if (it % nspool) != 0:
        return

    if (it-nspool) % ihfskip == 0:
        if data is not None:
            data.close()

        fileNum += 1

        ifile=(it-1)/ihfskip+1
        data = createNetCDFFileHeader('output_%i.nc' % fileNum, nNodes, nElements, nSides, nLayers)

    if data is None:
        return

    idxTime = (it // nspool) - 1

    # Geometric variables are only written once
    if idxTime == 0:
        view = monitor.getNodeVar("x")  # x is a 2D vector
        xy_cpu = view.get()
        SCHISM_hgrid_node_x_var = data.variables['SCHISM_hgrid_node_x']
        SCHISM_hgrid_node_x_var[:] = xy_cpu[0, :]
        
        SCHISM_hgrid_node_y_var = data.variables['SCHISM_hgrid_node_y']
        SCHISM_hgrid_node_y_var[:] = xy_cpu[1, :]
        del xy_cpu

        view = monitor.getSideVar("x")  # x is a 2D vector
        xy_cpu = view.get()
        SCHISM_hgrid_edge_x_var = data.variables['SCHISM_hgrid_edge_x']
        SCHISM_hgrid_edge_x_var[:] = xy_cpu[0, :]
        
        SCHISM_hgrid_edge_y_var = data.variables['SCHISM_hgrid_edge_y']
        SCHISM_hgrid_edge_y_var[:] = xy_cpu[1, :]
        del xy_cpu        

        view = monitor.getElementVar("x")  # x is a 2D vector
        xy_cpu = view.get()
        SCHISM_hgrid_face_x_var = data.variables['SCHISM_hgrid_face_x']
        SCHISM_hgrid_face_x_var[:] = xy_cpu[0, :]
        
        SCHISM_hgrid_face_y_var = data.variables['SCHISM_hgrid_face_y']
        SCHISM_hgrid_face_y_var[:] = xy_cpu[1, :]
        del xy_cpu

        view = monitor.getElementVar("FaceIndices")
        mask = (view[:,:] != 0)
        faceIndices_cpu = view.astype(cp.int64).get()
        firstPointer = view[mask].min().get()
        mask = mask.get()
        # To convert pointers to indices, we subtract the first pointer and then divide by the size of the node and add 1 (for 1-indexing)
        faceIndices_cpu[mask] = (faceIndices_cpu[mask] - firstPointer)//sizeNode + 1
        faceIndices_cpu[~mask] = -1
        nSCHISM_hgrid_face_nodes_var = data.variables['SCHISM_hgrid_face_nodes']
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

        nSCHISM_hgrid_edge_nodes_var = data.variables['SCHISM_hgrid_edge_nodes']
        nSCHISM_hgrid_edge_nodes_var[:, 0] = n1_cpu
        nSCHISM_hgrid_edge_nodes_var[:, 1] = n2_cpu
        del n1_cpu, n2_cpu

        


    # For each variable, get the data from GPU (using monitor) and store it in the NetCDF file.
    # Explicitly delete the data from CPU after storing it in the NetCDF file.
    view = monitor.getNodeVar("eta")
    eta_cpu = view.get()
    elevation_var = data.variables['elevation']
    elevation_var[idxTime, :] = eta_cpu
    del eta_cpu

    view = monitor.getNodeVar("isdry")
    isdry_cpu = view.astype(cp.int32).get()
    dryFlagNode_var = data.variables['dryFlagNode']
    dryFlagNode_var[idxTime, :] = isdry_cpu
    del isdry_cpu

    view = monitor.getSideVar("isdry")
    isdry_cpu = view.astype(cp.int32).get()
    dryFlagSide_var = data.variables['dryFlagSide']
    dryFlagSide_var[idxTime, :] = isdry_cpu
    del isdry_cpu    

    view = monitor.getElementVar("isdry")
    isdry_cpu = view.astype(cp.int32).get()
    dryFlagElement_var = data.variables['dryFlagElement']
    dryFlagElement_var[idxTime, :] = isdry_cpu
    del isdry_cpu    

   

    time_var = data.variables['time']
    time_var[idxTime] = it*dt









def closeOutputFiles():
    # Close the output files
    pass


def createNetCDFFileHeader(filename, nNodes, nElements, nSides, nLayers):
    """
    Create the NetCDF file "header", meaning the appropriate dimensions, variables, and attributes.
    The data variable is stored in the global variable data.
    """

    global data
    
    # Remove the file if it exists before creating the NetCDF file
    filepath = os.path.join(path, filename)
    if os.path.exists(filepath):
        os.remove(filepath)
    
    # Create the NetCDF file
    data = nc.Dataset(filepath, "w", format="NETCDF4")
    
    
    # Add dimensions to the file
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
    time_var.i23d = 0
    time_var.base_date = " 2017 12  1       0.00       0.00"
    time_var.units = "seconds since 2017-12-01T00:00:00+00"
    time_var.standard_name = "time"
    time_var.axis = "T"
    
    # Minimum depth
    min_depth_var = data.createVariable('minimum_depth', 'f8', ('one',))
    min_depth_var.units = "m"
    
    # SCHISM grid topology
    schism_grid_var = data.createVariable('SCHISM_hgrid', 'S1', ('one',))
    schism_grid_var.long_name = "Topology data of 2d unstructured mesh"
    schism_grid_var.topology_dimension = 2
    schism_grid_var.cf_role = "mesh_topology"
    schism_grid_var.node_coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    schism_grid_var.edge_coordinates = "SCHISM_hgrid_edge_x SCHISM_hgrid_edge_y"
    schism_grid_var.face_coordinates = "SCHISM_hgrid_face_x SCHISM_hgrid_face_y"
    schism_grid_var.edge_node_connectivity = "SCHISM_hgrid_edge_nodes"
    schism_grid_var.face_node_connectivity = "SCHISM_hgrid_face_nodes"
    
    # Node coordinates
    node_x_var = data.createVariable('SCHISM_hgrid_node_x', 'f8', ('nSCHISM_hgrid_node',))
    node_x_var.axis = "X"
    node_x_var.location = "node"
    node_x_var.mesh = "SCHISM_hgrid"
    node_x_var.units = "m"
    node_x_var.standard_name = "projection_x_coordinate"
    
    node_y_var = data.createVariable('SCHISM_hgrid_node_y', 'f8', ('nSCHISM_hgrid_node',))
    node_y_var.axis = "Y"
    node_y_var.location = "node"
    node_y_var.mesh = "SCHISM_hgrid"
    node_y_var.units = "m"
    node_y_var.standard_name = "projection_y_coordinate"
    
    # Depth and bottom index
    depth_var = data.createVariable('depth', 'f4', ('nSCHISM_hgrid_node',))
    depth_var.units = "m"
    depth_var.axis = "Z"
    depth_var.positive = "down"
    depth_var.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    depth_var.location = "node"
    depth_var.grid_mapping = "crs"
    depth_var.mesh = "SCHISM_hgrid"
    
    bottom_idx_var = data.createVariable('bottom_index_node', 'i4', ('nSCHISM_hgrid_node',))
    bottom_idx_var.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    bottom_idx_var.location = "node"
    bottom_idx_var.grid_mapping = "crs"
    bottom_idx_var.mesh = "SCHISM_hgrid"
    
    # Face coordinates
    face_x_var = data.createVariable('SCHISM_hgrid_face_x', 'f8', ('nSCHISM_hgrid_face',))
    face_x_var.axis = "X"
    face_x_var.location = "face"
    face_x_var.mesh = "SCHISM_hgrid"
    face_x_var.units = "m"
    face_x_var.standard_name = "projection_x_coordinate"
    
    face_y_var = data.createVariable('SCHISM_hgrid_face_y', 'f8', ('nSCHISM_hgrid_face',))
    face_y_var.axis = "Y"
    face_y_var.location = "face"
    face_y_var.mesh = "SCHISM_hgrid"
    face_y_var.units = "m"
    face_y_var.standard_name = "projection_y_coordinate"
    
    # Edge coordinates
    edge_x_var = data.createVariable('SCHISM_hgrid_edge_x', 'f8', ('nSCHISM_hgrid_edge',))
    edge_x_var.axis = "X"
    edge_x_var.location = "edge"
    edge_x_var.mesh = "SCHISM_hgrid"
    edge_x_var.units = "m"
    edge_x_var.standard_name = "projection_x_coordinate"
    
    edge_y_var = data.createVariable('SCHISM_hgrid_edge_y', 'f8', ('nSCHISM_hgrid_edge',))
    edge_y_var.axis = "Y"
    edge_y_var.location = "edge"
    edge_y_var.mesh = "SCHISM_hgrid"
    edge_y_var.units = "m"
    edge_y_var.standard_name = "projection_y_coordinate"
    
    # Connectivity arrays
    face_nodes_var = data.createVariable('SCHISM_hgrid_face_nodes', 'i4', ('nSCHISM_hgrid_face', 'nMaxSCHISM_hgrid_face_nodes'), fill_value=-1)
    face_nodes_var.start_index = 1
    face_nodes_var.cf_role = "face_node_connectivity"
    
    edge_nodes_var = data.createVariable('SCHISM_hgrid_edge_nodes', 'i4', ('nSCHISM_hgrid_edge', 'two'), fill_value=-1)
    edge_nodes_var.start_index = 1
    edge_nodes_var.cf_role = "edge_node_connectivity"
    
    # Time-dependent variables
    dry_flag_node_var = data.createVariable('dryFlagNode', 'f4', ('time', 'nSCHISM_hgrid_node'))
    dry_flag_node_var.i23d = 1
    dry_flag_node_var.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    dry_flag_node_var.location = "node"
    dry_flag_node_var.grid_mapping = "crs"
    dry_flag_node_var.mesh = "SCHISM_hgrid"
    
    elevation_var = data.createVariable('elevation', 'f4', ('time', 'nSCHISM_hgrid_node'))
    elevation_var.i23d = 1
    elevation_var.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    elevation_var.location = "node"
    elevation_var.grid_mapping = "crs"
    elevation_var.mesh = "SCHISM_hgrid"
    
    wind_speed_x_var = data.createVariable('windSpeedX', 'f4', ('time', 'nSCHISM_hgrid_node'))
    wind_speed_x_var.i23d = 1
    wind_speed_x_var.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    wind_speed_x_var.location = "node"
    wind_speed_x_var.grid_mapping = "crs"
    wind_speed_x_var.mesh = "SCHISM_hgrid"
    
    wind_speed_y_var = data.createVariable('windSpeedY', 'f4', ('time', 'nSCHISM_hgrid_node'))
    wind_speed_y_var.i23d = 1
    wind_speed_y_var.coordinates = "SCHISM_hgrid_node_x SCHISM_hgrid_node_y"
    wind_speed_y_var.location = "node"
    wind_speed_y_var.grid_mapping = "crs"
    wind_speed_y_var.mesh = "SCHISM_hgrid"
    
    dry_flag_element_var = data.createVariable('dryFlagElement', 'f4', ('time', 'nSCHISM_hgrid_face'))
    dry_flag_element_var.i23d = 4
    dry_flag_element_var.coordinates = "SCHISM_hgrid_face_x SCHISM_hgrid_face_y"
    dry_flag_element_var.location = "face"
    dry_flag_element_var.grid_mapping = "crs"
    dry_flag_element_var.mesh = "SCHISM_hgrid"
    
    dry_flag_side_var = data.createVariable('dryFlagSide', 'f4', ('time', 'nSCHISM_hgrid_edge'))
    dry_flag_side_var.i23d = 7
    dry_flag_side_var.coordinates = "SCHISM_hgrid_edge_x SCHISM_hgrid_edge_y"
    dry_flag_side_var.location = "edge"
    dry_flag_side_var.grid_mapping = "crs"
    dry_flag_side_var.mesh = "SCHISM_hgrid"
    
    return data

    




