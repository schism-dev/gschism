from os import path
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
from typing import Tuple
from iop.logger import *
from iop.readNamelist import readNamelist
from scipy.spatial import Delaunay

"""
This module handles meteorological data interpolation for the SCHISM model.

The module is organized around two main classes:
- SFlux: Main interface class that manages meteorological data interpolation
- DatasetInfo: Internal class that handles individual meteorological datasets

Main usage:
    1. Create an SFlux object with a namelist file:
       sflux = SFlux("path/to/sflux_inputs.nml", nodes)

       The object will construct the two datasets and compute the interpolation weights.
      
    2. Get interpolated meteorological data at any time:
       sflux.getWind(time_now, u_air_node, v_air_node, p_air_node, t_air_node, q_air_node)

The SFlux class can handle two datasets of different resolutions (typically global/regional) 
and combines them using weights specified in the namelist. It automatically handles:

    - Temporal interpolation to the current timestep
    - Spatial interpolation to the computational mesh
    - Dataset combination using specified weights

By default, dataset 2 has a relative weight of 99.0, which means that if it is used for the high-resolution
data, in the areas where it is available, it will use (mostly) the high-resolution data.
"""




t_freeze = 273.15


class DatasetInfo:
    """Represents a single sflux dataset and its interpolation information"""
    def __init__(self, in_dir: str,name: str, nodes: np.ndarray):
        
        self.name = name
        self.in_dir = in_dir


        self.out_x = nodes[0, :]
        self.out_y = nodes[1, :]

        self.max_file_times = 1000 # max. # of time records with each file
        self.max_files = 9999   # max. # of files. XXX Are these really needed in Python?


        # Default properties
        self.exist = False
        self.max_window_hours = 120.0
        self.fail_if_missing = False
        self.relative_weight = 1.0

        # Grid information
        self.nx = None
        self.ny = None
        self.lat = None
        self.lon = None
        
        # Time information
        self.times = None
        self.file_num_for_time = None  # mapping of times to file numbers
        self.time_idx_for_time = None  # mapping of times to local indices
        
        # Interpolation data
        self.weights = None
        self.interp_weights = None
        self.interp_indices = None

        # Obtain information from the first file of the dataset
        # -----------------------------------------------------

        # determine if the first file for this dataset exists
        file_name = self.getFileName(1)
        self.exist = path.exists(path.join(self.in_dir, file_name))
            
        if not self.exist:
            return

        # run should fail if dataset doesn't exist and fail_if_missing is set
        if (not self.exist) and (self.fail_if_missing):
            raise ValueError(f'missing dataset: {file_name}')

        
        infoMessage(f'getting additional info for: {self.name}')

        try:
            # Obtain information about the times of all the files in the dataset for which
            # there is data. 
            self.computeTimeInfo() 

        except Exception as e:
            if self.fail_if_missing:
                raise ValueError(f'Error getting times and files for {self.name}: {e}')            
            

        # read in lon and lat from the first file
        ds = Dataset(path.join(self.in_dir, file_name), 'r')
        try:
            self.lat = ds.variables['lat'][:,:].astype(np.float64)
            self.lon = ds.variables['lon'][:,:].astype(np.float64)
            self.lat = np.deg2rad(self.lat)
            self.lon = np.deg2rad(self.lon)

            # Save the shape and total number of nodes for the sflux grid
            self.nx, self.ny = self.lon.shape
            self.num_nodes = self.nx*self.ny


        except Exception as e:
            raise ValueError(f'Error reading lon and lat for {self.name}: {e}')
        
        ds.close()
            

    def getFileName(self, num: int) -> str:
        """
        Get the file name for a given file number.
        """

        prefix = 'sflux/'
        suffix = '.nc'
        
        return f'{prefix}{self.name}.{num:04d}{suffix}'


    def filter_triangles_by_grid(self, tri, nx, ny, i_orig, j_orig):
        """ 
        Remove triangles that are not inside the grid boundary.
        This is necessary because the Delaunay triangulation may create triangles that are
        outside the grid boundary. 
        """

        # Get i,j indices for all triangles at once
        i_indices = i_orig[tri.simplices]  # Shape: (n_triangles, 3)
        j_indices = j_orig[tri.simplices]  # Shape: (n_triangles, 3)
        
        # Calculate max differences for each triangle.
        # The idea is to find triangles where the vertices are immediate neighbors in the grid.
        i_diff = np.max(i_indices, axis=1) - np.min(i_indices, axis=1)
        j_diff = np.max(j_indices, axis=1) - np.min(j_indices, axis=1)
        
        # Create masks for boundary nodes
        # The idea is to find triangles where the vertices are *all* on the boundary of the grid.
        # These triangles are most likely to be outside of the grid boundary (created by
        # the Delaunay triangulation). However, we check that we are not excluding any
        # valid triangles.
        boundary_i = (i_indices == 0) | (i_indices == nx - 1)
        boundary_j = (j_indices == 0) | (j_indices == ny - 1)
        boundary_nodes = boundary_i | boundary_j
        
        # Triangle is valid if not all vertices are on the boundary of the grid.
        # Unless, i,j of the nodes have a difference of 1, in which case we assume
        # they are corner triangles.
        valid_mask = ((i_diff <= 1) & (j_diff <= 1) | ~(np.all(boundary_nodes, axis=1)))
        
        return tri.simplices[valid_mask], valid_mask
    

    def computeTriangleSignedArea(self, x1, y1, x2, y2, x3, y3):
        return 0.5 * ((x1 - x3)*(y2 - y3) - (y1 - y3)*(x2 - x3))
    
    
    def computeGridInterpWeights(self):
        """
        Calculate interpolation weights for unstructured grid interpolation using barycentric coordinates.

        This function is divided into two parts:
        1. Find the triangles that contain each output point
        2. Compute the interpolation weights for each output point
        
        1. Find the triangles that contain each output point:
            Performs a Delaunay triangulation of the input grid and filters out triangles that are not inside the grid boundary.
            This is necessary because the Delaunay triangulation may create triangles that are
            outside the grid boundary.
            Reorders the output points to increase locality and accelerate the search for the containing triangle.
            Finds the containing triangle for each output point.


        2. Compute the interpolation weights for each output point:
            Computes the barycentric coordinates as interpolation weights. These weights (w1,w2,w3) represent
            ratios of areas:
                w1 = Area(triangle_234)/Area(triangle_123)
                w2 = Area(triangle_134)/Area(triangle_123)
                w3 = Area(triangle_124)/Area(triangle_123)
            where 1,2,3 are the vertices of the input triangle and 4 is the output point.

            This makes sense geometrically: if point 4 is on edge 23, w1 = 0, since the interpolated 
            value should only depend on vertices 2 and 3.
            
            The sum of weights equals 1 when the point is inside the triangle. The weights
            can be used to interpolate any quantity from the input triangle vertices to the
            output point: value = w1*v1 + w2*v2 + w3*v3
        
        Args:
            nodes: (num_nodes_out, 2) array of (x,y) coordinates with the computational mesh node locations
        
        Returns:
            weights: Interpolation weights (3, num_nodes_out)
            indices: Indices of the containing triangles (3, num_nodes_out)
        """




        # 1. Find the triangles that contain each output point:

        # Get the coordinates of the nodes of the input grid and the indices for each node
        i, j = np.indices(self.lat.shape)
        x_grid = np.ma.filled(self.lon.flatten().astype(np.float64), fill_value=None)
        y_grid = np.ma.filled(self.lat.flatten().astype(np.float64), fill_value=None)
        
        # Create the Delaunay triangulation
        tri = Delaunay(np.column_stack((x_grid, y_grid)))


        # Filter out triangles that are outside the grid boundary
        valid_simplices, valid_simplices_mask = self.filter_triangles_by_grid(tri, self.nx, self.ny, i.flatten(), j.flatten())

        # Create virtual grid for sorting. The idea is to have an MxM grid where each cell
        # contains a number of points. This way, we can sort the points by their grid position
        # and then find which triangle contains each point. Since the points in each cell are
        # close to each other, finding the triangle is faster.
        M = 100  # Size of virtual grid (adjust as needed, different values have different performance)

        # Calculate virtual grid indices for each point
        x_min = self.out_x.min()
        x_max = self.out_x.max()
        y_min = self.out_y.min()
        y_max = self.out_y.max()
        tmp_x = np.floor((self.out_x - x_min) / (x_max - x_min) * M).astype(int)
        tmp_y = np.floor((self.out_y - y_min) / (y_max - y_min) * M).astype(int)

        # Create a single index for sorting (row-major order)
        grid_indices = tmp_y * M + tmp_x

        # Sort output points based on their grid position
        order = np.argsort(grid_indices)
        self.out_x = self.out_x[order]
        self.out_y = self.out_y[order]

        # Find which triangle contains each point
        # Note: find_simplex returns -1 for points outside the triangulation
        point_tri_indices = tri.find_simplex(np.column_stack((self.out_x, self.out_y)))

        # Filter out points that aren't in any triangle (-1) or the triangle is not valid (extra triangles
        # created by the Delaunay triangulation)
        valid_mask = (point_tri_indices != -1) & valid_simplices_mask[point_tri_indices]
        point_tri_indices[~valid_mask] = -1

        # At this point, point_tri_indices contains the indices of the triangles that contain each point.


        # 2. Compute the interpolation weights for each output point:
        
        tri_points = tri.simplices[point_tri_indices[valid_mask]]
        p1x = x_grid[tri_points[:,0]]
        p1y = y_grid[tri_points[:,0]]
        p2x = x_grid[tri_points[:,1]]
        p2y = y_grid[tri_points[:,1]]
        p3x = x_grid[tri_points[:,2]]
        p3y = y_grid[tri_points[:,2]]
        p4x = self.out_x[valid_mask]
        p4y = self.out_y[valid_mask]

        area = np.abs(self.computeTriangleSignedArea(p1x, p1y, p2x, p2y, p3x, p3y))

        # Signed areas of three subtriangles for each point
        a234 = self.computeTriangleSignedArea(p2x, p2y, p3x, p3y, p4x, p4y)
        a134 = self.computeTriangleSignedArea(p1x, p1y, p4x, p4y, p3x, p3y)
        a124 = self.computeTriangleSignedArea(p1x, p1y, p2x, p2y, p4x, p4y)

        if np.any(a234/area < 0) or np.any(a134/area < 0) or np.any(a124/area < 0):
            raise ValueError("Negative area found in subtriangles")

        # Compute and store the weights
        self.weights = np.zeros((3, valid_mask.size))
        self.weights[0, valid_mask] = a234/area
        self.weights[1, valid_mask] = a134/area
        self.weights[2, valid_mask] = a124/area
        
        
        # Recover the original order of the output points
        inverse_order = np.argsort(order)  # Create inverse permutation
        self.interp_weights = self.weights[:, inverse_order]
        self.interp_indices = tri.simplices[point_tri_indices[inverse_order]].T
        self.out_x = self.out_x[inverse_order]
        self.out_y = self.out_y[inverse_order]


    def interpGrid(self, data_in: np.ndarray):

        # XXX CHECK -1 indices!!!
        return (data_in[self.interp_indices]*self.interp_weights).sum(axis=0)


    def julianDay(self, year: int, month: int, day: int) -> int:
        """
        Convert Gregorian date to Julian day.
        Example: julian_day(1970,1,1) = 2440588
        """
        return datetime(year, month, day).toordinal() + 1721423
    
    
    def getFileTimes(self, file_name, file_idx):

        file_name = self.getFileName(file_idx + 1)

        ds = Dataset(path.join(self.in_dir, file_name),'r')

        times = ds.variables['time'][:].astype(np.float64)

        if (times.shape[0] > self.max_file_times):
            raise ValueError('sflux: num_file_times larger than max_file_times!')


        # Read the base_date. This is the time and date which corresponds with time zero
        # Note that only year,month,day of base_date are used (not 'hour').
        base_date = ds.variables['time'].getncattr('base_date')

        if (len(base_date) < 3):
            raise ValueError('sflux: base_date has less than 3 fields!')
        

        # convert to integer Julian date
        year, month, day = base_date[:3]

        ds.close()

        return self.julianDay(year, month, day) + times


    def computeTimeInfo(self):
        # Compute the times of all the files in the dataset for which there is data.
        # It also creates a mapping of global times to file numbers and local time indices.
        

        # Count how many data files there are
        num_files = 0
        while True:
            file_name = self.getFileName(num_files + 1)
            if not path.exists(path.join(self.in_dir, file_name)):
                break
            num_files += 1


        if (num_files > self.max_files):
            raise ValueError('num_files exceeds max_files!')


        # loop over files to get all the times in all the files
        files_times = []
        for file_idx in range(num_files):

            # get the times in this file
            file_julian_times = self.getFileTimes(file_name, file_idx)

            files_times.append(file_julian_times)


        # Creates a mapping of global times to file numbers and local time indices.
        # Upon return, the following variables are set:
        # - global_times: sorted list of unique global times
        # - file_num_for_time: file number for each global time
        # - time_num_for_time: local time index for each global time in its file
        
        # Step 1: Collect all times into a single array
        self.all_times = np.concatenate([times for times in files_times])
        
        # Step 2: Sort and deduplicate the times
        self.all_times = np.unique(self.all_times)

        # Step 3: Initialize file_num_for_time with zeros
        self.file_num_for_time = np.zeros_like(self.all_times, dtype=int)
        self.time_idx_for_time = np.zeros_like(self.all_times, dtype=int)

        # Step 4: Iterate over each file and update file_num_for_time
        for file_idx, times in enumerate(files_times):
            # Find the indices of the file's times in the global time array
            indices = np.searchsorted(self.all_times, times)
            
            # Update file_num_for_time with the current file number
            self.file_num_for_time[indices] = file_idx
            self.time_idx_for_time[indices] = np.arange(indices.size)
            

        log(f'num_files = {num_files}')
        log(f'num_times = {self.all_times.shape[0]}')
        log(f'first time = {self.all_times[ 0]}')
        log(f'last  time = {self.all_times[-1]}')


    def getTimeBracket(self, time: float) -> tuple[bool, int]:
        """
        Find the index where time falls between two consecutive times in input_times.

        Remember that self.all_times is a sorted array of unique times.
        
        Args:
            time: The time to find a bracket for
        
        Returns:
            tuple[bool, int]: (got_bracket, time_num_1) where:
                - got_bracket: True if time falls between two consecutive times
                - time_idx_1: Index of the lower bracketing time
        """
        time_idx_1 = np.searchsorted(self.all_times, time, side='right') - 1
        
        # Check if we found a valid bracket
        if 0 <= time_idx_1 < self.all_times.size - 1:
            got_bracket = self.all_times[time_idx_1] <= time <= self.all_times[time_idx_1 + 1]
        else:
            got_bracket = False
            
        return got_bracket, time_idx_1


    def interpInTime(self, time: np.float64, data_name: str, time_idx_1: int):
        """
        Interpolate data in time between two consecutive timesteps.

        NOTE: In the context of this function, _1 and _2 are the indices of the two consecutive timesteps
        that bracket the target time. It has nothing to do with dataset_1 and dataset_2.

        Args:
            time (float): Target time for interpolation
            data_name (str): Name of the variable to interpolate
            time_idx_1 (int): Index of the lower bracketing time

        Returns:
            np.ndarray: Interpolated data at the target time
            
        Notes:
            - Performs linear interpolation between two timesteps
            - Reads data from appropriate files based on time_idx_1 and time_idx_1 + 1
        """

        file_num_1 = self.file_num_for_time[time_idx_1]
        file_num_2 = self.file_num_for_time[time_idx_1 + 1]

        local_time_idx_1 = self.time_idx_for_time[time_idx_1]       # local time index in the file
        local_time_idx_2 = self.time_idx_for_time[time_idx_1 + 1]   # local time index in the file 

        data_1 = self.readData(file_num_1, data_name, local_time_idx_1).ravel()
        data_2 = self.readData(file_num_2, data_name, local_time_idx_2).ravel()

        weight = (time - self.all_times[time_idx_1]) / (self.all_times[time_idx_1 + 1] - self.all_times[time_idx_1])
            
        # Perform linear interpolation
        return (1 - weight)*data_1 + weight*data_2


    def readData(self, file_num: int, data_name: str, time_idx: int):
        """
        Read data from a netCDF file.
        NOTE: This is an internal function, used by the getData method via the interpInTime method.
        
        Args:
            file_num: Number of the netCDF file (starting with 1)
            data_name: Name of the variable to read
            time_num: Time index to read
        
        Returns:
            data: numpy array of shape (nx, ny)
        """
        try:
            # Open netCDF file in read-only mode
            file_name = self.getFileName(file_num)
            with Dataset(path.join(self.in_dir, file_name), 'r') as nc:
                # Read the data for the specified time
                # netCDF is 0-based, so we subtract 1 from time_num
                data_tmp = nc.variables[data_name][time_idx, :, :]

                # Check if the data has the correct shape
                if data_tmp.shape != (self.nx, self.ny):
                    raise ValueError(f"Data shape mismatch: expected ({self.nx}, {self.ny}), got {data_tmp.shape}")
                
                # Convert to double precision and return
                return data_tmp.astype(np.float64)
                
        except Exception as e:
            print(f"Error reading variable {data_name} from netCDF file: {file_name}\n{e}")
            return None 


    def getData(self, time_now,  data_name):

        hours_per_day = 24.0

        # find bracketing times for the current time
        got_bracket, time_idx_1 = self.getTimeBracket(time_now)
         
        # check to see whether we've exceeded the max_window
        window_hours = hours_per_day * (self.all_times[time_idx_1 + 1] - self.all_times[time_idx_1])
        got_suitable_bracket = got_bracket and (window_hours < self.max_window_hours)
             
        # fail if we don't have a suitable bracket and if fail_if_missing is set
        if not got_suitable_bracket:
            if self.fail_if_missing:
                warning(f'no appropriate time exists for: {self.name}, time_now = {time_now}, first time available = {self.all_times[0]}, last time available = {self.all_times[-1]}, got_bracket = {got_bracket}, got_suitable_bracket = {got_suitable_bracket}')
              
            if (got_bracket):
                infoMessage(f', window_hours = {window_hours}, max_window_hours = {self.max_window_hours}, time_1 = {self.all_times[time_idx_1]}, time_2 = {self.all_times[time_idx_1+1]}')
            
            error(f'No suitable bracket found for {self.name} at time {time_now}')


        # if we do have a suitable bracket, then read in the bracketing data
        # and interpolate the input data in time to the current time (but still on original grid)
        data_in_at_time = self.interpInTime(time_now, data_name, time_idx_1)
            
        # now the data needs to be interpolated in space from the original data
        # grid to the model grid
        data_out = self.interpGrid(data_in_at_time)


        return got_suitable_bracket, data_out
    

    
class SFlux:
    """Class that represents the meteorological data and allows for its interpolation 
       to the computational mesh.
    
    This data can be stored in two datasets, typically of different resolutions, and combined
    using weights.
    The information about the datasets used and the relative weights is stored in a namelist file.
    This is the namelist file name given to the constructor.
    
    This class also takes into account that:
    - The meteo data need to be interpolated in time to the current time step.
    - The meteo data need to be interpolated in space to the computational mesh.

    The interpolation weights, both time and space, are stored in the DatasetInfo objects.
    """

    def __init__(self, in_dir: str, namelist_file: str, nodes: np.ndarray):

        # Read the namelist file

        DEFAULTS = {
          'air_1_fail_if_missing': True,
          'air_2_fail_if_missing': False,
          'air_1_relative_weight': 1.0,
          'air_2_relative_weight': 99.0,
          'air_1_max_window_hours': 120.0,    
          'air_2_max_window_hours': 120.0,
          'air_1_file': 'sflux_air_1',
          'air_2_file': 'sflux_air_2',
          'uwind_name': 'uwind',
          'vwind_name': 'vwind',
          'prmsl_name': 'prmsl',
          'stmp_name': 'stmp',
          'spfh_name': 'spfh',
        }

        self.nodes = nodes

        try:
            # Read the namelist
            nml = readNamelist(path.join(in_dir, 'sflux', namelist_file))

            # Get the 'sflux_inputs' section
            self.inputs = nml.get('sflux_inputs', {})

            # Merge with defaults (only if the key is not already present in the namelist)
            for key, default in DEFAULTS.items():
                self.inputs.setdefault(key, default)


      
        except Exception as e:
            print(f"Error reading namelist: {e}")
            return None
        

        # output the information read from the namelist, including the defaults
        log('get_wind: sflux_inputs')
        for key, value in self.inputs.items():
            log(f'{key} = {value}')
        log('')
        

        # Initialize the two datasets
        self.dataset1 = DatasetInfo(in_dir, self.inputs['air_1_file'], nodes)
        self.dataset1.max_window_hours = self.inputs['air_1_max_window_hours']
        self.dataset1.fail_if_missing  = self.inputs['air_1_fail_if_missing']
        self.dataset1.relative_weight  = self.inputs['air_1_relative_weight']
        
        self.dataset2 = DatasetInfo(in_dir, self.inputs['air_2_file'], nodes)
        self.dataset2.max_window_hours = self.inputs['air_2_max_window_hours']
        self.dataset2.fail_if_missing  = self.inputs['air_2_fail_if_missing']
        self.dataset2.relative_weight  = self.inputs['air_2_relative_weight']


        # Initialize the spatial interpolation weights
        for dataset in [self.dataset1, self.dataset2]:
            if dataset.exist:
                dataset.computeGridInterpWeights()


    def getCombinedData(self, time_now, data_name, data_out):

        got_data_1 = False
        got_data_2 = False

        if (self.dataset1.exist):
            got_data_1, data_1 = self.dataset1.getData(time_now, data_name)


        if (self.dataset2.exist):
            got_data_2, data_2 = self.dataset2.getData(time_now, data_name)

        
        # depending on whether different datasets exist at this time, combine
        # data in different ways

        # if data_1 is missing, then fail
        if (not got_data_1):
            raise ValueError(f'missing data_1: {data_name}')

        # if data_1 exists, but not data_2, then simply use data_1 (for all)
        if (not got_data_2):
            data_out[:] = data_1
            return

        # if data_1 and data_2 exist, then combine
        if (got_data_1 and got_data_2):
            try:
                # Check if any weight for each node is bad (outside [0,1] range)
                # dataset_2.weights has shape (3, num_nodes)
                bad_node_2 = np.any((self.dataset2.weights < 0.0) | (self.dataset2.weights > 1.0), axis=0)
                
                # Set weights based on bad node mask
                local_weight_2 = np.where(bad_node_2, 0.0, self.dataset2.relative_weight)
                local_weight_1 = np.where(bad_node_2, 1.0, self.dataset1.relative_weight)

                data_out[:] = local_weight_1*data_1 + local_weight_2*data_2

            except Exception as e:
                raise ValueError(f'Error combining sflux data from two datasets: {e}')


    def rotate_winds(self, u, v):

 
        pass

        assert(False and 'IMPLEMENT THIS!!!!')
        # # XXX TODO: Do this using a rotation matrix (or complex numbers) instead of trigonometric functions
        # # calculate speed and direction (geographic)
        # dir = np.arctan2(-u,-v)
        # speed = np.sqrt(u*u + v*v)

        # # add rotation angle
        # dir += self.wind_rotate_angle

        # # calculate new u and v components
        # u = -speed * np.sin(dir)
        # v = -speed * np.cos(dir)


    def getWind(self, time_now: float, u_air_node: np.ndarray, v_air_node: np.ndarray, p_air_node: np.ndarray, t_air_node: np.ndarray, q_air_node: np.ndarray) -> tuple:

        self.getCombinedData(time_now, self.inputs['uwind_name'], u_air_node)
        self.getCombinedData(time_now, self.inputs['vwind_name'], v_air_node)
        self.getCombinedData(time_now, self.inputs['prmsl_name'], p_air_node)
        self.getCombinedData(time_now, self.inputs['stmp_name' ], t_air_node)
        self.getCombinedData(time_now, self.inputs['spfh_name' ], q_air_node)

           
        # convert air temperatures to Celsius
        t_air_node = t_air_node - t_freeze

        # rotate the winds from geographic to the grid's map projection
        self.rotate_winds(u_air_node, v_air_node)


    def getWinds(self, time: float):
        assert(False and 'IMPLEMENT THIS!!!!, should be similar to getWind()')