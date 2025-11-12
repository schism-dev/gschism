import os
import glob
import importlib.util
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import cupy as cp
import debug.monitor as monitor
import time
from python.pcgsolver import PCGSolver
from iop.writeOutput import writeOutput, initializeOutput

# pool = cp.cuda.MemoryPool()
# cp.cuda.set_allocator(pool.malloc)

# For profiling
from line_profiler import LineProfiler


VISTA = True
WINDOWS = not VISTA


# Enable verbose output
os.environ['DISTUTILS_DEBUG'] = '1'
os.environ['VERBOSE'] = '1'


# Print diagnostic information
print("Python version:", sys.version)
print("Include paths:")
for path in [
    r'C:\Users\jmgonzalez\AppData\Local\miniconda3\envs\gschism\Lib\site-packages\pybind11\include',
    r'C:\Users\jmgonzalez\apps\vcpkg\installed\x64-windows\include',
    'src/include'
]:
    print(f"  {path} exists: {os.path.exists(path)}")

import pybind11
import numpy as np
import timeit

import matplotlib
# matplotlib.use('TkAgg') 
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection


nws = 0

# sys.path.append(r'/home/jo.gonzalez/src/gschism/src')  
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'src'))

from iop.cudamodule import CudaSourceManager
from iop.processConstantsHeader import processConstants
from python.schism_init import readNamelist, schism_init
from python.schism_step import schism_step, solveWaveEquation
from python.schism_consts import *
import python.cuda_functions as cufunc




# # in_dir = '/media/jmg/Elements/src/tmp/schism_verification_tests/Test_Btrack_Gausshill'
# # in_dir = '/media/jmg/Elements/src/tmp/schism_verification_tests/Small_test/Small_grid'
# in_dir = r'C:\Users\jmgonzalez\src\schism_test_suit\Small_grid'
# in_dir = r'C:\Users\jmgonzalez\src\schism_test_suit\QuarterAnnulus\Pure_tri'
# # in_dir = r'C:\Users\jmgonzalez\src\schism_test_suit\RUN02a'


if WINDOWS:
    in_dir = r'C:\Users\jmgonzalez\src\schism_test_suit\RUN03b'
    out_dir = r'C:\Users\jmgonzalez\src\schism_test_suit\RUN03b\gpuoutput'
    # out_dir = r'c:\Users\jmgonzalez\src\temp'
elif VISTA:
    in_dir = r'/scratch/10922/jmglez/tests/RUN01b'
    out_dir = r'/scratch/10922/jmglez/tests/RUN01b/gpuoutput'
    in_dir = r'/home1/10922/jmglez/tests/RUN02a'
    out_dir = r'/home1/10922/jmglez/tests/RUN02a/gpuoutput'
else:
    raise Exception("Unsupported platform")




# import src.io.readSflux as readSflux
from iop.readSflux import SFlux


def compile_and_import(module_name, cpp_file):
    """Compile and import a C++ module dynamically."""

    cpp_file = os.path.abspath(cpp_file)
    
    # Ensure the C++ file exists
    if not os.path.exists(cpp_file):
        raise FileNotFoundError(f"C++ source file '{cpp_file}' does not exist.")

    # Build directory for intermediate files
    build_dir = os.path.abspath("build")
    os.makedirs(build_dir, exist_ok=True)

    # Output file name with platform-specific suffix
    output_files = glob.glob(f"{module_name}*.so")
    output_file = output_files[0] if output_files else f"{module_name}.so"

    # Check caching: Rebuild only if necessary
    if os.path.exists(output_file) and os.path.getmtime(output_file) >= os.path.getmtime(cpp_file):
        print(f"Using cached build: {output_file}")
    else:
        print("Compiling the C++ code...")

        # Get vcpkg root directory
        vcpkg_root = os.path.expanduser("~/apps/vcpkg")


        if VISTA:
            METIS_PREFIX = os.environ.get("TACC_METIS_DIR")
            METIS_INC    = os.path.join(METIS_PREFIX, "include")
            METIS_LIB    = os.path.join(METIS_PREFIX, "lib")
            METIS_LIB64  = os.path.join(METIS_PREFIX, "lib64")

            ext_modules = [
                Extension(
                    module_name,
                    sources=[cpp_file],
                    include_dirs=[
                        pybind11.get_include(),
                        pybind11.get_include(user=True),
                        "src/include",
                        METIS_INC,
                    ],
                    library_dirs=[METIS_LIB],
                    libraries=["metis"],                  # <-- drop "GKlib"
                    language="c++",
                    define_macros=[("IDXTYPEWIDTH", "32"), ("REALTYPEWIDTH", "64")],
                    extra_compile_args=["-std=c++17", "-O3", "-fPIC"],  # <-- GCC/Clang flags
                    extra_link_args=[f"-Wl,-rpath,{METIS_LIB}", f"-Wl,-rpath,{METIS_LIB64}"],        # ensure libmetis is found at runtime
                )
            ] 
        if WINDOWS:

            # Define the extension
            ext_modules = [
                Extension(
                    module_name,
                    sources=[cpp_file],
                    include_dirs=[
                        pybind11.get_include(),
                        os.path.join(vcpkg_root, "installed/x64-windows/include"),
                        "src/include"
                    ],
                    library_dirs=[os.path.join(vcpkg_root, "installed/x64-windows/lib")],
                    libraries=["metis", "GKlib"],
                    language="c++",
                    extra_compile_args=["/std:c++17", "/Zi", "/Od", "/EHsc"],
                )
            ]

        # ext_modules = [
        #     Extension(
        #         module_name,
        #         sources=[cpp_file],
        #         include_dirs=[
        #             pybind11.get_include(),
        #             os.path.join(vcpkg_root, "installed/x64-windows/include"),
        #             "src/include"
        #         ],
        #         library_dirs=[os.path.join(vcpkg_root, "installed/x64-windows/lib")],
        #         libraries=["metis"],
        #         language="c++",
        #         extra_compile_args=["/std:c++17", "/Zi", "/Od", "/EHsc"],
        #     )
        # ]        

        # Dynamically build the extension
        class BuildExt(build_ext):
            def build_extensions(self):
                for ext in self.extensions:
                    if WINDOWS:
                        ext.extra_compile_args = ["/std:c++17", "/Zi", "/Od", "/EHsc"]
                    elif VISTA:
                        ext.extra_compile_args = ["-std=c++17", "-O3", "-fPIC"]
                    else:
                        raise Exception("Unsupported platform")
                    #  ext.extra_compile_args = ["-std=c++17", "-g", "-O0", "-fno-omit-frame-pointer", "-march=native", "-DDEBUG"]
                super().build_extensions()

        setup(
            name=module_name,
            ext_modules=ext_modules,
            cmdclass={"build_ext": BuildExt},
            script_args=["build_ext", "--inplace"],
        )

    # Find the compiled shared object file
    if os.name == 'nt':
        module_files = glob.glob(f"{module_name}*.pyd")
    else:
        module_files = glob.glob(f"{module_name}*.so")

    if not module_files:
        raise ImportError(f"Compiled shared library for {module_name} not found.")
    module_path = module_files[0]

    # Load the module dynamically
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module


if __name__ == "__main__":
    # Get the project root directory (src folder)
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    # Set paths relative to project root
    cpp_file = os.path.join(project_root, "grid", "read_and_process_grid.cpp")

    print("Looking for file at:", os.path.abspath(cpp_file))
    module = compile_and_import("test", cpp_file)


    from iop.readTides import read_bctides

    namelist = readNamelist(in_dir)


    hgrid_file    = os.path.join(in_dir, 'hgrid.gr3')
    hgrid_ll_file = os.path.join(in_dir, 'hgrid.ll')
    vgrid_file    = os.path.join(in_dir, 'vgrid.in')

    
    numPartitions = 1
    i34, partElements, partNodes, partSides, elnode, elside, nodeCoords, \
        nodeCoords_ll, ic3, isdel, isidenode,isbnd, nond, nlnd, iond, ilnd, nodeNeighElems, nodeNeighNodes = \
            module.readHGrid(hgrid_file, hgrid_ll_file, numPartitions, namelist.opt.ncor)
    
    sizeMatrix = nodeCoords.shape[1] # The size of the matrix is equal to the number of nodes.
    totMatrixElements = nodeNeighNodes[0].size + sizeMatrix   # The total number of matrix elements is the number of node neighbors plus the number of nodes.

    # Create the solver
    solver = PCGSolver(sizeMatrix, totMatrixElements)
    
    print(f"Matrix size = {sizeMatrix}, total number of non-zero matrix elements = {totMatrixElements}")
    nope, mnond  = iond.shape
    nland, mnlnd = ilnd.shape
    bctides = read_bctides(in_dir, nope, mnond, [1,1,0,0,0,0,0,0,0,0,0,0], [1,2])




    [ivcor, nvrt, nzlev, nsiglev], ztot, [h_c, h_s, theta_b, theta_f], sigma_lcl, sigma, cs, dcs = \
        module.readVGrid(vgrid_file, nodeCoords.shape[1])

    if (nws == 2):
        wind_rotate_angle = module.readWindRotationAngle(os.path.join(in_dir, 'windrot_geo2proj.gr3'), nodeCoords.shape[1])
    else:
        wind_rotate_angle = np.zeros(nodeCoords.shape[1])



    # Before sending the mesh to the GPU, we need to do some preprocessing, including computing the 
    # interpolation weights.

    # Test the readSflux module
    # nodeCoords[:2,:] = np.deg2rad(nodeCoords[:2,:]) # XXX Why deg2rad?  XXX Why was I doing this? Is important in other contexts???

    sflux = SFlux(in_dir, r'sflux_inputs.txt', nodeCoords) 

    u_air_node = np.zeros(nodeCoords.shape[1])
    v_air_node = np.zeros(nodeCoords.shape[1])
    p_air_node = np.zeros(nodeCoords.shape[1])
    t_air_node = np.zeros(nodeCoords.shape[1])
    q_air_node = np.zeros(nodeCoords.shape[1])
    # sflux.getWind(2452424.0, u_air_node, v_air_node, p_air_node, t_air_node, q_air_node)

   

  
    # Read in bottom friction coefficients 
    nchi = namelist.opt.nchi
    drag      = np.zeros(1)
    roughness = np.zeros(1)
    rmanning  = np.zeros(1)
    if nchi == 0: # read in drag coefficients for 3D model
        drag    = np.zeros(nodeCoords.shape[1])
        drag[:] = module.readHGridData(os.path.join(in_dir, 'drag.gr3'), nodeCoords.shape[1], elnode.shape[1])

    elif nchi == 1:
        roughness    = np.zeros(nodeCoords.shape[1])
        roughness[:] = module.readHGridData(os.path.join(in_dir, 'rough.gr3'), nodeCoords.shape[1], elnode.shape[1])

    elif nchi == -1:
        rmanning   = np.zeros(nodeCoords.shape[1])
        rmanning[:] = module.readHGridData(os.path.join(in_dir, 'manning.gr3'), nodeCoords.shape[1], elnode.shape[1])

    else:
        raise Exception(f"nchi = {nchi} is not a valid option")

        
    # Generate the constants.h file that will be compiled with the CUDA kernels
    processConstants(namelist, bctides, ivcor, nvrt, nzlev, nsiglev, 
                     nope, mnond, nland, mnlnd,
                     templatePath="kernels/constants_template.h", 
                     outputPath="kernels/constants.h")
    
    # Initialize the CUDA functions
    cufunc.initializeCudaKernelManager(project_root)

    
    # convert arrays to cupy
    i34                  = cp.array(i34,               dtype=cp.int32, order='F')
    partElements         = cp.array(partElements,      dtype=cp.int32, order='F')
    partNodes            = cp.array(partNodes,         dtype=cp.int32, order='F')
    partSides            = cp.array(partSides,         dtype=cp.int32, order='F')
    elnode               = cp.array(elnode,            dtype=cp.int32, order='F')
    elside               = cp.array(elside,            dtype=cp.int32, order='F')
    ic3                  = cp.array(ic3,               dtype=cp.int32, order='F')
    isdel                = cp.array(isdel,             dtype=cp.int32, order='F')
    isidenode            = cp.array(isidenode,         dtype=cp.int32, order='F')
    isbnd                = cp.array(isbnd,             dtype=cp.int32, order='F')
    nond                 = cp.array(nond,              dtype=cp.int32, order='F')
    nlnd                 = cp.array(nlnd,              dtype=cp.int32, order='F')
    iond                 = cp.array(iond,              dtype=cp.int32, order='F')
    ilnd                 = cp.array(ilnd,              dtype=cp.int32, order='F')
    nodeNeighNodesStart  = cp.array(nodeNeighNodes[1], dtype=cp.int32, order='F')
    nodeNeighNodes       = cp.array(nodeNeighNodes[0], dtype=cp.int32, order='F')
    nodeNeighElemsStart  = cp.array(nodeNeighElems[1], dtype=cp.int32, order='F')
    nodeNeighElems       = cp.array(nodeNeighElems[0], dtype=cp.int32, order='F')
    nodeCoords           = cp.array(nodeCoords,        dtype=cp.float64, order='F')
    nodeCoords_ll        = cp.array(nodeCoords_ll,     dtype=cp.float64, order='F')
    wind_rotate_angle    = cp.array(wind_rotate_angle, dtype=cp.float64, order='F')

    ztot                 = cp.array(ztot,              dtype=cp.float64, order='F')
    sigma_lcl            = cp.array(sigma_lcl,         dtype=cp.float64, order='F')
    sigma                = cp.array(sigma,             dtype=cp.float64, order='F')
    cs                   = cp.array(cs,                dtype=cp.float64, order='F')
    dcs                  = cp.array(dcs,               dtype=cp.float64, order='F')

    drag                 = cp.array(drag,              dtype=cp.float64, order='F')
    roughness            = cp.array(roughness,         dtype=cp.float64, order='F')
    rmanning             = cp.array(rmanning,          dtype=cp.float64, order='F')


    if False:   

        fig, ax = plt.subplots(figsize=(12, 8))
        
        plotLabels = False

        # Plot all nodes
        plt.plot(nodeCoords[0,:].get(), nodeCoords[1,:].get(), '.', alpha=0.3, markersize=.3)

        if plotLabels:
            # Add node labels
            for i in range(nodeCoords.shape[1]):
                plt.annotate(f'N{i}', 
                            (nodeCoords[0,i].get(), nodeCoords[1,i].get()),
                            xytext=(2, 2),
                            textcoords='offset points',
                            color='black',
                            fontsize=6,
                            alpha=0.5)

        nodes = elnode[:3,:].get()
        # Append the first node to close the triangle
        closed_nodes = np.vstack((nodes, nodes[0,:]))
        
        # Get coordinates for all nodes
        x1 = nodeCoords[0, nodes[0, :]].get()
        y1 = nodeCoords[1, nodes[0, :]].get()
        x2 = nodeCoords[0, nodes[1, :]].get()
        y2 = nodeCoords[1, nodes[1, :]].get()
        x3 = nodeCoords[0, nodes[2, :]].get()
        y3 = nodeCoords[1, nodes[2, :]].get()
        
        # Plot element edges
        points = np.stack((np.column_stack((x1, y1)), np.column_stack((x2, y2))), axis=1)
        lc = LineCollection(points, colors='r', alpha=0.1)
        ax.add_collection(lc)
        points = np.stack((np.column_stack((x2, y2)), np.column_stack((x3, y3))), axis=1)
        lc = LineCollection(points, colors='r', alpha=0.1)
        ax.add_collection(lc)
        points = np.stack((np.column_stack((x3, y3)), np.column_stack((x1, y1))), axis=1)
        lc = LineCollection(points, colors='r', alpha=0.1)
        ax.add_collection(lc)



        if plotLabels:
            # Add element labels at centroids
            for i in range(elnode.shape[1]):
                # Get triangle nodes
                tri_nodes = elnode[:3,i].get()
                
                # Calculate triangle centroid
                centroid_x = np.mean(nodeCoords[0,tri_nodes].get())
                centroid_y = np.mean(nodeCoords[1,tri_nodes].get())
                
                # Element label at centroid
                plt.annotate(f'E{i}', 
                            (centroid_x, centroid_y),
                            xytext=(0, 0),
                            textcoords='offset points',
                            color='red',
                            fontsize=6,
                            alpha=0.6)

            # Add side labels
            for i in range(isidenode.shape[1]):
                # Get side nodes
                side_nodes = isidenode[:,i].get()
                
                # Calculate side midpoint
                mid_x = np.mean(nodeCoords[0,side_nodes].get())
                mid_y = np.mean(nodeCoords[1,side_nodes].get())
                
                # Side label at midpoint
                plt.annotate(f'S{i}',
                            (mid_x, mid_y),
                            xytext=(0, 3),
                            textcoords='offset points',
                            color='blue',
                            fontsize=6,
                            alpha=0.3)

            
        # Prevent plot distortion by setting equal aspect ratio
        plt.axis('equal')
        
        plt.show()
        plt.pause(0)

    # js         = cp.array(js)
    


    # Get number of GPUs
    n_gpus = cp.cuda.runtime.getDeviceCount()
    print(f"Found {n_gpus} GPUs")


    cufunc.initTidalHarmonicConstituents(bctides.amig, bctides.face, bctides.ff)

    # This list stores references to the device memory of the local arrays. This ensures that the memory is not freed until the mesh is destroyed.
    deviceMemRef = [ztot, sigma_lcl, sigma, cs, dcs]

    numNodes = nodeCoords.shape[1]
    numElements = elnode.shape[1]
    numSides = isdel.shape[1]  
 

    
    # Get the struct sizes
    sizes = cp.zeros(6, dtype=cp.int32)
    cufunc.getStructSizes(sizes)
    [nodesSize, elemsSize, sidesSize, boundaryNodesSize, sideBoundarySize, segmentSize] = sizes.tolist()



    
    # Compute the maximum number of nodes, elements, and sides per partition
    numPartElements = []
    numPartNodes = []
    numPartSides = []
    numPartBoundaryNodes = []
    numPartBoundarySides = []
    for ipart, device_id in enumerate(range(n_gpus)):

        # Count the number of elements in the partition
        numPartElements.append((partElements == ipart).sum().get())
        print(f"Partition {ipart} has {numPartElements} elements")

        # Count the number of nodes in the partition
        numPartNodes.append((partNodes == ipart).sum().get())
        print(f"Partition {ipart} has {numPartNodes} nodes")

        # Count the number of sides in the partition
        numPartSides.append((partSides == ipart).sum().get())
        print(f"Partition {ipart} has {numPartSides} sides")
       
        # Count the number of boundary nodes in the partition
        # idxBndNodes = cp.unique(cp.concatenate([ilnd.flatten(), iond.flatten()]))
        # XXXXX The whole computation of the number of boundary nodes and sides is wrong. 
        print("REMEMBER TO FIX THIS*****")
        idxBndNodes = cp.append(cp.unique(iond.flatten()), cp.unique(ilnd.flatten()))
        idxBndNodes = idxBndNodes[idxBndNodes >= 0]
        numPartBoundaryNodes.append((partNodes[idxBndNodes] == ipart).sum().get())
        print(f"Partition {ipart} has {numPartBoundaryNodes} boundary nodes")
       
        # Count the number of boundary sides in the partition
        numPartBoundarySides.append((partSides[idxBndNodes] == ipart).sum().get())
        # numPartBoundarySides.append((partSides[cp.unique(ilnd)] == ipart).sum().get())
        print(f"Partition {ipart} has {numPartBoundarySides} boundary sides")


    monitor.init(numNodes, numElements, numSides, numPartBoundarySides, nvrt) # XXX not valid for multiGPU
       
    nodesMemList = []
    elemsMemList = []
    sidesMemList = []
    boundaryNodesMemList = []
    boundarySidesMemList = []
    segmentMemList       = []
    for ipart, device_id in enumerate(range(n_gpus)):
        with cp.cuda.Device(device_id):

            # Create the local device arrays.
            i34_part    = cp.array(i34[partElements == ipart])
            elnode_part = cp.array(elnode[:, partElements == ipart])

            # Reserves GPU memory for the structs
            nodesMem = cp.zeros(max(numPartNodes   ).item()*nodesSize, dtype=cp.int8)
            elemsMem = cp.zeros(max(numPartElements).item()*elemsSize, dtype=cp.int8)
            sidesMem = cp.zeros(max(numPartSides   ).item()*sidesSize, dtype=cp.int8)

            # Memory for boundary nodes objects in each partition (We store in each device the BoundaryNodes objects 
            # corresponding to the nodes in the partition)
            boundaryNodesMem = cp.zeros(max(numPartBoundaryNodes)*boundaryNodesSize, dtype=cp.int8)

            # We don't know the number of boundary sides in each partition, so we reserve as much as the number of boundary nodes.
            boundarySidesMem  = cp.zeros((max(numPartBoundarySides)+1)*sideBoundarySize,  dtype=cp.int8)

            # We save all the segment info in all partitions.
            segmentMem = cp.zeros(ilnd.shape[0]*segmentSize, dtype=cp.int8)

            # Need to store the pointers to be used later in createMesh and populateMesh
            boundaryNodesMemList.append(boundaryNodesMem.data.ptr)
            boundarySidesMemList.append(boundarySidesMem.data.ptr)
            segmentMemList      .append(segmentMem.data.ptr)

            nodesMemList.append(nodesMem.data.ptr)
            elemsMemList.append(elemsMem.data.ptr)
            sidesMemList.append(sidesMem.data.ptr)

            # Add the arrays to the device memory reference list. This ensures that the memory is not freed until the mesh is destroyed.
            deviceMemRef.append([i34_part, elnode_part, nodesMem, elemsMem, sidesMem, boundaryNodesMem, boundarySidesMem, segmentMem])

    # Create the arrays of pointers to the start of the device memory
    nodesPtr = cp.array(nodesMemList, dtype=cp.uintp)
    elemsPtr = cp.array(elemsMemList, dtype=cp.uintp)
    sidesPtr = cp.array(sidesMemList, dtype=cp.uintp)

    partNodeNeighNodesMem = []
    partNodeNeighElemsMem = []
    partCoeffMatrixMem    = []

    # Create the local node and element neighbor arrays
    for ipart, device_id in enumerate(range(n_gpus)):
        with cp.cuda.Device(device_id):

            maskNodesPart = (partNodes == ipart)
            maskElemsPart = (partElements == ipart)

            partNodeNeighNodesSize = maskNodesPart.sum().item()
            partNodeNeighElemsSize = maskElemsPart.sum().item()

            # compute the local size of the node and element neighbor arrays  XXXX is the [maskNodesPart] below correct?
            partNodeNeighNodesCount = (nodeNeighNodesStart[1:][maskNodesPart] - nodeNeighNodesStart[:-1][maskNodesPart])
            partNodeNeighElemsCount = (nodeNeighElemsStart[1:][maskNodesPart] - nodeNeighElemsStart[:-1][maskNodesPart])

            pointerSize = cp.dtype(cp.uintp).itemsize 
            partNodeNeighNodesMem += [cp.zeros(partNodeNeighNodesCount.sum().item()*pointerSize, dtype=cp.int8)]
            partNodeNeighElemsMem += [cp.zeros(partNodeNeighElemsCount.sum().item()*pointerSize, dtype=cp.int8)]
            partCoeffMatrixMem    += [cp.zeros(partNodeNeighNodesCount.sum().item(), dtype=cp.float64)]
    print(f"SIZE OF partCoeffMatrixMem = {partCoeffMatrixMem[0].size}") # XXXXXX REMOVE
    print(f"zeros in partNodeNeighNodesCount = {(partNodeNeighNodesCount == 0).sum().get()}")
    print(f"max in partNodeNeighNodesCount = {partNodeNeighNodesCount.max().get()}")
    print(f"min in partNodeNeighNodesCount = {partNodeNeighNodesCount.min().get()}")
 


    
    # Create arrays that convert from global indices to (partition,localIdx)
    # These arrays will be released after building the mesh
    with cp.cuda.Device(0):
        localNodeIndices = cp.zeros(numNodes,    dtype=cp.int32)
        localElemIndices = cp.zeros(numElements, dtype=cp.int32)
        localSideIndices = cp.zeros(numSides,    dtype=cp.int32)
        localNodePart    = cp.zeros(numNodes,    dtype=cp.int32)
        localElemPart    = cp.zeros(numElements, dtype=cp.int32)
        localSidePart    = cp.zeros(numSides,    dtype=cp.int32)
        
        for ipart in range(n_gpus):
            
            maskNodesPart = (partNodes == ipart)
            maskElemsPart = (partElements == ipart)
            maskSidesPart = (partSides == ipart)

            # Given a global index, this array will return the local index
            localNodeIndices[maskNodesPart] = cp.arange(numPartNodes[ipart])
            localElemIndices[maskElemsPart] = cp.arange(numPartElements[ipart])
            localSideIndices[maskSidesPart] = cp.arange(numPartSides[ipart])

            # Given a global index, this array will return the partition index
            localNodePart   [maskNodesPart] = ipart
            localElemPart   [maskElemsPart] = ipart
            localSidePart   [maskSidesPart] = ipart

    for ipart in range(n_gpus):
        with cp.cuda.Device(ipart):

            maskNodesPart = (partNodes == ipart)
            maskElemsPart = (partElements == ipart)
            nodeIndicesPart = cp.where(maskNodesPart)[0].astype(cp.int32)

            i34Part    = i34[maskElemsPart].copy(order='F')
            elnodePart = elnode[:, maskElemsPart].copy(order='F')
            elsidePart = elside[:, maskElemsPart].copy(order='F')
            nodeCoordsPart      = nodeCoords[:, maskNodesPart].copy(order='F')
            windRotateAnglePart = wind_rotate_angle[maskNodesPart].copy(order='F')
            nodeCoords_llPart = cp.zeros((0, 0), dtype=cp.float64, order='F')
            if nodeCoords_ll.shape[0] > 0:
                nodeCoords_llPart  = nodeCoords_ll[:, maskNodesPart].copy(order='F')
                use_ll = True
            else:
                use_ll = False

            numSegments = ilnd.shape[0]

            # Creates the mesh object on each GPU. This function mostly deals with letting the GPU know
            # which memory is being used by the mesh, both locally and in other GPUs.
            cufunc.createMesh(ipart, numNodes, numElements, numSides, numSegments,
                              numPartNodes[ipart], numPartElements[ipart], numPartSides[ipart],  
                              max(numPartNodes), max(numPartElements), max(numPartSides),
                              nodesPtr, elemsPtr, sidesPtr,
                              nodesMemList[ipart], elemsMemList[ipart], sidesMemList[ipart],
                              partNodeNeighNodesMem[ipart], partNodeNeighElemsMem[ipart], 
                              partCoeffMatrixMem[ipart],
                              boundaryNodesMemList[ipart], boundarySidesMemList[ipart], segmentMemList[ipart],
                              h_c.item(), h_s.item(), theta_b.item(), theta_f.item(), sigma, sigma, ztot, cs, dcs,
                              monitor.debugPlot.data.ptr, monitor.debugPlotIdx.data.ptr)
            
            # Given a local index, this array will return the global index
            globalNodeIndices = cp.arange(numPartNodes   [ipart]).astype(cp.int32)[maskNodesPart]
            globalElemIndices = cp.arange(numPartElements[ipart]).astype(cp.int32)[maskElemsPart]
            globalSideIndices = cp.arange(numPartSides   [ipart]).astype(cp.int32)[maskSidesPart]

    sizes = {'global': 1, 'nodes': max(numPartNodes), 'elements': max(numPartElements), 'sides': max(numPartSides),
             'nodeBnd': max(numPartBoundaryNodes), 'sideBnd': max(numPartBoundarySides)}
    cufunc.initializeCudaFunctions(sizes)

    for ipart in range(n_gpus):
        with cp.cuda.Device(ipart):
            # Adds additional information to the mesh object. This function deals with the actual data
            # of the mesh, including the nodes, elements, sides, and their neighbors.
            sharedIdxNodeBoundary  = cp.zeros(1, dtype=cp.int32) - 1 # the first increment makes it 0
            sharedIdxSideBoundary  = cp.zeros(1, dtype=cp.int32) - 1 

            cufunc.populateNodes(i34Part, isbnd, sharedIdxNodeBoundary, 
                       nodeCoordsPart, nodeCoords_llPart, use_ll, windRotateAnglePart, 
                       nodeNeighNodesStart, nodeNeighNodes,
                       nodeNeighElemsStart, nodeNeighElems,
                       localNodeIndices, localElemIndices, localSideIndices,
                       localNodePart, localElemPart, localSidePart,
                       drag, roughness, rmanning)

            cufunc.populateSides(i34Part, sharedIdxSideBoundary,
                       nodeCoordsPart, nodeCoords_llPart, use_ll,
                       windRotateAnglePart,   
                       isdel, isidenode, isbnd,
                       nond, nlnd, iond, ilnd,                                  
                       nodeNeighNodesStart, nodeNeighNodes,
                       nodeNeighElemsStart, nodeNeighElems,
                       globalNodeIndices, globalElemIndices, globalSideIndices,
                       localNodeIndices, localElemIndices, localSideIndices,
                       localNodePart, localElemPart, localSidePart,
                       drag, roughness);        

            cufunc.populateElements(i34Part, elnodePart, elsidePart,
                       nodeCoordsPart, nodeCoords_llPart, use_ll,
                       windRotateAnglePart,   
                       isdel, isidenode, isbnd,
                       nond, nlnd, iond, ilnd,                                  
                       nodeNeighNodesStart, nodeNeighNodes,
                       nodeNeighElemsStart, nodeNeighElems,
                       globalNodeIndices, globalElemIndices, globalSideIndices,
                       localNodeIndices, localElemIndices, localSideIndices,
                       localNodePart, localElemPart, localSidePart,
                       drag, roughness);                                             


            cufunc.populateMesh(i34Part, elnodePart, elsidePart, nodeCoordsPart, nodeCoords_llPart, use_ll,
                                windRotateAnglePart, sharedIdxNodeBoundary, sharedIdxSideBoundary,
                                isdel, isidenode, isbnd,
                                nond, nlnd, iond, ilnd, 
                                nodeNeighNodesStart, nodeNeighNodes,
                                nodeNeighElemsStart, nodeNeighElems,
                                globalNodeIndices, globalElemIndices, globalSideIndices,
                                localNodeIndices, localElemIndices, localSideIndices,
                                localNodePart, localElemPart, localSidePart,
                                drag, roughness,
                                bctides.emo, bctides.efa, 
                                bctides.umo, bctides.ufa,
                                bctides.vmo, bctides.vfa)

            cufunc.populateBoundaryNodes(nond, iond,
                                         bctides.emo, bctides.efa, 
                                         bctides.umo, bctides.ufa,
                                         bctides.vmo, bctides.vfa)                                





    ntime, nrec = schism_init(namelist.core, namelist.opt)

    initializeOutput(out_dir, numNodes, numElements, numSides, nvrt, nodesSize, elemsSize, sidesSize)


    profiler = LineProfiler() 
    schism_step_wrapped = profiler(schism_step)
    profiler.add_function(solveWaveEquation)
    


    start_time = time.time()
    time1 = time.time()
    for it in range(1, ntime + 1, 1):
    # for it in range(1, 55, 1):
        print(f"> Step {it}, simulation time: {it*namelist.core.dt:.2f} s")
        # schism_step_wrapped(it, namelist, bctides, sflux, solver)

        schism_step(it, namelist, bctides, sflux, solver)

        writeOutput(it, namelist.core, namelist.opt)
        # if it % 100 == 0:
        #     monitor.plotSideVar("u", 0)

        time2 = time.time()
        print(f"Time per step: {time2 - time1:.3f} s")
        time1 = time2
        pass


    end_time = time.time()
    total_time = end_time - start_time
    print(f"\nTotal execution time: {total_time:.3f} s")
    print(f"Average time per step: {(total_time/ntime):.3f} s")

    profiler.print_stats()    
    pass


    

    