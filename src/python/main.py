import os
import glob
import importlib.util
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import cupy as cp
import debug.monitor as monitor
import time
import pybind11
from python.pcgsolver import PCGSolver, PCGSolverMultiGPU
from python.schism_consts import VISTA, WINDOWS
from iop.writeOutput import writeOutput, initializeOutput, getOutputTime
import iop.cudaKernelWrapper as ckw
import iop.deviceToHost.DeviceToHost as DeviceToHost
import iop.deviceToHost.getNodeVar as getNodeVar
import iop.deviceToHost.getElemVar as getElemVar
import iop.deviceToHost.getSideVar as getSideVar

pool = cp.cuda.MemoryPool()
cp.cuda.set_allocator(pool.malloc)

# For profiling
from line_profiler import LineProfiler





# Enable verbose output
os.environ['DISTUTILS_DEBUG'] = '1'
os.environ['VERBOSE'] = '1'


# # Print diagnostic information
# print("Python version:", sys.version)
# print("Include paths:")
# for path in [
#     r'C:\Users\jmgonzalez\AppData\Local\miniconda3\envs\gschism\Lib\site-packages\pybind11\include',
#     r'C:\Users\jmgonzalez\apps\vcpkg\installed\x64-windows\include',
#     'src/include'
# ]:
#     print(f"  {path} exists: {os.path.exists(path)}")


import numpy as np
import timeit

import matplotlib
# matplotlib.use('TkAgg') 
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection




sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'src'))

from iop.cudamodule import CudaSourceManager
from iop.processConstantsHeader import processConstants
from python.schism_init import readNamelist, schism_init
from python.schism_step import schism_step, solveWaveEquation, computeMeteorologicalForcing
from python.schism_consts import *
import python.cuda_functions as cufunc




# # in_dir = '/media/jmg/Elements/src/tmp/schism_verification_tests/Test_Btrack_Gausshill'
# # in_dir = '/media/jmg/Elements/src/tmp/schism_verification_tests/Small_test/Small_grid'
# in_dir = r'C:\Users\jmgonzalez\src\schism_test_suit\Small_grid'
# in_dir = r'C:\Users\jmgonzalez\src\schism_test_suit\QuarterAnnulus\Pure_tri'
# # in_dir = r'C:\Users\jmgonzalez\src\schism_test_suit\RUN02a'


if WINDOWS:
    # in_dir = r'C:\Users\jmgonzalez\src\schism_test_suit\RUN03b'
    # out_dir = r'C:\Users\jmgonzalez\src\schism_test_suit\RUN03b\gpuoutput'
    in_dir = r'C:\Users\jmgonzalez\src\schism_test_suit\RUN02e'
    out_dir = r'C:\Users\jmgonzalez\src\schism_test_suit\RUN02e\gpuoutput'
    # in_dir = r'C:\Users\jmgonzalez\src\schism_test_suit\Test_CORIE_LSC2\Test_CORIE_LSC2\GPU'
    # out_dir = r'C:\Users\jmgonzalez\src\schism_test_suit\Test_CORIE_LSC2\Test_CORIE_LSC2\GPU\gpuoutput'

    # out_dir = r'c:\Users\jmgonzalez\src\temp'
elif VISTA:
    in_dir = r'/scratch/10922/jmglez/tests/RUN02i'
    out_dir = r'/scratch/10922/jmglez/tests/RUN02i/gpuoutput'
    # in_dir = r'/scratch/10922/jmglez/tests/RUN02a'
    # out_dir = r'/scratch/10922/jmglez/tests/RUN02a/gpuoutput'
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

    # Initialize the CUDA functions

    configCuda = ckw.initCudaDeviceConfig()


    numPartitions = ckw.configCuda.deviceCount
    i34, partElements, partNodes, partSides, elnode, elside, nodeCoords, \
        nodeCoords_ll, ic3, isdel, isidenode,isbnd, nond, nlnd, iond, ilnd, nodeNeighElems, nodeNeighNodes, \
        numPartElements, numPartNodes, numPartSides, numPartBoundaryNodes, numPartBoundarySides = \
            module.readHGrid(hgrid_file, hgrid_ll_file, numPartitions, namelist.opt.ncor)
    

    # Initialize the solver
    if numPartitions == 1:
        N = nodeCoords.shape[1] # The size of the matrix is equal to the number of nodes.
        NNZ = nodeNeighNodes[0].size + N   # The total number of non-zero matrix elements is the number of node neighbors plus the number of nodes.
        solver = PCGSolver(N, NNZ)
        print(f"Matrix size = {N}, total number of non-zero matrix elements = {NNZ}")

    else:
        # Count nonzeros per row-block A_p (rows owned by partition p, columns = all nodes packed).
        idxStarts = nodeNeighNodes[1]
        idxNeigh = nodeNeighNodes[0]

        # ipart[k] = partition of the row node for the k-th off-diagonal neighbor entry
        ipart = np.repeat(partNodes, np.diff(idxStarts))

        NNZ = np.zeros(numPartitions, dtype=np.int64)
        np.add.at(NNZ, ipart, 1)  # off-diagonals whose row is in partition p
        del ipart

        N = [int(np.sum(partNodes == part)) for part in range(numPartitions)]
        NNZ += np.asarray(numPartNodes, dtype=np.int64)  # + diagonal per owned row

        if NNZ.sum() != idxNeigh.size + nodeCoords.shape[1]:
            raise Exception(
                f"The number of non-zero matrix elements is not correct. "
                f"Expected {idxNeigh.size + nodeCoords.shape[1]}, got {NNZ.sum()}"
            )

        solver = PCGSolverMultiGPU(ckw.configCuda, N, NNZ.tolist())
        print(f"Row-block NNZ = {NNZ.tolist()}, N = {N}")


    
    
    
    
    nope, mnond  = iond.shape
    nland, mnlnd = ilnd.shape
    bctides = read_bctides(in_dir, nope, mnond, [1,1,0,0,0,0,0,0,0,0,0,0], [1,2])




    [ivcor, nvrt, nzlev, nsiglev], ztot, kb, [h_c, h_s, theta_b, theta_f], sigma_lcl, sigma, cs, dcs = \
        module.readVGrid(vgrid_file, nodeCoords.shape[1])

    # If nw2 is 2, read the wind rotation angle from the windrot_geo2proj.gr3 file. Otherwise, set it to 0.
    if (namelist.opt.nws == 2):
        wind_rotate_angle = module.readWindRotationAngle(os.path.join(in_dir, 'windrot_geo2proj.gr3'), nodeCoords.shape[1])
    else:
        wind_rotate_angle = np.zeros(nodeCoords.shape[1])



    # Before sending the mesh to the GPU, we need to do some preprocessing, including computing the 
    # interpolation weights.

    # Test the readSflux module
    # nodeCoords[:2,:] = np.deg2rad(nodeCoords[:2,:]) # XXX Why deg2rad?  XXX Why was I doing this? Is important in other contexts???

    sflux = SFlux(in_dir, r'sflux_inputs.txt', nodeCoords_ll) 

    u_air_node = np.zeros(nodeCoords.shape[1])
    v_air_node = np.zeros(nodeCoords.shape[1])
    p_air_node = np.zeros(nodeCoords.shape[1])
    t_air_node = np.zeros(nodeCoords.shape[1])
    q_air_node = np.zeros(nodeCoords.shape[1])
  

  
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
                     numPartNodes,
                     templatePath="kernels/constants_template.h", 
                     outputPath="kernels/constants.h")
    
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
    kb                   = cp.array(kb,                dtype=cp.int32,   order='F')
    sigma                = cp.array(sigma,             dtype=cp.float64, order='F')
    cs                   = cp.array(cs,                dtype=cp.float64, order='F')
    dcs                  = cp.array(dcs,               dtype=cp.float64, order='F')

    drag                 = cp.array(drag,              dtype=cp.float64, order='F')
    roughness            = cp.array(roughness,         dtype=cp.float64, order='F')
    rmanning             = cp.array(rmanning,          dtype=cp.float64, order='F')



   

    print(f"Found {ckw.configCuda.deviceCount} GPUs")
    print(">>>> 1")
    cufunc.initTidalHarmonicConstituents(bctides.amig, bctides.face, bctides.ff)
    print(">>>> 2")
    # This list stores references to the device memory of the local arrays. This ensures that the memory is not freed until the mesh is destroyed.
    deviceMemRef = [ztot, sigma, cs, dcs]


    numNodes = nodeCoords.shape[1]
    numElements = elnode.shape[1]
    numSides = isdel.shape[1]  
    print(f">>>>> 3 numNodes = {numNodes}, numElements = {numElements}, numSides = {numSides}")

    
    # Get the struct sizes
    sizes = cp.zeros(6, dtype=cp.int32)
    cufunc.getStructSizes(sizes)
    [nodesSize, elemsSize, sidesSize, boundaryNodesSize, sideBoundarySize, segmentSize] = sizes.tolist()
    print(f">>>>> 4 nodesSize = {nodesSize}, elemsSize = {elemsSize}, sidesSize = {sidesSize}, boundaryNodesSize = {boundaryNodesSize}, sideBoundarySize = {sideBoundarySize}, segmentSize = {segmentSize}")


    # XXXXX Remember to remove
    # Compute the maximum number of nodes, elements, and sides per partition
    # numPartElements = []
    # numPartNodes = []
    # numPartSides = []
    # numPartBoundaryNodes = []
    # numPartBoundarySides = [] 
    # for ipart, device_id in enumerate(ckw.configCuda.devices):

    #     # Count the number of elements in the partition
    #     numPartElements.append((partElements == ipart).sum().get())
    #     print(f"Partition {ipart} has {numPartElements} elements")

    #     # Count the number of nodes in the partition
    #     numPartNodes.append((partNodes == ipart).sum().get())
    #     print(f"Partition {ipart} has {numPartNodes} nodes")

    #     # Count the number of sides in the partition
    #     numPartSides.append((partSides == ipart).sum().get())
    #     print(f"Partition {ipart} has {numPartSides} sides")
       
    #     # Count the number of boundary nodes in the partition
    #     # idxBndNodes = cp.unique(cp.concatenate([ilnd.flatten(), iond.flatten()]))
    #     # XXXXX The whole computation of the number of boundary nodes and sides is wrong. 
    #     print("REMEMBER TO FIX THIS*****")
    #     idxBndNodes = cp.append(cp.unique(iond.flatten()), cp.unique(ilnd.flatten()))
    #     idxBndNodes = idxBndNodes[idxBndNodes >= 0]
    #     numPartBoundaryNodes.append((partNodes[idxBndNodes] == ipart).sum().get())
    #     print(f"Partition {ipart} has {numPartBoundaryNodes} boundary nodes")
       
    #     # Count the number of boundary sides in the partition
    #     numPartBoundarySides.append((partSides[idxBndSides] == ipart).sum().get())
    #     # numPartBoundarySides.append((partSides[cp.unique(ilnd)] == ipart).sum().get())
    #     print(f"Partition {ipart} has {numPartBoundarySides} boundary sides")


    monitor.init(numNodes, numElements, numSides, numPartBoundarySides, nvrt) # XXX not valid for multiGPU
    print(f">>>>> 5 monitor initialized")
       
    nodesMemList = []
    elemsMemList = []
    sidesMemList = []
    boundaryNodesMemList = []
    boundarySidesMemList = []
    segmentMemList       = []
    for ipart, device_id in enumerate(ckw.configCuda.devices):
        with cp.cuda.Device(device_id):
            print(f">>>>> 6 device {device_id} initialized")
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

            # Memory for the land boundary nodes and segments.
            landNodeBndMem = cp.zeros(1 * boundaryNodesSize, dtype=cp.int8)
            landSegMem     = cp.zeros(1 * segmentSize, dtype=cp.int8)

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

    partNodeNeighNodesMem      = []
    partNodeNeighElemsMem      = []
    partNodePermutationSortMem = []
    partCoeffMatrixMem         = []

    # Create the local node and element neighbor arrays
    for ipart, device_id in enumerate(ckw.configCuda.devices):
        with cp.cuda.Device(device_id):

            maskNodesPart = (partNodes == ipart)
            maskElemsPart = (partElements == ipart)

            partNodeNeighNodesSize = maskNodesPart.sum().item()
            partNodeNeighElemsSize = maskElemsPart.sum().item()

            # compute the local size of the node and element neighbor arrays  XXXX is the [maskNodesPart] below correct?
            partNodeNeighNodesCount = (nodeNeighNodesStart[1:][maskNodesPart] - nodeNeighNodesStart[:-1][maskNodesPart])
            partNodeNeighElemsCount = (nodeNeighElemsStart[1:][maskNodesPart] - nodeNeighElemsStart[:-1][maskNodesPart])


            pointerSize = cp.dtype(cp.uintp).itemsize 
            partNodeNeighNodesMem      += [cp.zeros(partNodeNeighNodesCount.sum().item()*pointerSize, dtype=cp.int8)]
            partNodeNeighElemsMem      += [cp.zeros(partNodeNeighElemsCount.sum().item()*pointerSize, dtype=cp.int8)]
            partNodePermutationSortMem += [cp.zeros(partNodeNeighNodesCount.sum().item(), dtype=cp.int16)]
            partCoeffMatrixMem         += [cp.zeros(partNodeNeighNodesCount.sum().item(), dtype=cp.float64)]

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
        
        for ipart in range(ckw.configCuda.deviceCount):
            
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

    # Creates the mesh object on each GPU. This function mostly deals with letting the GPU know
    # which memory is being used by the mesh, both locally and in other GPUs.
    for ipart, device_id in enumerate(ckw.configCuda.devices):
        with cp.cuda.Device(device_id):
            numSegments = ilnd.shape[0]
            cufunc.createMesh(ipart, numNodes, numElements, numSides, numSegments,
                              numPartNodes[ipart], numPartElements[ipart], numPartSides[ipart],  
                              max(numPartNodes), max(numPartElements), max(numPartSides),
                              landNodeBndMem, landSegMem,
                              nodesPtr, elemsPtr, sidesPtr,
                              nodesMemList[ipart], elemsMemList[ipart], sidesMemList[ipart],
                              partNodeNeighNodesMem[ipart], partNodeNeighElemsMem[ipart], 
                              partNodePermutationSortMem[ipart],
                              partCoeffMatrixMem[ipart],
                              boundaryNodesMemList[ipart], boundarySidesMemList[ipart], segmentMemList[ipart],
                              kb, h_c.item(), h_s.item(), theta_b.item(), theta_f.item(), sigma, ztot, cs, dcs,
                              monitor.debugPlot.data.ptr, monitor.debugPlotIdx.data.ptr)



    sizes = {'global': 1, 'nodes': max(numPartNodes), 'elements': max(numPartElements), 'sides': max(numPartSides),
             'nodeBnd': max(numPartBoundaryNodes), 'sideBnd': max(numPartBoundarySides)}
    cufunc.initializeCudaFunctions(sizes)


    # Per-partition populate inputs (built once)
    partPopulate = []
    for ipart, device_id in enumerate(ckw.configCuda.devices):
        with cp.cuda.Device(device_id):
            maskNodesPart = (partNodes == ipart)
            maskElemsPart = (partElements == ipart)
            maskSidesPart = (partSides == ipart)
            nodeCoords_llPart = cp.zeros((0, 0), dtype=cp.float64, order='F')
            use_ll = nodeCoords_ll.shape[0] > 0
            if use_ll:
                nodeCoords_llPart = nodeCoords_ll[:, maskNodesPart].copy(order='F')
            partNodeNeighNodesCount = (
                nodeNeighNodesStart[1:][maskNodesPart] - nodeNeighNodesStart[:-1][maskNodesPart]
            )
            partNodeNeighElemsCount = (
                nodeNeighElemsStart[1:][maskNodesPart] - nodeNeighElemsStart[:-1][maskNodesPart]
            )
            localNodeNeighNodesStart = cp.zeros(partNodeNeighNodesCount.size + 1, dtype=cp.int32)
            localNodeNeighNodesStart[1:] = cp.cumsum(partNodeNeighNodesCount)
            localNodeNeighElemsStart = cp.zeros(partNodeNeighElemsCount.size + 1, dtype=cp.int32)
            localNodeNeighElemsStart[1:] = cp.cumsum(partNodeNeighElemsCount)
            globalNodeIndices = cp.arange(maskNodesPart.size, dtype=cp.int32)[maskNodesPart]
            globalElemIndices = cp.arange(maskElemsPart.size, dtype=cp.int32)[maskElemsPart]
            globalSideIndices = cp.arange(maskSidesPart.size, dtype=cp.int32)[maskSidesPart]
            if (globalNodeIndices.size != numPartNodes[ipart]
                    or globalElemIndices.size != numPartElements[ipart]
                    or globalSideIndices.size != numPartSides[ipart]):
                raise ValueError(f"Partition {ipart}: global index sizes do not match partition counts")
            partPopulate.append({
                'device_id': device_id,
                'sharedIdxNodeBoundary': cp.zeros(1, dtype=cp.int32) - 1,
                'sharedIdxSideBoundary': cp.zeros(1, dtype=cp.int32) - 1,
                'i34Part': i34[maskElemsPart].copy(order='F'),
                'elnodePart': elnode[:, maskElemsPart].copy(order='F'),
                'elsidePart': elside[:, maskElemsPart].copy(order='F'),
                'nodeCoordsPart': nodeCoords[:, maskNodesPart].copy(order='F'),
                'nodeCoords_llPart': nodeCoords_llPart,
                'use_ll': use_ll,
                'windRotateAnglePart': wind_rotate_angle[maskNodesPart].copy(order='F'),
                'localNodeNeighNodesStart': localNodeNeighNodesStart,
                'localNodeNeighElemsStart': localNodeNeighElemsStart,
                'globalNodeIndices': globalNodeIndices,
                'globalElemIndices': globalElemIndices,
                'globalSideIndices': globalSideIndices,
            })

    # Phase 1: all nodes on all GPUs
    for ipart, pd in enumerate(partPopulate):
        with cp.cuda.Device(pd['device_id']):
            print(f"HERE device {ipart} populateNodes")
            cufunc.populateNodes(
                pd['i34Part'], isbnd, pd['sharedIdxNodeBoundary'],
                pd['nodeCoordsPart'], pd['nodeCoords_llPart'], pd['use_ll'], pd['windRotateAnglePart'],
                nodeNeighNodesStart, nodeNeighNodes,
                nodeNeighElemsStart, nodeNeighElems,
                partNodePermutationSortMem[ipart], solver.csr_data(ipart),
                pd['localNodeNeighNodesStart'], pd['localNodeNeighElemsStart'],
                pd['globalNodeIndices'], pd['globalElemIndices'], pd['globalSideIndices'],
                localNodeIndices, localElemIndices, localSideIndices,
                localNodePart, localElemPart, localSidePart,
                drag, roughness, rmanning, kb, sigma_lcl)

    
    for ipart, pd in enumerate(partPopulate):
        print(f"HERE device {ipart} populateNodes2")
        with cp.cuda.Device(pd['device_id']):
            cufunc.populateNodes2(
                    pd['localNodeNeighNodesStart'],
                    solver.csr_indices(ipart),
                    solver.csr_indptr(ipart))

    # Phase 2: all sides
    for ipart, pd in enumerate(partPopulate):
        with cp.cuda.Device(pd['device_id']):
            print(f"HERE device {ipart} populateSides")
            cufunc.populateSides(
                pd['i34Part'], pd['sharedIdxSideBoundary'],
                pd['nodeCoordsPart'], pd['nodeCoords_llPart'], pd['use_ll'],
                pd['windRotateAnglePart'],
                isdel, isidenode, isbnd,
                nond, nlnd, iond, ilnd,
                nodeNeighNodesStart, nodeNeighNodes,
                nodeNeighElemsStart, nodeNeighElems,
                pd['globalNodeIndices'], pd['globalElemIndices'], pd['globalSideIndices'],
                localNodeIndices, localElemIndices, localSideIndices,
                localNodePart, localElemPart, localSidePart,
                drag, roughness)
    # Phase 3: all elements
    for ipart, pd in enumerate(partPopulate):
        with cp.cuda.Device(pd['device_id']):
            print(f"HERE device {ipart} populateElements")
            cufunc.populateElements(
                pd['i34Part'], pd['elnodePart'], pd['elsidePart'],
                pd['nodeCoordsPart'], pd['nodeCoords_llPart'], pd['use_ll'],
                pd['windRotateAnglePart'],
                isdel, isidenode, isbnd,
                nond, nlnd, iond, ilnd,
                nodeNeighNodesStart, nodeNeighNodes,
                nodeNeighElemsStart, nodeNeighElems,
                pd['globalNodeIndices'], pd['globalElemIndices'], pd['globalSideIndices'],
                localNodeIndices, localElemIndices, localSideIndices,
                localNodePart, localElemPart, localSidePart,
                drag, roughness)
    # Phase 4: finalize mesh (needs sharedIdx from phases 1 and 2)
    for ipart, pd in enumerate(partPopulate):
        with cp.cuda.Device(pd['device_id']):
            cufunc.populateMesh(
                pd['i34Part'], pd['elnodePart'], pd['elsidePart'],
                pd['nodeCoordsPart'], pd['nodeCoords_llPart'], pd['use_ll'],
                pd['windRotateAnglePart'], pd['sharedIdxNodeBoundary'], pd['sharedIdxSideBoundary'],
                isdel, isidenode, isbnd,
                nond, nlnd, iond, ilnd,
                nodeNeighNodesStart, nodeNeighNodes,
                nodeNeighElemsStart, nodeNeighElems,
                pd['globalNodeIndices'], pd['globalElemIndices'], pd['globalSideIndices'],
                localNodeIndices, localElemIndices, localSideIndices,
                localNodePart, localElemPart, localSidePart,
                drag, roughness,
                bctides.emo, bctides.efa,
                bctides.umo, bctides.ufa,
                bctides.vmo, bctides.vfa)
    # Phase 5: boundary nodes
    for ipart, pd in enumerate(partPopulate):
        with cp.cuda.Device(pd['device_id']):
            cufunc.populateBoundaryNodes(
                nond, iond,
                bctides.emo, bctides.efa,
                bctides.umo, bctides.ufa,
                bctides.vmo, bctides.vfa,
                bctides.ifltype, bctides.iettype)
    

    # initialize the module for bringing memory from the GPU to the CPU
    globalNodeIndices = []
    globalElemIndices = []
    globalSideIndices = []
    for pd in partPopulate:
        globalNodeIndices.append(pd['globalNodeIndices'])
        globalElemIndices.append(pd['globalElemIndices'])
        globalSideIndices.append(pd['globalSideIndices'])
    dtoh = DeviceToHost(numNodes, numElements, numSides, numPartBoundarySides, nvrt, 
                        globalNodeIndices, globalElemIndices, globalSideIndices, 
                        ckw.configCuda)
    


    ntime, nrec = schism_init(namelist.core, namelist.opt)

    initializeOutput(out_dir, numNodes, numElements, numSides, nvrt, nodesSize, elemsSize, sidesSize)


    profiler = LineProfiler() 
    schism_step_wrapped = profiler(schism_step)
    profiler.add_function(solveWaveEquation)
    profiler.add_function(solver.solve)

    profiler.add_function(computeMeteorologicalForcing)
    profiler.add_function(SFlux.getWinds)
    profiler.add_function(SFlux.getCombinedData)
    # profiler.add_function(SFlux.cacheGet)
    # profiler.add_function(SFlux.cachePut)


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

    print(f"Output time: {getOutputTime():.3f} s")
    end_time = time.time()
    total_time = end_time - start_time
    print(f"\nTotal execution time: {total_time:.3f} s")
    print(f"Average time per step: {(total_time/ntime):.5f} s")

    profiler.print_stats()    
    pass


    

    