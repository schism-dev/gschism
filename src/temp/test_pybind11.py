import os
import glob
import importlib.util
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import pybind11
import numpy as np
import timeit
import sys

sys.path.append(r'/blue/olabarrieta/jo.gonzalez/src/gschism/src')

# import src.io.readSflux as readSflux
from src.io.readSflux import SFlux

def launch_kernel(kernel_func, num_items, n_gpus, threads_per_block=256):
    """
    Launch a kernel across the GPUs.
    
    Args:
        kernel_func: CuPy kernel function to launch
        num_items: Number of items to process (e.g., numPartSides)
        threads_per_block: Number of threads per block (default 256)
    """
    num_blocks = (num_items + threads_per_block - 1) // threads_per_block
    
    for device_id in range(n_gpus):
        with cp.cuda.Device(device_id):
            kernel_func((num_blocks,), (threads_per_block,), ())


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

        # Define the extension
        ext_modules = [
            Extension(
                module_name,
                sources=[cpp_file],
                include_dirs=[pybind11.get_include(), "/apps/intel/2020.0.166/metis/5.2.1/include"],
                language="c++",
                # extra_compile_args=["-std=c++20", "-march=native"],
                extra_link_args=["-lmetis"],
                libraries=['metis'],  # Add METIS to libraries
                library_dirs=['/usr/lib/x86_64-linux-gnu', "/apps/intel/2020.0.166/metis/5.2.1/lib"],  # Add library directory
                runtime_library_dirs=['/usr/lib/x86_64-linux-gnu', "/apps/intel/2020.0.166/metis/5.2.1/lib"],
            )
        ]

        # Dynamically build the extension
        class BuildExt(build_ext):
            def build_extensions(self):
                for ext in self.extensions:
                    ext.extra_compile_args = ["-std=c++17", "-g", "-march=native", "-DDEBUG"]
                super().build_extensions()

        setup(
            name=module_name,
            ext_modules=ext_modules,
            cmdclass={"build_ext": BuildExt},
            script_args=["build_ext", "--inplace"],
        )

    # Find the compiled shared object file
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
    cpp_file = "../src/src/grid/read_and_process_grid.cpp"  # External C++ source file
    module = compile_and_import("test", cpp_file)



    # Test the module
    import cupy as cp

    # mesh_file = "/media/jmg/New450/src/gschism/test/hgrid2.gr3"
    # mesh_file = "/home/jo.gonzalez/src/gschism/test/hgrid.gr3"
    # mesh_file = "/media/jmg/Elements2/src/tmp/schism_verification_tests/Test_CORIE_LSC2/hgrid.gr3"
    # mesh_file = "/home/jo.gonzalez/src/gschism/test/schism_verification_tests/Test_CORIE_LSC2/hgrid.gr3"
    mesh_file = "/blue/olabarrieta/jo.gonzalez/projects/src/test/schism_verification_tests/Test_CORIE_LSC2/hgrid.gr3"

    
    numPartitions = 2
    i34, partElements, partNodes, partSides, elnode, elside, nodes, ic3, isdel = module.readMesh(mesh_file, numPartitions)


    # Before sending the mesh to the GPU, we need to do some preprocessing, including computing the 
    # interpolation weights.

    # Test the readSflux module
    # sflux_inputs = readSflux.readSFluxNamelist(r"/media/jmg/Elements2/src/tmp/schism_verification_tests/Test_CORIE_LSC2/sflux/sflux_inputs.txt")
    nodes = nodes.T # XXX Why .T?
    nodes[:2,:] = np.deg2rad(nodes[:2,:]) # XXX Why deg2rad?
    sflux = SFlux(r'/blue/olabarrieta/jo.gonzalez/projects/src/test/schism_verification_tests/Test_CORIE_LSC2/sflux/sflux_inputs.txt', nodes) 
    print(sflux)

    u_air_node = np.zeros(nodes.shape[1])
    v_air_node = np.zeros(nodes.shape[1])
    p_air_node = np.zeros(nodes.shape[1])
    t_air_node = np.zeros(nodes.shape[1])
    q_air_node = np.zeros(nodes.shape[1])
    sflux.getWind(2452424.0, u_air_node, v_air_node, p_air_node, t_air_node, q_air_node)


    # convert arrays to cupy
    i34             = cp.array(i34,          dtype=cp.int32)
    partElements    = cp.array(partElements, dtype=cp.int32)
    partNodes       = cp.array(partNodes,    dtype=cp.int32)
    partSides       = cp.array(partSides,    dtype=cp.int32)
    elnode          = cp.array(elnode,       dtype=cp.int32)
    nodes           = cp.array(nodes,        dtype=cp.float64)
    ic3             = cp.array(ic3,          dtype=cp.int32)
    isdel           = cp.array(isdel,        dtype=cp.int32)
    # js         = cp.array(js)
    # Read the CUDA source
    with open('../src/kernels/meshdata.cu', 'r') as f:
        cuda_source = f.read()

    # Create the kernel module
    mod = cp.RawModule(code=cuda_source, options=("-std=c++20", '-diag-suppress=64,1055')) # the diag-suppress is to suppress the warning about the anonymous union

    # Get the kernel functions
    createMesh          = mod.get_function('createMesh')
    getStructSizes      = mod.get_function('getStructSizes')
    computeFluxAdvHface = mod.get_function('computeFluxAdvHface')
    # destroy_mesh = mod.get_function('destroyMeshCPU')

    # Get number of GPUs
    n_gpus = cp.cuda.runtime.getDeviceCount()
    print(f"Found {n_gpus} GPUs")

    numNodes = nodes.shape[0]
    numElements = elnode.shape[1]
    numSides = isdel.shape[1]   




    # file_name = readSflux.getFileName('sflux_rad_1', 1)
    # times = readSflux.getFileTimes(file_name, 1000)
    # print(times)

    # times, files = readSflux.getTimesEtc('sflux_rad_1', 1000, 100)
    # print(times)
    # print(files)

    # info = 
    # print(info)




   
    # Get the struct sizes
    sizes = cp.zeros(3, dtype=cp.int32)
    getStructSizes((1,), (1,), sizes)
    [nodesSize, elemsSize, sidesSize] = sizes.tolist()

    # Compute the maximum number of nodes, elements, and sides per partition
    numPartElements = []
    numPartNodes = []
    numPartSides = []
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
       
    nodesMemList = []
    elemsMemList = []
    sidesMemList = []
    for ipart, device_id in enumerate(range(n_gpus)):
        with cp.cuda.Device(device_id):

            # Create the local device arrays.
            i34_part    = cp.array(i34[partElements == ipart])
            elnode_part = cp.array(elnode[:, partElements == ipart])

            # Reserves CPU memory for the structs
            nodesMem = cp.zeros(max(numPartNodes)   *nodesSize, dtype=cp.int8)
            elemsMem = cp.zeros(max(numPartElements)*elemsSize, dtype=cp.int8)
            sidesMem = cp.zeros(max(numPartSides)   *sidesSize,    dtype=cp.int8)

            nodesMemList.append(nodesMem.data.ptr)
            elemsMemList.append(elemsMem.data.ptr)
            sidesMemList.append(sidesMem.data.ptr)

    # Create the arrays of pointers to the start of the device memory
    nodesPtr = cp.zeros(numPartitions, dtype=cp.uintp)
    elemsPtr = cp.zeros(numPartitions, dtype=cp.uintp)
    sidesPtr = cp.zeros(numPartitions, dtype=cp.uintp)
    for ipart, device_id in enumerate(range(n_gpus)):
        nodesPtr[ipart] = nodesMemList[ipart]
        elemsPtr[ipart] = elemsMemList[ipart]
        sidesPtr[ipart] = sidesMemList[ipart]

    # Create the mesh objects, one per GPU
    for ipart, device_id in enumerate(range(n_gpus)):
        with cp.cuda.Device(device_id):

            createMesh((1,), (1,),  
                       (ipart, numNodes, numElements, numSides, 1,
                        max(numPartNodes), max(numPartElements), max(numPartSides),
                        numPartNodes[ipart], numPartElements[ipart], numPartSides[ipart],
                        nodesPtr, elemsPtr, sidesPtr,
                        nodesMemList[ipart], elemsMemList[ipart], sidesMemList[ipart],
                        # , partNodes, partSides,
                        i34[partElements], elnode[:, partElements == ipart], elside[:, partElements == ipart]) #, nodes, ic3, isdel)
            )


    # Launch kernel on each GPU
    launch_kernel(computeFluxAdvHface, max(numPartSides), n_gpus, 256)
    

    # launch_kernel(checkMassConservation, sumtr, max(numPartSides), n_gpus, 256)
    # sumtr.sum()




    time = timeit.timeit('module.readMesh("/media/jmg/New450/src/gschism/test/hgrid.gr3")',  globals={'module': module, 'array': array}, number=1)
    print(f"Time pybind11: {time:.6f} seconds")
    # time2 = timeit.timeit("result = np.sinh(abs(np.sinh(1/pow(array,  2.1))))",  globals={'array': array, 'np': np}, number=10)
    # print(f"Time pybind11: {time:.6f} seconds")
    # print(f"Time numpy: {time2:.6f} seconds")
    # print("Original:", array)
    # print("Doubled:", result)
    # print("numpy:", np.sinh(abs(np.sinh(1/pow(array,  2.1)))))