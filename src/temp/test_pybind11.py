import os
import glob
import importlib.util
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import pybind11
import numpy as np
import timeit
import sys
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection


sys.path.append(r'/home/jo.gonzalez/src/gschism/src')  
# sys.path.append(r'/media/jmg/New450/src/gschism/src')  

# in_dir = '/media/jmg/Elements2/src/tmp/schism_verification_tests/Test_CORIE_LSC2' # XXX
# in_dir = '/home/jo.gonzalez/src/gschism/test/schism_verification_tests/Test_CORIE_LSC2' # XXX
# in_dir = '/blue/olabarrieta/jo.gonzalez/projects/src/test/schism_verification_tests/Test_CORIE_LSC2' # XXX
# in_dir = '/mnt/c/Users/jondina/src/temp/schism_examples/schism_verification_tests/Test_CORIE'
in_dir = '/home/jo.gonzalez/src/gschism_old/test/schism_verification_tests/Test_CORIE'
# in_dir = '/media/jmg/Elements2/src/tmp/schism_verification_tests/Test_CORIE'

# import src.io.readSflux as readSflux
from iop.readSflux import SFlux

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

        # # Define the extension
        # ext_modules = [
        #     Extension(
        #         module_name,
        #         sources=[cpp_file],
        #         # include_dirs=[pybind11.get_include(), "/apps/intel/2020.0.166/metis/5.2.1/include"],
        #         include_dirs=[pybind11.get_include(), "/home/jmg/local/include", "/home/jmg/src/lib/metis/METIS/GKlib"],
        #         language="c++",
        #         # extra_compile_args=["-std=c++20", "-march=native"],
        #         extra_link_args=["-lmetis", '-lGKlib'],
        #         libraries=['metis', 'GKlib'],  # Add METIS to libraries
        #         # library_dirs=['/usr/lib/x86_64-linux-gnu', "/apps/intel/2020.0.166/metis/5.2.1/lib"],  # Add library directory
        #         library_dirs=['/usr/lib/x86_64-linux-gnu', '/home/jmg/local/lib', "/home/jmg/src/lib/metis/METIS/GKlib/build/Linux-x86_64"],  # Add library directory
        #         # runtime_library_dirs=['/usr/lib/x86_64-linux-gnu', "/apps/intel/2020.0.166/metis/5.2.1/lib"],
        #         runtime_library_dirs=['/usr/lib/x86_64-linux-gnu', '/home/jmg/local/lib', "/home/jmg/src/lib/metis/METIS/GKlib/build/Linux-x86_64"],
        #     )
        # ]


        # Define the extension
        ext_modules = [
            Extension(
                module_name,
                sources=[cpp_file],
                include_dirs=[pybind11.get_include(), "/apps/intel/2020.0.166/metis/5.2.1/include"],
                # include_dirs=[pybind11.get_include(), "/home/jmg/local/include", "/home/jmg/src/lib/metis/METIS/GKlib"],
                language="c++",
                # extra_compile_args=["-std=c++20", "-march=native"],
                extra_link_args=["-lmetis"],
                libraries=['metis'],  # Add METIS to libraries
                # library_dirs=['/usr/lib/x86_64-linux-gnu', "/apps/intel/2020.0.166/metis/5.2.1/lib"],  # Add library directory
                library_dirs=['/usr/lib/x86_64-linux-gnu', '/home/jmg/local/lib', "/home/jmg/src/lib/metis/METIS/GKlib/build/Linux-x86_64"],  # Add library directory
                # runtime_library_dirs=['/usr/lib/x86_64-linux-gnu', "/apps/intel/2020.0.166/metis/5.2.1/lib"],
                runtime_library_dirs=['/usr/lib/x86_64-linux-gnu', '/home/jmg/local/lib', "/home/jmg/src/lib/metis/METIS/GKlib/build/Linux-x86_64"],
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
    # Get the project root directory (src folder)
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    # Set paths relative to project root
    cpp_file = os.path.join(project_root, "grid", "read_and_process_grid.cpp")

    print("Looking for file at:", os.path.abspath(cpp_file))
    module = compile_and_import("test", cpp_file)



    # Test the module
    import cupy as cp

    # mesh_file = "/media/jmg/New450/src/gschism/test/hgrid2.gr3"
    # mesh_file = "/home/jo.gonzalez/src/gschism/test/hgrid.gr3"
    # mesh_file = "/media/jmg/Elements2/src/tmp/schism_verification_tests/Test_CORIE_LSC2/hgrid.gr3"
    # mesh_file = "/home/jo.gonzalez/src/gschism/test/schism_verification_tests/Test_CORIE_LSC2/hgrid.gr3"
    # mesh_file = "/blue/olabarrieta/jo.gonzalez/projects/src/test/schism_verification_tests/Test_CORIE_LSC2/hgrid.gr3"
    hgrid_file = os.path.join(in_dir, 'hgrid.ll')
    vgrid_file = os.path.join(in_dir, 'vgrid.in')

    

    numPartitions = 1
    i34, partElements, partNodes, partSides, elnode, elside, nodeCoords, \
        ic3, isdel, isidenode,isbnd, nodeNeighElems, nodeNeighNodes = \
            module.readHGrid(hgrid_file, numPartitions)
    
    print(nodeCoords.shape)
    print(nodeNeighNodes)


    [ivcor, nvrt, kz, nsiglev], ztot, [h_s, theta_b, theta_f], sigma_lcl, sigma, cs, dcs = \
        module.readVGrid(vgrid_file, nodeCoords.shape[1])

    wind_rotate_angle = module.readWindRotationAngle(os.path.join(in_dir, 'windrot_geo2proj.gr3'), elnode.shape[1])


    # Before sending the mesh to the GPU, we need to do some preprocessing, including computing the 
    # interpolation weights.

    # Test the readSflux module
    # sflux_inputs = readSflux.readSFluxNamelist(r"/media/jmg/Elements2/src/tmp/schism_verification_tests/Test_CORIE_LSC2/sflux/sflux_inputs.txt")
    nodeCoords[:2,:] = np.deg2rad(nodeCoords[:2,:]) # XXX Why deg2rad?
    # sflux = SFlux(r'/blue/olabarrieta/jo.gonzalez/projects/src/test/schism_verification_tests/Test_CORIE_LSC2/sflux/sflux_inputs.txt', nodes) 
    sflux = SFlux(in_dir, r'sflux_inputs.txt', nodeCoords) 
    print(sflux)

    u_air_node = np.zeros(nodeCoords.shape[1])
    v_air_node = np.zeros(nodeCoords.shape[1])
    p_air_node = np.zeros(nodeCoords.shape[1])
    t_air_node = np.zeros(nodeCoords.shape[1])
    q_air_node = np.zeros(nodeCoords.shape[1])
    sflux.getWind(2452424.0, u_air_node, v_air_node, p_air_node, t_air_node, q_air_node)


    
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
    nodeNeighNodesStart  = cp.array(nodeNeighNodes[1], dtype=cp.int32, order='F')
    nodeNeighNodes       = cp.array(nodeNeighNodes[0], dtype=cp.int32, order='F')
    nodeNeighElemsStart  = cp.array(nodeNeighElems[1], dtype=cp.int32, order='F')
    nodeNeighElems       = cp.array(nodeNeighElems[0], dtype=cp.int32, order='F')
    nodeCoords           = cp.array(nodeCoords,        dtype=cp.float64, order='F')
    wind_rotate_angle    = cp.array(wind_rotate_angle, dtype=cp.float64, order='F')

    ztot                 = cp.array(ztot,              dtype=cp.float64, order='F')
    sigma_lcl            = cp.array(sigma_lcl,         dtype=cp.float64, order='F')
    sigma                = cp.array(sigma,             dtype=cp.float64, order='F')
    cs                   = cp.array(cs,                dtype=cp.float64, order='F')
    dcs                  = cp.array(dcs,               dtype=cp.float64, order='F')

    if True:   

        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Plot all nodes
        plt.plot(nodeCoords[0,:].get(), nodeCoords[1,:].get(), '.')

        # Batch node labels
        # x = nodeCoords[0,:].get()
        # y = nodeCoords[1,:].get()
        # texts = []

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
        
        
        points = np.stack((np.column_stack((x1, y1)), np.column_stack((x2, y2))), axis=1)
        lc = LineCollection(points, colors='r', alpha=0.3)
        ax.add_collection(lc)
        points = np.stack((np.column_stack((x2, y2)), np.column_stack((x3, y3))), axis=1)
        lc = LineCollection(points, colors='r', alpha=0.3)
        ax.add_collection(lc)
        points = np.stack((np.column_stack((x3, y3)), np.column_stack((x1, y1))), axis=1)
        lc = LineCollection(points, colors='r', alpha=0.3)
        ax.add_collection(lc)

        
        # Function to display tooltip
        def on_hover(event):
            if event.inaxes == ax:
                # Check if the mouse is over a line
                cont, ind = lc.contains(event)
                if cont:
                    # Get the index of the line
                    line_index = ind['ind'][0]
                    # Display the tooltip
                    ax.set_title(f'Line between Node {nodes[0, line_index]} and Node {nodes[1, line_index]}')
                    fig.canvas.draw_idle()

        # for i in range(x.size):
        #     texts.append(plt.text(x[i] + 2, y[i] + 2, str(i), 
        #                         fontsize=8))
            
        # # Batch triangle plotting
        # nodes = elnode[:3,:].get()  # Get all triangle nodes at once
        # for i in range(elnode.shape[1]):
        #     # Get triangle nodes and close the triangle
        #     tri_nodes = np.append(nodes[:,i], nodes[0,i])
            
        #     # Get coordinates
        #     x = nodeCoords[0,tri_nodes].get()
        #     y = nodeCoords[1,tri_nodes].get()
            
        #     # Plot triangle
        #     plt.plot(x, y, 'r-', alpha=0.3)
            
        #     # Triangle label at centroid
        #     centroid_x = np.mean(x[:-1])
        #     centroid_y = np.mean(y[:-1])
        #     # texts.append(plt.text(centroid_x, centroid_y, f'T{i}', 
        #     #                     color='red', fontsize=8))
            
        #     # # Side labels
        #     # sides = elside[:3,i].get()
        #     # for j in range(3):
        #     #     mid_x = (x[j] + x[j+1]) / 2
        #     #     mid_y = (y[j] + y[j+1]) / 2
        #     #     texts.append(plt.text(mid_x, mid_y + 3, f'S{sides[(j + 2) % 3]}',
        #     #                         color='blue', fontsize=8))
        
        # Plot segments and intersections
        segment1 = [[-2.213340, 0.760580], [-2.203340, 0.760580]]
        intersection1 = [[-2.212069, 0.760580], [-2.210023, 0.760580], 
                        [-2.208977, 0.760580], [-2.207597, 0.760580], 
                        [-2.205915, 0.760580], [-2.204664, 0.760580]]
        plt.plot(np.array(segment1)[:,0], np.array(segment1)[:,1], 'g--')
        plt.plot(np.array(intersection1)[:,0], np.array(intersection1)[:,1], 'kx')

        segment2 = [[-2.213340, 0.760580], [-2.203340, 0.770580]]
        intersection2 = [[-2.212657, 0.761263], [-2.212017, 0.761903], 
                        [-2.210810, 0.763110], [-2.209846, 0.764074], 
                        [-2.209348, 0.764572], [-2.207661, 0.766259], 
                        [-2.206123, 0.767797], [-2.205463, 0.768457], 
                        [-2.203933, 0.769987]]
        plt.plot(np.array(segment2)[:,0], np.array(segment2)[:,1], 'g--')
        plt.plot(np.array(intersection2)[:,0], np.array(intersection2)[:,1], 'kx')
            

        # Connect the hover event
        fig.canvas.mpl_connect('motion_notify_event', on_hover)
        
        plt.show()
        plt.pause(0)

    # js         = cp.array(js)
    # Read the CUDA source
    cuda_source = ""
    source_files = [
        os.path.join(project_root, "kernels", "meshdata.h"),
        os.path.join(project_root, "kernels", "meshdata.cu")
    ]
    for file in source_files:
        with open(file, 'r') as f:
            cuda_source += f.read()

    # Create the kernel module
    mod = cp.RawModule(code=cuda_source, 
                       options=("-std=c++20", '-diag-suppress=64,1055')) # the diag-suppress is to suppress the warning about the anonymous union

    
    # Get the kernel functions
    createMesh          = mod.get_function('createMesh')
    populateMesh        = mod.get_function('populateMesh')
    getStructSizes      = mod.get_function('getStructSizes')
    computeFluxAdvHface = mod.get_function('computeFluxAdvHface')
    test                = mod.get_function('test')
    # destroy_mesh = mod.get_function('destroyMeshCPU')

    # Get number of GPUs
    n_gpus = cp.cuda.runtime.getDeviceCount()
    print(f"Found {n_gpus} GPUs")

    # This list stores references to the device memory of the local arrays. This ensures that the memory is not freed until the mesh is destroyed.
    deviceMemRef = [ztot, sigma_lcl, sigma, cs, dcs]

    numNodes = nodeCoords.shape[1]
    numElements = elnode.shape[1]
    numSides = isdel.shape[1]  
 



   
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
            sidesMem = cp.zeros(max(numPartSides)   *sidesSize, dtype=cp.int8)

            nodesMemList.append(nodesMem.data.ptr)
            elemsMemList.append(elemsMem.data.ptr)
            sidesMemList.append(sidesMem.data.ptr)

            # Add the arrays to the device memory reference list. This ensures that the memory is not freed until the mesh is destroyed.
            deviceMemRef.append([i34_part, elnode_part, nodesMem, elemsMem, sidesMem])

    # Create the arrays of pointers to the start of the device memory
    nodesPtr = cp.array(nodesMemList, dtype=cp.uintp)
    elemsPtr = cp.array(elemsMemList, dtype=cp.uintp)
    sidesPtr = cp.array(sidesMemList, dtype=cp.uintp)



    # Create the local node and element neighbor arrays
    for ipart, device_id in enumerate(range(n_gpus)):
        with cp.cuda.Device(device_id):

            maskNodesPart = (partNodes == ipart)
            maskElemsPart = (partElements == ipart)

            partNodeNeighNodesSize = maskNodesPart.sum().item()
            partNodeNeighElemsSize = maskElemsPart.sum().item()

            # compute the local size of the node and element neighbor arrays
            partNodeNeighNodesCount = (nodeNeighNodesStart[1:][maskNodesPart] - nodeNeighNodesStart[:-1][maskNodesPart])
            partNodeNeighElemsCount = (nodeNeighElemsStart[1:][maskNodesPart] - nodeNeighElemsStart[:-1][maskNodesPart])

            pointerSize = cp.dtype(cp.uintp).itemsize 
            partNodeNeighNodesMem = cp.zeros(partNodeNeighNodesCount.sum().item()*pointerSize, dtype=cp.int8)
            partNodeNeighElemsMem = cp.zeros(partNodeNeighElemsCount.sum().item()*pointerSize, dtype=cp.int8)


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
            windRotateAnglePart = wind_rotate_angle[maskElemsPart].copy(order='F')

            # Creates the mesh object on each GPU. This function mostly deals with letting the GPU know
            # which memory is being used by the mesh, both locally and in other GPUs.
            createMesh((1,), (1,),  
                       (ipart, numNodes, numElements, numSides, 
                        numPartNodes[ipart], numPartElements[ipart], numPartSides[ipart], nvrt, kz,
                        max(numPartNodes), max(numPartElements), max(numPartSides),
                        nodesPtr, elemsPtr, sidesPtr,
                        nodesMem, elemsMem, sidesMem,
                        partNodeNeighNodesMem, partNodeNeighElemsMem,
                        h_s.item(), theta_b.item(), theta_f.item()  , sigma, sigma, ztot, cs, dcs)
            )
            
            # Given a local index, this array will return the global index
            globalNodeIndices = cp.arange(numPartNodes   [ipart]).astype(cp.int32)[maskNodesPart]
            globalElemIndices = cp.arange(numPartElements[ipart]).astype(cp.int32)[maskElemsPart]
            globalSideIndices = cp.arange(numPartSides   [ipart]).astype(cp.int32)[maskSidesPart]


    for ipart in range(n_gpus):
        with cp.cuda.Device(ipart):
            # Adds additional information to the mesh object. This function deals with the actual data
            # of the mesh, including the nodes, elements, sides, and their neighbors.
            populateMesh((1,), (1,),  (i34Part, elnodePart, elsidePart, nodeCoordsPart, windRotateAnglePart, 
                                       isdel, isidenode, isbnd,
                                       nodeNeighNodesStart, nodeNeighNodes,
                                       nodeNeighElemsStart, nodeNeighElems,
                                       globalNodeIndices, globalElemIndices, globalSideIndices,
                                       localNodeIndices, localElemIndices, localSideIndices,
                                       localNodePart, localElemPart, localSidePart))

            test((1,), (1,), (-2.21334,  0.76058, -1.0))

    # Launch kernel on each GPU
    launch_kernel(computeFluxAdvHface, max(numPartSides), n_gpus, 256)

    print(len(deviceMemRef))
    

    # launch_kernel(checkMassConservation, sumtr, max(numPartSides), n_gpus, 256)
    # sumtr.sum()




    # time = timeit.timeit('module.readHGrid("/media/jmg/New450/src/gschism/test/hgrid.gr3")',  globals={'module': module, 'array': array}, number=1)
    # print(f"Time pybind11: {time:.6f} seconds")
    # time2 = timeit.timeit("result = np.sinh(abs(np.sinh(1/pow(array,  2.1))))",  globals={'array': array, 'np': np}, number=10)
    # print(f"Time pybind11: {time:.6f} seconds")
    # print(f"Time numpy: {time2:.6f} seconds")
    # print("Original:", array)
    # print("Doubled:", result)
    # print("numpy:", np.sinh(abs(np.sinh(1/pow(array,  2.1)))))