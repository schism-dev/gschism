import os
import sys

 
from iop.cudamodule import CudaSourceManager 


cuda_manager = None
mod = None
project_root = None

# Global CUDA function declarations
createMesh = None
populateMesh = None
getStructSizes = None
scaleGravity = None
computeTimeRamps = None
computeBedDeformation = None
computeTidalPotential = None
computeBottomDragVegetation = None
computeVegetationAndDragParameters = None
computeElevationBoundaryConditions = None
computeVerticallyIntegratedVelocity = None
computeDensityCubicSplines = None
computeDensityGradient = None
computeBaroclinicTerm = None
computeElemAverageGhat = None
computeWaveEquation = None
solveWaveEquation = None
computeGradients = None
test = None

def initializeCudaKernelManager(_project_root):
    global createMesh, populateMesh, getStructSizes, cuda_manager, project_root, mod

    project_root = _project_root

    source_file = os.path.join(project_root, "kernels", "meshdata.cu")

    # Create the CUDA source manager
    cuda_manager = CudaSourceManager(project_root)

    try:
        # Create the module with preprocessing
        mod, combined_source = cuda_manager.create_module(
            source_file, 
            options=("-std=c++20", '-diag-suppress=64,1055')
        )

        # Optionally save the combined source for debugging
        cuda_manager.save_combined_source(combined_source, "debug_cuda_source.cu")
        
        # Get kernel functions with better error handling
        try:
            sizes = {'global': 1}
            createMesh     = mod.get_function('createMesh',     sizes, kernel_type='global')
            populateMesh   = mod.get_function('populateMesh',   sizes, kernel_type='global')
            getStructSizes = mod.get_function('getStructSizes', sizes, kernel_type='global')

        except Exception as e:
            print(f"Error getting CUDA functions: {e}")
            # Continue with error handling or exit
            sys.exit(1)

    except Exception as e:
        print(f"Failed to compile CUDA code: {e}")
        sys.exit(1)


def initializeCudaFunctions(sizes):
    global mod
    global scaleGravity, computeTimeRamps
    global computeBedDeformation, computeTidalPotential
    global getBoundaryConditions, computeBottomDragVegetation, computeVegetationAndDragParameters
    global computeElevationBoundaryConditions, computeVerticallyIntegratedVelocity, computeDensityCubicSplines
    global computeDensityGradient, computeBaroclinicTerm, computeElemAverageGhat
    global computeWaveEquation, solveWaveEquation, computeGradients, backtrackAllNodes, test

 
    # Get kernel functions with better error handling
    try:
        scaleGravity                        = mod.get_function('scaleGravity',                        sizes, kernel_type='nodes')
        computeTimeRamps                    = mod.get_function('computeTimeRamps',                    sizes, kernel_type='global')
        computeBedDeformation               = mod.get_function('computeBedDeformation',               sizes, kernel_type='global')
        computeTidalPotential               = mod.get_function('computeTidalPotential',               sizes, kernel_type='nodeBnd')
        getBoundaryConditions               = mod.get_function('getBoundaryConditions',               sizes, kernel_type='nodeBnd')
        computeBottomDragVegetation         = mod.get_function('computeBottomDragVegetation',         sizes, kernel_type='nodes')
        computeVegetationAndDragParameters  = mod.get_function('computeVegetationAndDragParameters',  sizes, kernel_type='sides')
        computeElevationBoundaryConditions  = mod.get_function('computeElevationBoundaryConditions',  sizes, kernel_type='sides')
        computeVerticallyIntegratedVelocity = mod.get_function('computeVerticallyIntegratedVelocity', sizes, kernel_type='global')
        computeDensityCubicSplines          = mod.get_function('computeDensityCubicSplines',          sizes, kernel_type='global')
        computeDensityGradient              = mod.get_function('computeDensityGradient',              sizes, kernel_type='global')
        computeBaroclinicTerm               = mod.get_function('computeBaroclinicTerm',               sizes, kernel_type='global')
        computeElemAverageGhat              = mod.get_function('computeElemAverageGhat',              sizes, kernel_type='global')
        computeWaveEquation                 = mod.get_function('computeWaveEquation',                 sizes, kernel_type='global')
        solveWaveEquation                   = mod.get_function('solveWaveEquation',                   sizes, kernel_type='global')
        computeGradients                    = mod.get_function('computeGradients',                    sizes, kernel_type='global')
        backtrackAllNodes                   = mod.get_function('backtrackAllNodes',                   sizes, kernel_type='sides')

        test = mod.get_function('test', sizes, kernel_type='global')
    except Exception as e:
        print(f"Error getting CUDA functions: {e}")
        # Continue with error handling or exit
        sys.exit(1)
     