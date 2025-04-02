import os
import sys

 
from iop.cudamodule import CudaSourceManager 


# Global CUDA function declarations
createMesh = None
populateMesh = None
getStructSizes = None
scaleGravity = None
computeTimeRamps = None
computeBedDeformation = None
computeTidalPotential = None
computeBottomDragVegetation = None
computeElevationBoundaryConditions = None
computeDensityCubicSplines = None
computeDensityGradient = None
computeBaroclinicTerm = None
computeElemAverageGhat = None
computeWaveEquation = None
solveWaveEquation = None
computeGradients = None
test = None

def initialize_cuda_functions(project_root):
    global createMesh, populateMesh, getStructSizes, scaleGravity, computeTimeRamps
    global computeBedDeformation, computeTidalPotential, computeBottomDragVegetation
    global computeElevationBoundaryConditions, computeDensityCubicSplines
    global computeDensityGradient, computeBaroclinicTerm, computeElemAverageGhat
    global computeWaveEquation, solveWaveEquation, computeGradients, test

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
            createMesh     = mod.get_function('createMesh')
            populateMesh   = mod.get_function('populateMesh')
            getStructSizes = mod.get_function('getStructSizes')

            
            scaleGravity                       = mod.get_function('scaleGravity')
            computeTimeRamps                   = mod.get_function('computeTimeRamps')
            computeBedDeformation              = mod.get_function('computeBedDeformation')
            computeTidalPotential              = mod.get_function('computeTidalPotential')
            computeBottomDragVegetation        = mod.get_function('computeBottomDragVegetation')
            computeElevationBoundaryConditions = mod.get_function('computeElevationBoundaryConditions')
            computeDensityCubicSplines         = mod.get_function('computeDensityCubicSplines')
            computeDensityGradient             = mod.get_function('computeDensityGradient')
            computeBaroclinicTerm              = mod.get_function('computeBaroclinicTerm')
            computeElemAverageGhat             = mod.get_function('computeElemAverageGhat')
            computeWaveEquation                = mod.get_function('computeWaveEquation')
            solveWaveEquation                  = mod.get_function('solveWaveEquation')
            computeGradients                   = mod.get_function('computeGradients')

            test = mod.get_function('test')
        except Exception as e:
            print(f"Error getting CUDA functions: {e}")
            # Continue with error handling or exit
            sys.exit(1)

    except Exception as e:
        print(f"Failed to compile CUDA code: {e}")
        sys.exit(1)