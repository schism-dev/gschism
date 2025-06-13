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
initTidalHarmonicConstituents = None
printMatrix = None
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
setEssentialBC = None
updateElevation = None
setupCoefficientMatrixAndRHS = None
rollUpElevation = None
computeHorGradientAtSide = None
solveMomentumEquations = None
ul0UpdateZCoords = None
ul0ComputeDryElementsAndNodes = None
ul0UpdateRewettedNodes = None
ul0UpdateSides = None
zeroNodalVelocities = None
computeNodalVerticalVelocities1 = None
computeNodalVerticalVelocities2 = None
computeNodalHorizontalVelocities = None
resetVelocityOpenBndNodes = None
setVelocityOpenBndNodes = None
computeVerticalVelocity = None
computeMatrixProductAndDotProducts = None
updateSolutionAndDotProducts = None
updateSearchDirection = None
updateBedDeformation = None
shapiroFilter1 = None
shapiroFilter2 = None

getElementVar = None
getNodeVar = None
getSideVar = None

test = None

def initializeCudaKernelManager(_project_root):
    global createMesh, populateMesh, getStructSizes, initTidalHarmonicConstituents, cuda_manager, project_root, mod

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
            initTidalHarmonicConstituents = mod.get_function('initTidalHarmonicConstituents', sizes, kernel_type='global')

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
    global computeDensityGradient, computeBaroclinicTerm, computeElemAverageGhat, setEssentialBC, updateElevation
    global setupCoefficientMatrixAndRHS, rollUpElevation, computeHorGradientAtSide, solveMomentumEquations, backtrackAllNodes, test
    global updateBedDeformation, shapiroFilter1, shapiroFilter2
    global computeMatrixProductAndDotProducts, updateSolutionAndDotProducts, updateSearchDirection
    global initializePreconditionedResidual, printMatrix
    global ul0UpdateZCoords, ul0ComputeDryElementsAndNodes, ul0UpdateRewettedNodes, ul0UpdateSides
    global zeroNodalVelocities, computeNodalVerticalVelocities1, computeNodalVerticalVelocities2, computeNodalHorizontalVelocities
    global resetVelocityOpenBndNodes, setVelocityOpenBndNodes, computeVerticalVelocity
    global getElementVar, getNodeVar, getSideVar

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
        computeVerticallyIntegratedVelocity = mod.get_function('computeVerticallyIntegratedVelocity', sizes, kernel_type='sides')
        computeDensityCubicSplines          = mod.get_function('computeDensityCubicSplines',          sizes, kernel_type='elements')
        computeDensityGradient              = mod.get_function('computeDensityGradient',              sizes, kernel_type='elements')
        computeBaroclinicTerm               = mod.get_function('computeBaroclinicTerm',               sizes, kernel_type='sides')
        computeElemAverageGhat              = mod.get_function('computeElemAverageGhat',              sizes, kernel_type='elements')
        setEssentialBC                      = mod.get_function('setEssentialBC',                      sizes, kernel_type='nodes')
        updateElevation                     = mod.get_function('updateElevation',                     sizes, kernel_type='nodes')
        setupCoefficientMatrixAndRHS        = mod.get_function('setupCoefficientMatrixAndRHS',        sizes, kernel_type='nodes')
        rollUpElevation                     = mod.get_function('rollUpElevation',                     sizes, kernel_type='nodes')
        computeMatrixProductAndDotProducts  = mod.get_function('computeMatrixProductAndDotProducts',  sizes, kernel_type='nodes')
        updateSolutionAndDotProducts        = mod.get_function('updateSolutionAndDotProducts',        sizes, kernel_type='nodes')
        updateSearchDirection               = mod.get_function('updateSearchDirection',               sizes, kernel_type='nodes')
        initializePreconditionedResidual    = mod.get_function('initializePreconditionedResidual',    sizes, kernel_type='nodes')
        computeHorGradientAtSide            = mod.get_function('computeHorGradientAtSide',            sizes, kernel_type='sides')
        solveMomentumEquations              = mod.get_function('solveMomentumEquations',              sizes, kernel_type='nodes')
        # updateBedDeformation                = mod.get_function('updateBedDeformation',                sizes, kernel_type='nodes')
        shapiroFilter1                      = mod.get_function('shapiroFilter1',                      sizes, kernel_type='nodes')
        shapiroFilter2                      = mod.get_function('shapiroFilter2',                      sizes, kernel_type='nodes')
        ul0UpdateZCoords                    = mod.get_function('ul0UpdateZCoords',                    sizes, kernel_type='nodes')
        ul0ComputeDryElementsAndNodes       = mod.get_function('ul0ComputeDryElementsAndNodes',       sizes, kernel_type='elements')
        ul0UpdateRewettedNodes              = mod.get_function('ul0UpdateRewettedNodes',              sizes, kernel_type='nodes')
        ul0UpdateSides                      = mod.get_function('ul0UpdateSides',                      sizes, kernel_type='sides')
        zeroNodalVelocities                 = mod.get_function('zeroNodalVelocities',                 sizes, kernel_type='nodes')
        computeNodalVerticalVelocities1     = mod.get_function('computeNodalVerticalVelocities1',     sizes, kernel_type='elements')
        computeNodalVerticalVelocities2     = mod.get_function('computeNodalVerticalVelocities2',     sizes, kernel_type='nodes')
        computeNodalHorizontalVelocities    = mod.get_function('computeNodalHorizontalVelocities',    sizes, kernel_type='nodes')
        resetVelocityOpenBndNodes           = mod.get_function('resetVelocityOpenBndNodes',           sizes, kernel_type='nodeBnd')
        setVelocityOpenBndNodes             = mod.get_function('setVelocityOpenBndNodes',             sizes, kernel_type='nodeBnd')
        computeVerticalVelocity             = mod.get_function('computeVerticalVelocity',             sizes, kernel_type='elements')
        backtrackAllNodes                   = mod.get_function('backtrackAllNodes',                   sizes, kernel_type='sides')
        printMatrix                         = mod.get_function('printMatrix',                         sizes, kernel_type='nodes')
        getElementVar                       = mod.get_function('getElementVar',                       sizes, kernel_type='elements')
        getNodeVar                          = mod.get_function('getNodeVar',                          sizes, kernel_type='nodes')   
        getSideVar                          = mod.get_function('getSideVar',                          sizes, kernel_type='sides')   

        test = mod.get_function('test', sizes, kernel_type='global')
    except Exception as e:
        print(f"Error getting CUDA functions: {e}")
        # Continue with error handling or exit
        sys.exit(1)
     