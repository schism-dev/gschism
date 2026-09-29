import os
import sys


from iop.cudamodule import CudaSourceManager 


cuda_manager = None
mod = None
project_root = None

# Global CUDA function declarations
createMesh = None
populateMesh = None
populateNodes = None
populateNodes2 = None
populateSides = None
populateElements = None
populateBoundaryNodes = None
populateBoundarySides = None
populateBoundarySideLinks = None
getStructSizes = None
initTidalHarmonicConstituents = None
printMatrix = None
returnMatrix = None
returnX = None
returnB = None
scaleGravity = None
setWinds = None
computeTimeRamps = None
computeBedDeformation = None
computeTidalPotential = None
computeBottomDragVegetation = None
computeTurbulence = None
computeVegetationAndDragParameters = None
computeElevationBoundaryConditions = None
computeHorizontalViscosity = None
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
applyMomentumBC = None
ul0UpdateZCoords = None
ul0ComputeDryElements = None
ul0ComputeDryNodes = None
ul0UpdateRewettedNodes = None
ul0UpdateSides = None
ul0UpdateSides2 = None
zeroNodalVelocities = None
computeNodalVerticalVelocities1 = None
computeNodalVerticalVelocities2 = None
computeNodalHorizontalVelocities1 = None
computeNodalHorizontalVelocities2 = None
resetVelocityOpenBndNodes = None
setVelocityOpenBndNodes = None
computeVerticalVelocity = None
computeMatrixProductAndDotProducts = None
updateSolutionAndDotProducts = None
updateSearchDirection = None
updateBedDeformation = None
shapiroFilter1 = None
shapiroFilter2 = None


getSideBndVar = None
getElementVarAddr = None
getSideVarAddr = None
getNodeVarAddr = None
getElementVarFunc = None
getNodeVarFunc = None
getSideVarFunc = None
getNodeVar = None
getSideVar = None
getElementVar = None

setupCGSolver = None
# setupMultiGPUCGSolver = None XXX REMOVE?
solverCopyResult = None
gatherBFromMesh = None
scatterXToMesh = None
# test = None


def initializeCudaKernelManager(_project_root):
    global createMesh, getStructSizes
    global initTidalHarmonicConstituents, cuda_manager, project_root, mod


    project_root = _project_root

    source_file = os.path.join(project_root, "kernels", "meshdata.cu")

    # Create the CUDA source manager
    cuda_manager = CudaSourceManager(project_root)

    try:
        # Create the module with preprocessing
        mod, combined_source = cuda_manager.create_module(
            source_file, 
            options=("-std=c++20", '-diag-suppress=64,1055', '-D_FORCE_INLINES')
        )

        # Optionally save the combined source for debugging
        cuda_manager.save_combined_source(combined_source, "debug_cuda_source.cu")
        
        # Get kernel functions with better error handling
        try:
            sizes = {'global': 1}
            createMesh     = mod.get_function('createMesh',     sizes, kernel_type='global,onlycurrent')
            getStructSizes = mod.get_function('getStructSizes', sizes, kernel_type='global,onlycurrent')
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
    global populateMesh, populateNodes, populateNodes2, populateSides, populateElements, populateBoundaryNodes, populateBoundarySides, populateBoundarySideLinks
    global scaleGravity, computeTimeRamps, setWinds
    global computeBedDeformation, computeTidalPotential
    global getBoundaryConditions, computeFlowBcVelocities, resetBoundaryAreas, computeBoundaryAreas
    global computeBottomDragVegetation, computeVegetationAndDragParameters, computeTurbulence
    global computeElevationBoundaryConditions, computeVerticallyIntegratedVelocity, computeHorizontalViscosity, computeDensityCubicSplines
    global computeDensityGradient, computeBaroclinicTerm, computeElemAverageGhat, setEssentialBC, updateElevation
    global setupCoefficientMatrixAndRHS, rollUpElevation, computeHorGradientAtSide, solveMomentumEquations, applyMomentumBC, backtrackAllNodes, test
    global setupCGSolver, solverCopyResult, gatherBFromMesh, scatterXToMesh, ul0ComputeDryElements, ul0ComputeDryNodes
    global updateBedDeformation, shapiroFilter1, shapiroFilter2
    global computeMatrixProductAndDotProducts, updateSolutionAndDotProducts, updateSearchDirection
    global initializeScalars, initializePreconditionedResidual, printMatrix, returnMatrix, returnX, returnB
    global ul0UpdateZCoords, ul0ComputeDryElements, ul0ComputeDryNodes, ul0UpdateRewettedNodes, ul0UpdateSides, ul0UpdateSides2
    global zeroNodalVelocities, computeNodalVerticalVelocities1, computeNodalVerticalVelocities2, computeNodalHorizontalVelocities1, computeNodalHorizontalVelocities2
    global resetVelocityOpenBndNodes, setVelocityOpenBndNodes, computeVerticalVelocity
    global getSideBndVar
    global getElementVarAddr, getSideVarAddr, getNodeVarAddr
    global getElementVarFunc, getNodeVarFunc, getSideVarFunc
    global getNodeVar, getSideVar, getElementVar

    # Get kernel functions with better error handling
    try:
        populateMesh                        = mod.get_function('populateMesh',                        sizes, kernel_type='global,onlycurrent')
        populateNodes                       = mod.get_function('populateNodes',                       sizes, kernel_type='nodes,onlycurrent')
        populateNodes2                      = mod.get_function('populateNodes2',                      sizes, kernel_type='nodes,onlycurrent')
        populateSides                       = mod.get_function('populateSides',                       sizes, kernel_type='sides,onlycurrent')
        populateElements                    = mod.get_function('populateElements',                    sizes, kernel_type='elements,onlycurrent')
        populateBoundaryNodes               = mod.get_function('populateBoundaryNodes',               sizes, kernel_type='nodeBnd,onlycurrent')
        populateBoundarySides               = mod.get_function('populateBoundarySides',               sizes, kernel_type='sideBnd,onlycurrent')
        populateBoundarySideLinks           = mod.get_function('populateBoundarySideLinks',           sizes, kernel_type='sideBnd,onlycurrent')

        scaleGravity                        = mod.get_function('scaleGravity',                        sizes, kernel_type='nodes')
        setWinds                            = mod.get_function('setWinds',                            sizes, kernel_type='nodes')
        computeTimeRamps                    = mod.get_function('computeTimeRamps',                    sizes, kernel_type='global')
        computeBedDeformation               = mod.get_function('computeBedDeformation',               sizes, kernel_type='global')
        computeTidalPotential               = mod.get_function('computeTidalPotential',               sizes, kernel_type='nodes')
        getBoundaryConditions               = mod.get_function('getBoundaryConditions',               sizes, kernel_type='nodeBnd')
        computeFlowBcVelocities             = mod.get_function('computeFlowBcVelocities',             sizes, kernel_type='sideBnd')
        resetBoundaryAreas                  = mod.get_function('resetBoundaryAreas',                  sizes, kernel_type='sideBnd')
        computeBoundaryAreas                = mod.get_function('computeBoundaryAreas',                sizes, kernel_type='sideBnd')
        computeBottomDragVegetation         = mod.get_function('computeBottomDragVegetation',         sizes, kernel_type='nodes')
        computeTurbulence                   = mod.get_function('computeTurbulence',                   sizes, kernel_type='nodes')
        computeVegetationAndDragParameters  = mod.get_function('computeVegetationAndDragParameters',  sizes, kernel_type='sides')
        computeElevationBoundaryConditions  = mod.get_function('computeElevationBoundaryConditions',  sizes, kernel_type='nodeBnd')
        computeHorizontalViscosity          = mod.get_function('computeHorizontalViscosity',          sizes, kernel_type='sides')
        computeVerticallyIntegratedVelocity = mod.get_function('computeVerticallyIntegratedVelocity', sizes, kernel_type='sides')
        computeDensityCubicSplines          = mod.get_function('computeDensityCubicSplines',          sizes, kernel_type='elements')
        computeDensityGradient              = mod.get_function('computeDensityGradient',              sizes, kernel_type='elements')
        computeBaroclinicTerm               = mod.get_function('computeBaroclinicTerm',               sizes, kernel_type='sides')
        computeElemAverageGhat              = mod.get_function('computeElemAverageGhat',              sizes, kernel_type='elements')
        setEssentialBC                      = mod.get_function('setEssentialBC',                      sizes, kernel_type='nodes')
        updateElevation                     = mod.get_function('updateElevation',                     sizes, kernel_type='nodes')
        setupCoefficientMatrixAndRHS        = mod.get_function('setupCoefficientMatrixAndRHS',        sizes, kernel_type='nodes')
        rollUpElevation                     = mod.get_function('rollUpElevation',                     sizes, kernel_type='nodes')

        # computeMatrixProductAndDotProducts  = mod.get_function('computeMatrixProductAndDotProducts',  sizes, kernel_type='nodes')
        # updateSolutionAndDotProducts        = mod.get_function('updateSolutionAndDotProducts',        sizes, kernel_type='nodes')
        # updateSearchDirection               = mod.get_function('updateSearchDirection',               sizes, kernel_type='nodes')
        # initializeScalars                   = mod.get_function('initializeScalars',                   sizes, kernel_type='global')
        # initializePreconditionedResidual    = mod.get_function('initializePreconditionedResidual',    sizes, kernel_type='nodes')

        computeHorGradientAtSide            = mod.get_function('computeHorGradientAtSide',            sizes, kernel_type='sides')
        solveMomentumEquations              = mod.get_function('solveMomentumEquations',              sizes, kernel_type='sides')
        applyMomentumBC                     = mod.get_function('applyMomentumBC',                     sizes, kernel_type='sideBnd')
        # updateBedDeformation                = mod.get_function('updateBedDeformation',                sizes, kernel_type='nodes')
        shapiroFilter1                      = mod.get_function('shapiroFilter1',                      sizes, kernel_type='sides')
        shapiroFilter2                      = mod.get_function('shapiroFilter2',                      sizes, kernel_type='sides')
        ul0UpdateZCoords                    = mod.get_function('ul0UpdateZCoords',                    sizes, kernel_type='nodes')
        ul0ComputeDryElements               = mod.get_function('ul0ComputeDryElements',               sizes, kernel_type='elements')
        ul0ComputeDryNodes                  = mod.get_function('ul0ComputeDryNodes',                  sizes, kernel_type='nodes')
        ul0UpdateRewettedNodes              = mod.get_function('ul0UpdateRewettedNodes',              sizes, kernel_type='nodes')
        ul0UpdateSides                      = mod.get_function('ul0UpdateSides',                      sizes, kernel_type='sides')
        ul0UpdateSides2                     = mod.get_function('ul0UpdateSides2',                     sizes, kernel_type='sides')
        zeroNodalVelocities                 = mod.get_function('zeroNodalVelocities',                 sizes, kernel_type='nodes')
        computeNodalVerticalVelocities1     = mod.get_function('computeNodalVerticalVelocities1',     sizes, kernel_type='elements')
        computeNodalVerticalVelocities2     = mod.get_function('computeNodalVerticalVelocities2',     sizes, kernel_type='nodes')
        computeNodalHorizontalVelocities1   = mod.get_function('computeNodalHorizontalVelocities1',   sizes, kernel_type='elements')
        computeNodalHorizontalVelocities2   = mod.get_function('computeNodalHorizontalVelocities2',   sizes, kernel_type='nodes')
        resetVelocityOpenBndNodes           = mod.get_function('resetVelocityOpenBndNodes',           sizes, kernel_type='nodeBnd')
        setVelocityOpenBndNodes             = mod.get_function('setVelocityOpenBndNodes',             sizes, kernel_type='sideBnd')
        computeVerticalVelocity             = mod.get_function('computeVerticalVelocity',             sizes, kernel_type='elements')
        backtrackAllNodes                   = mod.get_function('backtrackAllNodes',                   sizes, kernel_type='sides')

        printMatrix                         = mod.get_function('printMatrix',                         sizes, kernel_type='global')
        returnMatrix                        = mod.get_function('returnMatrix',                        sizes, kernel_type='global')
        returnX                             = mod.get_function('returnX',                             sizes, kernel_type='global')
        returnB                             = mod.get_function('returnB',                             sizes, kernel_type='global')

        getSideBndVar                       = mod.get_function('getSideBndVar',                       sizes, kernel_type='sideBnd')
        getElementVarAddr                   = mod.get_function('getElementVarAddr',                   sizes, kernel_type='global,onlycurrent')
        getSideVarAddr                      = mod.get_function('getSideVarAddr',                      sizes, kernel_type='global,onlycurrent')
        getNodeVarAddr                      = mod.get_function('getNodeVarAddr',                      sizes, kernel_type='global,onlycurrent')
        getElementVarFunc                   = mod.get_function('getElementVarFunc',                   sizes, kernel_type='elements')
        getNodeVarFunc                      = mod.get_function('getNodeVarFunc',                      sizes, kernel_type='nodes')
        getSideVarFunc                      = mod.get_function('getSideVarFunc',                      sizes, kernel_type='sides')
        getNodeVar                          = mod.get_function('getNodeVar',                          sizes, kernel_type='nodes,onlycurrent')
        getSideVar                          = mod.get_function('getSideVar',                          sizes, kernel_type='sides,onlycurrent')
        getElementVar                       = mod.get_function('getElementVar',                       sizes, kernel_type='elements,onlycurrent')

        # test = mod.get_function('test', sizes, kernel_type='global')

        setupCGSolver            = mod.get_function('setupCGSolver', sizes, kernel_type='nodes')
        # setupMultiGPUCGSolver    = mod.get_function('setupMultiGPUCGSolver', sizes, kernel_type='nodes') XXX REMOVE?
        solverCopyResult         = mod.get_function('solverCopyResult', sizes, kernel_type='nodes')
        gatherBFromMesh          = mod.get_function('gatherBFromMesh', sizes, kernel_type='nodes,onlycurrent')
        scatterXToMesh           = mod.get_function('scatterXToMesh', sizes, kernel_type='nodes,onlycurrent')

    except Exception as e:
        print(f"Error getting CUDA functions: {e}")
        # Continue with error handling or exit
        sys.exit(1)
     