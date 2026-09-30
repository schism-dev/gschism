import python.cuda_functions as cufunc
import iop.readSflux as readSflux
import iop.writeOutput as writeOutput
import cupy as cp
import numpy as np
import debug.monitor as monitor
from cupyx.scipy.sparse.linalg import cg, gmres, minres, LinearOperator

from python.pcgsolver import PCGSolver
from python.schism_consts import VISTA
# from cupy import cudaMemcpyFromSymbol




import time

_starts = {}
_times = {}
def tic(name):
    cp.cuda.Device().synchronize()
    _starts[name] = time.perf_counter()
def toc(name):
    cp.cuda.Device().synchronize()
    _times[name] = _times.get(name, 0.0) + (time.perf_counter() - _starts[name])



def solveWaveEquation(solver, maxIter=1000, tol=1e-12, ): 
    """
    Solve wave equation using Conjugate Gradient method.
    Assumes all arrays are already stored in Node as doubles.

    NOTATION: rMr and similar variable names stand for the dot product (r^T)Mr.
    """
    global N, NN


    time1 = time.time()
  
    cp.cuda.profiler.start()
    its = solver.solve(rtol=tol, maxiter=maxIter)
    solver.checkSolution() # XXXXXXXX
    cp.cuda.profiler.stop()
    time2 = time.time()
    print(f"Time to setup and solve: {time2 - time1:.3f} , iterations: {its}")

    return its

    



def computeMeteorologicalForcing(opt, time, sflux):
    # Convert simulation time (seconds from start) to Julian day using namelist start and sflux's JD convention
    if opt.nws == 0:
        return

    jd_midnight = sflux.dataset1.julianDay(opt.start_year, opt.start_month, opt.start_day)
    julian_day = jd_midnight + opt.start_hour / 24.0 + time / 86400.0 + opt.utc_start/24.0

    u_air_node, v_air_node, p_air_node, t_air_node, q_air_node = sflux.getWinds(julian_day)

    cufunc.setWinds(u_air_node, v_air_node, p_air_node)

    print("WARNING: Sflux reading is not well optimized yet. I believe it does unnecesary computations every time step")  # XXXXXXXX





def schism_step(it: int, namelist, bctides, sflux, solver):

    core = namelist.core
    opt  = namelist.opt

    

    # SAL option: scale gravity
    cufunc.scaleGravity(opt.loadtide_coef)

    dt = core.dt

    time = it*dt

    # Broadcast to global module (XXX I think this is not needed)
    # time_stamp=time; it_main=it

    
    # Compute time ramps for boundary elevation forcing, wind and pressure forcing and tidal potential forcing
    cufunc.computeTimeRamps(time, opt.nws, opt.if_source, opt.drampbc, opt.drampwind, opt.drampwafo, opt.dramp_ss, opt.dramp)
    
    print(f"time = {time}, it = {it}")


    #  Compute new bed deformation
    cufunc.computeBedDeformation()


    # Earth tidal potential and loading tide at nodes: pre-compute to save time
    cufunc.computeTidalPotential(time, bctides.tfreq, bctides.tamp, bctides.tnf, 
                                 bctides.tear, bctides.jspc)


    tic("computeMeteorologicalForcing")
    computeMeteorologicalForcing(opt, time, sflux)  # Lines 458-1014
    toc("computeMeteorologicalForcing")


    # computeTracerStuff()            # Lines 1015-1145

    
    # # Compute hydraulic transfer blocks together with reading in flux values
    # computeHydraulicTransfer()      # Lines 1146-1330


    # Get boundary conditions from external files
    cufunc.getBoundaryConditions(time)          # Lines 1331-1584


    # # Read in volume/mass sources/sinks
    # readInSourcesSinks()              # Lines 1585-1610


    # # Calculation of cross-section areas and length for flow b.c.
    # computeFlowBcAreas()                   # Lines 1611-1758

    # Computes the cross-section areas and lengths for flow boundary conditions at each segment.
    cufunc.resetBoundaryAreas()
    cufunc.computeBoundaryAreas()


    # XXX I think I will combine all BCs into one function
    # Compute new vel. for flow b.c. (uth,vth)
    cufunc.computeFlowBcVelocities(time)            # Lines 1759-1800



    # Bottom drag coefficients for nchi=-1 or 1; Cd and Cdp for nchi=0 already read in
    ####cufunc.computeBottomDrag()                  # Lines 1801-1961


    # bypass solver for transport only XXXX
    # I haven't included this in the code yet


    # bottom drag and vegetation drag
    cufunc.computeBottomDragVegetation()      # Lines 2353-2662



    cufunc.computeTurbulence()                # Lines 2663-3839


    # Horizontal viscosity, implemented as a filter (also Shapiro filter)
    cufunc.computeHorizontalViscosity()       # Lines 3840-4160


    # Backtracking/upwind for momentum
    tic("backtrackAllNodes")
    try:
        monitor.debugPlotIdx[0] = 0   
             
        cufunc.backtrackAllNodes()              # Lines 4161-5029

        # if it >= 0 and not VISTA:
        #     import debug.newmonitor as newmonitor
        #     newmonitor.show()

        # if it >= 0:
        #     monitor.plotMesh()        

        cp.cuda.Device().synchronize() # XXXX needed? There are other synchronize calls inside backtrackAllNodes
        # if it >= 69:
        #     res = monitor.plotMesh()
        #     if not res:
        #         raise ValueError("Mesh is not valid")
    except Exception as e:
        print(f"Error in backtrackAllNodes: {e}")
    toc("backtrackAllNodes")
    print(_times)

    # Vertical advection 
    if opt.iupwind_mom != 0:
        raise ValueError("iupwind_mom != 0 is not implemented yet")
    # computeVerticalAdvection()        # Lines 5030-5122

    # if it >= 0:
    #     monitor.fortranDebugPlot(it)

    cufunc.computeElevationBoundaryConditions(time)


    cufunc.computeVerticallyIntegratedVelocity()       # Lines 5122-5320


    # compute vegetation and drag model parameters
    cufunc.computeVegetationAndDragParameters()


    # Density gradient at element centers
    cufunc.computeDensityGradient()      


    cufunc.computeBaroclinicTerm()


    # Elem. average ghat1 (in eframe if ics=2; \breve{G} in notes). Dimension: m^2/s
    cufunc.computeElemAverageGhat()                 # Lines 5533-5783


    # setup coefficient matrix, sparsem, for the wave equation
    #!...  No elevation essential b.c. are imposed yet but other b.c. is imposed
    cufunc.setupCoefficientMatrixAndRHS(time)       # Lines 5826-


    cufunc.rollUpElevation()


    cufunc.setEssentialBC()


    #  Solve the wave equation
    # cufunc.printMatrix(True)   
    solveWaveEquation(solver, maxIter = opt.mxitn0, tol = opt.rtol0)


    # cufunc.printMatrix(True)   
    # if it >= 2:
    #     N = 592
    #     A = cp.zeros((N, N))
    #     cufunc.returnMatrix(A)
    #     # Create array snz with the nonzero, non-diagonal elements of res
    #     # Extract the diagonal
    #     diagonal = cp.diag(A)
    #     # Create array with zeros on diagonal
    #     A_no_diag = A.copy()
    #     cp.fill_diagonal(A_no_diag, 0)

    #     # Count non-zeros per row (excluding diagonal)
    #     nnz_per_row = cp.sum(A_no_diag != 0, axis=1)

    #     # Extract non-zero, non-diagonal elements
    #     snz = A_no_diag[A_no_diag != 0]      

    #     X = cp.zeros(N)
    #     cufunc.returnX(X)

    #     B = cp.zeros(N)
    #     cufunc.returnB(B)    

    cufunc.updateElevation()


    # Momentum equations	
     
    # Precompute elevation gradient, atmo. pressure and earth tidal potential 
    cufunc.computeHorGradientAtSide()  # lines -7296


    cufunc.solveMomentumEquations(time)
    

    cufunc.applyMomentumBC()


    # Transport    lines 7388-8111

    if opt.ishapiro == 1:
        # The Shapiro filter is applied in two steps
        cufunc.shapiroFilter1()
        cufunc.shapiroFilter2()


    # cufunc.updateBedDeformation()  // XXXX CHECK IF updateZCoordData MUST BE CALLED HERE

    cufunc.computeVerticalVelocity()

    cufunc.ul0UpdateZCoords()
    
    cufunc.ul0ComputeDryElements()

    cufunc.ul0ComputeDryNodes()

    cufunc.ul0UpdateRewettedNodes(it, 0)  # iths XXXXXXX
    
    cufunc.ul0UpdateSides(it, 0)  # iths XXXXXXXX

    cufunc.ul0UpdateSides2()  # iths XXXXXXXX

    # Adjust mass after level change
    # XXXXXXXXX IMPLEMENT THIS


    # Compute nodal velocities from side velocities
    cufunc.zeroNodalVelocities()
    cufunc.computeNodalVerticalVelocities1()
    cufunc.computeNodalVerticalVelocities2()
    cufunc.computeNodalHorizontalVelocities1()
    cufunc.computeNodalHorizontalVelocities2()
    pass
    cufunc.resetVelocityOpenBndNodes()
    cufunc.setVelocityOpenBndNodes()



    # var = monitor.getElementVar("z", 0)
    # syncrhonize 


    # var = monitor.getNodeVar("eta")
    # eta[:, it] = var
    # cufunc.getElementVar(monitor.varCodeDict["x"], var, 0) # monitor.varCodeDict["x"]


    # Update bed deformation and depth info  lines 8156-
 

    # Optional calculations


    




    
