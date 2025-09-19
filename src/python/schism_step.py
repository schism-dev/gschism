import python.cuda_functions as cufunc
import iop.readSflux as readSflux
import iop.writeOutput as writeOutput
import cupy as cp
import numpy as np
import debug.monitor as monitor
from cupyx.scipy.sparse.linalg import cg, gmres, minres, LinearOperator

from python.pcgsolver import PCGSolver
# from cupy import cudaMemcpyFromSymbol

solver = None


import time


def solveWaveEquation(maxIter=1000, tol=1e-12): 
    """
    Solve wave equation using Conjugate Gradient method.
    Assumes all arrays are already stored in Node as doubles.

    NOTATION: rMr and similar variable names stand for the dot product (r^T)Mr.
    """
    global N, NN

    global solver

    time1 = time.time()
  
    if solver is None:
        solver = PCGSolver()
        

    its = solver.solve(rtol=tol, maxiter=maxIter)

    time2 = time.time()
    print(f"Time to setup and solve: {time2 - time1:.3f} , iterations: {its}")
    return

    time3 = time.time()
    iterCount = 0
    rMr0      = None

    # Checks that there's only one GPU device, otherwise it will fail
    if cp.cuda.runtime.getDeviceCount() > 1:
        raise ValueError("This function is not YET parallelized for multiple GPUs. Reductions need to be extended to all devices.")


    scalars = cp.zeros(3, dtype=cp.float64) # rMr_old, rMr_new, pAp
    cufunc.initializeScalars(scalars)


    # Compute initial preconditioned residual Mr
    cufunc.initializePreconditionedResidual()
    cp.cuda.Device().synchronize()

    rMr_new = scalars[1].item()
    if (rMr_new <= 1e-30):
        # If the initial residual is already small enough, return 0
        return 0
    
    rMr_0 = rMr_new

    while iterCount < maxIter:
        # compute Ap, rMr, pAp, α on GPU
        scalars[2] = 0.0 # pAp
        cp.cuda.Device().synchronize()
        cufunc.computeMatrixProductAndDotProducts()

        # Sync to make sure pAp reduction is finished
        cp.cuda.Device().synchronize()

        # update x,r,z and accumulate rr_new on GPU
        scalars[0] = scalars[1] # rMr_old = rMr_new
        scalars[1] = 0.0 # rMr_new
        cufunc.updateSolutionAndDotProducts()

        # Sync to make sure rr reduction is finished
        cp.cuda.Device().synchronize()

        # fetch rr_new and rr_old = <r, r> from the device
        rMr_new = scalars[1].item()
        rMr_old = scalars[0].item()

        if rMr_new <= tol*tol * rMr_0:
            break
        if (np.isnan(rMr_new)):
            raise ValueError("rMr_new is NaN")
        # update p ← z + beta*p on GPU
        cufunc.updateSearchDirection()

        iterCount += 1

    time4 = time.time()
    print(f"Time to solve: {time4 - time3:.3f} s")

    return iterCount


def solveWaveEquation2(maxIter=1000, tol=1e-12): 
    """
    Solve wave equation using Conjugate Gradient method.
    Assumes all arrays are already stored in Node as doubles.

    NOTATION: rMr and similar variable names stand for the dot product (r^T)Mr.
    """
    global N, NN, data, indices, indptr, b, diag, data3, indices3, indptr3, b3, diag3

    time1 = time.time()
    cufunc.setupCGSolver2(diag, data, indices, indptr, b)
    cufunc.setupCGSolver3(diag3, data3, indices3, indptr3, b3)

    indptr[N] = NN

    its = 1
    def cb(_xk):
        # called once per iteration
        nonlocal its
        its += 1
    A = csr_matrix((data, indices, indptr), shape=(N, N))
    # pcg = PCGFixed(A, diag)
    # x_out = cp.zeros(N, dtype=cp.float64)
    # pcg.solve(b, x_out, x0 = x00, rtol=1e-6, maxiter=maxIter)


    # cufunc.setupCGSolver(data, indices, indptr, b)

    A3 = csr_matrix((data3, indices3, indptr), shape=(N, N))
    # A.sort_indices()

    # A.sum_duplicates()

    # d = cp.where(d != 0, d, 1.0)
    M  = LinearOperator(A.shape,  matvec=lambda x: x * diag,  rmatvec=lambda x: x * diag, dtype=A.dtype)
    M3 = LinearOperator(A3.shape, matvec=lambda x: x * diag3, rmatvec=lambda x: x * diag3, dtype=A3.dtype)
    # def M_mv(x):
    #     cp.multiply(x, diag, out=buf_mv)   # y = D^{-1} x
    #     return buf_mv

    # def M_rmv(x):
    #     cp.multiply(x, diag, out=buf_rmv)  # y = (D^{-1})^T x == D^{-1} x (diagonal)
    #     return buf_rmv

    # M = LinearOperator(A.shape, matvec=M_mv, rmatvec=M_rmv, dtype=A.dtype)
  
    start = time.time()
    # x3, info = cg(A3, b3, x0 = 3*x00 - 3*x01 + x02, M=M3, tol=1e-6, maxiter=maxIter, callback=cb)
    x3, info = cg(A3, b3, x0 = 2*x00 - x01, M=M3, tol=1e-6, maxiter=maxIter, callback=cb)
    x02[:] = x01[:]
    x01[:] = x00[:]
    x00[:] = x3[:]
    # x, info = cg2(A, b, diag=diag, x0=x3.astype(A.dtype), rtol=tol, maxiter=maxIter)
    solver = PCGSolver(A)
    solver.diag[:] = diag
    # x, info = cg(A, b, M=M, x0=x3.astype(A.dtype), tol=tol, maxiter=maxIter)
    x, info = solver.solve(b, x0=x3.astype(A.dtype), rtol=tol, maxiter=maxIter)

    its += info
    end = time.time()
    print(f"GMRES solve time: {end - start:.3f} s, iterations: {its}, time per iteration: {(end - start)/its:.7f} s")

    cufunc.solverCopyResult(x)
    time2 = time.time()
    print(f"Time to setup and solve: {time2 - time1:.3f} s")
    return

    time3 = time.time()
    iterCount = 0
    rMr0      = None

    # Checks that there's only one GPU device, otherwise it will fail
    if cp.cuda.runtime.getDeviceCount() > 1:
        raise ValueError("This function is not YET parallelized for multiple GPUs. Reductions need to be extended to all devices.")


    scalars = cp.zeros(3, dtype=cp.float64) # rMr_old, rMr_new, pAp
    cufunc.initializeScalars(scalars)


    # Compute initial preconditioned residual Mr
    cufunc.initializePreconditionedResidual()
    cp.cuda.Device().synchronize()

    rMr_new = scalars[1].item()
    if (rMr_new <= 1e-30):
        # If the initial residual is already small enough, return 0
        return 0
    
    rMr_0 = rMr_new

    while iterCount < maxIter:
        # compute Ap, rMr, pAp, α on GPU
        scalars[2] = 0.0 # pAp
        cp.cuda.Device().synchronize()
        cufunc.computeMatrixProductAndDotProducts()

        # Sync to make sure pAp reduction is finished
        cp.cuda.Device().synchronize()

        # update x,r,z and accumulate rr_new on GPU
        scalars[0] = scalars[1] # rMr_old = rMr_new
        scalars[1] = 0.0 # rMr_new
        cufunc.updateSolutionAndDotProducts()

        # Sync to make sure rr reduction is finished
        cp.cuda.Device().synchronize()

        # fetch rr_new and rr_old = <r, r> from the device
        rMr_new = scalars[1].item()
        rMr_old = scalars[0].item()

        if rMr_new <= tol*tol * rMr_0:
            break
        if (np.isnan(rMr_new)):
            raise ValueError("rMr_new is NaN")
        # update p ← z + beta*p on GPU
        cufunc.updateSearchDirection()

        iterCount += 1

    time4 = time.time()
    print(f"Time to solve: {time4 - time3:.3f} s")

    return iterCount



def computeMeteorologicalForcing(time, sflux):
    sflux.getWinds(time)

    print("IMPLEMENT THIS: computeMeteorologicalForcing!!!!")
    # stop
#     # Re-scale wind
#     if (nws != 0 && iwindoff!=0) {
#         do i=1,npa
#           windx(i)=windx(i)*windfactor(i)
#           windy(i)=windy(i)*windfactor(i)
#         enddo !i
#       endif; endif

# !...  compute wind stress components (in lat/lon frame if ics=2; in map projection E-N direction if ics=1)
#       dragcmin=1.0d-3*(0.61d0+0.063d0*6.d0)
#       dragcmax=1.0d-3*(0.61d0+0.063d0*50.d0)

# !$OMP parallel default(shared) private(i,wmag,dragcoef,tmp,theta,z0_donelan)

# !$OMP workshare
#       tau=0.d0 !init.
# !$OMP end workshare

# !$OMP do
#       do i=1,npa
#         if(nws==0) then
#           tau(1,i)=0.d0
#           tau(2,i)=0.d0
#         else if(nws==2.and.ihconsv==1.and.iwind_form==0) then !tauxz and tauyz defined
#           if(idry(i)==1) then
#             tau(1,i)=0.d0
#             tau(2,i)=0.d0
#           else !rescale as well
#             tau(1,i)=-tauxz(i)/rho0*rampwind*windfactor(i)**2.d0 !sign and scale difference between stresses tauxz and tau
#             tau(2,i)=-tauyz(i)/rho0*rampwind*windfactor(i)**2.d0
#           endif
#         else !if(nws==1.or.nws>=4.or.nws>=2.and.ihconsv==0.or.iwind_form==-1) then
#           wmag=sqrt(windx(i)**2.d0+windy(i)**2.d0)
#           if(iwind_form==-1) then !P&P
#             dragcoef=1.0d-3*(0.61d0+0.063d0*wmag)
#             dragcoef=min(max(dragcoef,dragcmin),dragcmax)
#             tau(1,i)=dragcoef*0.001293d0*wmag*windx(i)*rampwind
#             tau(2,i)=dragcoef*0.001293d0*wmag*windy(i)*rampwind
#           else if(iwind_form==1) then !Hwang
#             if(wmag<=35.d0) then
#               dragcoef=1.0d-4*(-0.016d0*wmag*wmag+0.967d0*wmag+8.058d0)
#             else
#               dragcoef=2.23d-3*35.0d0/wmag
#             endif
#             dragcoef=max(dragcoef,4.d-4)
#             tau(1,i)=dragcoef*0.001293d0*wmag*windx(i)*rampwind
#             tau(2,i)=dragcoef*0.001293d0*wmag*windy(i)*rampwind
#           endif !iwind_form
#         endif !nws
#       enddo !i=1,npa
# !$OMP end do      



def schism_step(it: int, namelist, bctides, sflux):

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
    


    #  Compute new bed deformation
    cufunc.computeBedDeformation()


    # Earth tidal potential and loading tide at nodes: pre-compute to save time
    cufunc.computeTidalPotential(time, bctides.tfreq, bctides.tamp, bctides.tnf, 
                                 bctides.tear, bctides.jspc)


    # Meteorological forcing
    computeMeteorologicalForcing(time, sflux)  # Lines 458-1014


    # computeTracerStuff()            # Lines 1015-1145


    # # Compute hydraulic transfer blocks together with reading in flux values
    # computeHydraulicTransfer()      # Lines 1146-1330


    # Get boundary conditions from external files
    cufunc.getBoundaryConditions(time)          # Lines 1331-1584


    # # Read in volume/mass sources/sinks
    # readInSourcesSinks()              # Lines 1585-1610


    # # Calculation of cross-section areas and length for flow b.c.
    # computeFlowBcAreas()                   # Lines 1611-1758


    # XXX I think I will combine all BCs into one function
    # Compute new vel. for flow b.c. (uth,vth)
    cufunc.computeFlowBcVelocities(time)            # Lines 1759-1800



    # Bottom drag coefficients for nchi=-1 or 1; Cd and Cdp for nchi=0 already read in
    ####cufunc.computeBottomDrag()                  # Lines 1801-1961


    # bypass solver for transport only XXXX
    # I haven't included this in the code yet


    # bottom drag and vegetation drag
    cufunc.computeBottomDragVegetation()      # Lines 2353-2662



    # computeTurbulence()                # Lines 2663-3839


    # Horizontal viscosity, implemented as a filter (also Shapiro filter)
    cufunc.computeHorizontalViscosity()       # Lines 3840-4160


    # Backtracking/upwind for momentum
    try:
        monitor.debugPlotIdx[0] = 0        
        cufunc.backtrackAllNodes()              # Lines 4161-5029
    except Exception as e:
        print(f"Error in backtrackAllNodes: {e}")



    # Vertical advection 
    # computeVerticalAdvection()        # Lines 5030-5122


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


    # cufunc.printMatrix()    


    #  Solve the wave equation
    solveWaveEquation(maxIter = opt.mxitn0, tol = opt.rtol0)


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
    
    cufunc.ul0ComputeDryElementsAndNodes()
    
    cufunc.ul0UpdateRewettedNodes(0, 0)
    
    cufunc.ul0UpdateSides(0, 0)

    # Adjust mass after level change
    # XXXXXXXXX IMPLEMENT THIS


    # Compute nodal velocities from side velocities
    cufunc.zeroNodalVelocities()
    cufunc.computeNodalVerticalVelocities1()
    cufunc.computeNodalVerticalVelocities2()
    cufunc.computeNodalHorizontalVelocities1()
    cufunc.computeNodalHorizontalVelocities2()
    # cufunc.resetVelocityOpenBndNodes()
    # cufunc.setVelocityOpenBndNodes()



    # var = monitor.getElementVar("z", 0)
    # syncrhonize 
    # cp.cuda.Device().synchronize()
    # res = monitor.plotMesh()
    # if not res:
    #     raise ValueError("Mesh is not valid")

    # var = monitor.getNodeVar("eta")
    # eta[:, it] = var
    # cufunc.getElementVar(monitor.varCodeDict["x"], var, 0) # monitor.varCodeDict["x"]


    # Update bed deformation and depth info  lines 8156-
 

    # Optional calculations


    




    
