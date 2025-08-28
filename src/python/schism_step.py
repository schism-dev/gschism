import python.cuda_functions as cufunc
import iop.readSflux as readSflux
import iop.writeOutput as writeOutput
import cupy as cp
import debug.monitor as monitor
# from cupy import cudaMemcpyFromSymbol


def solveWaveEquation(maxIter=1000, tol=1e-12): 
    """
    Solve wave equation using Conjugate Gradient method.
    Assumes all arrays are already stored in Node as doubles.

    NOTATION: rMr and similar variable names stand for the dot product (r^T)Mr.
    """
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

        # fetch rr_new = <r, r> from the device
        rMr_new = scalars[1].item()
        rMr_old = scalars[0].item()

        if rMr_new <= tol*tol * rMr_0:
            break
        if (cp.isnan(rMr_new)):
            raise ValueError("rMr_new is NaN")
        # update p ← z + beta*p on GPU
        cufunc.updateSearchDirection()

        iterCount += 1

        

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

def checkU():
    #  REMOVE THIS!!!!!!
    import numpy as np
    # cp.cuda.Device().synchronize()
    maxVar = 0.0
    max_idx = 0
    max_level = 0
    var1 = monitor.getSideVar("u", 0)
    for i in range(1, 2):
        var2 = monitor.getSideVar("u", i)
        diff = np.abs(var1 - var2)
        max_diff = np.max(diff)
        if max_diff > maxVar:
            maxVar = max_diff
            max_level = i
            max_idx = np.argmax(diff)
        var1 = var2.copy()
    if maxVar > 0.4:
        print(f"Max difference: {maxVar} at node {max_idx}, level {max_level}")

    maxVar = 0.0
    max_idx = 0
    max_level = 0
    var1 = monitor.getNodeVar("u", 0)
    for i in range(1, 2):
        var2 = monitor.getNodeVar("u", i)
        diff = np.abs(var1 - var2)
        max_diff = np.max(diff)
        if max_diff > maxVar:
            maxVar = max_diff
            max_level = i
            max_idx = np.argmax(diff)
        var1 = var2.copy()
    if maxVar > 0.4:
        print(f"Max difference: {maxVar} at node {max_idx}, level {max_level}")

def schism_step(it: int, namelist, bctides, sflux, eta):

    core = namelist.core
    opt  = namelist.opt

    # SAL option: scale gravity
    cufunc.scaleGravity(opt.loadtide_coef)

    dt = core.dt

    time = (it + 1)*dt

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

    # checkU()
    cufunc.solveMomentumEquations(time)
    # checkU()
    

    cufunc.applyMomentumBC()
    # checkU()

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


    # I/O
    writeOutput.writeOutput(it)   




    
