
import python.cuda_functions as cufunc

def schism_step(it: int, namelist: dict):

    core = namelist.core
    opt  = namelist.opt

    # SAL option: scale gravity
    cufunc.scaleGravity((1,), (1,), ())

    dt = core.dt

    time = it*dt

    # Broadcast to global module (XXX I think this is not needed)
    # time_stamp=time; it_main=it

    
    # Compute time ramps for boundary elevation forcing, wind and pressure forcing and tidal potential forcing
    cufunc.computeTimeRamps((1,), (1,), (time, opt.nws, opt.if_source, opt.drampbc, opt.drampwind, opt.drampwafo, opt.dramp_ss, opt.dramp))
    


    #  Compute new bed deformation
    cufunc.computeBedDeformation()


    # Earth tidal potential and loading tide at nodes: pre-compute to save time
    computeTidalPotential()


    computeMeteorologicalForcing()  # Lines 458-1014


    # computeTracerStuff()            # Lines 1015-1145


    # # Compute hydraulic transfer blocks together with reading in flux values
    # computeHydraulicTransfer()      # Lines 1146-1330


    # Get boundary conditions from external files
    getBoundaryConditions()          # Lines 1331-1584


    # # Read in volume/mass sources/sinks
    # readInSourcesSinks()              # Lines 1585-1610


    # Calculation of cross-section areas and length for flow b.c.
    computeFlowBcAreas()                   # Lines 1611-1758


    # Compute new vel. for flow b.c. (uth,vth)
    computeFlowBcVelocities()            # Lines 1759-1800



    # Bottom drag coefficients for nchi=-1 or 1; Cd and Cdp for nchi=0 already read in
    computeBottomDrag()                  # Lines 1801-1961


    # bypass solver for transport only XXXX
    # I haven't included this in the code yet

    # bottom drag and vegetation drag
    computeBottomDragVegetation()      # Lines 2353-2662



    # computeTurbulence()                # Lines 2663-3839


    # Horizontal viscosity, implemented as a filter (also Shapiro filter)
    computeHorizontalViscosity()       # Lines 3840-4160


    # Backtracking/upwind for momentum
    backtrackAllNodes()             # Lines 4161-5029


    # Vertical advection 
    computeVerticalAdvection()        # Lines 5030-5122


    computeElevationBoundaryConditions(ramp, time)


    computeVerticallyIntegratedVelocity()       # Lines 5122-5320


    # compute vegetation and drag model parameters
    computeVegetationAndDragParameters()


    # Density gradient at element centers
    computeDensityGradient()      


    computeBaroclinicTerm(rampbc)


    # Elem. average ghat1 (in eframe if ics=2; \breve{G} in notes). Dimension: m^2/s
    computeElemAverageGhat()          # Lines 5533-5783



    # setup coefficient matrix, sparsem, for the wave equation
    #!...  No elevation essential b.c. are imposed yet but other b.c. is imposed
    computeWaveEquation()                # Lines 5826-



    #  Solve the wave equation
    solveWaveEquation()


    # Momentum equations	
     
    # Precompute elevation gradient, atmo. pressure and earth tidal potential 
    computeGradients()  # lines -7296


    # Transport    lines 7388-8111

    # Update bed deformation and depth info  lines 8156-

    # Optional calculations

    # I/O