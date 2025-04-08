
import python.cuda_functions as cufunc
import iop.readSflux as readSflux


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
    getBoundaryConditions(time)          # Lines 1331-1584


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