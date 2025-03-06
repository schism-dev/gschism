


def schism_step(it: int):


    # SAL option: scale gravity
    scale_gravity()


    time = it*dt

    # Broadcast to global module (XXX I think this is not needed)
    # time_stamp=time; it_main=it

    # Define ramp function for boundary elevation forcing, wind and pressure forcing and tidal potential forcing
    if ibc == 0:
        if drampbc > 0.0:
            rampbc = cp.tanh(2.0*time/86400.0/drampbc)
        else:
            rampbc=1.0

        if nws != 0 and drampwind > 0.0:
            rampwind = cp.tanh(2.0*time/86400.0/drampwind)
        else:
            rampwind = 1.0


    if drampwafo > 0.0:
        rampwafo = cp.tanh(2.0*time/86400.0/drampwafo)
    else:
        rampwafo = 1.0


    # For source/sinks
    if if_source != 0:
        if dramp_ss > 0.0:
            ramp_ss = cp.tanh(2.0*time/86400.0/dramp_ss)
        else:
            ramp_ss = 1.0


    if dramp > 0.0:
        ramp = cp.tanh(2.0*time/86400.0/dramp)
    else:
        ramp = 1.0


    #  Compute new bed deformation
    bdef2[:] = bdef[:]/real(ibdef,rkind)*real(min0(it,ibdef),rkind)




    # !...  Earth tidal potential and loading tide at nodes: pre-compute to save time
    #   do i=1,npa
    #     etp(i)=0.d0
    #     do j=1,ntip
    #       ncyc=int(tfreq(j)*time/2.d0/pi)
    #       arg=tfreq(j)*time-real(ncyc,rkind)*2.d0*pi+jspc(j)*xlon(i)+tear(j)
    #       etp(i)=etp(i)+0.69d0*ramp*tamp(j)*tnf(j)*fun_lat(jspc(j),i)*cos(arg)

    #       if(iloadtide==1) then !loading tide
    #         etp(i)=etp(i)+rloadtide(1,j,i)*cos(tfreq(j)*time-rloadtide(2,j,i))
    #       endif !iloadtide
    #     enddo !j
    #   enddo !i
