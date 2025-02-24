


def schism_step(it: int):

    print('SAL option: scale gravity [NOT IMPLEMENTED]')
    grav2=grav

    time = it*dt

    # Broadcast to global module (XXX I think this is not needed)
    # time_stamp=time; it_main=it

    # Define ramp function for boundary elevation forcing, wind and pressure forcing and tidal potential forcing
    # XXX IMPLEMENT THIS
    # if ibc == 0:
    #     if drampbc > 0.0:
    #         rampbc = cp.tanh(2.0*time/86400.0/drampbc)
    #     else:
    #         rampbc=1.0

    #       if(nws/=0.and.drampwind>0.d0) then
    #     rampwind=tanh(2.d0*time/86400.d0/drampwind)
    #   else
    #     rampwind=1.d0
    #   endif

    #   if(drampwafo>0.d0) then
    #     rampwafo=tanh(2.d0*time/86400.d0/drampwafo)
    #   else
    #     rampwafo=1.d0
    #   endif

    #   !For source/sinks
    #   if(if_source/=0) then
    #     if(dramp_ss>0.d0) then
    #       ramp_ss=tanh(2.d0*time/86400.d0/dramp_ss)
    #     else
    #       ramp_ss=1.d0
    #     endif
    #   endif

    #   if(dramp>0.d0) then
    #     ramp=tanh(2.d0*time/86400.d0/dramp)
    #   else
    #     ramp=1.d0
    #   endif



    #  Compute new bed deformation
    #   do i=1,npa
    #     bdef2(i)=bdef(i)/real(ibdef,rkind)*real(min0(it,ibdef),rkind)
    #   enddo !i



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
