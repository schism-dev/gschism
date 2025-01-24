/**
 * Tracks a particle through elements and vertical levels in the mesh
 * 
 * This function implements a "jump-and-walk" algorithm to find which element contains
 * a particle after it moves from a starting position to a target position.
 * 
 * @param mesh         Reference to mesh object containing grid information
 * @param x0,y0,z0     Starting position coordinates
 * @param xt,yt,zt     Target position coordinates
 * @param dtb          Time step for tracking
 * @param uuint,vvint,wwint Initial velocity components
 * @param jlev         Current/resulting vertical level
 * @param nextElem     Starting/resulting element index
 * 
 * @return Status code: 
 *         0 = Success, found containing element
 *         1 = Hit boundary or special case (dry element, death trap, etc)
 *         2 = Hit subdomain boundary upon entry
 *         3 = Hit subdomain boundary during tracking
 */
__device__ int trackThroughElementsAndLevels(Mesh& mesh,
                                             double x0, double y0, double z0,
                                             double xt, double yt, double zt, double dtb,
                                             double uuint, double vvint, double wwint, 
                                             int& jlev, Element *nextElem)
{
    Element *curElem = nextElem;

    // 1. Initial checks
    // Verify starting element is valid (wet)
    if (curElem->isDry) {
        return 1;  // Error code for invalid starting conditions
    }

    // Checks if track has a zero length.
    double trm = time; //time remaining


    // Nudging weights
    constexpr double nudgeW1 = 0.01;  // Nudging weights
    constexpr double nudgeW2 = 1.0 - nudgeW1;  // Nudging weights

    double xctr3, yctr3;
    if (ics == 1) {
        double xctr3 = curElem->xc;
        double yctr3 = curElem->yc;
    } else {
        // TODO: Handle lat/lon projection case
        // project_pt('g2l', xctr, yctr, zctr, gcor0, frame0, xctr3, yctr3, tmp)
    }
    
    // When close to an element side, nudge starting point towards element centroid to avoid ambiguity (nudgeW1 + nudgeW2 = 1.0)
    const double x0 = nudgeW2*x0 + nudgeW1*xctr3;
    const double y0 = nudgeW2*y0 + nudgeW1*yctr3;

    double xcg=x0, ycg=y0;
    double pathl2 = (xt - xcg)*(xt - xcg) + (yt - ycg)*(yt - ycg);
    
    if (pathl2==0.0 || trm==0.0) {
        return 2; // Error code for zero length track
    }

    // Check if target point is already in starting element. If it is not,
    // then we need to find the element that contains the target point.
    if (!curElem->containsPoint(xt, yt, zt)) {


        // Compute path from start to target
        // Find which edge the path intersects
        // Handle case where no intersection is found
        int intersectingEdge;
        double xIntersect, yIntersect;
        if (!curElem->findPathIntersection(x0, y0, xt, yt, 
                                        intersectingEdge, 
                                        xIntersect, yIntersect)) {
            return 3;  // No valid intersection found
        }
        
        // 3. Main tracking loop
        // - While not found target element:
        //   a. Update vertical position along path
        //   b. Check for boundary conditions
        //   c. Find next intersecting edge
        //   d. Move to neighboring element
        //   e. Check for special cases (dry elements, death traps)
        //   f. Handle coordinate system transformations if needed
        // Main tracking loop
        const int MAX_ITERATIONS = 1000;
        int iterations = 0;
        double xCurrent = x0, yCurrent = y0;

        
        while (iterations++ < MAX_ITERATIONS) {
            // Check if current point is in current element
            if (curElem->containsPoint(xt, yt)) {
                nextElem = curElem;
                return 0;
            }

            // Find intersection with current element
            int intersectingEdge;
            double xIntersect, yIntersect;
            if (!curElem->findPathIntersection(xCurrent, yCurrent, xt, yt, 
                                                intersectingEdge, nextElem,
                                                xIntersect, yIntersect)) {
                return 4;  // No intersection found
            }

            if (curElem->isEdgeBoundary(intersectingEdge)) {
                return 5;
            }

            if (nextElement->isDry) {
                return 6;
            }

            // Update position and element for next iteration
            xCurrent = xIntersect;
            yCurrent = yIntersect;
            curElem = nextElem;
        }

        if (iterations >= MAX_ITERATIONS) return 7;  // Hit iteration limit
    }

    // 4. Final position processing
    // - Compute final vertical level
    // - Calculate shape functions/area coordinates
    // - Interpolate z-coordinates
        // Final position processing
    double areaCoords[4];  // For both tri (3) and quad (4) cases
    if (!currentElem->computeShapeCoords(xt, yt, shapeCoords)) {
        return 8;  // Failed to compute area coordinates
    }

    // Get vertical coordinates at each node
    double zCoords[MAX_NVRT];  // Array for vertical coordinates
    int bottomLevel = currentElem->computeVerticalProfile(shapeCoords, zCoords);
    
    // Find vertical level containing zt
    if (!currentElem->findVerticalLevel(zt, zCoords, bottomLevel, jlev)) {
        return 9;  // Failed to find vertical level
    }


    // 5. Return appropriate status code
    // - Success if target found
    // - Various error/boundary conditions
        // Interpolate velocities at final position
    if (!currentElem->interpolateVelocities(shapeCoords, jlev, 
                                            uuint, vvint, wwint)) {
        return 10;
    }

    return 0;
}



      subroutine quicksearch(idx,itr,l_ns,ipsgb,gcor0,frame0,time,x0,y0,z0,nnel, &
     &jlev,xt,yt,zt,trm,nfl,kbpl,arco,zrat,ztmp,vis_coe,uuint0,vvint0,wwint0, &
     &uuint,vvint,wwint)
!
!********************************************************************************
!										*
!     Straightline search algorithm. For ics=2, this is done in the local frame of the starting pt. 
!     
!     Inputs: 
!       idx: ID identifying where quicksearch is called (info only);
!       itr: iteration # in btrack (info only);
!       l_ns: side if l_ns=2; element if l_ns=3; (info only)
!       ipsgb: global originating node or side or element index; (info only)
!       gcor0(3): global coord. of the originating pt (for ics=2);
!       frame0(3,3): frame tensor at originating pt (2nd index is axis id) (for ics=2);
!       vis_coe: weighting factor between continuous and discontinuous vel. - not used anymore
!       time: time step from (x0,y0,z0) to (xt,yt,zt);
!       x0,y0,z0:  starting pt. (x0,y0) must be inside nnel 
!                  (for ics=2, z0 is in vertical direction and (x0,y0,z0) are in frame0); 
!       nnel,jlev: starting element and level. nnel must be inside aug. domain.
!       xt,yt,zt: projected end pt; (for ics=2, zt is in vertical direction)
!       uuint0,vvint0,wwint0: vel. used in 'time' (info only). In frame0 if ics=2.
!       In addition, su2,sv2,ww2 are also used.
! 
!     Outputs:
!      nnel, jlev: end element and level. nnel must be inside aug. domain;
!      (xt,yt,zt):  the updated end pt (if so); 
!      trm: time remaining (trm<=time). trm=0 unless the path crosses aug. bnd (nfl=3);
!      nfl: a flag. nfl=1 if a bnd or dry element is hit and vel. there is small,              
!           or death trap is reached, or some extraordinary case like unable to
!           find intersecting side. In this case, all outputs are valid, 
!           and the time stepping in btrack is exited successfully. 
!           If nfl=2 (hit aug. bnd upon entry) or 3 (hit aug. bnd during iteration),
!           nnel, jlev, (xt,yt,zt) and trm are updated, and nnel should be inside 
!           one of the neighbor process. If nfl=2, (xt,yt,zt) are moved to (x0,y0,z0).

!      Following outputs are valid only if nfl<=1:
!      kbpl: the local bottom level at (xt,yt);
!      arco(4): shape functions of (xt,yt);
!      zrat: vertical ratio of zt (=0 when zt=ztmp(jlev));
!      ztmp(nvrt):  z-coords. at (xt,yt);
!      uuint,vvint,wwint: interpolated vel. at end pt. In frame0 if ics=2.
!********************************************************************************
!
      use schism_glbl
      use schism_msgp, only : myrank,parallel_abort
      use hydraulic_structures, only: nhtblocks
      implicit none

      integer, intent(in) :: idx,itr,l_ns,ipsgb
      real(rkind), intent(in) :: gcor0(3),frame0(3,3),time,x0,y0,z0,vis_coe,uuint0,vvint0,wwint0
      integer, intent(inout) :: nnel,jlev
      real(rkind), intent(inout) :: xt,yt,zt
      integer, intent(out) :: nfl,kbpl
      real(rkind), intent(out) :: trm,arco(4),zrat,ztmp(nvrt),uuint,vvint,wwint

      !Function
      real(rkind) :: signa1

      !Local
      logical :: qsearch_done
      integer :: jk(4)
      real(rkind) :: wild(10,2),wild2(10,2),wild3(4,2) !,xy_l(3,2)
      real(rkind) :: vxl(4,2),vyl(4,2),vzl(4,2),vxn(4),vyn(4),vzn(4),ztmp2(nvrt,4)

      integer :: nnel00,jlev00,nel,i,j,l,nel_j,jd1,jd2,iflag,it,md1,md2,lit, &
                 &isd,n1,n2,n3,kin,nd,lev,k,inside1,inside2
      real(rkind) :: xt00,yt00,zt00,xcg,ycg,zcg,xcg2,ycg2,zcg2,pathl,ar_min1, &
                     &ar_min2,xn1,yn1,zn1,xn2,yn2,zn2,ar1,ar2,xin,yin,zin,tt1, &
                     &tt2,xt2,yt2,zt2,xtmp,ytmp,eps,dist,tmp,xctr3,yctr3,vtan, &
                     &xvel,yvel,zvel,hvel,etal,dep,hmod2,uj,vj,uj1,vj1,uu,vv,uf,vf

!     Debug
!      fdb='qs_000000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!      open(98,file=out_dir(1:len_out_dir)//fdb,status='unknown')

!     Save for debug
      xt00=xt
      yt00=yt
      nnel00=nnel
      jlev00=jlev

!     Assumptions used in this routine:
!     (1) starting element and all crossing elements are wet
!     (2) starting pt is inside nnel; start and end pts are distinct;
!     (3) the trajectory intersects sides at 1 pt only, and when the
!         end pt (xt,yt) is near a side it's nudged into element to conclude
!         the process - this won't create CPU dependency

      if(idry_e(nnel)==1) then
        write(errmsg,*)'QUICKSEARCH: Starting element is dry:',idry_e(nnel)
        call parallel_abort(errmsg)
      endif

!     Initialize (moving) starting pt
      nfl=0
      trm=time !time remaining
      nel=nnel
      xcg=x0; ycg=y0
      pathl=sqrt((xt-xcg)*(xt-xcg)+(yt-ycg)*(yt-ycg))
      if(pathl==0._rkind.or.trm==0._rkind) then
!        write(12,*)'Last QUICKSEARCH: nodes'
!        do i=1,npa
!          do k=1,nvrt
!            write(12,*)iplg(i),k,uu2(k,i),vv2(k,i)
!          enddo !k
!        enddo !i
!        write(12,*)'Sides:'
!        do i=1,nsa
!          do k=1,nvrt
!            write(12,*)islg(i),iplg(isidenode(1:2,i)),k,su2(k,i),sv2(k,i)
!          enddo !k
!        enddo !i
        write(errmsg,*)'QUICKSEARCH: Zero path',idx,itr,l_ns,ipsgb,ielg(nel),jlev, &
     &x0,y0,xt,yt,xcg,ycg,time,uuint0,vvint0,wwint0
        call parallel_abort(errmsg)
      endif

!     Check start and end pts
      qsearch_done=.false.
      if(i34(nel)==3) then 
        call area_coord(0,nel,gcor0,frame0,xcg,ycg,arco)
        ar_min1=minval(arco(1:3)) !info for debug only
        call area_coord(0,nel,gcor0,frame0,xt,yt,arco)
        ar_min2=minval(arco(1:3))
        if(ar_min2>-small1) then !found & finish
          !Fix A.C. for 0/negative
          if(ar_min2<=0) call area_coord(1,nel,gcor0,frame0,xt,yt,arco)
          nnel=nel
          trm=0._rkind
          qsearch_done=.true.
        endif
      else !quad
        !Reproject pt in eframe of nel for ics=2
        if(ics==1) then
          xn1=xcg; yn1=ycg
          xn2=xt; yn2=yt
        else !ll
          call project_pt('l2g',xcg,ycg,0._rkind,gcor0,frame0,xcg2,ycg2,zcg2)
          call project_pt('g2l',xcg2,ycg2,zcg2,(/xctr(nel),yctr(nel),zctr(nel)/),eframe(:,:,nel),xn1,yn1,zn1)
          call project_pt('l2g',xt,yt,0._rkind,gcor0,frame0,xcg2,ycg2,zcg2)
          call project_pt('g2l',xcg2,ycg2,zcg2,(/xctr(nel),yctr(nel),zctr(nel)/),eframe(:,:,nel),xn2,yn2,zn2)
        endif !ics
        call quad_shape(0,1,nel,xn1,yn1,inside1,arco) !info only
        ar_min1=minval(arco) !info for debug only
        call quad_shape(0,2,nel,xn2,yn2,inside2,arco)
        ar_min2=minval(arco)
        if(inside2/=0) then !found & finish
          nnel=nel
          trm=0._rkind
          qsearch_done=.true.
        endif
      endif !i34

      if(.not.qsearch_done) then
!======================================================================
!     (xt,yt) not in nel, and thus (x0,y0) and (xt,yt) are distinctive
!     Find starting edge nel_j
!     Try this twice to account for underflow (e.g. inter-btrack etc),
!     the 2nd try may create CPU dependency but it occurs very rarely
      loop6: do i=1,2
        wild=0._rkind; wild2=0._rkind !initialize for debugging output
        nel_j=0
        do j=1,i34(nel) !sides
          jd1=elnode(nxq(1,j,i34(nel)),nel)
          jd2=elnode(nxq(2,j,i34(nel)),nel)
          if(ics==1) then
            xn1=xnd(jd1); yn1=ynd(jd1)
            xn2=xnd(jd2); yn2=ynd(jd2)
          else !lat/lon
            call project_pt('g2l',xnd(jd1),ynd(jd1),znd(jd1),gcor0,frame0,xn1,yn1,zn1)
            call project_pt('g2l',xnd(jd2),ynd(jd2),znd(jd2),gcor0,frame0,xn2,yn2,zn2)
          endif !ics
          wild3(j,1)=xn1; wild3(j,2)=yn1 !save for computing centroid and nudging later
          ar1=signa1(xcg,xn1,xt,ycg,yn1,yt)    
          ar2=signa1(xcg,xt,xn2,ycg,yt,yn2)    
          wild2(j,1)=ar1; wild2(j,2)=ar2
          if(ar1>0._rkind.and.ar2>0._rkind) then
            call intersect2(xcg,xt,xn1,xn2,ycg,yt,yn1,yn2,iflag,xin,yin,tt1,tt2)
            wild(j,1)=tt1; wild(j,2)=tt2; wild(3+j,1)=xin; wild(3+j,2)=yin
            if(iflag/=1) then !abnormal return
              if(ics==1) then
                xcg2=xcg; ycg2=ycg; zcg2=0._rkind; xt2=xt; yt2=yt; zt2=0._rkind
              else !lat/lon
                call project_pt('l2g',xcg,ycg,0._rkind,gcor0,frame0,xcg2,ycg2,zcg2)
                call project_pt('l2g',xt,yt,0._rkind,gcor0,frame0,xt2,yt2,zt2)
              endif !ics
              write(12,*)'QUICKSEARCH: Found no intersecting edges (1):',idx,itr, &
     &ielg(nel),xcg2,ycg2,zcg2,xt2,yt2,zt2,ar_min1,ar_min2,wild(1:3,1:2),wild(4:6,1:2),ar1,ar2, &
     &xcg,ycg,xt,yt,time,trm,jlev,uuint0,vvint0,wwint0,jlev00
              !Set final values and return
              nfl=1
              trm=0._rkind
              arco(:)=1.d0/i34(nnel)
              uuint=0._rkind
              vvint=0._rkind
              wwint=0._rkind
              zrat=1.d0
              kbpl=nvrt-1
              ztmp=znl(:,elnode(1,nnel)) !1st node
              return
            else !success
              nel_j=j; exit loop6
            endif
          endif !ar1>=0.and.ar2>=0
        enddo !j=1,i34

        if(nel_j==0) then
          if(ics==1) then
            xcg2=xcg; ycg2=ycg; zcg2=0._rkind; xt2=xt; yt2=yt; zt2=0._rkind
          else !lat/lon
            call project_pt('l2g',xcg,ycg,0._rkind,gcor0,frame0,xcg2,ycg2,zcg2)
            call project_pt('l2g',xt,yt,0._rkind,gcor0,frame0,xt2,yt2,zt2)
          endif !ics

          if(i==1) then !1st try
            write(12,*)'QUICKSEARCH: no intersecting edge; start ID (node/side/elem)=', &
     &l_ns,'; start gb. node/side/elem #=',ipsgb,'; start level=',jlev,'; current elem=',ielg(nel), &
     &'; cg (local) coord.=',xcg2,ycg2,zcg2,'; end coord.=',xt2,yt2,zt2, &
     &'; signed areas (cg,1,t)@ nodes followed by (cg,t,2)@ nodes=',wild2(1:i34(nel),1:2), &
     &'; xcg,ycg,xt,yt=',xcg,ycg,xt,yt, &
     &'; time step from cg to t=',time,'; time remaining=',trm, &
     &'; min. area coord. for cg, t=',ar_min1,ar_min2,'; input vel=',uuint0,vvint0,wwint0, &
     &idx,itr,jlev00
            !Nudge (xcg,ycg) to off centroid (to escape some tricky underflow)
            if(i34(nel)==3) then
              wild(1,1)=real(1./3.-1.12e-2,rkind); wild(2,1)=real(1./3.-1.09e-2,rkind); wild(3,1)=1._rkind-wild(1,1)-wild(2,1) !A.C.
            else
              wild(1,1)=real(0.25-1.12e-2,rkind); wild(2,1)=real(0.25-1.09e-2,rkind); 
              wild(3,1)=real(0.25+0.937e-2,rkind); wild(4,1)=1._rkind-sum(wild(1:3,1))
            endif !i34
            xtmp=dot_product(wild(1:i34(nel),1),wild3(1:i34(nel),1))
            ytmp=dot_product(wild(1:i34(nel),1),wild3(1:i34(nel),2))
            eps=real(1.019e-2,rkind)
            xcg=(1._rkind-eps)*xcg+eps*xtmp
            ycg=(1._rkind-eps)*ycg+eps*ytmp

            if(i34(nel)==3) then
              call area_coord(0,nel,gcor0,frame0,xcg,ycg,arco)
            else
              if(ics==1) call quad_shape(0,3,nel,xcg,ycg,inside1,arco)
            endif !i34
            ar_min1=minval(arco(1:i34(nel))) !info for debug only (undefined if ics=2 and quads)
          else !i=2; out of luck; return
            write(12,*)'QUICKSEARCH: no intersecting edge; start ID (node/side/elem)=', &
     &l_ns,'; start gb. node/side/elem #=',ipsgb,'; start level=',jlev,'; current elem=',ielg(nel), &
     &'; cg (local) coord.=',xcg2,ycg2,zcg2,'; end coord.=',xt2,yt2,zt2, &
     &'; signed areas (cg,1,t)@ nodes followed by (cg,t,2)@ nodes=',wild2(1:i34(nel),1:2), &
     &'; xcg,ycg,xt,yt=',xcg,ycg,xt,yt, &
     &'; time step from cg to t=',time,'; time remaining=',trm, &
     &'; min. area coord. for cg, t=',ar_min1,ar_min2,'; input vel=',uuint0,vvint0,wwint0, &
     &idx,itr,jlev00
!            call parallel_abort(errmsg)
            nfl=1
            trm=0._rkind
            arco(:)=1.d0/i34(nnel)
            uuint=0._rkind
            vvint=0._rkind
            wwint=0._rkind
            zrat=1.d0
            kbpl=nvrt-1
            ztmp=znl(:,elnode(1,nnel)) !1st node
            return
          endif !i
        endif !nel_j=0
      enddo loop6 !i: 2 tries

!     Check aug. exit
!     nnel, jlev, and trm are unchanged; (xt,yt,zt) moved to (x0,y0,z0)
!     to be ready for inter-subdomain tracking
      if(ic3(nel_j,nel)<0) then
        xt=x0; yt=y0; zt=z0; nnel=nel
        nfl=2; return
      endif

      zin=z0 !intialize
      it=0
      loop4: do
!----------------------------------------------------------------------------------------
      it=it+1

!     Exit loop if death trap is reached
      if(it>1000) then
!        if(ifort12(3)==0) then
!          ifort12(3)=1
        write(12,*)'QUICKSEARCH: Death trap reached'
!        endif
        nfl=1
        xt=xin
        yt=yin
        zt=zin
        nnel=nel
        trm=0._rkind !min(trm,time)
        exit loop4
      endif
      md1=elnode(nxq(1,nel_j,i34(nel)),nel)
      md2=elnode(nxq(2,nel_j,i34(nel)),nel)
      
!     Compute z position 
      dist=sqrt((xin-xt)*(xin-xt)+(yin-yt)*(yin-yt))
      tmp=min(1._rkind,dist/pathl)
      zin=zt-tmp*(zt-zin)
      trm=trm*tmp !time remaining
      pathl=dist !sqrt((xin-xt)**2+(yin-yt)**2)
      if(dist==0._rkind.or.trm==0._rkind) then !done
!        write(errmsg,*)'QUICKSEARCH: end pt on side:',idx,itr,l_ns,ipsgb,dist,it, &
!     &ielg(nnel00),ielg(nel),zin,time,x0,y0,xt00,yt00,xin,yin,xt,yt, &
!     &uuint0,vvint0,wwint0,ar_min2,jlev00,vis_coe
!        call parallel_abort(errmsg)

        if(i34(nel)==3) then
          call area_coord(1,nel,gcor0,frame0,xt,yt,arco) !'1' - fix A.C.
        else
          !Reproject pt in eframe of nel for ics=2
          if(ics==1) then
            xn2=xt; yn2=yt
          else !ll
            call project_pt('l2g',xt,yt,0._rkind,gcor0,frame0,xcg2,ycg2,zcg2)
            call project_pt('g2l',xcg2,ycg2,zcg2,(/xctr(nel),yctr(nel),zctr(nel)/),eframe(:,:,nel),xn2,yn2,zn2)
          endif !ics

          call quad_shape(1,4,nel,xn2,yn2,inside2,arco)
        endif
        nnel=nel
        trm=0._rkind
        exit loop4
      endif

!     Check for aug. exit
      if(ic3(nel_j,nel)<0) then
!       nnel is the last element inside aug. domain
!       xt,yt,zt, jlev, and trm are updated.
!       IMPORTANT: as new (xt,yt) is right on a side, make sure
!       that the vel. is pointing away from nnel!
        nfl=3
        xt=xin
        yt=yin
        zt=zin
        nnel=nel
        trm=min(trm,time) !>0
        nnel=nel
        return
      endif

!     Next element is inside aug. domain
      lit=0 !flag
!     For horizontal exit and dry elements, compute tangential vel.,
!     update target (xt,yt,zt) and continue.
!     max() added for bound check
      if(ic3(nel_j,nel)==0.or.idry_e(max(1,ic3(nel_j,nel)))==1) lit=1
      if(ihydraulics/=0.and.nhtblocks>0) then
        if(isblock_sd(1,elside(nel_j,nel))>0) lit=1 !active block
      endif
   
      if(lit==1) then
        isd=elside(nel_j,nel)
        if(isidenode(1,isd)+isidenode(2,isd)/=md1+md2) then
          write(errmsg,*)'QUICKSEARCH: Wrong side'
          call parallel_abort(errmsg)
        endif

!       Nudge intersect (xin,yin), and update starting pt
        eps=real(1.e-2,rkind) !100*small2
        if(ics==1) then
          xctr3=xctr(nel); yctr3=yctr(nel)
        else !lat/lon
!          eps=small2 !need more accuracy for south pole region
          call project_pt('g2l',xctr(nel),yctr(nel),zctr(nel),gcor0,frame0,xctr3,yctr3,tmp)
        endif !ics
        xin=(1._rkind-eps)*xin+eps*xctr3 !xctr(nel)
        yin=(1._rkind-eps)*yin+eps*yctr3 !yctr(nel)
        xcg=xin
        ycg=yin
   
        vtan=-su2(jlev,isd)*sny(isd)+sv2(jlev,isd)*snx(isd)
        !If open bnd is hit, optionally stop with nfl=1
        !if(ibtrack_openbnd/=0.and.isbs(isd)>0) vtan=0._rkind
        if(isbs(isd)>0) vtan=0._rkind
        xvel=-vtan*sny(isd)
        yvel=vtan*snx(isd)

        zvel=(ww2(jlev,md1)+ww2(jlev,md2))/2._rkind
        xt=xin-xvel*trm
        yt=yin-yvel*trm
        zt=zin-zvel*trm
        hvel=sqrt(xvel*xvel+yvel*yvel)
        if(hvel<=velmin_btrack) then
          nfl=1
          xt=xin
          yt=yin
          zt=zin
          nnel=nel
          trm=0._rkind
          exit loop4
        endif
        pathl=hvel*trm
      endif !abnormal cases

!     Search for nel's neighbor with edge nel_j, or in abnormal cases, the same element
      if(lit==0) nel=ic3(nel_j,nel) !next front element

      if(i34(nel)==3) then
        call area_coord(0,nel,gcor0,frame0,xt,yt,arco)
        ar_min1=minval(arco(1:3))

        if(ar_min1>-small1) then !found
          !arco will be fixed immediately outside loop4
          !if(ar_min1<=0) call area_coord(1,nel,gcor0,frame0,xt,yt,arco) !Fix
          nnel=nel
          trm=0._rkind
          exit loop4
        endif
      else !quad
        !Reproject pt in eframe of nel for ics=2
        if(ics==1) then
          xn2=xt; yn2=yt
        else !ll
          call project_pt('l2g',xt,yt,0._rkind,gcor0,frame0,xcg2,ycg2,zcg2)
          call project_pt('g2l',xcg2,ycg2,zcg2,(/xctr(nel),yctr(nel),zctr(nel)/),eframe(:,:,nel),xn2,yn2,zn2)
        endif !ics

        call quad_shape(0,5,nel,xn2,yn2,inside2,arco)
        ar_min1=minval(arco(1:4))
        if(inside2/=0) then !found
          !arco will be fixed immediately outside loop4
          !call quad_shape(1,?,nel,xt,yt,inside2,arco) !force the pt inside
          nnel=nel
          trm=0._rkind
          exit loop4
        endif
      endif !i34

!     Next intersecting edge
      wild=0._rkind; wild2=0._rkind !initialize for output
      nel_j=0
      do j=1,i34(nel)
        jd1=elnode(nxq(1,j,i34(nel)),nel)
        jd2=elnode(nxq(2,j,i34(nel)),nel)
!       For abnormal case, same side (border side) cannot be hit again
        if(jd1==md1.and.jd2==md2.or.jd2==md1.and.jd1==md2) cycle
        if(ics==1) then
          xn1=xnd(jd1); yn1=ynd(jd1)
          xn2=xnd(jd2); yn2=ynd(jd2)
        else !lat/lon
          call project_pt('g2l',xnd(jd1),ynd(jd1),znd(jd1),gcor0,frame0,xn1,yn1,tmp)
          call project_pt('g2l',xnd(jd2),ynd(jd2),znd(jd2),gcor0,frame0,xn2,yn2,tmp)
        endif !ics
        ar1=signa1(xcg,xn1,xt,ycg,yn1,yt)
        ar2=signa1(xcg,xt,xn2,ycg,yt,yn2)
        wild2(j,1)=ar1; wild2(j,2)=ar2
        if(ar1>0._rkind.and.ar2>0._rkind) then
          call intersect2(xcg,xt,xn1,xn2,ycg,yt,yn1,yn2,iflag,xin,yin,tt1,tt2)
          wild(j,1)=tt1; wild(j,2)=tt2; wild(3+j,1)=xin; wild(3+j,2)=yin
          if(iflag/=1) then !abort and return
            if(ics==1) then
              xcg2=xcg; ycg2=ycg; zcg2=0._rkind; xt2=xt; yt2=yt; zt2=0._rkind
            else !lat/lon
              call project_pt('l2g',xcg,ycg,0._rkind,gcor0,frame0,xcg2,ycg2,zcg2)
              call project_pt('l2g',xt,yt,0._rkind,gcor0,frame0,xt2,yt2,zt2)
            endif !ics      
            write(12,*)'QUICKSEARCH: Failed to find next edge (2):',lit,idx,itr,l_ns,ipsgb, &
     &xcg2,ycg2,zcg2,xt2,yt2,zt2,ielg(nel),iplg(md1),iplg(md2),ar_min1, &
     &wild(1:3,1:2),wild(4:6,1:2),ar1,ar2,xcg,ycg,xt,yt,time,trm,uuint0,vvint0,wwint0
!           call parallel_abort(errmsg)
            !Set final values and return
            nfl=1
            trm=0._rkind
            nnel=nel
            arco(:)=1.d0/i34(nnel)
            uuint=0._rkind
            vvint=0._rkind
            wwint=0._rkind
            zrat=1.d0
            kbpl=nvrt-1
            ztmp=znl(:,elnode(1,nnel)) !1st node
            return
          else !success
            nel_j=j; !next front edge
            cycle loop4
          endif
        endif !ar1>=0.and.ar2>=0
      enddo !j

      if(nel_j==0) then !return
        if(ics==1) then
          xcg2=xcg; ycg2=ycg; zcg2=0._rkind; xt2=xt; yt2=yt; zt2=0._rkind
        else !lat/lon
          call project_pt('l2g',xcg,ycg,0._rkind,gcor0,frame0,xcg2,ycg2,zcg2)
          call project_pt('l2g',xt,yt,0._rkind,gcor0,frame0,xt2,yt2,zt2)
        endif !ics   
        write(12,*)'QUICKSEARCH: no intersecting edge (2): ',idx,itr,l_ns,ipsgb,ielg(nel), &
     &xcg2,ycg2,zcg2,xt2,yt2,zt2,wild2(1:i34(nel),1:2),xcg,ycg,xt,yt,time,trm,uuint0,vvint0,wwint0,lit
!        call parallel_abort(errmsg)
        !Set final values and return
        nfl=1
        trm=0._rkind
        nnel=nel
        arco(:)=1.d0/i34(nnel)
        uuint=0._rkind
        vvint=0._rkind
        wwint=0._rkind
        zrat=1.d0
        kbpl=nvrt-1
        ztmp=znl(:,elnode(1,nnel)) !1st node
        return
      endif !nel_j

!----------------------------------------------------------------------------------------
      end do loop4 
!======================================================================
      endif !(.not.qsearch_done) then

!400   continue
!     No vertical exit from domain
      if(idry_e(nnel)==1) then
        write(errmsg,*)'QUICKSEARCH: Ending element is dry:',ielg(nnel)
        call parallel_abort(errmsg)
      endif

!     Compute area & sigma coord.
      if(i34(nnel)==3) then
        call area_coord(1,nnel,gcor0,frame0,xt,yt,arco)
      else
        !Reproject pt in eframe of nel for ics=2
        if(ics==1) then
          xn2=xt; yn2=yt
        else !ll
          call project_pt('l2g',xt,yt,0._rkind,gcor0,frame0,xcg2,ycg2,zcg2)
          call project_pt('g2l',xcg2,ycg2,zcg2,(/xctr(nnel),yctr(nnel),zctr(nnel)/),eframe(:,:,nnel),xn2,yn2,zn2)
        endif !ics

        call quad_shape(1,6,nnel,xn2,yn2,inside2,arco)
      endif !i34

!     Local z-coord.
      do j=1,i34(nnel) !nodes
        nd=elnode(j,nnel)
        call zcoor(2,nd,jk(j),ztmp2(:,j))
      enddo !j
      kbpl=minval(jk(1:i34(nnel))) !min. bottom index

      ztmp(kbpl:nvrt)=0._rkind
      do j=1,i34(nnel)
        do k=kbpl,nvrt
          ztmp(k)=ztmp(k)+ztmp2(max(k,jk(j)),j)*arco(j)
        enddo !k
      enddo !j

#ifdef DEBUG
      do k=kbpl+1,nvrt
        !todo: assert
        !Warning: can be 0 for degenerate case
        if(ztmp(k)-ztmp(k-1)<0._rkind) then
          write(errmsg,*)'QUICKSEARCH: Inverted z-level in quicksearch:', &
     &ielg(nnel),k,ztmp(k),ztmp(k-1),kbpl,jk(:),ztmp(kbpl:nvrt),';', &
     &arco(1:i34(nnel)),(ztmp2(jk(j):nvrt,j),j=1,i34(nnel)), &
     &eta2(elnode(1:i34(nnel),nnel)),dp(elnode(1:i34(nnel),nnel))
          call parallel_abort(errmsg)
        endif
      enddo !k
#endif

      if(zt<=ztmp(kbpl)) then
        zt=ztmp(kbpl)
        zrat=1._rkind
        jlev=kbpl+1
      else if(zt>=ztmp(nvrt)) then
        zt=ztmp(nvrt)
        zrat=0._rkind
        jlev=nvrt
      else
        jlev=0
        do k=kbpl,nvrt-1
          if(zt>=ztmp(k).and.zt<=ztmp(k+1)) then 
            jlev=k+1
            exit
          endif
        enddo !k
        !todo: assert
        if(jlev==0) then
          write(errmsg,*)'QUICKSEARCH: Cannot find a vert. level:',zt,etal,dep,(ztmp(k),k=kbpl,nvrt)
          call parallel_abort(errmsg)
        endif
        zrat=(ztmp(jlev)-zt)/(ztmp(jlev)-ztmp(jlev-1))
      endif

      !todo: assert
      !if(zrat<0.or.zrat>1) then
      !  write(errmsg,*)'QUICKSEARCH: Sigma coord. wrong (4):',jlev,zrat
      !  call parallel_abort(errmsg)
      !endif

!      if(kbpl==kz) then !in pure S region
!        ss=(1-zrat)*sigma(jlev-kz+1)+zrat*sigma(jlev-kz)
!      else
!        ss=-99
!      endif
!      if(ss<sigma(jlev-1).or.ss>sigma(jlev)) then
!        write(11,*)'Sigma coord. wrong (5):',jlev,ss,sigma(jlev-1),sigma(jlev)
!        stop
!      endif

!     Interpolate vel.
      if(indvel==-1) then !interpolated hvel using P_1^NC
        !Pure triangles only
!       Interpolate in vertical 
        do j=1,3 !sides and nodes
          nd=elnode(j,nnel)
          isd=elside(j,nnel)
!new37
          if(ics==1) then
            vxn(j)=su2(jlev,isd)*(1._rkind-zrat)+su2(jlev-1,isd)*zrat
            vyn(j)=sv2(jlev,isd)*(1._rkind-zrat)+sv2(jlev-1,isd)*zrat !side
          else !lat/lon
            call project_hvec(su2(jlev,isd),sv2(jlev,isd),sframe2(:,:,isd),frame0,uj,vj)
            call project_hvec(su2(jlev-1,isd),sv2(jlev-1,isd),sframe2(:,:,isd),frame0,uj1,vj1)
            vxn(j)=uj*(1-zrat)+uj1*zrat
            vyn(j)=vj*(1-zrat)+vj1*zrat
          endif !ics
          vzn(j)=ww2(jlev,nd)*(1._rkind-zrat)+ww2(jlev-1,nd)*zrat !node
        enddo !j=1,3

!       Interpolate in horizontal
        uuint=vxn(1)*(1._rkind-2._rkind*arco(1))+vxn(2)*(1._rkind-2._rkind*arco(2))+vxn(3)*(1._rkind-2._rkind*arco(3))
        vvint=vyn(1)*(1._rkind-2._rkind*arco(1))+vyn(2)*(1._rkind-2._rkind*arco(2))+vyn(3)*(1._rkind-2._rkind*arco(3))
        wwint=vzn(1)*arco(1)+vzn(2)*arco(2)+vzn(3)*arco(3)

      else !indvel>=0; interpolated hvel using P_1
!       No interpolate in time
        do j=1,i34(nnel) !nodes
          nd=elnode(j,nnel)
            do l=1,2 !levels
              lev=jlev+l-2
  !new37
              if(ics==1) then
                uu=uu2(lev,nd); vv=vv2(lev,nd)
              else !lat/lon
                call project_hvec(uu2(lev,nd),vv2(lev,nd),pframe(:,:,nd),frame0,uu,vv)
            endif !ics
            vxl(j,l)=uu !uu2(lev,nd) !+vis_coe*uf
            vyl(j,l)=vv !vv2(lev,nd) !vis_coe*vf
            vzl(j,l)=ww2(lev,nd)
          enddo !l
        enddo !j

!       Interpolate in vertical 
        do j=1,i34(nnel)
          vxn(j)=vxl(j,2)*(1._rkind-zrat)+vxl(j,1)*zrat
          vyn(j)=vyl(j,2)*(1._rkind-zrat)+vyl(j,1)*zrat
          vzn(j)=vzl(j,2)*(1._rkind-zrat)+vzl(j,1)*zrat
        enddo !j

!       Interpolate in horizontal
        uuint=dot_product(vxn(1:i34(nnel)),arco(1:i34(nnel))) !vxn(1)*arco(1)+vxn(2)*arco(2)+vxn(3)*arco(3)
        vvint=dot_product(vyn(1:i34(nnel)),arco(1:i34(nnel))) !vyn(1)*arco(1)+vyn(2)*arco(2)+vyn(3)*arco(3)
        wwint=dot_product(vzn(1:i34(nnel)),arco(1:i34(nnel))) !vzn(1)*arco(1)+vzn(2)*arco(2)+vzn(3)*arco(3)
      endif !indvel

      end subroutine quicksearch