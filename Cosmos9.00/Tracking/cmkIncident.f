!            make one incident particle, and make a copy of it
!       in IncientCopy.
!
      subroutine cmkIncident(incident, fin)
      implicit none
!
!     incident:  /track/.  output.   incident particl with track information.
!               copy of it is saved as 'IncidentCopy', and can be
!               inquired by call cqIncident(...)
!
#include  "Ztrack.h"
! #include  "Zmagfield.h"
#include  "Ztrackv.h"
#include  "Zmanagerp.h"
#include  "Zobs.h"
#include  "Zobsp.h"

      type(track)::incident
      integer fin  ! output. if 1, no more simulation
      type(coord)::angle

      integer icon
!

      icon = 1



      do while (icon .ne. 0)
!          sample energy, mass, code (mometum is not given)
         call csampPrimary(incident%p, fin)
!            DestEventNo is < 0; and cutoffed 1ry is counted, too
!           fin==1-->  all events generated
         if(fin .eq. 1) goto 10

!          If ObsPlane != spherical, fix angle at observation level
!          in detector system.
!          If ObsPlane = spherical, do the same tentatively.
!            (in this case Za1ry == 'is' or 'cos' guaranteed).

         call csPrimAng(angle)
         call cmkInc(incident, angle)
#if LABELING > 0
         incident%info = 0
         incident%label = 0
#endif
!           
         if(ObsPlane .eq. spherical) then
!              reset position and angle.

            call cresetPosAng(incident)
         endif    
         if(CutOffFile .ne. ' ') then
            call cifCutOff(icon)
         else
            icon =0
         endif
         if(icon .eq. 0) then
            if(Job .ne. 'newflesh') then
!
!                  for newflesh, next is managed by chookBgEvent
!
               EventsInTheRun = EventsInTheRun + 1
               if(Job .ne. 'flesh') then
                  EventNo = EventNo + 1
               endif
            endif
            call cupdtPrimC  ! update counter for each comp.
                             ! which is not rejected. 
         elseif( DestEventNo(2) .lt. 0 ) then
            if(Job .ne. 'newflesh') then
               EventsInTheRun = EventsInTheRun + 1
               if(Job .ne. 'flesh') then
                  EventNo = EventNo + 1
               endif
            endif
         endif
      enddo

 10   continue
      end
#if defined NEXT486
#define IMAG_P dimag
#elif defined PCLinux
#define IMAG_P dimag
#else
#define IMAG_P imag
#endif
!     *********************************
      subroutine cresetPosAng(incident)
      implicit none
!      After doing tentative business for energy, angle and
!      incident position for ObsPlnae != spherical,  this is
!      is used to reset incident position and angle for
!      spherical case.
!     
!      The incident position is uniform around a point given by the
!      (Latit, Longit, HeightOfInj)=PolarInjPos.  It will be  distributed
!      within  the half opnenig angle range given by Azimuth.
!      (if it is > 180, regarded as 180.  Hence, if Azimuth=(45,90),
!      incident position will be a ring like region on  a
!      sphere.)
!      As to the zenith angle at the incident position, 
!      it is determined isotropically from CosZenith.  Azimuth is
!      not used for this purpose.  Therefore, if zenith angle >= 90,
!      we will discard such particles. Hence, Imag(CosZenith) 
!      must be > 0.
!      
!
#include  "Zglobalc.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"      
#include  "Ztrackp.h"
#include  "Zobs.h"  
#include  "Zobsv.h"  
#include  "Zincidentp.h"

!
      type(track)::incident  ! input/output.
      type(coord)::incipos
      type(coord)::dir1
      type(coord)::dir2
      logical first
      real*8 len, cs, sn, sint, u, oa1, oa2

!     @@@@@@@@@@@@@For  bug correction to ObsPlane=3
      real*8 cosx, ux
!     @@@@@@@@@@@@@@

      data first/.true./
      save first
      
!         fix position 
      oa1 = real(Azimuth)
      oa2 = IMAG_P(Azimuth) 
!        if we don't use oa1, oa2 but use real(..) directly in
!        the subroutine call, Absoft compiler give always 0 for
!        oa1  
!
      cosx = 2.      ! @@@@@@@@@
      ux = 0.        ! @@@@@@@@@
      do while (cosx .gt. -ux)  !@@@@@@@@@
         if(first) then
            call cuonSphere(1, PolarInjPos%r(3), PolarInjPos%r(1),
     *           PolarInjPos%r(2), oa1, oa2,  incident%pos%xyz)
            first =.false.
         else
!               this is quicker a bit
            call cuonSphere(2, PolarInjPos%r(3), PolarInjPos%r(1),
     *      PolarInjPos%r(2),  oa1, oa2, incident%pos%xyz)
         endif
!           fix angle around zenith at incident.pos.xyz
!            convert PolarInjPos to xyz vector
         call ctransCoord2('xyz', incident%pos%xyz, incipos)
!            convert its direction to direction cos
         call c3DV2DDCos(incipos, dir1, len)
!            sample angle around dir1 (x,y axes arbitrary)
         call rndc(u)
         if(Za1ry .eq. 'is') then
            dir2%r(3) =  -(  (IMAG_P(CosZenith)- real(CosZenith) )*u +
     *             real(CosZenith) )    ! going down is negative
         elseif(Za1ry .eq. 'cos') then
!                  cos dcos
            call ksampLin(1.0d0, 0.0d0, 
     *       real(CosZenith), IMAG_P(CosZenith), dir2%r(3))
            dir2%r(3) = -dir2%r(3) ! going down is negative
         else
            call cerrorMsg('Za1ry =', 1)
            call cerrorMsg(Za1ry, 1)
            call cerrorMsg(
     *      ' for ObsPlane=3 is invalid; must be: "is" or "cos"',0)
         endif
         call rndc(u)
         call kcossn(cs, sn)
         sint = sqrt(1.d0-dir2%r(3)**2)
         dir2%r(1) = - sint*cs  ! - is needed for going down ptcl.
         dir2%r(2) = - sint*sn
!         convert it to xyz system
         call ctransVectZ(dir1, dir2, incident%vec%w)
!     @@@@@@@@@@@@@@@
         call cscalerProd(incident%vec%w, dir1, cosx)
         call rndc(ux)
      enddo   
!     @@@@@@@@@@@@@@@
!                uv 5.51
      if(Reverse .ne. 0) then
!          we must revert the angle
         dir2%r(1) = -dir2%r(1)
         dir2%r(2) = -dir2%r(2)
         dir2%r(3) = -dir2%r(3)
         call ctransVectZ(dir1, dir2, incident%vec%w)
      endif

!         reset others  
      call cresetPrim2(incident)
      end
!     *******************************
      subroutine cresetPrim(incidentp, angle)
      implicit none
!        reset primary.  This is typically used by the user,
!        at chookBgEvent to reset the primary which has been
!        set by the sytem so that the user can set own primray.
!
#include  "Ztrack.h"  
      type(ptcl)::incidentp ! input. must have E, mass, charge,subcode
      type(coord)::angle  ! input. direction cos at 'det' system
!
      type(track)::inc2
#if defined (KEKA) || defined (KEKB)
!  for IBM, we must  write as follows.
      inc2%p%charge=incidentp%charge
      inc2%p%subcode=incidentp%subcode
      inc2%p%fm%p(4) = incidentp%fm%p(4)
      inc2%p%mass = incidentp%mass
#else
      inc2%p = incidentp 
#endif
      call cmkInc(inc2, angle)
      call ciniTracking( inc2 )
      call cinitStack
      call cpush(inc2)
      end
!     *******************************
      subroutine cresetPrim2(incident)
      implicit none
!        reset primary.  This is typically used by the user,
!        at chookBgEvent to reset the primary which has been
!        set by the sytem so that the user can set own primray.
!     The difference from cresetPrim is that the parameter is
!     track, and  this is for
!     ObsPlane==spherical case where you can put very arbitrary
!     incdint injection point.
!          See cmkInc2 for what you must set for incident.
!
#include  "Ztrack.h"  
!
      type(track)::incident  ! input. you must give everything
!                                       about primary
      
      call cmkInc2(incident)
      call ciniTracking( incident ) 
      call cinitStack
      call cpush(incident)
      end
      
!     ********************************
      subroutine cmkInc(incident, angle)
      use modAtmosDef
      implicit none
!
#include  "Zglobalc.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"      
#include  "Ztrackv.h"
#include  "Ztrackp.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Zobsv.h"
#include  "Zincidentp.h"
#include  "Zincidentv.h"
! #include  "Zearth.h"
#include  "Zcode.h"
      type(track)::incident ! input/outut.
!             must have E, m, code, subcode, charge
      type(coord)::angle    ! input. must have direction cos in the det. sys.
      type(coord)::xyz
!     
!
      real*8 to100km, clenbetween2h, zenithcos, len
      integer i, j

      AngleAtObsCopy = angle
      if(Reverse .ne. 0 ) then
!           to see  cut off or to see if go out of Earth
         do i = 1, 3
            AngleAtObsCopy%r(i) =  - AngleAtObsCopy%r(i) 
         enddo
!             charge conjugate
         incident%p%charge = -incident%p%charge
         if(incident%p%code .ne. kgnuc) then
            incident%p%subcode = -incident%p%subcode
         endif
      endif
!           convert it to 'xyz' system
      call ctransVectZx(1, DetZaxis, DetXaxis, AngleAtObsCopy, 
     *    DcAtObsXyz)
!
      incident%pos%xyz%sys = 'xyz'  !  Exyz system 
!      --- fix injection point ---
!      get length from the deepest obs. place + OffsetHeight to 
!          HeightOfInj (=100 km); Normally OffsetHeight is 0.
!         if the primary is to be directed to different height
!         than the detector, make it non zero.
      zenithcos = - AngleAtObsCopy%r(3)
      if(ObsPlane .ne. spherical) then
         to100km = clenbetween2h(
     *    ObsSites(NoOfSites)%pos%radiallen + OffsetHeight, 
     *    Eradius + HeightOfInj,  zenithcos )  ! we need zenith angle here
      else
!           dummy 
         to100km = 10000.
      endif 

!         primary is going upward even if Reverse = 0.
      Upgoing = Reverse .eq. 0 .and. zenithcos .lt. 0. 
     *   .and. HeightOfInj .lt. ObsSites(NoOfSites)%pos%height

      if(( Reverse .eq. 0 .and. zenithcos .lt. 0. 
     *  .and.  HeightOfInj .gt. ObsSites(NoOfSites)%pos%height)
     *  .or. ( Reverse .ne. 0 .and. zenithcos .gt. 0.)) then
         if(ObsPlane .ne. spherical) then
!            distance to the conjugate point
            to100km = to100km -
     *       2*(ObsSites(NoOfSites)%pos%radiallen + OffsetHeight)*
     *       zenithcos
!             we should go reversed direction
            to100km = - to100km
         endif
      endif
!           injection point
      do i = 1, 3
         incident%pos%xyz%r(i) = ObsSites(NoOfSites)%pos%xyz%r(i) 
     *    + Offset%r(i)  + to100km * DcAtObsXyz%r(i)
      enddo
      call csetPos(incident%pos)
      call csetDirCos(DcAtObsXyz, incident)   ! set dc and coszenith in incident
!           momentum business
      call ce2p(incident)    !  px, py, pz   from direction cos

!          set time 0
      incident%t = 0.
      incident%wgt = 1.       ! weight is 1.


      do i = 1, NoOfSites
!               correction for Perpendicular : 2004/07/19
         if( ObsPlane .eq. Perpendicular ) then
!            fixing incident.where later.
         elseif( ObsPlane .ne. NotUsed ) then
            if( HeightOfInj .gt. ObsSites(i)%pos%height ) then
               incident%where = i
               goto 222
            endif
         endif
      enddo
      if(HeightOfInj .lt. BorderHeightL) then
         call cerrorMsg('Injection height is < BorderHeightL',0)
      endif
      incident%where = NoOfSites + 1
 222  continue
      if(HeightOfInj .gt. BorderHeightH) then
         call cerrorMsg('Injection height is > BorderHeightH',0)
      endif


      incident%asflag = 0
      incident%pos%colheight = Infty  ! latest nuc. collision height
      IncidentCopy = incident
!           shift the origin of detectors to be on the primary track
!           if OffsetHight=0
      if(OffsetHeight .eq. 0. .and. ObsPlane .ne. spherical) then
         if(zenithcos .ge. 0. .or. Upgoing) then
            do i = 1, NoOfSites-1
               len = clenbetween2h(
     *              ObsSites(NoOfSites)%pos%radiallen, 
     *              ObsSites(i)%pos%radiallen,
     *              zenithcos ) 

               do j = 1, 3
                  ObsSites(i)%pos%xyz%r(j) = 
     *            ObsSites(NoOfSites)%pos%xyz%r(j)
     *             + len * DcAtObsXyz%r(j)

               enddo
            enddo
            do i = 1, NoOfASSites-1
               len = clenbetween2h(
     *         ASObsSites(NoOfASSites)%pos%radiallen, 
     *         ASObsSites(i)%pos%radiallen,
     *         zenithcos ) 
               do j =1 , 3
                  ASObsSites(i)%pos%xyz%r(j) = 
     *                ASObsSites(NoOfSites)%pos%xyz%r(j)
     *                + len * DcAtObsXyz%r(j)
               enddo
            enddo
         endif
      endif
!          compute min time from  the injection point to each
!         obs level.
      if(ObsPlane .ne. spherical) then
         call csetMinTime(incident)
         if(HeightOfInj .lt. BorderHeightL) then
            call cerrorMsg('Injection height is < BorderHeightL',0)
         endif
      endif
      end
      
!     ****************************
      subroutine cmkInc2(incident)
!          this may be used when incident is ready (
!       even when   ObsPlane==spherical)
!      incident must have:
!          incident.p:   code, subcode, mass, energy
!          incident.pos: xyz.r(1), xyz.r(2), xyz.r(3) in E-xyz
!          incident.vec: w.r(1), w.r(2), w.r(3)   in E-xyz
!          incident.where:  from which height the incident particle
!                           crosses ?
!          Otheres are set here.
      use modAtmosdef
      implicit none
!
#include  "Zglobalc.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"      
#include  "Ztrackv.h"
#include  "Ztrackp.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Zobsv.h"
#include  "Zincidentp.h"
#include  "Zincidentv.h"
!  #include  "Zearth.h"

      type(track)::incident ! input
!
!
      incident%pos%xyz%sys = 'xyz'  !  Exyz system 
!         
      Upgoing = .false.

      call csetPos(incident%pos)
      call cgetZenith(incident, incident%vec%coszenith)
!      call csetDirCos(DcAtObsXyz, incident)   ! above is o.k 
!           momentum business
      call ce2p(incident)    !  px, py, pz   from direction cos

!          set time 0
      incident%t = 0.
      incident%wgt = 1.       ! weight is 1.


      incident%asflag = 0
      incident%pos%colheight = Infty  ! latest nuc. collision height
      IncidentCopy = incident
      end
!     *****************************
!        compute the minimum time the light needs to reach
!        each observation level from a given radial height from
!        along the primary direction
      subroutine csetMinTime(from)
!      use modAtmosDef   ! not needed ?
      implicit none
!
!
#include  "Ztrack.h"
! #include  "Zmagfield.h"      
#include  "Ztrackv.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Zobsv.h"
#include  "Zincidentp.h"
#include  "Zincidentv.h"
!#include  "Zearth.h"
      
      type(track)::from   ! input.  track to be origin
 
      real*8 leng
!        clenbetween2h, leng
      integer i, icon

      do i = 1, NoOfSites
!         ObsSites(i).minitime =
!     *       clenbetween2h(from.pos.radiallen, 
!     *       ObsSites(i).pos.radiallen,
!     *       from.vec.coszenith)    ! actually this is in m.
         call clenbetw2h(from%pos%radiallen, 
     *       ObsSites(i)%pos%radiallen, from%vec%coszenith,
     *       leng, icon)     ! actually leng is in m.
         if(icon .eq. 0) then
            ObsSites(i)%minitime = leng
         else
!              icon !=0 ==> light cannot come with this angle
            ObsSites(i)%minitime = 1.d10
         endif
      enddo

      from%t = 0.    ! reset time of incident track.
      end
!     *************************** inquire incident particle
      subroutine cqIncident(incident, AngleAtObs)
!     *************************** 
      implicit none
#include  "Ztrack.h"
#include  "Zincidentv.h"
      type(track)::incident
      type(coord)::AngleAtObs
      incident = IncidentCopy
      AngleAtObs = AngleAtObsCopy
      end
!     **********************************
      subroutine cuonSphere(ini, rin, teta, phi, oa1, oa2, pos)
      implicit none
#include "Zglobalc.h"
#include "Zptcl.h"
#include "Zcoord.h"
!         This is a modified version of epuonSphere in Epics
!         generate a random point uniformly distributed on the
!         surface of a sphere.  Points are distributed around
!         given polar angles (teta, phi) within a given opening angle
!        (oa1~oa2). 
!    By uniform  is meant that the points are uniformly distributed on
!    the surface of the sphere but not on a projected plane.
!
      integer ini     ! input
                      !  1-->  teta and phi are different from
                      !        previous call or this is the first call.
                      !  != 1 -->  teta, and phi are the same as
                      !        the previous call.
      real*8  rin           ! input.  radius of the sphere
      real*8  teta           ! input.  polar angle in degree
      real*8  phi            ! input.  azimutal angle in degree
      real*8  oa1            ! input.  starting half opnening angle in degree
      real*8  oa2            ! input.  ending half opnening angle in degree
      type(coord)::pos     ! output. an  obtained random point in Exyz

      type(fmom)::xyz
      type(fmom)::xyz2
      real*8  a(4, 4), b(4, 4), ba(4, 4)
      real*8  u, r
      real*8 fcos,  fsin
      save ba

      r = rin *0.999999999d0
      if(ini .eq. 1) then
         call cgetRotMat4(2, -teta*Torad, a)
         call cgetRotMat4(3, -phi*Torad, b)
         call cmultRotMat4(b, a, ba)
      endif

      call rndc(u)
      fcos = cos( min(oa2,180.d0) * Torad)
      fcos = ( cos( min(oa1,180.d0)*Torad ) - fcos) * u +  fcos
      fsin = sqrt(1.d0- fcos**2)
      call rndc(u)
      u = u*pi*2
      xyz%p(1) = r * (fsin * cos(u))
      xyz%p(2) = r * (fsin * sin(u))
      xyz%p(3) = r * fcos

      xyz%p(4) =  1. ! dummy
      call capplyRot4(ba, xyz, xyz2)
      pos%r(1) = xyz2%p(1)
      pos%r(2) = xyz2%p(2)
      pos%r(3) = xyz2%p(3)


      pos%sys = 'xyz'
      end
!       ************************* see if geomagnetic cut or not.
        subroutine cifCutOff(icon)
        implicit none
#include "Zglobalc.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zcode.h"
#include "Ztrack.h"
#include "ZrigCut.h"

         integer icon   !   output. 0 ==> not cut. 1 ==> cut.

         type(coord)::angleatOb
         type(track)::inc

         real*8 p, rig, zen, azm, theta,  prob, u

         call cqIncident(inc, angleatOb)

         if(inc%p%charge .eq. 0) then
            icon = 0
         else
            if(Rdatafmt .le. 4) then
!              angleatOb is down going ptcl's one, change sign
               angleatOb%r(1) = - angleatOb%r(1)
               angleatOb%r(2) = - angleatOb%r(2)
               angleatOb%r(3) = - angleatOb%r(3)
!                 convert to theta fai in deg
               call kdir2deg(angleatOb%r(1), angleatOb%r(2),
     *          angleatOb%r(3), theta, azm)
!
!                 azm is given from the current x-axis  (+ is counter
!                 clock wise) The x-axis is XaxisFromSouth
!                 degrees from the south in counter clockwise.
!                 convert azm so that measured from the south
!
               azm = mod(azm+ XaxisFromSouth, 360.d0)
               if(ZenValue .eq. 'cos') then
!                    table zenith is in cos
                  zen = angleatOb%r(3)
               else
                  zen = theta
               endif
            elseif(Rdatafmt .eq. 5) then
!                 in this case, azm is longitude; zen is latitude or cos(lati)
!       
!               zen =  atan2( inc.pos.xyz.z, 
!     *                sqrt(inc.pos.xyz.x**2+inc.pos.xyz.y**2) )
               zen =  atan2( inc%pos%xyz%r(3),
     *                sqrt(inc%pos%xyz%r(1)**2+inc%pos%xyz%r(2)**2) )
               if(ZenValue .eq. 'cos') then
                  zen = cos( zen-pi/2.0d0 )
               else
                  zen = zen*Todeg
               endif
!               azm = atan2(inc.pos.xyz.y, inc.pos.xyz.x)*Todeg
               azm = atan2(inc%pos%xyz%r(2), inc%pos%xyz%r(1))*Todeg
            else
               call cerrorMsg('dataformat for cut off invalid',0)
            endif

            p = sqrt(inc%p%fm%p(4)**2 - inc%p%mass**2)
            rig = p/abs(inc%p%charge)
            call crigCut(azm, zen, rig, prob)
            call rndc(u)
            if(u .lt. prob) then
               icon = 0
            else
               icon = 1
            endif
         endif
         end
