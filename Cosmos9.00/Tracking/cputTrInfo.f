!      ready made tracer
!          put trace info.
!     The output x,y,z values will be those defined as follows
!      depending on the  value of "Trace" 
!
!         0       This routine is not called.
!         
!           
!        <10       x, y, z in the primary sytem
!        <20       x, y, in the primary system. z in kg/m2
!
!        <30       x, y, z in the detector system
!        <40       x, y, in the detector system. z in kg/m2
!     
!        <50       x, y, z in 'xyz' system
!        <60       x, y in 'xyz' , z in kg/m2  
!
!    Primary system:  origin is the deepest detector.
!                     Z-axis is the primary direction
!                     X-aixs is Z x  Vertical axis
!                     X-Y plane is orthgonal to the primary
!    Detector system: origin is the deepest detector.
!                     Z-axis is the vertical one
!                     X-axis is directed to the magnetic east.
!                     X-Y palne is horizontal.
!    z in kg/m2 :     Vertical depth in kg/m2  above the 
!                     deepest detector to the particle.
!
!      60< < 100    ready made Cerenkov output.
!      100< <160      chookTrace is called.
!      160< <200      chookCeren is called with Cerenkov mode.
!        
!      
      subroutine cputTrInfo
      implicit none

#include  "Ztrack.h"
#include  "Ztrackp.h"
!  #include  "Ztrackv.h"

      
      type(coord)::f
      type(coord)::t

      if(Trace .ge. 100 .and. Trace .lt. 160) then
         call chookTrace
         return    !      ******************
      endif
      if(Trace .gt. 60 .or. Trace .gt. 160) then
         call cputCerenkov   ! Cerenkov output
         return    !      *****************
      else
!               convert coord.
         call ccoordForTr(Trace, f, t)
         call cwrtTrInfo(f, t)
      endif
      end
      subroutine cwrtTrInfo(f, t)
      implicit none

#include  "Ztrack.h"
! #include  "Zmagfield.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zobs.h"
#include  "Zobsv.h"
      
      type(coord)::f
      type(coord)::t

      real*8  xxx/-1.d37/, yyy/-1.d36/, zzz/1.d34/
      integer kkk/-1000/, chg/-1000/
      save xxx, yyy, zzz, kkk, chg

      if(kkk .ne. MovedTrack%p%code .or. f%r(1) .ne. xxx
     *    .or. f%r(2) .ne. yyy .or. f%r(3) .ne. zzz  .or.
     *     chg .ne. MovedTrack%p%charge ) then
         if(xxx .ne. -1.d37) then
            write(TraceDev, *) 
            write(TraceDev, *) 
         endif
         if(TimeStructure) then
!!            write(TraceDev, '(3g16.8, i4, g11.4, i4, g16.8)')
            write(TraceDev,
     *      '(1p,3g24.16,0p, i4,1p, g24.16,0p, i4, 1p, g24.16)')
     *        f%r(1), f%r(2), f%r(3), TrackBefMove%p%code,
     *        TrackBefMove%p%fm%p(4)-TrackBefMove%p%mass,
     *        TrackBefMove%p%charge,  TrackBefMove%t
!!!     *         TrackBefMove%p%charge,  TrackBefMove%pos%height
         else
            write(TraceDev, '(3g16.8, i4, g11.4, i4)')
     *        f%r(1), f%r(2), f%r(3), TrackBefMove%p%code,
     *        TrackBefMove%p%fm%p(4)-TrackBefMove%p%mass,
     *        TrackBefMove%p%charge
         endif
      endif
      if(TimeStructure) then
!!!         write(TraceDev,  '(3g16.8, i4, g11.4, i4, g16.8)')
         write(TraceDev,
     *      '(1p,3g24.16,0p, i4,1p, g24.16,0p, i4, 1p, g24.16)')         
     *      t%r(1), t%r(2), t%r(3), MovedTrack%p%code,
     *      MovedTrack%p%fm%p(4)-MovedTrack%p%mass,
     *      MovedTrack%p%charge,  MovedTrack%t
!!!     *      MovedTrack%p%charge,  MovedTrack%pos%height
      else
         write(TraceDev,  '(3g16.8, i4, g11.4, i4)')
     *      t%r(1), t%r(2), t%r(3), MovedTrack%p%code,
     *      MovedTrack%p%fm%p(4)-MovedTrack%p%mass,
     *      MovedTrack%p%charge
      endif
      xxx = t%r(1)
      yyy = t%r(2)
      zzz = t%r(3)
      kkk = MovedTrack%p%code
      chg = MovedTrack%p%charge
      end
!     ***************************************
!             convert coord for tracing.  This is called by
!          cputTrInfo and cputCerekov

      subroutine ccoordForTr(how, f, t)

      implicit none
!     ***************************************
#include  "Ztrack.h"
! #include  "Zmagfield.h"      
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zobs.h"
#include  "Zobsv.h"
      
      integer how   ! input. see Trace in cputTrInfo.
      type(coord)::f
      type(coord)::t  ! output. transformed 'from to' coord 

      if(how .le. 20) then
!              to 1ry system
         call cxyz2prim(ObsSites(NoOfSites)%pos%xyz,
     *                  TrackBefMove%pos%xyz, f)
         call cxyz2prim(ObsSites(NoOfSites)%pos%xyz,
     *                 MovedTrack%pos%xyz, t)
         if(how .ge. 10) then
            f%r(3) = ObsSites(NoOfSites)%pos%depth 
     *            - TrackBefMove%pos%depth 
            t%r(3) = ObsSites(NoOfSites)%pos%depth 
     *            - MovedTrack%pos%depth 
         endif
      elseif(how .le. 40) then
!              to detector system
         call cxyz2det(ObsSites(NoOfSites)%pos%xyz, 
     *                 TrackBefMove%pos%xyz, f)
         call cxyz2det(ObsSites(NoOfSites)%pos%xyz, 
     *                 MovedTrack%pos%xyz, t)
         if(how .ge. 30) then
            f%r(3) = ObsSites(NoOfSites)%pos%depth 
     *            - TrackBefMove%pos%depth 
            t%r(3) = ObsSites(NoOfSites)%pos%depth 
     *            - MovedTrack%pos%depth 
         endif
      elseif(how .le. 60) then
         f = TrackBefMove%pos%xyz
         t = MovedTrack%pos%xyz
         if(how .ge. 50) then
            f%r(3) = ObsSites(NoOfSites)%pos%depth 
     *            - TrackBefMove%pos%depth 
            t%r(3) = ObsSites(NoOfSites)%pos%depth 
     *            - MovedTrack%pos%depth 
         endif
      else
         call cerrorMsg('100>= how > 60 for ccoordForTr ', 0)
      endif
      end
