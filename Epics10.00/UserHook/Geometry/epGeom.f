!         this may be used from ephook.f
!         incident is assumed to be a neutrio
!      Suupose a newtrino passes thru  many componests
!      we see it runs in many different media. The path
!      is devided into many segments, each of which consists
!      of the same medium.  If medium changes or component
!      number changes during the path, segment becomes different.
!  make next defined by  #define STOPWATCH
#undef STOPWATCH
!  to measure the time for a primary to pass thru each component.
!  Usually, its small and 0.00000  sec is obtained. In some 
!  case, when the primary incident position is exactly on the
!  edge  or surface of a component, Epics might fall in some
!  kind of a basket case and might need order of few sec to even
!  min to escape that status.  One known example case could be 
!  generated when a particle is on such a point and the next
!  crossing point is made to be on the same point.
!  In such a case, the particle is made to go forward by
!  2*EpsLeng and go back -EpsLeng, and in efect 
!  go forward by EpsLeng (10^-6 cm in defaut). 
!  By making STOPWATCH defined, we can detect such component.
!  The output will be on error out.
!  
! Workaround: for +primary case.
!   1)  Instant workaroud is to shift the incident position
!       by 10^-5 cm or so.  (Use shiftxy.awk to change xy
!       in +primary and get new file. )
!     or
!       Use a sepicsfile  parameter, ShiftInciPos. E.G
!    a)   ShiftInciPos "1.e-5 "
!    b)   ShiftInciPos "1.e-5  xy"
!    c)   ShiftInciPos "2.e-5  yz"
!    d)   ShiftInciPos "1.e-5  zx"
!         a and b are the same. 'xy' means input position is
!         on the x-y plane, 'yz' on y-z and 'zx' on z-x plane
!         The world could be sphere but box is better.
!
!   2)  Fundamental fix is to find a cause such as mentioned
!       above, and repair the boundary search program for the
!       component.
!
      module  epModGeometry
      implicit none

      real(4),allocatable:: MatLen(:)  !cm-length of each segment
      real(4),allocatable:: MatGrm(:)  ! its g/cm2 equivalent
      real(4),allocatable:: MatRL(:)  ! its r.l equivalent
      real(4),allocatable:: MatPos1(:,:) ! (1:3,:) is coord of
               !  segment starting  point. each segment
      real(4),allocatable:: MatPos2(:,:)  ! ending point
      integer,allocatable:: MatCn(:)   !  component # of the segment
      character(12),allocatable:: MatName(:) ! media name
      real(8),save:: InciDir(3), InciPos(3)  ! Incident direction
                         ! cosine, and pos.  Always world coord.
      real(8),save:: tempPos(3)  ! working array to save position
      real(8),save:: ZeroLeng   ! EpsLeng(in epicsfile) * 1.01 
      integer:: NoOfSegments  ! # of segments for a neutrino passed
      integer,save:: SuppressBlank=0 ! suppress blank in the output(=1) or not(0)
      integer,save:: ToWorld = 1  ! coor. is written in world
                                   ! coordinate. if 0, local. 
#if defined (STOPWATCH)
      real(4):: start, goal, dt
#endif

      contains

      subroutine epGeom0
!          set countIO flag for all components.
!          To be called from uiaev
       implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepTrackp.h"
       integer n, m
       n = Det%nct
!
       write(0,*) ' # of comps=',n
!
       Det%cmp(1:n)%CountIO = 3  ! set boundary count flag for oall comp. 
       m = n*2  ! m could normally be much smaller if n is very large
       allocate(  MatLen(m) )
       allocate(  MatGrm(m) )
       allocate(  MatRL(m) )
       allocate(  MatPos1(3,m) )
       allocate(  MatPos2(3,m) )
       allocate(  MatName(m) )
       allocate(  MatCn(m) )
       ZeroLeng = EpsLeng*1.01d0

       call epqHooki(1, SuppressBlank)   ! read control parameters
       call epqHooki(2, ToWorld)

       end  subroutine epGeom0

      function epGeomPos(trk) result(pos)
!           track pos (is converted world, if
!           ToWorld == 1 and) is stored in pos 
!           as a function value. 
      implicit none
#include  "ZepTrack.h"
       type(epTrack)::  trk

      real(8):: pos(3)


      if( ToWorld == 1 ) then
         call epl2w(trk%cn, trk%pos, pos)
!/////////////
!          write(0,...) cannot be used, since
!          epGeomPos is used in write(0,*) epGeomPpos
!          and results is  recursive use of '0' device.
!         write(*,*) 'l pos=',
!     *   trk.pos.x, trk.pos.y, trk.pos.z
!         write(*,*) ' pos=',pos(:)
!///////////////
      else
!         since union/map is disabled now  in ZepPos.h  (why?)
!         we cannot use pos(:) = aTrack.pos.r(:)

         pos(1) = trk%pos%x
         pos(2) = trk%pos%y
         pos(3) = trk%pos%z
      endif
      end  function epGeomPos

      subroutine epGeom1
      implicit none
      !       init for 1 event; to be called from uafi1ev
#include  "ZepTrack.h"
       type(epTrack)::  incident
      NoOfSegments = 0

      call epqinc(incident) !  get incident info. local coord
      tempPos(:) = epGeomPos(incident)    
      call epl2wd( incident%cn, incident%w, InciDir)
      call epl2w( incident%cn, incident%pos, InciPos)

#if defined (STOPWATCH)
      call cpu_time(start)
#endif
      end      subroutine epGeom1


      subroutine epGeomB(info, aTrack, Move, media)
!         boundary  business; to be called from
!           userbd 
       implicit none
#include "ZepTrack.h"
#include "Zmove.h"
#include "Zmedia.h"
#include "Zcode.h"     
        integer,intent(in)::info 
             ! 0--> ptcl is exiting to void
             ! <0 --> ptcl is exiting to |info| comp.
             ! >0 --> ptcl is entering from info comp.
       type(epTrack):: aTrack 
             ! input. Current track before it is moved.
             !  If info <=0, the Move.Track's  position is 
             !  somewhere inside the component from which
             !  the track is existing.
             !  If info >0, aTrack.pos is the position 
             !  just before exiting the prvious component
       type(epmove)::  Move
             ! input/output. 
             ! Move.Track is the track
             ! infomation of the current particle
             ! moved to a new position.
             ! Say, Move.Track.cn is the  current
             ! comp. number.  For other details,
             ! see userde.
!
!                      comp.1    comp.2
!          info<0    |  *-----x|          |        
!                                  * is aTrack.pos 
!                             for nu * and x is same
!                                  x is Move.Track.pos
!          info>0    |        *|x         |       
!                          x is saved to tempPos.
       type(epmedia):: media   ! input.

      real(8):: dist
!///////////
       type(epPos)::  posw
!      write(0,*) ' info=',info
!!                next epGeomPos is in the write statement
!!                so if epGeomPos has write statement with 0
!!                dev#   inside, it will make recursive error
!!      write(0,'(a, 1p,3g15.5)') ' aTrack.pos=',epGeomPos(aTrack)
!!      write(0,'(a, 1p,3g15.5)')
!!     *  ' Move.pos=',epGeomPos(Move.Track)
!!      write(0,*) ' media=', media.name
!//////////////////
        if( info <= 0 ) then      ! exiting to void or other comp.
           dist =
     *       sqrt( sum( (epGeomPos(Move%Track) - tempPos(:))**2 ) )

#if defined (STOPWATCH)
           call cpu_time(goal)
           dt = goal - start
           call epl2w(Move%Track%cn, Move%Track%pos, posw)
           write(0,'(f9.3,1p,g14.5,1x, a8, i5, 6g14.5, i4 ) ')
     *      dt, dist, media%name, aTrack%cn, 
     *      Move%Track%pos%x, Move%Track%pos%y, Move%Track%pos%z
           write(0,'(37x, 1p, 3g14.5, i3)')
     *       posw%x,  posw%y,  posw%z,  NoOfSegments
#endif

         if( dist <= ZeroLeng ) then
                 !      <=  EpsLeng  ===> 0
         else
            NoOfSegments =  NoOfSegments  + 1

            MatName(NoOfSegments) = media%name
            MatCn(NoOfSegments) = aTrack%cn

            MatLen(NoOfSegments) = dist
            if( media%name == 'sp' .or.
     *           media%name == 'hollow' ) then
               MatGrm(NoOfSegments) = 0.
               MatRL(NoOfSegments) = 0.
            else
               MatGrm(NoOfSegments) = dist*media%rho * media%rhoc
               MatRL(NoOfSegments) =
     *          dist* media%rhoc/media%X0g/media%gtocm
            endif
            MatPos1(:,NoOfSegments) =tempPos(:)
            MatPos2(:,NoOfSegments) = epGeomPos(Move%Track)
         endif
         tempPos(:) = epGeomPos(Move%Track)

#if defined (STOPWATCH)
         start = goal
#endif

      else
         tempPos(:) = epGeomPos( Move%Track)

#if defined (STOPWATCH)
         call cpu_time(start)
#endif

      endif
      end      subroutine epGeomB


      subroutine epGeomE1ev
      ! to be called when 1 event ended from  ue1ev
      implicit none
      integer:: i, nev,icon
      integer:: nc  ! # of characters containe 'line' below
      character(200):: line
      real(8):: sumrl
      real(8):: opos(3)

      call  epqevn(nev)         ! # of events
      call  epProjPos(opos,icon)
      if(icon == 0)  then
         write(line,'(a,1p,3g15.5)') " projp ", opos(:)
         if( SuppressBlank == 1) then
            call ksupblank(line, nc) ! supress blank
         endif
         write(*, '(a)') trim(line)
      endif
      write(line,'(a, i7, i4, 1p,3g14.5, 3g17.9)')
     *  "h ", nev, NoOfSegments, InciPos(:), InciDir(:)
      if( SuppressBlank == 1) then
         call ksupblank(line, nc) ! supress blank
      endif
      write(*, '(a)') trim(line)
      sumrl = 0.
      do i = 1, NoOfSegments
!!      do i = 1, NoOfSegments ! ,10
         sumrl = sumrl +  MatRL(i)
         write(line,'(i3, i6, " ",  a8, 1p, 10g14.6)')
     *    i, MatCn(i), MatName(i), MatPos1(:,i), MatPos2(:,i),
     *    MatLen(i), MatGrm(i), MatRL(i), sumrl
         if(SuppressBlank == 1) then
            call ksupblank(line, nc)
         endif
         write(*, '(a)') trim(line)
      enddo
      write(line,'(i3, i6, " ",  a8, 1p, 10g14.6)')
     *  0, 0, "notmedia", 0., 0., 0., 0., 0., 0., 0., 0.,0., sumrl
      end      subroutine epGeomE1ev

      subroutine epProjPos(pos,icon)
      implicit none
#include  "ZepDirec.h"
#include  "ZepPos.h"      
#include "Zsparm.h"
!     type(epPos),intent(out):: pos  ! see below
      real(8),intent(out):: pos(3) ! see below
      integer,intent(out):: icon ! 0.  pos obtained
           ! InputP='u->gsp2' and InputA='fix"
           ! otherwise  1.  pos not obtained.

      if( InputP == "u->sph2" .and. InputA == "fix" ) then
!     get projected pos.  x-y coordinate  is on the surface
!     of projection plane. where points are uniform. Z is
!     (x,y,z) is on a spherical  surface.
         call epgonSphere2(pos)
         icon = 0
      else
         icon = 1
      endif
      end
      end module  epModGeometry            

