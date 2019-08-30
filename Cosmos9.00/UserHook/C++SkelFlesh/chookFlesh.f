#include "../cmain.f"
#include "chookHybAS.f"
#include "../ctemplCeren.f"
!
!  *************************************** hook for Beginning of a Run
!  * At this moment, all (system-level) initialization for this run
!  * has been ended.  After this routine is executed, the system goes into the
!  * event creation loop.
!  *
      subroutine chookBgRun
      implicit none
#include "Zmanagerp.h"
#include "Zprivate.h"

      real*8  temp
      character*100 msg
      integer klena
!     ==================================================
      
      integer seed(2)
!     ==================================================

      EventNo = 0
      RealBegin = .true.
      TopOfNode = .true.

!            namelist output
      call cwriteParam(ErrorOut, 0)
!            primary information
      call cprintPrim(ErrorOut)
!            observation level information
      call cprintObs(ErrorOut)

      call cqUHooki(1, Mdev)      ! get skeleton memo dev #
      call cqUHookc(1, msg)       ! get file name for sekelton memo
      call cgetfname(msg, Mskel)  ! add host name etc if needed

      open(Mdev, file=Mskel(1:klena(Mskel)), form='unformatted',
     *  status='old' )

      end

!     *********************************** hook for Beginning of  1 event
!     *  All system-level initialization for 1 event generation has been
!     *  eneded at this moment.
!     *  After this is executed, event generation starts.
!     *
      subroutine chookBgEvent
      implicit none
#include "Zprivate.h"


      integer nomore
      if( RealBegin ) then
         call cbegin1ev( nomore )
         if( nomore .eq. 1) then
            call cerrorMsg('all events are fleshed', 1)
            stop  !!!!!!!!!!!!  
         endif
         TopOfNode = .true.
      endif
      call c1by1


      end
      subroutine cbegin1ev(nomore)
      implicit none
#include "Zprivate.h"
#include "Ztrack.h"
#include "Ztrackv.h"
#include "Ztrackp.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"
#include "Zcode.h"
#include "Zmanager.h"
#include "Zmanagerp.h"
      
      integer nomore       !  output. 0 still there  are showers
                           !          1 no more skeleton showers to be fleshed
!          event number, primary      

      type(track)::incident
      type(track)::zsave
      type(coord)::angle
      
      integer i
      integer seed(2)
      integer cumnum, num, jeof, fin
      read( Mdev, end=1000 ) cumnum, num, SeedSave, Zfirst

      EventsInTheRun = EventsInTheRun + 1
      EventNo = EventNo + 1
!          get random seed at skelton making; this can work
!          if seed file is supplied
!      call creadSeed(SeedSave, EventNo, jeof)
!      if( jeof .ne. 0 ) goto 1000

!                 reset the seed.
      call rnd1r(SeedSave)

!         next incident; confirmed to be the same one as preserved one
      call cmkIncident(incident, fin)

      if(fin .ne. 0 ) goto 1000
      zsave = Zfirst    ! save;  this is reset in next 
      call ciniTracking( incident )   
!          set first interaction pos
      Zfirst = zsave
      call cresetTimer(Zfirst)

      RealBegin = .false.

!          do your own init for a one event here
!      ==========================================================
         call cqIncident( incident, angle)
      do i = 1, NoOfSites
         write(*, 999)
     *    sngl(ObsSites(i).pos.depth),
     *    EventNo, 
     *    incident.p.code,
     *    incident.p.subcode,
     *    incident.p.charge, 
     *    incident.p.fm.e,
     *    -angle.r(1),
     *    -angle.r(2),
     *    -angle.r(3)
 999      format(f10.3,i9,3i4,e15.5,3(1x,f12.8))
      enddo

!      ==========================================================
!

      call cgetHES(Mdev)  ! get high energy ptlcs
      call cobsHES        ! imitate their observation
      nomore = 0
      return

 1000 continue
      nomore = 1
      end

!     ************************************ hook for observation
!     *  One particel information is brought here by the system.
!     *  All information of the particle is in aTrack
!     *
      subroutine chookObs(aTrack, id)
!
!     Note that every real variable is in double  precision so
!     that you may output it in sigle precision to save the memory.
!     In some cases it is essential to put it in sigle (say,
!     for gnuplot).
! 
      implicit none
#include "Zcode.h"
#include "Ztrack.h"
#include "Zprivate.h"
      integer id  ! input.  2 ==> reached at an observation level
!                           1 ==> aTrack is going out from
!                                 outer boundery.
!                           2 ==> reached at an observation level
!                           3 ==> reached at inner boundery.
      type(track)::aTrack
!
!     For id =2, you need not output the z value, because it is always
!     0 (within the computational accuracy).
!
      if(id .eq. 2 .and. aTrack.p.code .ne. kneumu .and.
     *   aTrack.p.code .ne. kneue) then

!     ===================================================

       if( aTrack.p.code .le. 6 .and. aTrack.p.code .ne. 3 ) then
!         write(*, 959) 
!     *  aTrack.where,   
!     *  aTrack.p.code,   
!     *  aTrack.p.charge,  
!     *  sngl( aTrack.p.fm.p(4) - aTrack.p.mass ),  
!     *  sngl( aTrack.pos.xyz.r(1) ), 
!     *  sngl( aTrack.pos.xyz.r(2) ) ,  
!     *  sngl( aTrack.vec.w.r(1) ),   
!     *  sngl( aTrack.vec.w.r(2) ),  
!     *  sngl( aTrack.vec.w.r(3) ),  
!     *  sngl( aTrack.vec.coszenith )  
! 959      format(3i3,f12.3,2f16.6,4(1x,f12.8))
      endif

!     ===================================================

!         write(*,'(4i5, g15.4,g15.3)') 
!     *   aTrack.where, aTrack.p.code, aTrack.p.subcode,
!     *   aTrack.p.charge, sngl( aTrack.t ),
!     *   sngl( aTrack.p.fm.p(4) - aTrack.p.mass)
!     *   sngl( aTrack.pos.xyz.r(1) ), sngl( aTrack.pos.xyz.r(2) ),
!     *   sngl( aTrack.vec.w.r(1) ), sngl(aTrack.vec.w.r(2) ),
!     *   sngl(aTrack.vec.w.r(3) ), 
!     *   sngl(aTrack.vec.coszenith)

      endif
      end

!    *********************************** hook for end of 1 event
!    * At this moment, 1 event generation has been ended.
!    *
      subroutine chookEnEvent

      implicit none
#include "Zprivate.h"
#include "Ztrack.h"
#include "Ztrackv.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"

      integer i
      
      if(RealEnd) then
         call cfinTracking
!           real end of 1 event; if you need to do some here is
!           the place
!         ========================================================

         if(ObserveAS) then
!                   electron size in B approx.
            do i = 1, NoOfASSites
               write(*, *) ASObsSites(i).age, ASObsSites(i).esize
            enddo
         endif



!        ========================================================
!  
      else
!           there is still low energy skeleton ptcls
!           nothing to do here
      endif

      end


!     ********************************* hook for end of a run
!     *  all events have been created or time lacks
!     *
      subroutine chookEnRun
      implicit none
#include "Zprivate.h"
!     =========================================================

!     =========================================================
      call  cprintStatus   ! if don't like,  comment out
      end
!     ********************************* hook for trace
!     *  This is called only when trace > 100
!     *  User should manage the trace information here.
!     *  If you use this, you may need some output for trace
!     *  at the beginning of 1 event generatio and at the end of  1 event
!     *  generation so that you can identfy each event.
!     *
!     *
      subroutine chookTrace
            implicit none

#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Ztrackp.h"
#include  "Zobs.h"
#include  "Zobsv.h"

       real*4 h1,  h2
!
!    Every time a particle is moved in the atmosphere, this routine is called,
!    if trace > 100
!         For a one track segment,
!     TrackBefMove  has  track information at the beginning of the segment.
!     MoveTrack    has   track information at the end of the segment.
!   
!     You can know the  information a track contains in the 
!     chookObs routine. (Note however, no conversion of coordinate
!     has been done.  The values are in the Earth xyz system.)
!     Besides quantities explained there, you can use, for a  given 'track'
!
!     atrack.pos.xyz.x, atrack.pos.xyz.y, atrack.pos.xyz.z    (x,y.z)
!     atrack.pos.radiallen   (distance from the center of the earth)
!     atrack.pos.depth       (vertical depth)
!     atrack.pos.height      (vertical heigth from sea level)  
!

      h1 = TrackBefMove.pos.height- ObsSites(NoOfSites).pos.height
      h2 = MovedTrack.pos.height - ObsSites(NoOfSites).pos.height

      end
!     ********************* this is the hook called when
!       an electron made an interaction.
!
      subroutine chookEInt(never)
      implicit none
      integer never             ! input & output
      never = 1
      end

!     ********************* this is the hook called when
!       a gamma ray made an interaction.
!
      subroutine chookGInt(never)
      implicit none
      integer never   ! input & output
      never = 1
      end

!     ********************* this is the hook called when
!       non e-g particle made an interaction.
!
      subroutine chookNEPInt(never)
      implicit none
      integer never   ! input & output
      never = 1
      end


      subroutine cgetHES(from)
      implicit none
#include "Zprivate.h"
      integer from

      integer i

      read(from)  Np
      do i = 1, Np
         read(from) o(i)
      enddo
      end

      subroutine cobsHES
      implicit none
#include "Zprivate.h"
#include "Ztrack.h"
!
!           memorized high energy showers at the skeleton making
!     time is put into the chookObs as if they are really observed
      type(track)::aTrack

      integer i

      do i = 1, Np
         aTrack.where =  o(i).where 
         aTrack.p.code =  o(i).code 
         aTrack.p.subcode = o(i).subcode 
         aTrack.p.charge = o(i).charge 
         aTrack.t = o(i).atime 
         aTrack.p.fm.p(4) = o(i).erg
         aTrack.p.mass = o(i).mass 
         aTrack.pos.xyz.r(1) = o(i).x 
         aTrack.pos.xyz.r(2) = o(i).y 
         aTrack.vec.w.r(1) = o(i).wx
         aTrack.vec.w.r(2) = o(i).wy
         aTrack.vec.w.r(3) = o(i).wz
         aTrack.vec.coszenith = o(i).zenith 
         call chookObs(aTrack, 2)
      enddo
      end


!        process low energy partilces in the skeleton 1 by 1

      subroutine c1by1
      implicit none
#include "Zprivate.h"
#include "Ztrack.h"
#include "Ztrackv.h"

      character*100 msg

      call cinitStack  ! empty the stack

      if( TopOfNode ) then
         read(Mdev)  NoOfLowE,  p
         if( p.asflag .eq. -1 .and. ObserveAS )  then
            call embedAS
         endif
         NLowCounter = 0
         if( NoOfLowE .eq. -1 ) then
            RealEnd = .true.
            RealBegin = .true.
            return    ! ************
         endif
      endif

      RealBegin = .false.
      RealEnd = .false.


      if( NLowCounter .eq. NoOfLowE ) then
         TopOfNode =.true.
         return
      endif

      TopOfNode = .false.
!         still not the  end of 1 event

      read(Mdev)  c

      NLowCounter = NLowCounter + 1
      call cmove_c_stack   ! move c into stack

      end
!
      subroutine embedAS
      implicit none
#include "Zprivate.h"
#include "Ztrack.h"
#include "Zearth.h"

      type(track)::el
      
      el.pos.depth = p.depth
      el.vec.coszenith = p.coszenith
      el.pos.radiallen = p.height + Eradius
      el.pos.height = p.height
      el.p.fm.p(4) = p.erg
      el.wgt = 1.0
      call cobAS(el)
      end

      
      

      subroutine cmove_c_stack
      implicit none

#include "Zprivate.h"
#include "Ztrack.h"
#include "Zearth.h"

      type(track)::aTrack
!
!          a child of the current parent is moved to stack
!      as a track info.
!
      aTrack.pos.xyz.r(1) = p.posx
      aTrack.pos.xyz.r(2) = p.posy 
      aTrack.pos.xyz.r(3) = p.posz
      aTrack.pos.depth = p.depth 
      aTrack.pos.height = p.height
      aTrack.pos.colheight = p.colHeight 
      aTrack.t = p.atime 

      aTrack.where = p.where 

      aTrack.p.code = c.code 
      aTrack.p.subcode = c.subcode
      aTrack.p.charge  = c.charge
      aTrack.p.fm.p(1) = c.fm(1) 
      aTrack.p.fm.p(2) = c.fm(2) 
      aTrack.p.fm.p(3) = c.fm(3) 
      aTrack.p.fm.p(4) = c.fm(4) 
      aTrack.p.mass = c.mass

!        --------------- next must be compute here

      aTrack.pos.radiallen = Eradius +aTrack.pos.height
      aTrack.pos.xyz.sys = 'xyz'
      aTrack.vec.w.sys = 'xyz'
      aTrack.wgt = 1.0
      aTrack.asflag = 0

      call cresetDirec( aTrack )    ! set vec.w and coszenith 

      call cpush(aTrack)
      end
