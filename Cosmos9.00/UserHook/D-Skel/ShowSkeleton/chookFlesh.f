#include "cmain.f"
#include "chookHybAS.f"
#include "ctemplCeren.f"
#include "howmuch.f"
#include "asinfo.f"
#include "asdensity.f"
#include "crecprob.f"
!
!  *************************************** hook for Beginning of a Run
!  * At this moment, all (system-level) initialization for this run
!  * has been ended.  After this routine is executed, the system goes into the
!  * event creation loop.
!  *
      subroutine chookBgRun
      implicit none
#include "Zmanagerp.h"
#include "../../SkelFlesh/Zprivate.h"

      real*8  temp
      character*100 msg
      integer icon
      integer i
      EventNo = 0

!            namelist output
      call cwriteParam(ErrorOut, 0)
!            primary information
      call cprintPrim(ErrorOut)
!            observation level information
      call cprintObs(ErrorOut)

      call cqUHooki(1, Mdev)      ! get skeleton memo dev #
      call cqUHookc(1, msg)       ! get file name for sekelton data
      call cgetfname(msg, Mskel)  ! add host name etc if needed
      call copenfw2(Mdev, Mskel, 2, icon)
      if(icon .ne. 1) then
         call cerrorMsg(Mskel,1)
         call cerrorMsg(' could not be opened',0)
      endif

      call xBgRun
      call ihist                ! instanciate histogram
      end

!     *********************************** hook for Beginning of  1 event
!     *  All system-level initialization for 1 event generation has been
!     *  eneded at this moment.
!     *  After this is executed, event generation starts.
!     *
      subroutine chookBgEvent
      implicit none
#include "../../SkelFlesh/Zprivate.h"


      integer nomore
      call cbegin1ev( nomore )
      if( nomore .eq. 1) then
         call cerrorMsg('all events have been fleshed', 1)
         stop                   !!!!!!!!!!!!  
      endif
      call cpushInci
 
      call xBgEvent
      return
      end
      subroutine cbegin1ev(nomore)
      implicit none
#include "../../SkelFlesh/Zprivate.h"
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

      type(track):: incident, zsave
      type(coord):: angle
      
      integer i
      integer seed(2)
      integer cumnum, num, jeof, fin
      read( Mdev, end=1000, err=999 ) cumnum, num, SeedSave, Zfirst

      EventsInTheRun = EventsInTheRun + 1
      EventNo = EventNo + 1
!                 reset the seed.
      call rnd1r(SeedSave)
!         next incident; confirmed to be the same one as preserved one
      call cmkIncident(incident, fin)
      if(fin .ne. 0 ) goto 1000
      zsave = Zfirst    ! save;  this is reset in next 
      call ciniTracking( incident )   
!          set first interaction pos
      Zfirst = zsave
!      call cresetTimer(Zfirst)



!          do your own init for a one event here
!      ==========================================================


!      ==========================================================
!

      call cgetHES(Mdev)  ! get high energy ptlcs
      call cobsHES        ! imitate their observation
      nomore = 0
      return

 1000 continue
      nomore = 1
      return
 999  continue
      write(0,*) ' Mdev read err'
      stop 1111
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


      integer id  ! input.  2 ==> reached at an observation level
!                           1 ==> aTrack is going out from
!                                 outer boundery.
!                           2 ==> reached at an observation level
!                           3 ==> reached at inner boundery.
      type(track):: aTrack
!
!
      integer n, i
      real*8 eps, u

      if(aTrack.wgt .gt. 1.) then
         n=aTrack.wgt
         eps = aTrack.wgt - n
         call rndc(u)
         if(u .lt. eps) then
            n = n + 1
         endif
      else
         n = 1
      endif
      do i = 1, n 
         call xObs(aTrack, id)
      enddo

      end

!    *********************************** hook for end of 1 event
!    * At this moment, 1 event generation has been ended.
!    *
      subroutine chookEnEvent

      implicit none
#include "Ztrack.h"
#include "Ztrackv.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"

      integer i

!       for Job ='newflesh', we must call cfinTracking ourselves.
      call cfinTracking
!           end of 1 event; if you need to do some here is
!           the place

      call xEnEvent


      end


!     ********************************* hook for end of a run
!     *  all events have been created or time lacks
!     *
      subroutine chookEnRun
      implicit none
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
#include "../../SkelFlesh/Zprivate.h"
      integer from

      integer i

      read(from)  Np
      do i = 1, Np
         read(from) o(i)
      enddo
      end

      subroutine cobsHES
      implicit none
#include "../../SkelFlesh/Zprivate.h"
#include "Ztrack.h"
!
!           memorized high energy showers at the skeleton making
!     time is put into the chookObs as if they are really observed
      type(track):: aTrack

      integer i
      logical  HEobs            ! if T, currently observing 
      common /ZHEobs/ HEobs     !  particles those obsrved at skeelton making time

      HEobs = .true.
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
      HEobs = .false.
      end


!        push all  low energy partilces in the skeleton in the stack

      subroutine cpushInci
      implicit none
#include "../../SkelFlesh/Zprivate.h"
#include "Ztrack.h"
#include "Ztrackv.h"
      integer i

      type(track)::aTrack

      call cinitStack  ! empty the stack

      read(Mdev)  NoOfLowE
      do  i = 1, NoOfLowE
         read(Mdev) aTrack
!           aTrack is already complete track so push it directly.
         call cpush(aTrack)
         write(*,'(2i4,g13.4, g13.4)')
     *     aTrack.p.code, aTrack.p.charge, 
     *     aTrack.p.fm.p(4)-aTrack.p.mass, aTrack.pos.depth*0.1
      enddo
      stop
!           sort stack dscendent order
      call csortStack
      
      end
