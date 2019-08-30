#include "cmain.f"
#include "chookHybAS.f"
#include "ctemplCeren.f"
      subroutine time(xxx)
      integer xxx
      xxx = 1
      end
      
!  *************************************** hook for Beginning of a Run
!  * At this moment, all (system-level) initialization for this run
!  * has been ended.  After this routine is executed, the system goes into the
!  * event creation loop.
!  *
      subroutine chookBgRun
      implicit none
#include "Zmanagerp.h"

 
!         If you feel writing the parameters on stderr is
!         a bother, comment out the next or
!         use other device than ErrorOut.
!         Also you may comment out all output routines below.

!            namelist output
      call cwriteParam(ErrorOut, 0)
!            primary information
      call cprintPrim(ErrorOut)
!            observation level information
      call cprintObs(ErrorOut)
      end


!     *********************************** hook for Beginning of  1 event
!     *  All system-level initialization for 1 event generation has been
!     *  eneded at this moment.
!     *  After this is executed, event generation starts.
!     *
      subroutine chookBgEvent
      implicit none
#include "Ztrack.h"
#include "Zobs.h"
#include "Zobsv.h"
      include "Zprivate.h"

      integer i, j
      sume = 0
      sumepi0 = 0
      nc = 0
      ncpi0= 0 
      call cqUhookr(1, cut)
      end
  

!     ************************************ hook for observation
!     *  One particle information is brought here by the system.
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
#include "Ztrack.h"
#include  "Zcode.h"
      include "Zprivate.h"

      integer id  ! input.  1 ==> aTrack is going out from
!                                 outer boundery.
!                           2 ==> reached at an observation level
!                           3 ==> reached at inner boundery.
      type(track):: aTrack
!
!     For id =2, you need not output the z value, because it is always
!     0 (within the computational accuracy).
!
!      if(id .eq. 2 .and.  aTrack.vec.coszenith .gt. 0 ) then
!         if(aTrack.p.code .le. 3) then
!            nc(aTrack.p.code, aTrack.where) = 
!     *        nc(aTrack.p.code, aTrack.where) +1
!         endif
!      endif
!            output typical quantities.
!        write(*, '(3i5, ) 
!     *  aTrack.where,   !  observation level. integer*2.  1 is highest.
!     *  aTrack.p.code,    !  ptcl code.  integer*2.
!     *  aTrack.p.charge,  !  charge,  integer*2 
!     *  sngl(aTrack.t), !  relateive arrival time in nsec (NOT sec).
!                        !  if TimeStructure is F, nonsense.
!     *  sngl(aTrack.p.fm.p(4))  ! total energy in GeV.
!     *  sngl(aTrack.pos.xyz.r(1)), sngl(aTrack.pos.xyz.r(2)), !  x, y in m
!     *  sngl(aTrack.vec.w.r(1)),  ! direc. cos.x in the current detector system.
!     *  sngl(aTrack.vec.w.r(2)),  ! direc. cos.y
!     *  sngl(aTrack.vec.w.r(3)),  ! direc. cos.z
!     *  sngl(aTrack.vec.coszenith) ! cos of zenith angle
!      endif
!         you may need in some case other information such as
!       aTrack.p.subcode   ! sub code of the particle integer*2
!       aTrack.p.mass      ! mass 
!       aTrack.wgt         ! weight of the particle (may not be 1. if
                           ! ThinSampling =T)
!       aTrack.p.fm.p(1)      ! momentum x component.  Note. Momentum is
!                            given in the  Earth xyz system.

!       aTrack.p.fm.p(2)      !          y
!       aTrack.p.fm.p(3)      !          z

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
      include "Zprivate.h"

      integer i, j
!      do i = 1, NoOfSites
!         do j = 1, 3
!            write(*, *) i, j, nc(j, i)
!         enddo
!      enddo
      if(ObserveAS) then
!                   electron size in B approx.
!         write(*,*) 'e ', sngl(sume), nc, sngl(sumepi0), ncpi0
         do i = 1, NoOfASSites
            write(*, *) 's ',  sngl( ASDepthList(i)*0.1 ),
     *       sngl(ASObsSites(i).esize), sngl(ASObsSites(i).age)
         enddo
      endif
!      write(*,*)
      end


!     ********************************* hook for end of a run
!     *  all events have been created or time lacks
!     *
      subroutine chookEnRun

      implicit none
#include "Ztrack.h"
#include "Ztrackp.h"




      end


!     ********************* this is the hook called when
!       an electron made an interaction.
!
      subroutine chookEInt(never)
            implicit none

#include  "Ztrack.h"
#include  "Ztrackv.h"
!  #include  "Ztrackp.h"
      
      integer never   ! input & output
      
!         don't make never = 1, if you want to get
!         information after an electron made interaction
!         if this is made non zero, this routine will never be called.
!
!   MovedTrack is the electron that made interaction
!   Pwork contains produced particles.
!   Nproduced has the number of particles in Pwork
!   IntInfArray(ProcessNo) contains the type of interaction
!
!        default setting
      never = 1
!
!        IntInfArray(ProcessNo).process will have one of
!       'brems', 'mscat', 'bscat',or  'anihi'
!
      end

!     ********************* this is the hook called when
!       a gamma ray made an interaction.
!
      subroutine chookGInt(never)
            implicit none

#include  "Ztrack.h"
#include  "Ztrackv.h"
!  #include  "Ztrackp.h"
      
      integer never   ! input & output
      
!         don't make never = 1, if you want to get
!         information after a gamma ray made interaction
!         if this is made non zero, this routine will never be called.
!
!   MovedTrack is the gamma that made interaction
!   Pwork contains produced particles.
!   Nproduced has the number of particles in Pwork
!   IntInfArray(ProcessNo) contains the type of interaction
!
!        default setting
      never = 1
!         IntInfArray(ProcessNo).process will have one of
!        'pair', 'comp', 'photoe' or 'photop'
!       
      end

!     ********************* this is the hook called when
!       non e-g particle made an interaction.
!
      subroutine chookNEPInt(never)
            implicit none

#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Zcode.h"
      include "Zprivate.h"
      integer never   ! input & output
      integer i

!         don't make never = 1, if you want to get
!         information after a non-e-g particle  made interaction
!         if this is made non zero, this routine will never be called.
!
!   MovedTrack is the particle that made interaction
!   Pwork contains produced particles.
!   Nproduced has the number of particles in Pwork
!   IntInfArray(ProcessNo) contains the type of interaction
!
!        default setting is 1
      never = 1

!
!        IntInfArray(ProcessNo).process  will have
!             'coll' or 'decay'
      if(cut .lt. 1.) then
         if( IntInfArray(ProcessNo).process   .eq. 'coll' ) then
            do i = 1, Nproduced
!              for pions, kaons put E=mass if X > cut (for cut>0)
               if( Pwork(i).code .lt. knuc ) then
                  if( Pwork(i).fm.p(4) / MovedTrack.p.fm.p(4)
     *                 .gt. cut) then
                     sume = sume +  Pwork(i).fm.p(4)
                     nc = nc +1
                     if( Pwork(i).code .eq. kpion .and.
     *                   Pwork(i).charge .eq. 0 ) then
                        sumepi0 = sumepi0 + Pwork(i).fm.p(4)
                        ncpi0= ncpi0 + 1
                     endif
                     Pwork(i).fm.p(4) = Pwork(i).mass*1.1
                  endif
               endif
            enddo
         endif
      endif
      end


!     ********************************* hook for trace
!     *  This is called only when trace > 60
!     *  User should manage the trace information here.
!     *  If you use this, you may need some output for trace
!     *  at the beginning of 1 event generatio and at the end of  1 event
!     *  generation so that you can identfy each event.
!     *

      subroutine chookTrace
            implicit none

#include  "Ztrack.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zobs.h"
#include  "Zobsv.h"


!      h1 = TrackBefMove.pos.height- ObsSites(NoOfSites).pos.height
!      h2 = MovedTrack.pos.height - ObsSites(NoOfSites).pos.height
      
      type(coord)::  f, t



      real*8  xxx/-1.d37/, yyy/-1.d36/, zzz/1.d34/
      integer kkk/-1000/, chg/-1000/
      save xxx, yyy, zzz, kkk, chg

      if( MovedTrack.p.charge .eq. 0 ) return      
!               convert coord.
      call ccoordForTr( 21,  f,  t )
      
 
      if(kkk .ne. MovedTrack.p.code .or. f.r(1) .ne. xxx
     *    .or. f.r(2) .ne. yyy .or. f.r(3) .ne. zzz  .or.
     *     chg .ne. MovedTrack.p.charge ) then
         if(xxx .ne. -1.d37) then
!            write(TraceDev, *) 
!            write(TraceDev, *) 
            write(*, *) 
            write(*, *) 
         endif
!         write(TraceDev, '(4g16.8, i3,i3)')
         write(*, '(4g12.4, i3,i3)')
     *        f.r(1), f.r(2), f.r(3), TrackBefMove.t, 
     *        TrackBefMove.p.code, TrackBefMove.p.charge
         
      endif
!      write(TraceDev, '(4g16.8,i3,i3)') 
      write(*, '(4g12.4,i3,i3)') 
     *        t.r(1), t.r(2), t.r(3), MovedTrack.t,
     *        MovedTrack.p.code,   MovedTrack.p.charge
         
      xxx = t.r(1)
      yyy = t.r(2)
      zzz = t.r(3)
      kkk = MovedTrack.p.code
      chg = MovedTrack.p.charge
      end
