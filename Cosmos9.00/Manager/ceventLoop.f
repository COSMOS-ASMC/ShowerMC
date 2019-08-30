!              eventLoop
      subroutine ceventLoop
      use modEMcontrol
      implicit none
#include  "Zmanagerp.h"
#include  "Zmanager.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"      
#include  "Zprimary.h"
#include  "Zprimaryv.h"
#include  "Zcondc.h"
#include  "Ztrackp.h"


      type (track):: incident
      integer jold, twoir(2)
      integer fin, jeof


      do while (EventNo .lt. abs( DestEventNo(1)) .and. 
     *      EventsInTheRun .lt. abs( DestEventNo(2)))
!             switch to the random number generator 1
         call rndsw(jold, 1)
         if(Job .ne. 'newflesh') then
!              for newflesh, these are managed in chookBgEvent
!c            EventsInTheRun = EventsInTheRun + 1  ! moved to cmkIncident.f
!c            EventNo = EventNo + 1                ! from  v6.25
            if(RefreshIR) then
!               special case. init for ir of this event is
!               read from 14 which should have been opend by the
!               user in chookBrRun
               read(14, *, end=10) twoir
               call rnd1r(twoir)
            endif
            if(Job .eq. 'flesh') then
!                get random seed at skelton making
               call creadSeed(SeedSave, EventNo, jeof)
               if( jeof .ne. 0 ) goto 10
!                 reset the seed.
               call rnd1r(SeedSave)
            endif
!                 save the initial seed of random number
            call rnd1s(SeedSave)
            write(0,*) SeedSave, EventNo+1
            call cmkIncident(incident, fin)
            if(fin .ne. 0) goto 10 ! all including cut-offed end
            if(Job .eq. ' ' .and. SeedFile .ne. ' ') then
               call cwriteSeed
            endif

            call ciniTracking( incident ) ! init for tracking each event
            call cinitStack
            call cpush(incident)
         endif
!     if(Eabsorb(1) .ne.  0 ) then
         if(Eabsorb /=  0 ) then
            call chookEabsorbi(0)
         endif
         if(Job .eq. 'newflesh') then
!                  this may be used always; but next one is original
            call chookBgEvent   ! user hook for each event beginning
            call csetMolUnit    ! for newflesh, primary is fixed above
         else
            call csetMolUnit
            call chookBgEvent   ! user hook for each event beginning
         endif
         call ctrackingAll
         if(Job .ne. 'newflesh') then
!              next is managed by chookEnEven if new flesh
            call cfinTracking   ! enf of tracking of each event
         endif
         call chookEnEvent    ! user hook for each event  end

      enddo
 10   continue
      end
      subroutine  csetMolUnit
#include "Ztrack.h"
! #include "Zmagfield.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"

!         set  Moliere unit at each observtion depth (2r.l above along primary)
      type (track)::inci
      type (coord)::angle
      integer i
      real*8 mu
      call cqIncident(inci, angle)

      do  i = 1, NoOfSites
!         call cgetMoliereU(ObsSites(i).pos.depth, inci.vec.coszenith, ! bef v7.
          call cgetMoliereU(ObsSites(i)%pos%depth, -angle%r(3),
     *    ObsSites(i)%mu)
      enddo
      do i = 1, NoOfASSites
!         call cgetMoliereU(ASObsSites(i).pos.depth, inci.vec.coszenith, 
          call cgetMoliereU(ASObsSites(i)%pos%depth, -angle%r(3),
     *       ASObsSites(i)%mu)
      enddo
      end
