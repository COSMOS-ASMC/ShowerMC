       subroutine cpreSampEIntL(mediax )
!       use modXsecMedia
      use modMCScontrol
      use modEMcontrol
!*************
      use modSetIntInf
!************* 
!
!         brems, knock-on, and anihilation are considered
!
       implicit none
#include  "ZmediaLoft.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zelemagp.h"

!
       type(epmedia),intent(in):: mediax
!**************
       integer::i
!*************
       
       real*8  t
       real*8  cupsilon, cmBremMFP, dl, u
       real*8  ee, den, cvh2den, prob, syncmfp
!     
       ee =  TrackBefMove%p%fm%p(4)

!     density calc. needs some time, we skip if den is not needed.
!     in EPICS LpmBremEmin is not used       Sync .equiv. MagBrem

!     den = cthick2den(TrackBefMove.pos.depth)
!!          den = cvh2den(TrackBefMove%pos%height) ! kg/m3
                         ! --> g/cm3:  1000/10^6
!!          rrho = den*1.0d-3/mediax%rho ! 
!                   if LPMworks becomes F;
!                   but for Sync, always needed.
!!       else
!!          rrho = mediax%rhoc  !this is 1
!     !       endif
!     rrho =  mediax%rhoc  !  this should have been set at ptcl pop time
       
!     !!!  Ee > media%cnst%BremEeminLPM*Flpm/rrho   is tested inside
       !! next, if LPMeffect is T.  
       call cepSampEIntL(media, TrackBefMove%p)
          
!       if(Sync ==  2 .and.     ! this is also OK since Sync equv MagBrem
       if(MagBrem == 2 .and.
     *      mediax%X0/rrho .gt. 10.d5 ) then
          if(ee .gt. MagBremEmin) then  ! revived May. 15, 2019 
!             compute upsilon   ! Mag is already given in Ztrackv.h
!             Upsilon = cupsilon(TrackBefMove%p, Mag)
!             if( Upsilon .gt. UpsilonMin ) then
!                call rndc(u)                           ! must be 0 now
!                t = -log(u) *  cmBremMFP(ee, Upsilon, 0.d0)
!                call csetIntInf(t, .true., 'mbrem')
!             endif
!     endif
                        !  Upsilon is saved in modSetIntInf
             call epsyncp(TrackBefMove%p, Mag, Upsilon, syncmfp, dl)
!     call csetIntInf(dl, .true., 'mbrem')  ! save dl in m
             call csetIntInf(dl, .true., 'sync') ! save dl in m
          endif
       endif
       if( doNewMCS ) then
          call cfixMCSmodel( TrackBefMove%p )
       else
          ActiveMCS = 'Mol'
       endif
       if( ActiveMCS == 'El_con' )  then
          call cfixMixedConst(MediaNO, int(TrackBefMove%p%charge ) )
       elseif(ActiveMCS == 'El_hin' )  then
          write(0,*) ' El_hin  is not supported in Cosmos'
          stop
       endif
       end


!!!!        csampEIntL:  sample electron interaction length.
!!!       subroutine csampEIntL(media)
!!!!       use modXsecMedia
!!!       use modMCScontrol
!!!!
!!!!         brems, knock-on, and anihilation are considered
!!!!
!!!       implicit none
!!!#include  "Zmedia.h"
!!!#include  "Ztrack.h"
!!!#include  "Ztrackp.h"
!!!#include  "Ztrackv.h"
!!!#include  "Zelemagp.h"
!!!
!!!!
!!!       type(epmedia),intent(in):: media
!!!       
!!!       real*8  t
!!!       real*8  cupsilon, cmBremMFP, u
!!!       real*8  ee, den, cvh2den, prob
!!!!     
!!!       ee =  TrackBefMove%p%fm%p(4)
!!!       if(LpmEffect .and. ee .gt. LpmBremEmin) then
!!!!             den = cthick2den(TrackBefMove.pos.depth)
!!!             den = cvh2den(TrackBefMove%pos%height)
!!!             call cbremLPMXsec(ee, den, prob)
!!!             call rndc(u)
!!!             t = -log(u) / prob 
!!!       else
!!!          call cbremsPath(ee, t)
!!!       endif
!!!       call csetIntInf(t*X0, .false., 'brems')
!!!       if(Knockon) then
!!!          if(TrackBefMove%p%charge .eq. -1) then
!!!             call cmollerPath(ee,  RecoilKineMinE, prob, t)
!!!             call csetIntInf(t*X0, .false., 'mscat')
!!!          else
!!!             call cbhabhaPath(ee,  RecoilKineMinE, prob, t)
!!!             call csetIntInf(t*X0, .false., 'bscat')
!!!          endif
!!!       endif
!!!       if(TrackBefMove%p%charge .eq. 1 .and. 
!!!     *    ee .lt. AnihiE ) then
!!!!          if(ee - TrackBefMove.p.mass .lt. 100.e-6) then
!!!!             t = 0.
!!!!          else
!!!             call canihiPath(ee, prob, t)
!!!!          endif
!!!          call csetIntInf(t*X0, .false., 'anihi')
!!!       endif
!!!       if(MagBrem .eq. 2) then
!!!          if(ee .gt. MagBremEmin) then
!!!!             compute upsilon
!!!             Upsilon = cupsilon(TrackBefMove%p, Mag)
!!!             if( Upsilon .gt. UpsilonMin ) then
!!!                call rndc(u)
!!!                t = -log(u) *  cmBremMFP(ee, Upsilon, 0.d0)
!!!                call csetIntInf(t, .true., 'mbrem')
!!!             endif
!!!          endif
!!!       endif
!!!
!!!       if( doNewMCS ) then
!!!          call cfixMCSmodel( TrackBefMove%p )
!!!       else
!!!          ActiveMCS = 'Mol'
!!!       endif
!!!       if( ActiveMCS == 'El_con' )  then
!!!          call cfixMixedConst(1, int(TrackBefMove%p%charge ) )
!!!       elseif(ActiveMCS == 'El_hin' )  then
!!!          write(0,*) ' El_hin  is not supported in Cosmos'
!!!          stop
!!!       endif
!!!       end
