!        csampGIntL: renamed as below.>   sample photon interaction length in Cosmos
      subroutine cpreSampGIntL( media )
      use modEMcontrol
!       use modXsecMedia
!
!        pair, compton, photo-electric eff. photo-hadron
!
       implicit none
#include  "Zglobalc.h"
#include  "Zmedia.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zevhnv.h"
#include  "Zelemagp.h"

       type(epmedia),intent(in):: media

       real*8  t, tkgpm2, cxai, den, cvh2den, Eg
       real*8  prob, u, cmPairMFP, mfp, xs, path
       real*8  tcomp, tphot, tpair, tgp, tcoh
       real*8 xprob(5), txray(5)
            
       real*4  xsec(5)  !  coh, incoh,  P%E  1/(g/cm2)  
       real*8 xaimin/0.1/

!     real(8):: rrho  ; now in modEMcontrol
  
!       integer,save::nzero=0  not used now
!/////////////
!
       Eg = TrackBefMove%p%fm%p(4)
!         LPM will be treated in cepSampGIntL (its slave)

!     if(LpmEffect .and. Eg > LpmPairEmin) then
!     This LpmPairEmin .equv. LPMPairEmin. must be safe for the change of 
!     density below.          
!     LpmPairEmin is in 1e9 (GeV) in modEMcontrol.  resetable.
!     Pb ->10TeV.  10000*Pb-rho/Air-rho 10000*11.4/1.25e-3~ ~10^8 GeV.
!     At higher alt.  rho is 1/10 so this could be 10^9 GeV.
!     To be more safe, the value may be set by 1/3 :    3      
!                    epPrSampP
!     den = cthick2den(TrackBefMove.pos.depth)
!          den = cvh2den(TrackBefMove%pos%height) ! better than abvoe; unit is   kg/m3
!          rrho = den*1.0d-3/media%rho   ! kg/m3=1000g/10^6cm3 = 10^-3 g/cm3. rho is in g/cm3
!       else
!          rrho = media%rhoc    ! this is 1 for atmosphere
!     endif
!      rrho should have been set after cpop.

!       rrho will not be used if LPM is absecnt
       call cepSampGIntL( media, TrackBefMove%p)
!            paths set by above sub are in r.l
       
!           magnetic  pair.  min should be changed if non Earth
       if(MagPair .eq. 1 .and. Eg .gt.  MagPairEmin) then
!               MagPairEmin must be prepared  so that
!               Xai > 0.05 or so.     
          
!           Xai = cxai(TrackBefMove%p, Mag)
!     if(Xai .gt. xaimin) then
!         Mag has been fixed at tracking start time.
          if( Mag%sys /= 'xyz' ) then
             write(0,*) 'Mag ooordinate is not xyz in cpreSampGIntL'
             stop
          endif
!                             Xai will be saved in modSetIntInf
          call epmpairp(TrackBefMove%p, Mag, Xai, mfp, path)   
!              mfp, path : in m.
!           
          call csetIntInf(path, .true., 'mpair')
        endif

       end
