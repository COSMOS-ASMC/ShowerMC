      subroutine epPrSampP(media, Eg,  prob, path)
      use modEMcontrol
!          generic brems sampling routine
      implicit none
#include "Zmedia.h"
!  #include "ZepTrackp.h" 

      type(epmedia),intent(in):: media  ! input.
      real*8 Eg                 ! input. Gamma energy in GeV
!     real(8),intent(in):: rrho ! ratio= rho(true)/rho(norminal).
! rrho is now in modEMcontrol
!      
      !  rho(norminal) is media%rho. In EPICS, it is media%rhoc.
      !  In the case of Cosmos, density changes depending on the
      !  height. rrho=density(current location)/maeida%rho.    
      real*8 prob          ! output. pair prob. per r.l
      real*8 path          ! output. sampled path in r.l

      real*8 u
!      logical LPMworks   ! now in modEMcontrol

      LPMworks = .false.
      if( Eg <  media%cnst%PrScrE ) then
!          partial screeinig region
         call epPrLSampP(media, Eg,  prob)
      else
         
         LPMworks = LPMeffect .and.
!     *        ( Eg > media%cnst%PairEgmaxL/rrho)
     *        ( Eg > media%cnst%PairEgmaxL/media%rhoc )
         if( LPMworks ) then
            call epPrHSampP(media, Eg,  prob)
         else
!     complete screeing region
            call epPrCSampP(media, Eg, prob)
         endif
      endif
!         if(LPMeffect) then
!c            LPM region
!            if(media.rhoc .eq. 1.) then
!               call epPrHSampP(media, Eg, prob)
!            elseif(media.name(1:3) .eq. "Air") then
!               rhoin = media.rho*media.rhoc*1.e3  ! kg/m3
!               call cpairLPMXsec(Eg, rhoin, prob)
!            elseif( abs(media.rhoc-1.) .gt. 0.02 ) then
!               write(0,*)'  rhoc=', media.rhoc,' for ', media.name,
!     *         ' is not close to 1 (> 2%) with LPM effect; dangerous'
!               stop
!            else
!               call epPrHSampP(media, Eg, prob)
!            endif
!         else
!c           use  complete screeing region
!            call epPrCSampP(media, Eg, prob)
!         endif
!      endif
      call rndc(u)
      path = -log(u)/prob
      end
!     *********************************** 
      subroutine epPrSampE(media, Eg,  Ee)
      use modEMcontrol
      implicit none
#include "Zmedia.h"
!  #include "ZepTrackp.h"
       type(epmedia):: media  ! input.
      real*8 Eg    ! input. gamma energy in GeV
      real*8 Ee    ! output. sampled Ee in GeV. higher energy of pair

!      logical LPMworks  now in modEMcontrl

!      LPMworks = LPMeffect .and.
!     * ( Eg > media%cnst%PairEgmaxL/media%rhoc )
!     * ( Eg > media%cnst%PairEgmaxL/rrho )
!     * ( Eg > media%cnst%PairEgmaxL/min(media%rhoc, 1.d0) )
!  LPMworks has been set T/F in epPrSampP

      if( LPMworks ) then
         call epPrHSampE(media, Eg,  Ee)
      else
         if(Eg .le. media%cnst%PairNonSc + 2.d-3) then
!          table in epPrLSampE is  used >  PairNonSc + 2 MeV
!          because of a glitch  around PairNonSc
!         (for Nelson's case, there is no glitch; glitch
!          comes from diff.  of dsigma/dx at x~xmax.
!
            call epPrTSampE(Eg, Ee) ! media info not needed.
         elseif(Eg .le. media%cnst%PrScrE) then
!           partial screeinig region
            call epPrLSampE(media, Eg,  Ee)
         else
!            complete screeing region
            call epPrCSampE(media, Eg, Ee)
         endif
      endif
      end
