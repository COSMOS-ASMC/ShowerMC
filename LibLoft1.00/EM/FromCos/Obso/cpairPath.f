      subroutine cpairPath(Eg, path)
!        pair below LPM region
      implicit none
#include "ZbpCnst.h"
!      entry  cPrSampP(Eg, prob, path)


      real*8 Eg            ! input. Gamma energy in GeV

      real*8 path          ! output. sampled path in r%l

      real*8 prob          ! output. pair prob. per r%l
      real*8 u


      if( Eg .lt. PrScrE ) then
!          partial screeinig region
         call cPrLSampP( Eg,  prob)
      elseif(Eg .lt. PairEgmaxL) then
!            complete screeing region
         call cPrCSampP( Eg, prob)

      else
!            LPM region
!         call cPrHSampP( Eg, prob); this cannot be used due to rho change
!      branch to LPM routine is at outside.  See cbremsPath for more 
!      detail 

!          use  complete screeing
         call cPrCSampP( Eg, prob)
      endif
      call rndc(u)
      path = -log(u)/prob
      end
!     *********************************** 
      subroutine cpairEnergy( Eg, Ee)
      implicit none
#include "ZbpCnst.h"

      real*8 Eg            ! input. gamma energy in GeV
      real*8 Ee            ! output. sampled Ee in GeV. higher energy of pair


      if(Eg .le. PairNonSc + 2.d-3) then
!          table in cPrLSampE is  used >  PairNonSc + 2 MeV
!          because of a glitch  around PairNonSc
!         (for Nelson's case, there is no glitch; glitch
!          comes from diff.  of dsigma/dx at x~xmax.
!
         call cPrTSampE(Eg, Ee)  

      elseif(Eg .le. PrScrE) then
!           partial screeinig region
         call cPrLSampE(Eg,  Ee)

      elseif(Eg .le. PairEgmaxL) then
!            complete screeing region
         call cPrCSampE( Eg, Ee)

      else
!            LPM region; same as path routine; should not come
!         call cPrHSampE( Eg, Ee)
!
        call cPrCSampE( Eg, Ee)
      endif
      end

      
