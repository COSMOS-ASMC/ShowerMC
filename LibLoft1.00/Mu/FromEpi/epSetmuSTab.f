      subroutine epSetmuSTab(media, cnst)
      implicit none
#include "Zmedia.h"

       type(epmedia):: media  ! input. 
       type(SmpCnst)::  cnst   ! output. must be media.cnst

!      -   for muon interactions

!           nuclear interaction

!        media.Aeff  and media.Zeff may be  better to be 
!        replaced by medi.A and meia.Z in future when
!        all Tab is recreated.
!
      cnst%muNVmin = 1.d-3
      cnst%muNdETX = 0.1d0  ! log E step for total prob. dE/dx
      cnst%muNdE  = 0.2d0   ! log E step for sampling 
      cnst%muNEmin = 100./media%A**0.333  ! below this, contrib. < 1%
      cnst%muNLEmin = log10(cnst%muNEmin)
      cnst%muNEmax = 100.d3  ! upto 100 TeV. above this, use some scaling
                             ! this will be reset to an integral value
                             ! after  Esize is fixed
      cnst%muNdU = 1.d-2

      cnst%muNUsize = 1.00001d0/cnst%muNdU + 1
      cnst%muNTXT =
     *     (log10(cnst%muNEmax/cnst%muNEmin)+1.d-7)/cnst%muNdETX +1
!          this will be bit lower than 100 TeV.
      cnst%muNEmax1=
     *  cnst%muNEmin*10.d0**((cnst%muNTXT-1)*cnst%muNdETX)

      cnst%muNEsize =
     *     (log10(cnst%muNEmax/cnst%muNEmin)+1.d-7)/cnst%muNdE +1
      cnst%muNEmax =
     *   cnst%muNEmin*10.d0**( (cnst%muNEsize-1) *cnst%muNdE)

!
!            pair creation
!
      cnst%muPrVmin = 1.d-3
      cnst%muPrEmin = 70./media%Z**0.666
      cnst%muPrLEmin = log10(cnst%muPrEmin)
      cnst%muPrEmax = 100.d3
      cnst%muPrdETX = 0.1d0
      cnst%muPrdE = 0.2d0

      cnst%muPrTXT = 
     *     (log10(cnst%muPrEmax/cnst%muPrEmin)+1.d-7)/cnst%muPrdETX +1
      cnst%muPrEmax1 =
     *    cnst%muPrEmin*10.d0**((cnst%muPrTXT-1)*cnst%muPrdETX)

      cnst%muPrEsize =
     *     (log10(cnst%muPrEmax/cnst%muPrEmin)+1.d-7)/cnst%muPrdE +1
      cnst%muPrEmax =
     *    cnst%muPrEmin*10.d0**( (cnst%muPrEsize-1)*cnst%muPrdE)

      cnst%muPrdU = 0.01d0
      cnst%muPrUsize = 1.00001d0/cnst%muPrdU + 1

!
!            brems
! 
      cnst%muBrVmin = 1.d-3
      cnst%muBrdETX = 0.1d0
      cnst%muBrdE = 0.2d0
      cnst%muBrEmin = 100./media%Z**0.666
      cnst%muBrLEmin = log10(cnst%muBrEmin)
      cnst%muBrEmax = 100.d3
      cnst%muBrTXT = 
     *     (log10(cnst%muBrEmax/cnst%muBrEmin)+1.d-7)/cnst%muBrdETX +1
      cnst%muBrEmax1 =
     *      cnst%muBrEmin*10.d0**((cnst%muBrTXT-1)*cnst%muBrdETX)
      cnst%muBrEsize = 
     *     (log10(cnst%muBrEmax/cnst%muBrEmin)+1.d-7)/cnst%muBrdE +1
      cnst%muBrEmax =
     *    cnst%muBrEmin*10.d0**( (cnst%muBrEsize-1)*cnst%muBrdE)

      cnst%muBrdU = 0.01d0
      cnst%muBrUsize = 1.00001d0/cnst%muPrdU + 1
      end
