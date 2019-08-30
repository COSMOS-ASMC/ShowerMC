      subroutine epSetSTblCns(media, cnst)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"

       type(epmedia):: media  ! input. 
       type(SmpCnst)::  cnst   ! output. must be media.cnst

!   
!   Actual sampling table sizes are fixed here.
!   If they exceed the maximum availabe size, 
!   error will be issued to stop the execution. 
!   
!         Energy above which we can use complete screening 
!       cross-sections. evaluate at Eg/Ee=x= 0.99.
!       30 is marging.  This can be applied to pair creation, too.

      cnst%CompScrE =  !  30.0*( 0.05*0.99/0.01/media.Zeff**0.333)
!     *              150.0/media.Zeff**0.333
     *              150.0/media%Z**0.333
!
!      cnst%BremEgmin =  1.d-5   ! min Eg/Ee at  Ee > 100 MeV
       cnst%BrEgminS2 =  cnst%BremEgmin ! min Eg/Ee at  Ee > 100 MeV5.0d-6
!
!       For Seltzer data
!   
      cnst%BrEemaxS  = 100.d-3  ! max of the  lower Seltzer data 
                                !     region  (100MeV)
      cnst%BrEgminS = 1.0d-6    ! Eg  min for lower Seltzer region
                               ! 1 keV ***NOT RATIO****
                      ! so the ratio, Eg/Ee=0.0019569 ~ 1.e-5

      cnst%BrEeminS = masele+ 5.0d-6  ! below this, brems is neglected  (Ek=5keV)
      cnst%BrLEeminS = log10(cnst%BrEeminS-masele)  ! this is in K.E
!         table size for Energy;lower Seltzer region
      cnst%BrTXTS = mxBrTXS

      cnst%BrES = cnst%BrTXTS

      cnst%BrUminSA = 0.1d0
      cnst%BrUmaxSA = 1.0

      cnst%BrUszSA = 51
      cnst%BrdUSA =
     *    (cnst%BrUmaxSA-cnst%BrUminSA)/
     *      (cnst%BrUszSA-1)


      cnst%BrdETXS = 
     *  log10((cnst%BrEemaxS-masele)/(cnst%BrEeminS-masele))/
     *                (cnst%BrTXTS-1)

      cnst%BrdES =
     *   log10( (cnst%BrEemaxS-masele)/(cnst%BrEeminS-masele))/
     *    ( cnst%BrES-1)

      cnst%BrUszSB = 51
      cnst%BrUminSB = 0.
      cnst%BrUmaxSB = cnst%BrUminSA**0.25d0
      cnst%BrdUSB =
     *   (cnst%BrUmaxSB-cnst%BrUminSB)/
     *     (cnst%BrUszSB-1) 

! ----------------- 100MeV to 10 GeV Seltzer upper region
      cnst%BrEemaxS2  = 10.  ! max of the upper Seltzer data 
                                !     region  (10GeV)
      cnst%BrEgminS2 = 1.0d-5  ! Eg/Ee min for upper Seltzer region
                               ! *** ratio ***.
                               ! Eg=100keV to 1keV
      cnst%BrEeminS2 =  cnst%BrEemaxS ! 100. MeV
      cnst%BrLEeminS2 = log10(cnst%BrEeminS2-masele) ! this is in K.E

!         table size for Energy;upper Seltzer region
      cnst%BrTXTS2 = mxBrTXS2

      cnst%BrES2 = cnst%BrTXTS2

      cnst%BrUminSA2 = 0.1d0
      cnst%BrUmaxSA2 = 1.0

      cnst%BrUszSA2 = 51
      cnst%BrdUSA2 =
     *    (cnst%BrUmaxSA2-cnst%BrUminSA2)/
     *      (cnst%BrUszSA2-1)


      cnst%BrdETXS2 = 
     *  log10((cnst%BrEemaxS2-masele)/(cnst%BrEeminS2-masele))/
     *                (cnst%BrTXTS2-1)

      cnst%BrdES2 =
     *   log10( (cnst%BrEemaxS2-masele)/(cnst%BrEeminS2-masele))/
     *    ( cnst%BrES2-1)

      cnst%BrUszSB2 = 51
      cnst%BrUminSB2 = 0.
      cnst%BrUmaxSB2 = cnst%BrUminSA2**0.25d0
      cnst%BrdUSB2 =
     *   (cnst%BrUmaxSB2-cnst%BrUminSB2)/
     *     (cnst%BrUszSB2-1) 
!     Ee < Ee2H      ===> complete sc.+ LPM  table
!     Ee > Ee2H      ===> comp. + LPM  by rejection
!===================================================================
!


!           ordinary brems region
!
      cnst%BremEemin = 100.d-3  ! v9.14. 
!                                below this, don't make table for
!                                screened cross sec.
!                                Seltzer table will be created
!                               up to 10 GeV.
!                                  
      cnst%BremLEemin = log10(cnst%BremEemin)

!         LPM effect at small x=Eg/Ee appears at this Ee
!          next is BremEemaxL in older vesionsl   but 
!          100 times lower v9.15 Scales by X0(g/cm2)/rho(g/cm3)
!         = X0cm.  base is Pb (X0=0.56cm). So
!         for W, we consider LPM at Ee> 0.3 *0.35/0.56~ 190 MeV
!         for  Air Ee>16TeV. (X0~30000 cm)
      cnst%BremEeminLPM =
     *    max(0.1d0,  0.3d0 * media%X0/0.561)
!       At  100MeV, even in W, Eg affected by LPM is 10keV
!       (in ratio it is 1e-4) so we set absolute min for
!       LPM as above 
!       Flpm*cnst.BreEeminLPM is actual minimum for LPM (Flpm>=1)
!       Flpm canbe sepcified by the user.

!        For a while, we set const for
!        the case where LPM  is disabled

!      we use partial screened brems cross-section below BrScrE
      cnst%BrScrE = cnst%CompScrE
      cnst%BremTXTL = mxBrTXL
!     *   log10(cnst.CompScrE/cnst.BremEemin)*10
      cnst%BremEsize = cnst%BremTXTL
      cnst%BremUminLA = 0.1d0
      cnst%BremUmaxLA = 1.0
      cnst%BremUszLA = 51
      cnst%BremdULA =
     *    (cnst%BremUmaxLA-cnst%BremUminLA)/
     *      (cnst%BremUszLA-1)
      cnst%BremdETXL = 
     *  log10(cnst%BrScrE/cnst%BremEemin)/
     *                (cnst%BremTXTL-1)
      cnst%BremdEL =
     *   log10(cnst%BrScrE/cnst%BremEemin)/
     *    ( cnst%BremEsize-1)
      cnst%BremUszLB = 51
      cnst%BremUminLB = 0.
      cnst%BremUmaxLB = cnst%BremUminLA**0.25d0
      cnst%BremdULB =
     *   (cnst%BremUmaxLB-cnst%BremUminLB)/
     *     (cnst%BremUszLB-1) 


!       for pair creation 
!      ------------------------------------
      cnst%PairEgmin =  1.023d-3
      cnst%PairNonSc = 8.d-3

      cnst%PairLEgmin = log10(cnst%PairEgmin) 

!         LPM effect  appears at this Eg
      cnst%PairEgmaxL = 3000.*11.35/media%rho
!         we should use screened pair cross-section below this
!           probably always  CompScrE < PairEgmaxL
      cnst%PrScrE = min(cnst%CompScrE,  cnst%PairEgmaxL)

      cnst%PairTXTL =  mxPrTXL             
      cnst%PairEsize = cnst%PairTXTL
      cnst%PairUminLA  = 0.05d0
      cnst%PairUmaxLA = 1.
      cnst%PairUszLA = 51
      cnst%PairdULA = 
     *  (cnst%PairUmaxLA-cnst%PairUminLA)/
     *  (cnst%PairUszLA-1) 

      cnst%PairdETXL =
     *  log10(cnst%PrScrE/cnst%PairEgmin)/
     *   (cnst%PairTXTL-1)

      cnst%PairUminLB = 0.
      cnst%PairUmaxLB =   cnst%PairUminLA**0.25d0

      cnst%PairUszLB = 51
      cnst%PairdULB =
     * (cnst%PairUmaxLB-cnst%PairUminLB)/
     *   (cnst%PairUszLB-1)

      cnst%PairdELA = 
     * log10(cnst%PrScrE/ cnst%PairEgmin) /
     *  (cnst%PairEsize-1)

      cnst%PairdELB=
     *      sqrt( log10(cnst%PrScrE/cnst%PairEgmin) )
     *     /( cnst%PairEsize-1)

!      --------------------------------------------------

!           LPM treatment for Brems  For Ee >=EeminLPM (all region)

      cnst%BrEe1H =  cnst%BremEeminLPM
      cnst%BrEgminH =  cnst%BremEgmin   !  vmin  Eg/Ee  10^-5
      cnst%BrLEe1H = log10( cnst%BrEe1H )
      cnst%BrneH = 71  !  v9.15 (old value 51)
      cnst%BrdEH = .10d0  ! log E step
!              last 2  means  for 2D table
      cnst%BrneH2 = 71 ! v9.15 old value 51
      cnst%BrdEH2 = 0.05d0

!       cnst.BrdEH= log10(cnst.BrEe2H/cnst.BrEe1H)/(cnst.BrneH-1)
!       inverse of the above;  Pb: thie is 1000TeV, Air at sea level
!       3x10^19eV,   This is used when getting x-section
!          above this, 1/sqrt(rE), may be multiplied 
      cnst%BrEe2H =10.d0**( cnst%BrdEH *(cnst%BrneH-1)) *
     * cnst%BrEe1H

!       Next is used when sampling energy by rejection.
!       could be lower.
!       If Ee > this, rejection method will work
!       Pb: 3.2TeV; Air sea: 4x10^16eV 
!                   Air 20 km: 4x10^17eV
      cnst%BrEe2H2 =10.d0**( cnst%BrdEH2 *(cnst%BrneH2-1)) *
     * cnst%BrEe1H

      cnst%BrdU1H= 0.025d0
      cnst%BrU1H=0.2d0
      cnst%BrU2H=1.0d0
      cnst%Brnu1H=(cnst%BrU2H-cnst%BrU1H+1.d-9)/cnst%BrdU1H+1
!           ....................
!                                 
      cnst%BrPow = 4.0
      cnst%BrU3H = .0
      cnst%BrU4H = cnst%BrU1H**(1.d0/cnst%BrPow)
      cnst%Brnu2H = 51
      cnst%BrdVU2H = cnst%Brnu2H-1
      cnst%BrdU2H =  (cnst%BrU4H - cnst%BrU3H)/cnst%BrdVU2H
!         -----------

!           LPM region, pair

      cnst%PrEg1H =  cnst%PairEgmaxL 
      cnst%PrLEg1H = log10( cnst%PrEg1H )
!     cnst.PrneH = 25
      cnst%PrneH = 51
      cnst%PrdU1H= 0.025d0
      cnst%PrdEH = .1d0  ! log E step

      cnst%PrU1H=0.
      cnst%PrU2H=1.0
      cnst%Prnu1H=(cnst%PrU2H-cnst%PrU1H+1.d-9)/cnst%PrdU1H+1
      cnst%PrEg2H =10.d0**( cnst%PrdEH *(cnst%PrneH-1)) *  
     * cnst%PrEg1H
!            for muon 
      call epSetmuSTab(media, cnst)
      end
