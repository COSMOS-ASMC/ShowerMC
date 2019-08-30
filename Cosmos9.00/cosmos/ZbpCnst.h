!            constants used for pair/brem in Air

          real*8  mbtoPX0  ! mbtoPgrm x X0g.  If multiplied to sigma in mb, 
                             !  we obtain probability /radation length.

          real*8  X0g   ! radiation length.  in g/cm^2

          real*8  cScrC1     ! const which appears in the complete screening
                              ! crossection
          real*8  cScrC2     ! the other such one
          real*8  cScrMain   ! (4/3C1 + C2)


          real*8 CompScrE   !  Energy above which we can use complete screening
                            !  cross-sections. evaluate at Eg/Ee=x= 0.99.
                            !  10 is marging. 
                            !  This can be applied to pair creation, too.
                            !  10.0*( 0.05*0.99/0.01/media.Zeff**0.333)
                            !  ~ 50GeV/Z**0.33
          real*8 BrScrE     !  below this, screened cross-section is used.
          real*8 BremEgmin  ! Brems gamma min. energy in Gev (100 keV)

          real*8 BremEemin  ! Brems by e+/e- is considered at E > BremEemin 
                            ! (Ek ~ 300  keV )
          real*8 BremLEemin ! log10 of BremEemin
          real*8 BremEemaxL ! Energy of e+/e- below which LPM is negelected
                            ! for brems.   5*rho_Pb/rho GeV.
          integer BremTXTL  ! Size of the Brems total x-section table.
                            ! in the energy region BremEemin ~ BremEemaxL:
                            ! log10(BremEemaxL/BremEemin)*10
          integer BremEsize ! Size of log10 energy for 2D  table for
                            !  brems in the region A.  BremTXTL/2
          real*8 BremUminLA ! min of  uniform random number  in the region A
                            ! at energies BremEemin ~ BremEemaxL: 0.05
          real*8 BremUmaxLA ! max of  uniform random number in the region A
                            ! at energies BremEemin ~ BremEemaxL: 1.0
          integer BremUszLA ! Size  of uniform random nubmbers for 2D table
                            ! for brems in the region A.: 20
          real*8 BremdULA   ! step of u in region A  at Low energies
          real*8 BremdETXL   ! log10E  step for  brem total  cross secton
                            ! log10(BremEemaxL/BremEemin)/(BremTXTL-1)
          real*8 BremdEL    ! log10 step for 2D brem table at low energies
          real*8 BremUminLB ! min of uniform random number in the region B
                            !   0. sqrt(u)
          real*8  BremUmaxLB ! max. sqrt(BremUminLA)
          integer BremUszLB ! u talbe size in region B.  20
          real*8  BremdULB  ! step of u in B
!           --------------------------------------------- 

          real*8 PairEgmin  ! min. Eg above which pair cross section is
                            ! computed  1.1 MeV.   However,
                            ! at energies from  PairEgmin to PairNonSc, 
                            ! B.H original xsection is used as xsec=
                            ! Norm * B.H where
                            ! Norm * B.H(10MeV) = Pair(10MeV)
          real*8 PairNonSc  ! see above.  
          real*8 PairLEgmin ! log10 of PairEgmin  
          real*8 PairEgmaxL ! Eg where LPM effect starts to appear
          real*8 PrScrE     ! below this, screened cross-section is used.
          integer PairTXTL  ! Size of the Pair total x-section table.
                            ! in the energy region PairEgmin ~ PairEgmaxL;
                            ! log10(PairEgmaxL/PairEgemin)*10
          integer PairEsize ! Size of log10 energy for 2D  table for
                            ! pair in the region A,B.  PairTXTL/2
          real*8 PairUminLA !min of  uniform random number  in the region A
                            ! at energies PairEgmin ~ PairEgmaxL: 0.05
          real*8 PairUmaxLA ! max of  uniform random number in the region A
                            ! at energies PairEgmin ~ PairEgmaxL: 1.0
          integer PairUszLA ! Size  of uniform random nubmbers for 2D table
                            ! for pair in the region A.: 20

          real*8 PairdULA   ! step of u in region A  at Low energies
          real*8 PairdETXL  ! log10 step of total pair cross-sec at low E
          real*8 PairUminLB ! min of uniform random number in the region B
                            !   0. sqrt(u)
          real*8  PairUmaxLB ! max. sqrt(PairUminLA)
          integer PairUszLB ! u talbe size in region B.  20
          real*8  PairdULB  ! PairdULB=(PairUmaxLB-PairUminLB)/(PairUszLB-1))
          real*8  PairdELA  ! log10(PairEgmaxL/ PairEgmin) /( PairEsize-1)
          real*8  PairdELB  ! sqrt( log10(PairEgmaxL/ PairEgmin) )
                            !    /( PairEsize-1)
!
!            for Seltzer cross-section
          real*8  BrEeminS
          real*8  BrEgminS !  Eg min  for Seltzer Brems. 
          real*8  BrLEeminS
          real*8  BrEemaxS
          integer BrTXTS 
          integer BrES
          real*8  BrUminSA
          real*8  BrUmaxSA
          integer BrUszSA
          real*8  BrdUSA 
          real*8  BrdETXS
          real*8  BrdES
          integer BrUszSB
          real*8  BrUminSB
          real*8  BrUmaxSB 
          real*8  BrdUSB 



      common /Zrcnst/	
     *CompScrE, BrScrE, BremEgmin, BremEemin, BremLEemin, BremEemaxL, 
     *BremUminLA, BremUmaxLA, BremdULA, BremdETXL, BremdEL, BremUminLB, 
     *BremUmaxLB, BremdULB, PairEgmin, PairNonSc, PairLEgmin, 
     *PairEgmaxL, PrScrE, PairUminLA, PairUmaxLA, PairdULA, PairdETXL, 
     *PairUminLB, PairUmaxLB, PairdULB, PairdELA, PairdELB, BrEeminS, 
     *BrEgminS, BrLEeminS, BrEemaxS, BrUminSA, BrUmaxSA, BrdUSA, 
     *BrdETXS, BrdES, BrUminSB, BrUmaxSB, BrdUSB,
     *cScrC1, cScrC2, cScrMain, mbtoPX0, X0g 
!     * , muNVmin, muNdU, muNEmin, muNLEmin, muNEmax, 
!     *muNEmax1, muNdETX, muNdE, muNpwtx, muNpwdEdx0, muNpwdEdxt, 
!     *muBrVmin, muBrdU, muBrEmin, muBrLEmin, muBrEmax, muBrEmax1, 
!     *muBrdETX, muBrdE, muPrVmin, muPrdU, muPrEmin, muPrLEmin, 
!     *muPrEmax, muPrEmax1, muPrdETX, muPrdE

      common /Zicnst/
     *BremTXTL, BremEsize, BremUszLA, BremUszLB, PairTXTL, PairEsize, 
     *PairUszLA, PairUszLB, BrTXTS, BrES, BrUszSA, BrUszSB
!     * , muNTXT, muNUsize
!     *, muNEsize, muBrTXT, muBrUsize, muBrEsize, muPrTXT, muPrUsize, 
!     *muPrEsize

