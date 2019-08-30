!            constants used for the brems and pair creation
!     
       type SmpCnst 
       sequence

          real*8 CompScrE   !  Energy above which we can use complete screening
                            !  cross-sections. evaluate at Eg/Ee=x= 0.99.

          real*8 BrScrE     !  below this, partial screened cross-section is 
                         ! needed ( = ComScrE)
          real*8 BremEgmin  !  min. ratio  Eg/Ee for brems at high energy region

          real*8 BremEemin  ! Below this, partial screening brems x-section is
                      !  not made. (Seltzer table 1 is used)
          real*8 BremLEemin ! log10 of BremEemin
!          real*8 BremEemaxL ! Energy of e+/e- below which LPM is negelected
	    ! EemaxL is renamed to EeminLPM  (original naming is confusing) v9.15
          real*8 BremEeminLPM ! Min. energy of e+/e- above which LPM can be applied,
	    ! if  wanted.(  Acutal application will be done if, Ee > Flpm*this
            ! and LPMeffect=T. Flpm and LPMeffect can be controled by epicsfile
            !
          integer BremTXTL  ! Size of the Brems total x-section table.
                            ! in the energy region BremEemin ~ BrScrE:
                            ! ~log10(BrScrE/BremEemin)*10
          integer BremEsize ! Size of log10 energy for 2D  table for
                            !  brems in the region A.  BremTXTL/2
          real*8 BremUminLA ! min of  uniform random number  in the region A
                            ! at energies BremEemin ~ BrScrE: 0.1
          real*8 BremUmaxLA ! max of  uniform random number in the region A
                            ! at energies BremEemin ~ BrScrE: 1.0
          integer BremUszLA ! Size  of uniform random nubmbers for 2D table
                            ! for brems in the region A.: 20
          real*8 BremdULA   ! step of u in region A  at Low energies
          real*8 BremdETXL   ! log10E  step for  brem total  cross secton
                            ! log10(BrScrE/BremEemin)/(BremTXTL-1)
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
!            for Seltzer cross-section;  lower energy region
          real*8  BrEeminS
          real*8  BrEgminS !  Eg min  for Seltzer Brems.  (not ratio: 1keV)
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

!            for Seltzer cross-section;  higher energy region upto 10 GeV
          real*8  BrEeminS2
          real*8  BrEgminS2 !  Eg/Ee  min  for Seltzer Brems. 
          real*8  BrLEeminS2
          real*8  BrEemaxS2
          integer BrTXTS2 
          integer BrES2
          real*8  BrUminSA2
          real*8  BrUmaxSA2
          integer BrUszSA2
          real*8  BrdUSA2 
          real*8  BrdETXS2
          real*8  BrdES2
          integer BrUszSB2
          real*8  BrUminSB2
          real*8  BrUmaxSB2 
          real*8  BrdUSB2 

!          -------------- for LPM
          real*8  BrEgminH
          real*8  BrEe1H
          real*8  BrLEe1H  ! log10( BrEe1H)
          integer BrneH 
          real*8  BrdU1H
          real*8  BrdEH    ! log E step

!       cnst.BrdEH= log10(cnst.BrEe2H/cnst.BrEe1H)/(cnst.BrneH-1)
!       inverse of the above
          real*8  BrEe2H  ! max Ee where table is available
                          ! 

          real*8 BrU1H
          real*8 BrU2H
          integer Brnu1H    !=(cnst.BrU2H-cnst.BrU1H+0.00001d0)/cnst.BrdU1H+1
          integer BrneH2     ! for 2D table E size
          real*8 BrdEH2     ! //  E bin
          real*8 BrEe2H2    ! max E for 2D table 
!           ....................
!                                 
          real*8 BrU3H
          real*8 BrU4H    ! 
          integer Brnu2H  ! 
          integer BrdVU2H !  = cnst.Brnu2H-1
          real*8  BrdU2H  !  = (cnst.BrU4H - cnst.BrU3H)/cnst.BrdVU2H
          real*8  BrPow

          real*8 PrEg1H   ! minimum Eg above which LPM works
          real*8 PrLEg1H  ! log10 of PrEg1H
          integer PrneH    ! number of Eg bins
          real*8 PrdU1H   ! du
          real*8 PrdEH    ! dE in log10(Eg)
          real*8 PrU1H    ! minimum u= 0
          real*8 PrU2H    ! maximum u= 1 
          integer Prnu1H   ! numboer of u bins
          real*8 PrEg2H   ! max Eg where table is available.
!
!          -----------muons 
!              nuclear interaction
          real*8 muNVmin   ! min of Eg(virtual)/Emu by muon nuc.  int.
          real*8 muNdU     ! du for sampling table
          integer muNTXT  ! total xs, dEdx(v<vmin), dEdx(vall), tab size.
          real*8  muNEmin  ! above this, muon nuc. int. is treatable
          real*8  muNLEmin ! log10 of muNEmin
          real*8  muNEmax  ! above this, use some scaling(sampling)
          real*8  muNEmax1 ! max E of 1D table 
          real*8  muNdETX  ! log10 Energy step for total muon nuc. int prob. 
          real*8  muNdE    ! log10 Energy step for sampling table
          integer muNUsize ! sampling table size for u.
          integer muNEsize ! sampling table size for log10 E 
          real*8  muNpwtx  ! prob/X0 energy dependence; power.  set after 
                           ! table for total prob. is read
          real*8  muNpwdEdx0 ! dEdx(v<vmin)/Emu enery dependence; power
                             ! set after table is read
          real*8  muNpwdEdxt ! dEdXt(v<vmax)/Emu energy dependence: power
                             ! set after table  is read
!              brems
          real*8 muBrVmin !  min of Eg/Emu. for muon Brems
          real*8 muBrdU     ! du for sampling table
          integer muBrTXT  ! total xs, dEdx(v<vmin), dEdx(vall), tab size.
          real*8  muBrEmin  ! above this, muon brems is treatable
          real*8  muBrLEmin ! log10 of muBrEmin
          real*8  muBrEmax  ! above this, use some scaling
          real*8  muBrEmax1 ! max E of 1D table
          real*8  muBrdETX  ! log10 Energy step for total muon brems prob. 
          real*8  muBrdE    ! log10 Energy step for sampling table
          integer muBrUsize ! sampling table size for u.
          integer muBrEsize ! sampling table size for log10 E 
!                dependence can be neglected
!          real*8  muBpwtx  ! prob/X0 energy dependence; power.  set after 
                           ! table for total prob. is read
!          real*8  muBpwdEdx0 ! dEdx(v<vmin)/Emu enery dependence; power
                             ! set after table is read
!          real*8  muBpwdEdxt ! dEdXt(v<vmax)/Emu energy dependence: power

!                  pair creation
          real*8  muPrVmin  ! min of Eg(virtual)/Emu by muon pair cre.
          real*8  muPrdU     ! du for sampling table
          integer muPrTXT  ! total xs, dEdx(v<vmin), dEdx(vall), tab size.
          real*8  muPrEmin  ! above this, muon pair creation is treatable
          real*8  muPrLEmin ! log10 of muPrEmin
          real*8  muPrEmax  ! above this, use some scaling
          real*8  muPrEmax1 ! max E of 1D table
          real*8  muPrdETX  ! log10 Energy step for total muon pair prob. 
          real*8  muPrdE    ! log10 Energy step for sampling table
          integer muPrUsize ! sampling table size for u.
          integer muPrEsize ! sampling table size for log10 E 
!              dependence can be neglected
!          real*8  muPpwtx  ! prob/X0 energy dependence; power.  set after 
                           ! table for total prob. is read
!          real*8  muPpwdEdx0 ! dEdx(v<vmin)/Emu enery dependence; power
                             ! set after table is read
!          real*8  muPpwdEdxt ! dEdXt(v<vmax)/Emu energy dependence: power

!             normalization const
           integer::how
           real(8)::NormS
           real(8)::NormPS
           real(8)::NormCS
           real(8)::NormSH

       end type SmpCnst
!#if defined IBMAIX
!        character*500 cline
!#else
!        character*600 cline
!#endif
!        common /bpSamplec/ cline

