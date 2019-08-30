!
!     this is not the main program. for C++ interface.
!          this is called from C++ user interface.
	block data cblkElemag


!         parameters for Elemag process.
!	(->   ----------------------------------------------

	real*8  RecoilKineMinE  !2  Recoil Kinetic Min Energy above which the recoil (=knock-on process)
                                ! is treated. Below this energy, the effect is included as continuous
                                ! energy loss.  Used only if KnockOnRatio $>$ 1.
                                ! See also KnockOnRatio.
        real*8 KnockOnRatio   !2  KnockOnRatio* KEminoObs is used instead of RecoilKineMinE if KnockOnRatio $<$1.
	real*8  X0            !2  Radiation length in kg/m$^2$ for air. Normally the user should not touch this.
	real*8  Ecrit         !2  Critical energy in GeV. \newline
                              !   Employed only when calculating air shower size in the hybrid 
                              !  air shower generation.  The value would be dependent on the
                              !  experimental purpose.  The default  value, 81 MeV, is bit too
                              !  small in many applications (The air shower size is overestimated). 
                              !  Comparisons of sizes  by the hybrid method and by the  full Monte 
                              !  Carlo tell that \newline 
                              !  $N_e$ (full 3-D M.C) $ < N_e$ (hybrid AS with $E_c=81$  MeV ) $  < N_e$ (full 1-D M.C)
                              !  $ {\ \lower-1.2pt\vbox{\hbox{\rlap{$<$}\lower5pt\vbox{\hbox{$\sim$}}}}\ }
                              !  N_e$(hybrid AS with $E_c={76}$ MeV)  at around shower maximum.
                              !  Hybrid AS is always essentially 1-D. 

	logical Knockon       !2  Obsolete. Don't use this. See RecoilKineMinE
                              !    and KnockonRatio.
	real*8  AnihiE        !2  If E(positron) $<$ AnihiE, annihilation is  considered. 
	real*8  Es            !2  Modified scattering constant. 19.3d-3 GeV
	real*8  MaxComptonE   !2  Above this energy, Compton scattering is neglected. 
	real*8  MaxPhotoE     !2  Above this energy, photoelectric effect is neglected.
	real*8  MinPhotoProdE !1  Below this energy, no photo-prod of hadron.  See also PhotoProd.
	logical PhotoProd     !1  Switch. if .false., no photo prod. of hadron is considered at all. 
                              !   See also MinPhotoProdE.
        integer Moliere       !2  0$\rightarrow$ use Gaussian approx always (with air density change and
			      !    energy loss effect)\newline
                              !   1$\rightarrow$ use Moli\`ere scattering for non-electrons (default)\newline
                              !   2$\rightarrow$ use Moli\`ere scattering for all charged particles.\newline
                              !    If negative, anglular-correlated displacement is made to be 0 since Moli\`ere
                              !    theory cannot give it. (if $>0$, we use Gaussian approximation for correlation). 


!	<-)	----------------------------------------------

	common /Zelemagc/  RecoilKineMinE, KnockOnRatio,
     *  AnihiE, MaxComptonE,
     *  MaxPhotoE, MinPhotoProdE,   Es, X0, Ecrit, Knockon, 
     *  PhotoProd, Moliere



	data 
     *	RecoilKineMinE /0.2d-3/,     !
     *  KnockOnRatio  /1.0d0/,       ! 
     *  AnihiE /30.e-3/,            ! Eposi < 30 MeV, anihilation considered
     *  X0 /365.667/,                ! radiation length of air in kg/m2
     *  Ecrit /81.e-3/,             ! critical energy of air in GeV
     *  MaxComptonE /1./,       ! compton is considered below 75 MeV
     *  MaxPhotoE /1.e-3/,            ! above this, PhotoElectric effect neg.
     *  MinPhotoProdE /150.e-3/,     ! below 150 MeV, no gp --> hadrons
     *	Es /19.3d-3/,                 ! scattering const. Note, not 21 MeV.
     *  Knockon /.true./,            ! knockon is considered. Obsolete
     *  PhotoProd /.false./,        ! gp--> hadrons not considered.
     *  Moliere /2/                 ! Moliere scattering for non-electron 
	end
      block data cblkEvhnp

!            hadronic collision parameters 
!	(->  ---------------------------------------------

        real*8  Cepic0  !2  Obsolute
        real*8  Cekaon  !2  Obsolute
        integer SucInt  !2  The number of successive interactions inside A is affected by this parameter.\newline
                        !   If $0\rightarrow$ old formula (before uv3.500) is used, which give rather         
                        !   smaller number ($<Nsuc>$ in Air = 1.7 for 30 mb pp), \newline
                        !   if $1\rightarrow$ new one $<Nsuc>=2.2$ for 30 mb pp). \newline
                        !   Default is 0 (from V5.00 again).
        real*8  Ceneuc  !2 \verb| p -> n ; n-> p; p~->n~; n~->p~| prob.
        real*8  Mudirp  !2 \verb| DD~| \# enhancement factor. D is only for prompt muon.
     	real*8  Kpilog  !2  K$_{ch}/\pi_{ch}$=(Kpilog*log(ss+.069)+Kpicns)*exp(-8/s') where ss(GeV**2)= 
                        !   effective s. s'(GeV**2)=s - 4.63.  See also Kpicns.
        real*8  Kpicns  !2  See Kpilog. 0.077
        real*8  Elund   !2   obsolete (from v6.0)
!                  old def
!        real*8  Elund   !2   E $<$ Elund$\Rightarrow$  Lund M.C is used. in h-n. (IntModel='int1'). 
!                        !    If IntModel='int2' $\Rightarrow$    Gheisha is used. Elund can be as
!                        !    small as 4.99. (Default is 500. For NEXT 4.99).
        real*8  Elund2  !2  obsolete (from v6.0)
!              old def
!        real*8  Elund2  !2   Elund $<$ E $<$  Elund2 $\Rightarrow$   ad hoc model is used; Elund2 must
!                        !    be {\ \lower-1.2pt\vbox{\hbox{\rlap{$>$}\lower5pt\vbox{\hbox{$\sim$}}}}c                        !     }100, if Elund2 $<$     Elund3.   Defualt is Elund
         real*8 Elund3   !2 obsolete  (from v6.0).
!                old def
!        real*8  Elund3  !2  Elund2$<$ E $<$  Elund3 $\Rightarrow$  New Lund Fritiof is used.
!                        !   Elund3 can be $<$     Elund2. If E $>$  Elund3, ad hoc model is used.
!                        !   Default is Elund.
	real*8  Efermi  !2  If Kinetic E $<$ Efermi, Fermi Momentum is considered for Nucleus target.\\

        character*64  IntModel !1 Interaction model description. Usage was changed from v6.0. 
                      !    One may list  code name and upper energy bound  for the code.\newline
                      ! E.g.  IntModel = '"dpmjet3"' ; to specify the dpmjet3 in the entire energy region
                      !       IntModel = '"dpmjet3" 100 "nexus2" to specify dpmjet3 at $<$ 100 GeV and nexus2 
                      !        at E$>$100 GeV \newline
                      !       IntModel = '"nucrin" 5 "fritiof1.6" 500 "adhoc" to specify  Nucrin,
                      !               fritiof1.6, and ad-hoc model in the respective energy region. This 
                      !       corresponds to the old IntModel='int1'. \newline
                      !       IntModel = '"nucrin" 5 "fritiof1.6" 10 "fritiof7.02"  and \newline
                      !       IntModel = '"dpmjet3"' \newline
                      !        are most promissing models that fit the observed 
                      !       data (muons and gamma rays) for which the primary is well known by
                      !       BESS and AMS observations ($<$ 100 GeV).

!              next is old def.
!        character*16  IntModel   !1 interaction model.  If 'int1', Lund Fritiof/Nucrin/adhoc models are used
!                       ! for hadronic interactions.  If 'int2' Gheisha code is used for hadronic interactions,
!                       ! at low energies. (See, Elund?).
!                       !  If 'inclusive', an inclusive treatment is employed as given in Contrib/Inclusive.
!                       !  To make this work actually, INCLUSIVE in cosmos/Zcondc.h must be given a value of 1
!                       !  when making the cosmos library.  At execution time, InclusiveFile must be given
!                       !  the path to the xdist.d file. (normally in Contrib/Inclusive/xdist.d)
!                       !  The treatmen is one dimesional. May be used for neutrino calculations.
        character*100 InclusiveFile !2 The path to the inclusive data file, xdist.d. Default is
                       !   "../Contrib/Inclusive/xdist.d"
        real*8 SucPw   !2  In the 2nd, 3rd,.. collision of a nucleon inside a nucleus, the collision is 
                       !  made to be more elastic than normal one. The leading particle spectrum is
                       ! sampled from x**SucPw dx. SucPw should be in 1 to 2.
        real*8  Eta2Pi0 !2  eta/pi0 ratio. this is used to see the effect due to non-decay of pi0
                      ! at very high energies. Only source of h.e gamma can be eta and LPM may work
                      ! for them. default is 0.2
        integer MulLow !2 if 1, QCD predicted multiplicity law is used in the adhoc model else UA5
                       !  parametalization is used. Default is 1. (from v5), 
                       !  0.6135exp(23/18sqrt(2log(roots/0.3))) is QCD jet prediction. 
                       !  7.2roots**0.254 -7 is UA5 data. The branch point is set at roots = 900 GeV. 
                       !  (I have adjusted 0.6135 so that 900 GeV is the b.p)
        integer LundPara !2 To control Lund program. LundPara(1) is set to kfr(7); kfr(7)=1 is for Frititof
                       !  hard scattering. 2 is for Pythia H.S. 2 gives higher multiplicity but shape is 
                       !  strange.  Default is 1. LundPara(2) is set to kfr(12): 1 by for OPAL hard scattering
                       !   parameterization. 2 by DELPHI. Default is 2. (2 gives bit higher PT). LundPara(3)
                       !   $>$ 0 $\Rightarrow$    Pythia message will appear. LundPara(4) $>$ 0 Fritiof
                       !  message; both on ErrorOut. LundPara(5) =0 $\Rightarrow$    All kaons collisions
                       !  are treated as pi- in Fritiof, else they are treated by adhoc model as they are.

!	<-)     -----------------------------
!            next are not input parameters.
        integer nmdls
        parameter (nmdls = 8)
        character*16 RegMdls(nmdls)

        common /Zevhnp/ Cepic0, Cekaon, Ceneuc, Mudirp, Kpilog, Eta2Pi0,
     *         Kpicns,  Efermi, Elund, Elund2, Elund3, SucPw, MulLow, 
     *         LundPara(10), SucInt
        common /Zevhnc/RegMdls, InclusiveFile, IntModel 
!                 this should be in block common. 
!        data  Cepic0/0.3d0/, Cekaon/0.d0/, Ceneuc/0.35d0/,
!    *         Mudirp/1.d0/, Kpilog/0.0062/, Kpicns/0.077/,
!    *         Efermi/2.d0/

!            currently usable models
      data RegMdls/'dpmjet3', 'nexus2', 'fritiof7.02', 'fritiof1.6',
     *             'gheisha', 'nucrin', 'ad-hoc', 'incdpm3'/
      
       data  
     * Cepic0 /0.00d0/ , 
     * Cekaon /0.d0/ ,
     * Ceneuc /0.35d0/ ,
     * Mudirp /1.d0/ , 
     * Kpilog /0.0062/ ,
     * Kpicns /0.077/ ,

     * Elund /500./ , 
     * Elund2 /500./,
     * Elund3 /500./,

     * SucPw /1.5/,
     * MulLow /1/,      ! old one is 0. bef. v.5.0
     * LundPara /1, 2, 0, 0, 1, 5*0/,
     * IntModel /'"dpmjet3" 1.d8'/ ,
     * InclusiveFile /'../../Contrib/Inclusive/xdist.d'/ ,
     * Efermi /10.e-3/
       data 
     * Eta2Pi0/0.2/,

     * SucInt /1/     ! old one was 1 < v5.0; 0 <=v5.10   1 >= v5.11

      end
!     ********************
      block data cblkHeavy
!     ********************

!             ptcl kind code; kindmx is the no. of observable ptcls
!             klast; max ptcl code in the system.
!
        integer  kphoton, kelec, kmuon, kpion,  kkaon, knuc,
     1  kneue, kneumu, kindmx, knnb, kddb, kdmes, krho,
     2  komega, kphi, keta, kgnuc, kalfa, klibe, kcno, khvy, kvhvy,
     3  kiron, khvymax, klast, klambda, ksigma, kgzai, kbomega,
     4  ktriton, klambdac, krare
!            subcode
        integer
     1  regptcl, antip, k0s,  k0l, kneutron, 
     4  kneutronb, kd0, kd0b, kdirectg, kcasg

!
        parameter(
     1  kphoton=1, kelec=2, kmuon=3, kpion=4,  kkaon=5,
     2  knuc=6,
     3  kneue=7, kneumu=8, kgnuc=9, kalfa=10, klibe=11, kcno=12, 
     4  khvy=13, kvhvy=14, kiron=15, kdmes=16, 
!          next line added Nov. 17,'95. 
     5  ktriton=17, klambda=18, ksigma=19, kgzai=20, klambdac=21,
     6  kbomega=22,  kindmx=kbomega, krare = 0,

     7  knnb=kindmx+1, kddb=knnb+1,  krho=kddb+1,
     8  komega=krho+1, kphi=komega+1, keta=kphi+1,
     9  klast=keta, khvymax = kiron)
        parameter(
     1  regptcl=-1, antip=1,
     2  kdirectg=2, kcasg=3,
     3  k0s = 4,  k0l= 5,
     4  kneutron=regptcl,
     5  kneutronb = antip, kd0 =-8, kd0b =-kd0)
!
	integer maxHeavyMassN, maxHeavyCharge, maxHeavyG
	parameter (maxHeavyMassN = 56, maxHeavyCharge = 26, 
     1             maxHeavyG = 7)
!       kphoton: gamma ray 
!        kelec: electron, positiron
!        kmuon: muon
!        kpion: pion
!        kkaon: kaon
!        knuc: neucleon
!        kneue: electron neutrino
!       kneumu: muon neutrino
!        kgnuc: general nucleus(A>=2.)
!        kalfa: alpha  (heliunm)
!        klibe: Li, Be, B
!         kcno: C, N, O 
!         khvy: heavy such as, Na/Mg/Si
!        kvhvy: very heavy such as S/Cl/Ar
!        kiron: iron group
!        regptcl: particle index
!        antip: anti-particle index
!        krare:  used to set very rare particle code
!                which might come from imported soft.
!                such as tau. They are neglected in
!                Cosmos.





!
!	(->	------------------------------------------

	integer Charge2heavyG	!2  charge of heavy $\rightarrow$  heavy group index conversion array.
        integer HeavyG2massN    !2  heavy group index $\rightarrow$     mass number conversion array.
	integer HeavyG2charge	!2  heavy group index $\rightarrow$     charge of heavy conversion array.
        integer HeavyG2code     !2  heavy group index $\rightarrow$     particle code conversion array.
        integer Code2massN      !2  particle code $\rightarrow$     mass number conversion array.
        integer Code2heavyG	!2  particle code $\rightarrow$     heavy group index conversion array.
        real*8  FragmentTbl	!2  tbl(i,j)=$<$Number$>$  of frag. j when a heavy of heavy group index i
                                !    breaks up at air.
        real*8  PtAvNonInteNuc  !2  $<$Pt$>$  of non interacting nucleons.
	real*8  PtAvFrag        !2  $<$Pt$>$  of heavy fragments.
	character*4 HeavyG2symbol !2   heavy group index $\rightarrow$  'Fe' etc conversion array.
	 integer HowIntNuc       !2 If 0, the  number of interacting nucleons among a projectile heavy nucleus is 
                                 !  determined as the number of first collision of each interacting nucleon inside 
                                ! the  nucleus.  If 1, the number is determined as the total number of collisions 
                                !   including successive interactions. Default is 1. (There is uncertaninity in
                                !  interpretation of the formula; value 1 gives larger number of interacting
                                !  nucleons.)


 
!	<-)	--------------------------------------
        

	common /Zheavyc/
     *   PtAvNonInteNuc, PtAvFrag,
     *   FragmentTbl(maxHeavyG, maxHeavyG), 
     *	 Charge2heavyG(maxHeavyCharge),
     *   HeavyG2massN(maxHeavyG), HeavyG2charge(maxHeavyG),
     *   HeavyG2code(maxHeavyG), Code2massN(khvymax),
     *   Code2heavyG(khvymax), HowIntNuc
        common /Zheavycc/ HeavyG2symbol(maxHeavyG)



	
      integer i, j

      data
     *  HeavyG2massN /1,   4,    8,    14,    25,    35,     56/ ,
     *  HeavyG2charge /1,   2,    4,     7,    12,    17,     26/ ,
     *  HeavyG2symbol /'p', 'Alfa',  'L', 'CNO', 'H', 'VH', 'Fe'/ ,
     *  Charge2heavyG /1,   2,   3,3,  4,4,4,  5*5,   5*6,   9*7/

      data 
     *  HeavyG2code /knuc, kalfa, klibe, kcno, khvy, kvhvy, kiron/    
      data
     * Code2heavyG(knuc)/1/ ,
     * Code2heavyG(kalfa)/2/ ,
     * Code2heavyG(klibe)/3/ ,
     * Code2heavyG(kcno)/4/ , 
     * Code2heavyG(khvy)/5/ ,
     * Code2heavyG(kvhvy)/6/ ,
     * Code2heavyG(kiron)/7/
      data
     * Code2massN(knuc)/1/ ,
     * Code2massN(kalfa)/4/ , 
     * Code2massN(klibe)/8/ ,
     * Code2massN(kcno) /14/,
     * Code2massN(khvy)/25/,
     * Code2massN(kvhvy)/35/,
     * Code2massN(kiron)/56/

      data  HowIntNuc/1/
!
      data ( ( FragmentTbl(i,j), j=1,maxHeavyG), i=1,maxHeavyG)/
     1             1.,    0.,    0.,    0.,    0.,    0.,   0.,
     2             4.,    0.,    0.,    0.,    0.,    0.,   0.,
     3             5.81,  .61,   .11,   0.,    0.,    0.,   0.,
     4             6.62,  .72,   .24,   .17,   0.,    0.,   0.,
     5             10.5,  .77,   .21,   .39,   .16,   0.,   0.,
     6             12.8,  1.17,  .17,   .20,   .42,   .06,  0.,
     7             18.1,  1.71,  .24,   .17,   .22,   .20,  .17/

	data PtAvNonInteNuc/50.e-3/   ! 50 MeV
	data PtAVFrag/50.e-3/     ! 50 MeV
      end
       block data cblkIncident
!                parameters for primary angle sampling
!	(->	---------------------------------------------------

           complex*16  CosZenith   !1  Range of cos(zenith angle). Say, (0.5, 1.0). Used when Za2ry is 'is' 
                                   !   If ObsPlane=3 (spherical), real(CosZenith) must be $>$0, and means
                                   !   the zenith angle range at the incident point (not in Exyz system).
                                   !   In that case, azimuth is 0 to 2pi.
           complex*16  Azimuth     !1  Range of azimuthal angle in deg. Say, (0, 45). Default is (0,360).
                                   !   Can be such as (300., 390.).  Used when Za1ry is 'is'\newline
                                   !   If ObsPlane=3 (spherical), this is used to show the half opening angle
                                   !   range where the primary injection position is uniformly distributed 
                                   !   on a sphere.  The center of the opening angle is (Latit, Longit, HeightOfInj).
                                   !   In this case, for the upper opening angle,  min( Imag(Azimuth),180.) is used.
           character*4 Za1ry       !1  Specify the primary angle sampling method by one of 'is', 'ps' or 'aps'.\newline
                        !   "is" is isotropic. The range is by CosZenith.\newline
                        !   "ps" is for point source (See also SourceDec)\newline
                        !   "aps" is around point source (See also SourceDec and  Ddelta) \newline
                        !   If ObsPlane=3 (spherical), this must be "is".
           real*8  SourceDec    !1  Source declination of point source.(deg)
	   real*8  Ddelta       !1  SourceDec $\pm$ Ddelta is the region for 'aps' (deg).
           real*8  HeightOfInj  !1  The vertical height of primary injection point (m).
                                !   If this is $<$ deepest obs. level and zeinth angle of primary is $< 0$, 
                                !   the primary is  assumed to be upgoing even if Reverse =0.
	                        !   NOTE: BorderHeightH must be given explicitly in this case.
	   real*8  OffsetHeight !2  The vertical offset  height from  the deepest detector. 
                                ! The  primary is directed to this height above the detector.
                                ! If ObsPlane is  3 (spherical), not used.

!	<-)	----------------------------------------------------

           common /Zincident/ Azimuth, CosZenith, SourceDec, Ddelta,
     *      HeightOfInj,  OffsetHeight
           common /Zincidentc/ Za1ry
         data 
     *   Azimuth /(0.d0, 360.d0)/,
     *   Za1ry /'is'/ ,
     *   CosZenith /(0.6d0, 1.d0)/ ,
     *   SourceDec /30./ ,
     *   Ddelta /5./,
     *   HeightOfInj /100.d3/ ,
     *   OffsetHeight /0./
       end
      block data  cblkManager


















!   Parameters   needed  for the Launcher.
!
!	(->	------------------------------------

	integer ErrorOut    !2 Error output logical  dev number.
	character*70  PrimaryFile  !1  Primary Spectrum data file (full or relative path)
	character*70  CutOffFile   !1  Geomagnetic cut-off file
	character*70  ContFile     !1  Job continuation information file  (full or relative path).
                                   !   default is "ContInfo".  This will be created when job
                                   !   is finished normally.
	character*70  GeomagFile   !2  IGRF or WMM file path which contains geomagnetic field expansion 
                                   !   coefficients.  Their format is the same one as given in their web 
                                   !   page.  If ' ' (default), Cosmos/Data/Geomag/igrf is used.
	character*70  SkeletonFile !1   Skeleton information file (full or relative path). created if Job =
                                   !    'skeleton'. Default is 'skeletonParam'.  This is the Namelist data
                                   !     referred by Cosmos automatically  if Job='flesh' is specified. For
                                   !     Job='flesh', you have to modify some part of  this file.
        character*70  DpmFile      !2  control card to specify the dpmjet execution conditions. If ' ',
	                           !   Cosmos/Data/DPM/atmos.inp is assumed.  
	character*10  Job          !1  What kind of job you are going to do.\newline
                                   !   =' ' (default).  nothing special.\newline
                                   !   ='skeleton'.  Makes skeleton. \newline
                                   !   ='flesh'. Flesh skeleton events.  See manual.\newline
                                   !   ='newskel'   \newline
                                   !   ='newflesh'  see manual. \newline
        character*70  SeedFile     !1   File to  contain the initial random numbers for those events to 
                                   !    which you want to flesh. You can create the file by calling
                                   !     cwriteSeed in a user hook routine (say, in chookEnEvent) at 
                                   !     skeleton making time. Default is 'Seed'.  For a normal run with
                                   !      Job=' ', if SeedFile is not ' ',  two integer initial random numbers
                                   !      and the event number are  automatically output on the speicfied disk file.
        integer       SeedFileDev  !2   logical device number of SeedFile.
	logical       Cont         !1  If T, continuation from a previous job is assumed. Contfile content is used.
	integer       InitRN       !1  Initial random number seed. 2 integers. If InitRN(1) $<$ 0, file dev  \# 14
                                   !    is  assumed to have  pairs of IR in each row, and they are read to
                                   !    initialize each event.  This feature is ignored when Job = 'flesh' or 
                                   !    'newflesh'. The \# 14 file should be opened by the user routine
                                   !    (chookBgRun). This is almost debug purpose.\newline
                                   !   If InitRn(2)$<$0, timer, hostname and process number are used for the 
                                   !    initialization.
	integer       EventNo      !2  cumulative event number counter.(excluding discarded ones due to cutoff).
	integer       EventsInTheRun !2  Counter for event number in the run. Internal use.
                                     !          (excluding discarded ones due to cutoff).
	integer       DestEventNo    !1 2 integers: Final event no. to be generated and events to be generated
                                     !  in the current run.  If negative, their absolute is used and counting 
                                     !  includes discarded ones due to rigidity cutoff.
                                     !  If DestEventNo(2)=0, DestEventNo(1) is used. If it is negative, only
                                     !  DestEventNo(2) is checked to see events in the current run. For the
                                     !  flux calculation, negative ones are better.
	logical       Hidden         !1  Make T, if hidden parameters are to be written.
	integer	      TempDev	   !2  Logical Dev. number for temporary disk use.
	integer       PrevEventNo  !2  The event number already finished.  System use for Cont job.
                                   !        (excluding discarded ones due to cutoff).
	character*8   DeadLine     !1  The dead line before which the job should terminate.
                                   !   Should be given like '10.11.15' which means the nearest 10th, 11 O'clock,
                                   !   15 min.  Not used if Within has non zero value.  
        integer       Within       !1  The job should end within this minutes from now.  Default is 99999.
                                   !   If 0 is given,  DeadLine is used.
        real*8        BaseTime     !1  Rough cpu time needed for completing one event (say, for protons, or
                                   !   gamma rays) with energy BaseErg.  The cpu time estimation is based on 
                                   !   A * ( E1ry par nucleon )**BasePower / BaseErg * BaseTime, where A is mass number
                                   !  (for nucleus; otherwise 1).
        real*8        BaseErg      !2  See BaseTime.  The default is  1000 (GeV).
        real*8        BasePower    !2  See BaseTime.   Default is 1.0
        character*100 UserHookc    !2  array size is MAX\_USERHOOKC(=5). Usage is left for the user. To get the i-th
                                   !   component, the use may 'call  cqUHookc(i, cv)' in the userHook routine, 
                                   !   where cv is a character variable to receive the data.
        real*8        UserHookr    !2  array size is MAX\_USERHOOKR(=10). Usage is left for the user. To get the i-th
                                   !   component, the use may 'call cqUHookr(i, rv)' in the userHook routine,
                                   !   where rv is a real*8 variable to receive the data.
        integer       UserHooki    !2  array size is MAX\_USERHOOKI(=10). Usage is left for the user.  To get the i-th
                                   !   component, the use may 'call ccqUHooki(i, iv)' in the userHook rouitne,
                                   !   where iv is an integer varialbe to receive the data.
        character*100 AtmosFile    !2  path to the atmospheric data as in 'Cosmos/Data/Atmos/stdatmos2.d'

        character*32  AtEnv        !2  If this is non blank, an environmental variable with that name is
                                   !   assumed to exist and Cosmos tries to get the value of that env variable.
                                   !   If the value is obtained, the \verb/@/ in \verb/@_/ or \verb/@./
                                   !   expressing a part of a file name is replaced by that value. 
                                   !   (default is blank and in that case the \verb/@/ is replaced by
                                   !    the host name where  the job runs.)

        character*32 SharpEnv      !2  If this is non blank, an environmental variable with that name is
                                   !   assumed to exist and Cosmos tries to get the value of that env variable.
                                   !   If the value is obtained, the \verb/#/ in \verb/#_/ or \verb/#./ 
                                   !   expressing a  part of a file name is replaced by that value. 
                                   !   (default is blank and in that case the \verb/#/ is replaced by
                                   !    the process number of the run).

        character*32 PercentEnv    !2  If this is non blank, an environmental variable with that name is
                                   !   assumed to exist and Cosmos tries to get the value of that env variable.
                                   !   If the value is obtained, the \verb/%/ in \verb/%_/ or \verb/%./ 
                                   !   expressing a  part of a file name is replaced by that value. 
                                   !   (default is blank and in that case the \verb/%/ is replaced by
                                   !    the USER name).


!	<-)	-------------------------------------
	common /Zmanagerpc/
     *  BaseTime,  BaseErg, BasePower, Within, UserHookr(10),
     *  ErrorOut, Cont, InitRN(2), UserHooki(10),
     *  EventsInTheRun, DestEventNo(2), Hidden, TempDev, 
     *  PrevEventNo, SeedFileDev, EventNo


       common /Zmanagerpc2/
     * UserHookc(5), PrimaryFile,
     * CutOffFile,  Job, ContFile, AtmosFile, GeomagFile,
     * SkeletonFile, SeedFile, DpmFile, DeadLine, SharpEnv,
     * PercentEnv, AtEnv
 

       data 
     *  Errorout /0/   ! sun4 ; others.  if 0 is ng, give by namelist

       data  
     * DestEventNo /1,0/ ,
     * Hidden /.false./ ,
     * TempDev /11/ ,
     * PrimaryFile /' '/ ,
     * Job /' '/ ,
     * ContFile /'ContInfo'/ , 
     * SkeletonFile /'SkeletonParam'/ ,
     * SeedFile /'Seed'/ ,
     * PrevEventNo /0/ ,
     * EventNo /0/, 
     * SeedFileDev /22/ ,
     * DeadLine /' '/ ,
     * Within /99999/ , 
     * BaseTime /10./ ,
     * BaseErg /1000./ ,
     * BasePower /1./
       data
     * CutOffFile /' '/ 
       data
     * GeomagFile /' '/
       data
     * UserHookc /5*' '/
       data
     * UserHooki /10*99999999/
       data
     * UserHookr /10*1.d65/
       data
     * AtmosFile /' '/
       data
     * DpmFile /' '/
       data
     * PercentEnv /' '/
       data
     * SharpEnv /' '/
       data
     * AtEnv /' '/
      end
      block data cblkMuInt
      integer i
!         muon interaction related variables
!

!          
      real*8 Zeff, Zeff3  ! Zeff**(1/3)
      real*8 muPrVmin, muPrdETX, muPrdE, muPrEmin,
     *       muPrEmax, muPrdU,  muPrEmax1
      integer    muPrUsize, muPrEsize,  muPrTXT
!
      real*8 muBrVmin, muBrdETX, muBrdE, muBrEmin, muBrEmax,
     *       muBrdU, muBrEmax1
      integer muBrUsize, muBrEsize, muBrTXT
!
      real*8 muNVmin, muNdETX, muNdE, muNEmin,
     *  muNEmax, muNdU, muNEmax1
      integer muNUsize, muNEsize, muNTXT
!
      real*8 muNpwtx, muNpwdEdx0, muNpwdEdxt

      real*8 MuPrTX(38), MuPrdEdx0(38), MuPrdEdxt(38)
      real*8 MuBrTX(36), MuBrdEdx0(36), MuBrdEdxt(36)
      real*8 MuNTX(34), MuNdEdx0(34), MuNdEdxt(34)

      real*8 MuNTbl(101, 17)

      real*8  mupa, mura, mupb, muqb, muAk, muAkm, muAkm2,
     *        muPointLike, muShadow, mulogf0


      real*8  muNLEmin, muPrLEmin, muBrLEmin
      common /muintc/ MuNTbl,  MuPrTX, MuPrdEdx0, MuPrdEdxt,
     *    MuBrTX, MuBrdEdx0, MuBrdEdxt, MuNTX, MuNdEdx0, MuNdEdxt,
     *    muNpwtx, muNpwdEdx0, muNpwdEdxt, muPrVmin, muPrdETX, 
     *    muPrdE, muPrEmin, muPrEmax, muPrdU,  muPrEmax1,
     *    muBrVmin, muBrdETX, muBrdE, muBrEmin, muBrEmax,
     *    muBrdU, muBrEmax1, muNVmin, muNdETX, muNdE, muNEmin,
     *    muNEmax, muNdU, muNEmax1,  muNLEmin, muPrLEmin, muBrLEmin,
     *    mupa, mura, mupb, muqb, muAk, muAkm, muAkm2,
     *    muPointLike, muShadow, mulogf0, Zeff, Zeff3,
!
     *    muPrUsize, muPrEsize,  muPrTXT, muBrUsize, muBrEsize,
     *    muBrTXT, muNUsize, muNEsize, muNTXT

!          Pair total X-sec. 
          data ( MuPrTX(i), i=  1,   38)/
     1  0.431132E-02, 0.473175E-02, 0.514887E-02, 0.556010E-02,
     2  0.596286E-02, 0.635461E-02, 0.673285E-02, 0.709522E-02,
     3  0.743955E-02, 0.776396E-02, 0.806695E-02, 0.834742E-02,
     4  0.860479E-02, 0.883891E-02, 0.905012E-02, 0.923915E-02,
     5  0.940706E-02, 0.955517E-02, 0.968498E-02, 0.979806E-02,
     6  0.989601E-02, 0.998043E-02, 0.100528E-01, 0.101146E-01,
     7  0.101672E-01, 0.102116E-01, 0.102491E-01, 0.102805E-01,
     8  0.103068E-01, 0.103286E-01, 0.103467E-01, 0.103616E-01,
     9  0.103739E-01, 0.103839E-01, 0.103920E-01, 0.103987E-01,
     a  0.104040E-01, 0.104083E-01                                      
     * /   

!            Pair  dE/dx(v<vmin)/E /(g/cm2)
          data ( MuPrdEdx0(i), i=  1,   38)/
     1  0.415696E-07, 0.508098E-07, 0.606769E-07, 0.710744E-07,
     2  0.819072E-07, 0.930817E-07, 0.104506E-06, 0.116091E-06,
     3  0.127747E-06, 0.139387E-06, 0.150927E-06, 0.162285E-06,
     4  0.173381E-06, 0.184145E-06, 0.194509E-06, 0.204417E-06,
     5  0.213820E-06, 0.222682E-06, 0.230976E-06, 0.238687E-06,
     6  0.245811E-06, 0.252353E-06, 0.258325E-06, 0.263746E-06,
     7  0.268643E-06, 0.273044E-06, 0.276982E-06, 0.280489E-06,
     8  0.283600E-06, 0.286349E-06, 0.288768E-06, 0.290890E-06,
     9  0.292744E-06, 0.294359E-06, 0.295761E-06, 0.296974E-06,
     a  0.298020E-06, 0.298920E-06
     * /
!
!        Pair; dE/dx/E (v<all)  /(g/cm2)
          data ( MuPrdEdxt(i), i=  1,   38)/
     1  0.664623E-06, 0.725450E-06, 0.786305E-06, 0.846786E-06,
     2  0.906509E-06, 0.965081E-06, 0.102213E-05, 0.107730E-05,
     3  0.113026E-05, 0.118071E-05, 0.122843E-05, 0.127320E-05,
     4  0.131491E-05, 0.135348E-05, 0.138889E-05, 0.142119E-05,
     5  0.145044E-05, 0.147678E-05, 0.150034E-05, 0.152131E-05,
     6  0.153987E-05, 0.155620E-05, 0.157052E-05, 0.158301E-05,
     7  0.159386E-05, 0.160325E-05, 0.161134E-05, 0.161830E-05,
     8  0.162426E-05, 0.162936E-05, 0.163370E-05, 0.163739E-05,
     9  0.164053E-05, 0.164318E-05, 0.164542E-05, 0.164731E-05,
     a  0.164890E-05, 0.165023E-05
     * /
!
!
!       Brem; Total X-sec
          data ( MuBrTX(i), i=  1,   36)/
     1  0.348421E-03, 0.353857E-03, 0.358933E-03, 0.363712E-03,
     2  0.368194E-03, 0.372385E-03, 0.376290E-03, 0.379917E-03,
     3  0.383272E-03, 0.386366E-03, 0.389207E-03, 0.391806E-03,
     4  0.394175E-03, 0.396325E-03, 0.398269E-03, 0.400021E-03,
     5  0.401591E-03, 0.402995E-03, 0.404244E-03, 0.405351E-03,
     6  0.406329E-03, 0.407189E-03, 0.407942E-03, 0.408600E-03,
     7  0.409173E-03, 0.409669E-03, 0.410098E-03, 0.410467E-03,
     8  0.410784E-03, 0.411055E-03, 0.411286E-03, 0.411483E-03,
     9  0.411650E-03, 0.411791E-03, 0.411910E-03, 0.412010E-03
     * /

!        Brem:  dE/dx/E (v<vmin)  /(g/cm2)
          data ( MuBrdEdx0(i), i=  1,   36)/
     1  0.178842E-08, 0.178961E-08, 0.179057E-08, 0.179133E-08,
     2  0.179194E-08, 0.179242E-08, 0.179281E-08, 0.179312E-08,
     3  0.179336E-08, 0.179356E-08, 0.179371E-08, 0.179384E-08,
     4  0.179393E-08, 0.179401E-08, 0.179407E-08, 0.179412E-08,
     5  0.179416E-08, 0.179419E-08, 0.179422E-08, 0.179424E-08,
     6  0.179425E-08, 0.179426E-08, 0.179428E-08, 0.179428E-08,
     7  0.179429E-08, 0.179429E-08, 0.179430E-08, 0.179430E-08,
     8  0.179430E-08, 0.179431E-08, 0.179431E-08, 0.179431E-08,
     9  0.179431E-08, 0.179431E-08, 0.179431E-08, 0.179431E-08
     * /

!       Brem: dE/dx/E(v<all) /(g/cm2)
          data ( MuBrdEdxt(i), i=  1,   36)/
     1  0.785146E-06, 0.816296E-06, 0.845259E-06, 0.873845E-06,
     2  0.901994E-06, 0.929644E-06, 0.956731E-06, 0.983189E-06,
     3  0.100895E-05, 0.103394E-05, 0.105811E-05, 0.108137E-05,
     4  0.110368E-05, 0.112497E-05, 0.114520E-05, 0.116433E-05,
     5  0.118231E-05, 0.119914E-05, 0.121480E-05, 0.122929E-05,
     6  0.124263E-05, 0.125484E-05, 0.126595E-05, 0.127601E-05,
     7  0.128506E-05, 0.129317E-05, 0.130038E-05, 0.130677E-05,
     8  0.131240E-05, 0.131734E-05, 0.132165E-05, 0.132540E-05,
     9  0.132864E-05, 0.133143E-05, 0.133383E-05, 0.133588E-05
     * /

!       Nuc.Int Total X-sec
          data ( MuNTX(i), i=  1,   34)/
     1  0.131090E-03, 0.144395E-03, 0.158387E-03, 0.172382E-03,
     2  0.187685E-03, 0.202425E-03, 0.218363E-03, 0.234279E-03,
     3  0.251198E-03, 0.268982E-03, 0.286553E-03, 0.305978E-03,
     4  0.324348E-03, 0.335809E-03, 0.325788E-03, 0.318789E-03,
     5  0.316786E-03, 0.314891E-03, 0.313318E-03, 0.312124E-03,
     6  0.311322E-03, 0.310982E-03, 0.311065E-03, 0.311680E-03,
     7  0.312847E-03, 0.314774E-03, 0.317242E-03, 0.319977E-03,
     8  0.323055E-03, 0.326491E-03, 0.330301E-03, 0.334529E-03,
     9  0.339117E-03, 0.343873E-03
     * /

!        Nuc.Int. dE/dx/E (v<vmin) /(g/cm2)
          data ( MuNdEdx0(i), i=  1,   34)/
     1  0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00,
     2  0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00,
     3  0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00,
     4  0.000000E+00, 0.220165E-09, 0.915062E-09, 0.139798E-08,
     5  0.167378E-08, 0.189541E-08, 0.207196E-08, 0.220734E-08,
     6  0.230651E-08, 0.238182E-08, 0.243526E-08, 0.247029E-08,
     7  0.248771E-08, 0.248822E-08, 0.248032E-08, 0.247220E-08,
     8  0.246348E-08, 0.245385E-08, 0.244290E-08, 0.242993E-08,
     9  0.241769E-08, 0.241067E-08
     * /

!        Nuc.Int; dE/dx/E /(g/cm2)
          data ( MuNdEdxt(i), i=  1,   34)/
     1  0.402333E-06, 0.402855E-06, 0.402707E-06, 0.402319E-06,
     2  0.401983E-06, 0.401883E-06, 0.402184E-06, 0.403040E-06,
     3  0.404575E-06, 0.406782E-06, 0.409625E-06, 0.413063E-06,
     4  0.417039E-06, 0.421500E-06, 0.426368E-06, 0.431588E-06,
     5  0.437116E-06, 0.442897E-06, 0.448886E-06, 0.455043E-06,
     6  0.461337E-06, 0.467737E-06, 0.474219E-06, 0.480764E-06,
     7  0.487353E-06, 0.493973E-06, 0.500611E-06, 0.507257E-06,
     8  0.513904E-06, 0.520544E-06, 0.527173E-06, 0.533785E-06,
     9  0.540378E-06, 0.546948E-06
     * /

!
!      Nuc. Int. sampling table.
!      
          data ( MuNTbl(i,   1), i=  1,   68)/
     1 -0.172472E+01,-0.145269E+01,-0.125199E+01,-0.111430E+01,
     2 -0.101163E+01,-0.930870E+00,-0.865230E+00,-0.810390E+00,
     3 -0.763540E+00,-0.722800E+00,-0.686860E+00,-0.654790E+00,
     4 -0.625860E+00,-0.599560E+00,-0.575470E+00,-0.553250E+00,
     5 -0.532640E+00,-0.513430E+00,-0.495450E+00,-0.478620E+00,
     6 -0.462850E+00,-0.448010E+00,-0.434010E+00,-0.420760E+00,
     7 -0.408200E+00,-0.396260E+00,-0.384890E+00,-0.374040E+00,
     8 -0.363660E+00,-0.353720E+00,-0.344190E+00,-0.335030E+00,
     9 -0.326220E+00,-0.317730E+00,-0.309550E+00,-0.301640E+00,
     a -0.294000E+00,-0.286610E+00,-0.279440E+00,-0.272500E+00,
     b -0.265760E+00,-0.259210E+00,-0.252850E+00,-0.246650E+00,
     c -0.240630E+00,-0.234750E+00,-0.229020E+00,-0.223430E+00,
     d -0.217980E+00,-0.212650E+00,-0.207430E+00,-0.202330E+00,
     e -0.197340E+00,-0.192450E+00,-0.187660E+00,-0.182970E+00,
     f -0.178360E+00,-0.173830E+00,-0.169390E+00,-0.165030E+00,
     g -0.160740E+00,-0.156510E+00,-0.152360E+00,-0.148270E+00,
     h -0.144240E+00,-0.140260E+00,-0.136350E+00,-0.132480E+00          
     * /   
          data ( MuNTbl(i,   1), i=  69,   101)/
     1 -0.128660E+00,-0.124890E+00,-0.121160E+00,-0.117480E+00,
     2 -0.113830E+00,-0.110220E+00,-0.106640E+00,-0.103090E+00,
     3 -0.995700E-01,-0.960700E-01,-0.926000E-01,-0.891400E-01,
     4 -0.857000E-01,-0.822700E-01,-0.788500E-01,-0.754400E-01,
     5 -0.720200E-01,-0.686100E-01,-0.651800E-01,-0.617300E-01,
     6 -0.582600E-01,-0.547600E-01,-0.512200E-01,-0.476200E-01,
     7 -0.439500E-01,-0.401800E-01,-0.363100E-01,-0.321800E-01,
     8 -0.279000E-01,-0.225300E-01,-0.167700E-01,-0.132000E-01,
     9 -0.113000E-02                                                    
     * /   
          data ( MuNTbl(i,   2), i=  1,   68)/
     1 -0.192472E+01,-0.161718E+01,-0.139522E+01,-0.124252E+01,
     2 -0.112845E+01,-0.103907E+01,-0.966380E+00,-0.905490E+00,
     3 -0.853330E+00,-0.807830E+00,-0.767530E+00,-0.731410E+00,
     4 -0.698690E+00,-0.668970E+00,-0.641900E+00,-0.617080E+00,
     5 -0.594210E+00,-0.573010E+00,-0.553280E+00,-0.534850E+00,
     6 -0.517560E+00,-0.501300E+00,-0.485950E+00,-0.471420E+00,
     7 -0.457640E+00,-0.444530E+00,-0.432040E+00,-0.420110E+00,
     8 -0.408700E+00,-0.397770E+00,-0.387270E+00,-0.377180E+00,
     9 -0.367470E+00,-0.358110E+00,-0.349070E+00,-0.340340E+00,
     a -0.331890E+00,-0.323710E+00,-0.315770E+00,-0.308080E+00,
     b -0.300600E+00,-0.293330E+00,-0.286260E+00,-0.279370E+00,
     c -0.272660E+00,-0.266120E+00,-0.259730E+00,-0.253490E+00,
     d -0.247400E+00,-0.241430E+00,-0.235600E+00,-0.229890E+00,
     e -0.224290E+00,-0.218800E+00,-0.213420E+00,-0.208140E+00,
     f -0.202950E+00,-0.197860E+00,-0.192850E+00,-0.187930E+00,
     g -0.183080E+00,-0.178310E+00,-0.173600E+00,-0.168970E+00,
     h -0.164400E+00,-0.159900E+00,-0.155450E+00,-0.151050E+00          
     * /   
          data ( MuNTbl(i,   2), i=  69,   101)/
     1 -0.146710E+00,-0.142420E+00,-0.138170E+00,-0.133970E+00,
     2 -0.129810E+00,-0.125690E+00,-0.121600E+00,-0.117540E+00,
     3 -0.113510E+00,-0.109510E+00,-0.105530E+00,-0.101570E+00,
     4 -0.976200E-01,-0.936900E-01,-0.897600E-01,-0.858400E-01,
     5 -0.819100E-01,-0.779900E-01,-0.740500E-01,-0.701000E-01,
     6 -0.661200E-01,-0.621000E-01,-0.580400E-01,-0.539200E-01,
     7 -0.497200E-01,-0.454300E-01,-0.409800E-01,-0.363800E-01,
     8 -0.314100E-01,-0.262600E-01,-0.184200E-01,-0.141700E-01,
     9 -0.710000E-03                                                    
     * /   
          data ( MuNTbl(i,   3), i=  1,   68)/
     1 -0.212472E+01,-0.178196E+01,-0.153838E+01,-0.137024E+01,
     2 -0.124475E+01,-0.114634E+01,-0.106600E+01,-0.998400E+00,
     3 -0.940160E+00,-0.889050E+00,-0.844020E+00,-0.804090E+00,
     4 -0.768300E+00,-0.735960E+00,-0.706500E+00,-0.679500E+00,
     5 -0.654610E+00,-0.631550E+00,-0.610080E+00,-0.590020E+00,
     6 -0.571190E+00,-0.553470E+00,-0.536740E+00,-0.520910E+00,
     7 -0.505870E+00,-0.491570E+00,-0.477930E+00,-0.464890E+00,
     8 -0.452420E+00,-0.440450E+00,-0.428960E+00,-0.417900E+00,
     9 -0.407250E+00,-0.396970E+00,-0.387040E+00,-0.377430E+00,
     a -0.368130E+00,-0.359120E+00,-0.350370E+00,-0.341870E+00,
     b -0.333600E+00,-0.325560E+00,-0.317730E+00,-0.310090E+00,
     c -0.302640E+00,-0.295360E+00,-0.288250E+00,-0.281310E+00,
     d -0.274510E+00,-0.267870E+00,-0.261370E+00,-0.255010E+00,
     e -0.248780E+00,-0.242660E+00,-0.236670E+00,-0.230790E+00,
     f -0.225010E+00,-0.219330E+00,-0.213750E+00,-0.208260E+00,
     g -0.202860E+00,-0.197550E+00,-0.192310E+00,-0.187150E+00,
     h -0.182060E+00,-0.177040E+00,-0.172090E+00,-0.167190E+00          
     * /   
          data ( MuNTbl(i,   3), i=  69,   101)/
     1 -0.162360E+00,-0.157580E+00,-0.152850E+00,-0.148180E+00,
     2 -0.143540E+00,-0.138960E+00,-0.134410E+00,-0.129890E+00,
     3 -0.125410E+00,-0.120960E+00,-0.116540E+00,-0.112140E+00,
     4 -0.107760E+00,-0.103390E+00,-0.990400E-01,-0.946900E-01,
     5 -0.903400E-01,-0.859900E-01,-0.816200E-01,-0.772400E-01,
     6 -0.728300E-01,-0.683800E-01,-0.638800E-01,-0.593100E-01,
     7 -0.546700E-01,-0.499000E-01,-0.450000E-01,-0.398900E-01,
     8 -0.345400E-01,-0.286700E-01,-0.206700E-01,-0.147500E-01,
     9 -0.450000E-03                                                    
     * /   
          data ( MuNTbl(i,   4), i=  1,   68)/
     1 -0.232472E+01,-0.194684E+01,-0.168123E+01,-0.149722E+01,
     2 -0.135994E+01,-0.125171E+01,-0.116276E+01,-0.108731E+01,
     3 -0.102278E+01,-0.967050E+00,-0.918210E+00,-0.874890E+00,
     4 -0.836080E+00,-0.800980E+00,-0.769020E+00,-0.739710E+00,
     5 -0.712670E+00,-0.687610E+00,-0.664270E+00,-0.642440E+00,
     6 -0.621950E+00,-0.602640E+00,-0.584400E+00,-0.567110E+00,
     7 -0.550680E+00,-0.535040E+00,-0.520100E+00,-0.505820E+00,
     8 -0.492130E+00,-0.478980E+00,-0.466370E+00,-0.454250E+00,
     9 -0.442580E+00,-0.431330E+00,-0.420480E+00,-0.409990E+00,
     a -0.399850E+00,-0.390020E+00,-0.380490E+00,-0.371250E+00,
     b -0.362270E+00,-0.353540E+00,-0.345040E+00,-0.336760E+00,
     c -0.328690E+00,-0.320820E+00,-0.313140E+00,-0.305630E+00,
     d -0.298300E+00,-0.291120E+00,-0.284090E+00,-0.277210E+00,
     e -0.270470E+00,-0.263850E+00,-0.257360E+00,-0.250990E+00,
     f -0.244730E+00,-0.238580E+00,-0.232530E+00,-0.226580E+00,
     g -0.220720E+00,-0.214950E+00,-0.209260E+00,-0.203660E+00,
     h -0.198130E+00,-0.192670E+00,-0.187290E+00,-0.181960E+00          
     * /   
          data ( MuNTbl(i,   4), i=  69,   101)/
     1 -0.176700E+00,-0.171500E+00,-0.166360E+00,-0.161260E+00,
     2 -0.156220E+00,-0.151210E+00,-0.146260E+00,-0.141340E+00,
     3 -0.136450E+00,-0.131600E+00,-0.126770E+00,-0.121970E+00,
     4 -0.117190E+00,-0.112420E+00,-0.107660E+00,-0.102910E+00,
     5 -0.981700E-01,-0.934100E-01,-0.886400E-01,-0.838600E-01,
     6 -0.790400E-01,-0.741800E-01,-0.692700E-01,-0.642900E-01,
     7 -0.592200E-01,-0.540300E-01,-0.486900E-01,-0.431300E-01,
     8 -0.373000E-01,-0.309300E-01,-0.235700E-01,-0.151500E-01,
     9 -0.280000E-03                                                    
     * /   
          data ( MuNTbl(i,   5), i=  1,   68)/
     1 -0.252472E+01,-0.211161E+01,-0.182344E+01,-0.162318E+01,
     2 -0.147320E+01,-0.135396E+01,-0.125521E+01,-0.117290E+01,
     3 -0.110307E+01,-0.104274E+01,-0.989840E+00,-0.942900E+00,
     4 -0.900800E+00,-0.862710E+00,-0.827980E+00,-0.796100E+00,
     5 -0.766670E+00,-0.739340E+00,-0.713860E+00,-0.689990E+00,
     6 -0.667580E+00,-0.646520E+00,-0.626670E+00,-0.607900E+00,
     7 -0.590120E+00,-0.573220E+00,-0.557120E+00,-0.541770E+00,
     8 -0.527090E+00,-0.513020E+00,-0.499530E+00,-0.486560E+00,
     9 -0.474080E+00,-0.462060E+00,-0.450450E+00,-0.439240E+00,
     a -0.428390E+00,-0.417880E+00,-0.407700E+00,-0.397810E+00,
     b -0.388210E+00,-0.378870E+00,-0.369780E+00,-0.360930E+00,
     c -0.352300E+00,-0.343890E+00,-0.335670E+00,-0.327640E+00,
     d -0.319790E+00,-0.312110E+00,-0.304590E+00,-0.297220E+00,
     e -0.290000E+00,-0.282920E+00,-0.275970E+00,-0.269140E+00,
     f -0.262440E+00,-0.255850E+00,-0.249370E+00,-0.242990E+00,
     g -0.236710E+00,-0.230520E+00,-0.224420E+00,-0.218410E+00,
     h -0.212480E+00,-0.206620E+00,-0.200840E+00,-0.195130E+00          
     * /   
          data ( MuNTbl(i,   5), i=  69,   101)/
     1 -0.189480E+00,-0.183900E+00,-0.178370E+00,-0.172900E+00,
     2 -0.167470E+00,-0.162100E+00,-0.156770E+00,-0.151480E+00,
     3 -0.146230E+00,-0.141010E+00,-0.135830E+00,-0.130660E+00,
     4 -0.125520E+00,-0.120400E+00,-0.115280E+00,-0.110180E+00,
     5 -0.105070E+00,-0.999600E-01,-0.948400E-01,-0.896900E-01,
     6 -0.845200E-01,-0.793000E-01,-0.740300E-01,-0.686800E-01,
     7 -0.632400E-01,-0.576700E-01,-0.519400E-01,-0.459900E-01,
     8 -0.397200E-01,-0.329800E-01,-0.256100E-01,-0.154800E-01,
     9 -0.180000E-03                                                    
     * /   
          data ( MuNTbl(i,   6), i=  1,   68)/
     1 -0.272472E+01,-0.227609E+01,-0.196461E+01,-0.174754E+01,
     2 -0.158359E+01,-0.145190E+01,-0.134483E+01,-0.125594E+01,
     3 -0.118044E+01,-0.111513E+01,-0.105780E+01,-0.100684E+01,
     4 -0.961060E+00,-0.919560E+00,-0.881650E+00,-0.846920E+00,
     5 -0.814980E+00,-0.785470E+00,-0.758060E+00,-0.732510E+00,
     6 -0.708590E+00,-0.686130E+00,-0.664970E+00,-0.644970E+00,
     7 -0.626030E+00,-0.608040E+00,-0.590910E+00,-0.574580E+00,
     8 -0.558970E+00,-0.544030E+00,-0.529690E+00,-0.515920E+00,
     9 -0.502680E+00,-0.489910E+00,-0.477600E+00,-0.465710E+00,
     a -0.454210E+00,-0.443070E+00,-0.432280E+00,-0.421810E+00,
     b -0.411640E+00,-0.401750E+00,-0.392140E+00,-0.382770E+00,
     c -0.373640E+00,-0.364740E+00,-0.356040E+00,-0.347550E+00,
     d -0.339250E+00,-0.331140E+00,-0.323190E+00,-0.315410E+00,
     e -0.307780E+00,-0.300300E+00,-0.292960E+00,-0.285750E+00,
     f -0.278680E+00,-0.271720E+00,-0.264870E+00,-0.258130E+00,
     g -0.251500E+00,-0.244960E+00,-0.238520E+00,-0.232160E+00,
     h -0.225890E+00,-0.219690E+00,-0.213580E+00,-0.207530E+00          
     * /   
          data ( MuNTbl(i,   6), i=  69,   101)/
     1 -0.201550E+00,-0.195630E+00,-0.189780E+00,-0.183980E+00,
     2 -0.178230E+00,-0.172530E+00,-0.166870E+00,-0.161260E+00,
     3 -0.155680E+00,-0.150140E+00,-0.144620E+00,-0.139130E+00,
     4 -0.133660E+00,-0.128210E+00,-0.122770E+00,-0.117330E+00,
     5 -0.111890E+00,-0.106440E+00,-0.100980E+00,-0.954900E-01,
     6 -0.899700E-01,-0.844000E-01,-0.787700E-01,-0.730600E-01,
     7 -0.672400E-01,-0.613000E-01,-0.551800E-01,-0.488300E-01,
     8 -0.421600E-01,-0.350400E-01,-0.270700E-01,-0.158900E-01,
     9 -0.110000E-03                                                    
     * /   
          data ( MuNTbl(i,   7), i=  1,   68)/
     1 -0.292472E+01,-0.244015E+01,-0.210466E+01,-0.186977E+01,
     2 -0.169000E+01,-0.154712E+01,-0.143174E+01,-0.133577E+01,
     3 -0.125410E+01,-0.118330E+01,-0.112098E+01,-0.106546E+01,
     4 -0.101583E+01,-0.971160E+00,-0.930630E+00,-0.893620E+00,
     5 -0.859620E+00,-0.828230E+00,-0.799110E+00,-0.771980E+00,
     6 -0.746620E+00,-0.722810E+00,-0.700410E+00,-0.679260E+00,
     7 -0.659240E+00,-0.640240E+00,-0.622170E+00,-0.604960E+00,
     8 -0.588520E+00,-0.572790E+00,-0.557710E+00,-0.543250E+00,
     9 -0.529340E+00,-0.515950E+00,-0.503050E+00,-0.490590E+00,
     a -0.478550E+00,-0.466900E+00,-0.455620E+00,-0.444660E+00,
     b -0.434030E+00,-0.423690E+00,-0.413630E+00,-0.403840E+00,
     c -0.394290E+00,-0.384970E+00,-0.375880E+00,-0.366990E+00,
     d -0.358310E+00,-0.349800E+00,-0.341480E+00,-0.333320E+00,
     e -0.325330E+00,-0.317480E+00,-0.309780E+00,-0.302210E+00,
     f -0.294770E+00,-0.287460E+00,-0.280260E+00,-0.273180E+00,
     g -0.266200E+00,-0.259320E+00,-0.252530E+00,-0.245840E+00,
     h -0.239230E+00,-0.232700E+00,-0.226250E+00,-0.219870E+00          
     * /   
          data ( MuNTbl(i,   7), i=  69,   101)/
     1 -0.213560E+00,-0.207320E+00,-0.201130E+00,-0.195010E+00,
     2 -0.188930E+00,-0.182910E+00,-0.176930E+00,-0.170990E+00,
     3 -0.165080E+00,-0.159210E+00,-0.153370E+00,-0.147550E+00,
     4 -0.141750E+00,-0.135970E+00,-0.130200E+00,-0.124430E+00,
     5 -0.118650E+00,-0.112870E+00,-0.107060E+00,-0.101230E+00,
     6 -0.953700E-01,-0.894500E-01,-0.834600E-01,-0.774000E-01,
     7 -0.712200E-01,-0.649000E-01,-0.584000E-01,-0.516500E-01,
     8 -0.445800E-01,-0.370000E-01,-0.285400E-01,-0.164300E-01,
     9 -0.700000E-04                                                    
     * /   
          data ( MuNTbl(i,   8), i=  1,   68)/
     1 -0.300000E+01,-0.246289E+01,-0.211859E+01,-0.187312E+01,
     2 -0.168713E+01,-0.154212E+01,-0.142447E+01,-0.132610E+01,
     3 -0.124206E+01,-0.116974E+01,-0.110669E+01,-0.105107E+01,
     4 -0.100149E+01,-0.956920E+00,-0.916540E+00,-0.879710E+00,
     5 -0.845940E+00,-0.814790E+00,-0.785940E+00,-0.759110E+00,
     6 -0.734050E+00,-0.710570E+00,-0.688510E+00,-0.667700E+00,
     7 -0.648020E+00,-0.629360E+00,-0.611610E+00,-0.594710E+00,
     8 -0.578560E+00,-0.563120E+00,-0.548320E+00,-0.534110E+00,
     9 -0.520440E+00,-0.507280E+00,-0.494600E+00,-0.482340E+00,
     a -0.470500E+00,-0.459030E+00,-0.447920E+00,-0.437140E+00,
     b -0.426670E+00,-0.416490E+00,-0.406590E+00,-0.396950E+00,
     c -0.387550E+00,-0.378380E+00,-0.369430E+00,-0.360680E+00,
     d -0.352130E+00,-0.343760E+00,-0.335560E+00,-0.327540E+00,
     e -0.319660E+00,-0.311940E+00,-0.304360E+00,-0.296910E+00,
     f -0.289590E+00,-0.282400E+00,-0.275310E+00,-0.268340E+00,
     g -0.261470E+00,-0.254700E+00,-0.248030E+00,-0.241440E+00,
     h -0.234940E+00,-0.228520E+00,-0.222180E+00,-0.215910E+00          
     * /   
          data ( MuNTbl(i,   8), i=  69,   101)/
     1 -0.209710E+00,-0.203570E+00,-0.197490E+00,-0.191460E+00,
     2 -0.185490E+00,-0.179570E+00,-0.173690E+00,-0.167860E+00,
     3 -0.162060E+00,-0.156290E+00,-0.150550E+00,-0.144840E+00,
     4 -0.139140E+00,-0.133460E+00,-0.127790E+00,-0.122130E+00,
     5 -0.116460E+00,-0.110780E+00,-0.105090E+00,-0.993600E-01,
     6 -0.936100E-01,-0.878000E-01,-0.819200E-01,-0.759700E-01,
     7 -0.699100E-01,-0.637100E-01,-0.573300E-01,-0.507100E-01,
     8 -0.437600E-01,-0.363100E-01,-0.277500E-01,-0.187300E-01,
     9 -0.400000E-04                                                    
     * /   
          data ( MuNTbl(i,   9), i=  1,   68)/
     1 -0.300000E+01,-0.243876E+01,-0.208535E+01,-0.183439E+01,
     2 -0.164808E+01,-0.150166E+01,-0.138275E+01,-0.128472E+01,
     3 -0.120212E+01,-0.113132E+01,-0.106974E+01,-0.101555E+01,
     4 -0.967370E+00,-0.924170E+00,-0.885140E+00,-0.849610E+00,
     5 -0.817050E+00,-0.787040E+00,-0.759240E+00,-0.733370E+00,
     6 -0.709210E+00,-0.686550E+00,-0.665230E+00,-0.645120E+00,
     7 -0.626090E+00,-0.608030E+00,-0.590860E+00,-0.574500E+00,
     8 -0.558870E+00,-0.543920E+00,-0.529590E+00,-0.515830E+00,
     9 -0.502600E+00,-0.489850E+00,-0.477560E+00,-0.465690E+00,
     a -0.454220E+00,-0.443110E+00,-0.432340E+00,-0.421900E+00,
     b -0.411760E+00,-0.401900E+00,-0.392310E+00,-0.382960E+00,
     c -0.373860E+00,-0.364980E+00,-0.356310E+00,-0.347830E+00,
     d -0.339550E+00,-0.331450E+00,-0.323520E+00,-0.315740E+00,
     e -0.308120E+00,-0.300650E+00,-0.293310E+00,-0.286110E+00,
     f -0.279030E+00,-0.272060E+00,-0.265220E+00,-0.258470E+00,
     g -0.251830E+00,-0.245290E+00,-0.238840E+00,-0.232480E+00,
     h -0.226200E+00,-0.220000E+00,-0.213880E+00,-0.207820E+00          
     * /   
          data ( MuNTbl(i,   9), i=  69,   101)/
     1 -0.201830E+00,-0.195910E+00,-0.190050E+00,-0.184240E+00,
     2 -0.178480E+00,-0.172770E+00,-0.167100E+00,-0.161480E+00,
     3 -0.155890E+00,-0.150340E+00,-0.144810E+00,-0.139310E+00,
     4 -0.133830E+00,-0.128360E+00,-0.122910E+00,-0.117460E+00,
     5 -0.112010E+00,-0.106550E+00,-0.101070E+00,-0.955700E-01,
     6 -0.900400E-01,-0.844600E-01,-0.788100E-01,-0.730900E-01,
     7 -0.672600E-01,-0.613100E-01,-0.551800E-01,-0.488200E-01,
     8 -0.421400E-01,-0.349700E-01,-0.265800E-01,-0.181800E-01,
     9 -0.300000E-04                                                    
     * /   
          data ( MuNTbl(i,   10), i=  1,   68)/
     1 -0.300000E+01,-0.242609E+01,-0.205915E+01,-0.180517E+01,
     2 -0.161495E+01,-0.146767E+01,-0.135004E+01,-0.125346E+01,
     3 -0.117244E+01,-0.110329E+01,-0.104339E+01,-0.990760E+00,
     4 -0.943970E+00,-0.901990E+00,-0.864000E+00,-0.829390E+00,
     5 -0.797650E+00,-0.768380E+00,-0.741270E+00,-0.716030E+00,
     6 -0.692440E+00,-0.670320E+00,-0.649500E+00,-0.629860E+00,
     7 -0.611260E+00,-0.593620E+00,-0.576840E+00,-0.560840E+00,
     8 -0.545570E+00,-0.530950E+00,-0.516930E+00,-0.503480E+00,
     9 -0.490540E+00,-0.478070E+00,-0.466050E+00,-0.454450E+00,
     a -0.443220E+00,-0.432360E+00,-0.421820E+00,-0.411610E+00,
     b -0.401690E+00,-0.392040E+00,-0.382660E+00,-0.373520E+00,
     c -0.364620E+00,-0.355930E+00,-0.347450E+00,-0.339170E+00,
     d -0.331070E+00,-0.323150E+00,-0.315390E+00,-0.307790E+00,
     e -0.300340E+00,-0.293030E+00,-0.285860E+00,-0.278820E+00,
     f -0.271900E+00,-0.265100E+00,-0.258410E+00,-0.251830E+00,
     g -0.245340E+00,-0.238950E+00,-0.232650E+00,-0.226440E+00,
     h -0.220310E+00,-0.214260E+00,-0.208290E+00,-0.202380E+00          
     * /   
          data ( MuNTbl(i,   10), i=  69,   101)/
     1 -0.196540E+00,-0.190760E+00,-0.185040E+00,-0.179370E+00,
     2 -0.173760E+00,-0.168190E+00,-0.162670E+00,-0.157190E+00,
     3 -0.151750E+00,-0.146340E+00,-0.140950E+00,-0.135600E+00,
     4 -0.130260E+00,-0.124940E+00,-0.119630E+00,-0.114320E+00,
     5 -0.109010E+00,-0.103700E+00,-0.983700E-01,-0.930200E-01,
     6 -0.876300E-01,-0.822000E-01,-0.767200E-01,-0.711500E-01,
     7 -0.654800E-01,-0.596900E-01,-0.537300E-01,-0.475400E-01,
     8 -0.410500E-01,-0.340600E-01,-0.258200E-01,-0.178000E-01,
     9 -0.200000E-04                                                    
     * /   
          data ( MuNTbl(i,   11), i=  1,   68)/
     1 -0.300000E+01,-0.240477E+01,-0.203076E+01,-0.177108E+01,
     2 -0.158100E+01,-0.143556E+01,-0.132020E+01,-0.122609E+01,
     3 -0.114729E+01,-0.107998E+01,-0.102156E+01,-0.970180E+00,
     4 -0.924480E+00,-0.883450E+00,-0.846310E+00,-0.812460E+00,
     5 -0.781400E+00,-0.752760E+00,-0.726210E+00,-0.701490E+00,
     6 -0.678380E+00,-0.656710E+00,-0.636310E+00,-0.617050E+00,
     7 -0.598820E+00,-0.581530E+00,-0.565070E+00,-0.549380E+00,
     8 -0.534400E+00,-0.520060E+00,-0.506310E+00,-0.493110E+00,
     9 -0.480420E+00,-0.468190E+00,-0.456390E+00,-0.445000E+00,
     a -0.433990E+00,-0.423330E+00,-0.412990E+00,-0.402970E+00,
     b -0.393230E+00,-0.383770E+00,-0.374560E+00,-0.365600E+00,
     c -0.356860E+00,-0.348340E+00,-0.340020E+00,-0.331890E+00,
     d -0.323950E+00,-0.316170E+00,-0.308570E+00,-0.301110E+00,
     e -0.293810E+00,-0.286640E+00,-0.279610E+00,-0.272710E+00,
     f -0.265930E+00,-0.259260E+00,-0.252700E+00,-0.246250E+00,
     g -0.239890E+00,-0.233630E+00,-0.227460E+00,-0.221380E+00,
     h -0.215380E+00,-0.209450E+00,-0.203600E+00,-0.197810E+00          
     * /   
          data ( MuNTbl(i,   11), i=  69,   101)/
     1 -0.192090E+00,-0.186440E+00,-0.180840E+00,-0.175300E+00,
     2 -0.169800E+00,-0.164360E+00,-0.158960E+00,-0.153600E+00,
     3 -0.148270E+00,-0.142980E+00,-0.137720E+00,-0.132480E+00,
     4 -0.127260E+00,-0.122060E+00,-0.116870E+00,-0.111690E+00,
     5 -0.106500E+00,-0.101310E+00,-0.961100E-01,-0.908800E-01,
     6 -0.856200E-01,-0.803200E-01,-0.749600E-01,-0.695200E-01,
     7 -0.639900E-01,-0.583300E-01,-0.525100E-01,-0.464600E-01,
     8 -0.401200E-01,-0.333000E-01,-0.252000E-01,-0.174700E-01,
     9 -0.100000E-04                                                    
     * /   
          data ( MuNTbl(i,   12), i=  1,   68)/
     1 -0.300000E+01,-0.238211E+01,-0.199913E+01,-0.173882E+01,
     2 -0.155119E+01,-0.140901E+01,-0.129651E+01,-0.120452E+01,
     3 -0.112740E+01,-0.106148E+01,-0.100422E+01,-0.953830E+00,
     4 -0.909000E+00,-0.868730E+00,-0.832260E+00,-0.799010E+00,
     5 -0.768490E+00,-0.740340E+00,-0.714240E+00,-0.689940E+00,
     6 -0.667210E+00,-0.645890E+00,-0.625820E+00,-0.606880E+00,
     7 -0.588940E+00,-0.571910E+00,-0.555720E+00,-0.540280E+00,
     8 -0.525520E+00,-0.511410E+00,-0.497870E+00,-0.484870E+00,
     9 -0.472370E+00,-0.460330E+00,-0.448710E+00,-0.437500E+00,
     a -0.426650E+00,-0.416150E+00,-0.405970E+00,-0.396100E+00,
     b -0.386510E+00,-0.377190E+00,-0.368130E+00,-0.359300E+00,
     c -0.350690E+00,-0.342300E+00,-0.334110E+00,-0.326110E+00,
     d -0.318290E+00,-0.310630E+00,-0.303140E+00,-0.295810E+00,
     e -0.288620E+00,-0.281570E+00,-0.274650E+00,-0.267850E+00,
     f -0.261180E+00,-0.254620E+00,-0.248170E+00,-0.241820E+00,
     g -0.235570E+00,-0.229410E+00,-0.223340E+00,-0.217360E+00,
     h -0.211450E+00,-0.205630E+00,-0.199870E+00,-0.194190E+00          
     * /   
          data ( MuNTbl(i,   12), i=  69,   101)/
     1 -0.188570E+00,-0.183010E+00,-0.177510E+00,-0.172060E+00,
     2 -0.166660E+00,-0.161310E+00,-0.156010E+00,-0.150740E+00,
     3 -0.145510E+00,-0.140320E+00,-0.135150E+00,-0.130010E+00,
     4 -0.124880E+00,-0.119780E+00,-0.114680E+00,-0.109600E+00,
     5 -0.104510E+00,-0.994100E-01,-0.943100E-01,-0.891800E-01,
     6 -0.840200E-01,-0.788100E-01,-0.735600E-01,-0.682200E-01,
     7 -0.628000E-01,-0.572400E-01,-0.515300E-01,-0.456000E-01,
     8 -0.393800E-01,-0.327000E-01,-0.247300E-01,-0.172000E-01,
     9 -0.100000E-04                                                    
     * /   
          data ( MuNTbl(i,   13), i=  1,   68)/
     1 -0.300000E+01,-0.235981E+01,-0.197039E+01,-0.171262E+01,
     2 -0.152890E+01,-0.138951E+01,-0.127903E+01,-0.118860E+01,
     3 -0.111274E+01,-0.104785E+01,-0.991460E+00,-0.941810E+00,
     4 -0.897620E+00,-0.857910E+00,-0.821940E+00,-0.789130E+00,
     5 -0.759020E+00,-0.731230E+00,-0.705470E+00,-0.681470E+00,
     6 -0.659030E+00,-0.637970E+00,-0.618140E+00,-0.599420E+00,
     7 -0.581700E+00,-0.564870E+00,-0.548870E+00,-0.533610E+00,
     8 -0.519030E+00,-0.505070E+00,-0.491690E+00,-0.478840E+00,
     9 -0.466480E+00,-0.454580E+00,-0.443090E+00,-0.432000E+00,
     a -0.421280E+00,-0.410900E+00,-0.400840E+00,-0.391080E+00,
     b -0.381600E+00,-0.372380E+00,-0.363420E+00,-0.354690E+00,
     c -0.346180E+00,-0.337890E+00,-0.329790E+00,-0.321880E+00,
     d -0.314150E+00,-0.306580E+00,-0.299180E+00,-0.291930E+00,
     e -0.284820E+00,-0.277850E+00,-0.271020E+00,-0.264300E+00,
     f -0.257710E+00,-0.251230E+00,-0.244850E+00,-0.238580E+00,
     g -0.232400E+00,-0.226320E+00,-0.220330E+00,-0.214420E+00,
     h -0.208590E+00,-0.202830E+00,-0.197150E+00,-0.191540E+00          
     * /   
          data ( MuNTbl(i,   13), i=  69,   101)/
     1 -0.185990E+00,-0.180500E+00,-0.175070E+00,-0.169690E+00,
     2 -0.164360E+00,-0.159080E+00,-0.153850E+00,-0.148650E+00,
     3 -0.143490E+00,-0.138370E+00,-0.133270E+00,-0.128190E+00,
     4 -0.123140E+00,-0.118110E+00,-0.113080E+00,-0.108060E+00,
     5 -0.103050E+00,-0.980200E-01,-0.929900E-01,-0.879300E-01,
     6 -0.828400E-01,-0.777100E-01,-0.725300E-01,-0.672700E-01,
     7 -0.619200E-01,-0.564400E-01,-0.508100E-01,-0.449700E-01,
     8 -0.388200E-01,-0.322500E-01,-0.243800E-01,-0.169900E-01,
     9  0.000000E+00                                                    
     * /   
          data ( MuNTbl(i,   14), i=  1,   68)/
     1 -0.300000E+01,-0.234033E+01,-0.195058E+01,-0.169674E+01,
     2 -0.151557E+01,-0.137792E+01,-0.126871E+01,-0.117927E+01,
     3 -0.110419E+01,-0.103994E+01,-0.984080E+00,-0.934890E+00,
     4 -0.891080E+00,-0.851710E+00,-0.816040E+00,-0.783500E+00,
     5 -0.753630E+00,-0.726060E+00,-0.700490E+00,-0.676670E+00,
     6 -0.654390E+00,-0.633490E+00,-0.613810E+00,-0.595220E+00,
     7 -0.577620E+00,-0.560910E+00,-0.545020E+00,-0.529860E+00,
     8 -0.515380E+00,-0.501520E+00,-0.488230E+00,-0.475460E+00,
     9 -0.463180E+00,-0.451360E+00,-0.439950E+00,-0.428930E+00,
     a -0.418280E+00,-0.407960E+00,-0.397970E+00,-0.388270E+00,
     b -0.378850E+00,-0.369700E+00,-0.360790E+00,-0.352120E+00,
     c -0.343670E+00,-0.335430E+00,-0.327380E+00,-0.319520E+00,
     d -0.311840E+00,-0.304320E+00,-0.296970E+00,-0.289770E+00,
     e -0.282710E+00,-0.275780E+00,-0.268990E+00,-0.262320E+00,
     f -0.255770E+00,-0.249330E+00,-0.243000E+00,-0.236770E+00,
     g -0.230640E+00,-0.224600E+00,-0.218640E+00,-0.212770E+00,
     h -0.206980E+00,-0.201270E+00,-0.195630E+00,-0.190050E+00          
     * /   
          data ( MuNTbl(i,   14), i=  69,   101)/
     1 -0.184540E+00,-0.179090E+00,-0.173700E+00,-0.168360E+00,
     2 -0.163080E+00,-0.157830E+00,-0.152640E+00,-0.147480E+00,
     3 -0.142360E+00,-0.137270E+00,-0.132210E+00,-0.127170E+00,
     4 -0.122160E+00,-0.117160E+00,-0.112180E+00,-0.107200E+00,
     5 -0.102220E+00,-0.972300E-01,-0.922400E-01,-0.872200E-01,
     6 -0.821700E-01,-0.770800E-01,-0.719400E-01,-0.667200E-01,
     7 -0.614100E-01,-0.559800E-01,-0.504000E-01,-0.446000E-01,
     8 -0.385000E-01,-0.319800E-01,-0.241800E-01,-0.168600E-01,
     9  0.000000E+00                                                    
     * /   
          data ( MuNTbl(i,   15), i=  1,   68)/
     1 -0.300000E+01,-0.232414E+01,-0.193753E+01,-0.168662E+01,
     2 -0.150721E+01,-0.137075E+01,-0.126241E+01,-0.117362E+01,
     3 -0.109906E+01,-0.103522E+01,-0.979710E+00,-0.930800E+00,
     4 -0.887250E+00,-0.848090E+00,-0.812610E+00,-0.780240E+00,
     5 -0.750520E+00,-0.723080E+00,-0.697630E+00,-0.673920E+00,
     6 -0.651740E+00,-0.630930E+00,-0.611340E+00,-0.592830E+00,
     7 -0.575300E+00,-0.558670E+00,-0.542840E+00,-0.527740E+00,
     8 -0.513320E+00,-0.499510E+00,-0.486270E+00,-0.473560E+00,
     9 -0.461330E+00,-0.449550E+00,-0.438180E+00,-0.427200E+00,
     a -0.416590E+00,-0.406310E+00,-0.396360E+00,-0.386690E+00,
     b -0.377310E+00,-0.368190E+00,-0.359320E+00,-0.350680E+00,
     c -0.342260E+00,-0.334050E+00,-0.326030E+00,-0.318200E+00,
     d -0.310550E+00,-0.303060E+00,-0.295730E+00,-0.288550E+00,
     e -0.281520E+00,-0.274620E+00,-0.267860E+00,-0.261210E+00,
     f -0.254680E+00,-0.248270E+00,-0.241960E+00,-0.235750E+00,
     g -0.229640E+00,-0.223630E+00,-0.217700E+00,-0.211850E+00,
     h -0.206080E+00,-0.200390E+00,-0.194770E+00,-0.189220E+00          
     * /   
          data ( MuNTbl(i,   15), i=  69,   101)/
     1 -0.183730E+00,-0.178300E+00,-0.172930E+00,-0.167610E+00,
     2 -0.162350E+00,-0.157130E+00,-0.151950E+00,-0.146810E+00,
     3 -0.141710E+00,-0.136640E+00,-0.131610E+00,-0.126590E+00,
     4 -0.121600E+00,-0.116620E+00,-0.111660E+00,-0.106700E+00,
     5 -0.101740E+00,-0.967800E-01,-0.918000E-01,-0.868100E-01,
     6 -0.817800E-01,-0.767100E-01,-0.715900E-01,-0.664000E-01,
     7 -0.611100E-01,-0.557100E-01,-0.501500E-01,-0.443700E-01,
     8 -0.383000E-01,-0.318200E-01,-0.240600E-01,-0.167800E-01,
     9  0.000000E+00                                                    
     * /   
          data ( MuNTbl(i,   16), i=  1,   68)/
     1 -0.300000E+01,-0.231337E+01,-0.193063E+01,-0.168158E+01,
     2 -0.150326E+01,-0.136752E+01,-0.125968E+01,-0.117126E+01,
     3 -0.109697E+01,-0.103336E+01,-0.978030E+00,-0.929270E+00,
     4 -0.885830E+00,-0.846780E+00,-0.811390E+00,-0.779100E+00,
     5 -0.749440E+00,-0.722060E+00,-0.696670E+00,-0.673010E+00,
     6 -0.650870E+00,-0.630100E+00,-0.610540E+00,-0.592070E+00,
     7 -0.574570E+00,-0.557960E+00,-0.542150E+00,-0.527080E+00,
     8 -0.512680E+00,-0.498900E+00,-0.485680E+00,-0.472980E+00,
     9 -0.460770E+00,-0.449000E+00,-0.437650E+00,-0.426690E+00,
     a -0.416090E+00,-0.405820E+00,-0.395880E+00,-0.386230E+00,
     b -0.376860E+00,-0.367750E+00,-0.358880E+00,-0.350250E+00,
     c -0.341840E+00,-0.333640E+00,-0.325630E+00,-0.317810E+00,
     d -0.310170E+00,-0.302690E+00,-0.295370E+00,-0.288200E+00,
     e -0.281170E+00,-0.274280E+00,-0.267520E+00,-0.260880E+00,
     f -0.254360E+00,-0.247950E+00,-0.241650E+00,-0.235450E+00,
     g -0.229350E+00,-0.223340E+00,-0.217410E+00,-0.211570E+00,
     h -0.205810E+00,-0.200120E+00,-0.194510E+00,-0.188960E+00          
     * /   
          data ( MuNTbl(i,   16), i=  69,   101)/
     1 -0.183480E+00,-0.178060E+00,-0.172690E+00,-0.167380E+00,
     2 -0.162120E+00,-0.156900E+00,-0.151730E+00,-0.146600E+00,
     3 -0.141510E+00,-0.136440E+00,-0.131410E+00,-0.126400E+00,
     4 -0.121410E+00,-0.116440E+00,-0.111480E+00,-0.106530E+00,
     5 -0.101580E+00,-0.966200E-01,-0.916500E-01,-0.866600E-01,
     6 -0.816400E-01,-0.765800E-01,-0.714600E-01,-0.662800E-01,
     7 -0.610000E-01,-0.556000E-01,-0.500500E-01,-0.442900E-01,
     8 -0.382200E-01,-0.317500E-01,-0.240000E-01,-0.167400E-01,
     9  0.000000E+00                                                    
     * /   
          data ( MuNTbl(i,   17), i=  1,   68)/
     1 -0.300000E+01,-0.231061E+01,-0.192993E+01,-0.168175E+01,
     2 -0.150387E+01,-0.136836E+01,-0.126066E+01,-0.117231E+01,
     3 -0.109807E+01,-0.103448E+01,-0.979160E+00,-0.930400E+00,
     4 -0.886960E+00,-0.847900E+00,-0.812490E+00,-0.780180E+00,
     5 -0.750510E+00,-0.723120E+00,-0.697700E+00,-0.674030E+00,
     6 -0.651880E+00,-0.631090E+00,-0.611510E+00,-0.593020E+00,
     7 -0.575500E+00,-0.558880E+00,-0.543050E+00,-0.527970E+00,
     8 -0.513550E+00,-0.499750E+00,-0.486510E+00,-0.473800E+00,
     9 -0.461570E+00,-0.449790E+00,-0.438430E+00,-0.427450E+00,
     a -0.416830E+00,-0.406560E+00,-0.396600E+00,-0.386930E+00,
     b -0.377550E+00,-0.368420E+00,-0.359550E+00,-0.350900E+00,
     c -0.342480E+00,-0.334260E+00,-0.326240E+00,-0.318410E+00,
     d -0.310750E+00,-0.303260E+00,-0.295930E+00,-0.288750E+00,
     e -0.281710E+00,-0.274810E+00,-0.268030E+00,-0.261380E+00,
     f -0.254850E+00,-0.248430E+00,-0.242120E+00,-0.235910E+00,
     g -0.229790E+00,-0.223770E+00,-0.217830E+00,-0.211980E+00,
     h -0.206200E+00,-0.200510E+00,-0.194880E+00,-0.189320E+00          
     * /   
          data ( MuNTbl(i,   17), i=  69,   101)/
     1 -0.183830E+00,-0.178390E+00,-0.173020E+00,-0.167690E+00,
     2 -0.162420E+00,-0.157200E+00,-0.152010E+00,-0.146870E+00,
     3 -0.141770E+00,-0.136690E+00,-0.131650E+00,-0.126630E+00,
     4 -0.121630E+00,-0.116650E+00,-0.111680E+00,-0.106710E+00,
     5 -0.101750E+00,-0.967800E-01,-0.918000E-01,-0.868000E-01,
     6 -0.817700E-01,-0.767000E-01,-0.715700E-01,-0.663700E-01,
     7 -0.610800E-01,-0.556800E-01,-0.501100E-01,-0.443300E-01,
     8 -0.382600E-01,-0.317800E-01,-0.240200E-01,-0.167400E-01,
     9  0.000000E+00                                                    
     * /   

      end
!     *******************
      block data cblkObs
!     *******************
      implicit none

















!       Zobs.h     header file for observation sites definition
!
           integer maxNoOfSites, maxNoOfASSites, horizontal,
     *     perpendicular, notUsed, spherical
       parameter (

     *    maxNoOfSites = 20,


     *    maxNoOfASSites=20,

     *    notUsed = 0,           ! detector plane is not used
     *    horizontal = 1,        ! detector is horizontal
     *    perpendicular = 2,      ! detector is pependicular to 1ry.
     *    spherical = 3          ! detector is cocentric sphere as the earth
     *		      )

!       Zobsvp.h---parameters to be given by input.
!       This must be preceded by Zobs.h

!	(->	---------------------------------------------------

	real*8  HeightList  !1  Height of observation levels in m. This is  made from DepthList internally. 
                            ! I.e., this one is usually not an input. However, if the DepthList values are 
                            ! negative, this is used as input and corresponding DepthList is computed internally.
        real*8  DepthList   !1	Depth List of Observation level in kg/m$^2$. If $< 0$, HeightList has priority. 
                            !  (See HeightList)
        real*8  ASHeightList	!1  This is HeightList for Air Shower observ.  Used only if Generate contains
                            !  "as". See  HeightList.
        real*8  ASDepthList     !1  This is DepthList for AS observation.  Used only if Generate contains 
                            ! "as". See DepthList.
        real*8  LatitOfSite     !1  Latitude of the deepest observation level in degree.  East is positive.
        real*8  LongitOfSite    !1  Longitude of the deepest observation level in degree.  North is positive.

     	real*8  DtGMT           !1  Difference of the local time of the observation place from GMT (hour).
	real*8  YearOfGeomag    !1  Like 1999.5. Year when Geomagnetic field is to be calculated.
	integer ObsPlane        !1    How to observe particles. \newline
                                !    0$ \Rightarrow $ no detector plane is used for observation. BorderHeightL
                                !    and BorderHeightH are used to detect particles. This is for, say, neutrino
                                !    observation. See BorderHeight{L,H}. However, the primary is directed to
                                !    the deepest depth.  \newline
                                !    1,-1$ \Rightarrow $ detector at the observation place is horizontal. Note 
                                !    that the horizontal means not tangential plane, but rather a spherical surface \newline
                                !    2,-2$ \Rightarrow $ detector is perpendicular to the primary. \newline
                                !    3$ \Rightarrow $ spherical observation. See text. \newline
                                !    For ObsPlane={1,2}, the user observation routine will receive coordinate values in
                                !    the corresponding detector system. However, if it is 0, 3 or negative, Exyz values
                                !    are obtained.
        integer NoOfSites2    !2   No of Sites for particle observation; not to be touched; for skeleton/flesh use.
	real*8 XaxisFromSouth   !2 Angle between the horizontal detector X-axis and the south(deg). + is counter
                                ! clockwise.  If $|$XaxisFromSouth$| > 360$, it is computed so that the direction is
                                ! to the magnetic east at the deepest observation point. Default is 361.
!	<-)	--------------------------------------------

   

        common /Zobsc/
     *	 HeightList(maxNoOfSites),
     *   DepthList(maxNoOfSites),
     *   ASHeightList(maxNoOfASSites),
     *   ASDepthList(maxNoOfASSites),
     *   LatitOfSite, LongitOfSite, DtGMT,
     *   XaxisFromSouth, YearOfGeomag,
     *   ObsPlane, NoOfSites2

          
       data 
     * DepthList /maxNoOfSites*0./ ,
     * HeightList /maxNoOfSites*0./ ,
     * ASDepthList /maxNoOfASSites*0./ ,
     * ASHeightList /maxNoOfASSites*0./ ,
     * DtGMT /8.0/ ,
     * LatitOfSite /30.11/ ,
     * LongitOfSite /90.53/ ,
     * ObsPlane /1/ ,	
     * XaxisFromSouth /361.0/,
     * YearOfGeomag /2000.5/
      end
!     *************************
      block data  cblkTracking
!     *************************
      implicit none

!           Parameters used for Tracking.
!	(->	--------------------------------------------

	logical  ExactThick   !2  If T, a given length is converted into thickness with best accuracy even for very 
                              !  inclined trajectory by using numerical integration.
        logical IncMuonPolari  !1  if T, consider muon polarization
        logical Freec          !1  if F, the first interaction point is forced to be the injection point else
                               !   the interaction poin is randomly sampled.
	integer OneDim	       !1  If 0, 3 dimensional simulation. if $\ge$1, one
                               !  dimensional simulation is performed.  \newline
                               !  1: onedim without use of table. \newline
                               !  2: table is used for thickness $ \leftrightarrow$ length conversion. if cos $<$ .5 \newline
                               !  3: table is always used for any angle.
                               !  ( for height $>$ 30 km, table is not used in any case). 
        real*8  LamorDiv       !2  In the geomagnetic field, a charged particle can travel almost streight 
                               !  in (Lamor Radius)/LamorDiv.  Default is 5. For AMS like tracking 20 may be needed.
	real*8  Truncc         !2 coeff. for truncating path.
	real*8  Truncn         !2 coeff. for truncating path.
	real*8  Truncx         !2 coeff. for truncating path.
	real*8  KEminObs(8)      !1  The min kinetic energy of particles for observation.
	real*8  KEminObs2(8)     !2  Don't touch this. skeleton/flesh use.
	real*8  RatioToE0      !2 In the A.S generation,  hadronic interactions are followed down to at
                               !  least  RatioToE0 * E0/nucleon energy.
        real*8  WaitRatio      !1  Wait A.S generation until the electron energy, Ee, becomes $<$ WaitRatio* E0. 
                               !   This many be 1.0 for hadron origin case. But for gamma/electron primary, 
                               !   this should be as low as 0.01 to enjoy full fluctuation.
        integer EndLevel !2 Used for skeleton/flesh-out job. In a normal job, system default value 0 is reset by
            ! the system to be the max number of observation levels. (=NoOfSites).  Its real use is in such a
            ! skeleton/flesh-out job that one first follows the particles up to some high depth and later chooses
            ! events and flesh them out to deeper depths. In such a skeleton-making job, the user must give the
            ! depth list which is used flesh-out job, too.  In the skeleton job, particle tracking is terminated
            ! at the level specified by EndLevel.  In such a flesh-out job, the user must give a larger value 
            ! or 0 to EndLevel 
        integer EndLevel2 !2  Don't worry. This  is system use.
	integer Trace   !1  Flag for trace information output.\newline
            !   0 $\rightarrow$ no trace information is output.\newline
            !  $<$10$\rightarrow$ x, y, z in the primary system(say, 1)\newline
            !  $<$20 $\rightarrow$ x, y, in the primary sys. z in kg/m$^2$.(say,11)\newline
            !  $<$30 $\rightarrow$ x, y, z in the detector system\newline
            !  $<$40 $\rightarrow$ x, y, in the detector system. z in kg/m2\newline 
            !  $<$50 $\rightarrow$ x, y, z in 'xyz' system.\newline
            !  $<$60 $\rightarrow$ x,y, in 'xyz' and z in kg/m2\newline
            !  61-100 $\rightarrow$  for Cherenkov observation. For Coord system,  subtract 60.\newline
            !   if the value is even,  binary output is made on TraceDev.\newline
            !   if the last digit is 1 or 2, trace is always taken. if the last digit is 3 or 4, trace is taken
            !   only if the particle is located below the heighest observation depeth. 
            !  $>$ 101  $\rightarrow$ subtract 100 and apply the above, but chookTrace or chookCeren are used.\newline
            !  Primary system:  Origin is the deepest detector.  Z-axis is the primary direction. 
            !                   X-axis is Z x Vertical axis.  X-Y plane is orthogonal to the primary.\newline
            !  Detector system: origin is the deepest detector. Z-axis is the vertical one.  X-axis is 
            !                   directed to the magnetic east.  X-Y plane is horizontal.\newline
            !  z in kg/m$^2$ :     Vertical depth in kg/m$^2$  above the  deepest detector to the particle.
	integer TraceDev      !2  Logical dev \# for TraceDir/trace1,2,.... 
        character*70  TraceDir !1 Directory.  Default Trace information is put TraceDir/trace1, 2,..
                               ! for event 1, 2, ... The directory should exist.  Default is ' ' and in this case
                               ! /tmp/YourLoginName/ is employed. The last "/" should not be given.
                               ! *** NOTE that default Cherenkov output is made only using TraceDev,
                               !    TraceDir is not used.  You have to open the disk file at chookbgRun
                               !    It can by binary or ascii file depending on Trace value.
	logical ThinSampling  !1  if F, thinsampling is not tried. if T, alla Hillas thinning. Don't use with
                              !    the  skeleton/flesh method 
        real*8  EthinRatio    !2  if ThinsSamplig != F, thin sampling is performed if the energy of a particle is
                              !   $<$ EthinRatio * PrimaryEnergy(/nucleon) (=Ethin) ( EtinRatio$>$ 0).
                              !  If EthinRatio $<$ 0, Ethin will be |EthinRatio| (GeV). 
        logical TimeStructure !1  If T,  time information is computed
	integer HowGeomag     !2  if 1, no magnetic field until first coll. \newline
                              !      2, mag.f always exists. If Reverse not=0, use this. \newline
                              !     11, same as 1 but mag.f is const. \newline
                              !     12, same as 2 but mag.f is const. \newline
                              !     21, same as 1 but mag.f is const. \newline
                              !     22, same as 2 but mag.f is const.  \newline
                              !     31, same as 1 but mag.f is dependent on the position.  \newline
                              !  const value is the one at deepest observation plane. for 11,12  or should be given by
                              !  MagN, MagE, MagD for 21, 22. For normal applications,  11 is good.
                              !  If no magnetic field is applied, energy loss by dE/dx is considered.(bef.4.92,
                              !  and aft. 5.14)
        real*8  MagN          !2 See HowGeomag (in Tesla)
	real*8  MagE          !2 See HowGeomag (in Tesla)
        real*8  MagD          !2 See HowGeomag (in Tesla)

        real*8  MagChgDist    !2 Distance where mag. can be seen  as const.(m) at sea level
        integer UseRungeKutta !2 How to calculate deflection by the geomagnetic field. Let L be the distance
                              !   the  particle travels. \newline
                              ! 0$\rightarrow$Don't use RungeKutta method. Use the solution assuming the constant B, which
                              !     is exact if B is const.  Since the particle path is made  short, this is
                              !     enough for normal cases where particles are inside the atmosphere.(default) \newline
                              !     In every case below, if the particle height is $<$ 30km  
                              !    (= cheight in ccomPathEnd.f), the same method as 0 is used. \newline
                              ! 1$\rightarrow$ Use the Euler method.  Time needed is 20\% more than the 0 case.
                              !      As B, use the value at L/2 point obtained by using the current direction. \newline
                              ! 2$\rightarrow$ mixture of 1 and Runge-Kutta-Gill method. If gradient of B is large, RKG is
                              !      employed. This needs $\sim$4 times more cpu time than case of 1 when making a
                              !      cutoff table.  The step size of RKG is $\sim$1/10 of the Lamore radius. \newline
                              ! 3$\rightarrow$ The same as 2 but use the Runge-Kutta-Fehlberg method instead of RKG.
                              !      Step size is automattically adjusted ($\sim$1/20 $\sim$1/30 of Lamor radius) \newline
                              ! 4$\rightarrow$ As a middle point, use the point obtained by assuming the constant B at
                              !      initial point. If grad B is still large, use RKG. \newline
                              ! 5$\rightarrow$ The same as 4 but us RKF instead of RKG. \newline
                              ! 6$\rightarrow$ Use always RKG \newline
                              ! 7$\rightarrow$ Use always RKF.  This takes very long time.(50 times of 0). \newline
	real*8  BorderHeightH !2 If a particle goes higher than this, discard it.  This should be larger than 
                              !  HeightOfInj or 0.  
                              !  If 0, it is adjusted to be the same as HeightOfInj. NOTE: For upgoin primary cases, you have
                              !  to set this one explicitly. 
	real*8  BorderHeightL !2 If a particle reaches this hight, call observation routine. No further tracking is done. 
                              !  This is for neutrino observation.  See ObsPlane.
        real*8  BackAngLimit  !2 If the cosine of the angle between a particle and the primary becomes smaller than
                              !  this value, the particle is discarded. See also BorderHeighH. If you give a value 
                              !  less than -1.0, such rejection will never happen.   Default is -1.0
	character*16 Generate !1 specify what should be generated \newline
                              !   1)  Electro-magnetic cascade(em), \newline
                              !   2)  one dimensional  hybrid AS(as/qas) and/or \newline
                              !   3)  AS Lateral distribution(lat). \newline
                              !    If Generate= ' ', hadronic cascade shower is generated. \newline
                              !  For example, you may give as follows: \newline
                              !    Generate='em,as' or 'em/as' (order/case/separator insensitive) is to  generate EM-cascade and AS. \newline
                              !    Generate='as' will  generate AS with some  adequate EM cascade (EM cascade is automatically generated
                              !    so that hybrid A.S can be observed, but the minimum energy in EM cascade is independent of KEminObs). \newline
                              !   If 'qas' is given, quick generation of AS for heavy primaries is tried. See  chookASbyH.f 

        character*16 Generate2 !2  don't touch this.  for skeleton/flesh use.

	integer MagBrem       !2 If 0, no magnetic bremsstrahlung is considered. \newline
                              ! if 1 and Ee $>$ MagBremEmin, energy loss due to magnetic brems is considered \newline
                              ! if 2 and Ee $>$ MagBremEmin, real sampling of gamma is performed. \newline
                              ! (note, actually upsilon is referred further).
                              ! if generate='as' with really high energy primaries, WaitRatio
                              ! must be made small so that WaitRatio*E0 $\sim$ MagBremEmin 
        integer MagPair       !2 If 0, no magnetic pair creation is considered. \newline
                              !  if 1 and Eg > MagPairEmin, real sampling is tried.
                              ! (note, actually upsilon is referred further).  To see these magnetic effects,
                              !  HowGeoMag=2 and HightOfInj $\sim$ 5000 km are desirable.

        logical LpmEffect     !1 If t, the LPM effect is  considered when Ee $>$ LpmBremEmin for electrons and
                              !        Eg $>$ LpmPairEmin for gamma rays.

        real*8  MagBremEmin   !2  E $>$ this, magnetic bremsstrahlung by electrons may be considered. However, if
                              !   MagBrem = 0, not considered at all \newline
                              !   MagBrem = 1, total energy loss due to brems is considered. \newline
                              !   MagBrem = 2, gamma energy is sampled actually. \newline
                              !   If upsilon (Ee/m * B/Bcr) is small, the effective treatment will be
                              !   the same as MagBrem = 0 case.
        real*8  MagPairEmin   !2  E $>$ this, magnetic pair creation by gamma may be considered. However, if
                              !   MagPair = 0, not considered at all. \newline
                              !   MagPair = 1, pair creation is sampled.  \newline
                              !   However, again, actual occurrence will be dependent on the angle between
                              !   B and photon direction.
        real*8  UpsilonMin    !2  Magnetic bremsstralhung is considered only if upsilon $>$ UpsilonMin.
        real*8  LpmBremEmin   !2  The LPM effect is taken into account for bremsstrahlung when LpmEffect is .true. 
                              !   and the electron energy is higher than this.
        real*8  LpmPairEmin   !2  The LPM effect is taken into account for pair creation when LpmEffect is .true.
                              !    and the gamma energy is higher than this.
        integer Reverse       !2  0$\rightarrow$ Normal tracking. \newline
                              !   1$\rightarrow$ incident is tracked to a direction opposite to the given one.
                              !       the incident is charge-conjugated.
                              !       All interactions are ignored. (Use when to make cut-off table or to see
                              !       a given particle (say, observed anti proton) can go out of Earth. \newline
                              !   2$\rightarrow$ same as 1 but energy gain (not loss) is taken into account
                              !   TimeStructure should be T if Reverse != 0.  See BackAnglLimit.

        real*8 PathLimit      !2  If the sum of (path/beta) of a particle  exceeds this, it is judged as dead.
                              !   (to avoid infinite cyclotron loop).  However, for normal applications,
                              !   this will not be effective because of BackAnglLimit. See Reverse.
                              !   TimeStructure should be T if Reverse != 0 and PathLimit is to be effective.

       integer MuNI           !2 0$\rightarrow$ nuclear interaction of muon is completely neglected \newline
                              !  1$\rightarrow$ energy loss by n.i is subsumed in dE/dx of muons as a continuous energy loss.  Let v=
                              !     Etransfer/Emu,  the loss here is Int(vc:vmax) of (Emu vdsigma/dv).  (vc $\sim$0, vmax$\sim$1). \newline
                              !  2$\rightarrow$ (Default value). similar to 1 but as the continuous loss only v $<$ vmin=10$^{-3}$ of
                              !     fractional muon energy is subsumed (Int(vc: vmin) of (Emu vdsigma/dv)).  The portion
                              !     of loss by v$>$vmin is treated as a stocastic  process.  However, the product from the
                              !     n.i itself is neglected \newline
                              !  3$\rightarrow$ the same as 2, but the n.i is explicitly included to produce a number of particles.  
                              !     The n.i is treated as a photo-nucleus interaction.
      integer MuBr            !2  parameter similar to MuNI but for bremsstrahlung by muons.
       integer MuPr           !2  parameter similar to MuNI but for pair creation by muons.

!	 <-)	----------------------------------------------

	common /cZtracp/ Truncc, Truncn, Truncx, 
     *  KEminObs, KEminObs2, RatioToE0, PathLimit,
     *  WaitRatio,  EthinRatio, BackAngLimit, LamorDiv,
     *  BorderHeightH,  MagN, MagE, MagD, MagChgDist,
     *  BorderHeightL,  MuNI, MuBr, MuPr,
     *  MagBremEmin, MagPairEmin, UpsilonMin, LpmBremEmin, 
     *  LpmPairEmin, UseRungeKutta,
     *  ThinSampling, TimeStructure, HowGeomag,
     *  Trace, TraceDev,  ExactThick, OneDim, Reverse,
     *  Freec,  IncMuonPolari, MagBrem, MagPair, LpmEffect,
     *  EndLevel, EndLevel2
 
        common /cZtrackpc/ Generate, Generate2, TraceDir

      data 
     *  ExactThick /.false./ ,
     *  Freec /.true./ ,
     *  TraceDev /21/ , 
     *  TraceDir /' '/ ,
     *  BorderHeightH /0./ ,
     *  BorderHeightL /0./ ,
     *  HowGeomag /11/ ,
     *  MagN /0./ ,
     *  MagE /0./ ,
     *  MagD /0./ ,
     *  RatioToE0 /1.e-5/ ,
     *  UseRungeKutta /0/ ,
     *  MagChgDist /20.e3/ ,
     *  TimeStructure /.true./ ,
     *  Trace /0/ ,
     *  WaitRatio /1./ ,	
     *  Truncc/5.0/ ,
     *  Truncn /1.e-2/ ,
     *  Truncx /2./

      data
     * IncMuonPolari /.true./ ,
     * KEminObs /2*500d-6,6*50d-3/ ,
     * ThinSampling /.false./,
     * EthinRatio /3.e-5/ ,
     * Generate /'em'/ , 
     * Generate2 /' '/ ,
     * BackAngLimit  /-1.0/ ,
     * OneDim /0/,
     * MagBrem /2/,
     * MagPair /1/,
     * LpmEffect /.true./,
     * MagPairEmin /2.e10/,         ! Eg > this ==> magnetic pair considered.> 3x10^19 eV
     * MagBremEmin /3.e8/,          ! Ee > this ==> magnetic brems considered > 10^18 eV
     * UpsilonMin  /3.e-3/,         ! if upsilon is < this, no magnetic brem
     * LpmPairEmin /1.e9/,         ! Eg > LpmPairEmin && LpmEffect ==> LPM effect for g
     * LpmBremEmin /1.e8/           ! Ee > LpmBremEmin && LpmEffect ==> LPM effect for e

       data
     * Reverse /0/,
     * PathLimit /13000d4/,         !  ~ 20 x Eradius
     * EndLevel/0/,
     * MuNI /2/, 
     * MuBr /2/,
     * MuPr /2/,
     * LamorDiv /5./
      end
!           block data for Zxsectionp
      block data cblkXsec
!----      include 'Zxsectionp.h'
!           common parameters for tracking
!
!	(->  -------------------------------------

        real*8  Deltpp  !2  p-p xsection increases as $E^{Deltpp}$(E$>$ 100GeV)
        real*8  Deltpip !2  pi-p xsection increases as $E^{Deltpip}$ (E$>$ 100GeV)
        real*8  Deltkp  !2  k-p xsection increases as $E^{Deltkp}$ (E$>$ 100GeV)
        real*8  IncreaseXsec !2   how the xsection increases. 1.0$\rightarrow$  power of E

!	<-)  -----------------------------------

        common /Zxsectionp/
     *  Deltpp, Deltpip, Deltkp, IncreaseXsec
!

!
      data 
     *  Deltpp /0.08/ ,
     *  Deltpip /0.08/ ,
     *  Deltkp /0.08/ , 
     *  IncreaseXsec /1.0/
      end
      block data cblkdedx

      common /ZdedxAir/ stha, sthb, sthc, sthx0, sthx1, sthsa,
     *   wlg0, w0, betasq, tmax, norm, w0inMeV,
     *   Knckon, jdef

	real*8  stha, sthb, sthc, sthx0, sthx1, sthsa
	real*8  wlg0, w0, betasq, tmax, norm, w0inMeV
	integer jdef
        logical Knckon

!

      data jdef/0/
!
!        table of sternheimer's consts. by p.r.b vol.3 (1971)3681
!     for air
!                        
      data stha/7.68e-02/
      data sthb/ 17.86/
      data sthc/ -10.79/
      data sthx0/1.809/
      data sthx1/4.0/
      data sthsa/ 0.234/
      end

      block data  cblkSpecial

      integer modifyx
      real*8  modifyxpw1, modifyxpw2
      common /Zspecial/  modifyxpw1, modifyxpw2, modifyx


       data  
     * modifyx /0/ ,
     * modifyxpw1 /0.5/ ,
     * modifyxpw2 /0.5/
      end

      subroutine cmymain
!         this list block data names as external names so that the
!       block data is surely enabled.





!   dpmjet cannot be used on NEXTSTEP, so
!   you have to make the next 0. 



!            make DEBUG > 0 depending on the debug purpose. 


!
!   choose:    Old atmosphere or new segmented atmosphere
!            define 
!               old atmosphere --> 0
!           or  new with c-spline
!               new atmosphere --> 1
!           or  new with linear interp.
!               new atmosphere --> 2


!     if you want to put a lable on each particle to identify that
!     the one and the same particle crosses a given observation
!     plane more than once, make this 1 or 2.  Then the same particle
!     will have the  same label number in track record.
!     ( aTrack.label ).  If this is 0, aTrack.lable record dose not
!     exists. 
!     If 1; after any interaction (except for continuous energy
!     loss by dE/dx and deflection by B or scattering), label is
!     changed.
!     If 2: For knockon and Bremstrahlung, the survival particle
!     will have the same label. In the case of Moller scattring
!     higher enregy electrons are regarded as the survival one.
!

!     if you want to have a detailed info. for particle tracking
!     make the below >=1.  The user observation routine is called
!     with the following id  on the following  conditions:
!              chookobs(a, id)
!     1)  if it is >=1,  a particle is going to interact at a point given in
!         the track information, id=4
!     2)  if it is >=1,  a particle is going to die, id=5
!     3)  if it is >=2,  a particle is being discarded due to the large
!          angle (cos(angle relative to the parent) > BackAngLimit). id=6
!     4)  if it is >=3,  a particle makes a step. id=7
!        







      external cblkElemag
      external cblkdedx
      external cblkHeavy 
      external cblkManager 

      external cblkXsec
      external cblkEvhnp
      external cblkIncident
      external cblkObs
      external cblkTracking
      external cblkMuInt
!              For Lund related block common
      external blkhd1, blkhd2, blkhd3, blkhd4
      external blkdc1, blkdc2
      external ludataC
      external luedatC
      external luhdatC
      external cblkSpeical

         call cmanager   ! call Manager for Cosmos Simulation
      end
