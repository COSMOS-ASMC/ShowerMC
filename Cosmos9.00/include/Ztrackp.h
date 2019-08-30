/*
c           Parameters used for Tracking.
c	(->	--------------------------------------------

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
	real*8  KEminObs       !1  The min kinetic energy of particles for observation.
	real*8  KEminObs2      !2  Don't touch this. skeleton/flesh use.
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
        real*8  EthinRatio(2) !2  if ThinsSamplig != F, thin sampling is performed if the energy of a particle is
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

c	 <-)	----------------------------------------------
*/

//              not cztrackp;  due to spelling mistake in Fortran
extern  struct cztracp {
  double truncc;
  double truncn;
  double truncx;
  double keminobs;
  double keminobs2;
  double ratiotoe0;
  double pathlimit;
  double waitratio;
  double ethinratio[2];
  double backanglimit;
  double lamordiv;
  double borderheighth;
  double magn;
  double mage;
  double magd;
  double magchgdist;
  double borderheightl;
  int muni;
  int mubr;
  int mupr;
  double magbrememin;
  double magpairemin;
  double upsilonmin;
  double lpmbrememin;
  double lpmpairemin;
  int    userungekutta;
  logical    thinsampling;
  logical    timestructure;
  int    howgeomag;
  int    trace;
  int    tracedev;
  logical    exactthick;
  int    onedim;
  int    reverse;
  logical    freec;
  logical    incmuonpolari;
  int    magbrem;
  int    magpair;
  logical    lpmeffect;
  int    endlevel;
  int    endlevel2;
} cztracp_;


#define Truncc cztracp_.truncc
#define Truncn cztracp_.truncn
#define Truncx cztracp_.truncx
#define KEminObs cztracp_.keminobs
#define KEminObs2 cztracp_.keminobs2
#define RatioToE0 cztracp_.ratiotoe0
#define PathLimit cztracp_.pathlimit
#define WaitRatio cztracp_.waitratio
#define EthinRatio cztracp_.ethinratio
#define BackAngLimit cztracp_.backanglimit
#define LamorDiv cztracp_.lamordiv]
#define BorderHeightH cztracp_.borderheighth
#define MagN  cztracp_.magn
#define MagE  cztracp_.mage
#define MagD  cztracp_.magd]
#define MagChgDist cztracp_.magchgdist
#define BorderHeightL cztracp_.borderheightl
#define MuNI cztracp_.muni
#define MuBr cztracp_.mubr
#define MuPr cztracp_.mupr
#define MagBremEmin cztracp_.magbrememin
#define MagPairEmin cztracp_.magpairemin
#define UpsilonMin cztracp_.upsilonmin
#define LpmBremEmin cztracp_.lpmbrememin
#define LpmPairEmin cztracp_.lpmpairemin
#define UseRungeKutta cztracp_.userungekutta
#define ThinSampling cztracp_.thinsampling
#define TimeStructure cztracp_.timestructure
#define HowGeomag cztracp_.howgeomag

#define TraceDir cztrackpc_.tracedir

#define TraceDev cztracp_.tracedev
#define Trace cztracp_.trace
#define ExactThick cztracp_.exactthick
#define OneDim cztracp_.onedim
#define Reverse cztracp_.reverse
#define Freec cztracp_.freec
#define IncMuonPolari cztracp_.incmuonpolari
#define MagBrem cztracp_.magbrem
#define MagPair cztracp_.magpair
#define LpmEffect cztracp_.lpmeffect
#define EndLevel2 cztracp_.endlevel2
#define EndLevel cztracp_.endlevel

extern struct cztrackpc {
  char generate[16];
  char generate2[16];
  char tracedir[70];
  DUMMYCHAR
} cztrackpc_;

#define Generate2 cztrackpc_.generate2
#define Generate cztrackpc_.generate
