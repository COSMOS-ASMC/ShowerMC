/*
  c         parameters for Elemag process.
  c	(->   ----------------------------------------------
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
  real*8  Excom1        !2  (GeV).  If photon energy is <= Excom1, use XCOM data for
  !    compton/p.e/coherent scattering (must be < 100 GeV).
  real*8  Excom2        !2  (GeV).  If photon energy is <=Excom2,  use XCOM data for
  !           pair creation cross-section.  (must be< 100 GeV).    
  integer Moliere       !2  0$\rightarrow$ use Gaussian approx always (with air density change and
  !    energy loss effect)\newline
  !   1$\rightarrow$ use Moli\`ere scattering for non-electrons (default)\newline
  !   2$\rightarrow$ use Moli\`ere scattering for all charged particles.\newline
  !    If negative, anglular-correlated displacement is made to be 0 since Moli\`ere
  !    theory cannot give it. (if $>0$, we use Gaussian approximation for correlation). 
c	<-)	----------------------------------------------
*/

extern  struct zelemagc {
  double recoilkinemine ; 
  double knockonratio  ;
  double anihie;
  double maxcomptone;
  double maxphotoe;
  double minphotoprode;
  double es;
  double x0;
  double ecrit; 
  logical knockon; 
  logical photoprod;
  double Excom1;
  double Excom2;
  int moliere;
} zelemagc_;

#define RecoilKineMinE  zelemagc_.recoilkinemine
#define KnockonRatio    zelemagc_.knockonratio  
#define AnihiE          zelemagc_.anihie
#define MaxComtonE      zelemagc_.maxcomptone
#define MaxPhotoE       zelemagc_.maxphotoe
#define MinPhotoProdE   zelemagc_.minphotoprode
#define Es              zelemagc_.es
#define X0              zelemagc_.x0
#define Ecrit           zelemagc_.ecrit
#define KnockOn         zelemagc_.knockon
#define PhotoProd       zelemagc_.photoprod
#define Excom1          zelemagc_.Excom1
#define Excom2          zelemagc_.Excom2
#define Moliere         zelemagc_.moliere
