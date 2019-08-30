!         parameters for Elemag process.
!	(->   ----------------------------------------------

!-->mod         real*8  RecoilKineMinE  !2  Recoil Kinetic Min Energy above which the recoil (=knock-on process)
                                ! is treated. Below this energy, the effect is included as continuous
                 ! energy loss.  Used only if KnockOnRatio $>$ 1.
               	 ! If this is 0 or if KnockOnRatio =1,  KEminObs(gamma)=KEminObs(elec) is used. 
                                ! See also KnockOnRatio.
!->EM        real*8 KnockOnRatio   !2  KnockOnRatio* KEminoObs is used instead of RecoilKineMinE if KnockOnRatio $<$1.
!-->    real*8  X0            !2  Radiation length in kg/m$^2$ for air. Normally the user should not touch this.
!-->    real*8  Ecrit         !2  Critical energy in GeV. \newline
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

!         logical Knockon       !2  Obsolete. Don't use this. See RecoilKineMinE
                              !    and KnockonRatio.
!->-EMcon         real*8  AnihiE        !2  If E(positron) $<$ AnihiE, annihilation is  considered. 
!-->EMcon      real*8  Es            !2  Modified scattering constant. 19.3d-3 GeV
!-->EMcon     real*8  MaxComptonE   !2  Above this energy, Compton scattering is neglected. 
!-->EMcon         real*8  MaxPhotoE     !2  Above this energy, photoelectric effect is neglected.
!         real*8  MinPhotoProdE !1  Below this energy, no photo-prod of hadron.  See also PhotoProd.
!         logical PhotoProd     !1  Switch. if .false., no photo prod. of hadron is considered at all. 
                              !   See also MinPhotoProdE, HowPhotoP
!->EMcon        real*8  Excom1        !2  (GeV).  If photon energy is <= Excom1, use XCOM data for 
                              !    compton/p.e/coherent scattering (must be < 100 GeV).
!->EMcon        real*8  Excom2        !2  (GeV).  If photon energy is <=Excom2,  use XCOM data for
                              !           pair creation cross-section.  (must be< 100 GeV). 
!-->EMcon        integer Moliere       !2  2$\rightarrow$ use Moliere scat.\newline
                              !   0$\rightarrow$ use Gaussian scattrign. \newline 
                              !   1$\rightarrow$ use Moli\`ere scattering for non-electrons \newline
                              !   2$\rightarrow$ use Moli\`ere scattering for all charged 
                              !      particles.  But treatment is not so rigorous as case of 3.
                              !      \newline
                              !   3$\rightarrow$ use rigorus Moliere scattering.  Diff. from 2 is verysmall. May be some effect in the 
                              !   core region.  
!-->EMcon        integer ALateCor      !2  1$\rightarrow$ angular and lateral correlation is taken into account when Moliere=0 .\newline
                              !   t$\rightarrow$ Use angular-lateral correlation by Gaussian
                              !   approximation. No effect is seen if path length is short.

!	<-)	----------------------------------------------

!         common /Zelemagc/    KnockOnRatio,
! RecoilKineMinE,
!     *  AnihiE, MaxComptonE,
!     *  MaxPhotoE, MinPhotoProdE,   Es, 
!     *  MinPhotoProdE,   
!     *  PhotoProd
!      Moliere,  ALateCor 
!      Excom1, Excom2,
!    X0,
!     *  Ecrit, 
! Knockon, 




