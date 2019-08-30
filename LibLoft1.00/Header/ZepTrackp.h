!->EMcon          real*8  Bxu,  Byu,  Bzu  ! uniform mag field comp.
!->EMcon          real*8  Exu,  Eyu,  Ezu  ! uniform elec field comp.
!-->EMcon          integer MagField      ! how to give magnetic field
!-->EMcon          integer ElecField     ! //          electric field
!->EMcon        real*8  DtMax         ! max allowable path in r.l
!->EMcon         real*8  Tcoef         !
!->EMcon         real*8  Tmin          !
!->EMcon        real*8  Excom1        !  <100GeV;   use XCOM data for compton/p.e effect
!->EMcon       real*8  Excom2        !  < 100GeV;  use XCOM data for pair xsec.
!->EMcon       logical EdepdEdx      ! use energy dependent dE/dx
                                !  F is used to compare analytic calc
!-->EMcon        integer:: HowQuench   ! R: dE/dx restricted
	  !     F: dE/dx Full.  D delta ray contribution
          !                     O other // (knot known).
          !  cf: exprimetally obained quenching coef.
          ! ecf: quenching factor for R. D+O has no quenching
	  !  ecf*R + D + O = cf*F; D=F-R
          !  ecf + F/R -1 + O/R = cf*F/R; then
          !  ecf = 1 + (cf-1)F/R - O/R ; set O= 0 for first approx.
          !   if  HowQuench = 0  use cf for ecf (equiv. to  v <=9.154
          !   if            = 1  default:  use ecf as if old cf
!	  integer TargetElecBrems ! Flag to include target electron
   !    brems 0->not bit 1->mu  bit2 p,pi,K (bit3 heavy?)
! -->EMcon          logical LPMeffect     !  T(default). LPM is effective else off 
          logical TimeStruc     !  
!!!          logical Moliere       ! Moliere's scattering  or  Gaussian
!-->EMcon          integer Moliere       !  0--> Guass, 1-->Old Moliere 2-->
                                !   rigorous Moliere
!->EMcon       real*8  EminElec      ! e-'s E < EminElec is discarded
!->EMcon          real*8  EminGamma     ! Egamma < EminGamma is  discarded
	  real*8  EminH         ! Emin of Hadron/mu K.E/n 
                                ! (except for anti p,n )
                                ! if =0, max( 5 MeV, KEmin)
!->EMcon          real*8 Flpm         !  Ee > Flpm*BremEminLPM is the
                              !  actual LPM region
!->EMco         integer HowNormBrems ! -1 (default) : normalization base 
                 !is C.S
  		    ! 0: no normalization
    	            ! 1: base is Seltzer & Berger at 10GeV
!->EMcon          real*8 EpartialSC  ! default 10 (GeV)
               ! Seltzer to Tsai partial screening transition
  	       ! is made at this energy. (but normalization
	       ! is done always at 10GeV)
  	       ! Must be 100e-3 to 10 (GeV).  Around 1 GeV
               ! Selzter and Tsai gap is large for light media.
!          real*8 ElowerBndPair  ! Pair creation is neglected for Eg< this E
!          real*8 EupperBndCS    ! Compton scattering is neglected Eg> this E

!          integer IncGp  ! Include gp interactions or not.  This is copied from 
!            !       HowPhotoP which means
!            !   integer HowPhotoP      !2  if 0--> no photo hadron prod.
!            !   1--> Sofia at all E
!            !   2--> Exp. data < 2.5 GeV ; Sofia > 2.5GeV
!            !   3--> Sofia < 2.5GeV;  (rho, omega, phi) or pi0 or pi+/- at current model
!            !   4--> Exp. data < 2.5GeV   //

          logical Trace         ! how to take trace info.
          real*8 TraceErg(6)    ! a,b,c,d,e,f; energy range of trace particles 
                                ! (a, b), (c,d), (e,f) GeV are target region
          logical FreeC         ! first collision point is fixed or free
          real*8  EpsLeng       ! small distance to adjust boundary 
!          real*8 RecoilKEmin    ! Knock-on electron's E> this E is treated
!                                ! directly and not included dE/dx.
!->EMcon  real*8 Eanihi         ! Ee+ < Eanihi ==> anihilation is considered
          real*8 Ecut           ! cut off energy of brems. (obsolute)
          real*8 Escat          ! Modified scattering const. 
          integer Eabsorb       ! !0 --> photo electric absoption energy is
                                ! counted as energy loss in the media
!->EMcon          logical Photo         ! Photo-electric effect is considered or not
!          logical Knckon        ! obsolute
!          logical AngleB        ! T==>angle at bremsstrahlung is considered.
!-->EMCon  integer ALateCor      ! Angular-Lateral correlated scattering
                                ! is used or not.   0-> no cor. 1 & Moliere=f
                                ! ==> Gaussian cor.  2 ==> even if Moliere=f,
                                !   Gaussian cor. is used.
!          integer Sync          ! 0: No Synchrotron considered.
                                ! 1: energy loss by sync. for elec and 'sp'
                                !    if Magfield !=0.
                                ! 2; photon eimission is directly managed
                                !    note Eg is very low.
!->MuNucon        integer MuNI          ! 0: no muon n.i, 1: all loss --> dE/dx
                                ! 2: dE/dx(v<vmin) -->  dE/dx. v>vmin is
                                ! by stacastic process. 3: same as 2 but
                                ! n.i is treated explicitly.
!->MuNucon     integer MuBr          ! same as MuNI but for bremsstrahlung.
!->MuNucon     integer MuPr          ! same as MuNI but for pair creation
!          real*8  SyncLoop      ! electron may make a infinite loop in  'sp'.
!                                ! If it runs SyncLoop * r, we discard it,
                                ! where r is the cycrtron radius.
          real*8  KEmin         ! minimum kinetic energy of heavy particls
                                ! such as mu, pi, n, p...
                                ! if 0, appropriate value is set by ref. to
                                ! EminElec
!          integer MagPair       ! Magnetic pair production is considered or not
                                ! If 0, neglected. if 1 considered. 
                                ! It becomes effective, only if
                                ! Xai=B/Bc*Eg/me/2 > 0.05.
!->EMcon          integer Light         ! 0--> no light is treated at all  irres. of other setting
                                ! 11--> on fly ray tracing
                                ! 12--> on fly but scinti light tracing is at
                                !   end of each event
                                ! 21--> output data for Light=22
                                ! 22--> read data by Light=21 and do ray tracing 
         integer AutoEmin       ! 
                                ! 0-> no adjustment so given EminGamma,
                                ! EminElec, RecoilKEmin are used.
                                ! 1-> autoadjust Emin and recoil min energy
                                ! by using min thickness of media. 
                                ! 2--> same as above but minimum is little bit
                                !  higher. When ptcl energy becomes lower than
                                !  the min, it's residual range is computed
                                ! and if the length to the boundary is > the range,
                                ! particles is judged dead and kinetic energy is absorbed.
   		   	        ! (but pbar, e+ nbar are exeption; they are normally traced
				   ! 0 energy or anihilate) 
!->EMcon     integer StoppingPw    ! If positive value but Srim data is
                               ! available for a heavy ion,  Srim data
                               ! will be used for that heavy ion.
                               ! To force not to use the Srim data even
                               ! in such a case, make this vaule negative.
                               ! The absolue value of this means:
                               ! 1--> Effective charge with simple terms
                               !   by Pierce and Blann is used
                               !
                               ! 2--> Effective charge method with complex
                               !    terms is used for low E heavy dE/dx
                               ! Default is 1 and seems better than 2
                               ! at low energies.(Ek/n < 0.00* GeV)
                               ! However, the Srim data seems better at
                               ! such low energies.
        
!->EMcon        real(8)::SrimEmax      ! Srim data will not be used above this
                               ! energy (KE/n GeV) even if Srim data exists
                               ! Srim data dose not take into account the
                               ! restricted energy loss (always total 
                               ! energy loss including 
                               ! knockon > RecoilKEmin.  We explicitly
                               ! generate such high energy recoil.
                               ! default is 0.09 GeV/n 
!    PhitsXs is not used.  JamXs is managed by Cosmos
!        integer::PhitsXs 	! 0--> 
!                                ! 
!        integer::JamXs        ! 0--> only inelastic channel with Cosmos cross-section
!    		 	      ! 1--> use total cross-section with //

        logical::Zcorrec  ! (D=.true.)
              ! When multiple couloumb scatterign is
              ! teated with lareral displacement, distance travelled
              ! (Move.dl) is no more z-displacement (z is prticle
              ! direction.  Z-displacmentis < Move.dl. This 
              ! correction is considered, if this is .true.
        common /epTrackC/ TraceErg,
!     *  Bxu,  Byu,  Bzu, Exu, Eyu, Ezu,
     *  KEmin,
!    *  DtMax, Tcoef, Tmin,
     *  EminElec,  EminGamma,  EminH,
!	     EupperBndCS,  IncGp,
!     *  HowNormBrems, TargetElecBrems,
     *  Trace,  FreeC, EpsLeng,
!     *  Flpm, EpartialSC,
!     *       ElowerBndPair,
!        RecoilKEmin, SyncLoop,
!     *  Eanihi,
     *    Ecut, Escat, Eabsorb,
!	     Excom1, Excom2,
     *  MagField,  ElecField, TimeStruc, HowQuench,
!         EdepdEdx, 
     *  Zcorrec,
!     *  Moliere, 
						     
!        AngleB,
!!     *  ALateCor,
!         MuNI, MuBr, MuPr, 						     
!       Photo,  Sync, Knckon, 
!     * LPMeffect, MagPair,
     *   Light, StoppingPw, AutoEmin,  SrimEmax
!     *  PhitsXs, JamXs
