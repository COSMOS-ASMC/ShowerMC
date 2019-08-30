!          used by cEfield.f from Cosmos case.
      module modEfield
      implicit none
!         next one is only for printing
      namelist /Efparam/ HowEfield, DefofR, myEf

      type  Efield
      sequence
      real(8)::Ef(3)=0.
      real(8)::H1=0., H2=0.
      real(8)::R1=0., R2=0.
      real(8)::T1=0., T2=0.
      end type Efield

      integer,parameter:: nmaxEfield=5
!           HowEfield and myEf are namelist inputable
!           parameters
      integer,save:: HowEfield=0
! 0-->no Efield assumed
! 1-->simple Efield defined by myEf
! 2-->User supplies Efield subroutine cmyEfield.f

      character(len=1),save:: DefofR ='h' ! define distance to\
      the shower axis
!  h--> horizontal distance
!  p--> perpendicular to the shower axis.
      type(Efield),save:: myEf(nmaxEfield)
      integer,save:: nEfield = 0 ! acutal # of E flied


      logical,save:: checkH(nmaxEfield) = .false.
      logical,save:: checkR(nmaxEfield) = .false.
      logical,save:: checkT(nmaxEfield) = .false.
      logical,save:: useH=.false.
      logical,save:: useR=.false.
      logical,save:: useT=.false.

      real(8),save::Efxyz(3,nmaxEfield) !copy of
! myEf(:)%Ef(3) but in E-xyz coordinate

      end   module modEfield
      
      
!     csetIntInf:  set interaction information IntInfArray.
      module modSetIntInf
      implicit none
            ! Number of different kind of interactions
            ! considered for the current particle.
      integer,save::  NumberOfInte
!        Max number of kinds of interactions a particle can
!     take. (such as brems, knockon, anihilation)
      integer,save:: ProcessNo   ! The process really happend is the
! ProcessNo-th process in  IntInfArray.
!     To get this one, program must use modSetIntInf.
!     so there is no inteface routine to put this number 
      integer, parameter:: MaxInte = 6 

      
      type intinf          ! Interaction information
       sequence
       real(8):: thickness    ! in kg/m2 set if decay is F
       real(8):: length     ! in m, set if decay is T.  or eventually by cfixProc
       character(8):: process  ! process id string such as brems, pair
       logical:: decay          ! if decay, T, else F
!       logical:: dummy  ! to avoid warning by ifort (intinf must be 8*n bytes)
      end type intinf
!          define array of intinf
      type(intinf):: IntInfArray(MaxInte)

      end module modSetIntInf

      subroutine ciniIntInf
      use modSetIntInf
      implicit none
      processNo =  0
      NumberOfInte = 0 
      end
      subroutine cqintInfNo(num)
      use modSetIntInf
      implicit none
      integer,intent(out)::num  ! number of registered interactions.
      num =  NumberOfInte
      end

      subroutine cqintInfProcNo(num)
      use modSetIntInf
      implicit none
      integer,intent(out):: num !  process number selected among
      ! NumberOfInte interactions. This must be
      num = ProcessNo 
      end

      subroutine csetIntInf(lenOrThick, decay, procname)
      use modSetIntInf
      implicit none
!  #include  "Ztrack.h"
!  #include  "Ztrackv.h"
      
      real(8),intent(in)::  lenOrThick ! length (kg/m2) or (m)
             !  for the interaction 
      logical,intent(in):: decay   ! this acutully used not only for decay but
        ! for magnetic brems and pair cre.  for which
! length is in (m).  To convert this into  kg/m2, we
!     need medium info. which is dependent on the application.
!  (say uniform medium or earth atmosphere).      
      character(*),intent(in):: procname !  interaction process name
                  ! such as "brems" "pair" ...
!
      NumberOfInte = NumberOfInte + 1
      if(decay) then
         IntInfArray(NumberOfInte)%length = lenOrThick
      else
         IntInfArray(NumberOfInte)%thickness = lenOrThick
      endif
      IntInfArray(NumberOfInte)%decay = decay
      IntInfArray(NumberOfInte)%process = procname
      end
      subroutine cresetIntInf
!           this is intended to reset IntInfArray of the current
!       selected proces number.  Suppose, a pi+ of low
!       energy.  It can decay and can generate mu+ of some
!       kinetic enrgy.   So if the pi+ of such one comes to
!       stop, it is better to make it decay rather than
!       colide. This is to reset the present process and 
!     to force the decay.
      use modSetIntInf
      implicit none
!   #include  "Ztrack.h"
!  #include  "Ztrackv.h"
      
      IntInfArray(ProcessNo)%decay = .true.
      IntInfArray(ProcessNo)%process = 'decay'
      end

      module modMuNucontrol
      implicit none
      logical:: IncMuonPolari     !1  if T, consider muon polarization
      integer:: MuNI=3             !2  0$\rightarrow$ nuclear interaction of muon is completely
           !        neglected \newline
           !  1 $\rightarrow$ energy loss by n.i is subsumed in dE/dx of muons as a continuous
           !      energy loss.
           !      Let v=Etransfer/Emu,  the loss here is Int(vc:vmax) of (Emu vdsigma/dv).
           !      (vc $\sim$0, vmax$sim$1). \newline
           !  2 $\rightarrow$ (Default value). similar to 1 but as the continuous loss only
           !      v $<$vmin=10$^{-3}$ of fractional muon energy is subsumed (Int(vc: vmin) of
           !      (Emu vdsigma/dv)).  The portionof loss by v$>$vmin is treated as a stocastic 
           !       process.  However, the product from the n.i itself is neglected \newline
           !      3$\rightarrow$ the same as 2, but the n.i is explicitly included to produce 
           !      a number of particles.
           !     The n.i is treated as a photo-nucleus interaction.
      integer:: MuBr=3      !2  parameter similar to MuNI but for bremsstrahlung by muons.
      integer:: MuPr=3            !2  parameter similar to MuNI but for pair creation by muons.     
      
      real(8)::  MuonPolarization ! muon polarization value. not parameter

      end module modMuNucontrol

      
      module modEMcontrol
      ! many are moved from Ztrackp.h
      implicit none
      logical,save:: LPMeffct=.true. ! if this is f, LPM is negelected
      logical:: LpmEffect
      equivalence (LPMeffect, LpmEffect)
!               even if energy and density of medium is high enough.
!           if T, LPM possibilty is tested and if yes, LPMworks
!     will be set to T.   This is a member of namelist items.
      real(8):: Flpm =1.0d0     !  Ee > Flpm*BremEminLPM is the actual LPM region
      integer:: HowNormBrems=-1   ! -1 (default) : normalization base is C.S
                                ! 0: no normalization
                                ! 1: base is Seltzer & Berger at 10GeV
      logical,save:: LPMworks   ! if current e or g's energy is high
                 !  and media's rho is high, and LPM seems
!          to effective and LPMeffect is T, then
!     this becmoes T.   The value is set to T or F when the interaction
!     length sampling is tried. So if it is set T or F and, if e.g, brems 
!     process is selected for the electron interaction,  brems energy
!     sampleing should be perfomed referring to LPMworks (i.e,
!     with or without the LPM effect_
!     

!     There is some delicate but negligle effect due to energy loss during
!     tracking:   LPMworks may be set to T when brems interaction point is
!     sampled.  But after it is moved that point, it loses energy by dE/dx
!     and its energy may drop and escape from the LPM energy region. However,
!     in such a case, LPM effect is very weak so we may sefely neglect
!     the effect.

!
      real(8),save::LpmBremEmin = 0. ! GeV 10^16 eV.
!     If 0, after reading parameter, replaced by 10^7 (GeV.)         
!     This is used only in Cosmos (cepsampEIntL.f)  to  estimate quickly
!     that we should  consider the LPM  effect for electron brems.
!     More strict judgement is later done using air density at the current
!     electron position (--> LPMworks will be set to  T/F).
!     Gettign air density needs some computation time, we avoid such
!     calculations for  many low energy electrons.
!     If you change the air medium, this   might be better changed.
!     This is in the namelist items.

      real(8),save::LpmPairEmin = 0. ! GeV 10^18 eV.
!       same as LpmBremEmin.
      real(8),save::MagPairEmin = 2.0d10 !  GeV. 2x10^19 eV
!     magnetic  pari creation is considered above this energy.
!     namelist item.
      real(8),save::MagBremEmin = 2.0d10 !  GeV. 2x10^19 eV

      integer:: HowMagField=0     ! how to give magnetic field
      integer:: HowElecField=0         ! //          electric field
      real(8)::  Bxu=0.,  Byu=0.,  Bzu=0.   ! uniform mag field comp.
      real(8)::  Exu=0.,  Eyu=0.,  Ezu=0. ! uniform elec field comp.
      real(8)::  DtMax         ! max allowable path in r.l
      real(8)::  Tcoef=5.d0         !
      real(8)::  Tmin=1.d-3         !
      
      real(8),save:: rrho        !  actual(density)/media%rho
!        in Epics, this is media%rhoc.
!             when  a media is used, we can specify the ratio of (its
!            density)/(standard predefined density) as: e.g
!               Air*0.75
!            which specifies that the air density is 0.75 times of the
!            standard one.   and 0.75 becomes rhoc.      
!        in Cosmos, this is rho(@ presnet particle location)/media%rho
!         (rhoc should be default i.e, 1 in this case)
!     rrho  is calculated when LPMworks is calculated. for e/g.
      integer,save:: Sync=0            ! For Epics;  2-> Cosmos
      integer,save:: SyncLoop=10.d0    !  //
      integer:: MagBrem
      
      equivalence (Sync, MagBrem) 
      integer:: MagPair=0       ! For Epics, 1->Cosmos
!      real(8)::  MagChgDist =20.d3   !2 Distance where mag. can be seen  as const.(m) at sea level --> modAtmosDef

      real(8),save:: Xai !  output from path sampling (epmpairp)
!     and used in energy sampling (epmpaire)  Eg/Me *Bsin/Bcr/2
      real(8):: RecoilKineMinE
      real(8),save:: Upsilon    ! Ee/Me*Bsin/Bcr for syncrotron      　　
      real(8),save:: RecoilKEmin=100.0d-6 ! 100 keV
      equivalence (RecoilKEmin,RecoilKineMinE)
      logical,save:: AngleB=.true.   ! Brems Angle sampling is tried if 
      real(8),save:: EemaxBremAngle=10. !  Ee  <  this.  (GeV).
!      real(8)::  EminElec          ! e-'s E < EminElec is discarded 
!      real(8)::  EminGamma         ! Egamma < EminGamma is  discarded
      

      logical::Knockon=.true.
      real(8)::KnockOnRatio =1.0d0 
      real(8):: AnihiE
      real(8):: Eanihi=10.d-3
      equivalence ( AnihiE, Eanihi)
      real(8):: EpartialSC=10.0d0      ! default 10 (GeV)
!       Seltzer to Tsai partial screening transition
!       is made at this energy. (but normalization
!       is done always at 10GeV)
!       Must be 100e-3 to 10 (GeV).  Around 1 GeV
                     ! Selzter and Tsai gap is large for light media.

      real(8):: EupperBndCS=10  ! GeV ! Compton scattering is neglected Eg> this E
      real(8):: Ecut            ! cut off energy of brems. (obsolute)
      real(8):: Escat=19.3d-3   ! Modified scattering const.
      integer:: Eabsorb=6       ! In older versions of Cosmos Eabsorb(1).
         !   when a charged particle makes energy loss to Air, chookEabsorb is
         !   always called. When a particle dies (i.e, it K.E becomes < Emin) chookEabsorb2
         !   may be called dedpending on the particle and Eabsorb bit. (LSB is bit 1).
         !   calling cond (See Ztrackv.h, BitPhoton etc).
         !   bit 2  and photon.  bit 3 and e-/e+. bit 4 proton. bit 5 neutron
         !   bit 6 anti-N.  bit 7 decayable ptcl, bit8 others
         !   bit 1 is for photoelectric effect but is not used in Air.

         !  Comos Eabsorb(2) is now EabsorbL but not used in Epics
      
      integer:: EabsorbL=0      !  Used in Cosmos.  If 0, it will be made to be the
             ! # of Observation sited.  If > 0, and <= NoOfSites, it  is used.
             ! otherwise error stop.

      real(8):: MaxComptonE 
      equivalence (MaxComptonE, EupperBndCS)
      
      real(8):: ElowerBndPair=1.022d-3 ! GeV  Pair creation is neglected for Eg< this E
      logical:: Photo=.true.   ! photo electric effect
      real(8):: Excom1=1.d-3    ! if Eg <=Excom1,  XCOM data is used for
                             ! compton/p.e/coherent .  (<100)
      real(8):: Excom2=1.0d0    ! if Eg <= Excom2, XCOM data for Pair cre. (<100)
      logical:: EdepdEdx=.true.     ! use energy dependent dE/dx
             !  F is used to compare analytic calc
      integer:: HowQuench       ! R: dE/dx restricted
          !     F: dE/dx Full.  D delta ray contribution
          !                     O other // (knot known).
          !  cf: exprimetally obained quenching coef.
          ! ecf: quenching factor for R. D+O has no quenching
          !  ecf*R + D + O = cf*F; D=F-R
          !  ecf + F/R -1 + O/R = cf*F/R; then
          !  ecf = 1 + (cf-1)F/R - O/R ; set O= 0 for first approx.
          !   if  HowQuench = 0  use cf for ecf (equiv. to  v <=9.154
          !   if            = 1  default:  use ecf as if old cf
      integer:: TargetElecBrems=1  ! Flag to include target electron
                 !    brems 0->not bit 1->mu  bit2 p,pi,K (bit3 heavy?)

      integer:: IncGp=1         ! Include gp interactions or not.  This is copied from

      real(8):: MaxPhotoE=1.0d-3       !  Above this energy, photoelectric effect is neglected.

      
      integer:: HowPhotoP = 1
      real(8):: MinPhotoProdE = 153d-3
      equivalence (HowPhotoP, IncGp)
!            !       HowPhotoP which means
!            !   integer HowPhotoP      !2  if 0--> no photo hadron prod.
!            !   1--> Sofia at all E
!            !   2--> Exp. data < 2.5 GeV ; Sofia > 2.5GeV
!            !   3--> Sofia < 2.5GeV;  (rho, omega, phi) or pi0 or pi+/- at current model
!     !   4--> Exp. data < 2.5GeV   //
      integer:: Moliere=2
!           1==> simplified faster Moliere scattering
!           2==>  rigorous Moliere scattering. takes longer time than 1.
!           3==>  ?
      integer:: ALateCor=1      ! Angular-Lateral correlated scattering
         ! is used or not.
         ! 0==> no cor.      
         ! 1 & Moliere=1  ==> Gaussian cor.
         ! 2 ==> even if Moliere=f, Gaussian cor. is  forced.
      real(8):: Es = 19.3d-3    ! modified scattering constGeV.
      integer:: Light ! 0--> no light is treated at all  irres. of other setting
          ! 11--> on fly ray tracing
          ! 12--> on fly but scinti light tracing is at
          !   end of each event
          ! 21--> output data for Light=22
          ! 22--> read data by Light=21 and do ray tracing  
      integer:: StoppingPw=1        ! If positive value but Srim data is
          ! available for a heavy ion,  Srim data
          ! will be used for that heavy ion.
          ! To force not to use the Srim data even
          ! in such a case, make this vaule negative.
          ! The absolute value of this means:
          ! 1--> Effective charge with simple terms
          !   by Pierce and Blann is used
      real(8)::SrimEmax=0.09d0     ! Srim data will not be used above this
      ! energy (KE/n GeV) even if Srim data exists
      ! Srim data dose not take into account the
      ! restricted energy loss (always total
      ! energy loss including
      ! knockon > RecoilKEmin.  We explicitly
      ! generate such high energy recoil.
                                     ! default is 0.09 GeV/n

      logical::Zcorrec=.true.        ! (D=.true.)
        ! When multiple couloumb scatterign is
        ! teated with lareral displacement, distance travelled
        ! (Move.dl) is no more z-displacement (z is prticle
        ! direction.  Z-displacmentis < Move.dl. This
        ! correction is considered, if this is .true.
      
      end module modEMcontrol

      module modColInfo
      implicit none
! #include "Zmaxdef.h"
! #include "Zptcl.h"
!          media info
!     integer,save:: TargetMediaNo ! hadronic collision media is Media(TargetMeiaNo) =>media
!              this is same as MediaNo in ZmediaLoft.h      
      integer,save:: colElemNo !among them, media%elem(colElemNo) is the elelment.
      integer,save:: TargetNucleonNo ! its mass # A
      integer,save:: TargetProtonNo ! its charge # Z
      real(8),save:: TargetXs   !  x-section of (A,Z); mb

!!!!  below are moved to new Zpwork.h so two #include not needed now
!       info of  produced  particles
                    ! max # of ptcls producable in coll.
!      integer,parameter::  MaxPtcl = MAX_PTCL

    
!      type(ptcl):: Pwork(MaxPtcl) ! working array to store ptcls.
!      integer:: Nproduced         ! no. of ptcls produced and stored in Pwork.
!      integer:: Nstacked          ! no. of ptcls stacked. If ThinSampling=F,
! same as Nproduced.  Nstacked <= Nproduced
      end     module modColInfo


      
