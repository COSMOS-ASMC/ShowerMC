module modepLightPty
  use modepLightMaxDef
  use modcsampAF
  !    Component and Media Property for light. no procedure included
  !
  type property   ! common property of light related component
     real(8)::refracIndex(3) = 0.   !   / 2.2
     !      component  refracion index n1.
     !      If blank, data given in config file will be used. 
     !      if a,b:  (n-1)x10^8 = a + b f^2
     !      if a,b,c: (n-1)x10^8 = a + b/(c_f^2).  f= 1/w.l in 1/cm
     character(len=flen)::refracFile="dummy"  !  dummy / dummy 
     !  File name which is assumed to be in the dirctory shown in
     !  "LightDir" (i.e,  the directory given by LightDir in sepicsfile).
     !  The refraction index is usally a function of the wave
     !  length. If it is known, this must be given in this file: pairs
     !  of (wl, index) in each line. Wave length unit is nm.
     !  index is unitless. If dummy or blank, no file is
     !  assumed and the value given by refracIndex is used.

     real(8)::refMode(maxSurfaces)  !  1   0   1   0   1   1    / 
     !    see below.
     real(8)::wrapper(maxSurfaces)   ! 0.80 1.445 0.8 1.445 0.8 0.8  
     !  used when wrapperReflecFile="dummy" for reflection (refMode=1)
     ! or   when wrapperRefracFile="dummy" for refraction (refMode=0)
     !      wrapperRefracFile is different from RefracFile
     !  The values specify how the each surface of the compoent is
     !   surrounded.  Normally, we don't put very thin wrapper or grease
     !   surrounding the component in the config, but we have to consider
     !   such wrapper or grease for light propagation. The values in refMode
     !   and wrapper correspond to surface # 1, 2,...  of the component
     !   in local coordinate. (i.e canonical form: .see surfaceNo.pdf):
     !   box:  #1 xy@z=0, #2 yz@x=boxa, #3 xz@y=0, #4 xz@y=boxb
     !         #5 yz@x=0, #6 xy@z=boxc.
     !   cylinder: #1 bottom, #2 top, #3 column part.
     !      If the # of surfaces is n and  less than the # of values 
     !	 given above, only first n are used. 
     !
     !  refMode: 1--> wrapper is a reflector (such as Al, Ag, white sheet..)
     !            with the reflection coeff. given by the correspoinding 
     !	       to the value in wrapper.
     !           0--> wrapper is a grease or very thin transparent material
     !	           such as clad of a fiber. It's refractive index is given
     !                in the corresponding wrapper value.
     !           -1--> no warpper; the componet is adjacent to another componet
     !	           (or void) described in congig. wrapper is not consulted.
     !              If the refraction index of that adjacent component is
     !		    0, the light is instantly  absorbed.
     !
     !  wrapper. See above.  Normally, if the value is less than 1,
     !       it should be a reflector and if it is greater than 1, it should
     !	  be a transparent thin medium. So we may live without "refMode".
     !	  But in some case, we may use 1 for both so we introduced refMode.
     !	  (For low energy X-ray, refraction index may be < 1., though
     !  we don't treat such case at present).
     !	  If the index is given, we compute the reflection or refraction 
     !	  probability  and continue the ray tracing.
     !       If the component has multi-layer clads, we should put such 
     !	  clads in config file directly with wrapper=0.
     !    If refMode is 1 and r=wrapper > 0, r-fraction of light is reflected
     !    and (1-r)-fraction of light is absorbed.
     !    If refMode is 0 and n2=wrapper> 0, reflection or refraction happenes.
     !    Their ratio is computed automaticaly using n1 and n2 and light angle
     !    and ray tracing continues.  The refracted light may enter a sencor
     !    in some case.
     character(len=flen)::wrapperReflecFile="dummy"
     !           see above. wave length vs  wrapper's reflection coef. 
     character(len=flen)::wrapperRefracFile="dummy"
     !           see above. wave length vs wrapper's refraction coef. 
     real(8)::mirror(maxSurfaces)   ! 1.0 1.0 1.0 1.0 1.0 1.0 /
     !      (see FuzzyFac,too)
     !   The angle of the reflected light is dermined by this.
     !   The values correspond to the surface # as refMode or wrapper do.
     !   If mirror is 1, the reflected light has the same angle as
     !   the incident light with respect to the normal vector of the
     !   surface at the reflection point. (complete mirror reflection).
     !     It is  assumed that  (1.-"mirror") is diffusive refl. 
     !   Normally, if the wrapper is white, it is said that diffusive
     !   (random) reflection happens so that the reflection angle becomes
     !  isotropic.  Normally, mirror may be 0 or 1.
     real(8)::fuzzy(maxSurfaces) !  0. 0. 0. 0. 0. 0. 
     !   This introduces some fuzzy factor for the reflection/refraction angle.
     !   If 0, exact Snell's law is assumed. 
     !   If > 0., it is taken as a Gaussian sprNead sigma in deg around the
     !   exact Snell's law; light angle is fuzzified by the Gaussian random
     !   factor. (Not used for diffusive reflection).

     real(8)::maxPathFac   ! 100. / 50. etc　
     !   Let maxPath=maxPathFac*(max dimension of the component) 
     !   If the light travels longer than this maxPath, ray tracing
     !   is ceased and the light is assumed to be absorbed.
     real(8)::Rayleigh =0 
     real(8)::Mie      =0
     real(8)::attenL   !    0.   / 150 etc
     !   Light attenuation length.  If attenFile is not given, this is 
     !   is used for all light wave length.   (cm)
     !   (Note: 0 is special;  no attenuation in this case) 

     character(len=flen)::attenFile="dummy"  !  dummy 
     !    File name.  The file is assumed to be in "LightDir" (see refRacFile)
     !    referred when actual ray tracing is performed.
     !    If the light attenuation length is a function of wave length,
     !    and is known, give the list of (wl, attenL) in the file
     !    one pair in one line.  Unit is (nm, cm)
     real(8)::WLS(4)=(/0., 0., 0., 0./)      ! 0 / 
     !   If 0, the attenuated light is completely absorbed.
     !   If 1  wl1 wl2 wl are given, light with wave length
     !   in (wl1, wl2) is converted into a light of wave length wl.
     !   For other  integer n > 1,
     !   subroutine epLightWLS is called with argument n.
     !    See epLigthWLS in UserHook.  The user must
     !    modify the program and put a wave lenght shifted light into
     !    stack area by calling eppush.
     ! 
     real(8)::NpPerMeV   ! 200  / 3000 etc
     !    Absolute scintillation photon number generated by MeV energy deposit.
     !    Poisson fluctuation will be introduced if the number is small.


     character(len=flen)::waveLenFile="dummy"   ! dummy / dummy
     !      File name.  The file is assumed to be in "LightDir".
     !      This is referred when ray tracing is needed.
     !      If the wave length distribution of emitted light is known, give
     !      it in the file; pair of (wl, dN/d(lamda) must be given  in each
     !      line.  Unit (cm, arbitrary). If this is dummy,
     !      no file is assumed and peakWl is used.

     real(8)::peakWL     ! 550   / 500 
     !       Effective  Wave Length (nm) of the scinti. of emitted light. 
     !      If waveLenFile is not given, this is used, but
     !      even if  waveLenFile is given, this is used for the effective
     !      wave length (peak wave length) in some case.
     !      If the max value in the file and this one differ very much, 
     !      warning may be issued.
     real(8)::quench(3) =0   
     real(8)::minmaxWL(2)   !  200  800 /  200 800
     !      (nm).  When a Cerekov light is emitted, its  wave length 
     !      distribution is 1/L^2 dL and short wave length light is
     !      more probable. However, short wave length light is normally
     !      instantly absorbed in the media, or the sensor may not have
     !      quantum efficiency to that light. These factors might be not
     !      well known so we may cut the light outside of this range
     !      for Cerenkov light.  For low energy charged particles,
     !      the upper limit may be automatically lowered.
     real(8)::ElightH, ElightL ! computed from minimaxWL. corresponding
     !      photon energ in eV. minWL-->ElightH  maxWL-->ElightL
     !      so H > L 
     real(8)::NpSample   !  1000  / 2000
     !      Used when Light=12 or 22
     !      Maximum number of photons actually traced for scinti photons
     !      generated  in a given unit cell.  If the number of photons
     !      is larger than  this, we get the final result by increasing 
     !      the result so far by a factor of RNp/NpSample, 
     !       where RNp is the total number  of photons.  
     real(8)::CellSize(3)   !  1 10 10 / 1 10 10  (let 3 values be p, q, r)
     !       When Light is 12, 21 or 22, and scintillation is to be treated
     !      (i.e, B-i=1; see manual), we store energy deposit in Edepo(:::)
     !      which corresponds to a number of cells. The unit cell size
     !      (dx,dy,dz) is determined p,q,r(%) of the component size.
     !      If the component is a box, dx = p*boxa/100, dy=q*boxb/100,
     !      dz=r*boxc/100
     !      If the component is a cylinder, dx=p*H/100, dy=q*R/100, dz=r*R/100.
     !      For example, (p,q,r)=(1,10,10), and box size (20,2,2)cm
     !      (dx,dy,dz)=(2, 2, 2) mm.  One may think 10 % is too large, but
     !      for the direction transverse to the sensor direction, normally
     !      ~10 % is OK.
     !--------------------------------------- next is for sensor
     real(8)::Qeff !  0.85  / 0.8
     !     Quantum efficiency of the sensor. Used when QeffFile is dummy
     character(len=flen)::QeffFile="dummy"   ! dummy 
     !     Quantum efficienty of the sensor as a function of wavelength.
     !    Pairs of (wavelength,  efficiency) in unit of (nm, %) must be
     !    given  in each line. If dummy, Qeff is used for all wavelength. 
     !    Note that efficiency is expressed by percentage.
     integer::mnOfScinti=10  
     !      this sensor is attached to a scintillator of which property
     !      is specified by 'mn' of countDE part of the scintillator in config.
     !      (even indirectly via light guide).  This is used to get
     !      the scintillator's  peakWL when converting charged particle
     !      effect to the number of p.e. The charged particle's direct hit is
     !      nomally converted to the number of p.e by  dE*cF/Ex where dE
     !      is the energy deposited in the sensor, Ex the energy
     !      corresponding to peakWL of the scintillator and cF the 
     !      correction factor.
     real(8)::cf = 1.0    !  1 /  see above
     integer:: attenAF  ! id no for attenutation length dist. as func. of w.l
     integer:: refracAF ! id no for refraction index distribution as func of w.l;  obtained when csampAF0  is called
     integer:: wrapperReflecAF ! id  no for wrappers reflection index distribution as func of w.l

     integer:: wrapperRefracAF ! id  no for wrappers refraction index distribution as func of w.l

     integer:: waveAF   ! id no for wave length distribution for scinti light

     integer:: invN2intAF  ! id no for integral of 1/n^2 for E to Emax
     integer:: QeffAF   ! id no for quantum efficiency as func. of w.l

  end type property

  type eachInfo  ! info of each Lcomp
     ! ===================================================================
     !    Above: input variables
     !    Below: variables computed in prog.
     !    for scintillator
     integer::nx, ny, nz  ! number of cells in x,y,z direction
     real(8)::dx, dy, dz  ! cell size
     real(8)::dmin        ! min of dx, dy, dz
     !   for every component
!     real(8):: Edepo    !  energy deposit in GeV in the component

     !   for sensor 
!     real(8)::tempPE(5)  ! used to count p.e.  when a number of
     ! photons are generated from  a source.  The pointer, pointPE, 
     ! may points to this or CumPE, depending on "light's"  subcode. 
     !  see CimPEU why (5).
!     real(8):: CumPE(5) ! count p.e by scinti, Ceren, syn, tran, direct hit
     ! by charged ptcls.   CumPE(5:5) corresponds to Edepo
!     real(8),pointer::pointPE(:)  ! see above.  When bunch of photons is generated
         ! ray tracing will need lot of time. So if the # > NpSample, we make
         ! actual ray tracing only NpSample photons. and store the result
         ! in tempPE and at the  end of bunch, we multiply the effect by
         ! photons > NpSample to get CumPE.
         ! A photon (starter of bunch photons) contains digit 10 in subcode
         ! and stoper photon (at end of bunch) contains digit 100 in subcode
         ! others have normal subcode.  When light is poped up from the stack
         ! and if it is starter, we make the pointer to direct tempPE 
         ! and clear pointPE
         ! For other photons, we simply add p.e in pintPE.
         ! When stopper photon appears, we multipy the weight it contains
         ! to poinPE and add to CumPE.  Then pointPE is made to direct
         ! to CumPE.  *** For direct hit, we don't use this method, but
         ! use total dEeff for the component.
!     real(8):: CumPET   !  sum of above CumPE.
!     real(8):: tempPC(2,3) ! similar to tempPE. this is to count photons 
         ! at the exist of each "light" component. 2 is for scinti, Cerenkov
         ! 3 is to be able to count photons at a  max of  3 surfarces of the compoent.
         ! For box, smaller 2 faces and sum of other 4 faces. for (elliptic)cyl. 
         ! 2 cross-secions and side. 
!     real(8):: CumPC(2,3) ! similar to CumPE.
!     real(8),pointer::pointPC(:,:)  ! similar to pointPE.
!     real(8):: CumPCT(3)   !  sum of above; sum( CumPC,1)
     real(8) :: refracN  ! refractive index of current light in current medium
     integer:: compno   ! component # in the config.
     integer:: mn       ! mn 

     integer:: comInfoNo  ! info commonly used is in comInfo(comInfoNo)
  end type eachInfo
  !        for each  light component, we define.
  type(property), save,target :: comInfo(maxProperties)  ! commonly used info.
  type(eachInfo), save,target :: Lcomp(maxLightComp)  ! all light related components
end module modepLightPty
