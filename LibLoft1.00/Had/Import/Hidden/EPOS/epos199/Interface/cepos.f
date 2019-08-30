!     Make TUNE4M  1 to activate a workaround for the EPOS 4 momentum problem
!     !=1  ==>  output from EPOS will be accepted as it is; however, in this
!       case,  the relation,  E^2= p^2+ m^2 of each particls, is
!       sometimes broken non negligibly.    (see  eposWorkaround.pdf)
!     ==1  ==> The  energy of each particle is forced to be
!          sqrt(p^2 + m^2).  This workaround is used in CRMC.


#define TUNE4M  1
      module modEpos
      implicit none
#include "Zptcl.h"
      type(ptcl):: cms
      character*10,save:: from
      end module modEpos

      subroutine ceposIniAll
!        init of  EPOS  (once for all)

      implicit none
      call aaset(0)   ! default init for all events
!      call LHCparameters  ! needed for epos LHC version
      call ceposini00  ! some change from the default
 !          use ainit to specify r max target and projectile
      call cinieposmax
      end
      subroutine cinieposmax
      use modEpos
      include "../epos.inc"

      call cqfrom(from)  ! if called from Gencol, from = 'gencol' 

      idprojin = 1120
      idtargin = 1120
!      elab = 200.
      engy   = -1.
      pnll = 200.
      ekin = -1.
      ecms =  -1.
!      iframe=12                 ! 12=target frame (needed ?)
      if( from == 'gencol' ) then
         call cqGencolCMS(cms)
      endif
      if( from == 'epics'  .or. from ==  'gencol') then
         maproj = 56
         laproj = 26
         matarg = 14
         latarg = 1
      else
         maproj = 56
         laproj = 26
         matarg = 14 
         latarg = 1  
      endif
      call ainit
      end
      subroutine cqfrom(from)
#include "Zmanager.h"
      character(*),intent(out):: from 
      from = CosOrEpi
      end

      subroutine ceposIniOneEvent(pjin, tg, sig)
!        This is called from dummy init from ceposIniAll
!        and when getting xs for a given media to calculate mfp.
!        and when the target is fixed and event generation just starts.
!         (from within cepposGenOneEvent)
#include "Zcode.h"
#include "Zptcl.h"
      include "../epos.inc"
      type(ptcl):: pjin  ! input projectile particle
      type(ptcl):: tg  ! input target particle in lab.
      real(8),intent(out):: sig  ! inelastic cross section in mb

      type(ptcl)::pj
      real(8):: u

      pj = pjin

      elab = pj%fm%p(4)
      engy = -1.
      pnll = -1.
      ekin = -1.
      ecms =  -1.
      iframe=12            ! 12=target frame ; needed ?

      if( pj%code == kgnuc )  then
         maproj = pj%subcode           !proj A
         laproj = pj%charge
         idprojin = 1120
         elab = elab/pj%subcode   ! E(GeV)/n
      else
         if( pj%code == kphoton)  then
         !  replace it to pi0 or eta
            call rndc(u)
            if( u < 0.5d0 ) then
               call cmkptc(kpion, 1,  0, pj)
            else
               call cmkptc(keta, 1,  0, pj)
            endif
            ! adjust momentum
            cf = sqrt(1.d0 - (pj%mass/pj%fm%p(4))**2 )
            pj%fm%p(1:3) = cf * pj%fm%p(1:3)
         endif   
         call ccos2eposB( pj, idprojin)
                ! not idporj; before ainit "in"  must be added
         if( abs(idprojin) == 20  ) then ! k0s(20) or k0l(-20)
            call rndc(u)
            if(u < 0.5d0 ) then
               idprojin = 230   ! k0 
            else
               idprojin = -230  ! k0bar
            endif
         elseif(abs(idprojin) == 2130 ) then ! Lambda/lambda-bar
            idprojin = sign(1230, idprojin) ! use sigma/sigma-ba
         else
               !     use idprojin as it is;  
         endif
         laproj = -1
         maproj = 1
      endif
!/////////////
!      write(0,*) ' idprojin, laproj, maproj= ',
!     *        idprojin, laproj, maproj
!      write(0,*) ' idprojin has been changed to', idprojin
!//////////


      if( tg%code == kgnuc ) then
         idtargin = 1120
         latarg = tg%charge     !targ Z
         matarg = tg%subcode    !targ A
      else  ! p or n
         call ccos2eposB( tg, idtargin)
         matarg = 1           
         laproj = -1
      endif
!////////////
!      write(0,*) ' idtargin, latarg, matarg= ',
!     *  idtargin, latarg, matarg
!      write(0,*) ' entering  ainit'
!/////////////
      call ainit
      sig = sigineaa
      end
      
      subroutine ceposGenOneEvent(pj, ia, iz, a, n)
      use modEpos
!!!      implicit none ! cannot be used due to epos.inc
#include "Zcode.h"
! #include "Zptcl.h"

      type(ptcl):: pj   !  inp. projectile
      integer,intent(in):: ia  ! target A
      integer,intent(in):: iz  ! target Z
      type(ptcl):: a(*)  !  generated ptcl. must be >= n
      integer,intent(out):: n  ! # of  ptcls
      
      integer:: i, j
      include "../epos.inc"
      integer code, subcode, charge, kf
      type(ptcl):: tg
      real(8):: xs

      if( ia > 1 ) then
         call cmkptc(kgnuc, ia, iz, tg)
      else
         call cmkptc(knuc, ia, iz, tg)
      endif         
      tg%fm%p(1:3) = 0.
      tg%fm%p(4) = tg%mass
      call ceposIniOneEvent(pj, tg, xs)


      call aepos(-1)   ! generate 1 event -1 or 1 ?? 
              ! corsika use 1.
!       Fix final particles and some event parameters; epos code obtained
      if( from /= 'gencol') then
         call afinal   ! to the original system (NG at high E)
      endif
!       convert to HEP code 
      call ustore   !!!
      if(nhep .gt. nmxhep)then
         write(0,*) 'Error: produced number of particles=', nhep

         write(0,*) '>  nmxhep =',nmxhep
         stop
      endif
      
      n = 0
      do i = 1, nhep
#if  TUNE4M == 1
!          force E=sqrt( m^2+p^2) as in CRMC
         phep(4,i)=sqrt( dot_product(phep(1:3,i), phep(1:3,i))
     *      +   phep(5,i)**2  )
#endif    

         kf = idhep(i)
         call ckf2cos(kf, code, subcode, charge)
!///////////
!         write(0,
!     *  '(i3, a, i10, a, i3, a, i3, a, i4, a, i4, a, 1p,4g13.3)')
!     *  i, ' idhep=', idhep(i), ' code=',code, 
!     *  ' status=', isthep(i), ' moth=', jmohep(1,i), ' daug=',
!     *  jdahep(1,i), ' p(1:4)= ', phep(1:4, i)
!///////////
         if( isthep(i) /= 4 ) then
            if( code == krare .and. kf >= 1000000020 ) then
!                 Z= 0 A=subcode;  nucleus
               if( phep(4,i)/subcode > phep(5,i)*2 ) then
                  ! energy/n  > 2mass ; issue warning
                  write(0,*) 'from epos: idhep=',kf
                  write(0,*) ' subcode=',subcode,' charge=',charge
                  write(0,*) ' ia = ', ia, ' iz=',iz
                  write(0,*) ' 4mom=',phep(1:4,i) 
                  if( pj%code == kgnuc ) then
                  !  assign charge; not so bad method
                     call csetFragChg(
     *                    int(pj%subcode), subcode, charge)
                     code = kgnuc
                     n = n + 1
                     call cmkptc(code, subcode, charge, a(n))
                     a(n)%fm%p(1:4) = phep(1:4,i)
                     if( phep(4,i)/subcode > phep(5,i)*10 ) then
                        write(0,*) ' accepted as charge=',charge
                        write(0,*) ' nucleus'
                     endif
                  endif
               endif
            else
               n = n + 1
               call cmkptc(code, subcode, charge, a(n))
               a(n)%fm%p(1:4) = phep(1:4,i)
               !    use mass in phep
               a(n)%mass = phep(5,i)
            endif
         endif
      enddo
      if( from == 'gencol') then
         ! a is still in cms, so boost to lab system
         do i = 1, n
            call cibst1(i, cms, a(i), a(i))
         enddo
      endif
      call crot3mom(pj, a, n)  ! rotate result
      end


      subroutine ceposini00
!          change some of the parameters from default.  
!            (modification of IniEpos prepared for LHCf)
!      implicit none cannot be used because of epos.inc 
#include "Zkfcode.h"
      include "../epos.inc"

      integer:: idtrafo  ! pdg<--->epos code conversion
      

      integer:: kgetenv2 ! func to get Environmental variable
      integer:: leng     ! to get length of the string
      real(8):: u
      call rndc(u)
      seedi=u   !seed for random number generator: at start program
      call rndc(u)
      seedj=u   !seed for random number generator: for first event
      iwseed = 0   ! no record seed

! Initialize decay of particles
      nrnody=0       !number of particle types without decay (if 0 (default) : all unstable particles decay (at the end only (anti)nucleons, (anti)electrons and muons)
! Particle code is given as
!     id=+/-ijkl
!
!          mesons--
!          i=0, j<=k, +/- is sign for j
!          id=110 for pi0, id=220 for eta, etc.
!
!          baryons--
!          i<=j<=k in general
!          j<i<k for second state antisymmetric in (i,j), eg. l = 2130
!
!          other--
!          id=1,...,6 for quarks
!          id=9 for gluon
!          id=10 for photon
!          id=11,...,16 for leptons
!          i=17 for deuteron
!          i=18 for triton
!          i=19 for alpha
!          id=20 for ks, id=-20 for kl
!
!          i=21...26 for scalar quarks
!          i=29 for gluino
!
!          i=30 for h-dibaryon
!
!          i=31...36 for scalar leptons
!          i=39 for wino
!          i=40 for zino
!
!          id=80 for w+
!          id=81,...,83 for higgs mesons (h0, H0, A0, H+)
!          id=84,...,87 for excited bosons (Z'0, Z''0, W'+)
!          id=90 for z0
!
!          diquarks--
!          id=+/-ij00, i<j for diquark composed of i,j.
!
! Examples : 2130 = lambda, 1330=xi0, 2330=xi-, 3331=omega
!
! Conversion from epos to  pdg code can be done using
!      id_pdg=idtrafo('nxs','pdg',id_epos)

      nrnody=nrnody+1
      nody(nrnody)=120     !pi+

      nrnody=nrnody+1
      nody(nrnody)=-120    !pi-
      nrnody=nrnody+1

      nody(nrnody)=130     !K+
      nrnody=nrnody+1

      nody(nrnody)=-130    !K-
      nrnody=nrnody+1

      nody(nrnody)=-20     !Kl

      nrnody=nrnody+1
      nody(nrnody)=-14     !mu+
      nrnody=nrnody+1
      nody(nrnody)=14      !mu-

!         we inhibit decay of  pi0, Lambda0,  k0s, sigma, (gzai ?)
      nrnody=nrnody+1
      nody(nrnody)=idtrafo('pdg','nxs',kflambda)    !lambda using pdg code

      nrnody=nrnody+1
      nody(nrnody)=idtrafo('pdg','nxs',kfpi0)    !pi0 

      nrnody=nrnody+1
      nody(nrnody)=idtrafo('pdg','nxs',kfk0s)    !k0short 

      nrnody=nrnody+1
      nody(nrnody)=idtrafo('pdg','nxs', kfsigmam)    ! sigma-

      nrnody=nrnody+1
      nody(nrnody)=idtrafo('pdg','nxs',kfsigmap)    !sigma+

      nrnody=nrnody+1
      nody(nrnody)=idtrafo('pdg','nxs',kfsigma0)    !sigma0

      nrnody = nrnody + 1
      nody(nrnody) = idtrafo('pdg', 'nxs', kfgzai0) ! Xi0
      nrnody = nrnody + 1
      nody(nrnody) = idtrafo('pdg', 'nxs', kfgzai) ! Xi-
      nrnody = nrnody + 1
      nody(nrnody) = idtrafo('pdg', 'nxs', -kfgzai) ! Xi+
      nrnody=nrnody+1
      nody(nrnody)=1220     ! n
      nrnody=nrnody+1
      nody(nrnody)=-1220     ! an
      nrnody = nrnody + 1
      nody(nrnody) = idtrafo('pdg', 'nxs', kfeta) !  eta


      isigma=0   !do not print out the cross section on screen
!      ionudi=3   !count diffraction without excitation as elastic
      ionudi=1  ! this is used in  Corska.  include quasi elastic events but strict calculation of xs
      iorsce=0  ! used in corsika. color exchange turned on(1) or off(0)
      iorsdf=3  !corsika. droplet formation turned on(>0) or off(0)        
      iorshh=0  !corsika. other hadron-hadron int. turned on(1) or off(0) 
      istore = 0  !corsika  DO NOT STORE EVENTS ON zzz.data FILE
!      iframe=11                 !nucleon-nucleon frame (12=target)
      iframe=12                 ! 12=target frame
      iecho=0                     !"silent" reading mode

!    infragm= 0  ???????????

!      fnnx="./"                    ! path to main epos subdirectory
      leng=kgetenv2("LIBLOFT", fnnx)
      if( leng == 0 ) then
         write(0,*) 'env.  LIBLOFT not given '
         stop
      endif
!          file management
      fnnx = fnnx(1:leng)//"/Had/Import/EPOS/"
      nfnnx=len(trim(fnnx))  ! length of fnnx
!            files are opened and closed by epos prog.
!        with file # =1. 
!      nfnii=10                     ! epos tab file name length
!      fnii="epos.initl"            ! epos tab file name
!       files are  used in epos-sem except fnie which is in ep-qsh 
      fnii = trim(fnnx)//"epos.initl"
      nfnii = len(trim(fnii))

      fnid=trim(fnnx)//"epos.inidi"
      nfnid=len(trim(fnid))
!         used in ep-qsh
      fnie=trim(fnnx)//"epos.iniev"
      nfnie=len(trim(fnie))

      fnrj=trim(fnnx)//"epos.inirj"
      nfnrj=len(trim(fnrj))

      fncs=trim(fnnx)//"epos.inics"
      nfncs=len(trim(fncs))

! Debug
      ish=0       !debug level
      ifch=0      !debug output (screen); stderr
      ifmt = 0    ! r messages
!      ifch=31    !debug output (file)
!      fnch="epos.debug"
!      nfnch=index(fnch,' ')-1
!      open(ifch,file=fnch(1:nfnch),status='unknown')

!       These are postponed until init for each event.
!      nevent = 1  !number of events
!      modsho = 1  !printout every modsho events
!
!      ecms=14000  !center of mass energy in GeV/c2
!      
!      idproj = 1120   !proton
!      laproj = 1      !proj Z
!      maproj = 1      !proj A
!      idtarg = 1120   !proton
!      latarg = 1      !targ Z
!      matarg = 1      !targ A
!
      istmax = 0      !only final particles (istmax=1 includes mother particles)

! for main program
!      nevto  = nevent
!      isho   = ish

      end



!-----------------------------------------------------------------------     
      subroutine EposInput(nevto,isho)
!-----------------------------------------------------------------------     
! Read informations (new options or parameter change) in the file
! "epos.param". The unit "ifop" is used in aread. If not used, it will
! use the default value of all parameters.
!-----------------------------------------------------------------------     
! for TempDev
#include "Zmanagerp.h"   
      include "../epos.inc"
      nopen=0
!      ifop=35
      ifop = TempDev
      open(unit=ifop,file='example.param',status='old')
      call aread
      close(ifop)
! for main program
      nevto  = nevent
      isho   = ish
      end

!-----------------------------------------------------------------------
      function rangen()  result(rn)   ! single precision
!-----------------------------------------------------------------------
!     generates a random number
!         use Cosmos generator 
      implicit none
      real(8):: u
      real(4):: rn
      call rndc(u)
      rn = u
!-----------------------------------------------------------------------
!      include 'epos.inc'
!      double precision dranf
! 1    rangen=sngl(dranf(dble(rangen)))
!      if(rangen.le.0.)goto 1
!      if(rangen.ge.1.)goto 1
!      if(irandm.eq.1)write(ifch,*)'rangen()= ',rangen
!
!      return
      end

!-----------------------------------------------------------------------
!      double precision function drangen(dummy)
      function drangen(dummy)  result(u)
!-----------------------------------------------------------------------
!     generates a random number
!-----------------------------------------------------------------------
!      include 'epos.inc'
!      double precision dummy,dranf
!      drangen=dranf(dummy)
!      if(irandm.eq.1)write(ifch,*)'drangen()= ',drangen
!
      real(8):: dummy
      real(8):: u
      call rndc(u)
      
      end
!-----------------------------------------------------------------------
      function cxrangen(dummy)  result(rn)
!-----------------------------------------------------------------------
!     generates a random number
!-----------------------------------------------------------------------
!      include 'epos.inc'
!      double precision dummy,dranf
!      cxrangen=sngl(dranf(dummy))
!      if(irandm.eq.1)write(ifch,*)'cxrangen()= ',cxrangen
!
      real(8):: dummy
      real(4):: rn
      real(8):: u
      call rndc(u)
      rn = u
      end

! Random number generator from CORSIKA *********************************




!C=======================================================================
!
!      DOUBLE PRECISION FUNCTION DRANF(dummy)
!
!C-----------------------------------------------------------------------
!C  RAN(DOM  NUMBER) GEN(ERATOR) USED IN EPOS
!C  If calling this function within a DO-loop
!C  you should use an argument which prevents (dummy) to draw this function 
!C  outside the loop by an optimizing compiler.
!C-----------------------------------------------------------------------
!      implicit none
!      integer irndmseq
!      double precision uni,dummy
!C-----------------------------------------------------------------------
!
!      call RMMARD( uni,1,irndmseq)
!
!      DRANF = UNI
!      UNI = dummy        !to avoid warning
!
!      RETURN
!      END


!-----------------------------------------------------------------------
      subroutine ranfgt(seed)
!-----------------------------------------------------------------------
! Initialize seed in EPOS : read seed (output)
! Since original output seed and EPOS seed are different,
! define output seed as : seed=ISEED(3)*1E9+ISEED(2)
! but only for printing. Important values stored in /eporansto/
! Important : to be call before ranfst
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER          KSEQ
      PARAMETER        (KSEQ = 2)
      COMMON /CRRANMA3/CD,CINT,CM,TWOM24,TWOM48,MODCNS
      DOUBLE PRECISION CD,CINT,CM,TWOM24,TWOM48
      INTEGER          MODCNS
      COMMON /CRRANMA4/C,U,IJKL,I97,J97,NTOT,NTOT2,JSEQ
      DOUBLE PRECISION C(KSEQ),U(97,KSEQ)
      INTEGER          IJKL(KSEQ),I97(KSEQ),J97(KSEQ),
     *                 NTOT(KSEQ),NTOT2(KSEQ),JSEQ
      common/eporansto/diu0(100),iiseed(3)
      double precision    seed,diu0
      integer iiseed,i

      iiseed(1)=IJKL(1)
      iiseed(2)=NTOT(1)
      iiseed(3)=NTOT2(1)
      seed=dble(iiseed(3))*dble(MODCNS)+dble(iiseed(2))
      diu0(1)=C(1)
      do i=2,98
        diu0(i)=U(i-1,1)
      enddo
      diu0(99)=dble(I97(1))
      diu0(100)=dble(J97(1))
      return
      end

!-----------------------------------------------------------------------
      subroutine ranfst(seed)
!-----------------------------------------------------------------------
! Initialize seed in EPOS :  restore seed (input)
! Since original output seed and EPOS seed are different,
! define output seed as : seed=ISEED(3)*1E9+ISEED(2)
! but only for printing. Important values restored from /eporansto/
! Important : to be call after ranfgt
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER          KSEQ
      PARAMETER        (KSEQ = 2)
      COMMON /CRRANMA3/CD,CINT,CM,TWOM24,TWOM48,MODCNS
      DOUBLE PRECISION CD,CINT,CM,TWOM24,TWOM48
      INTEGER          MODCNS
      COMMON /CRRANMA4/C,U,IJKL,I97,J97,NTOT,NTOT2,JSEQ
      DOUBLE PRECISION C(KSEQ),U(97,KSEQ)
      INTEGER          IJKL(KSEQ),I97(KSEQ),J97(KSEQ),
     *                 NTOT(KSEQ),NTOT2(KSEQ),JSEQ
      common/eporansto/diu0(100),iiseed(3)
      double precision    seedi,seed,diu0
      integer i,iiseed

      seedi=seed
      IJKL(1)=iiseed(1)
      NTOT(1)=iiseed(2)
      NTOT2(1)=iiseed(3)
      C(1)=diu0(1)
      do i=2,98
        U(i-1,1)=diu0(i)
      enddo
      I97(1)=nint(diu0(99))
      J97(1)=nint(diu0(100))
      return
      end

!-----------------------------------------------------------------------
      subroutine ranflim(seed)
!-----------------------------------------------------------------------
      double precision seed
      if(seed .gt. 1d9)stop'seed larger than 1e9 not possible !'
      end

!-----------------------------------------------------------------------
      subroutine ranfcv(seed)
!-----------------------------------------------------------------------
! Convert input seed to EPOS random number seed
! Since input seed and EPOS (from Corsika) seed are different,
! define input seed as : seed=ISEED(3)*1E9+ISEED(2) 
!-----------------------------------------------------------------------
      IMPLICIT NONE
      COMMON /CRRANMA3/CD,CINT,CM,TWOM24,TWOM48,MODCNS
      DOUBLE PRECISION CD,CINT,CM,TWOM24,TWOM48
      INTEGER          MODCNS
      common/eporansto/diu0(100),iiseed(3)
      double precision    seed,diu0
      integer iiseed

      iiseed(3)=nint(seed/dble(MODCNS))
      iiseed(2)=nint(mod(seed,dble(MODCNS)))

      return
      end

!-----------------------------------------------------------------------
      subroutine ranfini(seed,iseq,iqq)
!-----------------------------------------------------------------------
! Initialize random number sequence iseq with seed
! if iqq=-1, run first ini
!    iqq=0 , set what sequence should be used
!    iqq=1 , initialize sequence for initialization
!    iqq=2 , initialize sequence for first event
!-----------------------------------------------------------------------
      IMPLICIT NONE
      COMMON /CRRANMA3/CD,CINT,CM,TWOM24,TWOM48,MODCNS
      DOUBLE PRECISION CD,CINT,CM,TWOM24,TWOM48
      INTEGER          MODCNS
      common/eporansto/diu0(100),iiseed(3)
      double precision    seed,diu0
      integer iiseed
      common/eporansto2/irndmseq
      integer irndmseq
      integer iseed(3),iseq,iqq,iseqdum
 
      if(iqq .eq.0)then
        irndmseq=iseq
      elseif(iqq .eq.-1)then
        iseqdum=0
        call RMMAQD(iseed,iseqdum,'R')   !first initialization
      elseif(iqq .eq.2)then
        irndmseq=iseq
        if(seed .ge. dble(MODCNS))then
           write(*,'(a,1p,e8.1)')'seedj larger than',dble(MODCNS)
           stop 'Forbidden !'
        endif
        iiseed(1)=nint(mod(seed,dble(MODCNS)))
! iiseed(2) and iiseed(3) defined in aread
        call RMMAQD(iiseed,iseq,'S') !initialize random number generator
      elseif(iqq .eq.1)then        !dummy sequence for EPOS initialization
        irndmseq=iseq
        if(seed .ge. dble(MODCNS))then
           write(*,'(a,1p,e8.1)')'seedi larger than',dble(MODCNS)
           stop 'Forbidden !'
        endif
        iseed(1)=nint(mod(seed,dble(MODCNS)))
        iseed(2)=0
        iseed(3)=0
        call RMMAQD(iseed,iseq,'S') !initialize random number generator
      endif
      return
      end

!=======================================================================

      SUBROUTINE RMMARD( RVEC,LENV,ISEQ )

!-----------------------------------------------------------------------
!  C(ONE)X
!  R(ANDO)M (NUMBER GENERATOR OF) MAR(SAGLIA TYPE) D(OUBLE PRECISION)
!
!  THESE ROUTINES (RMMARD,RMMAQD) ARE MODIFIED VERSIONS OF ROUTINES
!  FROM THE CERN LIBRARIES. DESCRIPTION OF ALGORITHM SEE:
!               http://consult.cern.ch/shortwrups/v113/top.html
!  IT HAS BEEN CHECKED THAT RESULTS ARE BIT-IDENTICAL WITH CERN
!  DOUBLE PRECISION RANDOM NUMBER GENERATOR RMM48, DESCRIBED IN
!               http://consult.cern.ch/shortwrups/v116/top.html
!  ARGUMENTS:
!   RVEC   = DOUBLE PREC. VECTOR FIELD TO BE FILLED WITH RANDOM NUMBERS
!   LENV   = LENGTH OF VECTOR (# OF RANDNUMBERS TO BE GENERATED)
!   ISEQ   = # OF RANDOM SEQUENCE
!
!  VERSION OF D. HECK FOR DOUBLE PRECISION RANDOM NUMBERS.
!  ADAPTATION  : T. PIEROG    IK  FZK KARLSRUHE FROM D. HECK VERSION
!  DATE     : Feb  17, 2009
!-----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER          KSEQ
      PARAMETER        (KSEQ = 2)
      COMMON /CRRANMA3/CD,CINT,CM,TWOM24,TWOM48,MODCNS
      DOUBLE PRECISION CD,CINT,CM,TWOM24,TWOM48
      INTEGER          MODCNS
      COMMON /CRRANMA4/C,U,IJKL,I97,J97,NTOT,NTOT2,JSEQ
      DOUBLE PRECISION C(KSEQ),U(97,KSEQ),UNI
      INTEGER          IJKL(KSEQ),I97(KSEQ),J97(KSEQ),
     *                 NTOT(KSEQ),NTOT2(KSEQ),JSEQ

      DOUBLE PRECISION RVEC(*)
      INTEGER          ISEQ,IVEC,LENV
      SAVE

!-----------------------------------------------------------------------

      IF ( ISEQ .GT. 0  .AND.  ISEQ .LE. KSEQ ) JSEQ = ISEQ
      
      DO   IVEC = 1, LENV
        UNI = U(I97(JSEQ),JSEQ) - U(J97(JSEQ),JSEQ)
        IF ( UNI .LT. 0.D0 ) UNI = UNI + 1.D0
        U(I97(JSEQ),JSEQ) = UNI
        I97(JSEQ)  = I97(JSEQ) - 1
        IF ( I97(JSEQ) .EQ. 0 ) I97(JSEQ) = 97
        J97(JSEQ)  = J97(JSEQ) - 1
        IF ( J97(JSEQ) .EQ. 0 ) J97(JSEQ) = 97
        C(JSEQ)    = C(JSEQ) - CD
        IF ( C(JSEQ) .LT. 0.D0 ) C(JSEQ)  = C(JSEQ) + CM
        UNI        = UNI - C(JSEQ)
        IF ( UNI .LT. 0.D0 ) UNI = UNI + 1.D0
!  AN EXACT ZERO HERE IS VERY UNLIKELY, BUT LET'S BE SAFE.
        IF ( UNI .EQ. 0.D0 ) UNI = TWOM48
        RVEC(IVEC) = UNI
      ENDDO

      NTOT(JSEQ) = NTOT(JSEQ) + LENV
      IF ( NTOT(JSEQ) .GE. MODCNS )  THEN
        NTOT2(JSEQ) = NTOT2(JSEQ) + 1
        NTOT(JSEQ)  = NTOT(JSEQ) - MODCNS
      ENDIF

      RETURN
      END

!=======================================================================

      SUBROUTINE RMMAQD( ISEED, ISEQ, CHOPT )

!-----------------------------------------------------------------------
!  R(ANDO)M (NUMBER GENERATOR OF) MA(RSAGLIA TYPE INITIALIZATION) DOUBLE
!
!  SUBROUTINE FOR INITIALIZATION OF RMMARD
!  THESE ROUTINE RMMAQD IS A MODIFIED VERSION OF ROUTINE RMMAQ FROM
!  THE CERN LIBRARIES. DESCRIPTION OF ALGORITHM SEE:
!               http://consult.cern.ch/shortwrups/v113/top.html
!  FURTHER DETAILS SEE SUBR. RMMARD
!  ARGUMENTS:
!   ISEED  = SEED TO INITIALIZE A SEQUENCE (3 INTEGERS)
!   ISEQ   = # OF RANDOM SEQUENCE
!   CHOPT  = CHARACTER TO STEER INITIALIZE OPTIONS
!
!  CERN PROGLIB# V113    RMMAQ           .VERSION KERNFOR  1.0
!  ORIG. 01/03/89 FCA + FJ
!  ADAPTATION  : T. PIEROG    IK  FZK KARLSRUHE FROM D. HECK VERSION
!  DATE     : Feb  17, 2009
!-----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER          KSEQ
      PARAMETER        (KSEQ = 2)
      COMMON /CRRANMA3/CD,CINT,CM,TWOM24,TWOM48,MODCNS
      DOUBLE PRECISION CD,CINT,CM,TWOM24,TWOM48
      INTEGER          MODCNS
      COMMON /CRRANMA4/C,U,IJKL,I97,J97,NTOT,NTOT2,JSEQ
      DOUBLE PRECISION C(KSEQ),U(97,KSEQ),UNI
      INTEGER          IJKL(KSEQ),I97(KSEQ),J97(KSEQ),
     *                 NTOT(KSEQ),NTOT2(KSEQ),JSEQ

      DOUBLE PRECISION CC,S,T,UU(97)
      INTEGER          ISEED(3),I,IDUM,II,II97,IJ,IJ97,IORNDM,
     *                 ISEQ,J,JJ,K,KL,L,LOOP2,M,NITER
      CHARACTER        CHOPT*(*), CCHOPT*12
      LOGICAL          FIRST
      SAVE
      DATA             FIRST / .TRUE. /, IORNDM/11/, JSEQ/1/

      
!-----------------------------------------------------------------------

      IF ( FIRST ) THEN
        TWOM24 = 2.D0**(-24)
        TWOM48 = 2.D0**(-48)
        CD     = 7654321.D0*TWOM24
        CM     = 16777213.D0*TWOM24
        CINT   = 362436.D0*TWOM24
        MODCNS = 1000000000
        FIRST  = .FALSE.
        JSEQ   = 1
      ENDIF
      CCHOPT = CHOPT
      IF ( CCHOPT .EQ. ' ' ) THEN
        ISEED(1) = 54217137
        ISEED(2) = 0
        ISEED(3) = 0
        CCHOPT   = 'S'
        JSEQ     = 1
      ENDIF

      IF     ( INDEX(CCHOPT,'S') .NE. 0 ) THEN
        IF ( ISEQ .GT. 0  .AND.  ISEQ .LE. KSEQ ) JSEQ = ISEQ
        IF ( INDEX(CCHOPT,'V') .NE. 0 ) THEN
          READ(IORNDM,'(3Z8)') IJKL(JSEQ),NTOT(JSEQ),NTOT2(JSEQ)
          READ(IORNDM,'(2Z8,Z16)') I97(JSEQ),J97(JSEQ),C(JSEQ)
          READ(IORNDM,'(24(4Z16,/),Z16)') U
          IJ = IJKL(JSEQ)/30082
          KL = IJKL(JSEQ) - 30082 * IJ
          I  = MOD(IJ/177, 177) + 2
          J  = MOD(IJ, 177)     + 2
          K  = MOD(KL/169, 178) + 1
          L  = MOD(KL, 169)
          CD =  7654321.D0 * TWOM24
          CM = 16777213.D0 * TWOM24
        ELSE
          IJKL(JSEQ)  = ISEED(1)
          NTOT(JSEQ)  = ISEED(2)
          NTOT2(JSEQ) = ISEED(3)
          IJ = IJKL(JSEQ) / 30082
          KL = IJKL(JSEQ) - 30082*IJ
          I  = MOD(IJ/177, 177) + 2
          J  = MOD(IJ, 177)     + 2
          K  = MOD(KL/169, 178) + 1
          L  = MOD(KL, 169)
          DO   II = 1, 97
            S = 0.D0
            T = 0.5D0
            DO   JJ = 1, 48
              M = MOD(MOD(I*J,179)*K, 179)
              I = J
              J = K
              K = M
              L = MOD(53*L+1, 169)
              IF ( MOD(L*M,64) .GE. 32 ) S = S + T
              T = 0.5D0 * T
            ENDDO
            UU(II) = S
          ENDDO
          CC    = CINT
          II97  = 97
          IJ97  = 33
!  COMPLETE INITIALIZATION BY SKIPPING (NTOT2*MODCNS+NTOT) RANDOMNUMBERS
          NITER = MODCNS
          DO   LOOP2 = 1, NTOT2(JSEQ)+1
            IF ( LOOP2 .GT. NTOT2(JSEQ) ) NITER = NTOT(JSEQ)
            DO   IDUM = 1, NITER
              UNI = UU(II97) - UU(IJ97)
              IF ( UNI .LT. 0.D0 ) UNI = UNI + 1.D0
              UU(II97) = UNI
              II97     = II97 - 1
              IF ( II97 .EQ. 0 ) II97 = 97
              IJ97     = IJ97 - 1
              IF ( IJ97 .EQ. 0 ) IJ97 = 97
              CC       = CC - CD
              IF ( CC .LT. 0.D0 ) CC  = CC + CM
            ENDDO
          ENDDO
          I97(JSEQ) = II97
          J97(JSEQ) = IJ97
          C(JSEQ)   = CC
          DO   JJ = 1, 97
            U(JJ,JSEQ) = UU(JJ)
          ENDDO
        ENDIF
      ELSEIF ( INDEX(CCHOPT,'R') .NE. 0 ) THEN
        IF ( ISEQ .GT. 0 ) THEN
          JSEQ = ISEQ
        ELSE
          ISEQ = JSEQ
        ENDIF
        IF ( INDEX(CCHOPT,'V') .NE. 0 ) THEN
          WRITE(IORNDM,'(3Z8)') IJKL(JSEQ),NTOT(JSEQ),NTOT2(JSEQ)
          WRITE(IORNDM,'(2Z8,Z16)') I97(JSEQ),J97(JSEQ),C(JSEQ)
          WRITE(IORNDM,'(24(4Z16,/),Z16)') U
        ELSE
          ISEED(1) = IJKL(JSEQ)
          ISEED(2) = NTOT(JSEQ)
          ISEED(3) = NTOT2(JSEQ)
        ENDIF
      ENDIF

      RETURN
      END