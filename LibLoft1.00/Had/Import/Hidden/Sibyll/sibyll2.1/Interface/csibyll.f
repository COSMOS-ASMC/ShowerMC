      subroutine csibyllinit
      implicit none
!
      COMMON /S_CLDIF/ LDIFF
      INTEGER          LDIFF
      COMMON /S_CSYDEC/IDB(49), CBR(102), KDEC(612), LBARP(49)
      REAL             CBR
      INTEGER          IDB,KDEC,LBARP
      integer  KODFRAG                  
      COMMON /CKFRAG/ KODFRAG                  

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN

      

      KODFRAG =0  ! complex fragmentation by Abrasion-ablation model
      LDIFF = 0       ! mix diff and non-diff

      LUN = 0         ! logical dev. # for (error) message  from sibyll
! sibyll default is 7.
      call SIBYLL_INI
      call SIGMA_INI   ! p-/pi-/K-Air  sigma and mfp
                   ! 
      call NUC_NUC_INI  ! A-A' xsec.  init.  seems false.
            ! NUC_NUC_INI calls sigma_ini agian; is it ok?
!      call INI_WRITE(0) !!!!! print pAir xsec.
!        don't let mu, pi, K decay;
!           for AA' col.
!       call SIGNUC_INI(...) !  this is almost useless
!                    so we make cgetAAXsSib for each
!     AA' col.
      IDB(4:12) = -IDB(4:12)

!            short life particle decay is managed by us
!
!       same for eta
      IDB(23)  = -IDB(23)
!       same for  SIGMA, XI, LAMBDA
!      do  i = 34, 39
      IDB(34:36) = -IDB(34:36)
      IDB(39) = -IDB(39)
      
!     decay control: can force not to decay for
!        some of short life ptcls s
!!!      call csibSetStblPtcl
      end
!     ****************
      subroutine  csibyllevent(pj, ia, iz, a, ntp)
!      use modXsecMedia
      implicit none
#include "Zmanager.h"      
#include "Zcode.h"
#include "Zptcl.h"
#include "Zmass.h"
#include "ZmediaLoft.h"      
!
      type(ptcl):: pj  ! input projectile particle
      integer,intent(in):: ia, iz  ! target A,Z (Z will not be
                  ! used)
      type(ptcl):: a(*)
      integer,intent(out)::ntp  ! # of ptcls produced in a.

      integer k, i, L
      integer::  icp, iat
      integer::  j
      real*8::  E0, Ecm, g, beta, Eg
      integer::icon

      real  sqs
!            not used
!      type(ptcl)::pjcms, tgcms, tgrest
!      data tgrest%fm%p(1:4)/0.d0,0.d0,0.d0,masp/
!      data tgrest.mass/masp/
!      type(ptcl)::Cmsp
!      type(ptcl)::pjpnlab  ! proj. /n in lab 
!      data pjpnlab%fm%p(1:2)/0.d0, 0.d0/


!       call csibyllIniEvent   ! no need at presnt
      call ccoscode2sibyll(pj, icp)  ! get sibyll code of pj
      iat =  ia                 !  target mass number
      if( CosOrEpi /= 'gencol') then
         if( Media(MediaNo)%name == "Air") then
            iat = 0
         endif
      endif   
      if(pj%code .eq. kgnuc) then
         E0 = pj%fm%p(4)/pj%subcode !  GeV/n in Lab
         Ecm = sqrt( ( E0 +masp )*2*masp ) ! nn system
      else
!         if( pj%code == klambda ) then   ! upto v8.02
!            icp = 14            ! n
!         endif
         E0 = pj%fm%p(4)
         Ecm = sqrt( E0*2*masp + pj%mass**2 + masp**2)
      endif

      sqs = Ecm
      g = (E0 + masp)/Ecm  ! cms gamma factor  !=sqs/2/mp for pi/K
      beta = sqrt( (g-1.d0) * (g + 1.d0)) /g
      if( pj%code /= kgnuc ) then
         call SIBYLL(icp, iat, sqs)  ! event generator
              ! icp = 7,8, 9,10, 11,12, 13,14, -13,-14       
              !   pi+-,K+-,  KL,KS, p,n,  pbar,nbar
              !  The output is contained in COMMON /S_PLIST/
              ! iat= 0--> Air . sqs = root(s) in GeV for pj-nucleon system.
         call DECSIB            !  all unfamilier ptcls should decay
      else
         call SIBNUC( icp, iat, sqs)  !A-A' col. icp  sqs is n-n system.
                            ! icp=A, iat= A', iat= 0 is for Air
      endif

      if( pj%code /= kgnuc ) then
         call csibhA2coscode(g, beta, a, ntp)
      else
         call csibAA2coscode(g, beta, a, ntp)
      endif
      call crot3mom(pj, a, ntp)
      end

      subroutine csibyllIniEvent
!        at present nothing to do 
      implicit none
      integer:: NP, LLIST
      real::   P
      COMMON /S_PLIST/ NP, P(8000,5), LLIST(8000)
      REAL             PA
      INTEGER          LLA,NPA
      COMMON /S_PLNUC/ NPA, PA(5,40000), LLA(40000)


      end      subroutine csibyllIniEvent


      subroutine csibhA2coscode(g, beta, a, ntp)
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
#include "Zmass.h"
!
      real(8),intent(in):: g  ! cms g factor
      real(8),intent(in):: beta  ! its beta 
      type(ptcl)::a(*)  ! output patlcs
      integer,intent(out):: ntp  ! # of ptcls put in a

      integer:: NP, LLIST
      real::   P
      COMMON /S_PLIST/ NP, P(8000,5), LLIST(8000)

      integer:: j, L
      
      ntp = 0
      do  j = 1, NP
         if(abs(LLIST(j)) .lt. 10000) then
            ntp = ntp + 1   
            a(ntp)%fm%p(1:2)=P(j, 1:2)
            a(ntp)%fm%p(3) = g* (P(j,3) +
     *              beta*P(j,4))
            a(ntp)%fm%p(4) = g* (P(j,4) +
     *              beta*P(j,3))
            L = mod(LLIST(j), 10000)
!            L = LLIST(j)
            call csibyllcode2cos(L, a(ntp))
         endif
      enddo
      end      subroutine csibhA2coscode
      subroutine csibAA2coscode(g, beta, a, ntp)
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
#include "Zmass.h"
!     
      real(8),intent(in):: g  ! cms gamma factor
      real(8),intent(in):: beta ! its beta
      type(ptcl)::a(*)  ! output. ptcls stored here
      integer,intent(out):: ntp  ! # of ptlcs in a

      REAL             PA
      INTEGER          LLA,NPA
      COMMON /S_PLNUC/ NPA, PA(5,40000), LLA(40000)

      integer:: L
      integer:: j
      
      ntp = 0
      do  j = 1, NPA
         if(abs(LLA(j)) .lt. 10000) then
            ntp = ntp + 1   
            a(ntp)%fm%p(1:2)=PA(1:2,j)
            a(ntp)%fm%p(3) = g* (PA(3,j) + beta*PA(4,j))
            a(ntp)%fm%p(4) = g* (PA(4,j) + beta*PA(3,j))
            L = mod(LLA(j), 10000)
            call csibyllcode2cos(L, a(ntp))
         endif
      enddo
      end      subroutine csibAA2coscode

       real  function  S_RNDM(X)
       integer  X  !  not used
       real*8 u
       call rndc(u)
       S_RNDM = u
       end
       subroutine ccoscode2sibyll(pj, ksib)
       implicit none
#include "Zptcl.h"
#include "Zcode.h"
      type(ptcl):: pj
      integer ksib
!
!      projectiel 7,8,  9,10, 11,12, 13,14, -13,-14
!                 pi+-,  K+-,  KL,KS, p,n,  pbar,nbar
      
      if( pj%code .eq. knuc ) then
         if( pj%charge .eq. 1) then
            ksib = 13
         elseif(pj%charge .eq. -1) then
            ksib = -13
         elseif( pj%charge .eq. 0) then
            if( pj%subcode .eq. antip ) then
               ksib = -14
            else
               ksib = 14
            endif
         endif
      elseif( pj%code .eq. kpion ) then
         if( pj%charge .eq. 1 ) then
            ksib = 7
         elseif( pj%charge .eq. -1 ) then
            ksib = 8
         else
            ksib = 6
         endif
      elseif( pj%code .eq. kkaon ) then
!   9 K+    10 K-  11  K0L   12 K0s
         if( pj%charge .eq. 1 ) then
            ksib= 9
         elseif( pj%charge .eq. -1 ) then
            ksib = 10
         elseif( pj%subcode .eq. k0s ) then
            ksib = 12
         else
            ksib = 11
         endif
      elseif( pj%code == kgnuc) then
         ksib = pj%subcode   ! if 0, assumed to be  Air
      elseif( pj%code == klambda) then
         ksib = 39
      else
         write(0,*) 
     *    ' code =',pj%code, ' not acceptable in sibyll'
         stop 12345
      endif
      end
      subroutine csibyllcode2cos(ksibin, pj)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
      integer,intent(in):: ksibin  ! sibyll code
      type(ptcl):: pj    !   output cos code is put here
!
      real*8 u 
      integer ksib, code, subcode, charge

      ksib = abs(ksibin)
      subcode = regptcl
      if(ksib .eq. 1) then
         code = kphoton
         charge = 0
      elseif(ksib .eq. 2) then
         code = kelec
         charge = 1
      elseif(ksib .eq. 3) then
         code = kelec
         charge = -1
      elseif(ksib .eq. 4) then
         code = kmuon
         charge = 1
      elseif(ksib .eq. 5) then
         code = kmuon
         charge = -1
      elseif(ksib .eq. 6) then
         code = kpion
         charge = 0
      elseif(ksib .eq. 7) then
         code = kpion
         charge = 1
      elseif(ksib .eq. 8) then
         code = kpion
         charge = -1
      elseif(ksib .eq. 9) then
         code = kkaon
         charge = 1
      elseif(ksib .eq. 10) then
         code = kkaon
         charge = -1
      elseif(ksib .eq. 11) then
         code = kkaon
         charge = 0
         subcode = k0l
      elseif(ksib .eq. 12) then
         code = kkaon
         charge = 0
         subcode = k0s
      elseif(ksib .eq. 13) then
         code = knuc
         subcode = regptcl
         charge = 1
      elseif(ksib .eq. 14) then
         code = knuc
         charge = 0
         subcode = regptcl
      elseif(ksib .eq. 15) then
         code = kneue
         subcode = regptcl
         charge = 0
      elseif(ksib .eq. 16) then
         code = kneue
         subcode = antip
         charge = 0
      elseif(ksib .eq. 17) then
         code = kneumu
         subcode = regptcl
         charge = 0
      elseif(ksib .eq. 18) then
         code = kneumu
         charge = 0
         subcode = antip
      elseif(ksib .eq. 19) then
         code = knuc
         charge = -1
         subcode = antip
      elseif(ksib .eq. 20) then
         code = knuc
         charge = 0
         subcode = antip
      elseif(ksib .eq. 21) then
         code = kkaon
         call rndc(u)
         charge = 0
         if(u .lt. 0.5) then 
            subcode = k0s
         else
            subcode = k0l
         endif
      elseif(ksib .eq. 22) then
         code = kkaon
         charge = 0
         if(u .lt. 0.5) then 
            subcode = k0s
         else
            subcode = k0l
         endif
      elseif(ksib .eq. 23) then
         code = keta
         charge = 0
      elseif(ksib .eq. 25) then
         code = krho
         charge = 1
      elseif(ksib .eq. 26) then
         code = krho
         charge = -1
      elseif(ksib .eq. 27) then
         code = krho
         charge = 0
      elseif(ksib .eq. 32) then
         code = komega
         charge = 0
         subcode = regptcl
      elseif(ksib .eq. 33) then
         code = kphi
         subcode = regptcl
         charge = 0
      elseif(ksib .eq. 34) then
         code = ksigma
         subcode = regptcl
         charge = 1
      elseif(ksib .eq. 35) then
         code = ksigma
         subcode = regptcl
         charge = 0
      elseif(ksib .eq. 36) then
         code = ksigma
         subcode = regptcl
         charge = -1
      elseif(ksib .eq. 37) then
         code = kgzai
         subcode = regptcl
         charge = 0
      elseif(ksib .eq. 38) then
         code = kgzai
         subcode = regptcl
         charge = -1
      elseif(ksib .eq. 39) then
         code = klambda
         subcode = regptcl
         charge = 0
      elseif( ksib > 1000) then
         subcode = ksib- 1000
         if( subcode > 1 ) then
            code = kgnuc
            charge = subcode/2.15 +0.7 ! some plausible Z
         else
            code  = knuc
            call rndc(u)
            if(u < 0.6) then
               charge = 0
            else
               charge = 1
            endif
            subcode = -1
         endif
      else
         write(0,*) ' ****************sibyllcode=',ksib, ksibin
         code = krare
         charge = 0
         subcode = regptcl
      endif
      if(ksibin .lt. 0) then
         subcode = antip
         charge = - charge
      endif
      call cmkptc(code, subcode, charge, pj)
      end
      subroutine csibylGetDiffCode(nwout, difcode)
      implicit none
!    in case of nuclei for each interaction]. The meaning is
!.       JDIF(JW)  = diffraction code    !!!! changed to field !!!!
!.                  (0 : non-diffractive interaction) 
!.                  (1 : forward diffraction)     
!.                  (2 : backward diffraction) 
!.                  (3 : double diffraction)  
      integer,intent(out):: nwout !
      !   NSD is 0 or 3
      integer,intent(out):: difcode(20)
!          for p/n projectile,nw =1, difcode(1) shows the                      
!          diffraction state.                                                  
      integer NW_max, NS_max, NH_max, NJ_max
      PARAMETER (NW_max = 20)
      PARAMETER (NS_max = 20, NH_max = 50)
      PARAMETER (NJ_max = (NS_max+NH_max)*NW_max)
      real(4):: X1J, X2J, X1JSUM, X2JSUM, PTJET, PHIJET
      integer:: NNPJET, NNPSTR, NNSOF, NNJET, JDIF, NW, NJET, NSOF

      COMMON /S_CHIST/ X1J(NJ_max),X2J(NJ_max),
     &    X1JSUM(NW_max),X2JSUM(NW_max),PTJET(NJ_max),PHIJET(NJ_max),
     &    NNPJET(NJ_max),NNPSTR(2*NW_max),NNSOF(NW_max),NNJET(NW_max),
     &    JDIF(NW_max),NW,NJET,NSOF

      nwout = NW
      difcode(1:nwout)=JDIF(1:nwout)
      end

