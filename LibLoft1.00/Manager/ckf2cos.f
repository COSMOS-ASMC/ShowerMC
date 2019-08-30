      subroutine ckf2cos(kf, code, subcode, chg)
!        kf code to cosmos code.
      implicit none
#include "Zcode.h"
#include "Zkfcode.h"

      integer kf   ! input
      integer code, subcode, chg  ! output. For spectator kf, code=krare
!                                           for tau neutrio,  code =krare
      character*80 msg
      real(8)::u
      integer kfabs
      kfabs = abs(kf)
!          special treatment for K0
      if( kfabs == kfk0 ) then  ! K0 : 50% K0s; 50 % K0L
         call rndc(u)
         if(u < 0.5 ) then
            kfabs = kfk0s
         else
            kfabs = kfk0l
         endif
      endif
!           
      if(kfabs .eq. kfpion) then
         code = kpion
         subcode = sign(1, -kf)
         chg = sign(1, kf)
      elseif(kfabs .eq. kfpi0) then
         code = kpion
         subcode = 0
         chg = 0
      elseif( kfabs .eq. kfkaon ) then
         code = kkaon
         subcode =sign(1, -kf)
         chg = sign(1, kf)
      elseif(kfabs .eq.  kfk0l) then
         code =kkaon
         subcode =sign( k0l, kf)
         chg = 0
      elseif(kfabs .eq.  kfk0s) then
         code = kkaon
         subcode =sign( k0s, kf)
         chg = 0
      elseif(kfabs .eq. kfneutron) then
         code =knuc
         if(kf .gt. 0) then
            subcode = kneutron
         else
            subcode = kneutronb
         endif
         chg = 0
      elseif(kfabs .eq. kfproton) then
         code = knuc
         if(kf .gt. 0) then
            subcode = regptcl
         else
            subcode = antip
         endif
         chg =sign(1, kf)
!cc      elseif(kfabs .ge. 10000) then
!cc         code = krare   ! target spectator. neglect
!         p(i,5) has mass. kfabs-10000=Z
!  &&&&&&&&&&&&&&&&&&
!         write(msg, *) 'kf code=',kf, ' not treatable'
!         call cerrorMsg(msg, 1)
!         call cerrorMsg('the particle is neglected',1)
!   &&&&&&&&&&&&&&&&
      elseif(kfabs  .eq. kfeta) then
         code= keta
         subcode = 0
         chg = 0
      elseif(kfabs .eq. kfelec) then
         code = kelec
         subcode = sign(1, -kf)
         chg = sign(1, -kf)
      elseif( kf .eq.  kfphoton ) then
         code = kphoton
         subcode = 0
         chg = 0
      elseif(kfabs .eq. kfmuon) then
         code = kmuon
         subcode = sign(1, -kf)
         chg = sign(1, -kf)
      elseif(kfabs .eq. kfneue) then
         code = kneue
         subcode =sign(1, -kf)
         chg = 0
      elseif(kfabs .eq. kfneumu) then
        code = kneumu
        subcode =sign(1, -kf)
        chg = 0
      elseif(kfabs .eq. kfdmes) then
         code = kdmes
         subcode = sign(1, -kf)
         chg = sign(1, kf)
      elseif(kfabs .eq. kfd0) then
         code =kdmes
         subcode = sign(1, -kf)
         chg = 0
      elseif(kfabs .eq. kflambda) then
         code = klambda
         subcode = sign(1, -kf)
         chg = 0
      elseif( kfabs .ge. 1000000020 .and. kfabs .le. 1000922350 ) then ! 2n system
!           nucleus is  10LZZZAAAI where L=0 (if not containing strangeness)
!           ZZZ is 3 digits of charge, AAA 3digits of mass no. I=0 for
!           gound state; so for our case, 1000010020 is minimum for deuteron
!           and 1000922350 is max for U(235)
         if(kf .lt. 0) then
            write(0,*) ' anti nuclues is not yet supported '
            code = krare
         else
            code=kgnuc
            subcode =( (kfabs/10)*10-(kfabs/10000)*10000 )/10
            chg =( (kfabs/10000)*10000 - (kfabs/1000000)*1000000 )
     *        /10000     
            if( kfabs < 1000010020 ) then  ! charge 0               
               code = krare
            endif
         endif
      elseif(kfabs == kfDelta0) then
         code = kDelta
         subcode = -1
         chg =0
      elseif(kfabs == kfDeltap) then
         code = kDelta
         subcode = -1
         chg = 1
      elseif(kfabs == kfDeltam) then
         code = kDelta
         subcode = 1
         chg = -1
      elseif(kfabs .eq. kfsigma0) then
         code = ksigma
         subcode = sign(1, -kf)
         chg = 0
      elseif(kfabs .eq. kfsigmap) then
         code = ksigma
         subcode = sign(1, -kf)
         chg = sign(1, kf)
      elseif(kfabs .eq. kfsigmam) then
         code = ksigma
         subcode = sign(1, -kf)
         chg = sign(1, -kf)
      elseif(kfabs .eq. kfgzai0 ) then
         code = kgzai
         subcode = sign(1, -kf)
         chg = 0
      elseif(kfabs .eq. kfgzai ) then
         code = kgzai
         subcode = sign(1, -kf)
         chg = sign(1, -kf)
      elseif(kfabs .eq. kflambdac) then
         code = klambdac
         subcode = sign(1, -kf)
         chg = sign(1, kf)
      elseif(kfabs .eq. kfbomega ) then
         code = kbomega
         subcode = sign(1, -kf)
         chg =  sign(1, -kf)
      elseif(kfabs .eq. kftau ) then
         code = ktau            
         subcode = sign(1, -kf)
         chg = sign(1, -kf)
      elseif( kfabs == kfneutau )then
         code = kfneutau
         subcode = sign(1,-kf)
         chg= 0
      elseif(kfabs .eq. kfrho) then  ! rho
         code = krho
         subcode = 0
         chg = 0
      elseif(kfabs == kfrhoc) then
         code = krho
         subcode = sign(1, -kf)
         chg = sign(1, kf)
      elseif(kfabs .eq. kfomega) then   ! omega
         code = komega
         subcode = 0
         chg = 0
      elseif( kfabs .eq. kfphi ) then  ! phi
         code = kphi
         subcode = 0
         chg = 0
      elseif( kfabs == kfds ) then ! Ds
         code = kds
         chg = sign(1, kf )
         subcode = sign(1,-kf)
      elseif( kfabs == kfXic ) then
         code = kXic
         chg = sign(1,kf)
         subcode = sign(1,-kf)
      elseif( kfabs == kfXic0 ) then
         code = kXic0
         chg =0
         subcode =sign(1,-kf)
      elseif( kfabs == kfomeC0 ) then
         code = komeC0
         chg = 0
         subcode =  sign(1,-kf)
      elseif( kfabs == kfetap ) then
         code = ketap
         chg = 0
         subcode = 0
      else
         code = krare
      endif
      if( code == krare ) then
         write(msg, *) 'not implemented kf code=', kf
         call cerrorMsg(msg, 1)
         call cerrorMsg('we neglect this particle',1)
      endif
      end

      subroutine ccos2kf(code, subcode, chg, kf)
!         cosmos code to kf code; 
      implicit none
#include "Zcode.h"
#include "Zkfcode.h"
      integer,intent(in):: code, subcode, chg ! cosmos code but not integer*2
      integer,intent(out):: kf  !  0 means code is not in PDG
!                              (say, kEdepo ...)          
      character*80 msg      
      if(code .eq. kelec) then
         kf = sign(kfelec, -chg)
      elseif(code .eq. kphoton) then
         kf = kfphoton
      elseif(code .eq. kpion) then
         if(chg .eq. 0) then
            kf = kfpi0
         else
            kf = sign(kfpion, chg)
         endif
      elseif(code .eq. kkaon) then
         if(chg .eq. 0) then
            if(abs(subcode) .eq. k0l) then
               kf = sign( kfk0l, subcode)
            else
               kf =sign( kfk0s,  subcode)
            endif
         else
            kf = sign(kfkaon, chg)
         endif
      elseif(code .eq. knuc) then
         if(chg .eq. 0) then
            kf = sign( kfneutron, -subcode)
         else
            kf = sign(kfproton, chg)
         endif
      elseif(code .eq. kmuon) then
         kf = sign(kfmuon, -chg)
      elseif(code .eq. kneue ) then
         kf = sign(kfneue, -subcode)
      elseif(code .eq. kneumu) then
         kf = sign(kfneumu, -subcode)
      elseif(code .eq. kdmes) then
         if(chg .eq. 0) then
            kf = sign(kfd0, -subcode)
         else
            kf = sign(kfdmes, -chg)
         endif
      elseif(code .eq. klambda) then
         kf = sign(kflambda, -subcode)
      elseif(code .eq. kgnuc ) then
         kf = 1000000000 + abs(chg)*10000 + subcode*10
         kf = sign(kf, chg)
      elseif( code == keta ) then
         kf = kfeta
      elseif( code == krho ) then
         if( chg == 0) then
            kf = kfrho
         else
            kf = sign(kfrhoc, chg)
         endif
      elseif( code == kDelta ) then
         if( chg == 0 ) then
            kf = kfDelta0
         elseif( chg == 1 ) then
            kf = kfDeltap
         else
            kf = kfDeltam
         endif
      elseif ( code == komega ) then
         kf = kfomega
      elseif( code == kphi ) then
         kf = kfphi
      elseif(code .eq. ksigma) then
         if(chg .eq. 1) then
            kf = sign(kfsigmap, -subcode)
         elseif(chg .eq. 0) then
            kf = sign(kfsigma0, -subcode)
         else
            kf = sign(kfsigmam, -subcode)
         endif
      elseif(code .eq. kgzai) then
         if(chg .eq. 0) then
            kf = sign(kfgzai0, -subcode)
         else
            kf = sign(kfgzai, -chg)
         endif
      elseif(code .eq. klambdac) then
         kf =sign(kflambdac, chg)
      elseif(code .eq. kbomega ) then
         kf = sign(kfbomega, -chg)
      elseif( code == kds ) then
         kf = sign(kfds, -chg)
      elseif( code ==  kXic ) then
         kf = sign(kfXic, -chg)
      elseif( code == komeC0 ) then
         kf = sign(kfomeC0, -subcode)
      elseif( code == ktau ) then
         kf = sign(kftau, -subcode)
      elseif( code == kneutau) then
         kf = sign(kfneutau, -subcode)
      elseif( code == ketap ) then
         kf = kfetap
      elseif( code <= 0 )then
         kf = 0
      else
         write(msg, *) 'code, subcode, chg to  ccos2kf=', 
     &      code, subcode, chg  
         call cerrorMsg(msg, 1)
         call cerrorMsg(' cannot be converted to kf code',1)
         kf=0
      endif


      end
      
      subroutine ckf2cosB(kf, p)
      implicit none
#include "Zptcl.h"
      integer,intent(in):: kf ! pdg ptcl code
      type(ptcl):: p   ! output p%code p%sucode and p%charge are set

      integer:: code, subcode, charge
      call ckf2cos(kf, code, subcode, charge)
      p%code = code
      p%subcode = subcode
      p%charge = charge
      end

      subroutine ccos2kfB(p, kf)  ! inverse of ckf2cosB
      implicit none
#include "Zptcl.h"
      type(ptcl),intent(in):: p  ! input; p%code, p%subcode and p%charge 
      integer,intent(out):: kf  ! pdg code for p%code...

      integer:: code, subcode, charge

      code = p%code
      subcode = p%subcode
      charge = p%charge
      call ccos2kf(code, subcode, charge, kf)
      end

       




