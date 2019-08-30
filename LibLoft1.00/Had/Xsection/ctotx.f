!     ctotx: x A collsion total cross section 
!      
      module modcosXsec
      logical,save:: nIsDiff=.false.
      end module modcosXsec
      subroutine ctotx(pjin, Ain, xs)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmass.h"
#include "Zevhnp.h"
      type(ptcl)::pjin   !input projectile particle
      real(8),intent(in)::  Ain          ! input. effective target mass number
      real(8),intent(out):: xs          ! output. in mb. total xsection. 
!           ctotx0 gives small Xs so we renormalize by
!           using cPDGsigmaTotpA(A) at @200 GeV 
      real(8),intent(in):: Zin
      real(8):: cPDGsigmaTotpA

      type(ptcl):: pj
      type(ptcl):: proton
      save proton
      real(8),save:: xst
      logical,save::first=.true. 
      real(8),save::Asave=-1.
      real(8),save::ratio
      real(8):: Z, A
      real(8),external:: cA2Z

      Z = cA2Z(Ain)
      goto 20


      entry ctotx2(pjin, Ain, Zin, xs)


      Z = Zin
 20   continue
      pj = pjin
      A = Ain      
      if( pjin%code == kgnuc .and. Ain == 1.0) then 
         call cmkptc(knuc, -1, int(Z), pj)
         pj%fm%p(4) = pjin%fm%p(4)/pjin%subcode
         A = pjin%subcode
         Z = pjin%charge
      endif


      call ctotx0(pj, A, Z,  xs)
      if( A /= 1.0d0 ) then 
         if(first) then
            call cmkptc(knuc, -1, 1, proton)
            proton%fm%p(4) =  200. + masp
            first = .false.
         endif
         if( A /= Asave ) then
            call ctotx0(proton, A, Z,  xst)
            ratio = cPDGsigmaTotpA(A)/xst
            Asave =  A
         endif
         xs = xs *ratio
      endif
      end

      subroutine ctotx0(pj, A, Z, xs)
      use modcosXsec
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zevhnp.h"
      type(ptcl)::pj   !input projectile particle
      real(8),intent(in):: A  ! effective target mass number
      real(8),intent(in):: Z  ! effective target charge
      real(8),intent(out)::xs  ! in mb. total xsection. 

      real*8 p

      real(8)::shp, shn
      if( pj%fm%p(4) .le. pj%mass) then
         if(pj%code .eq. knuc .and. pj%subcode .eq. antip) then
            xs = largexs
         elseif(pj%code .eq. kelec .and. pj%charge .eq. 1) then
            xs = largexs
         else
            xs  = smallxs
         endif
      else   
         p = sqrt(pj%fm%p(4)**2 - pj%mass**2)
!         p = max(p, 0.1d0) 
!         if( p .ge. 20.) then          
!            call cerrorMsg('Momentum is >20 GeV/c for ctotx',0)
!         endif

         if(pj%code .eq. knuc) then
            if(pj%charge .eq. 1) then
!            proton
!               shp = ctotpp1(p)
               if( A == 1.0 ) then
                  if( Z == 1.0) then  ! pp
                     call cppTotXs(p, xs)
                  elseif( Z== 0.0 ) then  ! pn
                     call cnpTotXs(p, xs)
                  else
                     write(0,*) 'pj=proton target A,Z=',A,Z
                     write(0,*) ' stragne for ctotx/ctotx2'
                     stop
                  endif
               elseif(A > 1.) then
                  call cppTotXs(p, shp)
                  if( nIsDiff ) then
                     call cnpTotXs(p, shn)
!               use average of pp,pn
                     shp = (shp + shn)/2.0
                  endif
                  call cxp2xAXsec(A, shp, xs)
               else
                  write(0,*) 'pj: p, target A,Z=',A,Z
                  write(0,*) 'strange for ctotx/ctotx2'
                  stop
               endif
            elseif(pj%charge .eq. -1) then
!                pbar  no n target
               call cpbarpTotXs(p, shp)
               if( A .gt. 1.) then
                  call cxp2xAXsec(A, shp, xs)
               else
                  xs = shp
               endif
            elseif(pj%subcode .eq. antip) then
!              anti-neutron; assume the same one as pbar
!               shp = ctotpbarp1(p)
               call cnbarpTotXs(p, shp)
               if( A .gt. 1.) then
                  call cxp2xAXsec(A, shp, xs)
               else
                  xs = shp
               endif
            else
!               neutron
!               shp = ctotnp1(p)
               call cnpTotXs(p, shp)
               if(A .gt. 1.) then
                  call cxp2xAXsec(A, shp, xs)
               else
                  xs = shp
               endif
            endif
         elseif(pj%code .eq. kpion) then
            if(pj%charge .eq. 1) then
               call cpippTotXs(p, shp)
            elseif(pj%charge .eq. -1) then
               call cpimpTotXs(p, shp) 
            else
!              at low energy, pi0 would not interact. any would be o.k
!               shp = ctotpiMp1(p)
               call cpimpTotXs(p, shp)
            endif
            if(A .ne. 1.0) then
               call cxp2xAXsec(A, shp, xs)
            else
               xs = shp
            endif
         elseif(pj%code .eq. kkaon) then
            if(pj%charge .eq. 1) then
               if( A == 1.0 ) then
                  if( Z == 1.0 ) then
                     call ckppTotXs(p, xs)
                  elseif( Z == 0. ) then
                     call ckpnTotXs(p, xs)
                  else
                     write(0,*) ' A,Z=',A,Z, 'for K+ pj'
                     write(0,*) ' strange to ctotx/ctotx2  '
                     stop
                  endif
               elseif(A > 1.) then
                  call ckppTotXs(p, shp)
                  if( nIsDiff ) then
                     call ckpnTotXs(p, shn)
                     shp = (shp + shn)/2.0
                  endif
                  call cxp2xAXsec(A, shp, xs)               
               else
                  write(0,*) ' A,Z=',A,Z, 'for K+ p'
                  write(0,*) ' strange to ctotx/ctotx2  '
                  stop
               endif
            elseif(pj%charge .eq. -1) then  ! K-
               call ckmpTotXs(p, shp)   ! kmp = kmn at present
               if(A .gt. 1.) then
                  if( nIsDiff ) then
                     call ckmnTotXs(p, shn)
                     shp = (shp + shn)/2.
                  endif
                  call cxp2xAXsec(A, shp, xs)               
               else
                  xs = shp
               endif
            else
!             k0; don't worry so much 
               call ckmpTotXs(p, shp)
               call ckmnTotXs(p, shn)
               shp = (shp+ shn)/2.
               if(A .gt. 1.) then
                  call cxp2xAXsec(A, shp, xs)               
               else
                  xs = shp
               endif
            endif
         elseif(pj%code .eq. kgzai .or. pj%code .eq. ksigma .or.
     *           pj%code .eq. kbomega .or. pj%code .eq. klambda .or.
     *           pj%code .eq. klambdac) then
!          don't worry, not used almost at all; use proton
!            shp = ctotpp1(p)
            call cppTotXs(p, shp)
            call cxp2xAXsec(A, shp, xs)               
         elseif( pj%code == kgnuc ) then
            write(0,*) 'Sorry: ctotx is not usable for heavy ions'
            write(0,*) ' only inelastic xs is used for heavy ions'
            write(0,*) ' so that you may use cinela or cAAXsec2'
            stop
         else
!          use pion
!            shp =ctotpiPp1(p)
            call cpippTotXs(p, shp)
            call cxp2xAXsec(A, shp, xs)               
         endif
      endif
      end

      function cA2Z(A) result(Z)
      implicit none
!      very rough charge assignment for mass # A nucleus
      real(8),intent(in):: A  ! mass #

      real(8)::Z
      Z = int( A/2.15 + 0.7 )
      end      function cA2Z


