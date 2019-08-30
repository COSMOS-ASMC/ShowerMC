!     ctotx: x A collsion total cross section 
!      
      subroutine ctotx(pj, A, xs)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmass.h"
#include "Zevhnp.h"
      type(ptcl)::pj   !input projectile particle
      real(8),intent(in)::  A          ! input. effective target mass number
      real(8),intent(out):: xs          ! output. in mb. total xsection. 
!           ctotx0 gives small Xs so we renormalize by
!           using cPDGsigmaTotpA(A) at @200 GeV 
      real(8):: cPDGsigmaTotpA
      type(ptcl):: proton
      save proton
      real(8),save:: xst
      logical,save::first=.true. 
      real(8),save::Asave=-1.
      real(8),save::ratio


      call ctotx0(pj, A, xs)
      if( A /= 1.0d0 ) then 
         if(first) then
            call cmkptc(knuc, -1, 1, proton)
            proton.fm.p(4) =  200. + masp
            first = .false.
         endif
         if( A /= Asave ) then
            call ctotx0(proton, A, xst)
            ratio = cPDGsigmaTotpA(A)/xst
            Asave =  A
         endif
         xs = xs *ratio
      endif
      end

      subroutine ctotx0(pj, A, xs)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zevhnp.h"
      type(ptcl)::pj   !input projectile particle
      real*8  A          ! input. effective target mass number
      real*8 xs          ! output. in mb. total xsection. 

      real*8 p

      real(8)::shp, shn
      if( pj.fm.p(4) .le. pj.mass) then
         if(pj.code .eq. knuc .and. pj.subcode .eq. antip) then
            xs = largexs
         elseif(pj.code .eq. kelec .and. pj.charge .eq. 1) then
            xs = largexs
         else
            xs  = smallxs
         endif
      else   
         p = sqrt(pj.fm.p(4)**2 - pj.mass**2)
!         p = max(p, 0.1d0) 
!         if( p .ge. 20.) then          
!            call cerrorMsg('Momentum is >20 GeV/c for ctotx',0)
!         endif

         if(pj.code .eq. knuc) then
            if(pj.charge .eq. 1) then
!            proton
!               shp = ctotpp1(p)
               call cppTotXs(p, shp)

               if(A .gt. 1.) then
                  call cnpTotXs(p, shn)
!               use average of pp,pn
!                  shp =( ctotpn1(p) + shp)/2.0
                  shp = (shp + shn)/2.0
                  call cxp2xAXsec(A, shp, xs)
               else
                  xs = shp
               endif
            elseif(pj.charge .eq. -1) then
!             pbar
!               shp = ctotpbarp1(p)
               call cpbarpTotXs(p, shp)
               if( A .gt. 1.) then
!                  shp = (shp + ctotpbarn1(p))/2.0
                  call cnbarpTotXs(p, shn)
                  shp = (shp + shn)/2.
                  call cxp2xAXsec(A, shp, xs)
               else
                  xs = shp
               endif
            elseif(pj.subcode .eq. antip) then
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
         elseif(pj.code .eq. kpion) then
            if(pj.charge .eq. 1) then
!               shp = ctotpiPp1(p)
               call cpippTotXs(p, shp)
            elseif(pj.charge .eq. -1) then
!               shp = ctotpiMp1(p)
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
         elseif(pj.code .eq. kkaon) then
            if(pj.charge .eq. 1) then
!               shp = ctotkPp1(p)
               call ckppTotXs(p, shp)
               if(A .gt. 1.) then
                  call ckpnTotXs(p, shn)
                  shp = (shp + shn)/2.0
                  call cxp2xAXsec(A, shp, xs)               
               else
                  xs = shp
               endif
            elseif(pj.charge .eq. -1) then
!               shp = ctotkMp1(p)
               call ckmpTotXs(p, shp)
               if(A .gt. 1.) then
                  call ckmnTotXs(p, shn)
                  shp = (shp + shn)/2.
                  call cxp2xAXsec(A, shp, xs)               
               else
                  xs = shp
               endif
            else
!             k0; don't worry so much 
!              xs =( ctotkMp1(p) +ctotkPp1(p))/2.0
               call ckmpTotXs(p, shp)
               call ckmnTotXs(p, shn)
               shp = (shp+ shn)/2.
               if(A .gt. 1.) then
                  call cxp2xAXsec(A, shp, xs)               
               else
                  xs = shp
               endif
            endif
         elseif(pj.code .eq. kgzai .or. pj.code .eq. ksigma .or.
     *           pj.code .eq. kbomega .or. pj.code .eq. klambda .or.
     *           pj.code .eq. klambdac) then
!          don't worry, not used almost at all; use proton
!            shp = ctotpp1(p)
            call cppTotXs(p, shp)
            call cxp2xAXsec(A, shp, xs)               
         elseif( pj.code == kgnuc ) then
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
