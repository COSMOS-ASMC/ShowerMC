!     cinelx: x A / AA' inelastic xsection
!        
      subroutine cinelx(pj, A, Z, xs)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmass.h"
#include "Zevhnp.h"
      type(ptcl)::pj   !input projectile particle
      real(8),intent(in)::  A     ! input. effective target mass number
      real(8),intent(in)::  Z     ! input. effective target charge number
      real(8),intent(out):: xs    ! output. in mb. inelastic xs.

      real(8):: gamma
      type(ptcl):: target
      integer:: tcharge
      real(8):: pjA, pjZ


      if( A == 1.0 .and. pj%code == kgnuc  ) then ! target is p/n, pj is heavy ion
               ! we use proj. rest system ;  if not, inside cinelx0 
               ! cAAXsec2 is called and it calles cinelx and recursive
               ! call happens.
         tcharge = Z
         call cmkptc(knuc, -1, tcharge, target)  ! p or n
         gamma = pj%fm%p(4)/pj%mass       
         pjA = pj%subcode
         pjZ = pj%charge   
         if( tcharge == 0 ) then
            target%fm%p(4) = gamma*masn
         else
            target%fm%p(4) = gamma*masp
         endif
         call cinelx0(target, pjA, pjZ, xs)   ! target is now proj. pjA is target

      else
         call cinelx0(pj, A, Z,  xs)      
      endif
      end

      subroutine cinelx0(pj, A, Z, xs)
      use modcosXsec  ! refer only nIsDiff
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zevhnp.h"
      type(ptcl)::pj   !input projectile particle
      real(8),intent(in)::  A          ! input. effective target mass number
      real(8),intent(in)::  Z     ! input. effective target charge number
      real(8),intent(out):: xs          ! output. in mb. total xsection. 


      real(8) p
      real(8),external:: cPDGsigmaInepA
      real(8),save:: Aprev=-100., norm
      real(8):: xs200, xs200A


      real(8)::shp, shn

      if( pj%fm%p(4) .le. pj%mass) then
         p =0.
         if(pj%code .eq. knuc .and. pj%subcode .eq. antip) then
            xs = largexs
         elseif(pj%code .eq. kelec .and. pj%charge .eq. 1) then
            xs = largexs
         else
            xs  = smallxs
         endif
      else   
         p = sqrt(pj%fm%p(4)**2 - pj%mass**2)
         if(pj%code .eq. knuc) then
            if(pj%charge .eq. 1) then
               if( A == 1.0 ) then
                  if( Z == 1.0 ) then
                     call cppInelaXs(p, xs)
                  elseif( Z == 0.) then
                     call cnpInelaXs(p, xs)
                  else
                     write(0,*) 'A,Z=',A,Z, 'for p pj'
                     write(0,*) 'in cinelx; invalid'
                     stop
                  endif
               elseif(A > 1.) then
                  call cppInelaXs(p, shp)
                  if( nIsDiff ) then 
                     call cnpInelaXs(p, shn)
                     shp = (shp + shn)/2.0
                  endif
                  call cxp2xAXsec(A, shp, xs)
               else
                  write(0,*) ' A,Z=', A, Z
                  write(0,*) ' for  cinelx invalid'
                  write(0,*) ' pj is p '
                  stop
               endif
            elseif(pj%charge .eq. -1) then
               call cpbarpInelaXs(p, shp)
               if( A .gt. 1.) then
                  if(nIsDiff) then
                     call cnbarpInelaXs(p, shn)
                     shp = (shp + shn)/2.
                  endif
                  call cxp2xAXsec(A, shp, xs)
               else
                  xs = shp
               endif
            elseif(pj%subcode .eq. antip) then
!              anti-neutron; assume the same one as pbar
               call cnbarpInelaXs(p, shp)
               if( A .gt. 1.) then
                  call cxp2xAXsec(A, shp, xs)
               else
                  xs = shp
               endif
            else
!               neutron
               call cnpInelaXs(p, shp)
               if(A .gt. 1.) then
                  call cxp2xAXsec(A, shp, xs)
               else
                  xs = shp
               endif
            endif
         elseif(pj%code == kpion .or. pj%code == keta) then
            if(pj%charge .eq. 1) then
               call cpippInelaXs(p, shp)
            elseif(pj%charge .eq. -1) then
               call cpimpInelaXs(p, shp) 
            else 
                  ! pi0 or eta
               call cpimpInelaXs(p, shp)
            endif
            if(A .ne. 1.0) then
               call cxp2xAXsec(A, shp, xs)
            else
               xs = shp
            endif
         elseif(pj%code .eq. kkaon) then
            if(pj%charge .eq. 1) then
               if(A == 1.0) then
                  if( Z==1.0) then
                     call ckppInelaXs(p, xs)
                  else
                     call ckpnInelaXs(p, xs)
                  endif
               elseif(A .gt. 1.) then
                  call ckppInelaXs(p, shp)
                  if( nIsDiff ) then
                     call ckpnInelaXs(p, shn)
                     shp = (shp + shn)/2.0
                  endif
                  call cxp2xAXsec(A, shp, xs)               
               else
                  write(0,*) 'A,Z=',A,Z, ' invalid for cinel%Kaon'
                  stop
               endif
            elseif(pj%charge .eq. -1) then
               call ckmpInelaXs(p, shp)
               if(A .gt. 1.) then
                  if( nIsDiff ) then
                     call ckmnInelaXs(p, shn)
                     shp = (shp + shn)/2.
                  endif
                  call cxp2xAXsec(A, shp, xs)               
               else
                  xs = shp
               endif
            else
!             k0; don't worry so much 
               call ckmpInelaXs(p, shp)
!!               call ckmnInelaXs(p, shn)
!!               shp = (shp+ shn)/2.
               if(A .gt. 1.) then
                  call cxp2xAXsec(A, shp, xs)               
               else
                  xs = shp
               endif
            endif
         elseif(pj%code .eq. kgnuc ) then
!                heavy ion 
            call cAAXsec2(pj, A, Z,  xs)
         elseif(pj%code .eq. kdmes) then
            call ckppInelaXs(p, shp)
            call cxp2xAXsec(A, shp, xs)               
         elseif(pj%code .eq. kgzai .or. pj%code .eq. ksigma .or.
     *           pj%code .eq. kbomega .or. pj%code .eq. klambda .or.
     *           pj%code .eq. klambdac) then
!          don't worry, not used almost at all; use proton
!            shp = ctotpp1(p)
            call cppInelaXs(p, shp)
            call cxp2xAXsec(A, shp, xs)
         elseif(pj%code == kphoton ) then
            call cgpXsec(A, pj%fm%p(4), xs)
         elseif(pj%code .eq. kneumu) then
            xs= smallxs
         elseif(pj%code .eq. kneue) then
            xs= smallxs
         elseif(pj%code == knnb .or. pj%code == kddb ) then
            xs = smallxs
         elseif(pj%code == kmuon) then
            write(0,*) 'cinelx: code=',pj%code, ' should not come'
            stop
         else
!          use pion
            call cpippInelaXs(p, shp)
            call cxp2xAXsec(A, shp, xs)               
         endif
      endif
!!!       normalize to pdg at 200 GeV. norm const is
!!!       applied to other collisions 
      if(A > 1. .and. A /=  Aprev ) then
         call cppInelaXs(200.d0, xs200)
         call cxp2xAXsec(A, xs200, xs200A)
!!         write(0,*) ' xs200A=',xs200A, ' pgd=',  cPDGsigmaInepA(A)
         norm =  cPDGsigmaInepA(A) /xs200A
!!         write(0,*) ' norm=', norm
         Aprev = A
      elseif( A == 1.)  then
         norm = 1.
         Aprev = A
      endif
      xs = xs*norm

      end
      function cinelCosByPdg(A)  result(ratio)
!             sigmaInela(pA)cosmos/sigmaInela(pA)PDG
!        @ 200 GeV.  SigmaInela by cosmos is the one
!             obtained by cinelx0  before  
!             cxp2xAXsec  included cxAbyxpXsec.
!       
      real(8),intent(in):: A
      real(8):: ratio

      real(8):: temp

      if( A > 48. ) then
         if(A >150.) then
            ratio = 1.0
         else
            ratio = (1.006+0.0268*A/100.)*0.957
         endif
      elseif(A > 20.) then
         ratio = 0.975
      elseif( abs(A - 16.) < 0.1 ) then
         ratio = 0.9427
      elseif( abs( A- 7.0) < 0.1 ) then
         ratio = 0.902
      elseif( A > 7. ) then 
         temp = log(A)
         ratio = (-0.100*temp + 0.57311)*temp +0.1578
      elseif( abs( A- 4.) < 0.1  ) then
         ratio = 0.975
      elseif( abs( A- 2.) < 0.1 ) then
         ratio = 1.137
      elseif (abs(A-1.0) < 0.1 ) then
         ratio =1.0138
      else
         ratio = -0.025*(A-4) + 0.975
      endif
      end

            
