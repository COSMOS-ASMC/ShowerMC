!         testing eta' cetapDecay  
!       implicit none
!#include  "Zptcl.h"
!#include  "Zcode.h"
!      record /ptcl/ etap
!      record /ptcl/a(10)
!      integer n, i, j
!c       make eta
!      call cmkptc(ketap, 0, 0, etap)
!      etap.fm.p(1) = 1.
!      etap.fm.p(2) = -1.
!      etap.fm.p(3) = 100.
!      etap.fm.p(4) = sqrt(etap.fm.p(3)**2 +
!     *      etap.fm.p(1)**2 + etap.fm.p(2)**2 +  etap.mass**2)
!      do i = 1, 10000
!         call cetapDecay(etap, a, n)
!         do j =1 , n
!             write(*, *) sngl(a(j).fm.p(4)), a(j).code
!          enddo
!      enddo         
!      end
!    ******************************************************************
!    *                                                                *
!    *   cetapDecay: eta' decay
!    *                                                                *
!    ******************************************************************
!  eta'--> pi+pi-eta 42.6 %                          cum 0.426
!     --> pi+pi-gamma 28.9 %  (inc. rho0 + gamma);       0.715
!     --> pi0 pi0 eta   22.8 %                           0.943
!     --> g g          2.22                              0.965
!     --> omega gamma  2.62                              ~1.000
       subroutine cetapDecay(pj,  a,  np)
       implicit none
#include  "Zptcl.h"
#include  "Zcode.h"

       integer np               !output. no. of ptcls produced
       type(ptcl),intent(in):: pj 
       type(ptcl),intent(out):: a(*)
       integer  i, icon
       real*8 u, w
!
       call rndc(u)
       if(u .lt. 0.426d0) then
!           pi+ pi- eta
          call cmkptc(kpion, -1, 1, a(1))
          call cmkptc(kpion, 1, -1, a(2))
          call cmkptc(keta, 0, 0,  a(3) )
          call cnbdcy(3, pj%mass, a, 0, w, icon)
          np=3
       elseif(u .lt. .715d0 ) then
!            pi+ pi- gamma
          call cmkptc(kpion, -1, 1, a(1))
          call cmkptc(kpion, 1, -1, a(2))
          call cmkptc(kphoton, 0, 0, a(3))
          call cnbdcy(3, pj%mass, a, 0, w, icon)
          np = 3
       elseif( u .lt. 0.943d0 ) then
!           pi0 pi0 eta
          call cmkptc(kpion, 0, 0, a(1))
          call cmkptc(kpion, 0, 0, a(2))
          call cmkptc(keta, 0, 0, a(3))
          call cnbdcy(3, pj%mass, a, 0, w, icon)
          np = 3
       elseif( u .lt. 0.965d0) then
!     g g
          call cmkptc(kphoton, 0, 0, a(1))
          call cmkptc(kphoton, 0, 0, a(2))
          call c2bdcy(pj, a(1), a(2))
          np = 2
       else
!         omega gamma  
          call cmkptc(komega, 0, 0, a(1))
          call cmkptc(kphoton, 0,0, a(2))
          call c2bdcy(pj, a(1), a(2))
          np = 2
       endif
       if( np > 2 ) then
          do   i=1, np
             call cibst1(i, pj, a(i), a(i))
          enddo
       endif          
      end
