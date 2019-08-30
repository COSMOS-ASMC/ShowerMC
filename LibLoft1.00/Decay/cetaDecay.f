!         testing cetaDecay
!       implicit none
!#include  "Zptcl.h"
!#include  "Zcode.h"
!      record /ptcl/ eta
!      record /ptcl/a(10)
!      integer n, i, j
!c       make eta
!      call cmkptc(keta, 0, 0, eta)
!      eta.fm.p(1) = 1.
!      eta.fm.p(2) = -1.
!      eta.fm.p(3) = 100.
!      eta.fm.p(4) = sqrt(eta.fm.p(3)**2 +
!     *      eta.fm.p(1)**2 + eta.fm.p(2)**2 +  eta.mass**2)
!      do i = 1, 10000
!         call cetaDecay(eta, a, n)
!         if(n .ne. 2) then
!            do j =1 , n
!               write(*, *) sngl(a(j).fm.p(4)), a(j).code
!            enddo
!         endif
!      enddo         
!      end
!    ******************************************************************
!    *                                                                *
!    *   cetaDecay: eta decay
!    *                                                                *
!    ******************************************************************
!  eta--> gg 39.41 %
!     --> 3 pi0 31.8 % --> 32.68 (2016)
!     --> pi+ pi- pi0 23.7 %-->22.92 (2016)
!     --> pi+pi- g     4.22%  (2016)
!     --> mu+ + mu- + gamma 3.1x10^-4 (not %)
!     --> mu+ + mu-         5.8 10^-6 (not %)     
      
       subroutine cetaDecay(pj,  a,  np)
       implicit none
#include  "Zptcl.h"
#include  "Zcode.h"

       integer np               !output. no. of ptcls produced
!     record /ptcl/ pj         ! input. eta
       type(ptcl),intent(in):: pj 
!     record /ptcl/ a(*)      ! output. produced ptcls
       type(ptcl),intent(out):: a(*)
       integer  i, icon
       real*8 u, w
!
       call rndc(u)
       if(u .lt. .3941d0) then
!           gg
          call cmkptc(kphoton, 0, 0, a(1))
          call cmkptc(kphoton, 0, 0, a(2))
          call c2bdcy(pj, a(1), a(2))
          np=2
       elseif(u .lt. .7209d0 ) then
!           3 pi0
          call cmkptc(kpion, 0, 0, a(1))
          call cmkptc(kpion, 0, 0, a(2))
          call cmkptc(kpion, 0, 0, a(3))
          call cnbdcy(3, pj%mass, a, 0, w, icon)
          np = 3
          do   i=1, np
             call cibst1(i, pj, a(i), a(i))
          enddo
       elseif( u .lt. 0.9501d0 ) then
!           pi+ pi- pi0
          call cmkptc(kpion, 0, 1, a(1))
          call cmkptc(kpion, 0, -1, a(2))
          call cmkptc(kpion, 0, 0, a(3))
          call cnbdcy(3, pj%mass, a, 0, w, icon)
          np = 3
          do   i=1, np
             call cibst1(i, pj, a(i), a(i))
          enddo
       elseif( u .lt. 0.9923d0) then
!     pi+ pi- + g
          call cmkptc(kpion, 0, 1, a(1))
          call cmkptc(kpion, 0, -1, a(2))
          call cmkptc(kphoton, 0, 0, a(3))
          call cnbdcy(3, pj%mass, a, 0, w, icon)
          np = 3
          do   i=1, np
             call cibst1(i, pj, a(i), a(i))
          enddo
       elseif(  u .lt. 0.99261d0) then
!     mu+ + mu- + g
          call cmkptc(kmuon, 0, 1, a(1))
          call cmkptc(kmuon, 0, -1, a(2))
          call cmkptc(kphoton, 0, 0, a(3))
          call cnbdcy(3, pj%mass, a, 0, w, icon)
          np = 3
          do   i=1, np
             call cibst1(i, pj, a(i), a(i))
          enddo
       elseif( u .lt. 0.9926158d0) then
!     mu+ + mu- 
          call cmkptc(kmuon, 0, 1, a(1))
          call cmkptc(kmuon, 0, -1, a(2))
          np = 2
          call c2bdcy(pj, a(1), a(2))
       else
! all rest : regards e+e-g  ~ 6.9e-3
          call cmkptc(kelec, 0, 1, a(1))
          call cmkptc(kelec, 0, -1, a(2))
          call cmkptc(kphoton, 0, 0, a(3))
          call cnbdcy(3, pj%mass, a, 0, w, icon)
          np = 3
          do   i=1, np
             call cibst1(i, pj, a(i), a(i))
          enddo
      endif
      end
