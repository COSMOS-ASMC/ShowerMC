!program main
!  implicit none
!  integer how, nr, ns
!  real(8)::a(0:3), u, v
!  complex(kind(0.d0))::x(3) 
!  integer i
!  how = 1
!  do while(.true.)
!!     write(0,*) 'Eneter a3,a2,a1,a0'
!!     read(*,*) a(3), a(2), a(1), a(0)
!     do i = 0, 3
!        call rndc(u)
!        call rndc(v)
!        u =int(u*1000.)
!        if(v < 0.5) then
!           u= -u
!        endif
!        a(i) = u
!     enddo
!!     call kcubicEq(a, how, x, nr, ns)
!     call kcubicEq(a,  x, nr, ns)
!     write(0,*) a(:)
!     write(0,*) nr
!     write(0,*) x(:)
!     do i = 1, 3
!        write(0,*) ( (a(3)*x(i) + a(2))*x(i) + a(1))* x(i) + a(0)
!     enddo
!     read(0,*)
!  enddo
!end program main
#include "Zcondc.h"
  subroutine kcubicEq(cf, x, nr, ns)
    implicit none
!                solve a cubic algebraic equation with real coef.
    real(8),intent(in)::cf(0:3)  ! cf(3)*x^3+cf(2)*x^2+cf(1)*x+cf(0) =0
                         !  cf(3) must not be 0.
!    integer,intent(in)::how  ! 0--> for real solution, Imag(x) is
                             ! made to be 0.  (normally use 0)
                             ! 1--> even for real solution, Imag(x)
                             ! is not cleared to see the accuracy of
                             ! numerial computation.
               ! not used now
    complex(kind(0d0)),intent(out)::x(3)  ! solution. 
                            ! first one is always real solution.
                            ! If nr=3, real x is sorted from small to large
    integer,intent(out):: nr  ! number of real solutions. (nr=1 or 3).
                              ! first nr of x are real solutions.
                      !  duplicated solutions are counted 2
    integer,intent(out):: ns ! nubmer of solutions. 
                             ! if cf(3)=0, ns may be < 3.
                             !  see kquadEq

    real(8)::p, q
#if defined MATHLOUSY
    real(8),parameter::sq3=1.73205080756887729353d0
#else
    real(8),parameter::sq3=sqrt(3.d0)
#endif
!!          the last (-1.0d0, -sq3)/2.d0 
!!   leads to error @ jaxa Fujitsu compiler; message is
!!   stragne; right parentheis is missing... 
!!   real reason :  -sq3 is NG.
!!   same notation for non complex is OK.
!!    complex(kind(0d0)),parameter::w(3) = (/(1.d0,0.d0), (-1.0d0,sq3)/2.d0, (-1.0d0, -sq3)/2.d0 /)
!!         next is workaround.  
    complex(kind(0d0)),parameter::w(3) = (/(1.d0,0.d0), (-1.0d0,sq3)/2.d0, -(1.0d0, sq3)/2.d0 /)

    complex(kind(0d0))::u(3), v(3), uu, vv
    complex(kind(0d0))::D
    complex(kind(0d0))::xx(3)  ! before sorted
    real(8)::aov3, a, b, c, temp, Delta, rx(3)
    integer::idx(3), i
    
    if( cf(3) == 0. ) then
       call  kquadEq(cf, x, nr, ns)
    else
   
       a = cf(2)/cf(3)
       b = cf(1)/cf(3)
       c = cf(0)/cf(3)   ! x^3 + ax^2 + bx +c = 0
   
       p = (b - a**2/3.d0)/3.d0
       q = (c + 2*a**3/27.d0 - a*b/3.d0)/2.d0
   
       Delta = -4*b**3 + b**2*a**2 - 4*c*a**3 + 18.0d0*c*b*a -27.0d0*c**2
       if(Delta > 0. ) then
          nr = 3
       elseif( Delta< 0. ) then
          nr = 1
       else
          nr= 3
       endif
       ns = 3
       D = q**2 + p**3
       if( real(D) < 0.d0 ) then
          uu = ( -q + sqrt(D) )**(1.d0/3.d0)
          vv = ( -q - sqrt(D) )**(1.d0/3.d0)
       else
          temp = -q + sqrt(D)
          uu = sign(abs(temp)**(1.d0/3.d0), temp)
          temp = -q - sqrt(D)
          vv = sign(abs(temp)**(1.d0/3.d0), temp)
       endif
       u(:) = w(:)*uu
       v(:) = w(:)*vv
   
       aov3 = a/3.d0
       xx(1) = u(1) + v(1) - aov3
   
       xx(2) = u(2) + v(3) - aov3
   
       xx(3) = u(3) + v(2) - aov3
   
       if( nr == 3 ) then
          rx(:)=real(xx(:))
   !             sort ascending order
          call kcsr1idx(rx,  3, idx, 'a')
          do i = 1, 3 
             x(i) = xx(idx(i))
             if( cf(0) == 0. ) then
                if( abs(x(i)) < 1.d-13 ) then
                   x(i) = cmplx(0.d0, 0.d0)
                endif
             endif
          enddo
          
!          if( how == 0 ) then
             x(1)=cmplx(real(x(1)), 0.d0)
             x(2)=cmplx(real(x(2)), 0.d0)
             x(3)=cmplx(real(x(3)), 0.d0)
!          endif
   
       else
          x(:)=xx(:)
          temp = min( abs(aimag(x(1))),abs(aimag(x(2))),abs(aimag(x(3))) )
          if( temp == abs(aimag(x(1)))) then
          elseif(  temp == abs(aimag(x(2)))) then
             x(1) = x(2)
             x(2) = conjg(x(3))
          else
             x(1)=x(3)
             x(3) = conjg(x(2))
          endif
!          if(how == 0) then
             x(1)=cmplx(real(x(1)), 0.d0)
!          endif
          if( cf(0) == 0. ) then
             x(1)=cmplx(0.d0, 0.d0)
          endif
       endif
    endif
  end subroutine kcubicEq
