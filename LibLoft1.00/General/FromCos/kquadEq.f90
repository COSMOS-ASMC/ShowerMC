!  program main
!    implicit none
!    real(8)::cf(0:2), u, v
!    complex(kind(0d0))::x(2) 
!    integer ns, i, nr, j
!
!    do while (.true.)
!!       write(0,*) 'Enter a b c'
!!       read(*,*) cf(2), cf(1), cf(0)
!       do j= 0, 2
!          call  rndc(u)
!          call  rndc(v)
!          u = int( u*1000 )/100.
!          if( v< 0.5) u= -u
!          cf(j) = u
!       enddo
!       call kquadEq(cf, x, nr, ns)
!       write(0,*) 'cf =', cf(:)
!       write(0,*) ' ns nr=', ns, nr
!       do i = 1, ns
!          write(0,*) 'x=',x(i)
!          write(0,*) (cf(2)*x(i) + cf(1))*x(i) + cf(0)
!       enddo
!       read(*,*)
!    enddo
!  end program main
  subroutine kquadEq(cf, x, nr, ns)
    implicit none
    real(8),intent(in)::cf(0:2)
    complex(kind(0d0)),intent(out)::x(2)  ! when nr=2, 
                                   ! real(x(1)) <= real(x(2))
    integer,intent(out)::nr  ! # of real solution
                  ! 1--> cf(2)=0 and the eq. is linear eqation
                  ! 2--> 2 real solutions (may be dup.)
                  ! 0--> 2 complex solutions
                  !  x(2) = conjugate of  x(1)  
    integer,intent(out)::ns ! # of solutions
                  !  nr=1-->ns=1
                  !  nr=0 or 2-->ns=2
                  !  ns = 0.  bad euation.
    real(8)::D
    complex(kind(0d0)):: temp

    if(cf(2) == 0.0d0) then
       call klinEq(cf, x, ns)
       if(ns > 0 ) then
          nr = 1
       else
          nr = 0
       endif
    else
       D = cf(1)**2 - 4.0d0*cf(2)*cf(0)
       if( D >= 0. )  then
          nr = 2
          temp = -( cf(1)+ sign(sqrt(D), cf(1)))/2.d0
          x(1) = temp/cf(2)
          x(2) = cf(0)/temp
          if(real(x(1)) > real(x(2)) ) then
             temp = x(1)
             x(1) = x(2)
             x(2) = temp
          endif
          if(cf(0) == 0.) then
             if( abs(x(1)) < abs(x(2)) ) then
                x(1) = 0.
             else
                x(2) = 0.
             endif
          endif
       else
          nr = 0
          
          temp = cmplx(0.d0, sqrt(-D) )
          x(1) = (-cf(1) + temp)/2/cf(2)
          x(2) = conjg(x(1))
       endif
       ns = 2
    endif
  end subroutine kquadEq
  subroutine  klinEq(cf, x, ns)
    implicit none
    real(8),intent(in)::cf(0:1)
    complex(kind(0d0)),intent(out)::x(1) 
    integer,intent(out)::ns ! ns=1. ok
                            ! ns=0. cf(1)=0.
    if(cf(1) == 0.d0) then
       ns = 0
    else
       x(1) = -cf(0)/cf(1)
       ns = 1
    endif
  end subroutine klinEq
