!  ifort kquarticEq.f90 kcubicEq.f90 kquadEq.f90 kcquadEq.f90 -L$COSMOSTOP/lib/PCLinuxIFC64 -lcosmos -L$EPICSTOP/lib/PCLinuxIFC64 -lepics
!  program main
!    implicit none
!    real(8)::cf(0:4), u, v
!    complex(kind(0d0))::x(4) 
!    integer ns, i, nr, j
!
!    do while (.true.)
!!       write(0,*) 'Enter a b c d e'
!!       read(*,*) cf(4), cf(3), cf(2), cf(1), cf(0)
!       do j= 0, 4
!          call  rndc(u)
!          call  rndc(v)
!          u = int( u*1000 )/100.
!          if( v< 0.5) u= -u
!          cf(j) = u
!       enddo
!!       cf(3) = 0.
!       call kquarticEq(cf, x, nr, ns)
!       write(0,'(5g14.5)') 'cf =', cf(:)
!       write(0,*) ' ns nr=', ns, nr
!       do i = 1, ns
!          write(0,*) 'x=',x(i)
!          write(0,*) (((cf(4)*x(i) + cf(3))*x(i) + cf(2))*x(i) + cf(1))*x(i)+cf(0)
!       enddo
!       read(*,*)
!    enddo
!  end program main

  subroutine kquarticEq(cf, x, nr, ns)
    implicit none
    real(8),intent(in)::cf(0:4)
    complex(kind(0.d0)),intent(out):: x(4)
    integer,intent(out):: nr  ! # of real solutions
    integer,intent(out):: ns  ! # of solutions. ns could be 
                              ! smaller than 4 if cf(4)=0 etc

    complex(kind(0d0)) :: y1(2), y2(2), y(4)
    complex(kind(0d0)) :: cfc(0:2)
    real(8)::yr(4)
    real(8)::cf2(0:2)
    real(8)::cf3(0:4) 
    integer:: i
    integer:: idx(4)
    complex(kind(0.d0))::u(3), temp
    integer nr1, nr2, ns1, ns2
    real(8)::B, p, q, r
    
    if( cf(4) == 0.d0) then
       call kcubicEq(cf, x, nr, ns)
    else
       if( cf(1) == 0.  .and. cf(3) == 0.) then
!            ax^4 + bx^2 + c = 0 type
          cf2(2) = cf(4) 
          cf2(1) = cf(2)
          cf2(0) = cf(0)
          call kbiquadEq(cf2, x, nr, ns)
       else
          cf3(0:4)=cf(0:4)/cf(4)
          B = cf3(3)/4.d0
          p = (cf3(2)-6.d0*B**2)
          q = (cf3(1)-2*cf3(2)*B + 8.d0*B**3)
          r = ( (-3.d0*B*B+cf3(2))*B -cf3(1))*B + cf3(0)
          if( abs(q) <= 1.0d-9 ) then
             cf2(0)=r
             cf2(1)=p
             cf2(2)=1.
             call kbiquadEq(cf2, y, nr, ns)
             do i = 1, ns
                x(i)= y(i) -B
             enddo
          else
             ! solve u(p+u)^2 -4ru=q^2
             ! i.e:   u^3 +2pu^2 + (p^2-4r)u -q^2 = 0
             cf3(0) = -q**2
             cf3(1) = (p**2 - 4*r)
             cf3(2) = 2*p
             cf3(3) = 1.
             call kcubicEq(cf3, u, nr, ns)
             ! first one should be real so we use it and solve
             ! y^2 -/+ sqrt(u)y +(p+u)/2 +/-sqrt(u)q/(2u)
             cfc(2) = 1.
             temp = sqrt(u(1)) 
             cfc(1) = -temp
             cfc(0) = (p+u(1))/2 + temp*q/2/u(1)
             call kcquadEq(cfc, y1, nr1, ns1)
             cfc(1) = temp
             cfc(0) = (p+u(1))/2 - temp*q/2/u(1)
             call kcquadEq(cfc, y2, nr2, ns2)
             nr = 0
             do i = 1, nr1
                nr = nr + 1
                yr(nr) = y1(i)
             enddo
             do i= 1, nr2
                nr = nr + 1
                yr(nr) = y2(i)
             enddo
             call kcsr1idx(yr,  nr, idx,  'a')
             ns = 0
             do i = 1, nr
                ns = ns + 1
                x(ns) = yr(idx(i)) -B
             enddo
             do i = nr1+1, ns1
                ns = ns + 1
                x(ns) = y1(i)-B
             enddo
             do i = nr2+1, ns2
                ns = ns + 1
                x(ns) = y2(i)-B
             enddo
          endif
       endif
    endif
   end subroutine kquarticEq
    

  subroutine kbiquadEq(cf, x, nr, ns)
!       solve  cf(2)x^4 + cf(1)x^2 + cf(0)= 0
    implicit none
    real(8),intent(in)::cf(0:2) ! real coef.
    complex(kind(0d0)),intent(out)::x(4)
    integer,intent(out)::nr  ! # of real solutions
    integer,intent(out)::ns  ! # of solutions.

    complex(kind(0d0)) :: y(2)
    integer:: nr2, ns2, i, j, np
    integer::cr(4)  ! if xx is complex(-->1  real -->0)
    integer::idx(4) ! index for cr
    integer::idx2(4) ! index for rx
    real(8)::rx(4), xx(4)

!       solve cf(2)y^2 + cf(1)y + cf(0)
    call kquadEq(cf, y, nr2, ns2)
    ns = 0
    nr = 0
!      here y=x^2
!      get 4 solutions by x=+/-sqrt(y(:));  
    do i = 1, ns2
       ns = ns + 1
       xx(ns) = sqrt(y(i))
       if(nr2 >= i .and. real(y(i)) >= 0.) then
          cr(ns) = 0
          nr = nr + 1   ! real solution
          ns = ns + 1
          cr(ns) = 0
          nr = nr + 1   ! real solution
          xx(ns) = -xx(ns-1)
       else
          cr(ns) = 1
          ns = ns + 1
          xx(ns) = -xx(ns-1)
       endif
    enddo

    call kqsorti(cr, idx, ns)
    np = 0
    do i = 1, ns
       if( cr(idx(i)) == 0 ) then
          np = np + 1
          rx(np)=real(xx(idx(i)))
       endif
    enddo
   !    there are np real solutions;   sort ascending order
    call kcsr1idx(rx,  np, idx2, 'a')

    do i = 1, np
       x(i) = rx(idx2(i))
    enddo
    j = np
    do i = 1, ns
       if(cr(i) /= 0 ) then
          j = j + 1
          x(j) = xx(idx(i))
       endif
    enddo
  end subroutine kbiquadEq
