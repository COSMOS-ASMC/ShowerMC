!  ifort  kcquadEq.f90 -L$COSMOSTOP/lib/PCLinuxIFC64 -lcosmos
!   program main
!     implicit none
!     real(8)::u, v, w
!     complex(kind(0d0))::x(2), cf(0:2) 
!     integer ns, i, nr, j
! 
!     do while (.true.)
!        write(0,*) 'Enter a b c'
!        read(*,*) cf(2), cf(1), cf(0)
!!        cf(2)=(1,0)
!!        cf(1)=-(1,1)
!!        cf(0)=(0,1)
!        do j= 0, 2
!           call  rndc(u)
!           call  rndc(v)
!           call  rndc(w)
!           u = int( u*1000 )/100.
!           v = int( v*1000 )/100.
!           if( w< 0.25) then
!              u= -u
!              v= 0.
!           elseif( w< 0.5 ) then
!              u = 0.
!           elseif( w< 0.75) then
!              v = -v
!           endif
!!           cf(j) = cmplx(u,v,8)
!        enddo
!        call kcquadEq(cf, x, nr, ns)
!        write(0,*) 'cf =', cf(:)
!        write(0,*) ' ns nr=', ns, nr
!        do i = 1, ns
!           write(0,*) 'x=',x(i)
!           write(0,*) (cf(2)*x(i) + cf(1))*x(i) + cf(0)
!        enddo
!        read(*,*)
!     enddo
!  end program main
  subroutine kcquadEq(cf, x, nr, ns)
!         coef. could be complex
    implicit none
    complex(kind(0d0)),intent(in)::cf(0:2)
    complex(kind(0d0)),intent(out)::x(2)  ! when nr=2, 
                                   ! real(x(1)) <= real(x(2))
    integer,intent(out)::nr  ! # of real solution
    integer,intent(out)::ns ! # of solutions
    real(8),parameter::zero=1.d-9
    complex(kind(0d0))::D, temp
    real(8)::rcf(0:2)

    if(abs(cf(2)) == 0.0d0) then
       if(abs(cf(1)) == 0.d0)  then
          ns = 0
          nr = 0
       else
          ns = 1
          x(1) = -cf(0)/cf(1)
          if( abs(aimag(x(1))) < zero) then
             nr = 1
          else
             nr = 0
          endif
       endif
    else
       D = cf(1)**2 - 4.0d0*cf(2)*cf(0)
       if( aimag(cf(0)) == 0.d0 .and.  &
            aimag(cf(1)) == 0.d0 .and. &
            aimag(cf(2)) == 0.d0 ) then
          rcf(:) = real(cf(:))
          call kquadEq(rcf, x, nr, ns)
       else
          temp = sqrt(D)
          x(1) = (-cf(1) + temp)/2/cf(2)
          x(2) = (-cf(1) - temp)/2/cf(2)
          ns = 2
          nr = 0
          if(abs(aimag(x(1))) < zero) then
             x(1) = cmplx(real(x(1)), 0.d0, 8) 
             nr = nr + 1
          endif
          
          if(abs(aimag(x(2))) < zero) then
             x(2) = cmplx(real(x(2)), 0.d0, 8) 
             nr = nr + 1
          endif
          if( nr == 2 ) then
             if(real( x(1) ) > real( x(2) ) ) then
                temp = x(1)
                x(1) = x(2)
                x(2) = temp
             endif
          elseif( nr == 1) then
             if( aimag(x(2))< zero ) then
                temp = x(1)
                x(1) = x(2)
                x(2) = temp
             endif
          endif
       endif
    endif
  end subroutine kcquadEq
