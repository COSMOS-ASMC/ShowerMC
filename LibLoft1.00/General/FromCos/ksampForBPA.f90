!        sampling from functions
!   ksampForBPA1:   1/(1+x)**2 dx
!   ksampForBPA2:   x/(1+x)**4 dx
!
!  These are related to brems photon emission angle or
!  pair creation electron angle  
!
  subroutine ksampForBPA2(x)
    implicit none
    real(8),intent(out):: x  ! sampled random number
    
    real(8):: coef(0:3)
    complex(kind(0d0)):: sol(3)
    real(8):: u
    integer:: nr, ns
    
    call rndc(u)
    coef(0) = u -1.
    coef(2) = 3*u
    coef(1) = coef(2) - 3.0d0
    coef(3) = u
    
    call  kcubicEq(coef, sol, nr, ns)
    if(nr == 3 ) then  ! always here. 
       if(real(sol(3)) < 0. ) then   ! not happen
          write(0,*) ' u ',u, 'nr=',nr
          write(0,*) ' sol=',sol
          stop
       else
          x =real( sol(3) )
       endif
    elseif( nr == 1 ) then
       if(real(sol(1)) < 0.) then
          write(0,*) ' u ',u, 'nr=',nr
          write(0,*) ' sol=',sol
          stop
       else
          x =real( sol(1) )
       endif
    else
       write(0,*) ' u ',u, 'nr=',nr
       write(0,*) ' sol=',sol
       stop
    endif
!!    write(0,'(a,i2,1p,6g15.6)') 'nr ', nr, sol
  end subroutine ksampForBPA2
  subroutine ksampForBPA1(x)
!   samples random x from 1/(1+x)**2 dx
    implicit none
    real(8),intent(out):: x
    real(8):: u

    call rndc(u)
    x = 1.0d0/u - 1.
  end subroutine ksampForBPA1
!  program main
!    implicit none
!    integer:: i
!    real(8):: x
!
!    do i = 1, 100000
!!       call ksampForBPA2(x)
!       call ksampForBPA1(x)
!       write(*,*) x
!    enddo
!  end program main
       
