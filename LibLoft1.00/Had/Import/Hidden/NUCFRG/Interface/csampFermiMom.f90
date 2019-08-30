  subroutine csampFermiMom(A,Z, chg, p)
    implicit none
    integer,intent(in):: A, Z !  A and Z of nucleus
    integer,intent(in):: chg  ! 0 or 1 for n or p
    real(8),intent(out):: p(3) ! sampled momentum in GeV/c
    
    real(8)::u1,u2, u3, pabs, cosu
    
    call rndc(u1)
    call rndc(u2)
    call rndc(u3)
    pabs = max(u1, u2, u3)* 0.272  ! sample from x^2dx
    if(chg == 1 ) then
       pabs = pabs * ( dble(Z)/A )**0.333333
    else
       pabs = pabs * (dble(A-Z)/A)**0.333333
    endif
    call kcossn(u1, u2)  ! cosfai,sinfai 
    call rndc(u3)  
    cosu = 2*u3 - 1    !  for cos dcos
    p(3) = pabs*cosu
    u3 = sqrt(1.d0-cosu**2)   ! sin teta
    p(1) = pabs*u3*u1
    p(2) = pabs*u3*u2
  end subroutine csampFermiMom
!  program main
!    implicit none
!    integer i, A, Z, chg
!    real(8):: p(3)
!    write(0,*) 'Enter A,Z, chg'
!    read(*,*)  A, Z, chg
!    do i = 1, 100000
!       call csampFermiMom(A, Z, chg, p)
!       write(*,'(1p,4g14.5)') p(:),sqrt(dot_product(p,p))
!    enddo
!  end program main
