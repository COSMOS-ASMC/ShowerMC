  subroutine csampHFMom(A,Z, A1,Z1, p1, p2)
    use NucMassFormula
!    samples heavy fragmen momentum (fission)
    implicit none
    integer,intent(in):: A, Z ! spectator A,Z
          ! this spectator is assumed to break into A1,Z1 and  (A2,Z2)
          !            where   A2= A-A1, Z2=Z-Z1
    integer,intent(in):: A1,Z1 ! A,Z of fragment escaping from spectator
    real(8),intent(out):: p1(3) ! sampled momentum in GeV/c of fragment 1
    real(8),intent(out):: p2(3) ! sampled momentum in GeV/c of residual nucleus
    
    real(8),parameter:: alfa=1./137., hbarc=0.200  
    real(8),parameter:: Ef= hbarc*alfa   ! GeV Fermi
    real(8):: r1, r2, KE, m1, m2,  T
    integer:: A2, Z2
    real(8)::u,  pabs, cosu, sinu, cosf, sinf
    A2 = A - A1
    Z2 = Z - Z1
    if(Z2 < 1 .or. A2 < 1) then
       pabs = 0.001
    else
       r1 = 1.2* dble(A1)**0.3333
       r2 = 1.2* dble(A2)**0.3333
       KE = Z1*Z2* Ef/(r1+r2)*2   ! GeV
       m1 = cNucMass(A1, Z1)/1000.  ! GeV
       m2 = cNucMass(A2, Z2)/1000. 
       T = KE + m1 + m2
       pabs =sqrt( ( (T**2 - (m1+m2)**2) * (T**2 - (m1-m2)**2) )/(4*T**2))
    endif

    call rndc(cosu)
    cosu = 2*cosu - 1
    sinu = sqrt(1.d0 - cosu**2)
    call kcossn(cosf, sinf)
    p1(1) = pabs*sinu*cosf
    p1(2) = pabs*sinu*sinf
    p1(3) = pabs*cosu
    p2(1:3) =-p1(1:3)
  end subroutine csampHFMom

!  program main
!    use NucMassFormula
!    implicit none
!    integer:: A, Z ! spectator A,Z
!    integer:: A1,Z1 ! A,Z of fragment escaping from spectator
!
!    real(8):: p1(3) ! sampled momentum in GeV/c of fragment 1
!    real(8):: p2(3) ! sampled momentum in GeV/c of residual nucleus
!
!    real(8):: m1, m2
!    integer::i
!
!    write(0,*) 'Enter A, Z, A1, Z1'
!    read(*,*) A,Z,A1,Z1
!    write(0,*) ' A,Z =', A,Z
!    write(0,*) ' A1,Z1 =', A1,Z1
!    m1 = cNucMass(A1, Z1)/1000.  ! GeV
!    m2 = cNucMass(A-A1, Z-Z1)/1000.
!    write(0,*) ' m1,m2=',m1,m2
!
!    do i = 1, 100000 
!       call csampHFMom(A,Z, A1,Z1, p1, p2)
!       write(*,'(1p,8g14.4)') &
!         p1(:), p2(:), sqrt(dot_product(p1,p1)+m1**2)-m1, &
!         sqrt(dot_product(p2,p2)+m2**2) - m2 
!    enddo
!  end program main
