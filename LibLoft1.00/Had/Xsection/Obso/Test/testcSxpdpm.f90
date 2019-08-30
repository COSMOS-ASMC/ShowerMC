  program main
    use modDPMXsec
    ! test cSxpdpm.f
    implicit none
    real(8):: E, s, roots, p, xs, m1, m2
    real(8),parameter:: dE=10.d0**0.1d0
    real(8),parameter:: mpi=0.13957
    real(8),parameter:: mp=0.93827
    real(8),parameter:: mK=0.4936
    
    E = 50.
    m1 = mp  !  change this for p,pi,K
    m2 = mp
    do while (E < 1.e11)
       p = sqrt(E**2 - m1**2)
       s = m1**2 + m2**2 + 2*m2*E
       roots = sqrt(s)
       call cSppdpm(E, xs)   ! change name for p,pi,K
       write(*,'(1p, 5g14.4)')  E, xs, p, s, roots
       E = E *dE
    enddo
  end program main
