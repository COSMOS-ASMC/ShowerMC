!           to test ksbwig: sampling of breight-wigner type function
!       implicit none
!       integer i
!       real*8 e  
!       do i=1, 50000
!          call ksbwig(770.d0, 150.d0, e)
!          write(*, *) sngl(e)
!       enddo
!       end
!       *************************************************************
!       *
!       * ksbwig: give a random variable followng the breight-wigner
!       *         type distribution
!       *
!       *    note: no cut on the both sides of the distribution is made.
!       *         The user should do it for physical applications
!       *
!       *************************************************************
!
!       ksbwig(e0, g, e)
! distribution funcition is
!      de/( (e-e0)**2 + (g/2)**2)
!
        subroutine ksbwig(e0, g, e)
        implicit none
        real*8 e0, g, e
        real*8  u, pi

        parameter (pi=3.141592653589793238)
        call rndc(u)
        e=g/2 * tan( (u-.5)*pi ) + e0
        end
