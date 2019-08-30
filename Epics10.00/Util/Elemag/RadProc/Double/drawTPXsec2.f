!         This uses the fact that probability f(v, Eg, rho)dv
!     scales as f(v, Eg * rho)dv
!     The range of Eg is 10^8 GeV to 10^13 GeV
!     and rho= 1.e-3 to 1.e-8 g/cm3
!     The X = Eg * rho is
!          1 to 10^10
!     We fix the energy; 10^19 eV and variate rho
!       
      implicit none
      real*8 rho,  z, a,  ans, v(300)
      real*8 e
      common /landu1/ e
      real*8 x
      integer nx

      z = 7.25
      a = 14.5


      nx = 0
      x = 1.
      e = 1.e19/1.e9
      do while (x  .lt. 1.2e10)
         rho = x / e
         call zpart( z, a,  rho)
         call totcp(0.d0, 1.d0, ans)
         write(*,*) sngl(x), sngl(ans)
         nx = nx + 1
         v(nx) = log10(ans)
         x = x* 10.**.1
      enddo

      end

 
.1
      enddo

      end

 
