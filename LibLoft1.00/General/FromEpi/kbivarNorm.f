!        test  kbivarNorm
!      real*8 m1, s1, m2, s2, rho
!      integer n
!      real*8 x1, x2
!
!      integer i
!
!      write(0,*) ' enter m1, s1, m2, s2, rho, N'
!      read(*, *) m1, s1, m2, s2, rho,n
!
!      do i = 1, n
!         call kbivarNorm(m1, s1, m2, s2, rho, x1, x2)
!         write(*,*)  sngl(x1),  sngl(x2)
!      enddo
!      end
!
      subroutine kbivarNorm(m1, s1, m2, s2, rho, x1, x2)
!       Bivariate normal distribution;
      real*8  m1   ! input meand of X1
      real*8  s1   ! input variance of X1
      real*8  m2   ! input meand of X2
      real*8  s2   ! input variance of X2
      real*8  rho  ! input correlaion between X1 and X2
      real*8  x1   ! output  random variate for X1 
      real*8  x2   ! output  random variate for X2

      real*8 sqrtrho, rhosave, y1, y2
      save sqrtrho, rhosave
      data  rhosave/31245671d0/

      


      if(abs(rho) .gt. 1.0d0) then
         write(0, *) ' |rho| > 1 for kbivarNorm'
         stop 99999
      elseif(rho .ne. rhosave) then
         rhosave = rho
         sqrtrho = sqrt(1.0d0 -rhosave**2)
      endif
!       first get x1 and x2 with mean 0 and unit variance      
      call kgauss2(0.0d0, 1.d0, y1, y2)      
      
      x1 = y1*s1 + m1

      x2 = ( rhosave*y1 + sqrtrho*y2 ) * s2 + m2

      end

