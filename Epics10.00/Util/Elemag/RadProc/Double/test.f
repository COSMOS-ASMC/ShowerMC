!
      implicit none
      real*8 rho,  z, a,  ans
      real*8 e, cvh2den
      common /landu1/ e
      real*8 v(100), h, x0
      integer ne, nrho
      x0 = 364.
      z = 7.25
      a = 14.5
      h = 10.d3
      h = 100.
      e = 1.e20/1.e9
      do while (h .lt. 100.d3)
         rho = cvh2den(h)
         call zpart(z, a, rho)
         call totcp(0.d0, 1.d0, ans)
         write(*,*)' h=',h/1000., '  km  rho=', rho,
     *  ' ans= ', ans, ' mfp=',1./(ans/x0*rho)/1000.,
     *   '  km'
         h = h * 10.**.01
      enddo
      end

 
