!c       teste ksampAF
!c      (x,y) is wave length vs light yield for  Li2B4O7 
!      integer n
!      parameter (n=33)
!      real*8 x(n), y(n)
!c      data x/300., 320., 340., 360., 380., 390., 400.,
!      data x/ 320., 340., 360., 380., 390., 400.,
!     *       410., 420., 430., 440.,450.,460.,470., 480., 490, 500.,
!     *       510., 520., 530., 540.,550.,560.,570., 580., 590, 600.,
!     *       610., 620., 630., 640., 650., 660., 670./
!c      data y/0.,  0.1, 1.,4., 15.,20.,32.,
!      data y/ 0.0, 1., 4., 15.,20.,32.,
!     *      39., 56., 70., 78., 90., 95., 99., 105.,111.,118.,
!     *      117., 99., 79., 63.,56., 47., 43., 40.,37.,31.,23.,
!     *      20., 17., 16., 10., 5., 0./
!      integer  nc
!      parameter (nc=n+1)
!      real*8 coef(nc, 3)
!
!      real*8 yi(n)
!      real*8 coef2(nc,3)
!      real*8 a, b, dx, xx, u, yy
!
!      real*8 xs
!      integer i
!      call ksampAF0(x, y, n, coef, nc,  yi, coef2)
!        to see how well interpolation was done
!      do i = 1, n
!         write(*,'(a, 3g15.4)') 'o ', x(i), y(i),  yi(i)
!      enddo
!      a= x(1)
!      b= x(n)
!      dx = (b-a)/(5*n)
!      do i = 1, 5*n+1
!        xx = a + (i-1)*dx
!         u =  (i-1)/(5.0*n)
!         call  kcsplIntp(x, y, n, coef, nc,   xx,  xs)
!         call  kcsplIntp(yi, x, n, coef2, nc,  u,  yy)
!         write(*,'(a, 4e15.4)') 'i ',  xx, xs, u, yy
!      enddo
!
!c       generate random number
!      do i = 1, 100000
!         call ksampAF(x, yi, n, coef2, nc, xs)
!         write(*,*) xs
!      enddo
!      end
!
!    sample a random number following an arbitrary function
!    which is expressed as a table of n-point  (x,y)
!    This is initialization routine
!
      subroutine ksampAF0(x, y, n, coef, nc,  yi, total,
     *      coef2)
      implicit none
      integer n  !  input number of data points
      real*8 x(n), y(n) ! input data points (x,y)
!                          zero is assumed at x<x(1) or x>x(n)
      integer nc  ! input. must be >=n-1
      real*8 coef(nc, 3) ! output. to keep the spline coefficents
      real*8 yi(n) ! output. yi(i)= Normalized integral from x=x(1) to x(i)
!                   (i = 1, n)      (yia(1)=0. yi(n)= 1.0) 
      real*8  total  ! integral valeu of the function
!
      real*8 coef2(nc, 3)  !ouptut. coef2 to be used for (yi, x)
!                   spline interploation      
!



      real*8 dx, a, b, sum
      integer i
      real*8 x1, x2
 
      if( nc .lt. n-1 ) then
         write(0,*) ' nc< n-1 in ksamplAF0'
         stop 1234
      endif

      call kcsplCoef(x, y, n, coef, nc)
      a = x(1)
      b = x(n)
      call kcsplInteg(x, y, n, coef, nc, a, b, total)

      yi(1) = 0.
      do i = 2, n
         b = x(i)
         call kcsplInteg(x, y, n, coef, nc, a, b, sum)
         yi(i) = sum/total
      enddo
!          to assure yi(n) = 1.
      yi(n) = 1.
      call kcsplCoef(yi, x, n, coef2, nc)
      end


      subroutine ksampAF0_b(x, y, n, coef, nc)
      implicit none
!              this is only for getting interplation
      integer n  !  input number of data points
      real*8 x(n), y(n) ! input data points (x,y)
!                          zero is assumed at x<x(1) or x>x(n)
      integer nc  ! input. must be >=n-1
      real*8 coef(nc, 3) ! output. to keep the spline coefficents


      real*8 dx, a, b
      integer i
      real*8 x1, x2
 
      if( nc .lt. n-1 ) then
         write(0,*) ' nc< n-1 in ksamplAF0'
         stop 1234
      endif

      call kcsplCoef(x, y, n, coef, nc)
      end

      subroutine ksampAF(x, yi, n, coef2, nc, xs)
      implicit none
      integer n
      real*8 x(n), yi(n)
      integer nc
      real*8 coef2(nc, 3)
      real*8 xs  !output.  sampled x
      real*8 u

      call rndc(u)
      call kcsplIntp(yi, x, n, coef2, nc, u,  xs)
      end
