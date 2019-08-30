      implicit none
!      next 4 lines  for fortran 90 dynamic memory allocation
#include "Z90histc.f"
#include "Z90hist.f"
      type(histogram1)  h
      type(histogram1)  k
!      next 2 is for static old forran
!      include 'Zhist.f'
!      type(histogram):: h, k
      save h
      real*8 x
      integer i
      
      call kwhisti(h,  0., 0.1, 300,  b'0010')
      call kwhistc(h)
      call kwhisti(k, 1.5, 0.1, 30, b'1111' )
      call kwhistc(k)

      do i = 1, 1000000
         call kgauss(10.0d0, 1.d0, x)
         call kwhist(h, sngl(x), 1.0 )
         call rndc(x)
         x = x**(-0.5)
         call kwhist( k, sngl(x), 1.0 )
!         write(*,*) sngl(x), 1.0
      enddo
      
      call kwhists(h, 1000000.0 )
      call kwhistp(h, 'a' )
      write(*,*) 
      call kwhists(k, 1000000.0)
      call kwhistp(k, 'b' )
      end

