!          Trapezoidal integral of a given table data.
      subroutine ktrpzIntegT(t, intv, n, x0, dx, x, ans)
      implicit none
      integer intv  ! input. see below
      integer n     ! input. number of data values
      real*8 t(intv, n)  ! input.  t(1, 1), t(1, 2), .. t(1, n) are used.
                         ! function values at x0, x0+dx, ..
      real*8 x0  ! input  the leftmost value of x.
      real*8 dx  ! input. step of x.
      real*8 x   ! input. integral is from x0 to x. 
      real*8 ans ! output.  integral value
      real*8 nx
      integer i, m

      m = (x - x0)/(0.9999999d0*dx) + 1
      nx = x0 + (m-1) * dx


      ans = 0.
      do i = 2, m-1
         ans = ans + t(1, i)
      enddo

      if(m .gt. 1) then
         ans = (ans + (t(1,1) + t(1, m))/2)* dx
      endif

      if(m .lt. n) then
         ans = ans + 
     *   ((t(1, m+1) - t(1,m))*(x - nx)/dx + t(1,m))*(x-nx)/2
      endif
      end
