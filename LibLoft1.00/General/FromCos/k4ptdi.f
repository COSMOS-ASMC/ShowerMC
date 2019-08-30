!     ****************************************************************
!     *                                                              *
!     * k4ptdi:  4-point two dimensional interpolation               *
!     *                                                              *
!     ****************************************************************
!
!   /usage/
!
!       call
!              k4ptdi(f, im, jm, iadj,  x0, y0, hx, hy, x, y, ans)
!
!     f is a 2-dimensional table of some function with 2 arguments.
!     f containes the function values at (x0,y0), (x0+hx, y0+hy),...
!     (x0+(im-1)*hx, y0+(jm-1)*hy).
!     iadj is the adjustable dimension.
!     ans gets the value of the funtion at (x,y)
!
!
!
      subroutine k4ptdi(f, im, jm, iadj, x0, y0, hx, hy, x, y, ans)
      implicit none
      integer im, jm, iadj
      real*8 f(iadj,jm),  x0, y0, hx, hy, x, y, ans
!
      integer i, j
      real*8 a, b, p, q, p1, q1

      a=(x-x0)/hx
      b=(y-y0)/hy
      i=a
      j=b
      i=min(max(i,0)+1,im-1)
      j=min(max(j,0)+1,jm-1)
      p=a+1.-i
      q=b+1.-j
      p1=1.-p
      q1=1.-q
      ans=( f(i,j)*p1 + f(i+1,j)*p ) * q1 +
     *                  ( f(i,j+1)*p1 + f(i+1,j+1)*p ) * q
      end
