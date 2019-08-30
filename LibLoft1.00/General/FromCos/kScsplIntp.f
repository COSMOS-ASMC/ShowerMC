      subroutine kScsplIntp(x, y, n, coef, nc, v,  f)
      implicit none
!        interpolation by the cubic spline.
      integer n    ! input. size of x, y
      real(4):: x(n)  ! input.
      real(4):: y(n)  ! input.
      integer  nc  ! input. size of coef(nc,3)
      real(4):: coef(nc,3)  ! input. which is an output from
                         !     kcsplCoef
      real(4):: v    ! input. argument where the interpolated
                  !        value is requested.
      real(4):: f    ! output. obtained interpolated values
!
      integer i
      real(4):: v1

      call kwhereis(v, n, x, 1, i)
      if(i .eq. 0) then 
         i = 1
      elseif(i .eq. n) then
         i = n - 1
      endif
      v1 = v - x(i)
      f = y(i) + 
     *    v1 * (coef(i,1) + v1 * (coef(i,2) + v1 * coef(i,3)))
      end
