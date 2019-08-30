      subroutine kcsplDif(x, n, coef, nc, v, d1, d2)
      implicit none
      integer n  ! input size of x
      real*8 x(n)  ! input
      integer nc   ! input. size of coef(nc, 3)
      real*8 coef(nc, 3) ! input.  which is output from
                         !       kcsplCoeff
      real*8 v           ! input. value where the interpolation is 
                         !      to be made to obtain the numerical
                         !      differenctiation value of y
      real*8 d1          ! output. 1st order differentiation
      real*8 d2          ! output. 2nd ...


      integer i
      real*8 v1, t

      call kdwhereis(v, n, x, 1, i)
      if(i .eq. 0) then
         i = 1
      endif
      v1 = v- x(i)
      t = coef(i,2) + 3.0d0 * coef(i,3) * v1
      d1 = coef(i,1) + (t + coef(i,2)) * v1
      d2 = t + t
      end

