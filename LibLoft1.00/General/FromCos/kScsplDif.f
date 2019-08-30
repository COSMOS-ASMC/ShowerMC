      subroutine kScsplDif(x, n, coef, nc, v, d1, d2)
      implicit none
      integer n  ! input size of x
      real(4):: x(n)  ! input
      integer nc   ! input. size of coef(nc, 3)
      real(4):: coef(nc, 3) ! input.  which is output from
                         !       kcsplCoeff
      real(4):: v           ! input. value where the interpolation is 
                         !      to be made to obtain the numerical
                         !      differenctiation value of y
      real(4):: d1          ! output. 1st order differentiation
      real(4):: d2          ! output. 2nd ...


      integer i
      real(4):: v1, t

      call kwhereis(v, n, x, 1, i)
      if(i .eq. 0) then
         i = 1
      endif
      v1 = v- x(i)
      t = coef(i,2) + 3.0 * coef(i,3) * v1
      d1 = coef(i,1) + (t + coef(i,2)) * v1
      d2 = t + t
      end
