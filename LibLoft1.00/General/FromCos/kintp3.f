!     ****************************************************************
!     *                                                              *
!     * kintp3:  lagrange's 3-point interpolation                    *
!     *                                                              *
!     ************************ tested 87.03.07 ***************k.k ****
!
!    /usage/   call kintp3(f, intv, n, x1, h, x, ans)
!
!     f:  table of some function values to be interpolated
!  intv:  f is used with step intv
!     n:  # of available f
!    x1:  f containes the function values at x1, x1+h,...x1+(n-1)*h
!     h:  interval of argument
!     x:  argument at which the value of the function is to be obtained
!
!   ans:  interpolated value
!
      subroutine kintp3(f, intv, n, x1, h, x, ans)
      implicit none
!
      integer intv, n
      real*8 f(intv, n), x1, h, x, ans
!
      integer i
      real*8 p, ta, tb
!
      if(n .lt. 3) then
          ans=0.
      else
          i=(x-x1)/h
          if(i .lt. 0) then
              i=0
          elseif(i .gt. 0) then
              if(i+3 .gt. n) then
                 i=n-3
              endif
          endif
          p=(x-x1-h*float(i+1))/h
          ta=p-1.
          tb=p+1.
          ans=0.5*p*(ta*f(1, i+1)+tb*f(1, i+3)) - ta*tb*f(1, i+2)
      endif
      end


