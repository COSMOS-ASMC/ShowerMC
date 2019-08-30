!      external ff
!      real*8 ff
!      real*8 x1, x2, x, eps, ans
!      integer icon
!      x1 = 0.
!      x2 = 1.8
!      x = 1.
!      eps= 1.d-11
!      call kbinChop(ff, x1, x2, x, eps, ans, icon)
!      write(*, *) icon, ans, ff(ans)
!      end
!      real*8 function ff(x)
!      real*8 x
!      ff = sin(x) - 0.5d0
!      end
!      
!      Binary Chop for getting a solution of  f(x) = 0.
!
      subroutine kbinChop(f, x1, x2, x, eps, ans, icon)
      implicit none
!
      real*8 f  ! input. function name. to be used as f(x)
                !   f(x) = 0 is solved.
      real*8 x1 ! input. lower bound of solution
      real*8 x2  ! input.  upper bound of solution
      real*8 x   ! input. initial guess of  solution.  
      real*8 eps  ! input. relative error of solution.
      real*8 ans  ! output. obtained solution
      integer  icon ! output. condition code. 0--> ok.
!               1--> unconvergence after 45 iterations
!               2--> x not in the range
      real*8 xa, xb, fa, fb, xt, ft
      integer  n
!
      if(x .lt. x1 .or. x .gt. x2) then
         icon = 2
      else
         xa = x1
         xb = x2
         fa = f(xa)
         fb = f(xb)
         icon = 1
         do n = 0, 45
            if(fa * fb .gt. 0.) then
               icon = 1
               goto 100
            else
               xt = (xa + xb)/ 2
               ft = f(xt)
               if( ft * fa .gt. 0.) then
                  xa = xt
                  fa = ft
               else
                  xb = xt
                  fb = ft
               endif
               if(abs(xt) .gt. 1.) then
                  if(abs( (xa-xb) / xt ) .lt. eps) then
                     icon = 0
                     goto 100
                  endif
               else
                  if(abs(xa-xb) .lt. eps) then
                     icon = 0
                     goto 100
                  endif
               endif
            endif
         enddo
 100     continue
         ans = xt
      endif
      end


                     
