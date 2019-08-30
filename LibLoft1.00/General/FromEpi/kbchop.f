!     ****************************************************************
!     *                                                              *
!     * kbchop:  primitive binary chop to find a root of a function  *
!     *                                                              *
!     ****************************************************************
!
!   /usage/
!            call kbchop(f, x1, x2, er, x, n)
! all must be double precision
!
!     f:  a double precistion function name with one argument
!    x1:  lower boundary of root to be found
!    x2:  upper //
!    er:  specifies relative accuracy of root to be found
!     x:  root found
!     n:  condition code:  no. of iteration needed to find the root
!         (>=0)  or error code (<0).  n = -1 means unconvergence after
!         30 iterations.  x may be errorneous.   n=-2 means f(x1)
!         and f(x2) have same sign.  x becomes undef.
!
!  *** note ***
!         must x1 < x2.   f(x1) and f(x2) must have different sign
!
!
!
!
      subroutine kbchop(f, x1, x2, er, x, n)
      implicit none
      external f
      real*8 f, x1, x2, er, x
      integer n

      character*160 msg

!
!
      real*8 a, b, fa, fb, t, ft

      real*8  td
      integer j
      save t, ft
!
      a = x1
      b = x2
      fa=f(a)
      fb=f(b)
      td = 1.d50
      if(fa * fb .gt. 0.) then
         write(msg,*) ' a, b, fa, fb=',a, b, fa, fb, 'in kbchop'
         call cerrorMsg(msg, 1)
         write(0,*) ' t, ft=',t, ft 
         n=-2
         if( abs(fa) < abs(fb) ) then
            x = a
         else
            x = b
         endif
         return
      else       
         do j = 1, 100
            if(fa*fb .le. 0.) then
               t=(a+b)*.5
               if(abs(td) .gt. er) then
                  if(abs((td-t)/td) .le. er) then
                     n = j
                     goto 100
                  endif
               else
                  if(abs(td - t) .le. er) then
                     n =  j
                     goto  100
                  endif
               endif

               td = t
               ft = f(t)
               if(fa * ft .lt. 0.)then
                  b=t
                  fb = ft
               else   
                  a = t
                  fa = ft
               endif
            else
               n = j
               write(0,*) 'warning;kbchop. at iteration=',n
               write(0,*)'a, fa=', a,fa
               write(0,*)'b, fb=', b,fb
               if(abs(fa) .lt. abs(fb) ) then
                  t = a
               else
                  t = b
               endif
               write(0,*) t, ' is used as solution'
               goto 200
            endif    
         enddo
         n = j
      endif
 100  continue
 200  continue
      x=t
      end
