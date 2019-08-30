      subroutine ksmooth(x,intvx,  y, intvy, n, jin, repeat,
     * cgap, icon)
!        very simple smoothing of data
!         This is intended to smooth arrival time distribution
!         which is integrated from smaller time side  and normalized
!         to be 1 at the largest time.  
!       In many cases, the distribution has long tail approaching to 1.
!       The smooting is tried until the maximum gap of current y values
!       and those after smoothing becomes < cgap/Neff. ( y is updated
!       at each smoothing).  Where Neff
!       is the effective number of data points. Neff is defined as
!       the data number which lie below y=0.95. (This is to avoid the
!       effect of long tail close to y=1.0).
      implicit none
      integer intvx  ! input see below
      integer intvy  ! input see below. Normally these two are  1
      integer n      ! input number of data point
      real*8  x(intvx, n)  ! input. data points are given at
                           ! x(1, i), i=1,...n
      real*8  y(intvy, n)  ! input/output. data to be smoothed at x. 
      integer jin  !  input. 0.  no smoothing for the 1 st and last points
                 !         1.  smoothing is tried also to the 1st and last
                 !         For time distribution, 0 is ok

      integer repeat  !  input.  The maximum number of  smoothing.
                      ! smooting is ceased at this number even cgap/Neff
                      ! condtion is not satisfied.  Normally 500~1000 is ok
!
!               if repeat is large, this procedure gives a straight line.
!
      real*8  cgap  ! input.   critical value to judge the max gap.  As explained
!                      above. 0.03~0.2.
      integer icon  ! output.  number of smoothing actually tried.
                    !          -1 means no success.

!          take 3 data points from smallest x, and get average of
!          3  y's and make it to be the value of 2nd y. 
!          Then, extract next 3 x's until the 3rd point becomes then
!          n-th point.  if jin==1,   

      integer i, j
      real*8 sum, y2, y3, yn1, yn2, r, dy1, dy2, dy3
      real*8 temp
      real*8 gap
      integer neff
      neff = 0
      do i = 1, n
         if(y(1,i) .gt. 0.95 ) exit
         neff = neff+ 1
      enddo
      if(neff .lt. 3) then
         neff= n
      endif

      if( n .le. 2) then
         icon = -1
      else
         icon = 0
         do i = 1, repeat
            icon = icon + 1
            gap = 0.
            if(jin .eq. 1) then
               y2 = y(1, 2)
               y3 = y(1, 3)
               yn1 = y(1, n-1)
               yn2 = y(1, n-2)
            endif

            do j = 2, n-1
               r =(x(1,j)-x(1,j-1))/(x(1,j+1)-x(1,j))
               sum = y(1, j-1) + y(1,j) +
     *          y(1,j+1)*r
               if( j .eq. 2) then
                  temp = sum/(2.0d0 + r)
               else
                  if(gap .lt. abs(temp-y(1,j-1)) ) then
                     gap = abs(temp-y(1,j-1))
                  endif
                  y(1,j-1) = temp
                  temp = sum/(2.0d0 + r)
               endif
            enddo
            y(1,n-1) = temp
            if( jin .eq. 1) then
               dy3= y3-y(1,3)
               dy2 = y2-y(1,2)
               y(1, 1) = (dy3-dy2)/(x(1,3)-x(1,2)) *
     *                   (x(1,1)-x(1,2)) + dy2 + y(1, 1)

               dy1 = yn1 - y(1, n-1)
               dy2 = yn2 - y(1, n-2)
               y(1,n) =(dy2-dy1)/(x(1,n-2)-x(1,n-1))*
     *                 (x(1,n)-x(1,n-1)) + dy1 + y(1,n)
            endif
            if(gap .lt. cgap/neff) exit
         enddo  
      endif
      end
!c        test program
!c         ifort ksmooth.f
!c          ./a.out cgapvalue < data > outdta
! 
!      implicit none
!      integer n, an, i
!      parameter (n=2000)
!      integer repeat, icon
!      real*8 x(n), y(n)
!      real*8 cgap 
!      integer count,  status
!      character*80 buf
!
!
!      count = NARGS()
!      if(count .ne. 2) then
!         write(0,*) ' must give cgap'
!         stop
!      endif
!      call getarg(1, buf, status)
!      write(0,*) status, buf
!      read(buf, *) cgap
!      write(0, *) ' cgap =',cgap
!      an = 0
!      repeat = 1000
!      do while(.true.)
!c         write(0,*) an
!         read(*,*, end=1000) x(an+1), y(an+1)
!         if(an .gt. n) then
!            write(0,*) ' too many data'
!            stop
!         endif
!         an = an +1
!      enddo
! 1000 continue
!      write(0,*) ' data read: n=',an
!
!      call ksmooth(x, 1,  y, 1, an, 0, repeat, cgap, icon)
!      write(0,*) ' icon=',icon
!      do i = 1, an
!         write(*,*) x(i), y(i)
!      enddo
!      end
