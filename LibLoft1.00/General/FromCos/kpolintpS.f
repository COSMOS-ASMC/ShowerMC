!         single precision version
!         test kpolintpS, kpolintpSFE
!
!       implicit none
!       integer n
!       parameter (n = 20)
!       real(4):: xa(n), ya(n), x, y, error
!       integer i,  m
! 
!      do i = 1, n
!          xa(i) = i/3.d0
!          ya(i) = exp(xa(i))
!       enddo
!       m = 5   ! use m points around x.
!       do i =1, n
!          x = xa(i) - 0.2
!          call kpolintpSFE(xa, 1, ya, 1, n, m, x, y, error)
!          write(*, *) x, exp(x), y, error
!       enddo
!       end
      subroutine kpolintpSLogxyFE(xa, xstep, ya, ystep, nt, m,
     *  logxy,  x, y, error)
      implicit none
!        This is a front end for kpolintpS for which we must give
!     some few to several points around x. This manages such 
!     business automatically.  This version takes log of x and or y before
!     kpolintpS is called.  
!
      integer xstep ! input.   see below
      integer nt    ! input.   total number of points 
      integer m     ! input.   the number of points to be used
                    !          for interpolation. must be <=10.
      real(4):: xa(xstep, nt)  ! input. values of x-coordinate at xa(1, i)
                           !        (i=1, nt) are valid x data.
      integer ystep ! input.  see below
      real(4):: ya(ystep, nt)  ! input. values of y-coordinate at ya(1, i)
                           !       (i=1, nt) are valid y data.
      integer logxy   !   input.  1-->log(x) 2-->log(y) 3-->log(x),log(y)
      real(4):: x             ! input. x-value where an interpolated y
                           !        value is wanted
      real(4):: y             ! output. see above

      real(4):: error         ! output. estimated error
! -----------------------------------------------

      real(4):: logx(10), logy(10)      ! working array.
      real(4):: xx, yy

      integer  loc, k, i
      logical  kbitest
!          find location of  x  in xa
      call kwhereis(x, nt, xa, xstep,  loc)
      k = min(max(loc - (m-1)/2,1), nt+1-m) ! max of m points from k

      do  i = k, m+k-1
         if( kbitest(logxy,  1)) then
            if( xa(1, i) .gt. 0.) then
               logx(i-k+1) = log(xa(1,i))
               xx = log(x)
            else
!               if some of x is <= 0; we don't use log 
               call kpolintpSFE(xa, xstep, ya, ystep, nt, m,
     *           x, y, error)
               return           !  **************
            endif
         else
            logx(i-k+1) = xa(1,i)
            xx = x
         endif
         if( kbitest(logxy, 2)) then
            if( ya(1, i) .gt. 0.) then
               logy(i-k+1) = log(ya(1,i))
            else
!               if some of y is <= 0; we don't use log 
               call kpolintpSFE(xa, xstep, ya, ystep, nt, m,
     *           x, y, error)
               return           !  **************
            endif
         else
            logy(i-k+1) = ya(1,i)
         endif
      enddo
      call kpolintpS(logx, 1, logy, 1, m, xx, yy, error)
      if( kbitest(logxy, 2) ) then
         y = exp(yy)
      else
         y = yy
      endif
      end
      subroutine kpolintpSLogFE(xa, xstep, ya, ystep, nt, m,
     *     x, y, error)
      implicit none
!        This is a front end for kpolintpS for which we must give
!     some few to several points around x. This manages such 
!     business automatically.  This version takes log of y before
!     kpolintpS is called.  
!
      integer xstep ! input.   see below
      integer nt    ! input.   total number of points 
      integer m     ! input.   the number of points to be used
                    !          for interpolation. must be <=10.
      real(4):: xa(xstep, nt)  ! input. values of x-coordinate at xa(1, i)
                           !        (i=1, nt) are valid x data.
      integer ystep ! input.  see below
      real(4):: ya(ystep, nt)  ! input. values of y-coordinate at ya(1, i)
                           !       (i=1, nt) are valid y data.
      real(4):: x             ! input. x-value where an interpolated y
                           !        value is wanted
      real(4):: y             ! output. see above
      real(4):: error         ! output. estimated error
! -----------------------------------------------

      real(4):: logy(10)      ! working array.
      real(4):: yy

      integer  loc, k, i
!          find location of  x  in xa
      call kwhereis(x, nt, xa, xstep,  loc)
      k = min(max(loc - (m-1)/2,1), nt+1-m) ! max of m points from k

      do  i = k, m+k-1
         if(ya(1,i) .gt. 0.) then
            logy(i-k+1) = log(ya(1,i))
         else
!            if some of y is <= 0; we don't use log 
            call kpolintpSFE(xa, xstep, ya, ystep, nt, m,
     *           x, y, error)
            return   !  **************
         endif
      enddo
      call kpolintpS(xa(1, k), xstep, logy, 1, m, x, yy, error)
      y = exp(yy)
      end
      subroutine kpolintpSFE(xa, xstep, ya, ystep, nt, m,  x, y, error)
      implicit none
!        This is a front end for kpolintpS for which we must give
!     some few to several points around x. This manages such 
!     business automatically.
!
      integer xstep ! input.   see below
      integer nt    ! input.   total number of points 
      integer m     ! input.   the number of points to be used
                    !          for interpolation. must be <=10.
      real(4):: xa(xstep, nt)  ! input. values of x-coordinate at xa(1, i)
                           !        (i=1, nt) are valid x data.
      integer ystep ! input.  see below
      real(4):: ya(ystep, nt)  ! input. values of y-coordinate at ya(1, i)
                           !       (i=1, nt) are valid y data.
      real(4):: x             ! input. x-value where an interpolated y
                           !        value is wanted
      real(4):: y             ! output. see above
      real(4):: error         ! output. estimated error


      integer  loc, k
!          find location of  x  in xa
      call kwhereis(x, nt, xa, xstep,  loc)
      k = min(max(loc - (m-1)/2,1), nt+1-m) ! max of m points from k
      call kpolintpS(xa(1, k), xstep, ya(1, k), ystep, m, x, y, error)
      end


      subroutine kpolintpS(xa, xstep, ya, ystep, n,  x, y, error)
!   
!      integer   n. input. number of points.   
!      real(4)::   xa(xstep, n). input.
!      real(4)::   ya(ystep, n). input.  function values at xa.
!      real(4)::   x.  input.
!      real(4)::   y.  output.  interepolated functon value at x.
!      real(4)::  error. output. estiamted rough error.
!
      implicit none
      integer n, xstep, ystep
      real(4):: xa(xstep, n), ya(ystep, n), x, y, error

      integer i, maxm, j
      parameter (maxm = 10)
      real(4)::  c(maxm), d(maxm), diff, difft

      integer ns,  m
      real(4):: h0, hp, w, den

      if(n .gt. maxm) then
         write(*, *) ' kpolintpS: use lesser number of points'
         stop
      endif
      ns = 1      
      diff = abs(x - xa(1, 1))
      do i = 1, n
         difft= abs(x - xa(1, i))
         if(difft .le. diff) then
            ns = i
            diff = difft
         endif
         c(i) = ya(1, i)
         d(i) = ya(1, i)
      enddo
      y = ya(1, ns)
      ns = ns-1
      do m = 1, n-1
         do i=1, n-m
            h0 = xa(1,i) -x
            hp = xa(1, i+m) - x
            w = c(i+1) - d(i)
            den = h0- hp
            if(den .ne. 0.) then
               den = w/den
               d(i) = hp*den
               c(i) = h0*den
            else
               write(0,*)  ' error in kpolintpS'
               write(0,*) 'x=',x
               write(0,'(10G12.4)' ) ' xa=', (xa(j,1), j=1, n)
               write(0,'(10G12.4)' ) ' ya=', (ya(j,1), j=1, n)
               stop
            endif
         enddo
         if(2*ns .le. n-m) then
            error = c(ns+1)
         else
            error = d(ns)
            ns = ns-1
         endif
         y = y + error
      enddo
      end

      subroutine kpolintpSeqs(x0, dx, ya, ystep, n,  x, y, error)
!   
!      integer   n. input. number of points.   
!      real(4)::   x0. input. x0, x0+dx, x0+2dx, ...x0+(n-1)dx
!                          are given data points
!      real(4)::   ya(ystep,n)  input.  function values at x0,..
!      real(4)::   x.  input.
!      real(4)::   y.  output.  interepolated functon value at x.
!      real(4)::  error. output. estiamted rough error.
!
      implicit none
      integer n,  ystep
      real(4):: x0,  ya(ystep, n), x, y, error, dx

      integer i, maxm
      parameter (maxm = 10)
      real(4)::  c(maxm), d(maxm), diff, difft

      integer ns,  m
      real(4):: h0, hp, w, den
      integer p, q
      real(4):: xa
      xa(p, q) = (q-1)*dx + x0

      if(n .gt. maxm) then
         write(*, *)
     *    ' kpolintpSeqs: use lesser number of points'
         stop
      endif
      ns = 1      
      diff = abs(x - xa(1, 1))
      do i = 1, n
         difft= abs(x - xa(1, i))
         if(difft .le. diff) then
            ns = i
            diff = difft
         endif
         c(i) = ya(1, i)
         d(i) = ya(1, i)
      enddo
      y = ya(1, ns)
      ns = ns-1
      do m = 1, n-1
         do i=1, n-m
            h0 = xa(1,i) -x
            hp = xa(1, i+m) - x
            w = c(i+1) - d(i)
            den = h0- hp
            if(den .eq. 0.) then
               write(0,*) ' error in kpolintpSeqs'
               stop
            endif
            den = w/den
            d(i) = hp*den
            c(i) = h0*den
         enddo
         if(2*ns .le. n-m) then
            error = c(ns+1)
         else
            error = d(ns)
            ns = ns-1
         endif
         y = y + error
      enddo
      end
