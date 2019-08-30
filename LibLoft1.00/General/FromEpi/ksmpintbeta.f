!      implicit none
!      real*8 x
!      integer i, m,n
!      call cerrorMsg('Enter m, n', 1)
!      read(*, *) m, n
!      do i = 1, 100000
!         call ksmpintbeta(m, n, x)
!         write(*, *) sngl(x)
!      enddo
!      end

      subroutine ksmpintbeta(m, n, x)
      implicit none
!         samples random variable x of which density is
!          x**m * (1-x)** n dx 
!      where m >= 0 and n>=0 are integers.
!
!        This can be sampled generating n+m+1 uniform 
!     random variables in 0~1 and take the n+1-th largest one
!
      integer m  ! input. >=0
      integer n  ! input. >=0.  m+n must be =< 100
      real*8 x   ! output. sampled variable
!
      integer maxsize
      parameter (maxsize=100)
      real*8 ua(maxsize), u
      integer i

      if( m+n .gt. maxsize) then
         write(*,*)' ksmpintbeta:  too large m+n; m=',m, ' n=',n
         stop 9999
      endif
      if(m .lt. 0  .or.  n .lt. 0) then
         write(*, *) 
     *   ' ksmpintbeta: too small m or n; m=',m, ' n=',n
         stop 9999
      endif
      do i = 1, m+n+1
         call rndc(u)
         ua(i) = u
      enddo
!           sort
      call kcsr1(ua, m+n+1, 'd')
      x = ua(n+1)
      end

      
