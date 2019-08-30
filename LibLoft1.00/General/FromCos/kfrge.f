!c        test kfrge.
!      program test_kfrge
!      implicit none
!      integer l, icon
!c
!      real*8 x(6)/4., 3., -1., 9., 2., 4./
!      call kfrge(x, 1, 6, 3.d0, l, icon)
!      write(*, *) l, icon, x(l)
!      end
!     ****************************************************************
!     *                                                              *
!     *  kfrge: find real data (position) .ge. given value           *
!     *                                                              *
!     ***********************  tested 87.06.07  **********************
!
!   /usage/
!          call kfrge(x, intvx, n, c, l, icon)
!
!     x:  real*8. input.  data array
! intvx:  integer. input. interval of data in x
!     n:  integer. input. |n| is no. of data in x
!     c:  real*8. input.  given value.  x  .ge.  c is sought for.
!     l:  integer. output. position of xf in x.  xf=x(1,m)
!  icon:  0 if found else 1
!
! *** note ***
!         if n>0 search is made for from 1st, else from last
!         if icon=1 resluts, m will be n+1 or 0 dependingon n>0 or
!         n<0.
!
!
      subroutine kfrge(x, intvx, n, c, m, icon)
      implicit none
      integer intvx, n, m, icon
      real*8  x(intvx, *), c
!
      integer i
!
      if( n .gt. 0 ) then
          do i=1, n
             if(x(1,i) .ge. c) then
                 m=i
                 icon=0
                 goto 100
             endif
          enddo   
          icon=1
          m=n+1
      elseif(n .lt. 0) then
          do i=-n, 1, -1
              if(x(1,i) .ge. c) then
                  m=i
                  icon=0
                  goto 100
              endif
          enddo    
          m=0
          icon=1
      else
          icon=1
      endif
  100 continue
      end
