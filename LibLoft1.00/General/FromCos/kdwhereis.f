!              test kdwhereis
!  
!     
!      implicit none
!c      real*8 xa(5)/-1., 0., 1.33, 1.33, 9./
!      real*8 xa(5)/9.d0, 1.33d0, 1.33d0, 0.d0, -1.d0/
!      real*8 x
!      integer ios, loc
!      do while(1)
!        write(*,*) xa
!        write(*, *) 'enter x'
!        read(*, *, iostat=ios)  x
!        if(ios .ne. 0) stop
!        call kdwhereis(x, 5, xa, 1, loc)
!        if(loc .eq. 0 .or. loc .eq. 5) then
!           write(*, *) " out of range"
!        endif
!        write(*, *) " loc=", loc, " xa =", xa(loc), xa(loc+1) 
!      enddo
!      end
! --------------------------------------------------------------------
!       find location of a given double value in a sorted given double
!     array.
!
!     *********************************
      subroutine kdwhereis(x, in, a, step, loc)
!     *********************************
!      x : real*8. input. given double value.
!     in : integer. input. number of data in a.
!      a : real*8 a(step, in). input array.
!                 a(1, 1), a(1, 2), a(1, 3)... are examined.
!    step: integer. input. see above. give 1 for one dim. array.
!     loc: integer. ouput.  a(1,loc) <= x < a(1,loc+1)  if a is ascending
!                                (if a(1,in) = x, loc= in)
!                           a(1, loc) > x >= a(1, loc+1) if a is dscending
!c                              (if a(1,1) = x, loc = 0)
!                   if loc=0 or loc =in; x is out of range.

      implicit none
      integer in, loc, step
      real*8  x, a(step, in)

      logical  ascending
      integer i1, i2, im
      
      i1 = 0             ! lower and
      i2 = in + 1        ! upper bound
      
      ascending = a(1, in) .gt. a(1, 1) 
      do while (i2 - i1 .gt. 1)
            im = (i1 + i2)/2
            if( ascending .eqv. x .ge. a(1, im) ) then
               i1 = im
            else
               i2 = im
            endif
      enddo       
      loc = i1
      end
!
!     *********************************
      subroutine kwhereis(x, in, a, step, loc)
!     *********************************
!      x : real*4. input. given  value.
!     in : integer. input. number of data in a.
!      a : real*4 a(step, in). input array.
!                 a(1, 1), a(1, 2), a(1, 3)... are examined.
!    step: integer. input. see above. give 1 for one dim. array.
!     loc: integer. ouput.  a(1,loc) <= x < a(1,loc+1)  if a is ascending
!                                (if a(1,in) = x, loc= in)
!                           a(1, loc) > x >= a(1, loc+1) if a is dscending
!c                              (if a(1,1) = x, loc = 0)
!                   if loc=0 or loc =in; x is out of range.

      implicit none
      integer in, loc, step
      real*4  x, a(step, in)

      logical  ascending
      integer i1, i2, im
      
      i1 = 0             ! lower and
      i2 = in + 1        ! upper bound
      
      ascending = a(1, in) .gt. a(1, 1) 
      do while (i2 - i1 .gt. 1)
            im = (i1 + i2)/2
            if( ascending .eqv. x .ge. a(1, im) ) then
               i1 = im
            else
               i2 = im
            endif
      enddo       
      loc = i1
      end
