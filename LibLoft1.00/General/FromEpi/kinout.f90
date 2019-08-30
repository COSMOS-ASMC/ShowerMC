!     kxplsl is use.
!          test kinout
!  implicit none
!!      concave case test
!!  integer,parameter::n=11
!!  real(8) xa(n), ya(n)
!!  real(8):: eps, x0, y0
!!  integer icon, i
!!  real(8)::  xmin, ymin, xmax, ymax, u
!!  data xa/0.,     3.0d0, 3.d0, 6.d0, 0.d0, 0.d0, -3.0d0, -9.d0, -2.d0, -7.d0, 0.d0/
!!  data ya/-5.0d0, -2.d0, 1.d0, 3.d0, 3.d0, 8.d0,   3.d0, 5.d0,  -1.d0, -5.5d0, -1.5d0/
!
!!     convex case test
!  integer,parameter::n=6
!  real(8) xa(n), ya(n)
!  real(8):: eps, x0, y0
!  integer icon, i, irsave(2)
!  real(8)::  xmin, ymin, xmax, ymax, u
!  data xa/0.,     3.0d0,  6.d0,  0.d0,  -9.d0,  -7.d0 /
!  data ya/-5.0d0, -2.d0,  3.d0,  8.d0,   5.d0,  -5.5d0/
!
!  eps=1.e-10     !  10^-6 is too large
!  xmin = minval(xa)*1.6
!  xmax = maxval(xa)*1.6
!  ymin = minval(ya)*1.6
!  ymax = maxval(ya)*1.6
!!  irsave=(/527214161,  1512820383/)
!  call rnd1i(irsave)
!  do i = 1, 10000000
!!     call rnd1s(irsave)
!!     write(0,*) irsave
!     call rndc(u)
!     x0 = (xmax-xmin)*u + xmin
!     call rndc(u)
!     y0 = (ymax-ymin)*u + ymin
!     call kinout(xa, 1, ya, 1, n, x0, y0, eps, 'xxxx', icon)
!!     write(*,*) icon, x0, y0
!  end do
!end program
!     ****************************************************************
!     *  kinout: judge whether a given point is in inner or outer    *
!     *          region of a given polygon.                          *
!     ****************************************************************
!
! /usage/ 
!        call kinout(xa, intvx, ya, intvy, n, x0, y0, eps,  convex, icon)
!  ***   kxplsl needed.   ***
!!
  subroutine  kinout(xa, intvx, ya, intvy, n, x0, y0, eps, convex, icon)
    implicit none
    integer,intent(in) :: intvx ! interval of data in xa
    integer,intent(in) :: intvy !  //                 ya
    integer,intent(in) :: n    ! # of points composing the polygon.
                        ! line segment between  first and n-th point is 
                        ! one of the edges.
    real(8),intent(in) :: xa(intvx,n)  !  array containe x values of a polygon
    real(8),intent(in) :: ya(intvy,n)  !  array containe y values of a polygon
    real(8),intent(in) :: x0, y0  ! given point
    real(8),intent(in) :: eps  !  accuracy to judge //ity or coincident lines.
                         ! ~ 10^-10 is ok
    character(*),intent(in) :: convex ! Now  better not to put 'conv'
           !                   even for convex polygon.
           !    'conv'->polygon is assumed to be convex, else it
           !     may be convex or concave.  Even for convex cae,
           !     shown in test, convex="xxxx" is little bit faster
           !     than convec="conv". (10^7 points, 5.23 s vs 5.32 s)
           !    (without any io by 3GHz machine). 
    integer,intent(out) :: icon 
              !   icon:   1 if the point is in outer region
              !          -1                 in inner region
              !           0  on the edge of the polygon
              !          -2 input is invalid (n<=2)
!
! method:  draw a line at x0 and get crossing point y(i) with each
!         edge.  Sort y(i) and see if y0 is inbetween
!         some of y(1)-y(2), y(3)-y(4), y(5)-y(6)...
!         If yes, (x0,y0) is inside. 
!         If y0 coinsides with y(i), the point is on the line.
!  ***    When x0  overlappes with some of the edges,
!            see if (x0,y0) is on the edge.
!            if no, take any of end point y
!     
!
!        if two consecutive points of the polygon are same, one of them
!        will be ignored.  if resultant vertex points become <=2,
!        icon=-2 will result.
!

    real(8):: y(n) ! dynamic allocation
    integer:: idx(n)  !//

    integer:: ny
    complex(kind(0d0))::z0, expa, z1, z2
    real(8)::p, q, u
    integer :: jcon, i, j
!////////////
!    write(0, *) n
!    write(0, *) xa(1,:)
!    write(0, *) ya(1,:)
!    write(0, *) 'x0,y0=', x0, y0
!/////////////
    if(n .le. 2 .or. intvx .le. 0 .or. intvy .le. 0)then
       icon=-2
    else
       ny = 0
!         find xing point of a line, x=x0,  and each edge
       z0=cmplx( x0, y0, 8)
       call rndc(u)
       u = u*6.28318530717958647692d0  ! 0< u <2pi
       expa= exp( (0.d0, 1.d0)*u)   ! exp(i u)
!       expa = cmplx(0.d0, 1.d0, 8)   ! exp(i pi/2)
       do i = 1, n
          z1 = cmplx( xa(1,i), ya(1,i), 8 )
          if( i == n ) then
             z2 = cmplx( xa(1,1), ya(1,1), 8)
          else
             z2 = cmplx( xa(1,i+1), ya(1,i+1), 8)
          endif
          call kxplsl(z0, expa, z1,z2, eps, p, q, jcon) 
!/////////////
!          write(0,*) ' jcon =',jcon, ' for i=',i
!          write(0,*) ' p, q=', p, q
!///////////
            !   p--distance from z1 to the crossing point (singned)
            !      this is the portion of the segment length.
            !   q--distance from z0 to    //
            !   jcon=0  crossing point found on the segment (0 <= p <= 1)
            !        1  the line overlaps the segment.  the position of
            !           z0 may be judged by q.  if q < 0, z0 is 'left'
            !           to z1,  if q > 1, z0 is 'right' to z2
            !           if q=0, z0 is on the segment
            !        2  crossing point is outside of the segment
            !           if p < 0, it is 'left' to z1, if p > 1, it is
            !           'ritht' to z2
            !        3  no crossing point at all, i.e., they are parallel
            !           and distance between two lines > eps
            !           in this case q becomes distance between two lines
          if(jcon == 0 ) then
             if( abs(q) <= eps  ) then
                icon = 0
                return ! ************
             endif
             ny = ny + 1
             y(ny) = q + y0
          elseif ( jcon == 1 )  then
             if( q == 0.d0 ) then
                ! z0 on the segment
                icon = 0
                return ! ************
             endif
             ny = ny + 1
             y(ny) = ya(1,i)
          endif
          if(convex == "conv" .and. ny == 2 ) then
             exit
          endif
       enddo
!/////////
!       write(0,*) ' after loop ny=',ny
!       write(0,*) ' y=',y(1:ny)
!//////////

       if( ny == 2 ) then
!///////////n
!          write(0, *) '(y0 - y(1)) * (y0 - y(2)) = ', (y0 - y(1)) * (y0-y(2)) 
!//////////////
          if( ( y0 - y(1)) * (y0- y(2) )  <= 0.d0 ) then
             icon = -1
          else
             icon = 1
          endif
          return !******************
       elseif( ny == 0 ) then
          icon = 1
       elseif( ny > 1) then
          if( 2*(ny/2) /= ny ) then
             write(0,*) ' error in kinout:ny= ', ny
             stop
          endif
          !  sort y : ascending order
          call kqsortd(y, idx, ny)
          do j = 1, ny-1, 2
             if( (y(idx(j))- y0) * (y(idx(j+1))-y0) <= 0.d0 ) then
                icon = -1
                return !************
             endif
          enddo
          icon = 1
       else
          write(0,*) ' strange ny=',ny, ' in kinout'
          write(0,*) ' x0, y0=',x0, y0
          stop
       endif
    endif
  end subroutine kinout





