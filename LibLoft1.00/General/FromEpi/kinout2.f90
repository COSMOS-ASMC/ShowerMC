!!          test kinout2
!  implicit none
!!!      concave case test
!   integer,parameter::n=11   ! 10^7 -->2.462 s
!   real(8) xa(n), ya(n)
!   real(8):: eps, x0, y0
!   integer icon, i
!   real(8)::  xmin, ymin, xmax, ymax, u
!   data xa/0.,     3.0d0, 3.d0, 6.d0, 0.d0, 0.d0, -3.0d0, -9.d0, -2.d0, -7.d0, 0.d0/
!   data ya/-5.0d0, -2.d0, 1.d0, 3.d0, 3.d0, 8.d0,   3.d0, 5.d0,  -1.d0, -5.5d0, -1.5d0/
!!
!!!     convex case test:  10^7 test 1.792 sec
!!  integer,parameter::n=6
!!  real(8) xa(n), ya(n)
!!!  real(8):: eps, x0, y0
!!  real(8)::  x0, y0
!!  integer icon, i, irsave(2)
!!  real(8)::  xmin, ymin, xmax, ymax, u
!!  data xa/0.,     3.0d0,  6.d0,  0.d0,  -9.d0,  -7.d0 /
!!  data ya/-5.0d0, -2.d0,  3.d0,  8.d0,   5.d0,  -5.5d0/
!!
!!  eps=1.e-10     !  10^-6 is too large
!  xmin = minval(xa)*1.6
!  xmax = maxval(xa)*1.6
!  ymin = minval(ya)*1.6
!  ymax = maxval(ya)*1.6
!!  irsave=(/527214161,  1512820383/)
!!  call rnd1i(irsave)
!  do i = 1, 10000000
!!     call rnd1s(irsave)
!!     write(0,*) irsave
!     call rndc(u)
!     x0 = (xmax-xmin)*u + xmin
!     call rndc(u)
!     y0 = (ymax-ymin)*u + ymin
!!     x0 = 0
!!     y0 = 0
!!     write(0,*) 'enter x0,y0'
!!     read(*,*) x0,y0
!     call kinout2(xa, 1, ya, 1, n, x0, y0,  icon)
!!     write(*,*) icon, x0, y0
!  end do
!end program
!     ****************************************************************
!     *  kinout2: judge whether a given point is in inner or outer    *
!     *          region of a given polygon.                          *
!     ****************************************************************
!
! /usage/ 
!  call kinout2(xa, intvx, ya, intvy, n, x0, y0,  icon)
!!
  subroutine  kinout2(xa, intvx, ya, intvy, n, x0, y0,  inside)
    implicit none
    integer,intent(in) :: intvx ! interval of data in xa
    integer,intent(in) :: intvy !  //                 ya
    integer,intent(in) :: n    ! # of points composing the polygon.
                        ! line segment between  first and n-th point is 
                        ! one of the edges.
    real(8),intent(in) :: xa(intvx,n)  !  array containe x values of a polygon
    real(8),intent(in) :: ya(intvy,n)  !  array containe y values of a polygon
    real(8),intent(in) :: x0, y0  ! given point
!    real(8),intent(in) :: eps  !  accuracy to judge //ity or coincident lines.
                         ! ~ 10^-10 is ok
    integer,intent(out) :: inside
              !  insdie:   1 if the point is in outer region
              !          -1                 in inner region
              !           0  on the edge of the polygon
              !          -2 input is invalid (n<=2)
!
! method: make a triangle with a given point  and 
!         each edge, and comput angle at the given point
!         The sum of the angle is 2 pi if the given 
!         point is inside otherwise, 0.

    real(8):: sint, teta, sumteta, absr1, absr2, r1xr2
    real(8):: cost
    real(8):: r1(2), r2(2)
    integer:: i
!////////////
!    real(8)::pi
!    pi=asin(1.d0)*2
!    write(0, *) n
!    write(0, *) xa(1,:)
!    write(0, *) ya(1,:)
!    write(0, *) 'x0,y0=', x0, y0
!/////////////
    if(n .le. 2 .or. intvx .le. 0 .or. intvy .le. 0)then
       inside = -2
    else
       sumteta = 0.
       do i = 1, n
          r1(1) =  xa(1,i) - x0
          r1(2) =  ya(1,i) - y0
          if(i == n ) then
             r2(1) = xa(1,1) - x0
             r2(2) = ya(1,1) - y0
          else
             r2(1) = xa(1,i+1) - x0
             r2(2) = ya(1,i+1) - y0
          endif
          absr1 = sqrt(  dot_product(r1(:), r1(:)) )
          if(absr1 == 0.d0 ) then
             ! (x0,y0) is p1
             inside = 0 
             return ! ************************
          endif
          absr2 = sqrt(  dot_product(r2(:), r2(:)) )
          if( absr2 == 0.d0 ) then
             inside = 0
             return ! ************************
          endif

          r1xr2 = r1(1)*r2(2) - r1(2)*r2(1)
          sint = r1xr2/(absr1*absr2)
          cost = dot_product( r1(:), r2(:))/absr1/absr2
!          teta = asin(sint)
          teta = atan2(sint,cost)
          sumteta = sumteta + teta
!///////////
!          write(0,*) 'i=',i
!          write(0,*) 'absr1,2=',absr1,  absr2
!          write(0,*) 'r1xr2=', r1xr2
!          write(0,*) 'sint=',sint
!          write(0,*) 'teta=',teta, teta*180./pi
!          write(0,*) 'sumteta=',sumteta*180./pi
!          write(0,*) ' '
!////////////
       enddo
    endif
    if( abs(sumteta) > 1.) then
       inside = -1
    else
       inside = 1
    endif
!    write(*,'(1p,3g16.7)') x0, y0, sumteta
  end subroutine kinout2





