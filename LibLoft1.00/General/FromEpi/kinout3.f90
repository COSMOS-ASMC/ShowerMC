!!!          test kinout3
!     ****************************************************************
!     *  kinout3: judge whether a given point is in inner or outer    *
!     *          region of a given polygon.                          *
!     * The polygon is defined in 3D space but the points are coplaner
!     ****************************************************************
!
!
  subroutine  kinout3(p, n,  r, nv, sumteta, cond)
!    judge if a given point on a plane formed by a polygon is inside
!    of the polygon or not.
    implicit none
    integer,intent(in) :: n    ! # of points composing the polygon.
                        ! line segment between  first and n-th point is 
                        ! one of the edges.
    real(8),intent(in) :: p(3,n)  !  array containing x,y,z values of a polygon
    ! Points p(:,i), i=1,2,...are globary given so  that the normal of
    ! the surface is directed to the moving direction of  right-handed screw 
    real(8),intent(in) :: r(3)    ! given point
    real(8),intent(in) :: nv(3)   ! normal vector of the surface
    real(8),intent(out):: sumteta  ! if the point is inside, and
           !  if normal vector is proper direction,  sumteta is 2pi 
           !  if opposite direction, -2pi. 
           !  If the point is outside 0 (normally +/- 10^-15 ~10^-10) 
           !  Proper direction means right-handed srew rule.
    character(*),intent(out) :: cond
           ! at least length is 4.  
           !  'in'   the point is in the inner region
           ! 'out'   //           in the outer region
           ! 'vtx'   the point x is on a vertex point
           ! 'edge'  on the edge of the polygon (vertex is not included)
           ! 'ng'   input is invalid (n<=2)

    real(8),parameter::pi=3.14159265358979323846d0 ! asin(1.d0)*2
    real(8),parameter::twopi=2*pi
    real(8)::  xmin, xmax,  ymin, ymax, zmin, zmax

! method: make a triangle by connecting a given point and the
!         verteces of the  polygon.  Then, compute the angle
!         between the two verteces at the given point. 
!         The sum of the angle is 2 pi if the given 
!         point is inside otherwise, 0.

    real(8),parameter:: eps=1.d-12

    real(8):: teta, absr1, absr2, para
    real(8):: cost
    real(8):: r1(3), r2(3),  r1xr2(3)
    integer:: i
    if(n <= 2  ) then
       write(0,*) 'invalid input to kinout3'
       cond='ng'
    else
       xmin = minval(p(1,:)) 
       xmax = maxval(p(1,:)) 
       ymin = minval(p(2,:)) 
       ymax = maxval(p(2,:)) 
       zmin = minval(p(3,:)) 
       zmax = maxval(p(3,:)) 
       sumteta = 0.
       if( r(1) < xmin .or. r(1) > xmax) then
          cond='out'
          return
       elseif(r(2) < ymin .or. r(2) > ymax) then
          cond='out'
          return
       elseif(r(3) < zmin .or. r(3) > zmax) then
          cond='out'
          return
       endif

       do i = 1, n
          r1(:) = p(:,i)  - r(:)
          if(i == n ) then
             r2(:) = p(:,1) - r(:)
          else
             r2(:) = p(:,i+1) - r(:)
          endif
          absr1 = sqrt(  dot_product(r1(:), r1(:)) )
          if(absr1 == 0.d0 ) then
             cond='vtx'
             return ! ************************
          endif
          absr2 = sqrt(  dot_product(r2(:), r2(:)) )
          if( absr2 == 0.d0 ) then
             cond='vtx'
             return ! ************************
          endif
          r1(:) = r1(:)/absr1
          r2(:) = r2(:)/absr2
          cost = dot_product(r1(:), r2(:))
          call epvectorProd(r1, r2, r1xr2)
          para = dot_product(r1xr2, r1xr2)
          if( abs(para) < eps )  then
             if( cost < 0. ) then
                cond="edge"
                sumteta= twopi
                return
             endif
          endif
          cost=max( min( cost, 1.d0), -1.d0)
          teta = acos(cost)
          if( dot_product( r1xr2, nv)  < 0. ) then
             teta = -teta
          endif
          sumteta = sumteta + teta
       enddo
    endif
!       sumteta is 2pi if normal vector is correct direction
!     if opposite, -2pi
!    write(0,*) ' sum=',sumteta

    if( abs(sumteta) > pi) then
       cond='in'
    else
       cond='out'
    endif
  end subroutine kinout3
!
!  program main
!  implicit none
!
!!   ifort kinout3.f90 -L$COSMOSTOP/lib/PCLinuxIFC64 -lcosmos -L$EPICSTOP/lib/PCLinuxIFC64 -lepics
!!!!      concave case test
!   real(8),parameter:: digitize = 10.  ! make 0 if not wanted.
!   logical,parameter:: rotate=.true.
!   integer,parameter::n=8   !! 10^7 -->2.462 s
!!   real(8):: r(3)=(/-1.5, 1.5, 0./), r2(3)
!   real(8):: r(3)=(/1.0d0, 0.999999999999d0, 0.d0/), r2(3)
!   integer icon, i
!
!   real(8),parameter::pi=asin(1.d0)*2
!   real(8):: sumteta
!   
!   character(4)::cond
!   real(8):: cost, sint, q(3,n)
!   real(8)::nv(3) =(/0., 0., -1./) , nv2(3) 
!   real(8):: p(3,n)=(/  0.,   0.,  0.,  &
!                        1.,   1.,  0.,  &
!                        2.,   1.,  0.,  &
!                        2.,   4.,  0.,  &
!                        1.,   2.,  0.,  &
!                        1.,   5.,  0.,  &
!                        0.,   5.,  0.,  &
!                        0.5,  1.,  0.  /)
!   real(8)::rm(3,3), rm1(3,3), rm2(3,3)
!   real(8):: xmin, ymin, xmax, ymax, u
!  xmin = minval(p(1,:)) - 0.3
!  xmax = maxval(p(1,:)) + 0.3
!  ymin = minval(p(2,:)) - 0.3
!  ymax = maxval(p(2,:)) + 0.3
!  if(rotate) then
!     cost=cos(pi/4.)
!     sint=sin(pi/4.)
!     call cgetRotMat3(1, cost, sint, rm)
!     cost = cos(pi/3.)
!     sint = sin(pi/3.)
!     call cgetRotMat3(3, cost, sint, rm1)
!     call cmultRotMat3(rm, rm1, rm2)
!     cost = cos(-pi/5.)
!     sint = sin(-pi/5.)
!     call cgetRotMat3(2, cost, sint, rm1)
!     call cmultRotMat3(rm2, rm1, rm)
!
!     call capplyRot3(rm, nv,    nv2)
!     do i = 1, n
!       call capplyRot3(rm, p(1,i), q(1,i))
!     enddo
!   endif
!
!   if(rotate) then
!     p(:,:) = q(:,:)
!     nv(:) = nv2(:)
!   endif
!
!   do i = 1, n
!      write(0,*) p(:,i)
!   enddo
!   write(0,*) p(:,1)
!
!   do i = 1, 1000000
!      call rndc(u)
!      r(1) = (xmax-xmin)*u + xmin
!      call rndc(u)
!      r(2) = (ymax-ymin)*u + ymin
!      r(3) = 0.
!      if(digitize > 0) then
!         r(:) = int(digitize*r(:))/digitize
!      endif
!      if(rotate) then
!         call capplyRot3(rm, r, r2)
!         r(:) = r2(:)
!      endif
!      call kinout3(p, n, r, nv, sumteta, cond)
!      write(*,'(a,1p 4g16.7)') cond, r(:), sumteta
!   enddo
! end program main
