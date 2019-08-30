  subroutine kCheckCoPlan(p, n, eps, nn, icon)
!        check coplanarity of given points
    implicit none
    integer,intent(in):: n  ! number of points
    real(8),intent(in):: p(3,n)  ! poitns x,y,z
    real(8),intent(in):: eps
    real(8),intent(out):: nn  ! max diff of normal vecotr scaler prod
    integer,intent(out)::icon ! =0. coplaner
    
    integer i
    real(8):: r1(3), r2(3), v1(3), v2(3), v1s, v2s, temp
    if(n <= 3) then
       icon = 0
       nn = 0.
       return
    endif
    r2(:) = p(:,2) - p(:,1)
    nn = 0.d0
    do i = 1, n-3
       r1(:) = r2(:)
       r2(:) = p(:,i+2) - p(:,1)
       call epvectorProd(r1,r2, v1)
       v1s = dot_product(v1,v1)
       r1(:) = r2(:)
       r2(:) = p(:,i+3) - p(:,1)
       call epvectorProd(r1,r2, v2)
       v2s = dot_product(v2,v2)
       temp = dot_product(v1,v2)
       if( temp /= 0.d0 ) then
       !   write(0,*) 'temp=', temp/sqrt(v1s*v2s)
          nn = max(abs(abs(temp/sqrt(v1s*v2s))-1.0d0), nn)
       endif
    enddo
    if(nn < eps) then
       icon = 0
    else
       icon = -1
    endif
  end subroutine kCheckCoPlan
!
!  program main
!  implicit none
!   integer,parameter::n=8   ! 10^7 -->2.462 s
!   real(8)::r(3, 1000), r2(3,1000)
!   integer icon, i
!   real(8),parameter::pi=asin(1.d0)*2
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
!   real(8)::rm(3,3), rm1(3,3), rm2(3,3), nn
!   real(8)::  eps = 1.d-10
!   real(8)::  xmin, ymin, xmax, ymax, u
!
!  p(3,3)= 1.d-4
!
!  xmin = minval(p(1,:)) - 1.
!  xmax = maxval(p(1,:)) + 1
!  ymin = minval(p(2,:))  -1
!  ymax = maxval(p(2,:)) + 1.
!
!  cost=cos(pi/5.)
!  sint=sin(pi/5.)
!  call cgetRotMat3(1, cost, sint, rm)
!  cost = cos(pi/3.)
!  sint = sin(pi/3.)
!  call cgetRotMat3(3, cost, sint, rm1)
!  call cmultRotMat3(rm, rm1, rm2)
!  do i = 1, n
!     call capplyRot3(rm2, p(1,i), q(1,i))
!     call capplyRot3(rm2, nv,    nv2)
!  enddo
!
!  call epCheckCopla(q, n, eps, nn,  icon)
!  write(0,*) ' nn=',nn, ' icon=',icon
!  write(0,*) '--------------'
!  do i = 1, 1000
!     call rndc(u)
!     r(1,i) = (xmax-xmin)*u + xmin
!     call rndc(u)
!     r(2,i) = (ymax-ymin)*u + ymin
!     r(3,i) = 0.
!     call capplyRot3(rm2, r(1,i), r2(1,i))
!  enddo
!  call kCheckCoplan(r2, 1000, eps, nn,  icon)
!  write(0,*) ' nn=',nn, ' icon=',icon
!end program


