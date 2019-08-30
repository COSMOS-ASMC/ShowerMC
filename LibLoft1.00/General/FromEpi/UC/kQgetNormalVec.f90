subroutine  kQgetNormalVec(p, n, req, nv, k, icon)
    implicit none
    integer,intent(in)::n  ! # of points n>=3
    real(16),intent(in)::p(3,n) !x,y,z of coplanaer points
          ! complanarity is not checked. To check, use kCheckCoplan
    integer,intent(in):: req ! 0:  nv does  not  not consider surface side
                           !   1:  nv should be  directed to the
         ! right-handed screw moving direction, when points are ordered
         ! globally as screw rotaion.
    real(16),intent(out):: nv(3) ! normal vector
    real(16),intent(out):: k  ! origin to the plane distace
    integer,intent(out):: icon ! 0. nv, k obtained
                             ! -1   p is colinear or n<3.
    real(16),parameter:: pi=3.14159265358979323846264338327950q0
    integer:: i, isave, inside
    real(16):: r1(3), r2(3), r1xr2(3), r1abs, r2abs, temp, r(3)
    real(16):: cost, costmin, sumteta

    character(5):: condi

    if( n< 3 ) then
       icon = -1
       return
    endif
    r1(:) = p(:,2) - p(:,1)
    r1abs = dot_product(r1,r1)
    costmin = 1.0
    do i = 3, n
       r2(:) = p(:,i) - p(:,1)
       r2abs = dot_product(r2, r2)
       if(r1abs > 0 .and. r2abs> 0 ) then
          temp =  sqrt(r1abs*r2abs)
          cost = dot_product(r1,r2)/temp
          if( abs(cost) < 0.99619q0 ) then
             !  opening angle is > 5 deg. so stable nv will be obtained
             call epQvectorProd(r1, r2, r1xr2)
             !  normalize
             nv(:) = r1xr2(:)/sqrt( dot_product(r1xr2,r1xr2) )
             goto 100
          elseif( costmin > abs(cost) ) then
             costmin = abs(cost)
             isave = i
          endif
       endif
    enddo
    r2(:) = p(:,isave) - p(:,1)
    r2abs = dot_product(r2, r2)
    temp =  sqrt(r1abs*r2abs)
    call epQvectorProd(r1, r2, r1xr2)
    nv(:) = r1xr2(:)/sqrt( dot_product(r1xr2,r1xr2) )
100 continue
    if( req == 1 ) then
       ! find an inner point. This algorism seems not perfect 
       ! though very difficult to think exception
       do i = 1, n-1
          if( i /= n-1 ) then
             r(:) = (p(:,i+2) + p(:,i))*0.5q0
          else
             r(:) =(p(:,1) + p(:,i))*0.5q0
          endif
          ! see if r is inside
          call  kQinout3(p, n, r, nv, sumteta, condi)
          if( condi == 'in' ) then
             if( sumteta < pi ) then
                nv(:) = -nv(:)
             endif
             exit
          endif
       enddo
    endif
    k = dot_product(nv(:), p(:,1))
    icon = 0 
  end subroutine kQgetNormalVec
!  program main
!!!          test kgetNormalVec
!  implicit none
!!!!!      concave case test
!
!   integer icon, i
!   real(8),parameter::pi=asin(1.d0)*2
!   real(8):: sumteta, k
!   
!   real(8)::nv(3) =(/0., 0., 1./) , nv2(3) 
!!   integer,parameter::n=8   
!!   real(8):: p(3,n)=(/  0.9,   1.00001,  0.,  &
!!                        1.,   1.,  0.,  &
!!                        2.,   1.,  0.,  &
!!                        2.,   4.,  0.,  &
!!                        1.,   2.,  0.,  &
!!                        1.,   5.,  0.,  &
!!                        0.,   5.,  0.,  &
!!                        0.5,  1.,  0.  /)
!   integer,parameter::n=9   
!   real(8):: p(3,n) =(/ 0., 0., 0., &
!                       2.,  2., 0., &
!                       3.,  0., 0., &
!                       2.5, 3., 0.,  &
!                       4.,  3.2, 0., &
!                       3.0, 3.3, 0., &
!                       1.5, 3.2, 0., &
!                       2.4, 2.4, 0., &
!                       1.3, 1.9, 0. &
!                       /)
!   real(8):: cost, sint, q(3,n)
!   real(8)::rm(3,3), rm1(3,3), rm2(3,3)
!  cost=cos(pi/4.)
!  sint=sin(pi/4.)
!  call cgetRotMat3(1, cost, sint, rm)
!  cost = cos(pi/3.)
!  sint = sin(pi/3.)
!  call cgetRotMat3(3, cost, sint, rm1)
!  call cmultRotMat3(rm, rm1, rm2)
!  cost = cos(-pi/5.)
!  sint = sin(-pi/5.)
!  call cgetRotMat3(2, cost, sint, rm1)
!  call cmultRotMat3(rm2, rm1, rm)
!  do i = 1, n
!     call capplyRot3(rm, p(1,i), q(1,i))
!     call capplyRot3(rm, nv,    nv2)
!  enddo
!  call  kgetNormalVec(q, n, 1, nv,  k, icon)
!!  call  kgetNormalVec(p, n, 1, nv2,  k, icon)
!  write(0, *) ' icon =', icon
!  write(0,*) '  nv=',nv(:)
!  write(0,*) ' nv2=',nv2(:)
!  write(0,*) 'k=', k
!end program main
