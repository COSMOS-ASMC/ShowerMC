!  real(8):: ps(3), p1(3), p2(3), dir(3), h, el, teta, fai
!  real(8),parameter:: pi = asin(1.d0)*2
!  real(8),parameter:: torad=pi/180.d0
!  integer:: icon
!  ps = (/0., 0., 1./)
!  p1 = (/5., 1., 0.0/)
!  p2 = (/1., 2., 0.0/)
!  teta =45.d0 *torad
!  fai = 35.d0 *torad 
!  h = 10.d0
!  dir = (/sin(teta)*cos(fai), sin(teta)*sin(fai), cos(teta)/)
!  write(0,*) ' dir=', dir
!  call kxplvsq(ps, dir, p1, p2, h, el, icon)
!  write(0,*) icon, el
!  if(icon == 0) then
!     write(0,*) ps(:)+el*dir(:)
!  endif
!end program
!
  subroutine kxplvsq(ps, dir, p1, p2, h, el, icon)
!   get crossing point of a line and vertical square
!      (z-directin is vertical)

    implicit none
    real(8),intent(in)::ps(3)  ! a line passes this point
    real(8),intent(in)::dir(3) ! it's direction cosines
    real(8),intent(in)::p1(3) ! point of the vertical square
    real(8),intent(in)::p2(3) ! another // the square.  p1(3)=p2(3)
    real(8),intent(in)::h ! height of the square
    real(8),intent(out):: el ! length to the crossing point with
                       ! the square plane el > 0. from ps.
    integer,intent(out):: icon ! 0: xp obtained 
                               ! -1: xp not exists 
    ! ==============
!      first consider by projection onto the x-y plane
!     
    complex(kind(0d0))::expa, z1, z2, z0
    real(8)::cost, sint, p, q, xp(3)
    real(8),parameter::eps=1.d-10, eps2=1.d-6
    integer::jcon

    cost = dir(3)
    if( abs(cost) >= 1.d0) then
       ! no x-point
       icon = -1
    else

       sint=sqrt(1.d0 - cost**2)
       expa = cmplx( dir(1), dir(2), 8)/sint

       z0 = cmplx(ps(1), ps(2), 8)
       z1 = cmplx(p1(1), p1(2), 8)
       z2 = cmplx(p2(1), p2(2), 8)

       call kxplsl(z0, expa, z1,z2, eps, p,q, jcon)
!         jcon 0: xp found on (z1->z2) segment. 0<=p<=1
       if( jcon /= 0 ) then
          icon = -1
       else
          ! back to 3D
          el = q/sint
          xp(:) = ps(:)+ el*dir(:)

          if( q <=0. ) then
             icon = -1
          elseif( xp(3) >= p1(3)-h*eps2 .and. &
               xp(3) <= p1(3)+ h +h*eps2 ) then
             icon = 0
          else
             icon = -1
          endif
       endif
    endif
  end subroutine kxplvsq


    
