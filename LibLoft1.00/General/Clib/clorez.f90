!  for test see Test/testLorez/test.f
!       *
!       * clorez: Lorentz transformation in z direction.
!       *         relative accuracy is always ~16 digits.
!       *  (24,Jun.2016)
!       *  This one is simpler  and
!       *  better than old one (of which accuracy is
!       *  sometimes 12~13 digits)   
!       *********************************************************
!
! /usage/  call clorez(gb, p,  po)
!
!        Suppose two systems K and K'.  K' is moving with a
!        constant velocity relative to K (4 velocity is 
!        in gb).  The axises in the both system are parallel.
!        K' is moving along the direction of the z axis.
!        (i.e., gb=(0, 0, gb(3), gamma)).
!        p is a 4 momentum given in K'.  This routine
!        transforms p into po seen from K.
!
subroutine clorez(gb, p, po)
  implicit none
#include  "Zptcl.h"
  type(fmom),intent(in)::  gb
  type(ptcl),intent(in)::  p
  type(ptcl),intent(out):: po  ! po can be p

  type(ptcl):: qo        ! output to be copied to po
!
  real(8):: tmass2       ! transverse mass^2
  real(8):: gbPz, beta   ! gamma* beta*Pz, beta
  real(8):: E, Pz, g
  
  !
  E= p%fm%p(4)
  Pz = p%fm%p(3)
  qo = p
  gbPz = gb%p(3)*Pz
  g = gb%p(4)
  if( gbPz >= 0.0d0 ) then
     qo%fm%p(4)=g*p%fm%p(4) + gbPz
     qo%fm%p(3)=g*p%fm%p(3) + gb%p(3)*p%fm%p(4)
  else
     beta = gb%p(3)/g
     tmass2 = dot_product(p%fm%p(1:2), p%fm%p(1:2)) + p%mass**2

     qo%fm%p(4) = g*((Pz/g)**2 + tmass2 )/ (E-beta*Pz)
     qo%fm%p(3) = g*(( E/g)**2 - tmass2 )/ (Pz-beta*E)
  endif
  po = qo
end subroutine clorez
