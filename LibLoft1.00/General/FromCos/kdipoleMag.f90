subroutine kdipoleMag(r, teta, Br, Bt, B)
  ! Compute dipole magnetic field at a given radial distance and
  ! polar angle. The radial distnace is normalized (by r0, say, such as
  ! Earth radius).  The field strength is also normalized by the value
  ! at the r0 and polar angle at 90 deg. 
  ! 
  implicit none
  real(8),intent(in):: r  ! normalized distance from the center
  real(8),intent(in):: teta ! polar angle (rad) of the position
  real(8),intent(out):: Br  ! normalized radial direction  component
  real(8),intent(out):: Bt  ! normalized polar angle component
  real(8),intent(out):: B   ! normalzied field strength

  real(8):: cost, sint

  cost = cos(t)
  sint = sin(t)

  call  kdipoleMag2(r, cost, sint, Br, Bt, B)
end subroutine kdipoleMag

subroutine kdipoleMag2(r, cost, sint, Br, Bt, B)
  implicit none
  real(8),intent(in):: r  ! normalized distance from the center
  real(8),intent(in):: cost, sint
  real(8),intent(out):: Br  ! normalized radial direction  component
  real(8),intent(out):: Bt  ! normalized polar angle component
  real(8),intent(out):: B   ! normalzied field strength

  real(8):: r3

  r3 = 1.d0 / r**3

  Br = -r3 * cost * 2
  Bt = -r3 * sint
  B =   r3 * sqrt( cost**2*3.0d0 * 1.0d0)
end subroutine kdipoleMag2

! subroutine kdipoleMagExyz(R, Lat, Long, B0, atxyz, Bxyz)
subroutine kdipoleMagExyz(R, npole, B0, atxyz, Bxyz)
  implicit  none
  real(8),intent(in):: R
!  real(8),intent(in):: Lat 
  !  real(8),intent(in):: Long
  real(8),intent(in):: nploe  ! north pole pos.
  real(8),intent(in):: B0
  real(8),intent(in):: axxyz(3)

  real(8),intent(out):: Bxyz(3)

  real(8):: rn, Bn
  
!  real(8),parameter:: torad=asin(1.0d0)/90.d0
 

  rpos = sqrt( dot_product(atxyx(:),atxyx(:)) )
  rn = rpos/R
  rsurf = rpos - R
!  npole(1) = R * sin( (90.d0-Lat)*torad ) * cos(torad*Long)
!  npole(2) = R * sin( (90.d0-Lat)*torad ) * sin(torad*Long)
!  npole(3) = R * cos( (90.d0-Lat)*torad )
  cost = dot_product(atxyz, atxyz)/rpos/R
  sint = sqrt( 1.0d0 - cost**2)
  call kdipoleMag2(rn, cost, sint, Br, Bt, B)  
  
  
  
  
