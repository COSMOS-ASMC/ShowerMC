function epGSf(lmax, coef1, coef2,  cost)
  implicit none 
!         compute Eq 4. of GS
  integer,intent(in)::lmax  ! max l for smmation
  real(8),intent(in):: coef1, coef2 ! see epGSGl
  real(8),intent(in)::cost  ! cos of scattering angle 
  real(8):: epGSf


  real(4):: kpnx    ! Legendre polynomial

  integer  l
  real(8):: sum, Gl
  real(8),parameter:: fpi =4* asin(1.d0)*2
  
  sum = 0.
  do l = 0, lmax
     call epGSGl(l, coef1, coef2, Gl)
     sum = sum + (2*l+1)*Gl*kpnx( sngl(cost) )
  enddo
  epGSf = sum /fpi
  end
