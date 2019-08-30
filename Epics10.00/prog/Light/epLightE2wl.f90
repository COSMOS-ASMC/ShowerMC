module  modepLightE2wl
  !   hneu = E ;  l0 neu = c; l0E/h=c; l0 =hc/E
  !   neu l0= c ;  nue l = c/n;   l/l0 = 1/n; l = l0/n 
  real(8),parameter::hbarc=197.0e6*1e-15/1.e-9 ! hbarc in eV nm
  real(8),parameter::hc=hbarc*2*3.14159
end module modepLightE2wl
! ???  we may use always wl0; since we cannot mesure wl;
!   E.g wave length distribution from scintillation is based on
!   ~wl0.   ???
!   ***** use alwasy n=1.0d0  *****
subroutine epLightE2wl(E, n,  wl0, wl)
  use modepLightE2wl
  implicit none
  real(8),intent(in)::E  ! energy of light in eV
  real(8),intent(in)::n  ! refraction index of the medium
  real(8),intent(out)::wl0  ! wave length (nm) in the vacuumn
  real(8),intent(out)::wl ! wave length (nm) in the medium
  wl0 = hc/E
  wl = wl0/n 
end subroutine epLightE2wl

subroutine epLightwl2E(wl, n, wl0, E)
  use modepLightE2wl
  implicit none
  real(8),intent(in)::wl ! wave length (nm) in the medium
  real(8),intent(in)::n  ! refraction index of the medium
  real(8),intent(out)::wl0  ! wave length (nm) in the vacuumn
  real(8),intent(out)::E  ! energy of light in eV
    
  wl0 = wl*n
  E = hc/wl   ! we assume wl is always that in vacuum


end subroutine epLightwl2E

