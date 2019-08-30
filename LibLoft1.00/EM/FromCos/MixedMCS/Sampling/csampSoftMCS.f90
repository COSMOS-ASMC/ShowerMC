!   sample soft multiple scattring angle between
! two hard scattering separated by 's'
!  Before using this
! 1)  cfixMixedConst must be called
!  to fix e- or e+ and media (required in cgetAveMul12)
! 2)  ciniFsot must be called
! 3) Then, angle sampling is possible by calling
!       csampSoftMCSang(mu, cosa)
! Test program is at bottom. This one dose not use
! cgetAveMul12 and asfot and bsoft parameters are
! fixed by hand.
!
module modsoftMCS
  real(8),save:: asoft, bsoft
  real(8),save:: avemu, avemu2 
end module modsoftMCS

subroutine ciniFsoft(Ein, s)
  use modsoftMCS
  implicit none
  real(8),intent(in):: Ein ! e-/e+ energy in eV
  real(8),intent(in):: s   ! dist between two hard col. in g/cm2
!          get <mu>  <mu2>
  call cgetAveMu12(Ein, s, avemu, avemu2)
  bsoft =( 2*avemu - 3*avemu2 )/(1-2*avemu)
  asoft = (1-2*avemu) + bsoft
end subroutine ciniFsoft

function  cFsoft(mu) result(ans)
  use modsoftMCS
  implicit none
  real(8),intent(in):: mu
  real(8):: ans

  real(8),external:: cUsoft

  ans = asoft*cUsoft(0.d0, bsoft, mu) + (1-asoft)*cUsoft(bsoft,1.0d0,mu)

end function cFsoft

function cUsoft(a,b,mu) result(ans)
  implicit none
  real(8),intent(in):: a, b
  real(8),intent(in):: mu
  real(8):: ans
  if(mu >= a .and. mu <= b) then
     ans = 1./(b-a)
  else
     ans = 0.
  endif
end function cUsoft

subroutine csampSoftMCSang(mu, cosa)
  use modsoftMCS 
  implicit none
  real(8),intent(out):: mu   ! sampled polar  angle  mu=(1-cosa)/2
  real(8),intent(out):: cosa   ! sampled polar  angle 
  
  real(8):: u

  call rndc(u)
  if(u < asoft) then
     mu = bsoft/asoft* u
  else
     mu = (u-asoft)/(1-asoft) *(1-bsoft) + bsoft
  endif
  cosa = (1-2*mu)
end subroutine csampSoftMCSang

!program  test
!!  fort csampSoftMCS.f90 -L$COSMOSTOP/lib/PCLinuxIFC64 -lcosmos
!  use modsoftMCS 
!  implicit none
!  real(8),external:: cFsoft
!  integer:: i
!  real(8):: mu, cosa
!
!  write(0,*) 'Enter a(0.3), b(0.1)'
!  asoft = 0.3d0
!  bsoft = 0.1d0
!  read(*,*) asoft, bsoft
!  write(0,*)  asoft, bsoft
!  
!  do i = 1, 100000
!     call csampSoftMCSang(mu, cosa)
!     write(*,*) mu
!  enddo
!  mu = 0.
!  do while (mu< 1.) 
!     write(*,*) '# ', mu, cFsoft(mu)
!     mu = mu + 0.01d0
!  enddo
!end program test

     


     

