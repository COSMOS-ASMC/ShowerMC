! ifort integral.f90 -L$COSMOSTOP/lib/PCLinuxIFC64 -lcosmos
module gbetaAve
!        intergral of (1/gbeta^2)2 dt and 1/beta^2 dt
!        with constnt energy loss  is compared with approximate average
!           

  real(8):: E, mass, t,  dEdx

contains
  real*8 function igbeta2(x)
!           1/( g beta^2) ^2
!         =  1/(g - 1/g)^2
    implicit none

    real(8),intent(in):: x  !  path in g/cm2
    real(8):: g

    g = (E-x*dEdx)/mass
    igbeta2 = 1./(g-1./g)**2
  end function igbeta2
  real*8 function ibeta2(x)
!          1/beta^2 = 1/(1- 1/g^2)
    implicit none

    real(8),intent(in):: x  !  path in g/cm2
    real(8):: g
    g = (E-x*dEdx)/mass
    ibeta2 =1./( 1. - 1/g**2)
  end function ibeta2
end module gbetaAve

program main
  use integral
  implicit none
  integer icon
  real(8):: avibetasq, avigbetasqsq, avibetasqp
  real(8):: eps, ans1, ans2, error
  real(8):: g1, beta1, gbeta1
  real(8):: g2, beta2, gbeta2
  
  mass = 0.511
  E = 10.
  t = 0.1
  dEdx = 2
  do while(.true.)
     write(0,'(a)') 'Enter E, t, dEdx (MeV, g/cm2, MeV/(g/cm2)'
     write(0,'(a,1p,3g12.3)') 'default ', E, t, dEdx
     read(*, *, END=100 ) E, t, dEdx
     g1 = E/mass
     g2 = (E-t*dEdx)/mass
     beta1 =  sqrt(1.-1./g1/g1)
     beta2 =  sqrt(1.-1./g2/g2)
  
     avibetasq =1./( 1.- 1./g1/g2) *t
!           above is better
     avibetasqp =1./beta1/beta2 * t

     avigbetasqsq =1./( (g1-1./g1)* (g2- 1./g2) ) *t

     eps = 1.e-5
     call kdexpIntF(igbeta2, 0.d0, t, eps, ans1, error, icon)  
     write(*,*) 'icon ', icon
     call kdexpIntF(ibeta2, 0.d0, t, eps, ans2, error, icon)  
     write(*,*) 'icon ', icon
     write(*,'(a,1p,3g13.4)') 'E t dEdx', E, t, dEdx

     write(*,'(a,1p,2g13.4)') 'igbetasqsq ', ans1, avigbetasqsq
     write(*,'(a,1p,3g13.4)') 'ibetasq ',  ans2,  avibetasq, avibetasqp
  enddo
100 continue
end program main


