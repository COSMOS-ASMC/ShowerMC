program main
  use gammafunc
  implicit  none
  integer::iz
   complex(8):: z, ans
   !   complex(8),external:: kcgamma, kgamma
   !  complex(8),external:: kcgamma, kgamma
!   z=cmplx(-3.0d0, -7.d0, kind=8)
   z=cmplx(3.0d0, 7.d0, kind=8) ! Abramowits & Steguns
                                 ! give wrong imaginary part
   !   ans = kcgamma(z)
      ans = kgamma(z)
   write(0,*) ' by kgamma'
   write(0,*) ' z =', z, ' ans =',ans
   write(0,*) ' log ans =', log(ans)
   iz = 100.0
   !   ans = kcgamma(z)
   ans = kgamma(iz)
   write(0,*) ' z=',iz, ' ans=',ans
   z = iz
   ans = kgamma(z)
   write(0,*) ' z=',z, ' ans=',ans
 end program main
