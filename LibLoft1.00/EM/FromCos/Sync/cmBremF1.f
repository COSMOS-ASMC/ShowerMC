!      compute (1-v)/v K1(2zeta) which is the first term of
!      F in Brainerd and Petrosian's Eq(5). v is fractional energy
!      of emitted gamma ray.
!     
      real*8 function cmBremF1(up,  v)
      implicit none
      real*8 up  !  input. upsilon
      real*8 v   !  input. faractional energy.
      real*8 zeta2, cmBremF11
      
      zeta2 = v/(1.0d0-v) /3.0d0/up * 2.0d0
      cmBremF1 = (1.0d0-v)/v * cmBremF11(zeta2)
      end
!       test cmBremF11
!      program testcmBremF11
!      implicit none
!
!      real*8 x, y, cmBremF11
!
!      x = 0.0001
!      do while (x .lt. 30.)
!         y = cmBremF11(x)
!         write(*, *) sngl(x), sngl(y)
!         x = x * 10.**0.025
!      enddo
!      end
!   ************************************************************************
!    This is an approximaiton of bremExact.f by polynomial of degree 13
!    Eq. 2.9 of Erber ;  i.e., k(z)=z Inte(z,inf)dxK5/3(x).
!       At least 5 digit accuracy at  0.001 < z < 10.
!       At z < 0.001 or z> 10, assmptotic expansion is used.
!       There is some small gap at z=10. Don't worry.
!   ************************************************************************
!
      real*8 function cmBremF11(z)
      implicit none
      real*8 z
      integer pol
      parameter (pol = 13)
      real*8 coef(pol + 1)
      real*8 x, y
      integer i
      data coef/
     *  -0.4285957, -0.6851291, -0.4785896, -0.1615990, -0.4218448E-01,
     *  -0.8632375E-02, -0.1389391E-02, -0.1832289E-03,
     *  -0.2284316E-04, -0.3119965E-05,
     *  -0.4083642E-06, -0.3924065E-07, -0.2244432E-08, -0.5629939E-10/
      save coef
!
      if(z .lt. 0.001) then
         cmBremF11 = 2.1495 * z**0.3333333333 - 1.8138* z
      elseif(z .gt. 10.) then
         cmBremF11 = 1.2533 * sqrt(z) * exp(-z)
      else
         y = coef(pol + 1)
         x = log(z)
         do i = pol, 1, -1
            y = y * x + coef(i)
         enddo
         cmBremF11 = exp(y)
      endif
      end
