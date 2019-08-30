      implicit none
!     The result (ans) of this calculation can be fitted by the
!     following 4th order polynomial.  
!
!         log(ans) = sum c_i p**i  (i=0 to 3)
!
!     c0        .77099  
!     c1        1.3470  
!     c2        .12049   
!     c3       -.57001E-02 
!
!  ans  gives integral 1 to 5 of 1/(x- sqrt(x**2-1))**p/sqrt(x**2-1)
!      for p=0.1 to 10. 
!     the result is used in the decay in flight.
!
      real*8 p
      common /zzz/p
      integer i, icon
      real*8 func
      external func
      real*8  ans, error, eps

      eps = 1.d-5
      do i = 1, 21
         p= 0.1d0*10.**((i-1)*0.1d0)
         call kdexpIntF(func, 1.0d0, 5.0d0, eps,
     *        ans, error, icon)
         write(*,*) sngl(p), sngl(ans), log(sngl(ans))
      enddo
      end
!     *********************************
      real*8 function func(xa)
      real*8 xa(2)
      real*8 p
      common /zzz/p
      real*8 x
!      write(*,*) xa
      x = xa(1)
      if(xa(2) .lt. 0.) then
!         x = 1.d0 - xa(2) so x-1=-xa(2)
!
         func =
     *   1.d0/(x-sqrt(-xa(2)*(x+1.d0)))**p/
     *   sqrt(-xa(2)*(x+1.d0))
      else
         func = 1.d0/(x- sqrt(x**2-1.d0))**p
     *          /sqrt(x**2-1.d0)
      endif

      end

