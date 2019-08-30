!      implicit none
!      real*8 a, b, alpha, beta, x
!      integer i
!
!      a = -10.
!      b = 1.
!      alpha = -2.
!      beta =-1.
!      do i = 1, 50000
!         call ksampLin(a, b, alpha, beta, x)
!         write(*,*) sngl(x)
!      enddo
!      end
      subroutine ksampLin(a, b, alpha, beta, x)
!        samples a random variable with density 
!        f(x)dx = ax + b in x=(alpha, beta); alpha  <=  beta
!        f(x) need not be normalized.
!        f(x) > 0  in the region since it is prob. density.
!
      implicit none
      real*8 a ! input.
      real*8 b ! input.
      real*8 alpha ! input
      real*8 beta  ! input.  beta > alpha.
      real*8 x     ! output.
!
      real*8 u1, u2, u
      real*8 fa     ! f(alpha)
      real*8 fb     ! f(beta)
      real*8 ba     ! beta -alpha
!
      ba = beta - alpha
      if(ba .eq. 0.) then
         x = beta
      elseif(ba .lt. 0.) then
         stop 'error input to ksampLin'
      else
         fa = a*alpha + b
         fb = a*beta  + b
         if(fa .lt. 0 .or. fb .lt. 0.) then
            stop 'bad region to ksampLin'
         endif

         call rndc(u)
         if(u .lt. abs(fb - fa)/(fb + fa)) then
            call rndc(u1)
            call rndc(u2)
            if(fa .lt. fb) then
               x = max(u1, u2)
            else
               x = min(u1, u2)
            endif
            x = ba * x + alpha
         else
            call rndc(u)
            x = ba * u + alpha
         endif
      endif
      end
