!       test kgamma
!     real*8  kgamma, y, x
!     do x=-0.01d0, -6.d0,-0.001d0
!          y=kgamma(x)
!          write(*,*) sngl(x), sngl(y)
!     enddo
!     end
!     *************************************************************
!     *                                                           *
!     * kgamma: gamma function in real domain. 
!     *         although real*8 must be specified, this gives
!     *         single precision accuracy.
!     *                                                           *
!     *************************************************************
!
!  Usage:  y=kgamma(x).   x
!
!     Computes gamma(x) with 6 significant digit.  gamma(x)=factorial of
!     (x-1).  gamma(1)=gamma(2)=1 .........................................
!
!
!
!
      real*8 function kgamma(x)
!
      implicit none
      real*8 x
!
      real*8 pi, z, f, t
!
!
      parameter (pi=3.141592653)
!
      if(abs(x).gt.15.d0) then
          z=x
          if(z .le. 0.d0) then
               f=pi/sin(pi*z)
               z=1.d0-z
          endif     
          kgamma=2.506628274d0*exp(-z)*z**(z-0.5d0)*
     *         ((3.47222222222d-3/z+8.3333333133d-2)/z+1.d0)
          if(x .lt. 0.d0) then
               kgamma=f/kgamma
          endif
      else
           f=1.
           z=x
           do while (z .gt. 3.0)
                z=z-1.
                f=f*z
           enddo   
           do while (z .lt. 2.0) 
               f=f*z
               z=z+1.
           enddo    
           z=z-2.0
           t = (((((1.08298598d-2*z - 3.42705226d-3)*z + 7.7549276d-2)*z
     *     + 8.01782477d-2)*z + 4.12102903d-1)*z +4.22766368d-1)* z +
     *          1.0000002d0
           if(x .lt. 2.0d0) f=1.0d0/f
           kgamma=t*f
        endif   
      end
