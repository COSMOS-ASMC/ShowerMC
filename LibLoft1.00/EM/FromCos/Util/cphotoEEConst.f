!    
      implicit none
      real*8  z, a, x0inkgpm2
      read(*, *) z, a, x0inkgpm2
      call cphotoEEConst(z, a, x0inkgpm2)
      end
!          compute Photo Electric Effect constants.
!      ************
       subroutine cphotoEEConst(z, a, x0inkgpm2)
!      ************
!
       implicit none
       real*8  z, a, x0inkgpm2, az, b0, b1, b2, fa, 
     * cnsta, cnstp,  ek, cnstl
       real*8 emass/.511e-3/
!                  make table
       az=z/137.
       b0 = 1./(  ( .9663*az + 5.023) *az + .9211)
       b1 = (2.56*az - 2.632)*az + 1.90
       b2 = ( (-6.563*az+ 8.25)*az - 5.616)*az + 2.097
       fa = ( .2762*az - .0288) *az + 1.083
       cnsta = emass - 13.5e-6*z**2
       cnstp = .06 * az**4 *z/a* x0inkgpm2
!                  eg>the  then no photoe
!             the(i)= 10.e-3*(z/82.)**4.5
       ek = 13.6e-6*z**2
       cnstl= (8. - z*2.5/80.)
       write(*, *) b0, b1, b2, fa, cnsta, cnstp, ek, cnstl
       end
