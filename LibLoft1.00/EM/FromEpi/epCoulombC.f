      real*8 function epCoulombC(z)
      implicit none
      real*8 z  ! input.  z= (Z/137)**2
      
      epCoulombC =
     1   z * (1./(1. + z) +0.20206 -0.0369*z
     2    + 0.0083*z**2 - 0.002*z**3)
!
!      next is by Tsai.  good but little bit smaller
!      than correct values at z > 0.6
!      = 1.202*z - 1.0369*z**2 + 1.008*z**3/(1.+z)
!
      end
