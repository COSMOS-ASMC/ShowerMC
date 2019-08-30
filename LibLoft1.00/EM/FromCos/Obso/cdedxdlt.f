!        density effect correction 
!
      subroutine cdedxdlt(rhoin, gin, delta)
      implicit none
      real*8  rhoin  ! input. air density in kg/m^3   
      real*8 gin  ! input  ! gamma factor of the particle
      real*8 delta  ! output density correcton foactor

      real*8 g
      real*8 rho0
      data  rho0/1.205/      ! standard density in kg/m3                        
      save rho0

      g = sqrt( (gin**2-1.0) * rhoin/rho0 + 1.d0 )
      call cdedxDenC0(g, delta)
      end


      subroutine cdedxDenC0(g, delta)
      implicit none
      real*8 g ! input  ! gamma factor of the particle
      real*8 delta  ! output density correcton foactor
!      (Z/A)       I[eV]   a       k      x0    x1     Cbar  delta0            
!      0.49919   85.7  0.1091  3.3994  1.7418  4.2759  10.5961 0.00            
!                       | this is sh.sa not sh.a                    
      real*8   shc, shx0, shx1, shsa, shk
      data  shc/-10.5961/
      data shx0/1.7418/, shx1/4.2759/, shk/3.3994/
      data shsa/0.1091/
      save   shc, shx0, shx1, shsa,  shk

      
      real*8 x, cbar
      real*8 tln10/4.60517/
      save tln10

      x=log10( (g- 1.)*(g+1.) ) / 2 ! = log10(gbeta) = 0.4343log(gbeta)        


      if(x .lt. shx0) then
         delta = 0.
      else
         cbar = - shc
         delta = tln10*x - cbar 
         if(x .lt. shx1) then
            delta = delta + shsa*(shx1-x)**shk
         endif
      endif
      end
