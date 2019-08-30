!        density effect correction 
!
      subroutine epdEdxDenC(media, gin, delta)
      implicit none
#include "Zmedia.h"
       type(epmedia):: media     ! input 
      real*8 gin  ! input  ! gamma factor of the particle
      real*8 delta  ! output density correcton foactor

      real*8 g
      if( media%rhoc .eq. 1.d0 ) then
         g = gin
      else
         g = sqrt( (gin**2-1.0d0) * media%rhoc + 1.d0 )
      endif

      call epdEdxDenC0(media%sh, g, delta)
      end


      subroutine epdEdxDenC0(sh, g, delta)
      implicit none
#include "Zstern.h"
       type(sternh):: sh     ! input 
      real*8 g ! input  ! gamma factor of the particle
      real*8 delta  ! output density correcton foactor
      
      real*8 x, cbar
      real*8 tln10/4.60517/
      save tln10

      x=log10( (g- 1.)*(g+1.) ) / 2 ! = log10(gbeta) = 0.4343log(gbeta)        


      if(x .lt. sh%x0) then
         if( sh%delta0 .eq. 0) then
            delta = 0.
         else
            delta = sh%delta0 *10.0d0**(2*(x-sh%x0))
         endif
      else
         cbar = - sh%c
         delta = tln10*x - cbar 
         if(x .lt. sh%x1) then
            delta = delta + sh%sa*(sh%x1-x)**sh%k
         endif
      endif
      end
