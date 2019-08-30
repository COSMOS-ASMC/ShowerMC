      subroutine getConvFac( g2rl, g2cm)
      implicit none
#include "ZepTrackv.h"
      real(8),intent(out):: g2rl ! if x is g/cm2, x*g2rl is in r.l 
      real(8),intent(out):: g2cm ! if x is g/cm2, x*g2cm is in cm
      integer::md=1
      g2rl = 1./Media(md)%X0g 
      g2cm =  Media(md)%gtocm / Media(md)%rhoc      
      end

      
      
      
