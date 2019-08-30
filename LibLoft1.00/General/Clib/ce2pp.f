!     ce2pp: When energy of a partilce is changed without chaning the
!     direction, we may need to adjust  px, py, pz.  ce2p uses direciton
!     cosine but this one use current direction (px,py,pz) so that
!     input can be only type(ptcl)      
!     direction vector and mass.
!        same  call cadjm(atpcl, aptcl)
      subroutine ce2pp(aptcl)
      implicit none
#include  "Zptcl.h"
      type(ptcl),intent(inout)::aptcl
!
      real(8):: p, pc
      p  = aptcl%fm%p(4)**2 - aptcl%mass**2

      if(p < 0.) then
         aptcl%fm%p(1:3) = 0.
         aptcl%fm%p(4) =  aptcl%mass
      else
         pc = sqrt( dot_product( aptcl%fm%p(1:3),aptcl%fm%p(1:3)) )         
         p = sqrt(p)
         if(pc > 0.) then
            aptcl%fm%p(1:3) =  aptcl%fm%p(1:3)*p/pc
         else
            aptcl%fm%p(1:3) = 0.
            aptcl%fm%p(4) =  aptcl%mass
         endif
      endif
      end
