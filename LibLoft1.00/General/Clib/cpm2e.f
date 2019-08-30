      subroutine cpm2e(p, q)
!          get energy from 3 momentum and mass
!         p: /ptcl/ Input. 3 momentum and mass must be given
!         q: /ptcl/ Outpu. q.fm.p(4) is given.  
!          q may be p
      implicit none
!----      include '../Zptcl.h'
#include  "Zptcl.h"
!
      type(ptcl):: p, q
!  
      q%fm%p(4) =  sqrt( p%fm%p(1)**2 + p%fm%p(2)**2 + p%fm%p(3)**2
     *   + p%mass**2)
      end
