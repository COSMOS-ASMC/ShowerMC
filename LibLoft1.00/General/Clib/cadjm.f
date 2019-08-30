!       *****************************************************
!       * 
!       *  cadjm: adjust momentum when energy has changed
!       *          or mass has changed. 
!       *  but keep the direction of \vec p unchanged
!       ***************************************************
!  call cadjm(p, q)
!  
!  p:  /ptcl/  Input.   4 momeuntum and mass must be given.
!  q:  /ptcl/  Output.  Momentum of p is adjusted and stored in q
!              q can be p.;  In this case, call ce2pp(p) can be used.
!
      subroutine cadjm(p, q)
!----      include '../Zptcl.h'
#include  "Zptcl.h"
      type(ptcl):: p, q
!      
      real*8  pabs, cpabs, r
!        |p| from 3 momentum
      call cpxyzp(p%fm, pabs)
!        true |p|      
      cpabs = p%fm%p(4)**2 - p%mass**2
      if(cpabs .gt. 0.d0 .and. pabs .gt. 0.d0) then
           cpabs = sqrt(cpabs)
           r = cpabs/pabs
           q%fm%p(1) = p%fm%p(1) * r
           q%fm%p(2) = p%fm%p(2) * r
           q%fm%p(3) = p%fm%p(3) * r
      else
           q%fm%p(1) =0.
           q%fm%p(2) =0.
           q%fm%p(3) =0.
      endif
      q%fm%p(4) = p%fm%p(4)
      q%mass = p%mass
      end
