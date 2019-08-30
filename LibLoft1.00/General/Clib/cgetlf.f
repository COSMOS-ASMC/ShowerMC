      subroutine cgetlf(p,  gb)
!         get Lorentz factor of a particle or 
!         system of particles of which 4 momentum
!         is given in p.
!      p: type ptcl, Input.  4 momentum and mass 
!                           must be given
!     gb: type fmom, Output. (g*beta, g)
!
         implicit none
!----         include '../Zptcl.h'
#include  "Zptcl.h"
         type(ptcl):: p
         type(fmom):: gb
!
         if(p%mass .gt. 0.d0) then
             gb%p(1)=p%fm%p(1)/p%mass
             gb%p(2)=p%fm%p(2)/p%mass
             gb%p(3)=p%fm%p(3)/p%mass
             gb%p(4)=max(p%fm%p(4)/p%mass, 1.d0)
         else
             write(0, *) ' mass=', p%mass, ' invalid to cgetlf'
             call cbackTrace(p)
             stop 9999
         endif
        end   
      
