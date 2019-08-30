        subroutine crot3mom(p, a, n)
!           a(n)'s momentum is given in a system whose z axis
!         coinsides with momentum direction of p.
!         This subroutine converts 'a' so that 
!         it is seen from the system where p is defined.
!
        implicit none
#include "Zptcl.h"
        integer n
        type(ptcl)::p, a(n)
        integer i

        if( (p%fm%p(1) .ne. 0.d0 .or. p%fm%p(2) .ne. 0.d0 )
     *    .or. p%fm%p(3) < 0.d0)  then
           do i = 1, n
              call crot3vec(p%fm, a(i)%fm, a(i)%fm)
           enddo
        endif
        end

