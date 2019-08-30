       subroutine ctestOnShell(com, p)
#include "Zptcl.h"
        character*(*) com
        type(ptcl):: p
        real*8 erg
!        write(*, *) '---------',com
        erg =sqrt(
     *  p%fm%p(1)**2 +
     *  p%fm%p(2)**2 +
     *  p%fm%p(3)**2 +
     *  p%mass**2  )
!        write(*, *) ' sqrt(p^2 + m^2)=', erg,  ' E=', p.fm.p(4)
        if(abs(p%fm%p(4)/erg -1.) .gt. 1.d-6) then
           write(*, *) '---------',com
           write(*, *) ' sqrt(p^2 + m^2)=', erg,  ' E=', p%fm%p(4)
           write(*, *) ' code, mass =',p%code, p%mass
        endif
        end
