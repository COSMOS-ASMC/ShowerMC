#include "ZcosmosBD.h"
#include "ZepicsBD.h"


#include "jamdummy.h"
#include "phitsdummy.h"

c       ****************
        program stdepics
        implicit none
#include  "ZepManager.h"

c          block data names
         external epblksepi
#include  "ZcosmosExt.h"

         integer ngen


!//////////
!         call mydummy  !  to see jam common block setting
!         call myPhitsDummy  ! // phits
!/////////////
         call sepics(ngen)

         if(ngen .eq. 0) then
             stop  100
         endif
         end
c         following 3 are to bypass the problem of unresolved
c         external references (in connection with Cosmos).
         subroutine chookCeren
         end
         subroutine chookCerenS
         end
         subroutine chookCerenE
         end

