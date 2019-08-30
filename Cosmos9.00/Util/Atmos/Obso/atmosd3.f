#include "BlockData/cblkGene.h"
      implicit none
!  #include "BlockData/cextGene.h"
#include "Zmanagerp.h"
#include "Zcondc.h"
#include "Zatmos.h"


      integer i

      real*8 cthick2h
      real*8  h
      real*8 t
!      real*8 sh2
#if ATMOSPHERE == 1
      call creadParam(5)
      call creadAtmosD
!
      call catmosCnst1
      call catmosCnst2
#endif
      t = 10300.
      write(*,'(a,a)') '#   h     press  depth '
      do i = 1, 1000000
          h = cthick2h(t)
          write(*, *) sngl(t),  sngl(h)
          t = t -50.d0
          if(t .lt. 10.d0) goto 10
      enddo
 10   continue
      end

