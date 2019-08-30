#include "BlockData/cblkGene.h"
      implicit none
#include "Zcondc.h"
#include "Zmanagerp.h"
#include "Zatmos.h"
#include "ZcosmosExt.h"

      integer i

      real*8 cvh2den, cvh2denp, cvh2den2p, cvh2temp, cthick2h
      real*8 cvh2scaleh, cvh2thick, cthick2den, h

      real  rho, rhop, rhopp, temp,  sh, h1, den
      real*8 t
!      real*8 sh2
#if ATMOSPHERE == 1
      call creadParam(5)
      call creadAtmosD
!
      call catmosCnst1
      call catmosCnst2
#elif ATMOSPHERE == 2
!          read segmented atmosphere data
      call creadAtmosD
!          manipulate data
      call catmosCnst1
#endif
      h = -1000.
      write(*,'(a,a)') '#   h       rho    rhop     rhopp   Temp ',
     *                 '        depth      H        h   rho '
      do i = 1, 1000000
          rho= cvh2den(h)
          rhop = cvh2denp(h)
          rhopp = cvh2den2p(h)
          temp = cvh2temp(h)
          t = cvh2thick(h)
          sh = cvh2scaleh(h)
          h1 = cthick2h(t)
          den = cthick2den(t)
!          sh2 = - cvh2den(h)/cvh2denp(h)
!          write(*, '(1pE12.4,)')  h, rho, rhop, rhopp, temp, t, sh, 
!     *       h1, den
          write(*, '(1p3E12.4)')  h,   t, rho
          h = h + 25.d0
          if(h .gt. 100.d3) goto 10
      enddo
 10   continue
      end


