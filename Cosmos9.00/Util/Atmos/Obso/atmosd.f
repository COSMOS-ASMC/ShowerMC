#include "BlockData/cblkGene.h"
      implicit none
#include "Zmanagerp.h"
#include "Zatmos.h"
      integer i

      call creadParam(5)
      call cerrorMsg(AtmosFile, 1)
      call creadAtmosD
!
      call catmosCnst1

      write(*,'(a,a)')
     * '#    Height        Temp         Press          rho  ',
     * '            a             d0           cumd            H'
      write(*,'(a)') "#-----------------------------------------"

      do i = 1, atmos%nodes-1
         write(*,
     *   '(8g14.5)')
     *   atmos%z(i), atmos%T(i),  atmos%P(i),
     *   atmos%rho(i), atmos%a(i),  atmos%d0(i),
     *   atmos%cumd(i), atmos%H(i)
      enddo

      end

