      implicit none
#include "Zmanagerp.h"
#include "ZrigCut.h"

      character*80 file
      real*8 rigcut, azm, zen, maxZ
!
      file = '../Contrib/Data/Cutoff/Kamioka'
      write(ErrorOut, *) ' enter file name=',file
      read(*, *) file
      call crigCut0(file)
      maxZ = MinZen + DZen * ( ZenSize -1 )
      do azm = MinAzm, 360.00001d0, 5.d0
         do zen = MinZen, maxZ, DZen/2.
            call cgetrigCut(azm, zen, rigcut)
            write(*, *) sngl(azm), sngl(zen), sngl(rigcut)
         enddo
         write(*, *)
      enddo
      end
