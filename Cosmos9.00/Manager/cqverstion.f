      subroutine cqversion(cosv)
      implicit none
#include "Zmanagerp.h"
      character*8::cosv ! output  cosmos version such as 7.58; left justified
      character*64::COSMOSTOP
      character*128:: filen
      integer kgetenv2, icon

      if( kgetenv2("COSMOSTOP", COSMOSTOP) == 0 ) then
         write(0,*) " Environmental variable "
         write(0,*) "COSMOSTOP cannot be obtained in cqversion"
         stop
      endif
      filen = trim(COSMOSTOP)//"/Version/version" 
      call copenf(TempDev, filen, icon)
      read(TempDev, '(a)') cosv
      close(TempDev)
      end
