      subroutine fort1(intvalue)
      integer  intvalue
      write(*,*)  intvalue
      end
      subroutine fort2(intvalue)
      integer  intvalue
      write(*,*)  intvalue
      end
      subroutine confirm
      include "fcom.h"
      write(*,*) ' from fort common ',  xxx
      end
