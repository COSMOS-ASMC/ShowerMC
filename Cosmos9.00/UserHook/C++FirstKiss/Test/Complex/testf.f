      subroutine fort1(intvalue)
      integer  intvalue
      write(*,*)  intvalue
      end
      subroutine fort2(intvalue)
      integer  intvalue
      write(*,*)  intvalue
      end
      subroutine dummy
      include "fcom.h"
      write(*,*) instabc(1).i8, xxx
      end
