      subroutine Fort1(intvalue)
      integer  intvalue
      write(*,*)  intvalue
      end
      subroutine fort2(intvalue, i)
      integer  intvalue, i
      write(*,*)  intvalue, i
      end
      subroutine dummy
      include "fcom.h"
      write(*,*) instabc(1).i8, xxx
      end
      subroutine fort3(array)
      integer array(3)
      array(1) = 1
      array(2) = 2
      array(3) = 3
      end
