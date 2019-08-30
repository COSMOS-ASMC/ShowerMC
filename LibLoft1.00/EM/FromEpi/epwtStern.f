      subroutine epwtStern(sh)
      implicit none
#include "Zstern.h"
       type(sternh)::  sh
      write(*,*) sh%a, sh%b, sh%c, sh%x0, sh%x1,  sh%sa
!      write(*,*) sh.tcut, sh.w0, sh.wlg0
      write(*,*) sh%tcut, sh%w0   ! these two may change
                          ! for each component v9.153
      end
