      subroutine episoDirec(w)
      implicit none
#include "ZepDirec.h"
       type(epDirec)::  w  ! output. sampled isotropic direction cosine

      real*8 cosu, cs, sn, sinu

      call rndc(cosu)
      cosu = 2.*cosu - 1.
      call kcossn(cs, sn)
      sinu = sqrt(1.d0 - cosu**2)
      w%x = sinu * cs
      w%y = sinu * sn
      w%z = cosu
      end
