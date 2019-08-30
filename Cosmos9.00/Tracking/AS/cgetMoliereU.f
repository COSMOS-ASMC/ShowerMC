!      implicit none
!      real*8 dep, cosz, rmu
!      dep = 9000.
!      cosz = 0.5
!      call cgetMoliereU(dep, cosz, rmu)
!      write(*, *) ' rmu=',rmu
!      end
!            get moliere unit (in m) at 2 r.l above dep along cosz
!            direction.
      subroutine cgetMoliereU(dep, cosz, rmu)
      use modEMcontrol
       implicit none
! #include "Zelemagp.h"  ! only Es is used;   X0 and Ecrit there, no more used
#include "ZmediaLoft.h"       
       real*8 dep  ! input.  depth in kg/m2
       real*8 cosz ! input.  zenith angle
       real*8 rmu  ! output. Moliere length in m
!
       real*8  cthick2den, tmp
       real(8):: X0

       X0 = Media(MediaNo)%X0g *10.d0 ! X0g  r.l in g/cm2; X0  kg/m2

       tmp=max( dep - 2*X0*cosz, 100.d0) ! 2 r%l above
       rmu = Es/Media(MediaNo)%Ecrit * X0 / cthick2den(tmp)
       end
      
