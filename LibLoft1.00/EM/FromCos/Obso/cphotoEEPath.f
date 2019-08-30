!c      cphotoEEPath, Energy
!      real*8 eg, egom, em/0.511e-3/,  path, prob
!      integer i
!      egom=.01
!      do  i=1,150
!          eg=egom*em
!          call cphotoEEPath(eg, prob,  path)
!          write(*, *) eg
!          egom=egom*10.**(.03333)
!      enddo
!      end
!     ****************************************************************
!     *
!     * cphotoEEPath: samples photo electric effect path in r.l
!     * cphotoEEEnergy: gives energy of photo-electron in GeV
!     *
!     ****************************************************************
!
!   /usage/  call cphotEEPath(eg, prob, path)
!            call cphotoEEEnergy(eg, e)
!
!   --input--
!    eg: photon energy in GeV
!
!   -- output --
! prob: prob/r.l 
! path: sampled path in r.l
!

!
!
!
!
       subroutine cphotoEEPath(eg, tp, path)
       implicit none
!        use XCOM result by fitting it 
#include "ZbpCnst.h"
!
!
       real*8  eg            ! input energy in GeV
       real*8  path          !  output sample path in r%l 
       real*8  tp             ! prob./r%l  
!       
!  XCOM result ; Energy in MeV.  Xsection in b
       real*8  x, u

       x = log10(eg) +3.0     ! log10(Eg/MeV)
       if(x .lt. -3.0) then
          tp = -3.0*(x+3.0) + 3.55376
       elseif(x .lt. -0.5) then
          tp = (-0.11136*x -3.5363)*x - 6.0163
       elseif( x .lt. 1.5) then
          tp =((-0.25541*x + 0.86430)*x-2.0364)*x - 5.4876
       else
          tp = -x- 5.964
       endif
       tp = 10.0**tp
!               tp in /(g/cm2) --> /r.l
       tp = tp*X0g
       call rndc(u)
       path=-log(u) / tp

       end

!      ************
       subroutine cphotoEEEnergy(egin, eout)
!      ************
       implicit none
#include "Zmass.h"

       real*8 egin
       real*8 eout
       
!
       real*8 cnsta
       data cnsta/ -1.9859375d-04 /
       save cnsta

       eout=egin + cnsta
       if(eout .le. masele) then
          eout=masele*1.0000000001d0
       endif
       end
