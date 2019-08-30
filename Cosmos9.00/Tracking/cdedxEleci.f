       subroutine cdedxEleci(w0in,  knck)
       implicit none
!      ***********
#include  "ZdedxAir.h"
      real*8  w0in ! input.
!          kinetic energy of recoil electron in GeV.  k.e>w is not
!          included as energy loss so that it should be treated by
!          producing the recoil electron actually.  
!          if w0in is < KEminObs/10 or so,
!          it should be adjusted to be KEminObs/10 or so.
!      
      logical knck ! input. obsolete now

!
      jdef = 1

       w0=w0in
       w0inMeV = w0 *1000.0d0
       Knckon=knck
!          log(w*1000.)
       wlg0=log(w0) + 6.907d0
      end

