!          integral of muon pair creation function v* ds/dvdr by dr
!     so that the result is v*ds/dv. (in mb)
!
!     .  v=E_pair/Emu.  r=(E+-E-)/E_pair
!
      real*8 function epmudsdv(v)
      implicit none
#include "Zmedia.h"
#include "ZmuBPNgene.h"
      real*8  v  ! input. E_pair/Emu

!           muon Energy must be given in ZmuBPNgene.h

      real*8  vmin, vmax, rhomax, ans

      real*8  epmuPairRmax, epmuvmax, epmuPairVmn

      real*8   epmudsdvdrv
      external epmudsdvdrv
!         to be  common  with epmudsdvdrv
      real*8 vv, rr
      common /mupairint/ vv, rr
!

      
      vmax = epmuvmax(Emu)
      vmin = epmuPairVmn(Emu)
      if(v .le. vmin .or. v .ge. vmax) then
         ans = 0.
      else
         rhomax = epmuPairRmax(Emu, v)
         vv = v
         call k16pGaussLeg(epmudsdvdrv, 0.d0, rhomax, 16,  ans)
      endif
      
      epmudsdv  = ans *2
      end
!     ***********************
      real*8 function epmudsdvdrv(r)
      implicit none
#include "Zmedia.h"
#include "ZmuBPNgene.h"
      real*8 r  !  integration variable (rho=(E+- E-)/(E+ +E-)

!         to be  common  with epmudsdvdrx
      real*8 vv, rr
      common /mupairint/ vv, rr
      real*8  epmudsdvdr
!       note that   epmudsdvdr is v* ds/dvdr
      epmudsdvdrv = epmudsdvdr(Emu, vv, r)
      end


