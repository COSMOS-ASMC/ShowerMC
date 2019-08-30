!          integral of muon pair creation function ds/dvdr by v
!  so that  the result is ds/dr (in mb)
!           r = assymetry factor, v= E_pair/Emu
!
      real*8 function epmudsdr(r)
      implicit none
#include "ZmuBPNgene.h"
      real*8  r  ! input. (E+- E-)/E_pair  

!           muon Energy must be given in ZmuBPNgene.h

      real*8  vmin, vmax, rhomax, ans

      real*8  epmuvmax, epmuPairVmn

      real*8   epmudsdvdrr
      external epmudsdvdrr
!         to be  common  with epmudsdvdrr
      real*8 vv, rr
      common /mupairint/ vv, rr
!
      real*8  v1,  v2,  ans1

      vmax = epmuvmax(Emu)
      vmin = epmuPairVmn(Emu)

      rr =  r
      ans  = 0.
      if(abs(r) .lt.  rhomax) then
         v2 =  vmax
         do while(v1  .ne.  vmin)
            v1 = max(vmin, v2/10.d0)
            call k16pGaussLeg(epmudsdvdrr, v1, v2, 16,  ans1)
            ans =  ans  +  ans1
            v2 = v1
         enddo
      endif
      epmudsdr  = ans 
      end
!     ***********************
      real*8 function epmudsdvdrr(v)
      implicit none
#include "ZmuBPNgene.h"
      real*8 v  !  integration variable (rho=E++ E-)/Emu

!         to be  common  with epmudsdvdrx
      real*8 vv, rr
      common /mupairint/ vv,rr
      real*8  epmudsdvdr  ! is v* ds/dvdr
      epmudsdvdrr = epmudsdvdr(Emu, v, rr)/v
      end

