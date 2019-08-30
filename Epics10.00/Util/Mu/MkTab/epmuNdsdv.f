!      v*ds/dv of muon nuclear interaction for an atom (in mb)
!
!         v = E_t/Emu
!
      real*8 function epmuNdsdv(v)
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "ZmuBPNgene.h"
#include "Zmuint.h"
#include "Zmass.h"
      real*8  v  ! input. E_t/Emu : E_t is the transfer energy.

!           muon Energy must be given in ZmuBPNgene.h
      real*8 mu2, m12, m22
      parameter (mu2 = masmu*masmu, m12 = 0.54, m22 = 1.8)

      real*8  x, Gx, k, t, epmuG
      
      real*8 Eg, xs, f1, f2, f3

      Eg = Emu * v
      call cgpxs1(Eg, xs)       ! cosmos
      if(xs .le. 0. .or. v .gt. 1.0) then
         epmuNdsdv =  0.
      else
         x =  A3*0.00282d0 *xs*1.d3  ! xs to be in micro barn
         Gx = epmuG(x)
     
         k = 1 + 2*(1/v -1)/v
         t = mu2*v*v/(1-v)
      
         f1 = 0.75d0*Gx*(k*log(1+m12/t) - k*m12/(m12+t)
     *         - 2*mu2/t)
         f2 = 0.25d0*(k*log(1+m22/t) - 2*mu2/t)
         f3 = Gzai*mu2*2/t*( 0.75d0*Gx* m12/(m12+t) 
     *          + 0.25d0*m22/t*log(1+t/m22))
      
         epmuNdsdv = alpha/2/pi*A*xs * v**2* (f1 + f2 + f3)
      endif

      end
      real*8 function epmuG(x)
      implicit none
      real*8 x
      epmuG = 3./x**3*(x*x/2 -1.0 + exp(-x)*(1.0+x))
      end



