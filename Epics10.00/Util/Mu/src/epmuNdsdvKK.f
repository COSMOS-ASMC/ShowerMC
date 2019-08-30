!  This is the cross-section used by Kobayakawa
!
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
      real*8 mu2
      parameter (mu2 = masmu*masmu)


      
      real*8 Eg,  xs,  epvmhd1, epvmhd2

      Eg = Emu * v
      call cgpxs1(Eg, xs)       ! cosmos
      if(xs .le. 0. .or. v .gt. 1.0) then
         epmuNdsdv =  0.
      else
         epmuNdsdv = alpha/pi/2*A*xs *
     *          ( epvmhd1(v)  + epvmhd2(v) )
     
      endif

      end
      real*8  function epvmhd1(v)
      implicit none
!         point like part 
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "ZmuBPNgene.h"
#include "Zmuint.h"
#include "Zmass.h"
!                 point like
       real*8 v

       real*8 g

!             v*fai1(v) for nuclear interaction  of mu

       g = Emu/masmu

         epvmhd1 = ((1.d0+(1.d0-v)**2)*log(2* g*8.88d0*(1.-v)/v)
     *       + 2*(v-1.)
!    *       + ( 1./g**2*v**2/(1.-v) - 2*8.88/g*v)/2
!    *       + v/g/8.88
     *                )
         epvmhd1 = epvmhd1 * PointLike
        end
!       ************************
      real*8 function epvmhd2(v)
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "ZmuBPNgene.h"
#include "Zmuint.h"
#include "Zmass.h"

!             v*fai2(v) for nuclear interaction  of mu
!        hadron like
      real*8 v
      real*8 g, alm2
      parameter (alm2=0.365d-6)
      
      g = Emu/masmu

      epvmhd2=   (v-1.)
!    *    + v*g/8.88/2
!    *           +alm2/Emu**2/4*log(  ( 2*v+alm2/g/8.88)/
!    *          (v**2/(1.-v)*g/8.88+ alm2/g/8.88) )
     *       +   (1.-v+ v**2/2*(1.+2* masmu**2/alm2) ) *
     *           log( (1.+ (1.-v)/v**2 * alm2/masmu**2) )
!    *  log( ( v+ (1.-v)/v*alm2/masmu**2 )/(v + alm2/2/masp/Emu) )
      epvmhd2 = epvmhd2* (1.-PointLike)*Shadow
      end
