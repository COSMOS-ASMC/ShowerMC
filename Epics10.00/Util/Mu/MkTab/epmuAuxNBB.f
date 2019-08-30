!     This one uses cross-section given by Bezrukov and Bugaev (
!     Sov. J. N.P (vol.33, May 1981, p.635)
!     This should also be usable for Kobayakawa's formula.
!
!     ************************************
      subroutine eptotcmuN(vmin, vmax, ans)
!     ************************************
      implicit none
      
      real*8 vmin, vmax, ans
!
!
!        integration of mu Nuc int. function from vmin to vmax.
!        result is in mb
!
      external  epmuNSdsdv
      real*8    epmuNSdsdv


      real*8 v1, v2

      real*8 ans1

!      real*8  epsa, err, icon
!      data epsa/1.d-7/

      ans=0.
      v2 = vmax
      v1 = 0.
!      v1 = vmin
      do while (v1 .ne. vmin)
         v1=max( v2/5.d0,  vmin)
!         call kdexpIntFb(epmuNSdsdv, v1, v2, epsa, ans1, err, icon)

         call k16pGausLeg2(epmuNSdsdv, v1, v2, 10,  ans1)
!
!
         ans=ans+ans1
         v2=v1
      enddo
      end
!    *************************
      real*8 function epmuNS(v)
!    *************************
      implicit none
#include "Zmedia.h"
#include "ZmuBPNgene.h"
!
!       sum of (   v dsigma/dv ) for a target
! 
      real*8 v         ! input. E_t/Emu (E_t transferred energy)

      integer i
      real*8 sum
      real*8 epmuNdsdv  !  v* ds/dv

      sum = 0.
      do i = 1,  media%noOfElem
         call epmuSetCnst(media%elem(i)%Z, media%elem(i)%A)
         sum = sum + epmuNdsdv(v) * media%No(i)
      enddo
      epmuNS = sum
      end
!     ************************************
      subroutine epmuElossN(vmin, vmax, ans)
!     ************************************
      implicit none
      real*8 vmin, vmax, ans
!
!      Int(v=vmin  to vmax) of v*ds/dv of  nuc.int function
!      result is in mb.
!   Note:      v* ds/dv = epmuNS (=sum of epmuNdsdv)
!
      external  epmuNS
      real*8    epmuNS
      real*8    v1,  ans1, v2    

      real*8  epsa, err,  icon



      data epsa/1.d-7/

      v1 =  vmin
      v2 = 0.
!      v2 = vmax
      ans  = 0.
      do while (v2 .ne.  vmax)
         v2 =  min(v1*10, vmax)
!         call k16pGausLeg2(epmuNS, v1, v2, 16,  ans1)
         call kdexpIntFb(epmuNS, v1, v2, epsa, ans1, err, icon)

         ans  =  ans +  ans1
         v1=v2 
      enddo
      end
!     *********************
      real*8 function epmuNSdsdv(v)
      implicit none
#include "Zmedia.h"
#include "ZmuBPNgene.h"
!
!        sum of dsigma/dv  of a target in mb
!
      real*8 v
      real*8 epmuNS

         
      epmuNSdsdv = epmuNS(v)/v

      end
