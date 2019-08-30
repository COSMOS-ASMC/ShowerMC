!     ************************************
      subroutine eptotcmuB(vmin, vmax, ans)
!     ************************************
      implicit none
      
      real*8 vmin, vmax, ans
!
!
!        integration of brems function from vmin to vmax.
!        result is in mb
!        brems function is  ds/dv= epmuBdsdv(v)/v
!
      external  epmuBrSdsdv
      real*8    epmuBrSdsdv

      real*8 v1, v2,  ans1

      ans=0.
      v2=vmax
      v1 = 0.
      do while (v1 .ne. vmin)
         v1=max( v2/10.d0,  vmin)
         call k16pGausLeg2(epmuBrSdsdv, v1, v2, 16,  ans1)
         ans=ans+ans1
         v2=v1
      enddo
      end
!    *************************
      real*8 function epmuBrS(v)
!    *************************
      implicit none
#include "Zmedia.h"
#include "ZmuBPNgene.h"
!
!       sum of (   v dsigma/dv ) for a target
! 
      real*8 v         ! input. Eg/Emu

      integer i
      real*8 sum, vmax,  epmuvmax
      real*8 epmuBrem  !  v* ds/dv

      sum = 0.
      do i = 1,  media%noOfElem
         call epmuSetCnst(media%elem(i)%Z, media%elem(i)%A)
         vmax = epmuvmax(Emu)
         if(v .lt. vmax)  then
            sum = sum + epmuBrem(Emu, v) * media%No(i)
         endif
      enddo
      epmuBrS = sum
      end
!     ************************************
      subroutine epmuElossB(vmin, vmax, ans)
!     ************************************
      implicit none

      real*8 vmin, vmax, ans
!
!      Int(v=vmin  to vmax) of v*ds/dv of  brems function
!      result is in mb.
!   Note:      v* ds/dv = epmuBrS (=sum of epmuBrem)
!
      external  epmuBrS
      real*8    epmuBrS
      real*8    v1, v2, ans1
      v1 = vmin
      v2 = vmax
      ans  = 0.
!      do while (v2 .ne.  vmax)
!         v2 =  min(v1*10.d0**0.5, vmax)
         call k16pGausLeg2(epmuBrS, v1, v2, 16,  ans1)
         ans  =  ans +  ans1
!         v1=v2 
!      enddo
      end
!     *********************
      real*8 function epmuBrSdsdv(v)
      implicit none
#include "Zmedia.h"
#include "ZmuBPNgene.h"
!
!        sum of dsigma/dv  for a target in mb
!
      real*8 v
      real*8 epmuBrS

         
      epmuBrSdsdv = epmuBrS(v)/v

      end
