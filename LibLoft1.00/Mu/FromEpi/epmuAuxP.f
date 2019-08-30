!     ************************************
      subroutine eptotcmuP(vmin, vmax, ans)
!     ************************************
      implicit none
      
      real*8 vmin, vmax, ans
!
!
!        integration of pair function from vmin to vmax.
!        result is in mb
!       pair function is  ds/dv= epmudsdv(v)/v
!
      external  epmuPrSdsdv
      real*8    epmuPrSdsdv

      real*8 v1, v2,  ans1

      real*8  epsa,  err, icon
      data  epsa/1.d-8/



      ans=0.
      v2=vmax
      v1 = 0.
!      do while (v1 .ne. vmin)
!         v1=max( v2/10.d0,  vmin)
!         call k16pGausLeg2(epmuPrSdsdv, v1, v2, 16,  ans1)
         call kdexpIntFb(epmuPrSdsdv, vmin, vmax, epsa,
     *        ans1, err,icon)

         ans=ans+ans1
!         v2=v1
!      enddo
      end
!    *************************
      real*8 function epmuPrS(v)
!    *************************
      implicit none
#include "Zmedia.h"
#include "ZmuBPNgene.h"
!
!       sum of (   v dsigma/dv ) for a target
! 
      real*8 v         ! input. E_pair/Emu

      integer i
      real*8 sum, vmax,  epmuvmax
      real*8 epmudsdv  !  v* ds/dv

      sum = 0.
      do i = 1,  media%noOfElem
         call epmuSetCnst(media%elem(i)%Z, media%elem(i)%A)
         vmax = epmuvmax( Emu)
         if(v .lt. vmax)  then
            sum = sum + epmudsdv(v) * media%No(i)
         endif
      enddo
      epmuPrS = sum
      end
!     ************************************
      subroutine epmuElossP(vmin, vmax, ans)
!     ************************************
      implicit none
!             E * dsigma/dv
      real*8 vmin, vmax, ans
!
!      Int(v=vmin  to vmax) of v*ds/dv of  pair function
!      result is in mb.
!   Note:      v* ds/dv = epmuPrS (=sum of epmudsdv)
!
      external  epmuPrS
      real*8    epmuPrS
      real*8    v1, v2, ans1, error, eps
      integer icon

      data eps/1.d-9/


      v1 =  vmin
      ans  = 0.
      v2 = vmax
!      do while (v2 .ne.  vmax)
!         v2 =  min(v1*10.d0**0.5, vmax)
!         call k16pGausLeg2(epmuPrS, v1, v2, 16,  ans1)
      if(vmax > vmin )  then
          call kdexpIntFb(epmuPrS, vmin, vmax, eps,  ans1,
     *         error, icon)
         ans  =  ans +  ans1
      endif
!         v1=v2 
!      enddo
      end
!     *********************
      real*8 function epmuPrSdsdv(v)
      implicit none
#include "Zmedia.h"
#include "ZmuBPNgene.h"
!
!        sum of dsigma/dv  of a target in mb
!
      real*8 v
      real*8 epmuPrS

         
      epmuPrSdsdv = epmuPrS(v)/v

      end
