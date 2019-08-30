      subroutine epBrSSampP(media, Ee, prob)
!     brems prob. sampling by Seltzer data.
      use modEMcontrol
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
!  #include "ZepTrackp.h"

       type(epmedia):: media
      real*8 Ee
      real*8 prob  ! output probability of Brems / X0

      real*8 ale, Ek


      if(Ee .le. media%cnst%BrEeminS) then
         prob= 1.d-40
      elseif( Ee <= media%cnst%BrEemaxS) then
         Ek = Ee - masele
         ale=log10(Ek)

         call kintp3(media%tbl%BrTXS(:,1),
     *   1, media%cnst%BrTXTS, media%cnst%BrLEeminS,
     *   media%cnst%BrdETXS, ale, prob) 
      elseif(Ee <=  media%cnst%BrEemaxS2 ) then
         Ek = Ee - masele
         ale=log10(Ek)
         call kintp3(media%tbl%BrTXS2(:,1), 
     *   1, media%cnst%BrTXTS2, media%cnst%BrLEeminS2,
     *   media%cnst%BrdETXS2, ale, prob) 
      else
         write(0,*) ' Ee=',Ee,' too large for Seltzer'
         stop
      endif
      if( HowNormBrems == -1 ) then
         ! nothing to do
      else
         if( HowNormBrems == 1 ) then
            prob =prob / media%cnst%NormS
         elseif( HowNormBrems == 0 ) then
            prob =prob / media%cnst%NormS
         else
            write(0,*) 'HowNormBrems=',HowNormBrems
            write(0,*) ' invalide in epBrSSampP'
            stop
         endif
      endif
      end
!     ************
      subroutine epBrSSampE(media, Eein, Eg)
!     ************
!         brems energy by Seltzer
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
       type(epmedia):: media
      real*8 Eein,  Eg

      real*8 Ee, Ek 

      real*8 u, ale, us, ans, vmax
      integer count
      logical ok
!&&&&&&&
      real*8 dummy(1)
      real*8 error
      data dummy/0.d0/
      save
!&&&&&&&&&&&


      Ek = Eein-masele
      if(Eein .lt. media%cnst%BrEeminS) then
!         this  happne when electron loses energy    
!         after sampling of processes
         Ee = media%cnst%BrEeminS
      else
         Ee = Eein
      endif
      if(Ee <= media%cnst%BrEemaxS) then
         call epBrSSampE1(media, Ee, Eg)
      elseif( Ee <= media%cnst%BrEemaxS2) then
         call epBrSSampE2(media, Ee, Eg)
      else
         write(0,*)
     *   ' Ee=',Ee, ' too large for Seltzer'
         stop
      endif
      end

      subroutine epBrSSampE1(media, Eein, Eg)
!         brems energy by Seltzer
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
       type(epmedia):: media
      real*8 Eein,  Eg

      real*8 Ee, Ek 

      real*8 u, ale, us, ans, vmax
      integer count
      logical ok
!&&&&&&&
      real*8 dummy(1)
      real*8 error
      data dummy/0.d0/
      save
!&&&&&&&&&&&
      Ee = Eein
      Ek = Ee - masele
      vmax = 1.-masele/Ee
      ale = log10(Ek)

      count = 0      
 10   continue
      call rndc(u)
      if(u .gt. media%cnst%BrUminSA) then
!          region A
!&&&&&&&&&&&&&
!         call k4ptdi(media.tbl.BrSTSA, 
!     *        media.cnst.BrUszSA, 
!     *        media.cnst.BrES,
!     *        media.cnst.BrUszSA, 
!     *        media.cnst.BrUminSA,
!     *        media.cnst.BrLEeminS,
!     *        media.cnst.BrdUSA,
!     *        media.cnst.BrdES, u,  ale,  ans)  
!         Eg= exp( ans*(1.-u))*media.cnst.BrEgminS
         call kpolintp2(media%cnst%BrUminSA, 1,  media%cnst%BrdUSA,
     *            media%cnst%BrLEeminS, 1, media%cnst%BrdES,
     *            media%tbl%BrSTSA,  media%cnst%BrUszSA, 
     *            media%cnst%BrUszSA,  media%cnst%BrES,
     *             3, 3,  u, ale, ans, error)
!&&&&&&&&
!         Eg= exp( ans*(1.-u))*media.cnst.BrEgminS*Ee
         Eg= exp( ans*(1.-u))*media%cnst%BrEgminS   !  BrEgminS is not ratio v9.135
      else
!         region B
         us = u**0.25d0
         call kpolintp2(dummy, 1, media%cnst%BrdUSB,  
     *            media%cnst%BrLEeminS, 1, media%cnst%BrdES,
     *            media%tbl%BrSTSB,  media%cnst%BrUszSB, 
     *            media%cnst%BrUszSB,  media%cnst%BrES,
     *             3, 3,  us, ale, ans, error)
         Eg = exp(-ans*u)*Ee*vmax
      endif

      ok = Eein-Eg .ge. masele

      if(.not. ok) then
         count  = count + 1
         if(count .gt. 100 ) then
            call rndc(u)
            Eg = (Eein-masele)*u
            return  !  ******
         endif
         goto 10
      endif
      end   subroutine epBrSSampE1

      subroutine epBrSSampE2(media, Eein, Eg)
!         brems energy by Seltzer, upper region
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
       type(epmedia):: media
      real*8 Eein,  Eg

      real*8 Ee, Ek 

      real*8 u, ale, us, ans, vmax
      integer count
      logical ok
!&&&&&&&
      real*8 dummy(1)
      real*8 error
      data dummy/0.d0/
      save
!&&&&&&&&&&&
      Ee= Eein
      Ek = Ee - masele
      vmax = 1.-masele/Ee
      ale = log10(Ek)

      count = 0      
 10   continue
      call rndc(u)
      if(u .gt. media%cnst%BrUminSA2) then
!          region A
         call kpolintp2(media%cnst%BrUminSA2, 1, 
     *             media%cnst%BrdUSA2,
     *            media%cnst%BrLEeminS2, 1, media%cnst%BrdES2,
     *            media%tbl%BrSTSA2,  media%cnst%BrUszSA2, 
     *            media%cnst%BrUszSA2,  media%cnst%BrES2,
     *             3, 3,  u, ale, ans, error)
!&&&&&&&&
         Eg= exp( ans*(1.-u))*media%cnst%BrEgminS2*Ee
      else
!         region B
         us = u**0.25d0
         call kpolintp2(dummy, 1, media%cnst%BrdUSB2,  
     *        media%cnst%BrLEeminS2, 1, media%cnst%BrdES2,
     *        media%tbl%BrSTSB2,  media%cnst%BrUszSB2, 
     *        media%cnst%BrUszSB2,  media%cnst%BrES2,
     *             3, 3,  us, ale, ans, error)
         Eg = exp(-ans*u)*Ee*vmax
      endif

      ok = Eein-Eg .ge. masele

      if(.not. ok) then
         count  = count + 1
         if(count .gt. 100 ) then
            call rndc(u)
            Eg = (Eein-masele)*u
            return  !  ******
         endif
         goto 10
      endif
      end subroutine epBrSSampE2
