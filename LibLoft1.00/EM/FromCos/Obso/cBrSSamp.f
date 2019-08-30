      subroutine cBrSSampP(Ee, prob)
!          brems prob. sampling by Seltzer data.
      implicit none
#include "ZbpCnst.h"
#include "ZbpTable.h"

      real*8 Ee
      real*8 prob  ! output probability of Brems / X0

      real*8 ale


      if(Ee .le. BrEeminS) then
         prob= 1.d-40
      else
         ale=log10(Ee)
         call kintp3(BrTXS,
     *   1, BrTXTS, BrLEeminS,
     *   BrdETXS, ale, prob) 
      endif
      end
!     ************
      subroutine cBrSSampE( Eein, Eg)
!     ************
!         brems energy by Seltzer
      implicit none
#include "Zmass.h"
#include "ZbpCnst.h"
#include "ZbpTable.h"

      real*8 Eein,  Eg

      real*8 Ee
      real*8 u, ale, us, ans
      integer count
      logical ok
      
      count = 0
      if(Eein-masele .le. BrEgminS ) then
!         simply avoid loop;  ( Eek < 10 keV)
         call rndc(u)
         Eg =( Eein-masele ) * u
         return  ! ************
      elseif(Eein .lt. BrEeminS) then
!         this and above  happne when electron loses energy
!         after sampling of processes 
         Ee =  BrEeminS
      else
         Ee = Eein
      endif
      ale = log10(Ee)

 10   continue
      call rndc(u)
      if(u .gt. BrUminSA) then
!          region A

         call k4ptdi(BrSTSA, 
     *        BrUszSA, 
     *        BrES,
     *        BrUszSA, 
     *        BrUminSA,
     *        BrLEeminS,
     *        BrdUSA,
     *        BrdES, u,  ale,  ans)  
         Eg= exp( ans*(1.-u))*BrEgminS
      else
!         region B
         us = u**0.25d0

         call k4ptdi(BrSTSB, 
     *        BrUszSB, 
     *        BrES,
     *        BrUszSB, 
     *        0.d0,
     *        BrLEeminS,
     *        BrdUSB,
     *        BrdES, us,  ale,  ans)  
         Eg = exp(-ans*u)*(Ee-masele)
!         Eg = exp(-ans*u)*Ee
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
      end

