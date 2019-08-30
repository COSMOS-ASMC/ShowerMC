      subroutine cBrLSampP(Ee, prob)
      implicit none
#include "ZbpCnst.h"
#include "ZbpTable.h"

      real*8 Ee
      real*8 prob  ! output probability of Brems / X0

      real*8 ale
!        Ee
!           BremEemin  BremEeCS  BremEemaxL  BremEemaxH
!     ----------------------------------------------------
!                |         |            |          |
!     no brems       p.s        c.s        lpm       lpm by
!                                                    rejection
!

      if(Ee .le. BremEemin) then
         prob= 1.d-40
      else
         ale=log10(Ee)
         call kintp3(BrTXL,
     *   1, BremTXTL, BremLEemin,
     *   BremdETXL, ale, prob) 
      endif
      end
!     ************
      subroutine cBrLSampE(Ee, Eg)
!     ************
      implicit none
#include "ZbpCnst.h"
#include "ZbpTable.h"

#include "Zmass.h"

      real*8 Ee,  Eg

      real*8 u, ale, us, ans

      call rndc(u)
      ale = log10(Ee)
      if(u .gt. BremUminLA) then
!          region A

         call k4ptdi(BrSTLA, 
     *        BremUszLA, 
     *        BremEsize,
     *        BremUszLA, 
     *        BremUminLA,
     *        BremLEemin,
     *        BremdULA,
     *        BremdEL, u,  ale,  ans)  
         Eg= exp( ans*(1.-u))*BremEgmin
      else
!         region B
         us = u**0.25d0

         call k4ptdi(BrSTLB, 
     *        BremUszLB, 
     *        BremEsize,
     *        BremUszLB, 
     *        0.d0,
     *        BremLEemin,
     *        BremdULB,
     *        BremdEL, us,  ale,  ans)  
         Eg = exp(-ans*u)*(Ee - masele)
      endif
      end

