      subroutine cPrLSampP(Eg, prob)
      implicit none
#include "ZbpCnst.h"
#include "ZbpTable.h"
      real*8 Eg
      real*8 prob  ! output probability of Pair / X0


      real*8 ale


      if(Eg .lt. PairEgmin) then
         prob= 1.d-40
      else
         ale=log10(Eg)
         call kintp3(PrTXL,
     *   1, PairTXTL, PairLEgmin,
     *   PairdETXL, ale, prob) 
      endif
      end
!     ************
      subroutine cPrLSampE(Eg, Ee)
!     ************
!          samples higher energy pair electron
      implicit none
#include "ZbpCnst.h"
#include "ZbpTable.h"

      real*8 Ee,  Eg

      real*8 u, ale, us, ans, ex

      call rndc(u)
      ale = log10(Eg)
      if(u .gt. PairUminLA) then
!          region A
         call k4ptdi(PrSTLA, 
     *        PairUszLA, 
     *        PairEsize,
     *        PairUszLA, 
     *        PairUminLA,
     *        PairLEgmin,
     *        PairdULA,
     *        PairdELA, u,  ale,  ans)  
         Ee = (ans*(1.-u) + 0.5d0)*Eg
      else
!         region B
         us = u**0.25d0
         ex = sqrt(ale - PairLEgmin)
         
         call k4ptdi(PrSTLB, 
     *        PairUszLB, 
     *        PairEsize,
     *        PairUszLB, 
     *        0.d0,
     *        0.d0,
     *        PairdULB,
     *        PairdELB, us,  ex,  ans)  
         Ee = ans* Eg
      endif
      end

