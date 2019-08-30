#define USELogTabOnlyLowE
!            !#undef USELogTabOnlyLowE
      subroutine epPrLSampP(xmedia, Eg, prob)
      implicit none
!     #include "Zmedia.h"
!     #include "ZepTrackv.h"
#include "ZmediaLoft.h"      
!
#ifdef  USELogTabOnlyLowE      
      real(8),parameter:: UseLogBelow=1.20d-3 ! below this, use log table
      integer,parameter:: maxIdxLog=6 !make log table below this index
!          see next table for index.      
!    idx    Eg  Gev     prob/r.l (for  W)
!     1  1.02300E-03  1.36858E-10
!     2  1.22149E-03  6.69279E-04
!     3  1.45848E-03  4.34835E-03
!     4  1.74146E-03  1.21478E-02
!     5  2.07935E-03  2.41705E-02
!     6  2.48279E-03  4.01101E-02
!     7  2.96451E-03  5.95280E-02
!     8  3.53970E-03  8.19728E-02
!     9  4.22648E-03  0.10703
!     10 5.04652E-03  0.13431    
!
#endif      
      type(epmedia):: xmedia
      logical,save::first(Maxmedia)=.true.
      real*8 Eg
      real*8 prob  ! output probability of Pair / X0
!     real*8,allocatable::logtbl(mxPrTXL,:) ! log table  not used now  2018/Dec/7
      real*8,allocatable::logtbl(:,:) ! log table 
      real*8 ale
      integer i
      save  logtbl
      logical,save::onlyonce=.true.

      if( onlyonce ) then
#ifdef  USELogTabOnlyLowE         
         allocate( logtbl(maxIdxLog, NoOfMedia)) ! log table
#else
         allocate( logtbl(mxPrTXL, NoOfMedia)) ! log table for entire region
#endif         
         onlyonce =.false.
      endif
      if( first(MediaNo))  then
!     take log of the total crross-section
#ifdef  USELogTabOnlyLowE         
         do i = 1, maxIdxLog
#else
         do i = 1, xmedia%cnst%PairTXTL         
#endif                     
            logtbl(i,MediaNo)=log( xmedia%tbl%PrTXL(i) )
         enddo
         first(MediaNo) = .false.
      endif
      if(Eg .lt. xmedia%cnst%PairEgmin) then
         prob= 1.d-40
      else
         ale=log10(Eg)
!     without using log table, accuracy is OK.
!     but for very small Eg, prob < 10^-8 and sometimes it becomes
!     ~  -10^-4 and leads to boundary error. This error can be
!     aovided with ElowerBndPair = 1.2 MeV in epics file.
!     This version is bit slower than the no-log version by ~7% or so
!     but can be use ElowerBndPair = 1.022 MeV . if this is 1.2MeV
!     the value becomes ~2%
#ifdef  USELogTabOnlyLowE                  
         if(Eg < UseLogBelow ) then
               ! PairTXTL need not be maxIdxLog 
            call kintp3(logtbl(1,MediaNo),
     *      1, xmedia%cnst%PairTXTL, xmedia%cnst%PairLEgmin,
     *      xmedia%cnst%PairdETXL, ale, prob) 
            prob = exp(prob)
         else
            call kintp3( xmedia%tbl%PrTXL,
     *      1, xmedia%cnst%PairTXTL, xmedia%cnst%PairLEgmin,
     *      xmedia%cnst%PairdETXL, ale, prob) 
         endif
#else
         call kintp3(logtbl(1,MediaNo),
     *    1, xmedia%cnst%PairTXTL, xmedia%cnst%PairLEgmin,
     *    xmedia%cnst%PairdETXL, ale, prob) 
          prob = exp(prob)
#endif         
      endif
      end
!     ************
      subroutine epPrLSampE(media, Eg, Ee)
!     ************
!          samples higher energy pair electron
      implicit none
#include "Zmedia.h"
       type(epmedia):: media
      real*8 Ee,  Eg

      real*8 u, ale, us, ans, ex

      call rndc(u)
      ale = log10(Eg)
      if(u .gt. media%cnst%PairUminLA) then
!          region A
         call k4ptdi(media%tbl%PrSTLA, 
     *        media%cnst%PairUszLA, 
     *        media%cnst%PairEsize,
     *        media%cnst%PairUszLA, 
     *        media%cnst%PairUminLA,
     *        media%cnst%PairLEgmin,
     *        media%cnst%PairdULA,
     *        media%cnst%PairdELA, u,  ale,  ans)  
         Ee = (ans*(1.-u) + 0.5d0)*Eg
      else
!         region B
         us = u**0.25d0
         ex = sqrt(ale - media%cnst%PairLEgmin)
         
         call k4ptdi(media%tbl%PrSTLB, 
     *        media%cnst%PairUszLB, 
     *        media%cnst%PairEsize,
     *        media%cnst%PairUszLB, 
     *        0.d0,
     *        0.d0,
     *        media%cnst%PairdULB,
     *        media%cnst%PairdELB, us,  ex,  ans)  
         Ee = ans* Eg
      endif
      end

