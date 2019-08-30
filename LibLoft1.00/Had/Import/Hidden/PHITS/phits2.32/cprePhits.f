      subroutine cprePhits
      implicit none
!        this is init qmd before qdm00 is called;
!        do the business which is done in read00 of PHITS
!   
      real(8):: bplus
      integer::  icrhi, ijudg, imadj, iqmax
      common /crshi/  bplus, icrhi, ijudg, imadj, iqmax
      integer:: ielas,icasc,iqstep,lvlopt,igamma
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      integer:: icugn
      
      real(8):: dtg, fric
      common /cugnon/ icugn
      real(8)::andt
      integer:: jevap,npidk
      common /bparm/  andt,jevap,npidk

      ijudg = 1  ! = mstz( 56 )      
      imadj = 1
      icrhi  = 1 ! =  mstz( 54  ) (1--> nasa xsection; 0-->Shen's Xs)
      iqmax = 150
      bplus = 0.  ! parz(168)

      icasc = 1  ! mstz(9) + 1
      iqstep = 6 ! mstz(6)
      ielas = 2  ! mstz(8)
      lvlopt = 3 ! mstz(11) 
      igamma = 0 ! mstz(12)
      icugn = 1  ! mstz(13)

      andt  = 0.  ! parz(22)
      jevap = 3 ! mstz(6)
      
      call cprePhits2

      end subroutine

      subroutine  cprePhits2
      implicit real*8 (a-h,o-z)
#include "param01.inc"      
      common /input1/ mstq1(mxpa1), parq1(mxpa1)
!           ntry1  # dist is prop exp(-x/10.712).
!       so x=128  makes prob. be ~ 1/1000.
!       over this, discard event
!       k+ 1 GeV.
!       128:    12230:   52 min  no Error
!        20:     5000:    3:56   576 Errors
!        50:     5000:   21:03   no Error
!        40:      //     18:45   //
!        30:      //     20:28   15 Errors
      mstq1(91) = 50   ! mntry (default is 1000)
      end subroutine

