      subroutine epRdBrSTbl1(io, cnst, tbl)
      implicit none
#include "ZbpSample.h"
#include "ZbpTbl.h"

!
      integer io
       type(bpTbl):: tbl   !    must be media.tbl
       type(SmpCnst)::  cnst  ! must be media.cnst

!

      call epRdBrCnst(io, cnst)

!            total x-sec.
      if(cnst%BremTXTL .gt. mxBrTXL) then
         write(0,*) 'cnst%BremTXTL=', cnst%BremTXTL, 
     *     " >  mxBrTXL=",mxBrTXL
         stop 1111
      endif
      call epRd1dTbl(io,  tbl%BrTXL, cnst%BremTXTL)

!           sampling table in A region
      if( cnst%BremUszLA * cnst%BremEsize .gt. mxBrTblLA) then
         write(0,*) "cnst%BremUszLA*cnst%BremEsize=",
     *     cnst%BremUszLA,"*", cnst%BremEsize, "=",
     *     cnst%BremUszLA*cnst%BremEsize,
     *     "> mxBrTblLA=",mxBrTblLA
         stop 1112
      endif
      call  epRd2dTbl(io,  tbl%BrSTLA, 
     *                    cnst%BremUszLA, cnst%BremEsize)
!           sampling table in B region
      if( cnst%BremUszLB * cnst%BremEsize .gt. mxBrTblLB) then
         write(0,*) 'cnst%BremUszLB * cnst%BremEsize=',
     *   cnst%BremUszLB,"*",  cnst%BremEsize,"=", 
     *    cnst%BremUszLB * cnst%BremEsize,
     *  " >  mxBrTblLB= ",mxBrTblLB
         stop 1113
      endif

      call  epRd2dTbl(io,  tbl%BrSTLB,
     *             cnst%BremUszLB, cnst%BremEsize)
      end

!     ****************************************
      subroutine epRd1dTbl(io,  tbl, size)
!     ****************************************
      implicit none
#include "ZepManager.h"

      integer size, io
      real*8  tbl(size,2), dummy

      character*160  msg
      character(40)::almostdummy(2)
      integer i, nf

      read(io, '(a)') msg
      if(MsgLevel-1 .ge. 0) then
         call cerrorMsg(msg, 1)
      endif
      read(io, '(a)') msg
      if(MsgLevel-1 .ge. 0) then
         call cerrorMsg(msg, 1)
      endif
      msg = ' '
      read(io,'(a)') msg
!      call kgetField(trim(msg), msg(80:), 2, nf)
      call kgetField(trim(msg), almostdummy, 2, nf)
      
      do i = 1, size
         if(nf == 2 ) then
            read(msg, *) dummy, tbl(i,1)
            tbl(i,2) = 0
         elseif( nf == 3 ) then
            read(msg, *) dummy, tbl(i,1),tbl(i,2)
         else
            write(0,*) trim(msg)
            write(0,*)
     *      'above data:cross section table has fields > 3=',nf
            write(0,*) ' in  epRd1dTbl'
            stop
         endif
         if( i < size) then
            read(io, '(a)') msg
         endif
      enddo
      read(io,*)
      end

!     **********************************************
      subroutine epRd2dTbl(io, bla, sizeu, sizee)
      implicit none
#include "ZepManager.h"

      integer sizeu, sizee, io
      real*8  bla(sizeu,  sizee)

      character*160 msg
      integer iu, i

      read(io,'(a)') msg
      if(MsgLevel-1 .ge. 0) then
         call cerrorMsg(msg, 1)
      endif

      do i = 1, sizee
         read(io, '(7f11.5)') ( bla(iu, i), iu = 1, sizeu )
         read(io, *) 
      enddo
      read(io,*)
      end

      subroutine epRdBrCnst(io, cnst)
      implicit none
#include  "ZbpSample.h"
#include  "ZepManager.h"
       type(SmpCnst)::  cnst

      integer io
      character*160 msg
 
      read(io, '(a)') msg
      if(MsgLevel-1 .ge. 0) then
         call cerrorMsg(msg, 1)
      endif
      read(io, *)  cnst%CompScrE,
     * cnst%BremEgmin, cnst%BremEemin, cnst%BremLEemin,
     *  cnst%BremEeminLPM, cnst%BrScrE,
     *  cnst%BremUminLA, cnst%BremUmaxLA, cnst%BremTXTL,
     *  cnst%BremdULA, cnst%BremdETXL, cnst%BremdEL, cnst%BremUminLB,
     *  cnst%BremUmaxLB, cnst%BremdULB,  cnst%BremEsize,
     *  cnst%BremUszLA, cnst%BremUszLB
      read(io, *)
      end
      subroutine epRdBrSTblH(io, cnst, tbl)
!          read brems sampling table in LPM region
      implicit none
#include "ZbpSample.h"
#include "ZbpTbl.h"

!
      integer io
       type(bpTbl):: tbl   !    must be media.tbl
       type(SmpCnst)::  cnst  ! must be media.cnst

!

      call epRdBrCnstH(io, cnst)

!            total x-sec.
      if(cnst%BrneH .gt. mxBrTXH) then
         write(0,*)" LPM brem table:  cnst%BrneH=",cnst%BrneH,
     *   "> mxBrTXH =",mxBrTXH 
         stop 2221
      endif
      call epRd1dTbl(io,  tbl%BrTXH, cnst%BrneH)

!           sampling table in A region
      if( cnst%Brnu1H * cnst%BrneH2 .gt. mxBrTblHA) then
         write(0,*) " LPM brem table A: cnst%Brnu1H * cnst%BrneH2=",
     *   cnst%Brnu1H, "*",  cnst%BrneH2, "=", cnst%Brnu1H*cnst%BrneH2,
     *   "> mxBrTblHA=" , mxBrTblHA
         stop 2222
      endif
      call  epRd2dTbl(io,  tbl%BrSTHA, 
     *                    cnst%Brnu1H, cnst%BrneH2)

!           sampling table in B region
      if( cnst%Brnu2H * cnst%BrneH2 .gt. mxBrTblHB) then
         write(0,*) " LPM brem table B: cnst%Brnu2H*cnst%BrneH2=",
     *     cnst%Brnu2H, "*", cnst%BrneH2, "=",cnst%Brnu2H*cnst%BrneH2,
     *   ">mxBrTblHB=",mxBrTblHB
         stop 2223
      endif
      call  epRd2dTbl(io,  tbl%BrSTHB,
     *             cnst%Brnu2H, cnst%BrneH2)
      end
      subroutine epRdBrCnstH(io, cnst)
      implicit none
#include  "ZbpSample.h"
#include  "ZepManager.h"

       type(SmpCnst)::  cnst

      integer io
      character*160 msg
 
      read(io, '(a)') msg
      if(MsgLevel-1 .ge. 0) then
         call cerrorMsg(msg, 1)
      endif

      read(io, *)   cnst%BrEgminH, cnst%BrEe1H, cnst%BrLEe1H,
     *   cnst%BrneH, cnst%BrdU1H, cnst%BrdEH,
     *   cnst%BrEe2H, cnst%BrdU1H,  cnst%BrU1H,
     *   cnst%BrU2H,  cnst%Brnu1H,  cnst%BrU3H,
     *   cnst%BrU4H,  cnst%Brnu2H,  cnst%BrdVU2H,
     *   cnst%BrdU2H, cnst%BrneH2,  cnst%BrdEH2,
     *   cnst%BrEe2H2, cnst%BrPow
      read(io, *)
      end
      subroutine epRdBrSTblS(io, cnst, tbl)
      implicit none
#include "ZbpSample.h"
#include "ZbpTbl.h"

!
      integer io
       type(bpTbl):: tbl   !    must be media.tbl
       type(SmpCnst)::  cnst  ! must be media.cnst

!

      call epRdBrCnstS(io, cnst)

!            total x-sec.
      if(cnst%BrTXTS .gt. mxBrTXS) then
         write(0,*) " cnst%BrTXTS=",cnst%BrTXTS, ">mxBrTXS=", 
     *      mxBrTXS
         stop 3331
      endif
      call epRd1dTbl(io,  tbl%BrTXS, cnst%BrTXTS)

!           sampling table in A region
      if( cnst%BrUszSA * cnst%BrES .gt. mxBrTblSA) then
         write(0,*)
     *     "cnst%BrUszSA*cnst%BrES=", cnst%BrUszSA,"*",cnst%BrES,
     *    " =",  cnst%BrUszSA * cnst%BrES, " > mxBrTblSA=",
     *     mxBrTblSA
         stop 3332
      endif
      call  epRd2dTbl(io,  tbl%BrSTSA, 
     *                    cnst%BrUszSA, cnst%BrES)
!           sampling table in B region
      if( cnst%BrUszSB * cnst%BrES .gt. mxBrTblSB) then
         write(0,*)
     *     "cnst%BrUszSB*cnst%BrES=", cnst%BrUszSB,"*",cnst%BrES,
     *    " =",  cnst%BrUszSB * cnst%BrES, " > mxBrTblSB=",
     *     mxBrTblSB
         stop 3333
      endif
      call  epRd2dTbl(io,  tbl%BrSTSB,
     *             cnst%BrUszSB, cnst%BrES)


!           upper Selzer region
!            total x-sec.
      if(cnst%BrTXTS2 .gt. mxBrTXS2) then
         write(0,*) " cnst%BrTXTS2=",cnst%BrTXTS2, ">mxBrTXS2=", 
     *      mxBrTXS2
         stop 3331
      endif
      call epRd1dTbl(io,  tbl%BrTXS2, cnst%BrTXTS2)
      write(0,*) ' xsection read'

!           sampling table in A region
      if( cnst%BrUszSA2 * cnst%BrES2 .gt. mxBrTblSA2) then
         write(0,*)
     *     "cnst%BrUszSA2*cnst%BrES2=", cnst%BrUszSA2,"*",cnst%BrES2,
     *    " =",  cnst%BrUszSA2 * cnst%BrES2, " > mxBrTblSA2=",
     *     mxBrTblSA2
         stop 3332
      endif
      call  epRd2dTbl(io,  tbl%BrSTSA2, 
     *                    cnst%BrUszSA2, cnst%BrES2)
      write(0,*) ' tbl read'
!           sampling table in B region
      if( cnst%BrUszSB2 * cnst%BrES2 .gt. mxBrTblSB2) then
         write(0,*)
     *     "cnst%BrUszSB2*cnst%BrES2=", cnst%BrUszSB2,"*",cnst%BrES2,
     *    " =",  cnst%BrUszSB2 * cnst%BrES2, " > mxBrTblSB2=",
     *     mxBrTblSB2
         stop 3333
      endif
      call  epRd2dTbl(io,  tbl%BrSTSB2,
     *             cnst%BrUszSB2, cnst%BrES2)
      write(0,*) ' tbl read'
      end
      subroutine epRdBrCnstS(io, cnst)
      implicit none
#include  "ZbpSample.h"
#include  "ZepManager.h"
       type(SmpCnst)::  cnst

      integer io
      character*160 msg
 
      read(io, '(a)') msg
      if(MsgLevel -1 .ge. 0) then
         call cerrorMsg(msg, 1)
      endif
      read(io, *)  
     *  cnst%BrEeminS, cnst%BrEgminS,  cnst%BrLEeminS,
     *  cnst%BrEemaxS, 
     *  cnst%BrUminSA, cnst%BrUmaxSA, cnst%BrTXTS,
     *  cnst%BrdUSA, cnst%BrdETXS, cnst%BrdES, cnst%BrUminSB,
     *  cnst%BrUmaxSB, cnst%BrdUSB,  cnst%BrES,
     *  cnst%BrUszSA, cnst%BrUszSB, cnst%how,
     *  cnst%NormS, cnst%NormPS, cnst%NormCS, cnst%NormSH
!    The last 5 are specail for brems normalization 
      read(io, '(a)') msg
      if(MsgLevel -1 .ge. 0) then
         call cerrorMsg(msg, 1)
      endif
      read(io, *)  
     *  cnst%BrEeminS2, cnst%BrEgminS2,  cnst%BrLEeminS2,
     *  cnst%BrEemaxS2, 
     *  cnst%BrUminSA2, cnst%BrUmaxSA2, cnst%BrTXTS2,
     *  cnst%BrdUSA2, cnst%BrdETXS2, cnst%BrdES2, cnst%BrUminSB2,
     *  cnst%BrUmaxSB2, cnst%BrdUSB2,  cnst%BrES2,
     *  cnst%BrUszSA2, cnst%BrUszSB2
      read(io, *)
      end
      subroutine epRdPrSTbl1(io, cnst, tbl)
!        read  pair sampling table
      implicit none
#include "ZbpSample.h"
#include "ZbpTbl.h"

!
      integer io
       type(bpTbl):: tbl   !    must be media.tbl
       type(SmpCnst)::  cnst  ! must be media.cnst

!

      call epRdPrCnst(io, cnst)

!            total x-sec.
      if(cnst%PairTXTL .gt. mxPrTXL) then
         write(0,*) "cnst%PairTXTL=", cnst%PairTXTL,
     *  ">  mxPrTXL=",  mxPrTXL
         stop 4441
      endif
      call epRd1dTbl(io,  tbl%PrTXL, cnst%PairTXTL)

!           sampling table in A region
      if( cnst%PairUszLA * cnst%PairEsize .gt. mxPrTblLA) then
         write(0,*) "cnst%PairUszLA*cnst%PairEsize=",
     *    cnst%PairUszLA, "*", cnst%PairEsize, "=",
     *    cnst%PairUszLA * cnst%PairEsize,
     *   "> mxPrTblLA=",mxPrTblLA
         stop 4442
      endif
      call  epRd2dTbl(io,  tbl%PrSTLA, 
     *                    cnst%PairUszLA, cnst%PairEsize)
!           sampling table in B region
      if( cnst%PairUszLB * cnst%PairEsize .gt. mxPrTblLB) then
         write(0,*) "cnst%PairUszLB * cnst%PairEsize=",
     *   cnst%PairUszLB,"*", cnst%PairEsize, "=",
     *   cnst%PairUszLB*cnst%PairEsize, ">",
     *   "mxPrTblLB=",  mxPrTblLB
         stop 4443
      endif
      call  epRd2dTbl(io,  tbl%PrSTLB,
     *             cnst%PairUszLB, cnst%PairEsize)
      end

      subroutine epRdPrCnst(io, cnst)
      implicit none
#include  "ZbpSample.h"
#include  "ZepManager.h" 
       type(SmpCnst)::  cnst

      integer io
      character*160 msg
 
      read(io, '(a)') msg
      if(MsgLevel-1 .ge. 0) then
         call cerrorMsg(msg, 1)
      endif
      read(io, *)  cnst%PairEgmin, cnst%PairLEgmin,
     * cnst%PairNonSc, cnst%PrScrE,
     * cnst%PairEgmaxL,  cnst%PairTXTL, cnst%PairEsize,
     * cnst%PairUminLA,  cnst%PairUmaxLA, cnst%PairUszLA,
     * cnst%PairdULA, cnst%PairdETXL,  cnst%PairUminLB,
     * cnst%PairUmaxLB, cnst%PairUszLB, cnst%PairdULB,
     * cnst%PairdELA, cnst%PairdELB
      read(io, *)
      end
      subroutine epRdPrSTblH(io, cnst, tbl)
!        read  pair sampling table with LPM
      implicit none
#include "ZbpSample.h"
#include "ZbpTbl.h"

!
      integer io
       type(bpTbl):: tbl   !    must be media.tbl
       type(SmpCnst)::  cnst  ! must be media.cnst

!

      call epRdPrCnstH(io, cnst)

!            total x-sec.
      if(cnst%PrneH .gt. mxPrTXH) then
         write(0,*) "LPM Pair: cnst%PrneH=",cnst%PrneH, " >",
     *   " mxPrTXH=", mxPrTXH
         stop 5551
      endif
      call epRd1dTbl(io,  tbl%PrTXH, cnst%PrneH)

      if( cnst%Prnu1H * cnst%PrneH .gt. mxPrTblH) then
         write(0,*) "cnst%Prnu1H * cnst%PrneH=",
     *     cnst%Prnu1H,"*", cnst%PrneH, "=",
     *     cnst%Prnu1H*cnst%PrneH, ">",
     *     mxPrTblH
         stop 5552
      endif
      call  epRd2dTbl(io,  tbl%PrSTH,
     *             cnst%Prnu1H, cnst%PrneH)
      end

      subroutine epRdPrCnstH(io, cnst)
      implicit none
#include  "ZbpSample.h"
#include  "ZepManager.h"
       type(SmpCnst)::  cnst

      integer io
      character*160 msg
 
      read(io, '(a)') msg
      if(MsgLevel -1 .ge. 0) then
         call cerrorMsg(msg, 1)
      endif
      read(io, *)  
     *  cnst%PrEg1H, cnst%PrneH, cnst%PrdU1H, cnst%PrdEH,
     *  cnst%PrU1H,  cnst%PrU2H, cnst%Prnu1H, cnst%PrLEg1H,
     *  cnst%PrEg2H
      read(io, *)
      end
