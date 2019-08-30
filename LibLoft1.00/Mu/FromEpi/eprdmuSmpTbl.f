      subroutine epRdmuPrSTbl(io, cnst, tbl)
      implicit none
#include "ZbpSample.h"
#include "ZbpTbl.h"
!     read muon  pair ceation sampling table.
!     This table must be read first (earlier than 
!     Br, N).  If the table is not found,
!     muon energy loss is only by dE/dx due to
!     ionization.
!
      integer io
       type(bpTbl):: tbl   !    must be media.tbl
       type(SmpCnst)::  cnst  ! must be media.cnst
!

      call epRdmuPrCnst(io, cnst)
!            total x-sec. v> vmin
      if(cnst%muPrTXT .gt. mxMuPrTX) then
         call cerrorMsg(
     *      'mu pair cre. total X-sec table too large',0)
      else
         write(0,*) 'muon Pr cnst has been read'
      endif
!        total x-sec
      call epRd1dTbl(io,  tbl%MuPrTX, cnst%muPrTXT)
      write(0,*) ' pr total xs read'
!          dE/dx/E (/(g/cm2)
      call epRd1dTbl(io,  tbl%MuPrdEdx0, cnst%muPrTXT)
      write(0,*) ' pr dE/dx/E 0 read'
!          dE/dx/E (/(g/cm2))
      call epRd1dTbl(io,  tbl%MuPrdEdxt, cnst%muPrTXT)
      write(0,*) ' pr dE/dx/E t read'
!

!           sampling table
!      if( cnst.muPrUsize * cnst.muPrEsize .gt. mxMuPrTbl) then
!         call cerrorMsg(
!     *     'muon pair cre. sampling table too large', 0)
!      endif
!      call  epRd2dTbl(io,  tbl.MuPrTbl, 
!     *     cnst.muPrUsize, cnst.muPrEsize)
      end

      subroutine epRdmuBrSTbl(io, cnst, tbl)
      implicit none
#include "ZbpSample.h"
#include "ZbpTbl.h"

!
      integer io
       type(bpTbl):: tbl   !    must be media.tbl
       type(SmpCnst)::  cnst  ! must be media.cnst

!

      call epRdmuBrCnst(io, cnst)

!            total x-sec. v> vmin
      if(cnst%muBrTXT .gt. mxMuBrTX) then
         call cerrorMsg(
     *   'mu pair cre. total X-sec table too large',0)
      endif
!        total x-sec
      call epRd1dTbl(io,  tbl%MuBrTX, cnst%muBrTXT)
!        dE/dx/E (/(g/cm2)
      call epRd1dTbl(io,  tbl%MuBrdEdx0, cnst%muBrTXT)
!        dE/dx/E (/(g/cm2))
      call epRd1dTbl(io,  tbl%MuBrdEdxt, cnst%muBrTXT)
!
!         sampling table
!      if( cnst.muBrUsize * cnst.muBrEsize .gt. mxMuBrTbl) then
!         call cerrorMsg(
!     *   'muon pair cre. sampling table too large', 0)
!      endif
!      call  epRd2dTbl(io,  tbl.MuBrTbl, 
!     *       cnst.muBrUsize, cnst.muBrEsize)
!
      end
      subroutine epRdmuNSTbl(io, cnst, tbl)
      implicit none
#include "ZbpSample.h"
#include "ZbpTbl.h"

!
      integer io
       type(bpTbl):: tbl   !    must be media.tbl
       type(SmpCnst)::  cnst  ! must be media.cnst

!

      call epRdmuNCnst(io, cnst)

!            total x-sec. v> vmin
      if(cnst%muNTXT .gt. mxMuNTX) then
!///////
!         
!         write(0, *)  'cnst.muNTXT=', cnst.muNTXT, ' mx=', mxMuNTX
!////////////
         call cerrorMsg(
     *    'mu Nuc%int Total X-sec table too large',0)
      endif
!        total x-sec
      call epRd1dTbl(io,  tbl%MuNTX, cnst%muNTXT)
!        energy dependence
      read(io, *) cnst%muNpwtx
!        dE/dx/E (/(g/cm2)
      call epRd1dTbl(io,  tbl%MuNdEdx0, cnst%muNTXT)
!         energy dependence
      read(io, *) cnst%muNpwdEdx0
!        dE/dx/E (/(g/cm2))
      call epRd1dTbl(io,  tbl%MuNdEdxt, cnst%muNTXT)
!         energy dependence
      read(io, *) cnst%muNpwdEdxt
!
!         sampling table
      if( cnst%muNUsize * cnst%muNEsize .gt. mxMuNTbl) then
         call cerrorMsg(
     *   'muon nuc. int sampling table too large', 0)
      endif
      call  epRd2dTbl(io,  tbl%MuNTbl, 
     *       cnst%muNUsize, cnst%muNEsize)

      end

      subroutine epRdmuPrCnst(io, cnst)
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
     *  cnst%muPrVmin, cnst%muPrdETX, cnst%muPrdE, cnst%muPrEmin,
     *  cnst%muPrEmax, cnst%muPrdU, cnst%muPrUsize, cnst%muPrEsize,
     *  cnst%muPrTXT, cnst%muPrEmax1

      read(io, *)
      end

      subroutine epRdmuBrCnst(io, cnst)
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
     *   cnst%muBrVmin, cnst%muBrdETX, cnst%muBrdE, cnst%muBrEmin,
     *   cnst%muBrEmax,  cnst%muBrdU, cnst%muBrUsize, cnst%muBrEsize,
     *   cnst%muBrTXT, cnst%muBrEmax1

      read(io, *)
      end


      subroutine epRdmuNCnst(io, cnst)
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
!/////
!         call cerrorMsg(msg, 1)
!////
      read(io, *)  
     *   cnst%muNVmin,  cnst%muNdETX, cnst%muNdE,  cnst%muNEmin,
     *   cnst%muNEmax,  cnst%muNdU, cnst%muNUsize, cnst%muNEsize,
     *   cnst%muNTXT, cnst%muNEmax1
!/////
!      write(0, *)  
!     *   cnst.muNVmin,  cnst.muNdETX, cnst.muNdE,  cnst.muNEmin,
!     *   cnst.muNEmax,  cnst.muNdU, cnst.muNUsize, cnst.muNEsize,
!     *   cnst.muNTXT, cnst.muNEmax1
!////
      read(io, *)
      end
