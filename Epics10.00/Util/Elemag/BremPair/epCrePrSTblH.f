!     ****************************************************************
!     *                                                              *
!     *   create sampling table for pair
!     *            with landau effect at high energies               *
!     *                                                              *
!     ****************************************************************
      subroutine epCrePrSTblH(mediain, cnst)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

!
       type(epmedia):: mediain  !input.  media
       type(SmpCnst)::  cnst  !input.  must be media.cnst
!

!
      integer ws
      parameter (ws = 10000)
      real*8 work(ws), work2(ws)
!        cp mediain into common area
      media = mediain
!         print pair const. at high energies
      call epwtPrCnstH(cnst)
!            total x-sec.
      if(cnst%PrneH .gt. mxPrTXH) then
         call cerrorMsg('LPM pair total X-sec. table too large',0)
      endif
      call epCrePrTXTH(cnst, work, work2, cnst%PrneH)
!           sampling table
      if( cnst%Prnu1H*cnst%PrneH .gt. mxPrTblH) then
         call cerrorMsg('Pair table for LPM is too large', 0)
      endif
      call  epCrePrH(cnst, work, cnst%Prnu1H, cnst%PrneH)
      end
!     ****************************************
      subroutine epCrePrTXTH(cnst, erg,  tbl, size)
!     ****************************************
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

       type(SmpCnst)::  cnst  ! must be media.cnst
      integer size
      real*8  tbl(size), erg(size)

      real*8 E,  vmax,   tcp, de
      character*160 msg
      integer i


      write(msg,*) ' Matter=', media%name,
     *  ': Computing the total X-sec of Pair from Eg= ',
     *    cnst%PrEg1H, ' GeV with LPM effect'
      call cerrorMsg(msg, 1)

      E = cnst%PrEg1H*1.000001d0

      de = 10.**cnst%PrdEH
      do i = 1, size
         Egme = E/masele
         vmax=1.d0 - masele/ E
!         call eptotcpH(.5d0, vmax, tcp)
         call epPrgeneTX(0.5d0, vmax, tcp)
         tbl(i) = 2* tcp* media%mbtoPX0  ! prob. per radiation length.
!                                        ! don't use mbtoPX02. 
         erg(i) = E
         E = E * de
      enddo
      write(msg, *) 
     * 'Pair table with LPM has been made:E=', cnst%PrEg1H,
     *  ' to  E=', E/de, ' GeV'
      call cerrorMsg(msg, 1)
      call epwt1dTbl(
     * 'Pair total x-sec. with LPM', erg, tbl, size, media%name)
      end


!     **********************************************
      subroutine epCrePrH(cnst, blb, sizeu,  sizee)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

!

       type(SmpCnst)::  cnst  ! must be media.cnst
      integer sizeu, sizee
      real*8  blb(sizeu, sizee)

      real*8 de1, E, vl,  u, tcp,
     *   eps, v

      external epPrgeneSolv
      real*8  epPrgeneSolv

      integer  iu, j, ie


      common/upsic/upsi, vmax
      real*8 vmax, upsi

      data eps/1.d-5/

      E = cnst%PrEg1H*1.000001d0
      de1=10.**cnst%PrdEH
!
      call cerrorMsg(
     *'Creating pair sampling table in the LPM region',1)
      do   ie = 1, sizee
         Egme = E/masele
         vmax = 1.d0 - masele/E
!         call eptotcpH(.5d0, vmax,  tcp)
         call epPrgeneTX(0.5d0, vmax,  tcp)
!         tcp = tcp*2
         u = cnst%PrU1H
         vl=.5d0
         do   iu = 2, sizeu-1
            u = u + cnst%PrdU1H
            upsi = u*tcp
            call kbchop(epPrgeneSolv, vl, vmax, eps, v, j)
            if(j .le. 0) then
               write(0,*) ' error in epCrePrH, E=',E,' u=',u
               write(0, *) ' tcp=',tcp, ' vl=',vl, ' vmax=',vmax
            endif
            blb(iu,ie)=v
         enddo
         blb(1,  ie)=1.
         blb(sizeu, ie)=0.5
         E = E*de1
       enddo

       call epwt2dTbl('Pair table with LPM', blb, sizeu, sizee)
       call cerrorMsg('the table has been created',1)
      end

