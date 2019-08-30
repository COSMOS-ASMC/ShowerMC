!     **********************
      subroutine epCrePrSTbl1(mediain, cnst)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

!     *                                                             
!     *            create sampling table for 
!     *            Pair creation  at low energies where
*     *            LPM effect can be neglected
!     *            all other correction included                    
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

      call epwtPrCnst(cnst)
!

!           total cross-section
      if(cnst%PairTXTL .gt. mxPrTXL) then
         call cerrorMsg(
     *     'Low energy pair total X-sec%tbl too large',0) 
      endif
      call epCrePrTXTL(cnst, work, work2, cnst%PairTXTL)

!           sampling table in A region
      if( cnst%PairUszLA * cnst%PairEsize .gt. mxPrTblLA) then
         call cerrorMsg('pair table for LA is too large', 0)
      endif
      call  epCrePrLA(cnst, work, cnst%PairUszLA, cnst%PairEsize)

!           sampling table in B region
      if( cnst%PairUszLB * cnst%PairEsize .gt. mxPrTblLB) then
         call cerrorMsg('pair table for LB is too large', 0)
      endif
      call  epCrePrLB(cnst, work, cnst%PairUszLB, cnst%PairEsize)
      end
!     ****************************************
      subroutine epCrePrTXTL(cnst,  erg,  tbl, size)
!     ****************************************
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

       type(SmpCnst)::  cnst  ! must be media.cnst
      integer size
      real*8  tbl(size), erg(size)

      real*8 E,  vmax, vmin,  tcp, de1
      character*160 msg
      integer i


      write(msg,*) ' Matter=', media%name,
     *  ': Computing the total Xsec of Pair from Eg= ',
     *    cnst%PairEgmin, ' GeV'
      call cerrorMsg(msg, 1)


      de1=10.d0**cnst%PairdETXL

      E = cnst%PairEgmin*1.0000001d0     ! Eg. 1.000..1d0 is for IFC

      do   i = 1, size

         call epPrChkEndValL(E,i, size)

         Egme = E/masele
         vmin=.5
         vmax= 1.d0 - masele/E
         if(vmin .ge. vmax) then
            call cerrorMsg('E/Eg <= vmin', 0)
         endif
!         call eptotcp(vmin, vmax, tcp)
         call epPrgeneTX(vmin, vmax, tcp)
         tbl(i)=tcp*2.*media%mbtoPX0 ! prob. per rad. length
         erg(i) = E
         E=E*de1
       enddo
       
       write(msg, *) 'Table has been made up to  E=',
     *    E/de1,' GeV'
       call cerrorMsg(msg, 1)

       call  epwt1dTbl('total pair x-sec. at low E',
     *  erg,  tbl, size, media%name)
       end

!     **********************************************
      subroutine epCrePrLA(cnst, pla, sizeu,  sizee)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

       type(SmpCnst)::  cnst  ! must be media.cnst
      integer sizeu, sizee
      real*8  pla(sizeu,  sizee)

      real*8 de1, E, vmin, vl, vr, u, tcp, 
     * eps, v

      external epPrgeneSolv
      real*8  epPrgeneSolv

      real*8  epPrgenex

      character*160 msg
      integer ie, iu, j


      common/upsic/upsi, vmax
      real*8 vmax, upsi

      data eps/1.d-5/
!
!
!          region a
!
      call cerrorMsg(
     * 'Creating Pair sampling table at low Eg for small Ee',1)
      vmin=.5
      de1=10.0**cnst%PairdELA

      E = cnst%PairEgmin * 1.0000001d0    ! ifc needs 1.00...d0

      do   ie = 1, sizee
          vmax=1.- masele/E
          call epPrChkEndValL(E,ie, sizee)
          Egme = E/masele
!          call eptotcp(vmin,vmax,tcp)
          call epPrgeneTX(vmin,vmax,tcp)
          vl=vmin
          vr=vmax
          u=cnst%PairUminLA
          do   iu = 1, sizeu - 1
             upsi=u*tcp
             call kbchop(epPrgeneSolv, vl, vr, eps, v, j)
             if(j .le. 0) then
                  write(msg,*)' Err in epCrePrLA  E,u=', E, u
                  call cerrorMsg(msg, 0)
             endif
             pla(iu,ie) = (v-.5d0)/(1.d0-u)
             u = u + cnst%PairdULA
          enddo
          pla(sizeu,ie)=tcp/epPrgenex(0.5d0)
          E = E*de1
       enddo

       call  epwt2dTbl('Pair table (v-.5)/(1-u) at low E',
     *                       pla, sizeu,sizee)
       call cerrorMsg('the table has been created',1)
       end

!     ***********************************************
      subroutine epCrePrLB(cnst, plb, sizeu,  sizee)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

       type(SmpCnst)::  cnst  ! must be media.cnst
      integer sizeu, sizee
      real*8  plb(sizeu, sizee)

      real*8 de1, E, vmin, vl, vr, u, tcp,
     *   eps, us, v, ex

      external epPrgeneSolv
      real*8  epPrgeneSolv


      integer ie, iu, j

      character*160 msg  



      common/upsic/upsi, vmax
      real*8 vmax, upsi

      data eps/1.d-5/

      call cerrorMsg(
     * 'Creating pair sampling table at low Eg for large Ee',1)
      vmin=.5d0
      de1=cnst%PairdELB
      ex=0.
      do   ie = 1, sizee

         E = 10.d0**(ex**2) * cnst%PairEgmin *1.00000001d0  ! ifc needs 1.00.01d0
         call epPrChkEndValL(E, ie, sizee)
         vmax = 1.d0 - masele/E 
         Egme = E/masele
!         call eptotcp(vmin,vmax,tcp)
         call epPrgeneTX(vmin,vmax,tcp)
         vl = vmin
         vr = vmax
         us = cnst%PairUminLB + cnst%PairdULB
         do   iu=2, sizeu
            u= us**4
            upsi=u*tcp
            call kbchop(epPrgeneSolv, vl, vr, eps, v, j)
            if(j .le. 0) then
               write(msg, *)'Err in epCrePrLB,  E, u=', E, u
               call cerrorMsg(msg, 0)
            endif
            plb(iu, ie) = v
            us = us + cnst%PairdULB
         enddo
         plb( 1,  ie)=vmax
         ex=ex+de1
      enddo
      
      call  epwt2dTbl('Pair table at low E; v',
     *             plb, sizeu, sizee)
      call cerrorMsg('the table has been created',1)
      end

      subroutine  epPrChkEndValL(Eg, i, size)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

      real*8 Eg  ! in/out
      integer i, size

!            table end value should not exceed the region end value
      if(i .eq. size) then 
         if( Eg .ge. media%cnst%PrScrE ) then
            if(Eg .gt. media%cnst%PrScrE*1.0001d0) then
               write(0,*) '*********** too large end value *** '
               write(0,*) 'Eg=',Eg, 
     *           '  > media%cnst%PrScrE*1.0001d0=',
     *           media%cnst%PrScrE*1.0001d0
               stop
            endif
            Eg = media%cnst%PrScrE*0.999999999d0
         endif
      endif
      end






