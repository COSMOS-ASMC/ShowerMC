!     ****************************************************************
!     *                                                              *
!     * epCreBrSTbH:   create sampling table for brems
!     *            with landau effect at high energies               *
!     *                                                              *
!     ****************************************************************
      subroutine epCreBrSTbH(mediain, cnst)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"
!     *
!
       type(epmedia):: mediain  !  media
       type(SmpCnst)::  cnst  ! must be media.cnst
!
!
      integer ws
      parameter (ws = 10000)
      real*8 work(ws), work2(ws)
!        cp mediain into common area
      media = mediain
!         print brem const. at high energies
      call epwtBrCnstH(cnst)

!            total x-sec.
      if( cnst%BrneH .gt. mxBrTXH) then
         call cerrorMsg(
     *    'too large LPM Brem total X-section tbl',0)
      endif
      call epCreBrTXTH(cnst, work, work2, cnst%BrneH)
!           sampling table in A region

      if( cnst%Brnu1H*cnst%BrneH2 .gt. mxBrTblHA) then
         call cerrorMsg('brem table for HA is too large', 0)
      endif
      call  epCreBrHA(cnst, work,  cnst%Brnu1H, cnst%BrneH2)


!           sampling table in B region
      if( cnst%Brnu2H*cnst%BrneH2 .gt. mxBrTblHB) then
         call cerrorMsg('brem table for HB is too large', 0)
      endif
      call  epCreBrHB(cnst, work, cnst%Brnu2H, cnst%BrneH2)
      end
!     ****************************************
      subroutine epCreBrTXTH(cnst, erg,  tbl, size)
!     ****************************************
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

       type(SmpCnst)::  cnst  ! must be media.cnst
      integer size
      real*8  tbl(size), erg(size)

      real*8 E,  vmax, vmin,  tcb, de
      character*160 msg
      integer i


      E = cnst%BrEe1H*1.00001d0  !  normally few x 0.1 GeV

      de = 10.**cnst%BrdEH
      do i = 1, size
         Eeme = E/masele
!           need not chk
!         call epChkEndValH(E)

         vmax= (1.d0 - masele/ E)
         vmin=cnst%BrEgminH
         if(vmin .ge.  vmax) then
            call cerrorMsg('Eg/E <= vmin', 0)
         endif
!         call eptotcbH(vmin,  vmax,  tcb)
         call epBrgeneTX(vmin,  vmax,  tcb)
!         tbl(i) = tcb* media.mbtoPX02  ! this is wrong since tcb has been
!                                  normalized already
         tbl(i) = tcb* media%mbtoPX0  ! prob. per radiation length.
         erg(i) = E
         E = E * de
      enddo
      write(msg, *) 
     * 'Matter=',media%name,
     * ' Brems total X-sec. table with LPM has'//
     * ' been made up to  E=', E/de, ' GeV'
      call cerrorMsg(msg, 1)

      write(msg, *) 'Brems total X-sec. table with LPM upto E=',
     *  E/de, ' GeV' 
      call epwt1dTbl(msg,  erg, tbl, size, media%name)
      end
!     **********************************************
      subroutine epCreBrHA(cnst, bla, sizeu,  sizee)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

       type(SmpCnst)::  cnst  ! must be media.cnst
      integer sizeu, sizee
      real*8  bla(sizeu,  sizee)

      real*8 de1, E, vmin, vl, vr, u, tcb, vx, tcbx,
     *   eps, v

      external epBrgeneSolv
      real*8  epBrgeneSolv


      character*160 msg
      integer i, iu, j


      common/upsic/upsi, vmax
      real*8 vmax, upsi

      real*8 sqrtv
      real*8  epBrgenex

      data eps/1.d-8/


      de1 = 10.**cnst%BrdEH2

      E = cnst%BrEe1H*1.00001d0 

      write(msg, *) ' Creating Brems sampling table with LPM: E>=',
     *    cnst%BrEe1H
      call cerrorMsg(msg, 1)
      call cerrorMsg('for small Eg region', 1)
      
      do  i = 1, sizee
         Eeme = E/masele
!           need not chk         
!         call epChkEndValH(E)

         vmin=cnst%BremEgmin
         vmax=1.d0 - masele/E
         sqrtv = sqrt(vmin)
!         call eptotcbH(vmin, vmax, tcb)
         call epBrgeneTX(vmin, vmax, tcb)
         vl=vmin
         vr=vmax

         u=cnst%BrU1H
         do   iu = 1, cnst%Brnu1H-1
            upsi=u*tcb
            call kbchop(epBrgeneSolv, vl, vr, eps, v, j)
            if(j .le. 0) then

               write(0, *) ' cond, E, u=', j, E, u
               do vx = vmin, vmax, (vmax-vmin)/100.d0
!                  call eptotcbH(vx, vmax, tcbx)
                  call epBrgeneTX(vx, vmax, tcbx)
                  write(0, *) vx, tcbx
               enddo

            endif
            bla(iu,i)=( sqrt(v) - sqrtv)/(1.d0 -u)
            u=u+cnst%BrdU1H
!            vr = v
         enddo
         bla(sizeu ,i)=  tcb/sqrtv/epBrgenex(vmin)/2.0d0
         E =E* de1
      enddo

      call epwt2dTbl('Brems with LPM sampling tbl for small Eg',
     *                bla, sizeu, sizee)
      call cerrorMsg('the table has been created',1)
      end
!     ***********************************************
      subroutine epCreBrHB(cnst, blb, sizeu,  sizee)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

       type(SmpCnst)::  cnst  ! must be media.cnst
      integer sizeu, sizee
      real*8  blb(sizeu, sizee)

      real*8 de1, E, vmin, vl, vr, u, tcb,
     *   eps, us, v

      external epBrgeneSolv
      real*8  epBrgeneSolv
      

      integer  iu, j, ie


      common/upsic/upsi, vmax
      real*8 vmax, upsi          
      real*8 aa/4.d0/

      real*8  epBrgenex

      call cerrorMsg(
     * 'Creating LPM brems table for large Eg',1)

      E = cnst%BrEe1H*1.0001d0

      de1 = 10.d0**cnst%BrdEH2
      do ie = 1, sizee
          Eeme = E/masele
!           need not check
         call epChkEndValH(E)

          vmin=cnst%BrEgminH
          vmax=(1.d0 - masele/E)*0.999d0
!          call eptotcbH(vmin, vmax, tcb)
          call epBrgeneTX(vmin, vmax, tcb)
          vl=vmin
          vr=vmax
          us= cnst%BrU3H
          do  iu = 2, sizeu
             us = us + cnst%BrdU2H
             u = us**aa
             upsi=u*tcb
             if(iu .le. 3) then
                eps = 1.d-8
             else
                 eps = 1.d-8
             endif
             call kbchop(epBrgeneSolv, vl, vr, eps, v, j)
!             call  kbinChop(epBrgeneSolv, vl, vr, guess,
!     *            eps, v, j)
             if(j .le. 0) then
                write(0,'("back from bchop: E,u,v=",1p,3g18.7)') 
     *         E, u, v
             endif
             blb(iu,ie) = log(vmax/v)/u
!             vr =  v
          enddo
          blb(1, ie) = tcb/epBrgenex(vmax)/vmax
          E= E*de1
       enddo

       call epwt2dTbl(
     * 'Brems with LPM smp tbl; (sqrt(v)-sqrt(vmin)/(1-u)',
     *  blb, sizeu, sizee)
       call cerrorMsg('the talbe has been created', 1)
       end

      subroutine  epChkEndValH(Ee)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

      real*8 Ee  ! in/out
!            table end value should not exceed the region end value
      if(Ee .gt. media%cnst%BrEe2H) then
         if(Ee .gt. media%cnst%BrEe2H*1.0001d0) then
            write(0,*) '*********** too large end value *** '
            write(0,*) 'Ee=',Ee, 
     *           '  >media%cnst%BrEe2H*1.0001d0=',
     *          media%cnst%BrEe2H*1.0001d0
            stop
         else
            Ee =media%cnst%BrEe2H
         endif
      endif
      end
