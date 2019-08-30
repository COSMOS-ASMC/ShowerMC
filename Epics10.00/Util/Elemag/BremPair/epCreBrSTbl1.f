      subroutine epCreBrSTbl1(mediain, cnst)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

!     *                                                             
!     *       create sampling table for Bremstrahlung 
!     *       at low energies where
!     *       LPM effect can be neglected
!     *       all other correction included                    
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
      call epwtBrCnst(cnst)

!            total x-sec.
      if( cnst%BremTXTL .gt. mxBrTXL) then
         call cerrorMsg(
     *     'too large total X-sec%table for low E bresm', 0)
      endif
      call epCreBrTXTL(cnst, work, work2,  cnst%BremTXTL)

!           sampling table in A region
      if( cnst%BremUszLA * cnst%BremEsize .gt. mxBrTblLA) then
         call cerrorMsg('brem table for LA is too large', 0)
      endif
      call  epCreBrLA(cnst, work, cnst%BremUszLA, cnst%BremEsize)

!           sampling table in B region
      if( cnst%BremUszLB * cnst%BremEsize .gt. mxBrTblLB) then
         call cerrorMsg('brem table for LB is too large', 0)
      endif
      call  epCreBrLB(cnst, work, cnst%BremUszLB, cnst%BremEsize)
      end

!     ****************************************
      subroutine epCreBrTXTL(cnst, erg,  tbl, size)
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



      write(msg,*) ' Matter=', media%name,
     *  ': Computing the total Xsec of Brems from Ee= ',
     *    cnst%BremEemin, ' GeV'
      call cerrorMsg(msg, 1)
      write(0,*) ' log step =',cnst%BremdETXL, ' size=',size
      write(0,*)

      E = cnst%BremEemin*1.000001d0
      de = 10.d0**cnst%BremdETXL
      do   i=1, size

         call epChkEndValSc(E, i, size)

         Eeme = E/masele
         vmax=1.d0 - 1.0001d0 * masele/ E
         vmin=cnst%BremEgmin
         if(vmin .ge.  vmax) then
            call cerrorMsg('Eg/E <= vmin', 0)
         endif
!         call eptotcb(vmin,  vmax,  tcb)
         call epBrgeneTX(vmin,  vmax,  tcb)
         tbl(i) = tcb* media%mbtoPX0  ! prob. per radiation length.
         erg(i) = E
         E = E * de
      enddo
      write(msg, *) 'Table has been made up to  E=', E/de,' GeV'
      call cerrorMsg(msg, 1)
      write(msg, *) 'Brems total X-sec. table upto E=',
     *     E/de,' GeV'
      call epwt1dTbl(msg, erg, tbl, size, media%name)

      end

!     **********************************************
      subroutine epCreBrLA(cnst, bla, sizeu,  sizee)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

       type(SmpCnst)::  cnst  ! must be media.cnst
      integer sizeu, sizee 
      real*8  bla(sizeu,  sizee)

      real*8 de1, E, vmin, vl, vr, u, tcb, vx, tcbx,
     *  tcbleft, tcbright, eps, v

      external epBrgeneSolv
      real*8  epBrgeneSolv

      real*8  epBrgenex

      character*160 msg
      integer i, iu, j


      common/upsic/upsi, vmax
      real*8 vmax, upsi

      data eps/1.d-8/


      de1 = 10.d0**cnst%BremdEL

      E = cnst%BremEemin*1.000001d0

      write(msg, *) ' Creating Brem sampling table P%S region: E>=',
     *    cnst%BremEemin
      call cerrorMsg(msg, 1)
      write(msg, *) ' for small Eg region; log step= ',cnst%BremdEL
      call cerrorMsg(msg, 1)
      write(0, *) ' size=', sizee



      do  i = 1, sizee
         call epChkEndValSc(E, i, sizee)
         Eeme = E/masele
         vmin=cnst%BremEgmin
         vmax=1.d0 -  masele/E
         
!         call eptotcb(vmin, vmax, tcb)
         call epBrgeneTX(vmin, vmax, tcb)
         vl=vmin
         vr=vmax
         
         u=cnst%BremUminLA
         
         do  iu=1, sizeu -1
            
            upsi=u*tcb
            call kbchop(epBrgeneSolv, vl, vr, eps, v, j)
            
            if(j .le. 0) then
               
               write(0, *) ' a-- cond, E, u=', j, E, u
               do vx = vmin, vmax, (vmax-vmin)/100.d0
!                  call eptotcb(vx, vmax, tcbx)
                  call epBrgeneTX(vx, vmax, tcbx)
                  
                  write(0, *) vx, tcbx
               enddo
               write(0,*) 'u upsi, vl, vr', u,
     *              upsi, vl, vr,
     *              ' vmin=', vmin, ' tcb=', tcb
!               call eptotcb(vl,vmax,tcbleft)
!               call eptotcb(vr,vmax,tcbright)
               call epBrgeneTX(vl,vmax,tcbleft)
               call epBrgeneTX(vr,vmax,tcbright)
               
               write(0, *) ' tcb at vl, vr=', tcbleft,
     *              tcbright
            endif
            bla(iu,i)= log(v/vmin)/(1.-u)
#if defined (DEBUG1)
            write(0,'(a, i4, 1p, 4g14.5)')
     *           'tt1 ', iu, E-masele, u, v, bla(iu, i)
#endif            
            u=u + cnst%BremdULA
!            vr = v
         enddo
!         bla(sizeu, i)=tcb/vmin/epBremS(vmin)
         bla(sizeu, i)=tcb/vmin/epBrgenex(vmin)
#if defined (DEBUG1)
            write(0,'(a, i4, 1p, 4g14.5)')
     *           'tt1 ', sizeu, E-masele, u, vmin, bla(sizeu, i)
#endif            
         E =E* de1
      enddo

      call epwt2dTbl(
     *  'Brems sampling table at low E for Small Eg',
     *   bla, sizeu, sizee)
      call cerrorMsg('the table has been created', 1)

      end

!     ***********************************************
      subroutine epCreBrLB(cnst, blb, sizeu,  sizee)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

       type(SmpCnst)::  cnst  ! must be media.cnst
      integer sizeu, sizee 
      real*8  blb(sizeu, sizee)

      real*8 de1, E, vmin, vl, vr, u, tcb, 
     *   eps, us, v, temp

      external epBrgeneSolv
      real*8  epBrgeneSolv
      real*8  epBrgenex

      integer  iu, j, ie


      common/upsic/upsi, vmax
      real*8 vmax, upsi

!      data eps/1.d-5/


      call cerrorMsg(
     * 'Creating brems sampling table P%S for large Eg', 1)

      E = cnst%BremEemin*1.000001d0

      de1 = 10.d0**cnst%BremdEL
!      do ie = 1, cnst.BremEsize
      do ie = 1, sizee
          call epChkEndValSc(E, ie, sizee)
          Eeme = E/masele


          vmin=cnst%BremEgmin

          vmax=1.d0 -  masele/E
!          call eptotcb(vmin, vmax, tcb)
          call epBrgeneTX(vmin, vmax, tcb)
          vl=vmin
          vr=vmax
          us=cnst%BremUminLB + cnst%BremdULB
!          do   iu = 2,  cnst.BremUszLB
          do   iu = 2,  sizeu
             u=us**4
             upsi=u*tcb
             if(iu .le .5) then
                eps = 1.d-8
             else
                eps = 1.d-8
             endif

             call kbchop(epBrgeneSolv, vl, vr, eps, v, j)
             if(j .le. 0) then
                write(0,'(''b-- E,u='',2g12.3)') E, u
             endif
             blb(iu,ie) = log(vmax/v)/u
#if defined (DEBUG1)
            write(0,'(a, i4, 1p, 5g14.5)')
     *           'tt2 ', iu,  E-masele, us, v, bla(sizeu, ie), u
#endif            
             us=us + cnst%BremdULB
!             vr = v
          enddo
!          temp = epBremS(vmax)
          temp = epBrgenex(vmax)
          if(temp .ne. 0.) then
             blb(1, ie) = tcb/vmax/temp
          else
             blb(1, ie) = blb(2, ie)
          endif
#if defined (DEBUG1)
            write(0,'(a, i4, 1p, 5g14.5)')
     *         'tt2 ', 1,  E-masele, cnst%BremUminLB,
     *         vmax, blb(1, ie), cnst%BremUminLB**4
#endif            

          E= E*de1
       enddo
       
       call epwt2dTbl('Brems sampling tbl at low E for large Eg',
     *  blb, sizeu, sizee)
       call cerrorMsg('the table has been created', 1)
       end
      subroutine  epChkEndValSc(Ee, ie, sizee)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

      real*8 Ee  ! in/out

      integer ie, sizee
!            table end value should not exceed the region end value
      if(ie .eq. sizee) then
         if(Ee .ge. media%cnst%BrScrE) then
            if(Ee .gt. media%cnst%BrScrE*1.0001d0) then
               write(0,*) '*********** too large end value *** '
               write(0,*) 'Ee=',Ee, 
     *           '  > media%cnst%BrScrE*1.001d0=',
     *           media%cnst%BrScrE*1.0001d0
               stop
            endif
            Ee = media%cnst%BrScrE*0.99999999d0
         endif
      endif
      end
