      subroutine epCreMuBTab(mediain, cnst)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZmuBPNgene.h"

!     *                                                             
!     *       create sampling table for brems
!     *       of muons 
!
       type(epmedia):: mediain  !  media
       type(SmpCnst)::  cnst  ! must be media.cnst
!
      integer ws, ws1, ws2
      parameter (ws = 3000, ws1=500, ws2=500)
!      real*8 work(ws1, ws2),
       real*8  work1(ws),  work2(ws),  work3(ws), work4(ws)
!        cp mediain into common area
      media = mediain
      call epwtmuBrCnst(cnst)

!            total x-sec.
      if( cnst%muBrTXT .gt. mxMuBrTX) then
         call cerrorMsg(
     *     'too large total X-sec.table for mu Brems ', 0)
      endif
      call epCreMuBrTXT(
     *       cnst, work1, work2,  work3, work4, cnst%muBrTXT)

!        sampling will be performed by rejection method
!     so that we don't make table. 
!           sampling table
!      if( cnst.muBrUsize * cnst.muBrEsize .gt. mxMuBrTbl ) then
!         call cerrorMsg('Too large Mu Brems table requested', 0)
!      endif
!      call  epCreMuBr(cnst, work, cnst.muBrUsize, cnst.muBrEsize)

      end

!     ****************************************
      subroutine epCreMuBrTXT(
     *       cnst, erg,  txs, tdEdx0, tdEdxt, size)
!     ****************************************
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZmuBPNgene.h"


       type(SmpCnst)::  cnst  ! must be media.cnst
      integer size
      real*8  txs(size), erg(size), tdEdx0(size), tdEdxt(size)

      real*8  xmax, xmin,  de, tprob
      character*160 msg
      integer i
!      integer klena

      real*8  dEdx0,  dEdxt

      write(msg,*) ' Matter=', media%name,
     *  ': Computing the total Xsec of Mu Brems. from E= ',
     *    cnst%muBrEmin, ' GeV'
      call cerrorMsg(msg, 1)
      Emu = cnst%muBrEmin
      de = 10.d0**cnst%muBrdETX


      do   i=1, size
         xmin = cnst%muBrVmin
         xmax = 1.d0 - masmu/Emu
         call eptotcmuB(xmin, xmax, tprob)
         call epmuElossB(0.d0, xmin, dEdx0)
!          total loss
         call epmuElossB(0.d0, xmax, dEdxt)
         erg(i) = Emu
         txs(i) =tprob * media%mbtoPX0 ! prob/X0
!            dE/dx(v<vmin)/Emu
         tdEdx0(i) = dEdx0 * media%mbtoPgrm !  /(g/cm^2)
!            dE/dx/(all v)/Emu
         tdEdxt(i) = dEdxt * media%mbtoPgrm !  /(g/cm^2)

         Emu = Emu * de
      enddo
      write(msg, *) 'Table has been made up to  Emu=',
     *       Emu/de,' GeV'
      call cerrorMsg(msg, 1)
      write(msg, *) 
     * 'Muon brems total X-sec. table (Prob./X0) upto E=',
     *     Emu/de,' GeV'
      call epwt1dTbl(msg, erg, txs, size, media%name)
!        at higher energies, const.
!
!      pw1 =( log10(txs(size-1)/txs(size-2))/cnst.muNdETX +
!     *       log10(txs(size)/txs(size-1))/cnst.muNdETX +
!     *       log10(txs(size)/txs(size-2))/cnst.muNdETX/2)/ 3.d0
!      write(msg, *) pw1,
!     * ' = power of energy dependence at higher energies'
!      call cerrorMsg(msg, 1)
!      write(*,*)  msg(1:klena(msg))
      
      msg= 'dE/dx(v<vmin)/Emu (/(g/cm2))by muon brems'
      call epwt1dTbl(msg, erg, tdEdx0, size, media%name)
!      pw2 =( log10(tdEdx0(size-1)/tdEdx0(size-2))/cnst.muNdETX +
!     *       log10(tdEdx0(size)/tdEdx0(size-1))/cnst.muNdETX +
!     *       log10(tdEdx0(size)/tdEdx0(size-2))/cnst.muNdETX/2)/3.d0
!      write(msg, *) pw2,
!     * ' = power of energy dependence at higher energies'
!      call cerrorMsg(msg, 1)
!      write(*,*)  msg(1:klena(msg))

      msg= 'dE/dx(v<vmax)/Emu (/(g/cm2))by muon brems'
      call epwt1dTbl(msg, erg, tdEdxt, size, media%name)
!      pw3 =( log10(tdEdxt(size-1)/tdEdxt(size-2))/cnst.muNdETX +
!     *       log10(tdEdxt(size)/tdEdxt(size-1))/cnst.muNdETX +
!     *       log10(tdEdxt(size)/tdEdxt(size-2))/cnst.muNdETX/2)/3.d0
!      write(msg, *) pw3,
!     * ' = power of energy dependence at higher energies'
!      call cerrorMsg(msg, 1)
!      write(*,*) msg(1:klena(msg))


      end

!     **********************************************
      subroutine epCreMuBr(cnst, tbl, sizeu,  sizee)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZmuBPNgene.h"
!         make 2-D sampling table for muon nuclear interaction
 
       type(SmpCnst)::  cnst  ! must be media.cnst
      integer sizeu, sizee 
      real*8  tbl(sizeu,  sizee)

      real*8 de, vmin,  u, tcn,  tcnx,  v, dv, vmax, y, error
      integer nvmax
      parameter (nvmax=1000)
      real*8 uv(nvmax), va(nvmax)

      real*8 uu

      character*160 msg
      integer i, j, jmax, k, nerr

!      data a/0.02d0/

      de = 10.d0**cnst%muBrdE

      Emu = cnst%muBrEmin    

      write(msg, *)
     *  ' Creating Mu Brems sampling table  E>=',
     *    cnst%muBrEmin
      call cerrorMsg(msg, 1)
      nerr = 0
      do  i = 1, sizee
         vmin= cnst%muBrVmin
         vmax=1.d0 -  masmu/Emu
         call eptotcmuB(vmin, vmax, tcn)
         v = vmin
         uv(1) = 0.
         va(1) = vmin
         dv =  10.d0**0.01d0
         j= 2
         do  while(.true.)
            v = v * dv
            va(j) = v
            if(v .ge. vmax) then
               jmax = j
               goto 10
            endif
            call eptotcmuB(vmin, v, tcnx)
            uv(j) = min(tcnx/tcn, 1.d0)
            if(uv(j-1) .ge. uv(j)) then
               nerr = nerr + 1
!               numerical error. should be uv(j-1) < uv(j) 
!               neglect this j-th value
            else
               if(uv(j) .eq. 1.) then
                  uv(j) = (uv(j) + uv(j-1))/2
               endif
!/////////////
!         write(*,*) sngl(va(j)), sngl(uv(j))
!/////////////
               j = j + 1
            endif
         enddo

 10      continue
         va(jmax) = vmax
         uv(jmax) = 1.0
!
!            u -> v table
!         uniform in u

!         tbl(1, i) = vmin
!         do k = 2, sizeu-1
!            u = u + cnst.muNdU
!            if(u .gt. 0.50001d0) then
!               kmax = k
!               goto 20
!            endif
!            call kpolintpFE(uv, 1,  va,  1, jmax, 4,  u, y, error)
!            tbl(k, i) = y
!         enddo
!
! 20      continue
!        
!            uniform in uu= vm**(1-u)

         uu = 0.
         tbl(1, i) = vmin
         do k = 2, sizeu-1
            uu = uu + cnst%muBrdU
!            u = 1.d0 - log( uu )* vmin 
            u = uu
            call kpolintpFE(uv, 1,  va,  1, jmax, 4,  u, y, error)
            tbl(k, i) = y
         enddo
         tbl(sizeu, i) = vmax
         Emu =Emu * de
      enddo

      call epwt2dTbl(
     *  'muon brems  sampling table ',
     *   tbl, sizeu, sizee)
      call cerrorMsg('the table has been created', 1)
      if(nerr .gt. 0) then
         call  cerrorMsg(
     *   '              ********************************* ', 1)
         call  cerrorMsg(
     *  'warning: at making sampling tab for muon brems',
     *   1)
         write(msg, *) 
     *  'The number of numerical precision errors occurred=',nerr
         call cerrorMsg(msg, 1)
         call cerrorMsg(
     *   'You should check the 2-D table if there are some'//
     *   ' irregurer part such as > 1',1)
         call cerrorMsg(
     *   'If there is, make smoothing of the tab. or'//
     *   ' try a smaller error bound in epmuAuxB.f(epsa, epsr)', 1)
      endif
      end
