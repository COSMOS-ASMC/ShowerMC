      subroutine epCreMuPTab(mediain, cnst)
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
!       real*8 work(ws1, ws2)
       real*8  work1(ws),  work2(ws),  work3(ws), work4(ws)
!        cp mediain into common area
      media = mediain
      call epwtmuPrCnst(cnst)

!            total x-sec.
      if( cnst%muPrTXT .gt. mxMuPrTX) then
         call cerrorMsg(
     *     'too large total X-sec.table for mu pair creation', 0)
      endif
      call epCreMuPrTXT(
     *       cnst, work1, work2,  work3, work4, cnst%muPrTXT)

!        sampling can be performed by rejection method
!        and we don't make table here
!           sampling table
!      if( cnst.muPrUsize * cnst.muPrEsize .gt. mxMuPrTbl ) then
!         call cerrorMsg('Too large Mu pair creation tab requested', 0)
!      endif
!      call  epCreMuPr(cnst, work, cnst.muPrUsize, cnst.muPrEsize)

      end

!     ****************************************
      subroutine epCreMuPrTXT(
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

      real*8  xmax, xmin,  de, tprob, epmuPairVmn, vcut
      real*8  epmuvmax2
      character*160 msg
      integer i
!      integer klena

      real*8  dEdx0,  dEdxt
      real*8 pw1, pw2, pw3

      write(msg,*) ' Matter=', media%name,
     *  ': Computing the total Xsec of Mu pair cre. from E= ',
     *    cnst%muPrEmin, ' GeV'
      call cerrorMsg(msg, 1)
      Emu = cnst%muPrEmin
      de = 10.d0**cnst%muPrdETX

      do   i=1, size
         xmin = cnst%muPrVmin  ! > 2masele/Emu
         vcut = epmuPairVmn(Emu)

         xmin = max(xmin, vcut)
!         xmax = 1.d0 - masmu/Emu
         xmax =  epmuvmax2(media, Emu)
         call eptotcmuP(xmin, xmax, tprob)  ! prob in mb

         call epmuElossP(vcut, xmin, dEdx0)  ! dE/dx/Emu *mb  

!          total loss
         call epmuElossP(vcut, xmax, dEdxt)  ! dE/dx/Emu *mb
         erg(i) = Emu
         txs(i) =tprob * media%mbtoPX0 ! prob/X0
!            dE/dx(v<vmin)/Emu
         tdEdx0(i) = dEdx0 * media%mbtoPgrm ! dE/dx/Emu /(g/cm^2)
!            dE/dx/(all v)/Emu
         tdEdxt(i) = dEdxt * media%mbtoPgrm ! dE/dx/Emu (g/cm^2)

         Emu = Emu * de
      enddo
      write(msg, *) 'Table has been made up to  Emu=',
     *       Emu/de,' GeV'
      call cerrorMsg(msg, 1)
      write(msg, *) 
     * 'Muon pair cre. total X-sec. table (Prob./X0) upto E=',
     *     Emu/de,' GeV'
      call epwt1dTbl(msg, erg, txs, size, media%name)
!
      pw1 =( log10(txs(size-1)/txs(size-2))/cnst%muPrdETX +
     *       log10(txs(size)/txs(size-1))/cnst%muPrdETX +
     *       log10(txs(size)/txs(size-2))/cnst%muPrdETX/2)/ 3.d0
      write(msg, *) pw1,
     * ' = power of energy dependence at higher energies'
      call cerrorMsg(msg, 1)
!      write(*,*)  msg(1:klena(msg))
      
      msg= 'dE/dx(v<vmin)/Emu (/(g/cm2))by muon pair cre.'
      call epwt1dTbl(msg, erg, tdEdx0, size, media%name)
      pw2 =( log10(tdEdx0(size-1)/tdEdx0(size-2))/cnst%muPrdETX +
     *       log10(tdEdx0(size)/tdEdx0(size-1))/cnst%muPrdETX +
     *       log10(tdEdx0(size)/tdEdx0(size-2))/cnst%muPrdETX/2)/3.d0
      write(msg, *) pw2,
     * ' = power of energy dependence at higher energies'
      call cerrorMsg(msg, 1)
!      write(*,*)  msg(1:klena(msg))

      msg= 'dE/dx(v<vmax)/Emu (/(g/cm2))by muon pair cre.'
      call epwt1dTbl(msg, erg, tdEdxt, size, media%name)
      pw3 =( log10(tdEdxt(size-1)/tdEdxt(size-2))/cnst%muPrdETX +
     *       log10(tdEdxt(size)/tdEdxt(size-1))/cnst%muPrdETX +
     *       log10(tdEdxt(size)/tdEdxt(size-2))/cnst%muPrdETX/2)/3.d0
      write(msg, *) pw3,
     * ' = power of energy dependence at higher energies'
      call cerrorMsg(msg, 1)
!      write(*,*) msg(1:klena(msg))


      end

!     **********************************************
      subroutine epCreMuPr(cnst, tbl, sizeu,  sizee)
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



      de = 10.d0**cnst%muPrdE

      Emu = cnst%muPrEmin    

      write(msg, *)
     *  ' Creating Mu Pair cre. sampling table  E>=',
     *    cnst%muPrEmin
      call cerrorMsg(msg, 1)
      nerr = 0
      do  i = 1, sizee
         vmin= cnst%muPrVmin
         vmax=1.d0 -  masmu/Emu
         call eptotcmuP(vmin, vmax, tcn)
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
            call eptotcmuP(vmin, v, tcnx)
            uv(j) = min(tcnx/tcn, 1.d0)
            if(uv(j-1) .ge. uv(j)) then
               nerr = nerr + 1
!               numerical error. should be uv(j-1) < uv(j) 
!               neglect this j-th value
            else
               if(uv(j) .eq. 1.) then
                  uv(j) = (uv(j) + uv(j-1))/2
               endif
!//////////
!               write(*,*) sngl(uv(j)), sngl(va(j))
!//////////
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
            uu = uu + cnst%muPrdU
!            u = 1.d0 - log( uu )* vmin 
            u = uu
            call kpolintpFE(uv, 1,  va,  1, jmax, 4,  u, y, error)
            tbl(k, i) = y
         enddo
         tbl(sizeu, i) = vmax
         Emu =Emu * de
      enddo

      call epwt2dTbl(
     *  'muon cre.  sampling table ',
     *   tbl, sizeu, sizee)
      call cerrorMsg('the table has been created', 1)
      if(nerr .gt. 0) then
         call  cerrorMsg(
     *   '              ********************************* ', 1)
         call  cerrorMsg(
     *  'warning: at making sampling tab for muon pair cre.',
     *   1)
         write(msg, *) 
     *  'The number of numerical precision errors occurred=',nerr
         call cerrorMsg(msg, 1)
         call cerrorMsg(
     *   'You should check the 2-D table if there are some'//
     *   ' irregurer part such as > 1',1)
         call cerrorMsg(
     *   'If there is, make smoothing of the tab. or'//
     *   ' try a smaller error bound in epmuAuxP.f(epsa, epsr)', 1)
      endif
      end
