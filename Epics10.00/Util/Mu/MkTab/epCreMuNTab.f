      subroutine epCreMuNTab(mediain, cnst)
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
      real*8 work(ws1, ws2),
     * work1(ws),  work2(ws),  work3(ws), work4(ws)
!        cp mediain into common area
      media = mediain
      call epwtmuNCnst(cnst)

!            total x-sec.
      if( cnst%muNTXT .gt. mxMuNTX) then
         call cerrorMsg(
     *     'too large total X-sec.table for mu nuc int', 0)
      endif
      call epCreMuNTXT(
     *       cnst, work1, work2,  work3, work4, cnst%muNTXT)

!           sampling table
      if( cnst%muNUsize * cnst%muNEsize .gt. mxMuNTbl ) then
         call cerrorMsg('Too large Mu N-int table requested', 0)
      endif
      call  epCreMuN(cnst, work, cnst%muNUsize, cnst%muNEsize)



      end

!     ****************************************
      subroutine epCreMuNTXT(
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
      integer i, klena

      real*8  dEdx0,  dEdxt, pw1, pw2, pw3

      write(msg,*) ' Matter=', media%name,
     *  ': Computing the total Xsec of Mu Nuclear Int. from E= ',
     *    cnst%muNEmin, ' GeV'
      call cerrorMsg(msg, 1)
      Emu = cnst%muNEmin
      de = 10.d0**cnst%muNdETX



      do   i=1, size
         xmin = max ( cnst%muNVmin, masrho/Emu)
         xmax = 1.d0 - masmu/Emu
         call eptotcmuN(xmin, xmax, tprob)
         if(xmin .le. masrho/Emu) then
            dEdx0 = 0.
         else
            call epmuElossN(masrho/Emu, xmin, dEdx0)
         endif
!          total loss
         call epmuElossN(masrho/Emu, xmax, dEdxt)
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
     * 'Muon Nuc. Int. total X-sec. table (Prob./X0) upto E=',
     *     Emu/de,' GeV'
      call epwt1dTbl(msg, erg, txs, size, media%name)
!
      pw1 =( log10(txs(size-1)/txs(size-2))/cnst%muNdETX +
     *       log10(txs(size)/txs(size-1))/cnst%muNdETX +
     *       log10(txs(size)/txs(size-2))/cnst%muNdETX/2)/ 3.d0
      write(msg, *) pw1,
     * ' = power of energy dependence at higher energies'
      call cerrorMsg(msg, 1)
      write(*,*)  msg(1:klena(msg))
      
      msg= 'dE/dx(v<vmin)/Emu (/(g/cm2))by muon Nuc. int'
      call epwt1dTbl(msg, erg, tdEdx0, size, media%name)
      pw2 =( log10(tdEdx0(size-1)/tdEdx0(size-2))/cnst%muNdETX +
     *       log10(tdEdx0(size)/tdEdx0(size-1))/cnst%muNdETX +
     *       log10(tdEdx0(size)/tdEdx0(size-2))/cnst%muNdETX/2)/3.d0
      write(msg, *) pw2,
     * ' = power of energy dependence at higher energies'
      call cerrorMsg(msg, 1)
      write(*,*)  msg(1:klena(msg))

      msg= 'dE/dx(v<vmax)/Emu (/(g/cm2))by muon Nuc. int'
      call epwt1dTbl(msg, erg, tdEdxt, size, media%name)
      pw3 =( log10(tdEdxt(size-1)/tdEdxt(size-2))/cnst%muNdETX +
     *       log10(tdEdxt(size)/tdEdxt(size-1))/cnst%muNdETX +
     *       log10(tdEdxt(size)/tdEdxt(size-2))/cnst%muNdETX/2)/3.d0
      write(msg, *) pw3,
     * ' = power of energy dependence at higher energies'
      call cerrorMsg(msg, 1)
      write(*,*) msg(1:klena(msg))


      end

!     **********************************************
      subroutine epCreMuN(cnst, tbl, sizeu,  sizee)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZmuBPNgene.h"
!         make 2-D sampling table for muon nuclear interaction
 
       type(SmpCnst)::  cnst  ! must be media.cnst
      integer sizeu, sizee 
      real*8  tbl(sizeu,  sizee)

      real*8 de, vmin,  u,  tcnx, dv, vmax, y, error
      real*8  v1, v2, temp
      integer nvmax
      parameter (nvmax=1000)
      real*8 uv(nvmax), va(nvmax)

      real*8 uu, a

      character*160 msg
      integer i, j, jmax, k

      data a/0.02d0/  ! if change, must change also epmuNsmp.f

      de = 10.d0**cnst%muNdE

      Emu = cnst%muNEmin    

      write(msg, *)
     *  ' Creating Mu Nuc. Int. sampling table  E>=',
     *    cnst%muNEmin
      call cerrorMsg(msg, 1)
      do  i = 1, sizee
         vmin=max( cnst%muNVmin, masrho/Emu)
         vmax=1.d0 -  masmu/Emu
!         call eptotcmuN(vmin, vmax, tcn)
         v2 = vmin
         uv(1) = 0.
         va(1) = vmin
         dv =  10.d0**0.01d0
         j= 2
         tcnx = 0.
         do  while( v2 .ne. vmax )
            v1 = v2
            v2 = min(v2 * dv, vmax)
            va(j) = v2 
            call eptotcmuN(v1, v2, temp)
            tcnx = tcnx + temp
            uv(j) =  tcnx
            j = j + 1
         enddo
         jmax = j - 1 
         do  j = 1, jmax
            uv(j) = uv(j) / tcnx
         enddo
!        
!           uniform in 0.01*u/(1.01-u**0.5);  good  but inverse cannot be
!                                             obtained.
!                      0.001/(1.001-u**0.05)  good but non 0 at u=0
!                      0.02u/(1.02- u)       good and solvable
!
!                       ( a=0.02 )
!c         uu = a*u/(1+a - u)
         uu = 0.
         tbl(1, i) = log10( vmin )
         do k = 2, sizeu-1
            uu = uu + cnst%muNdU
            u = uu*(1+a)/(uu+a)
            call kpolintpFE(uv, 1,  va,  1, jmax, 6,  u, y, error)
            y = log10(y)
            tbl(k, i) = y
         enddo
         tbl(sizeu, i) =  log10( vmax )
         Emu =Emu * de
      enddo

      call epwt2dTbl(
     *  'muon nuc. int.  sampling table ',
     *   tbl, sizeu, sizee)
      call cerrorMsg('the table has been created', 1)
      end

