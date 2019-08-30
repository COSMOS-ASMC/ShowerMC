#include "ZepicsBD.h"
      implicit none
!
!         to write data needed to draw brems function at arbitrary energies
!   This is the same as DrawSelzterBrems.f except for
!   this is only for Selzter brems function and disregards
!   media.cnst.BrEemaxS and make it 10 GeV.
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"
#include "ZepTrackp.h"

      character*50  file
      character*16  basemedia
      integer nerg
      integer io, i, norm
      real*8  Ee, x, f, epBrgenex, xmax
      real*8  xmin, tprob, xx, fxx, Nc
      real*8  E1

      io = 10
      write(0,*)
     *  'echo norm media Ek LPMeffect | drawbrem*.out'
      call cerrorMsg('1) In unit of / r%l (1)',  1)
      call cerrorMsg('2) In unit of mb/ingredient', 1) 
      call cerrorMsg('3) In unit of /(g/cm^2)', 1)
      call cerrorMsg('4) In unit of /cm', 1)
      call cerrorMsg(
     * '5) Area-normalization '//
     * ' for comparison with M.C data such as k vs xdN/dx/Nt or '//
     * ' k vs dN/dx/Nt',   1)

      read(*,*)  norm, basemedia, E1


      file ="../../../Data/BaseM/"//basemedia

      open(io, file=file, action ='read')
      call epBPgeneini(io, media)
      media%cnst%BrEemaxS = 10.01  ! force

      Ee = E1 + masele

      Eeme = Ee/masele
      if(Ee .lt. media%cnst%BrEemaxS) then
!            xmin = media.cnst.BrEgminS/Ee
         xmin = media%cnst%BrEgminS
      else
!            xmin = media.cnst.BremEgmin/Ee
         xmin = media%cnst%BremEgmin
         write(0,*) 'should not come'
         stop
      endif
      xmax = 1. - masele/Ee
      write(0,*) ' xmin =', xmin
      write(0,*) ' xmax =', xmax
      call epBrgeneTX(xmin, xmax, tprob)
      if(norm .eq. 5) then
         Nc = 1./tprob
      elseif(norm .eq. 1) then 
         Nc=media%mbtoPX0  ! prob/X0
      elseif( norm .eq.  2) then
         Nc = 1
      elseif(norm .eq.  3) then
         Nc = media%mbtoPgrm
      elseif(norm .eq. 4) then
         Nc = media%mbtoPcm
      else
         call cerrorMsg('input error for norm',0)
      endif
      write(0,*) Ee- masele, tprob*NC
      x = xmin
             
      do while (x .lt.  xmax)
         f = epBrgenex(x) *NC
         if(f .gt. 0.) then
!           f = ds/dx so fxx= ds/dxx is (Ek/Ee)* ds/dx
            xx = x*Ee/E1  ! xx = Eg/Ek
            fxx = f*(E1/Ee)
            write(*, '(1p, 5g12.4)')
     *           xx, xx*fxx, Ee, Nc, tprob
         else
            write(0,*) 'f =', f, '<0 at ',x, Ee,
     *      ' in DrawBrams'
             stop 1235 
         endif
         if(x .gt.  xmax-0.5d0 ) then
            x = x + 1.d-3
         else
            x = x * 10.**0.01
         endif
      enddo
      write(*,*)
      end
