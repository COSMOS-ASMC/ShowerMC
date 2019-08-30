#include "ZepicsBD.h"
      implicit none
!
!         to write data needed to draw pair function at arbitrary energies
!
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
!   #include "ZBPgene.h"
#include "ZepTrackp.h"

      character*120  file
      character*8  basemedia
      integer io, i, norm, nerg
      real*8  Eg, x, f, xmax, xmin, xx, fxx, Egme
      real*8  tprob, epPrgene, dx, Nc
       type(epmedia)::  media
      character(8):: force
      integer how


      io = 10
!             next will result  in error for intel compiler
!          with some option
!      call cerrorMsg(' ', 1)
!       write(0,*)
!     * 'echo norm media Eg LPMeffect | drawpair2.out'
!      call cerrorMsg('1) In unit of / r.l (1)',  1)
!      call cerrorMsg('2) In unit of mb/ingredient', 1)
!      call cerrorMsg('3) In unit of /(g/cm^2)', 1)
!      call cerrorMsg('4) In unit of /cm', 1)
!      call cerrorMsg(
!     * '5) Area normalization for comarison with M.C'//
!     * ' result such as x  vs ds/dx ',    1)


      read(*,*)  norm, basemedia, Eg, LPMeffect

      write(0,*) 'norm, basemedia, Eg, LPMeffect ='
      write(0,*) norm, basemedia, Eg, LPMeffect 


      file ="$EPICSTOP/Data/BaseM/"//trim(basemedia)

!     open(io, file=file, action ='read')
      how =  0  ! dummy for pair

      if(Eg .le. 2*masele) stop
      Egme = Eg/masele

      call epBPgeneini(file,  media, how)
      call epPrgenePreInte(media,  Egme)

      xmax = 1. - masele/Eg
!       if only write x>0.5
!     xmin = 0.5d0
!       Ee/Eg min = me/Eg
      xmin = masele/Eg
      x = xmin
      dx =min( (xmax - xmin)/100., 1.d-3)
      call epPrgeneTX(xmin, xmax, tprob)
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

      
      do while (x .lt.  xmax)
         f = epPrgene(media, Eg, x)*Nc
!              f=ds/dx x=Ee/Eg
!              xx= (Ee-masele)/(Eg-2*m)=(x-m/Eg)/(1-2m/Eg)
!         hence dxx = dx/(1-2m/Eg)
!             fxx= ds/dxx = (1-2m/Eg)* (ds/dx) = (1-2m/Eg)*f
         xx = (x*Eg-masele)/(Eg-2*masele)
         fxx = (1.-2.0*masele/Eg)*f
         write(*,'(1p,5g12.4)')
     *        xx, fxx, Eg, Nc, tprob
         x = x + dx
      enddo
      write(*,*)
      end
