      implicit none
!
!         to write data needed to draw muon pair creation function 
!      ds/dr (= ds/dvdr integrated over v)
!
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZmuBPNgene.h"

      character*50  file
      character*8  basemedia

      integer io, i, norm, ie, icon
      real*8  f, xmin, xmax, x, E1
      real*8  tprob,  epmuPrS, Nc
      integer kgetenv2, leng


      io = 10
      call cerrorMsg(' ', 1)
      call cerrorMsg('1) In unit of / r.l (1)',  1)
      call cerrorMsg('2) In unit of mb/ingredient', 1) 
      call cerrorMsg('3) In unit of /(g/cm^2)', 1)
      call cerrorMsg('4) In unit of /cm', 1)
      call cerrorMsg(
     * '5) Area-normalization '//
     * ' for comparison with M.C data such as v vs vdN/dv/Nt  ',
     * 1 )

      
      norm = 1
      write(0,*)'Enter norm, basemedia, E1'
      read(*,*) norm, basemedia, E1

      leng =  kgetenv2("EPICSTOP", file)
      if( leng  == 0 )  then
         write(0,*) "EPICSTOP is not given yet"
         stop
      endif
      file =file(1:leng)//"/Data/Media/"//basemedia
      call copenf(io, file, icon)
      if(icon .ne. 0) then
         call cerrorMsg('media file path invalid',  0)
      endif
      Emu = E1
!       actual xmax may be smaller than 1. and if so,
!       cross section is automatiaally made to be 0.
      call epReadTab(io, media)
      xmin = max( 4*masele/Emu, media%cnst%muPrVmin)
!///////////
!      xmin = 4*masele/Emu  ! this is to drawo entire region
!/////////////
      xmax = 1.  

      call eptotcmuP(xmin, xmax, tprob)   ! in mb
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
      x = xmin 
      do while(x < xmax)
         f = epmuPrS(x) * Nc    ! mb to prob/X0, etc  xds/dx
         if(f .gt. 0.) then
            write(*, '(1p,6g13.4)')
     *         x, f, Emu, Nc, tprob, Nc*tprob
         elseif( x> xmin+1.d-3) then
            exit
         endif
         if( x > 0.1) then
            x = x + 5.d-3
         else
            x = x*10.0**0.025
         endif
      enddo
      write(0,*) 'v vds/dx(mb) Emu mbto??  sigma(mb)  prob/?? '
      end
