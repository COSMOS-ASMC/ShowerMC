      implicit none
!
!                Int(vmin:v) of  v* ds/dv is obtained 
!        for muon pair creation.
!          vm  .vs.  Int v*ds/dv is written  on  stdout.
!
!             
!
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZmuBPNgene.h"

      character*50  file
      integer io,  norm,  icon, frac
      real*8  vmax
      real*8  vmin, total, v, v1, ans,  ans1
      real*8  dEdx, epmuPairVmn



      io = 10
      call cerrorMsg(' ', 1)
      call cerrorMsg('1) In unit of / r.l (1)',  1)
      call cerrorMsg('2) In unit of mb/ingredient', 1) 
      call cerrorMsg('3) In unit of /(g/cm^2)', 1)
      call cerrorMsg('4) In unit of /cm', 1)
      call cerrorMsg('5) Area normalization', 1)
      
      norm = 1
      read(*, *)  norm
      call cerrorMsg(
     *  "Enter a basic media file path ( such as"//
     *    "  '../../Data/BaseM/Pb')",  1)
      file ='../../Data/BaseM/Pb'
      read(*, *) file

      call copenf(io, file, icon)
      if(icon .ne. 0) then
         call cerrorMsg('media file path invalid',  0)
      endif

      call cerrorMsg('1) dE/dx in absolute scale OR', 1)
      call cerrorMsg('2) dE/dx/E (defalut)',1)
      frac = 2
      read(*,*) frac
      

      call epmuBPNgeneI(io, media)

!      vmin = media.cnst.muPrVmin

      call cerrorMsg('Enter muon energy (GeV)', 1)
      read(*, *) Emu
      vmin = epmuPairVmn(Emu)
      if(norm .eq. 5) then
         call epmuElossP(vmin, 1.d0, total)
!             this is fractional energy loss
      endif
      write(*,*)
     * '# v,  dE/dx(v<vm)'
      v = vmin
      ans = 0.
      do while(v .lt. 0.9999d0)
         v1 = v
         v = v*10.d0**0.025d0
         vmax = min(1.d0, v)
         call epmuElossP(v1, vmax, ans1)
         ans = ans + ans1
         dEdx = ans
         if(frac  .eq. 1 .and. norm  .ne. 5) then
            dEdx = dEdx* Emu
         endif

         if(norm .eq. 1) then 
            dEdx = dEdx * media%mbtoPX0 ! /X0
         elseif( norm .eq.  2) then
!            nothing to do  . in mb
         elseif(norm .eq.  3) then
            dEdx = dEdx *media%mbtoPgrm
         elseif(norm .eq. 4) then
            dEdx = dEdx *media%mbtoPcm
         elseif(norm .eq.  5) then
            dEdx  = dEdx/total
         else
            call cerrorMsg('input error for norm',0)
         endif
         write(*,'(3G13.3)') sngl(v), sngl(dEdx),
     *           sngl(Emu)
      enddo
      end
