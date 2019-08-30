      implicit none
!
!         to write data needed to draw muon brems total Xsect.
!         dE/dx/, dEdx/E
!
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZmuBPNgene.h"

      character*50  file
      integer io,  norm, ie, icon, frac, nerg
      real*8  xmax
      real*8  xmin, tprob
      real*8  dEdx, dEdx0, dEdx1


      nerg = 75
      io = 10
      call cerrorMsg(' ', 1)
      call cerrorMsg('1) In unit of / r.l (1)',  1)
      call cerrorMsg('2) In unit of mb/ingredient', 1) 
      call cerrorMsg('3) In unit of /(g/cm^2)', 1)
      call cerrorMsg('4) In unit of /cm', 1)
      
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

      xmin = media%cnst%muBrVmin
      Emu = 1.
      write(*,*)
     * '# Emu(GeV), dE/dx, dE/dx(v<vm), dE/dx(v>vm) ToalXS'
      do ie = 1, nerg
!         xmax = epmuvmax(Emu)             
!          cannot be determined here since it is dependent
!          on each element. So we set it  1.0 tentatively
!          actual value is set in  epmuBremS
!          and  if x> xmax, 0 is returned
!
         xmax = 1.d0
         call eptotcmuB(xmin, xmax, tprob)
         call epmuElossB(xmin, xmax, dEdx1)
         call epmuElossB(0.d0, xmin, dEdx0)
         if(frac  .eq. 1) then
            dEdx0 = dEdx0* Emu
            dEdx1 = dEdx1* Emu
         endif

         if(norm .eq. 1) then 
            tprob = tprob * media%mbtoPX0 ! prob/X0
            dEdx0 = dEdx0 * media%mbtoPX0 ! /X0
            dEdx1 = dEdx1 * media%mbtoPX0 ! /X0
         elseif( norm .eq.  2) then
!            nothing to do  . in mb
         elseif(norm .eq.  3) then
            tprob = tprob* media%mbtoPgrm
            dEdx0 = dEdx0 *media%mbtoPgrm
            dEdx1 = dEdx1 *media%mbtoPgrm
         elseif(norm .eq. 4) then
            tprob = tprob* media%mbtoPcm
            dEdx0 = dEdx0 *media%mbtoPcm
            dEdx1 = dEdx1 *media%mbtoPcm
         else
            call cerrorMsg('input error for norm',0)
         endif
         dEdx = dEdx0 +  dEdx1
         write(*,'(5G13.3)') sngl(Emu), sngl(dEdx),
     *           sngl(dEdx0), sngl(dEdx1), tprob
         Emu = Emu * 10.d0**0.1
      enddo
      end
