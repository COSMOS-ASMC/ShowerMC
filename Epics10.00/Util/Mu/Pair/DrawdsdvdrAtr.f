      implicit none
!
!         to write data needed to draw muon pari creation function 
!
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZmuBPNgene.h"


      character*50  file
      integer nerg
      integer io, i, norm, ie, icon
      real*8  x, f, xmax
      real*8  xmin, tprob, rho, epmuPairVmn, epmudsdvdr

      parameter (nerg=10)
      real*8  EmuA(nerg)

      io = 10
      call cerrorMsg(' ', 1)
      call cerrorMsg('1) In unit of / r.l (1)',  1)
      call cerrorMsg('2) In unit of mb/ingredient', 1) 
      call cerrorMsg('3) In unit of /(g/cm^2)', 1)
      call cerrorMsg('4) In unit of /cm', 1)
!      call cerrorMsg(
!     * '5) Area-normalization '//
!     * ' for comparison with M.C data such as v vs vdN/dv/Nt  ',
!     * 1 )

      
      norm = 1
      read(*, *)  norm
      do i = 1, nerg
         EmuA(i) = 0.
      enddo
      call cerrorMsg(
     *  "Enter a basic media file path ( such as"//
     *    "  '../../Data/BaseM/Pb')",  1)
      file ='../../Data/BaseM/Pb'
      read(*, *) file
      call cerrorMsg(
     * '  Energy of the muon in (GeV) upto 10 (with /)', 1)
      read(*, *) EmuA

      call copenf(io, file, icon)
      if(icon .ne. 0) then
         call cerrorMsg('media file path invalid',  0)
      endif
!
      call cerrorMsg('Enter assymetry factor ', 1)
      read(*,*) rho
      call epmuBPNgeneI(io, media)

!     **********
      call epmuSetCnst(media%elem(1)%Z, media%elem(1)%A)
!     *********


      do ie = 1, nerg
         Emu = EmuA(ie)
         if(Emu .eq. 0.) goto 10
!         xmin =max(media.cnst.muPrVmin, epmuPairVmn(Emu))
         xmin = epmuPairVmn(Emu)

         if(norm .eq. 5) then
!            call eptotcmuB(xmin, xmax, tprob)
            tprob =  1.
         else
            tprob =  1.
         endif

         x = xmin

!         xmax = epmuvmax(Emu)             
!          cannot be determined here since it is dependent
!          on each element. So we set it  1.0 tentatively
!          actual value is set in  epmuBremS
!          and  if x> xmax, 0 is returned
!
         xmax = 1.d0

         do while (x .lt.  xmax)
            f = epmudsdvdr(Emu, x, rho)/x

            if(norm .eq. 5) then
               f = f/tprob
            else
               if(norm .eq. 1) then 
                  f = f * media%mbtoPX0 ! prob/X0
               elseif( norm .eq.  2) then
!                 nothing to do  . in mb
               elseif(norm .eq.  3) then
                  f = f* media%mbtoPgrm
               elseif(norm .eq. 4) then
                  f = f* media%mbtoPcm
               else
                  call cerrorMsg('input error',0)
               endif
            endif
!
            if(f .gt. 0.) then
               write(*,'(5g13.3)')
     *         sngl(x), sngl(f*x),  sngl(f), sngl(Emu), sngl(rho)
            elseif( x .ge. .8) then
               goto  5
            endif
            if(x .gt.  xmax-0.5d0 ) then
               x = x + 1.d-3
            else
               x = x * 10.**0.02
            endif
         enddo
 5       continue
         write(*,*)
      enddo
 10   continue
      end