      implicit none
#include "Zmedia.h"
#include "Zmass.h"
!
!          test Bhabha scattering cross-section
!     At fixed e+ energy, see the cross-section
!     change as a function of cut off recoil electron
!     energy.
! 
       type(epmedia):: media
      integer i, io, icon
      real*8 Ee, prob,  w, path, v, vm
      real*8 xm ,   g
      character*130 file
      character(len=12):: name
      real*8  epbhabhatx,  tx, txg, t0, beta2

      io = 10
      
      call  cerrorMsg(
     * "Enter  e+ energy(T.E 0.1 GeV),"//
     * " minimum cutoff  recoil energy(100e-6 GeV)"//  
     * " and media name(BGO)", 1)
      Ee = 1.
      w = 100.0d-6
      name = "BGO"

      
      read(*, *)  Ee, w, name
      call cerrorMsg(name,  1)
      file = '$EPICSTOP/Data/Media/'//trim(name)
      call copenf(io, file, icon)
      call epReadTab(io, media)
      close(io)
      call epGetEffZA(media)
      call epStern(media)  ! only 1 arg
      g = Ee/masele
      t0 = Ee-masele
      beta2 = 1.-1/g/g
      write(*,'(a)')
     *  "# E-min     prob/r.l     tx   txg/r.l"
      do while ( (Ee-masele) .gt. w )
         vm =  w/t0
         call epbhabhap(media, Ee, w,  prob, path)
         tx = epbhabhatx(g, vm)
         txg = tx * masele/t0/beta2 * media%basearea*2.0  ! /r.l
         write(*,'(1p,4g12.4)')  w, prob, tx, txg
         w = w *10.**0.01
      enddo
      end
