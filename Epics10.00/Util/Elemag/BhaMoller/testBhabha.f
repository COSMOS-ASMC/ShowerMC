      implicit none
#include "Zmedia.h"
#include "Zmass.h"
!
!          test Bhabha scattering
! 
       type(epmedia):: media
      integer i, io, icon
      real*8 Ep, prob, Epp, cosp, cose, Ee, w, path
      character*120 file
      character*12 name

      io = 10
      
      call  cerrorMsg(
     * 'Enter positron energy(T.E, 1GeV), cutoff kinetic'//
     * ' energy(100.d-6 GeV)  and media name(BGO)', 1)
      name = "W"
      Ep = 1.
      w = 100.d-6
      read(*, *) Ep, w, name
      call cerrorMsg(name,  1)
      file = "$EPICSTOP/Data/Media/"//trim(name)
      call copenf(io, file, icon)
      call epReadTab(io, media)
      close(io)
      call epGetEffZA(media)
      call epStern(media)  !  argument is now 1. 

      write(*,'(a)')
     * "# prob/r.l sampled:  E+'   E-/E+    cos+   cos-    path(r.l)"
      do i = 1,  100000
         call epbhabhap(media, Ep, w,  prob, path)
         call epbhabhae(Ep, w, Epp, Ee, cosp, cose)
         write(*,'(1p,6g12.4)') prob, Epp, (Ee-masele)/(Ep-masele),
     *     cosp, cose, path
      enddo
      end
