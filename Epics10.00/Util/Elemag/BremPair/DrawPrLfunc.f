#include "ZepicsBD.h"
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZepTrackp.h"
!
      type(epmedia):: media

      real*8 Ee,  Eg
      character*80 file
      character*24 name
      real(8),external:: epPairLowE
      real(4):: rhoc
      real*8 Egme, xs               ! input  Eg/me
      integer  icon, io
      real(8):: Z    ! Atom's Z
      real*8 x        ! input  Ee/Eg.   me/Eg =< x <= 1.-me/Eg   
      io = 10
      write(0,*) 'Enter media name, Z'
      read(*, *)  name
      call epgetRhoc(name, name, rhoc)
      file = "$EPICSTOP/Data/Media/"//trim(name) 
      call cerrorMsg(file,  1)
      call copenf(io, file, icon)
      call epReadTab(io, media)
      Z = media%Z
      media%rhoc = rhoc

      Eg= 2.e-3
      do while (Eg < 10000.)
         Egme = Eg/masele
         x=masele/Eg
         Ee = x*Eg
         do while (x< 1.0-masele/Eg)
            xs =epPairLowE(Z, Egme, x)
            write(*,'(1p,4g13.3)') x, xs, Ee, Eg
            x = x + 0.002d0
         enddo
         Eg = Eg*10.
      enddo
      end
