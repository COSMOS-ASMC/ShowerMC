      implicit none
#include "Zmedia.h"
#include "Zmass.h"
!
!          test Moller scattering
! 
       type(epmedia):: media
      integer i, io, icon
      real*8 Ee, prob, Eep, cose, cosr, Er, w, path
      character*130 file
      character(len=12):: name

      io = 10
      
      call  cerrorMsg(
     * 'Enter electron energy(T.E GeV),'//
     * ' cutoff energy(K.E 100.d-6 GeV)'//
     * ' and media nameh(BGO)', 1)
      w = 100.d-6
      name = "BGO"

      
      read(*, *) Ee, w, name
      call cerrorMsg(name,  1)
      file = '$EPICSTOP/Data/Media/'//trim(name)
      call copenf(io, file, icon)
      call epReadTab(io, media)
      close(io)
      call epGetEffZA(media)
      call epStern(media)       ! one argm. now
      write(*,'(a)')
     *   "#  prob   sampled:  Es    KEr/KEe  coss  cosr  path (r.l)"
      
      do i = 1,  100000 
         call epmollerp(media, Ee, w,  prob, path)
         call epmollerea(Ee, w, Eep, Er, cose, cosr)
         write(*,'(1p,6g12.4)') prob, Eep, (Er-masele)/(Ee-masele),
     *    cose, cosr, path
      enddo
      end
