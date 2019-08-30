      implicit none
#include "Zmedia.h"
#include "Zmass.h"
!
!          test compton
! 
       type(epmedia):: media
      integer i, io, imax, icon
      real*8 Ee, prob, Eg, cosg, cose, Egp, path
      character*120 file
      character(len=8):: name

      io = 10
      imax = 50000
      call cerrorMsg('compton sampling test', 1)
      call cerrorMsg(
     * " prob/r.l  Eg'  Ee   cosg  cose  sampled_path", 1)


      call  cerrorMsg(
     * "Enter gamma energy, events(50000) and media name"//
     * "BGO",  1)

      read(*, *) Eg, imax, name
      call cerrorMsg(name,  1)
      file = '$EPICSTOP/Data/BaseM/'//trim(name)
      call copenf(io, file, icon)

      call epReadMTbl(io, media)
      close(io)
      call epGetEffZA(media)
      write(*,*) '# "compton test"'
      write(*,*) 
     * '# "prob/r.l"    "Eg_out"  "Ee"   "cosg"  "cose"  "path"'
      do i = 1,  imax
         call epcompp(media, Eg, prob, path)
         call epcompea(Eg, Egp, Ee, cosg, cose)
         write(*,'(1p, 6g13.4)') prob, Egp, Ee, cosg, cose, path
      enddo
      end

