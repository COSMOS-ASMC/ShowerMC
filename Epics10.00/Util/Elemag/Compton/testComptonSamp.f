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
      character(len=8)::name

      io = 10
      write(0,*) ' Enter Eg(GeV), # of samplig, media name'
      read(*, *) Eg, imax, name


      file  = "$EPICSTOP/Data/Media/"//trim(name)

      call copenf(io, file, icon)
      if(icon /= 0 ) then
         write(0,*) ' open file=',trim(file), 'error'
         stop
      endif

      call epReadMTbl(io, media)
      write(0,*) ' media%name=', media%name, ' read'
      close(io)
      call epGetEffZA(media)
      write(0,*)
     * 'compton test; see also data in /tmp/$USER/Work'
      write(0,*) 'comptonX.dat which contains data for'
      write(0,*) 
     * '# "Egp/Eg"   "Ee"  "cosg"   "cose"  "path"  "Eg"'
      write(*,*)
     * '# "Egp/Eg"   "Ee"  "cosg"   "cose"  "path"  "Eg"'
      do i = 1,  imax
         call epcompp(media, Eg, prob, path)
         call epcompea(Eg, Egp, Ee, cosg, cose)
         write(*,'(1p,2g12.4, 0p, 2f10.6,1p,3g12.4)')
     *    Egp/Eg, Ee, cosg, cose, path, prob, Eg
      enddo
      end
