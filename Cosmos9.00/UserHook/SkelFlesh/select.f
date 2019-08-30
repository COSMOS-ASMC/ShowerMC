      implicit none
#include "Zprivate.h"
#include "Ztrack.h"



      integer num, cumnum, irevent(2), i,  nr, firstev, lastev
      type(track)::Zfirst
      character*160 str
      character*80  skelin, stro(4)

      real*8  temp

      integer klena, nlow

      Mdev = 71
      Rdev = 70

      str = ' '
      read(*,'(a)')  str
!
      do i = 1, 4
         stro(i) = ' '
      enddo
!
      call ksplit(str,  80, 4, stro,  nr)
      if(nr .lt. 4) stop 'must give 2 files and 2 event #'

      skelin = stro(1)
      read( stro(2), * ) firstev
      read( stro(3), * ) lastev
      Mskel = stro(4)


      open(Rdev, file=skelin(1:klena(skelin)), form='unformatted',
     *       status='old')
      open(Mdev, file=Mskel(1:klena(Mskel)), form='unformatted',
     *       status='unknown')

      do while(.true.)
         read(Rdev, end= 1000) cumnum, num, irevent, Zfirst

         if(cumnum .gt. lastev) goto 1000

         Copy =cumnum .ge. firstev .and. cumnum .le. lastev
         if(Copy) then
            write(Mdev)  cumnum, num, irevent, Zfirst
         endif
         call cgetHES(Rdev)
         nlow = 1
         do while (nlow .ne. -1)
            read(Rdev)  nlow,  p
            if(Copy) then
               write(Mdev)  nlow,  p
            endif
            do i = 1, nlow
               read(Rdev) c
               if(Copy) then
                  write(Mdev) c
               endif
            enddo
         enddo
      enddo
 1000 continue
      if(lastev .gt. cumnum) then
         call cerrorMsg(
     *   'EOF before the specified last event reached; ok?', 1)
      else
         call cerrorMsg(
     *   ' all events copied  successfully', 1)
      endif

      end
      

      subroutine cgetHES(from)
      implicit none
#include "Zprivate.h"
      integer from
!

      integer i

      read(from) Np
      if(Copy) then
         write(Mdev) Np
      endif
      do i = 1, Np
         read(from) o(i)
         if(Copy) then
            write(Mdev) o(i)
         endif
      enddo

      end
