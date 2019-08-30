      implicit none
#include "Zprivate.h"
#include "Ztrack.h"



      integer num, cumnum, irevent(2), i,  nr, firstev, lastev
      type(track)::Zfirst
      character*160 str
      character*80  skelin, stro(3)

      real*8  temp

      integer klena, nlow

      Mdev = 71
      Rdev = 70

      str = ' '
      read(*,'(a)')  str
!
      do i = 1, 3
         stro(i) = ' '
      enddo
!
      call ksplit(str,  80, 3, stro,  nr)
      if(nr .lt. 3) stop 'must give 1 file and 2 event #'

      skelin = stro(1)
      read( stro(2), * ) firstev
      read( stro(3), * ) lastev


      open(Rdev, file=skelin(1:klena(skelin)), form='unformatted',
     *       status='old')

      do while(.true.)
         read(Rdev, end= 1000) cumnum, num, irevent, Zfirst

         if(cumnum .gt. lastev) goto 1000

         Copy =cumnum .ge. firstev .and. cumnum .le. lastev
         if(Copy) then
            call toasciiH(cumnum, num, irevent, Zfirst)
         endif
         call cgetHES(Rdev)
         nlow = 1
         do while (nlow .ne. -1)
            read(Rdev)  nlow,  p
            if(Copy) then
               call toasciiN( nlow,  p)
            endif
            do i = 1, nlow
               read(Rdev) c
               if(Copy) then
                  call toasciiC(c)
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
     *   ' all events processed successfully', 1)
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
         call toasciiNp(Np)
      endif
      do i = 1, Np
         read(from) o(i)
         if(Copy) then
            call toasciiHE(o(i))
         endif
      enddo

      end
