      implicit none
#include "Zprivate.h"
#include "Ztrack.h"



      integer num, cumnum, irevent(2), i,  nr, ndev, withir, nev
      type(track)::Zfirst
      character*160 str, msg
      character*80  skelin, stro(3), numbin


      integer klena, nlow

      Mdev = 71
      Rdev = 70
      ndev = 69

      str = ' '
      read(*,'(a)')  str
!
      do i = 1, 3
         stro(i) = ' '
      enddo
!
      call ksplit(str,  80, 3, stro,  nr)
      if(nr .lt. 3) stop 'must give 2 files and flag '

      skelin = stro(1)
      numbin = stro(2)
      read( stro(3), * )  withir


      open(Rdev, file=skelin(1:klena(skelin)), form='unformatted',
     *       status='old')

      open(ndev, file=numbin(1:klena(numbin)), form='formatted',
     *       status='old')


      do while(.true.)

         if(withir .eq. 0) then
            read(ndev, *, end= 1000) nev
         else
            read(ndev, *, end= 1000) irevent, nev
         endif
         Copy = .false.
         do while (.not. Copy)
            read(Rdev, end= 1000) cumnum, num, irevent, Zfirst
            if(cumnum .gt. nev) then
               write(msg,*) 'specified event #=',nev,
     *             ' not exist in skeleton-node file'
               call cerrorMsg(msg, 1)
               write(msg,*) 'skip to event=',cumnum
               call cerrorMsg(msg, 1)

               do while(nev .lt. cumnum)
                  if(withir .eq. 0) then
                     read(ndev, *, end= 1000) nev
                  else
                     read(ndev, *, end= 1000) irevent, nev
                  endif
               enddo
            endif
            Copy =cumnum .eq. nev
            if(Copy) then
               call toasciiH( cumnum, num, irevent, Zfirst)
            endif
            call cgetHES(Rdev)
            nlow = 1
            do while (nlow .ne. -1)
               read(Rdev)  nlow,  p
               if(Copy) then
                  call toasciiN(nlow,  p)
               endif
               do i = 1, nlow
                  read(Rdev) c
                  if(Copy) then
                     call toasciiC(c)
                  endif
               enddo
            enddo
         enddo
      enddo
 1000 continue
      call cerrorMsg(
     *   ' all events finished', 1)
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
