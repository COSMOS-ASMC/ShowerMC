      implicit none
#include "Zprivate.h"
#include "Ztrack.h"



      integer num, cumnum, irevent(2), i,  nr, ndev, withir, nev
      type(track):: Zfirst
      character*160 str, msg
      character*80  skelin, stro(4), numbin


      integer klena, nlow

      Mdev = 71
      Rdev = 70
      ndev = 69

      str = ' '
      read(*,'(a)')  str
!
      do i = 1, 4
         stro(i) = ' '
      enddo
!
      call ksplit(str,  80, 4, stro,  nr)
      if(nr .lt. 4)  then
         write(0,*) 'give name of ' 
         write(0,*) '        input SkelFile'
         write(0,*) '        event num file'
         write(0,*) '        ouput SkelFile'
         write(0,*) '        1/0 depending on event num file has  IR'
         write(0,*)
     *    'e.g echo  "skelfile Seed outskel 0" | randomselPCLinuxIFC'
         stop  111
      endif

      skelin = stro(1)
      numbin = stro(2)
      Mskel = stro(3)
      read( stro(4), * )  withir
      open(Rdev, file=skelin(1:klena(skelin)), form='unformatted',
     *       status='old')

      open(ndev, file=numbin(1:klena(numbin)), form='formatted',
     *       status='old')

      open(Mdev, file=Mskel(1:klena(Mskel)), form='unformatted',
     *       status='unknown')


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
      enddo
 1000 continue
      call cerrorMsg(
     *   ' all events copied', 1)
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
