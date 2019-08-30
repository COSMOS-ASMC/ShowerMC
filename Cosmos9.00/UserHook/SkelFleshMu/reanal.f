      implicit none
#include "Zprivate.h"
#include "Ztrack.h"


      integer num, i,  nr, firstev, lastev
      type(track):: Zfirst
      character*160 str
      character*80  skelin, stro(3)
      
      real*8  temp

      integer klena, nlow

      integer  cumnum, irevent(2)
      common /cselev/ cumnum, irevent

      Rdev = 70

      str = ' '
      read(*,'(a)')  str
!
      do i = 1, 3
         stro(i) = ' '
      enddo
!
      call ksplit(str,  80, 3, stro,  nr)
      if(nr .lt. 3) stop 'must give 1 files and 2 event #'

      skelin = stro(1)
      read( stro(2), * ) firstev
      read( stro(3), * ) lastev

      open(Rdev, file=skelin(1:klena(skelin)), form='unformatted',
     *       status='old')

      do while(.true.)
         read(Rdev, end= 1000) cumnum, num, irevent, Zfirst

         if(cumnum .gt. lastev) goto 1000

         call cgetHES(Rdev)
         nlow = 1
         do while (nlow .ne. -1)
            read(Rdev)  nlow,  p
            do i = 1, nlow
               read(Rdev) c
            enddo
         enddo
      enddo
 1000 continue
      if(lastev .gt. cumnum) then
         call cerrorMsg(
     *   'EOF before the specified last event reached; ok?', 1)
      else
         call cerrorMsg(
     *   ' all events processed  successfully', 1)
      endif

      end
      

      subroutine cgetHES(from)
      implicit none
#include "Zprivate.h"
      integer from

      integer  cumnum, irevent(2)
      common /cselev/ cumnum, irevent

!
      integer i
      read(from) Np
      do i = 1, Np
         read(from) o(i)
      enddo
!     ==========================================
!         do your analysis and if you want to
!         memorize the it for flesh, make 
!         Copy = .true.  Then, the event number
!         (or seed  and event number) is written on stdout.
!       You can use
!     o(i).xxx  where  xxx is
!
!       where: obsv. location
!       code, subcode, charge:  code=9 is Nucleus, A=subcode,Z=charge
!       atime: arrival time ns
!       erg:   energy GeV
!       mass:  mass in GeV
!       x,y:   in m
!       wx,wy,wx:  direction cos.
!       zenith:    cos of zenith angle.
!

      Copy = .true.


!     ==========================================
      if(Copy) then
          write(*,*) cumnum  
!            or 
!         write(*,*) irevent, cumnum
       endif
      end
