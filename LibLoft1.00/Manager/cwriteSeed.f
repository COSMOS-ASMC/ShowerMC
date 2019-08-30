!     ****************************************************
!          inquire the initial random seed
      subroutine  cqIniRn(ir)
      implicit none
#include  "Zmanager.h"
      integer ir(2)
      
      ir(1) = SeedSave(1)
      ir(2) = SeedSave(2)
      end
!     ****************************************************
!          write iformation for initial random number for flesh
!        at skelton making time
      subroutine cwriteSeed
      implicit none
#include "Zmanagerp.h"      
#include "Zmanager.h"      

#undef  F90
#undef FLUSH

#if defined PCLinux
#define FLUSH
#define F90
#elif defined PCLinuxIFC
#define FLUSH
#define F90
#elif defined Solaris
#define FLUSH
#elif defined NEXT486
#define FLUSH
#elif defined DECALPHA
#define FLUSH
#elif defined CF_AlphaLinux
#define FLUSH
#elif defined SGI
#define FLUSH
#elif defined IBMAIX
#define FLUSH
#endif

#if defined FLUSH
      character*120 fn
      logical first
      integer klena
      save first, fn
      data first/.true./
      
      if(first) then
         first = .false.
         fn = ' '
         call cgetfname(SeedFile, fn)
      else
#if defined F90
         open(SeedFileDev, file=fn(1:klena(fn)), 
     *    form='formatted', position='append')
#elif defined IBMAIX
!               ibm should apppend if status old
         open(SeedFileDev, file=fn(1:klena(fn)), 
     *    form='formatted', status='old')
#else
         open(SeedFileDev, file=fn(1:klena(fn)), 
     *    form='formatted', access='append')
#endif
      endif
#endif

      write(SeedFileDev, *)  SeedSave, EventNo+1

#if defined FLUSH
      close(SeedFileDev)
#endif
      end
!     ****************************************************
!          read iformation for initial random number at flesh time
      subroutine creadSeed(ir, no, jeof)
      implicit none
#include "Zmanagerp.h"      

      integer ir(2), no, jeof

      integer nn
      integer count
      data count/0/
      save count

      jeof = 0
      do while(.true.)
         read(SeedFileDev, *, end=100) ir, nn
         if(nn .ge. EventNo) goto 200  ! skip already processed ones
      enddo

 100  continue
      jeof =1
      if(count .eq. 0) then
         count = count + 1
         call cerrorMsg(
     *        'No more events to be fleshed',1)
         call cerrorMsg(
     *        'EOF of the seed file reached', 1)
      else
         call cerrorMsg(
     *  '2nd time to read EOF of the  seed file',0)
      endif
 200  continue
      no = nn
      end
