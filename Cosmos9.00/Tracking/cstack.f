!       stacking tracking information
!
      subroutine cpush(a)
      implicit none
#include  "Zmaxdef.h"
#include  "Ztrack.h"
#include  "Zstackv.h"
#include  "Zevhnv.h"
      type(track)::a
      character*70  msg

!

      if( a%p%code /=  0) then
         if( Stack_pos .ge. Max_stack_size) then
            write(msg,*) 'stack area full=',Max_stack_size
            call cerrorMsg(msg, 0)
         else
            Stack_pos = Stack_pos + 1
            Stack(Stack_pos) = a
         endif
      else
!           although very rare, 0 code appears;  neglect it (once /5 days)

         write(0,*) ' code=0  appeared subcode=',a%p%subcode
         write(0,*) ' charge=',a%p%charge
         write(0,*) ' px,py, pz, E, mass=',a%p%fm%p(1:4), a%p%mass
         write(0,*) ' neglected for stacking in cstack.f'
         write(0,*) 'ActiveMdl=',ActiveMdl
!c         call checkstat("in cpush")
      endif
      end
!------------------------------------------------------------
      subroutine cpop(a, remain)
      implicit none
#include  "Zmaxdef.h"
#include  "Ztrack.h"
#include  "Zstackv.h"
      type(track)::a

      integer remain
!
!         remain: int. is the number of ptcls remaining unprocessed
!                 including the current one to be processed now.
      if( Stack_pos .le. 0) then
         remain = 0
      else
         remain =  Stack_pos 
         a = Stack(Stack_pos)
         Stack_pos = Stack_pos -1
      endif
      end
      subroutine cgetStacked(stackpos, aTrack, icon)
      implicit none
#include  "Zmaxdef.h"
#include  "Ztrack.h"
#include  "Zstackv.h"
      integer stackpos  !  input.  stack pos from which track is
                        ! to be extracted.
      type(track)::aTrack  ! output extracted track.
      integer icon   ! output. 0 OK 1--> no. track at stackpos 
!          diff. from cpop;  Stack_pos is not affected.      
      if(stackpos .le. Stack_pos .and. stackpos .ge. 1) then
         aTrack = Stack(stackpos)
         icon = 0
      else
         icon = 1
      endif
      end

!---------------------------------------------------------
      subroutine cinitStack
!         initialize stack.
      implicit none
#include  "Zmaxdef.h"
#include  "Ztrack.h"
#include  "Zstackv.h"
      Stack_pos = 0
      end
      subroutine cgetCurrentStackPos(stackpos)
      implicit none
#include  "Zmaxdef.h"
#include  "Ztrack.h"
#include  "Zstackv.h"
      integer stackpos
      stackpos=Stack_pos
      end
      subroutine cresetStackPos(stackpos)
      implicit none
#include  "Zmaxdef.h"
#include  "Ztrack.h"
#include  "Zstackv.h"
      integer stackpos
      Stack_pos=stackpos
      end


      subroutine csortStack
!        sort stack dscending order 
      implicit none
#include  "Zmaxdef.h"
#include  "Ztrack.h"
#include  "Zstackv.h"
      real*8 erg(Max_stack_size)
      type(track)::sorted(Max_stack_size)
      integer  idx(Max_stack_size)
      
      integer i, j, k
      type(track)::temp
      do  i = 1,  Stack_pos
         erg(i) = Stack(i)%p%fm%p(4)
      enddo
      
      call kqsortd(erg, idx, Stack_pos)
      call ksortinv(idx, Stack_pos)
      
      do i = 1, Stack_pos
         sorted(i) = Stack(idx(i))
      enddo

      do i = 1, Stack_pos
         Stack(i) = sorted(i)
      enddo
      end

      

      
