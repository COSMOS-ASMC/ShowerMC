      subroutine eppush(aTrack)
      implicit none
!#include  "ZepTrack.h"
#include  "ZepTrackv.h"
#include  "ZepStack.h"
#include  "ZsepManager.h"
#include  "Zcode.h"
      
       type(epTrack)::  aTrack ! input ptcl to be stacked
      integer:: size
      integer i
      character(len=192)::stackdiskpath

       type(epTrack)::  bTrack ! copy of aTrack
#if defined (jaxa) || (jaxaflat) 
       type(epTrack):: dummyTrack(2)   ! to get size
#endif

      bTrack = aTrack
      if( Making1ry == 1 ) then
         ! now 1ry is being stacked
         Total1ryE =Total1ryE + bTrack%p%fm%p(4) -bTrack%p%mass
         bTrack%inciflag = 1
      else
         if( Move%Track%inciflag == 0 )  then
            bTrack%inciflag = 0 ! bulk of cases come here
!                  so separated from  the next subrotuine
         else
            call epManageInciflag(aTrack, bTrack) ! in epSetupV1ry
         endif
      endif

      if(Stack_pos .ge. MaxStackSize) then
         if( DdiskNo == 0 ) then
           ! open direct access file in binary 
            DdiskNo = StackDiskFileNo
#if defined (jaxa) || (jaxaflat) 
!  this is not used by  ifort;  loc is non standard but
!  loc seems to be usable by ifort too. 
            size = loc( dummyTrack(2)) -loc( dummyTrack(1))
#else
            size = sizeof(bTrack)      ! in bytes      
#endif
            if(StackDiskFile == "scratch" ) then
               open(DdiskNo,  status="scratch",
     *              form="unformatted", access="direct",
#if defined (Solaris) || (jaxa) || (jaxaflat) || (MacGfort) || (LinuxGfort)
     *           RECL=size)   
#else
     *           recordsize=size)   
#endif
            else
               call cgetfname(StackDiskFile, stackdiskpath)
               open(DdiskNo,  file=stackdiskpath,
     *           form="unformatted", access="direct",
#if defined Solaris || (jaxa) || (jaxaflat) || (MacGfort) || (LinuxGfort)
     *           RECL=size)   
#else
     *           recordsize=size)   
#endif
            endif
            write(0,*) 'memory stack area=',MaxStackSize
            write(0,*) ' is not enought so that direct '
            write(0,*) ' access file has been opened: '
            if(StackDiskFile == "scratch" ) then
               write(0,*) 
     *          ' it is a scratch file and the system created it'
               write(0,*)
     *          ' somewhere (probaby in /tmp/$USER/) and will '
               write(0,*)
     *          ' delete it at end of job even if abnormal end '
            else
               write(0,*) trim(stackdiskpath)
               write(0,*) 'This will be deleted at normal end of job'
               write(0,*)
     *           'for abnormal end, you must delete it yourself'
            endif
         endif
         Diskpos = Diskpos + 1
!//////////////////////
!         if( bTrack.cn == 0 ) then
!            write(0,*) ' push cn=0. code =',bTrack.p.code
!            write(0,*) ' pos =',bTrack.pos.x,bTrack.pos.y,bTrack.pos.z
!         endif
!///////////////////
         write(DdiskNo, rec=Diskpos) bTrack
      else
         Stack_pos = Stack_pos + 1
         Stack(Stack_pos) = bTrack
      endif
      end  subroutine eppush


   !       ************************
   !        put 1ptcl in the stack area (for unknown cn)
   !       if stack overflows, disk is automatically used.
      subroutine epputTrack(aTrack)
      implicit none
#include "ZepTrack.h"

#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(epTrack)::  aTrack   ! input. 

!                 find component # for pos (world coord.)
!              and  set it in aTrack.cn
      call eppos2cn(0,  aTrack,  aTrack%cn)
!         convert pos/w into  l.c
      if(aTrack%cn .le. Det%nct ) then

         call epw2l(aTrack%cn, aTrack%pos, aTrack%pos)
!c           call epw2ld(aTrack.cn, aTrack.w,  aTrack.w)
         call epw2ldm(aTrack%cn, aTrack%w,  aTrack%w, aTrack%p)
         call eppush(aTrack)
      else
         write(0,*) ' current comp. #=',aTrack%cn,
     *        ' total comp.#=', Det%nct
         write(0,*) ' epputTrack: ptcl is not in detector'
         stop
      endif
      end subroutine epputTrack
!       ************************
!         put 1 ptcl in the stack disk (for unknown cn)
!     
      subroutine epputTrk2(io, aTrack)
      implicit none
#include "ZepTrack.h"

      integer io               !  disk file number to stack the track
       type(epTrack)::  aTrack   ! input. 
!                 find component # for pos (world coord.)
!              and  set it in aTrack.cn
      write(0,*) " You called  epputTrk2 "
      write(0,*) " please use 'call epputTrack(aTrack)' "
      write(0,*) " you need not open the disk "
      write(0,*) " However, you must define the file name path in "
      write(0,*) " sepicsfile thru StackDiskFile and it logical #"
      write(0,*) " by StackDiskFileNo (default is 13) "
      stop
      end subroutine epputTrk2
!      ******************* 
      subroutine epempty
      implicit none
!             clear the stack
#include  "ZepTrack.h"
#include  "ZepStack.h"
      logical,save::first=.true.

      Stack_pos = 0
      Diskpos = 0
      if(first) then
         DdiskNo = 0
         first =.false.
      endif

      end  subroutine epempty
!     ******************************
      subroutine eppop(aTrack, icon)
      implicit none
!             get 1 particle from the stack area
!              icon=0:  1 ptcl gotten
!                  -1:  no more ptcl
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"
#include  "ZepStack.h"       
!
       type(epTrack)::  aTrack
       integer::icon
!/////////////
!       integer onlyonce
!       common /dddebug/ onlyonce
!//////////

       icon=0
       do while( .true. )
          if( Diskpos > 0 ) then
             read(DdiskNo, rec=Diskpos ) aTrack
             Diskpos = Diskpos - 1
          elseif(Stack_pos > 0)  then
!                vectors are in local.c
             aTrack = Stack(Stack_pos)     
             Stack_pos = Stack_pos - 1
          else
             icon  = -1
             exit
          endif
!                make cosmos ptcl; e may not be informed. g may  not be,
!                                  if IncGp=0.
          if( aTrack%p%code > 0 ) then !   <<<<>>>> light
             call ep2cosPtcl(aTrack%p)
          endif                 !   <<<<>>>> light

!              do business when a new component appear
          call epnewComp(aTrack)
          if( Light > 0 .and. aTrack%p%code < 0 ) then
                ! special routine for light related ptlcs.
             call epLightAtPop(aTrack, icon)
             if( icon /=  0 ) cycle
          endif
!***        call epchckE(aTrack, icon)
!***        if( icon /= 0 ) cycle
          exit
       enddo
       end  subroutine eppop
!       *************************
      subroutine epqsTrack(n, sTrack)
      implicit none
#include  "ZepTrack.h"
      integer,intent(in):: n  ! stack pos
       type(epTrack)::  sTrack   ! output

      integer:: icon
      call epgetTrack(n, sTrack, icon)
      if(icon /= 0 ) then
         sTrack%p%fm%p(4) = 0.
      endif
      end subroutine epqsTrack

      subroutine epgetTrack(n, aTrack, icon)
      implicit none
  !           returns a specified track in the stack
  !           if n > memory stack size, disk stack is
  !           sought for.
#include  "ZepTrack.h"
#include  "ZepStack.h"

!
      integer n                 ! input. track in n-th stack position
       type(epTrack)::  aTrack   ! output
      integer icon              ! output. 0--> ok
                                 !        non 0--> no n-th stack
      if(n >= 1 .and. n <= Stack_pos) then                   
         aTrack = Stack(n)
         icon = 0
      elseif(Stack_pos == 0 ) then
         !  no ptcl stacked
         icon =1 
      elseif( DdiskNo > 0 .and. Diskpos > 0 ) then
         if( Diskpos > n-MaxStackSize  .and. 
     *       n-MaxStackSize > 0 ) then 
            read(DdiskNo, rec = n-MaxStackSize) aTrack
            icon = 0
         else
            icon = 1 
         endif
      else
         icon  = 1
      endif
      end subroutine epgetTrack
!       *************************
      subroutine epretTrack(n, aTrack, icon)
      implicit none
!            inverse  of epgetTrack
!            returns aTrack  to the n-th pos in the stack
!                 
#include  "ZepTrack.h"
#include  "ZepStack.h"
!
      integer n                 ! input. track in n-th stack position
       type(epTrack)::  aTrack   ! input. 
      integer icon              ! output.  0--> ok
!          1--> n is outside of the region
      if(n >=  1 .and. n <= Stack_pos) then                   
         Stack(n)= aTrack 
         icon = 0
      elseif( Diskpos > n-MaxStackSize ) then
         write(DdiskNo, rec = n-MaxStackSize) aTrack
         icon = 0
      else
         icon = 1
      endif

      end   subroutine epretTrack


      subroutine epqstn(n)
      implicit none
#include  "ZepTrack.h"
#include  "ZepStack.h"
      integer,intent(out)::n  ! current stack #
      n = Stack_pos
      end  subroutine epqstn
      subroutine epqstn2(ns,nd)
      implicit none
#include  "ZepTrack.h"
#include  "ZepStack.h"
      integer,intent(out)::ns  ! current stack #
      integer,intent(out)::nd  ! current disk pos, if ns==MaxStackSize else 0
      ns = Stack_pos
      nd = Diskpos
      end  subroutine epqstn2

      subroutine epcloseStackDisk
      implicit none
#include  "ZepTrack.h"
#include  "ZepStack.h"
      if( DdiskNo > 0 ) then
         close(DdiskNo, status="delete")
      endif
      end
