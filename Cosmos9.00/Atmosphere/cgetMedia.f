      subroutine cgetMedia
!     This is almost same as epgetMedia
      use modAtmosDef
      implicit none
#include "ZmediaLoft.h"
      
      character(30):: msg
      logical fatal
      integer i, j, icon

      character*8  missingMedia(Maxmedia)
      integer  missing, k

      character(8)::  tempc

      missing = 0
      MediaNo = 0
      fatal = .false.

      do i = 1,  atmos%nodes
         tempc = atmos%matter(i)
         if( tempc =="sp2" ) tempc="sp"
!            for the moment don't use this;  (see eprcnf.f)
!     call epSeeIfAlias(tempc) ! if alias, replace it by true
         do j = 1,  MediaNo         
            if(tempc .eq. Media(j)%name) then ! v9.14
!               aleady read
!                i-th component has media index j.
               atmos%node2mediaNo(i) = j
               goto 10
            endif
         enddo
!            this is new matter
         if(MediaNo .lt. Maxmedia) then
            MediaNo = MediaNo + 1
            j = MediaNo
         else
            write(msg, *) Maxmedia
            call cerrorMsg(
     *       'too many different media > '//msg, 0)
         endif
!            see if the matter file exists, and if exist, read it
         call eprdMFile(tempc, icon) ! file is in Preinte
!     inside above;execute  epReadTab (epReadMTbl,epGetEffZA,epExpot
!            reading-many-tables, epRdmuTab, epSetMu)
!            epStern, epSetUrban; epReadXXsec   epSrimChk         
!                       
         if(icon .ne. 0) then
!            see if already appeared
            do k = 1, missing
!!!               if(missingMedia(k) .eq. Det.cmp(i).matter) then
               if(missingMedia(k) .eq. tempc ) then  ! v9.14
                  goto 5
               endif
            enddo
            missing = missing + 1
            if(missing .le. Maxmedia) then
               write(msg, *) 'Media file:',
!!     *          Det.cmp(i).matter, ' not exist'
     *          tempc, ' not exist'   ! v9.14
               call cerrorMsg(msg, 1)
!!               missingMedia(missing) = Det.cmp(i).matter
               missingMedia(missing) = tempc   !  v9.14
            else
               call cerrorMsg(
     *              'Seems too many missing media names',1)
               missing = Maxmedia
            endif
 5          continue
            
            MediaNo = MediaNo - 1
            fatal = .true.
         else
!                i-th node has media index j.
            atmos%node2mediaNo(i) = j
         endif
 10      continue
      enddo

      if(fatal) then
         call cerrorMsg(
     *    'Media files shown above are missing in any of the'//
     *    ' following directories: ',  1)
         do i = 1, MaxMediaDir
            if(MediaDir(i) .ne. ' ') then
               call cerrorMsg(MediaDir(i), 1)
            endif
         enddo
         call cerrorMsg('Fatal error, good bye', 0)
      endif

      NoOfMedia = MediaNo       ! memorize the # of media
!      For the moment, don't use next 2 two 
!   If  RangeMkTbl is to be used, next must be called 
!     call epRangeAlloc( NoOfMedia )  ! allocate memory for Range calc.
      
!      do i = 1, NoOfMedia
!         call epRangeMkTbl(Media(i), i)   ! make Range tbl for i-th media
!      enddo

      MediaNo = 1               ! may be not needed ; for safety
!       for setting minimum for Knock-on.
!       see cdedexEleci. in cbeginRun.f
      end
