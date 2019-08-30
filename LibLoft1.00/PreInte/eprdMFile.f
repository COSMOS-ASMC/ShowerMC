!    **********************
      subroutine eprdMFile(name, icon)
      use srimdata
      implicit none
#include "ZepManager.h"      
#include "ZmediaLoft.h"
!#include "ZepTrackp.h"
!#include "ZepTrackv.h"
!#include "Zcnfig.h"


      character*(*) name  ! input. media file name
                          ! Absoft Fortran cannot get correct length
                          ! if *(*) is used and 
                          !  call e...(Det.cmp(i).matter,   )
                          
 
      integer icon        ! outpu. 0, if file exists  and read
                          !        1, if file does not exist

      character*100 msg
      character*150 mdpath
      integer  i, result, klena, jcon

      icon = 0
      do i = 1, MaxMediaDir
!!!!!!!!!!!!
!         write(0,*) 'rmMfile:  i=',i, ' Dir ', trim( MediaDir(i) ),
!     *     ' name ',   trim(name)
!!!!!!!!!!!         
         if( klena(MediaDir(i)) .gt. 0 ) then
            mdpath = ' '
            mdpath = MediaDir(i)(1:klena(MediaDir(i)))
     *            // '/'  //name(1:klena(name))
            call copenf(iowk, mdpath, result)

            if(result .eq. 0) then
!              in the case  of  dir, result=0 and error will happen

               call epReadTab(iowk, Media(MediaNo))
!                 some  correction
!                 knockon cut  energy; this has not been called
!                 in ReadeTab.
!!!!               call epStern(RecoilKEmin, Media(MediaNo))
               call epStern( Media(MediaNo) )   ! v9.154
!     prepare for Urban model; later epResetUrban must be called
!         after fixing Tcut for recoil.
               call epsetUrban(Media(MediaNo), Media(MediaNo)%urb)
!                  cross-section precalculation; not needed
!               call epixsec(Media(MediaNo)) 
               icon =0 
               close(iowk)
!     
!                   read XCOM file
               call epReadXXsec( Media(MediaNo), jcon)
!                  examin if there is SRIM data in the current
!                  media dir v9.09l

               call epSrimChk(iowk, trim(MediaDir(i)),
     *           trim(name),  Media(MediaNo)%srim)

               if(Media(MediaNo)%name .ne. name ) then
!                 file name and name  contained 
!                 there conflict; issue warning ---> change to stop
!                 from v7.25 
!                  if(MsgLevel .ge. 1) then
                     write(msg, *) ' ^^^^^ Warning: Media file:',
     *               name, ' contains media name:',
     *               Media(MediaNo)%name
                     call cerrorMsg(msg, 0)
!                  endif
               endif
               goto 10
            endif
         endif
      enddo
!        name not exists
      icon = 1
 10   continue
!!!!!!
!      write(0,*) '  ret eprdMfile; icon =',icon
!!!!!!!!!!      
      end

      subroutine epwriteMedia
      implicit none
!          print media info
#include "ZmediaLoft.h"
      

      integer:: i
      write(0,*) 'list of meida used'
      do i = 1, NoOfMedia
         write(0, *) Media(i)%name
      enddo
      end       subroutine epwriteMedia
      
