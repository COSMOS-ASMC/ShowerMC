      subroutine mediaCount(cmp, ncmp, media, nmedia)
!       count number of different media in 'cmp' component 
      implicit none
#include "ZepDraw.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      integer ncmp                ! input.  number of components in cmp
       type(Component)::  cmp(ncmp)   ! input.  
      integer nmedia    ! output.  No. of different media used.
      character*8 media(*)    ! output.  media(i) (i=1, nmedia) is media id.


      integer i, j
      character*8 mat
!
!      if(NumComp .gt. 0) then
      NoOfCompsEachMedia(:)=0
      if(ncmp .gt. 0) then
         nmedia = 1
         compnToMn(1) = 1
         media(1) = cmp(1)%matter
         NoOfCompsEachMedia(1) = 1
!         do i = 2, NumComp
         do i = 2, ncmp
            mat = cmp(i)%matter
            do j = 1, nmedia
               if(mat .eq. media(j)) then
                  compnToMn(i) =  j
                  NoOfCompsEachMedia(j) = 
     *               NoOfCompsEachMedia(j) + 1
                  goto 100
               endif
            enddo
            if( nmedia == maxmat ) then
               write(0,*) ' error in mediaCount '
               write(0,*) ' increase maxmat in ZepDraw%h'
               stop
            endif
            nmedia = nmedia + 1
            compnToMn(i) = nmedia
            media(nmedia) = mat
            NoOfCompsEachMedia(nmedia) = 1
 100        continue
         enddo
      else
         nmedia = 0
      endif
      end
