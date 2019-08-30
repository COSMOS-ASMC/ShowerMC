!!!!!!!!!!!! Now not needed. !!!!!!!!!!!!!!!
!     calling seq.
!        0) after call eprcnf
!        1) call  epiniXsMedia which calles 
!               call cAllocXsecMediaArray( NoOfMedia )  !alloc Media
!               call cAllocXsecElemArray(i, Media(i).NoOfElem) ! Elem
!               call epsetXsMedia(i,  Media(i)) ! set values
!        2) In ccosIntF. instead of epgetxs
!           use cgetXsec
      subroutine epiniXsMedia

      use modXsecMedia, only: cAllocXsMediaArray, cAllocXsElemArray
      implicit none
   ! NoOfMedia, Media,  here  are used
#include "ZepTrackv.h"
      integer::i

      call cAllocXsMediaArray( NoOfMedia ) ! alloc # of media
      do i = 1, NoOfMedia     ! this is in modXsecMedia
         call cAllocXsElemArray(i, Media(i)%NoOfElem) ! alloc # of
                ! elements in each medium
         call epsetXsMedia(i,  Media(i)) ! give some
              ! values not chainging
      enddo
      end      subroutine epiniXsMedia


      subroutine epsetXsMedia(n, media)
!                         to avoid name collisions
      use modXsecMedia, xmedia=>media, xelement=>element
      implicit none
#include "Zmedia.h"
      integer,intent(in)::n ! n-th medium
       type(epmedia)::  media ! input it's ingredient info.

      integer:: i

      xmedia(n)%name = media%name
      xmedia(n)%noOfElem = media%noOfElem

      xmedia(n)%elem(1:media%noOfElem)%A = 
     *         media%elem(1:media%noOfElem)%A
      xmedia(n)%elem(1:media%noOfElem)%Z = 
     *         media%elem(1:media%noOfElem)%Z
      xmedia(n)%elem(1:media%noOfElem)%No = 
     *         media%No(1:media%noOfElem)
      xmedia(n)%mbtoPkgrm = media%mbtoPgrm/10.d0
      end      subroutine epsetXsMedia
