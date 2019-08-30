      module modepTPXS
      use modXsecMedia, xmedia=>media, xelement=>element
      use modTPXS
      implicit none
      type(TPXSconst),pointer::TPXSem(:) ! e-
      type(TPXSconst),pointer::TPXSep(:) ! E+
      end   module modepTPXS

      subroutine epmkTPXSconst
      use modepTPXS
      implicit none
#include "ZepTrackv.h"
      integer:: i
      allocate( TPXSem(NoOfMedia) )
      do i = 1, NoOfMedia
         call cMWconstForMedia(Media(i), -1, TPXSem(i))
      enddo

      allocate( TPXSep(NoOfMedia) )
      do i = 1, NoOfMedia
         call cMWconstForMedia(Media(i), 1, TPXSep(i))
      endif
      end   subroutine epmkTPXSconst
