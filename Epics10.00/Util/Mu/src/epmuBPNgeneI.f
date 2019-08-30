      subroutine epmuBPNgeneI(io, media)
      implicit none

#include "Zmedia.h"
       type(epmedia)::  media  ! ouput. media data is read

!          initialization before using 
!                                      
!         Besides this,
!         you have to give Emu 
!         in ZmuBPgene.h
!
      integer io  ! input. logical  device number for
!                          basic media file which must have been
!                          opened with this number.
!


      call epReadMTbl(io, media)
      call epGetEffZA(media)
      call epSetSTblCns(media, media%cnst)
      end
