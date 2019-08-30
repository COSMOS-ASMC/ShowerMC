      subroutine readMedia
      implicit none
#include "ZepTrackv.h"
#include "Zmass.h"
!      record /epmedia/media
      integer  io,  icon
      character*120 file
      character(len=8):: name
      io = 10
      write(0,*)
     *  'Enter base Media file name in $EPICSTOP/Data/BaseM/'
      write(0,*) ' such as PWO'
      read(*, *)  name
      call cerrorMsg(name,  1)
      file = '$EPICSTOP/Data/BaseM/'//trim(name)
      call copenf(io, file, icon)
      call epReadMTbl(io, media)
      close(io)
      call epGetEffZA(media)
      write(0,*) ' exiting readMedia'
      NoOfMedia =1
      end       subroutine readMedia




