#include "Zmedia.h"
     
      integer,parameter:: Maxmedia=MAX_MEDIA
      type(epmedia):: Media(Maxmedia)

       !     current particle media no = Cn2media(cn)
      integer:: MediaNo
      integer:: NoOfMedia   !  No of media

        ! max number of directories where
      integer,parameter:: MaxMediaDir=MAX_MEDIA_DIR
 
      character(128)::MediaDir(MaxMediaDir)
      common /ZmediaLoft/ Media, MediaNo,NoOfMedia
      common /ZmediaLoftc/ MediaDir
