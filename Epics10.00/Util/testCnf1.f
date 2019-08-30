!     expand config data.
#include "ZepicsBD.h"
      implicit none
!      
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
#include  "ZepManager.h"
      character*100  dsn1
      MediaDir(1) = '$EPICSTOP/Data/Media'
      write(0,*) ' Enter path to the config file'
      read(*,'(a)') dsn1
      write(0,*) ' reading config'
      call eprcnf(dsn1)
!  from level 1 to 1000 --> all data output-->6 stdout
      call epOutCnf(0, 1000, 6)
      end

