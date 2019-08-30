#include "ZepicsBD.h"
#include "ZcosmosBD.h"
      use modsubdTree
      implicit none
!      
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
#include  "ZepManager.h"
      character*100  dsn1
      MediaDir(1) = '$EPICSTOP/Data/Media'
      read(*,'(a)') dsn1
      read(*,*) option
      write(0,*) ' reading config'
      call eprcnf(dsn1)
      call epsubdTree(6)

      end

