      program showParam
!         to read namelist and write them
!         This  may be used to get all the parameter list.
!         
      implicit none
#include "Zmanagerp.h"
!          testing creadParm
      call cerrorMsg(
     * 'You can get Defaults by putting &Param &end', 1)
      call creadParam(5)
      call cwriteParam(ErrorOut, 1)
      end
#include "BlockData/cblkGene.h"
