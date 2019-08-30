#include "ZepMaxdef.h"
        integer NPreDefName   ! number of predefined structures
        parameter ( NPreDefName = 5 )
  
        integer MaxNewStruc   ! max number of new structures constructable
        parameter ( MaxNewStruc = MAX_NEW_STRUC )

        integer nattrib(NPreDefName + MaxNewStruc) 
                                     ! each structure has nattrib(i)
                                     ! attributes for a volume
                                     ! box--> i=1, cyl--> i=2
                                     ! for new-i, it is determined at run time


        character*8  PreDefName(NPreDefName)

        
        common /Zcnfig2/  nattrib

        common /Zcnfic2/ PreDefName


