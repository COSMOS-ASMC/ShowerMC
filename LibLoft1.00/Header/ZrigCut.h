#include "Zmaxdef.h"

!         rigidity cut related commons
!     *******************************
      integer maxAzm, maxZen, maxDim, maxRig, maxDim2
      parameter (
#ifdef MAX_AZM
     * maxAzm = MAX_AZM,
#else
     * maxAzm = 37, 
#endif
#ifdef MAX_ZEN
     * maxZen=MAX_ZEN,
#else
     * maxZen=21,
#endif
     * maxRig = 40, 		 
     * maxDim=maxAzm* maxZen, maxDim2 = 12*20*maxRig)

      real*4  RigTbl2(maxDim2)  !!!!!!!  *4
      real*8  RigCutTbl(maxDim), DAzm, DZen, Longi, Latit,
     *        LogDRig, DRig
      integer	AzmSize, ZenSize, RigSize, Rdatafmt
      real*8  MagDec, MinAzm, MinZen, MinRig, ZenMax
      common /Zrigcut/ RigCutTbl, Longi, Latit, MagDec,
     *  DAzm, DZen, MinAzm, MinZen, LogDRig, MinRig,
     *  ZenMax, DRig,  RigTbl2,
     *  AzmSize, ZenSize, RigSize

      character*16 Place
      character*3  AzmValue, ZenValue
      common /Zrigcutc/ Place, AzmValue, ZenValue, Rdatafmt
!     *****************************


