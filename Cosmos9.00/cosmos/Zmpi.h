! #include "Zmaxdef.h"  /*  use this   if not used yet */
!  #include "Ztrack.h"  /*  use this   if not used yet */
#include "Zmpibasic.h"
         integer maxLowEptcl, fnoLowE, nFiles, maxFiles
!           actual logical file #  fnoLowE+nFiles
!           it will be used during mpi-skeleton making time
!           so when the program starts to normal run,
!           it will not be used.
          parameter (maxLowEptcl=MAX_SKELPTCL, fnoLowE=15, maxFiles=15)
!         record /track/, allocatable ::  Loft(:)

         integer nLowEptcl
         character*8 skelfleshstat
         character*120 skelworkbase
         integer lengOfBase
          character*124 skelworkfile
          real*8  KEminObsSave  
         integer EndLevelSave
         character*16 GenerateSave  ! this save is non-sense
!             to send/recv track structure
         integer  block(10), byteloc(10), btype(10)

         common /commmpi2/ KEminObsSave, EndLevelSave,
     *   nLowEptcl, nFiles,  lengOfBase, block, byteloc,
     *   btype 
         common /commmpic/ skelfleshstat, GenerateSave, skelworkbase,
     *   skelworkfile      
