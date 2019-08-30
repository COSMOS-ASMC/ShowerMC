!         this list block data names as external names so that the
!       block data is surely enabled.
#ifdef  ZCONDC
#else
#include  "Zcondc.h"
#endif
      external cblkElemag
      external cblkdedx
      external cblkHeavy 
      external cblkManager 
#if  ATMOSPHERE == 0 
      external cblkStdAtmos
#endif
      external cblkXsec
      external cblkEvhnp
      external cblkIncident
      external cblkObs
      external cblkTracking
      external cblkMuInt
!              For Lund related block common
      external blkhd1, blkhd2, blkhd3, blkhd4
      external blkdc1, blkdc2
      external ludataC
      external luedatC
      external luhdatC
      external cblkSpecial


