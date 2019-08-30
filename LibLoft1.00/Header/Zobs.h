#include "Zmaxdef.h"
!       Zobs.h     header file for observation sites definition
!
           integer maxNoOfSites, maxNoOfASSites, horizontal,
     *     perpendicular, notUsed, spherical
       parameter (
#ifdef MAX_NO_OF_SITES
     *    maxNoOfSites = MAX_NO_OF_SITES,
#else		  
     *    maxNoOfSites = 15,
#endif
#ifdef MAX_NO_OF_AS_SITES
     *    maxNoOfASSites=MAX_NO_OF_AS_SITES,
#else
     *    maxNoOfASSites=20,
#endif
     *    notUsed = 0,           ! detector plane is not used
     *    horizontal = 1,        ! detector is horizontal
     *    perpendicular = 2,      ! detector is pependicular to 1ry.
     *    spherical = 3          ! detector is cocentric sphere as the earth
     *		      )

