#include "Zmaxdef.h"
/*
c       Zobs.h     header file for observation sites definition
c
*/
#ifdef MAX_NO_OF_SITES
const int  maxnoofsites = MAX_NO_OF_SITES;
#else		  
const int  maxnoofsites = 15;
#endif
#ifdef MAX_NO_OF_AS_SITES
const int  maxnoofassites = MAX_NO_OF_AS_SITES;
#else
const int  maxnoofassites = 20;
#endif
const int  notused = 0; // detector plane is not used
const int  horizontal = 1;   // detector is horizontal
const int  perpendicular = 2;  // detector is pependicular to 1ry.
const int  spherical = 3;  // detector is cocentric sphere as the earth
