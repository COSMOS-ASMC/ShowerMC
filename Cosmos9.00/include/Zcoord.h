/*

#include "Zunionmap.h"

c        sys="xyz":  origin is center of the Earth.
c              x:  directed to longitude 0, latitude 0
c              y:  directed to longitude 90 deg, latitude 0
c              z:  center to the North pole
c  ****************************************************************
c  *       During the paticle tracking, this system is used.      *
c  ****************************************************************
c


c
      structure /coord/
#ifdef  UNIONMAP
          union
              map
#endif
                  real*8 r(3)
#ifdef UNIONMAP
              endmap
              map
                  real*8 x, y, z  !  x,y,z in m
              endmap          ! 'xyz'
              map 
                  real*8 lat,    ! latitude in deg.  + is to the north.
     *                  long,    ! longitude in deg. + is to the east.
     *                     h     ! height in m       
              endmap      !  'llh'
              map
                  real*8 theta,      ! polar angle
     *                   phi,        ! azimuthal angle
     *                   radius      ! radial distance
              endmap         !   'sph'
          endunion
#endif
          character*4 sys  ! which system. 'xyz', 'llh', 'sph'
      end structure
*/

struct coord {
  double r[3];
  char sys[4]; 
  DUMMYCHAR  
};

