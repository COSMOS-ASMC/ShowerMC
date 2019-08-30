!  #ifndef Zcoord_
!  #define Zcoord_
#include "Zunionmap.h"
!        sys="xyz":  origin is center of the Earth.
!              x:  directed to longitude 0, latitude 0
!              y:  directed to longitude 90 deg, latitude 0
!              z:  center to the North pole
!  ****************************************************************
!  *       During the paticle tracking, this system is used.      *
!  ****************************************************************
!


!
      type coord
        sequence
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
      end type coord
!  #endif /* Zcoord.h */
