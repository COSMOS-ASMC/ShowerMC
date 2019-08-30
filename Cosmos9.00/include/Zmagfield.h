/*
#include "Zunionmap.h"
c
       structure /magfield/
c          Note that position vector where the magnetic field is given
c          is not included here.
c          unit of field strength is  in T (1 gauss = 10**-4 T)
#ifdef  UNIONMAP
         union
              map
#endif
                  real*8 x,       ! in earth_center coordinate
     *                   y,       !
     *                   z        !

#ifdef UNIONMAP
              endmap              ! 'hva'.  this is for geomag.
              map
                  real*8 n,      ! north com.
     *                   e,      ! east comp.
     *                   d       ! down com.
              endmap             ! 'ned'   this is for geomag
              map
                  real*8 h,       ! horizontal comp.
     *                   v,       ! vertical comp.(down is +)
     *                   a        ! deflection angle (deg. east is +)

              endmap 
          endunion     
#endif
c  
          character*4 sys  ! which system. 'xyz',  'ned',  'hva'
        end structure

*/



struct magfield {
  double x;
  double y; 
  double z;
  char sys[4];
  DUMMYCHAR
};


