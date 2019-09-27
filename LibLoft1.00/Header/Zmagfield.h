! #ifndef  Zmagfield_
! #define  Zmagfield_
#include "Zunionmap.h"
!
       type magfield
         sequence
!          Note that position vector where the magnetic field is given
!          is not included here.
!          unit of field strength is  in T (1 gauss = 10**-4 T)
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
	      map
		real(8):: r(3)
	      endmap
          endunion     
#endif
!  
          character*4 sys  ! which system. 'xyz',  'ned',  'hva'
        end type magfield
!  #endif  /* Zmagfield.h */
