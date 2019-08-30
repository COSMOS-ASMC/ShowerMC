!  #ifndef Zdirec_
!  #define Zdirec_
      type direc
        sequence
          type (coord) w
           real*8  coszenith   ! cos of the zenith angle.  
!               it is defined as follows:
!                   Let's assume w and position are given
!                   in xyz sytem.
!                  
!                   coszenith = -( x*w.x + y * w.y + z * w.z )/
!                                (length of (x,y,z)) 
!                   This should be computed whenever w is
!                   updated.
      end type direc
!  #endif /* Zdirec.h */
