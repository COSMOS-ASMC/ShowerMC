! #ifndef Zpos_
! #define Zpos_
!          location of a ptcl 
!       Zcoord.h must be  preceeded
!
         type  position
           sequence
           type(coord):: xyz   ! in xyz
            real*8  radiallen    ! in m . radial length
            real*8  depth       ! in kg/m2   depth.
            real*8  height      ! in m.  vertical height(from sea level
           real*8  colheight   ! in m.  //  where the  
!                           latest nuclear collision took place.
!                           (iniitial value is very large value).
         end  type position
!  #endif  /* Zpos.h */
