/*
c          location of a ptcl 
c       Zcoord.h must be  preceeded
c
	structure  /position/
	   record  /coord/ xyz   ! in xyz
	   real*8  radiallen    ! in m . radial length
	   real*8  depth       ! in kg/m2   depth.
	   real*8  height      ! in m.  vertical height(from sea level
           real*8  colheight   ! in m.  //  where the  
c                           latest nuclear collision took place.
c                           (iniitial value is very large value).
	end  structure
*/

struct position {
  struct coord xyz;
  double radiallen;
  double depth;
  double height; 
  double colheight;
};

