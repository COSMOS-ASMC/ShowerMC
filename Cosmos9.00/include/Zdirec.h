/*
      structure /direc/
          record /coord/ w
	  real*8  coszenith   ! cos of the zenith angle.  
c               it is defined as follows:
c                   Let's assume w and position are given
c                   in xyz sytem.
c                  
c                   coszenith = -( x*w.x + y * w.y + z * w.z )/
c                                (length of (x,y,z)) 
c                   This should be computed whenever w is
c                   updated.
      end structure

*/
struct direc {
  struct coord w;
  double coszenith;
};




