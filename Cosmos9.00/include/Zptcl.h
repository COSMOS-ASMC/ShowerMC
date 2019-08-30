/*
#include "Zunionmap.h"
c    structure defining a particle at production
c         Basic idea of what is to be contained in 
c         the particle structue is that
c        1) dynamical ones should be included
c        2) those derivable from the particle code
c           is not included 
c     ******************************************************
      structure /fmom/     ! 4 momentum
#ifdef UNIONMAP
          union
              map
#endif
                  real*8 p(4)
#ifdef UNIONMAP
              endmap    
              map
                  real*8 px, py, pz, e
              endmap
              map
                  real*8  x, y, z, t
              endmap
              map
c                         pt before pz is set
                  real*8  dummy1, dummy2, pt, rap
              endmap
c                                tm: transverse mass
              map
                  real*8  dummy3, dummy4, tm
              endmap
          endunion   
#endif
      end structure   
c     \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
c       Important note:   Bug in sun fortran
c           If we define, say,
c                 record /fmom/ p1
c           and set
c                 p1.e = some value (or p1.p(4)= ...)
c           where some value is a constant or arithmetic
c           expression which results in a value > 1.d37
c           then overflow message comes out on SUN fortran
c           although the result is correct.
c           Setting the same into, say, p1.px does not
c           cause such. (as of 1993/08/14)
c     \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
c     ******************************************************
      structure /ptcl/       ! particle at production
c                   4 momentum. 
          record /fmom/ fm 
c
          real*8 mass
          integer*2 code, subcode  
          integer*2 charge
c       code: ptcl code
c    subcode:used mainly to identify paticle/antiparticle
c            if the difference is important.
c            To set particle, "ptcl" is used.
c                   anti-partilce, 'antip" is used for particles
c           For particles of which partilce/antiparticle nature
c            can be judded by its code and charge, the user 
c            need not specify it when using cmkptc subroutine.
c            give 0.
c            subcode for gamma ray may be used to identify
c            brems gamma and direct gamma by kdiretg, kcasg
      end structure   
c     ******************************************************

*/

struct fmom {
  double p[4];
};
struct ptcl {
  struct fmom fm; 
  double  mass;
  short int code;
  short int subcode;  
  short int charge;
};




