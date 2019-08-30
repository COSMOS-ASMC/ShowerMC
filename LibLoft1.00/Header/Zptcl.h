! #pragma once
!       both of above and below NG if a file contains many subroutines
! #ifndef Zptcl_
! #define Zptcl_
#include "Zunionmap.h"
!    structure defining a particle at production
!         Basic idea of what is to be contained in 
!         the particle structue is that
!        1) dynamical ones should be included
!        2) those derivable from the particle code
!           is not included 
!     ******************************************************
      type fmom     ! 4 momentum
	sequence
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
!                         pt before pz is set
                  real*8  dummy1, dummy2, pt, rap
              endmap
!                                tm: transverse mass
              map
                  real*8  dummy3, dummy4, tm
              endmap
          endunion   
#endif
      end type fmom
!     \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
!       Important note:   Bug in sun fortran
!           If we define, say,
!                 record /fmom/ p1
!           and set
!                 p1.e = some value (or p1.p(4)= ...)
!           where some value is a constant or arithmetic
!           expression which results in a value > 1.d37
!           then overflow message comes out on SUN fortran
!           although the result is correct.
!           Setting the same into, say, p1.px does not
!           cause such. (as of 1993/08/14)
!     \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
!     ******************************************************
      type ptcl       ! particle at production
        sequence
!                   4 momentum. 

      type(fmom):: fm 
!
          real*8 mass
          integer*2 code, subcode  
          integer*2 charge
!       code: ptcl code
!    subcode:used mainly to identify paticle/antiparticle
!            if the difference is important.
!            To set particle, "ptcl" is used.
!                   anti-partilce, 'antip" is used for particles
!           For particles of which partilce/antiparticle nature
!            can be judded by its code and charge, the user 
!            need not specify it when using cmkptc subroutine.
!            give 0.
!            subcode for gamma ray may be used to identify
!            brems gamma and direct gamma by kdiretg, kcasg
      end type ptcl
!     ******************************************************
!  #endif  /* Zptcl.h */
