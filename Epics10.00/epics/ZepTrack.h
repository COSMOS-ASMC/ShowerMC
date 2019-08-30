!    structure used when tracking a particle in epics
!    *************************
#ifdef   PTCL
#else
#include  "Zptcl.h"
#endif
#include  "ZepPos.h"
#include  "ZepDirec.h"
!     ---------------------
       type epTrack       ! full particle attributes in epics
       sequence

       type(ptcl)::  p    ! basic ptcl attributes
!               position and time
       type(epPos)::  pos
          real*8 t           ! time in length/beta (m)
       type(epDirec)::  w  ! direction cos.
          real*4  wgt      ! weight for thin sampling or scintillation light
                           ! from a cell
	  real*4  pol      ! polarization of photons (light)
	  real*8  wl       ! wave length for light or path length for
                           !  charged particle for Cerenkov light emmision
			   ! (nm or cm).
	  integer cn    ! current componet number where the particle is.
          integer inciflag  ! for incident particle, 1 is set
                     ! decendents are 0
	  real*8  user  ! user use.
       end type epTrack


