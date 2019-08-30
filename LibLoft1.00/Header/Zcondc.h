#ifndef Zcondc_
#define  Zcondc_
!    Description about MYEFIELD is obso; Now it should be
!         defined always. (Aug. 2019)

!      next is probably used only for Cosmos but
!      the Zcondc.h is placed in LibLoft/Header
!     (not harmfule).

!         The simple cmyEfield.f is now in $COSMOSTOP/cosmos
!         and referred from every chook.f. It might be
!         used actually depending on HowEfield.
!         In such a case, the user would need to copy the
!         subroutine to the applicaiton area to make it more
!         complex. Then, the user must 
!  define MYEFIELD  if Electric field is to be supplied 
!  by the user using cmyEfield.f of which template is
!  in UserHook/.  The user may copy it to the users
!  application area, modify it and may add  cmyEfield.o in
!  the chook.mk like:
!      objs =  chook.o cmyEfield.o
!  Also the user must give a value of
!  >1 to the 'HowEfield' parameter  in the namelist ($HPARAM).
!  Note simple Electric field can be specified
!  without using this but by giving HowEfield=1 in the
!  namelist parameter.  Then simple electric
!  field can be specified (together with other parameters)
! #undef MYEFIELD 
#define MYEFIELD
!   dpmjet cannot be used on NEXTSTEP, so
!   you have to make the next 0. 
#ifdef NEXT486
#define  USEDPMJET    0
#else
#define  USEDPMJET    1
#endif
#undef  MATHLOUSY
!   if parameter statement does not permit to use math such as 
!      parameter::pi=asin(1.d0)*2, define MATHLOUSY
#if defined  (Solaris) || (jaxa) || (jaxaflat)
#define MATHLOUSY
#endif

!            make DEBUG > 0 depending on the debug purpose. 

#define  DEBUG 0
!

!     if you want to put a lable on each particle to identify that
!     the one and the same particle crosses a given observation
!     plane more than once, make this 1 or 2.  Then the same particle
!     will have the  same label number in track record.
!     ( aTrack.label ).  If this is 0, aTrack.lable record dose not
!     exists. 
!     If 1; after any interaction (except for continuous energy
!     loss by dE/dx and deflection by B or scattering), label is
!     changed.
!     If 2: For knockon and Bremstrahlung, the survival particle
!     will have the same label. In the case of Moller scattring
!     higher enregy electrons are regarded as the survival one.
!
#define  LABELING    0
!     if you want to have a detailed info. for particle tracking
!     make the below >=1.  The user observation routine is called
!     with the following id  on the following  conditions:
!              chookobs(a, id)
!     1)  if it is >=1,  a particle is going to interact at a point given in
!         the track information, id=4
!     2)  if it is >=1,  a particle is going to die, id=5
!     3)  if it is >=2,  a particle is being discarded due to the large
!          angle (cos(angle relative to the parent) > BackAngLimit). id=6
!     4)  if it is >=3,  a particle makes a step. id=7
!        
#define  DETAILED_TRACKING 0
#endif  /* Zcondc%h */
