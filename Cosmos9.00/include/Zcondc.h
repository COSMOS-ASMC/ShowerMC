#ifndef ZCONDC
#define  ZCONDC

/*
c   dpmjet cannot be used on NEXTSTEP, so
c   you have to make the next 0. 
*/

#ifdef NEXT486
#define  USEDPMJET    0
#else
#define  USEDPMJET    1
#endif

// c            make DEBUG > 0 depending on the debug purpose. 

#define  DEBUG 0
/*
c
c   choose:    Old atmosphere or new segmented atmosphere
c            define 
c               old atmosphere --> 0
c           or  new with c-spline
c               new atmosphere --> 1
c           or  new with linear interp.
c               new atmosphere --> 2
*/

#define  ATMOSPHERE  2
/*
c     if you want to put a lable on each particle to identify that
c     the one and the same particle crosses a given observation
c     plane more than once, make this 1 or 2.  Then the same particle
c     will have the  same label number in track record.
c     ( aTrack.label ).  If this is 0, aTrack.lable record dose not
c     exists. 
c     If 1; after any interaction (except for continuous energy
c     loss by dE/dx and deflection by B or scattering), label is
c     changed.
c     If 2: For knockon and Bremstrahlung, the survival particle
c     will have the same label. In the case of Moller scattring
c     higher enregy electrons are regarded as the survival one.
c
*/
#define  LABELING    0
/*
c     if you want to have a detailed info. for particle tracking
c     make the below >=1.  The user observation routine is called
c     with the following id  on the following  conditions:
c              chookobs(a, id)
c     1)  if it is >=1,  a particle is going to interact at a point given in
c         the track information, id=4
c     2)  if it is >=1,  a particle is going to die, id=5
c     3)  if it is >=2,  a particle is being discarded due to the large
c          angle (cos(angle relative to the parent) > BackAngLimit). id=6
c     4)  if it is >=3,  a particle makes a step. id=7
c        
*/

#define  DETAILED_TRACKING 0
#endif




