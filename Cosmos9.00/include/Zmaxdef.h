#define MAX_NO_OF_SITES 50
#define MAX_NO_OF_AS_SITES  50
#define MAX_SEGMENTS 60
#define MAX_NO_OF_COMPS 8
#define MAX_AZM 37
#define MAX_ZEN 21
#define MAX_STACK 50000
#define MAX_PTCL  35000
#define MAX_USERHOOKC 5
#define MAX_USERHOOKI 10
#define MAX_USERHOOKR 10

#undef SHOWTHIS
#ifdef SHOWTHIS
/*
c    The user may change the definition of the max
c    size of arrays by the above variables. It is better not 
c    to touch the parameter definition directly such 
c    as in Zobs.h 
c   
c    MAX_NO_OF_SITES:  Max number of observation sites for each ptcl
c    MAX_NO_OF_AS_SITES:  Max number of observation sites for hybrid A.S
c    MAX_SEGMENTS      :  Max number of segments when expressing primary
c                         energy spectrum by a segments of lines for one 
c                         type of primary.
c    MAX_NO_OF_COMPS 8:   Max number of different primary components
c                         definable at a time.
c    MAX_AZM 37:          Max number of azimuthal angle bins for a one rigidity
c                         cut table can contain.
c    MAX_ZEN 21:          Max number of zenith angle bins for a one rigidity 
c                         cut table can contain.
c    MAX_STACK:           Max stack size for particle tracking. If there is
c                         no multi-particle produciton, this may be as small
c                         as ~30, but if multipicity of n must be treated,
c                         this must be as large as a few times * n.
c    MAX_PTCL 8000:       Max number of partilces at a particle production.
c
c    MAX_USERHOOKC/I/R:      size of UserHookc, UserHooki, UserHookr
*/
#endif



