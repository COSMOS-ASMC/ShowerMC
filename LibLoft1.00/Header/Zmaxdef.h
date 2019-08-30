#define MAX_NO_OF_SITES 50
#define MAX_NO_OF_AS_SITES  50
#define MAX_SEGMENTS 60
#define MAX_NO_OF_COMPS 11
#define MAX_AZM 37
#define MAX_ZEN 21
#define MAX_STACK 50000
#define MAX_PTCL  35000
#define MAX_USERHOOKC 5
#define MAX_USERHOOKI 10
#define MAX_USERHOOKR 10
!   make next as large as 1~2 milion for practical case
!  (for general MPI application)
#define MAX_SKELPTCL  1
#undef SHOWTHIS
#ifdef SHOWTHIS
!    The user may change the definition of the max
!    size of arrays by the above variables. It is better not 
!    to touch the parameter definition directly such 
!    as in Zobs.h 
!   
!    MAX_NO_OF_SITES:  Max number of observation sites for each ptcl
!    MAX_NO_OF_AS_SITES:  Max number of observation sites for hybrid A.S
!    MAX_SEGMENTS      :  Max number of segments when expressing primary
!                         energy spectrum by a segments of lines for one 
!                         type of primary.
!    MAX_NO_OF_COMPS 8:   Max number of different primary components
!                         definable at a time.
!    MAX_AZM 37:          Max number of azimuthal angle bins for a one rigidity
!                         cut table can contain.
!    MAX_ZEN 21:          Max number of zenith angle bins for a one rigidity 
!                         cut table can contain.
!    MAX_STACK:           Max stack size for particle tracking. If there is
!                         no multi-particle produciton, this may be as small
!                         as ~30, but if multipicity of n must be treated,
!                         this must be as large as a few times * n.
!    MAX_PTCL 8000:       Max number of partilces at a particle production.
!
!    MAX_USERHOOKC/I/R:      size of UserHookc, UserHooki, UserHookr
#endif



