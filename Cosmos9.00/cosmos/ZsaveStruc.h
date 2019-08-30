!  this is to treat the save statement for the structure construct.
!    if it cannot be recognized by the compiler, but static compile
!    is possible, define BYSTATC
!    if it cannot be recognized totally,  define DONTSAVE
! 
#undef BYSTATIC
#undef DONTSAVE

#ifdef DECALPHA
#define BYSTATIC
#endif

#ifdef CONVEX
#define BYSTATIC
#endif

#define USESAVE
#ifdef BYSTATIC
#undef USESAVE
#endif

#ifdef DONTSAVE
#undef USESAVE
#endif

