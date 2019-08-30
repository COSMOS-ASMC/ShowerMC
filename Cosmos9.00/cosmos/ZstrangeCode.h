!  this is to treat the stange code generated in Lund Fritiof.
!   If you find that you can get reliable events by discarding
!   events with such strange code, then you may define
!   the machine definition here. 
! 
!     If you want to discard the strange code events but
!     want to be notified,  define below.
#undef DEBUG_STRANGECODE

#ifdef DECALPHA
#define DEBUG_STRANGECODE
#endif

#ifdef Solaris
#define DEBUG_STRANGECODE
#endif
!  if you want discard the strage code events without
!  any message, define below. If  DEBUG_STRANGECODE is
!  defined, you should not give the below.
#undef DISCARD_STRANGECODE
#ifdef SOME_MACHINE
#define DISCARD_STRANGECODE
#endif

!  if none of above two is defined, messages will be issued and
!  execution will be halted.
!
!  you should not touch below.
! ----------------------------------------
!
#define NOTIFY_STRANGECODE

#ifdef DISCARD_STRANGECODE
#undef NOTIFY_STRANGECODE
#endif

#undef USE_STRANGEFLAG
#ifdef DEBUG_STRANGECODE
#define USE_STRANGEFLAG
#endif

#ifdef DISCARD_STRANGECODE
#define USE_STRANGEFLAG
#endif


