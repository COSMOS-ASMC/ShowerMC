#define  UNIONMAP

!#    for gfortran UNIONMAP must be disabled.  
#if defined  (IBMAIX) || (LinuxGfort) || (MacGfort)
#undef UNIONMAP
#endif

