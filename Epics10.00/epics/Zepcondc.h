#ifndef EPCONDC
!       angle in the config file is in degree; (e.g  sin(90.)  is 1.0)
#define DEGANGLE
!       if f90 is available, USEFORMULA will permit to use
!        mathmatical formula in config file
#undef  USEFORMULA
!

!       change the next #undef INTINFO to
!       #define INTINFO, then the system will inform the user
!       the interaction information for each interaction by calling epUI.
!       To make this works
!         1)   The user must copy epUI.f in UserHook to the user's
!              hook directory (say in UserHook/myAplication)
!         2)   change the following part in ephook.mk.
!              objs =  ephook.o eppos2B.o 
!           as
!              objs =  ephook.o eppos2B.o epUI.o
!         3)  To keep INTINFO defined, do the same thing 
!            for all your applications. 

#undef  INTINFO
#define  INTINFO
#endif
#define EPCONDC




