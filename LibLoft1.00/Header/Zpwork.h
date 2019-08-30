#include "Zmaxdef.h"
      !  #include "Zptcl.h" must preceed this one.
!       info of  produced  particles
                    ! max # of ptcls producable in coll.
      integer,parameter::  MaxPtcl = MAX_PTCL
      type(ptcl):: Pwork(MaxPtcl) ! working array to store ptcls.
      integer:: Nproduced         ! no. of ptcls produced and stored in Pwork.
      integer:: Nstacked          ! no. of ptcls stacked. If ThinSampling=F,
! same as Nproduced.  Nstacked <= Nproduced
      common /Zpwork/ Pwork, Nproduced, Nstacked 


      
