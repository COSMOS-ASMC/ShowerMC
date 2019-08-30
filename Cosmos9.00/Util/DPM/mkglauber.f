#include "ZcosmosBD.h"
      PROGRAM DPMJET
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "Zmanagerp.h"
#include "Zmanager.h"
#include  "Zevhnp.h"      
! event flag
      integer  NEVENT, ICASCA
      COMMON /DTEVNO/ NEVENT,ICASCA
!-----------------------------------------------------------------------
!       DpmFile is the basic input file for dpmjet.
!       Here, we have to use NO "satart" control word
!       and only make Glauber data.
!
      read(*,'(a)') DpmFile
      write(ErrorOut,*) DpmFile

      IntModel ='"dpmjet3"'

!  
!      CALL DT_DTUINI(NEVTS,EPN,NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,IEMU)
      call cintModels('cosmos')
      END




