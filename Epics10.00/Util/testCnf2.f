#include "ZepicsBD.h"
!           This is to print Emin for each compoent
      implicit none
!       to print configuration with Emin 
!       unlike tesntCnf1, you need epicsfile path
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
#include  "ZepManager.h"
      character*100  dsn1, dsn2 
      integer length, i, j

      read(*, '(a)') dsn1
      read(*, '(a)') dsn2
      write(0,*) trim(dsn1)
      write(0,*) trim(dsn2)

      call epprmr(dsn2)
      call epcmp1  ! compute some (EminGsave etc)
      write(0,*) 'reading config'

      call eprcnf(dsn1)
      call epOutCnfWithEmin(0, 1000, 6)

      end





