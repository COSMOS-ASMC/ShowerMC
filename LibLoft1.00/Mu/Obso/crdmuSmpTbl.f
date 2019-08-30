      subroutine cRdmuPrSTbl
      implicit none
!     set  muon  pair ceation cnst
      call cRdmuPrCnst
      end
!     ************************
      subroutine cRdmuBrSTbl
      implicit none
      call cRdmuBrCnst
      end
!     ************************
      subroutine cRdmuNSTbl
      implicit none
#include  "Zcmuint.h"
      call cRdmuNCnst
!        energy dependence
      muNpwtx = 5.982385835778501d-2 
!        energy dependence      
      muNpwdEdx0 =  -1.727630385862580d-2
!         energy dependence
      muNpwdEdxt = 5.289883443973254d-2
      end
!     ************************
      subroutine cRdmuPrCnst
      implicit none
#include  "Zcmuint.h"

      muPrVmin =1.d-3 
      muPrdETX = 0.1d0
      muPrdE = 0.2d0
      muPrEmin = 18.4826551756518d0
      muPrEmax = 73580.7755629467d0
      muPrdU = 1.d-2
      muPrUsize = 101
      muPrEsize = 19
      muPrTXT = 38
      muPrEmax1 = 92632.7081757169
      end
!     ************************
      subroutine cRdmuBrCnst
      implicit none
#include "Zcmuint.h"
      muBrVmin = 1.d-3
      muBrdETX = 0.1d0
      muBrdE = 0.2d0
      muBrEmin = 26.4037931080739d0
      muBrEmax = 66323.3296485571d0
      muBrdU = 1.d-2
      muBrUsize = 101
      muBrEsize = 18
      muBrTXT = 36
      muBrEmax1 = 83496.1250893700d0
      end
!     ************************
      subroutine cRdmuNCnst
      implicit none
#include  "Zcmuint.h"
      muNVmin  = 1.d-3
      muNdETX = 0.1d0
      muNdE  = 0.2d0
      muNEmin = 40.7451833474073d0
      muNEmax = 64576.7637128858d0  
      muNdU = 1.d-2
      muNUsize = 101
      muNEsize = 17
      muNTXT = 34
      muNEmax1 = 81297.3288495795d0
      end
