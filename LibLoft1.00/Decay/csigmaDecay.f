!    ******************************************************************
!    *                                                                *
!    *   csigmaDecay: sigma +0- and their antiptcl
!    *                                                                *
!    ******************************************************************
!
       subroutine csigmaDecay(pj,  a,  np)
       implicit none
#include  "Zptcl.h"
#include  "Zcode.h"

       integer np               !output. no. of ptcls produced
       type(ptcl):: pj         ! input. kaon
       type(ptcl):: a(*)      ! output. produced ptcls
!
!
      if(pj%charge .eq. 1) then
!           Sigma +
         call csigmaPDcy(pj, a, np)
      elseif(pj%charge .eq. 0) then
!           sigma 0
         call csigma0Dcy(pj, a, np)
      elseif(pj%charge .eq. -1) then
!           sigma -1         
         call csigmaMDcy(pj, a, np)
      endif
      end
      subroutine csigmaPDcy(pj,  a, np)
!
!         sigma+ decay
!         1) Sigma---->p pi0      (51.6%0
!         2)      ---->n pi+       48.4 
      implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
      integer np
      type(ptcl):: pj
      type(ptcl):: a(*)
      integer subcode, charge

      real*8 u
      call rndc(u)
      subcode = pj%subcode
      charge = -pj%subcode
      if(u  .lt. .516) then
!          p+pi0  or conjugae       
         call cmkptc(kpion, 0, 0, a(1))
         call cmkptc(knuc, subcode, charge, a(2))
         call c2bdcy(pj, a(1), a(2))
         np=2
      else
!          n pi+
         call cmkptc(kpion, regptcl, -charge, a(1))
         call cmkptc(knuc,  subcode,  0, a(2))
         call c2bdcy(pj, a(1), a(2))
         np=2
      endif
      end
      subroutine csigma0Dcy(pj, a, np)
!
!         sigma0 decay
!         1) Sigma---->Lamda gamma  100 %
      implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
      integer np
      type(ptcl):: pj
      type(ptcl):: a(*)
      integer subcode

      subcode = pj%subcode
      call cmkptc(klambda, subcode, 0, a(1) )
      call cmkptc(kphoton, 0, 0, a(2) )
      call c2bdcy(pj, a(1), a(2))
      np=2
      end
      subroutine csigmaMDcy(pj,  a, np)
!
!         sigma- decay
!         1) Sigma---->n pi-   100 %
      implicit none
#include  "Zptcl.h"
#include  "Zcode.h"      
      integer np
      type(ptcl):: pj
      type(ptcl):: a(*)
      integer subcode, charge

      subcode = pj%subcode
      charge = pj%subcode
      
!          n+pi-
         call cmkptc(kpion, 0, charge, a(1))
         call cmkptc(knuc, subcode, 0 , a(2))
         call c2bdcy(pj, a(1), a(2))
         np=2
      end
