!    ******************************************************************
!    *                                                                *
!    *   cgzaiDecay: gzai and anti gzai decay
!    *                                                                *
!    ******************************************************************
!
       subroutine cgzaiDecay(pj,  a,  np)
       implicit none
#include  "Zptcl.h"
#include  "Zcode.h"

       integer np               !output. no. of ptcls produced
       type(ptcl):: pj         ! input. kaon
       type(ptcl):: a(*)      ! output. produced ptcls
!
       integer subcode  ! of lambda
       integer charge   ! of pi
!           gzai- --> lambda + pi-
!           gzai-Bar --> lambda_bar + pi+
       
       charge = pj%subcode
       subcode = pj%subcode
!     
       call cmkptc(kpion, regptcl, charge, a(1))
       call cmkptc(klambda, subcode, 0, a(2))
       call c2bdcy(pj, a(1), a(2))
       np=2
      end

