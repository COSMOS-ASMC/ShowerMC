!    ******************************************************************
!    *                                                                *
!    *   clambdaDcy: Lamda0 decay
!    *                                                                *
!    ******************************************************************
!
       subroutine clambdaDcy(pj,  a,  np)
       implicit none
#include  "Zptcl.h"
#include  "Zcode.h"

       integer np               !output. no. of ptcls produced
       type(ptcl):: pj         ! input. kaon
       type(ptcl):: a(*)      ! output. produced ptcls
       real*8 u
       integer subcode, chargepi, chargen
!           lambda --> p pi- (64.2%)
!                 ---> n pi0 (35.8%)
!
        call rndc(u)
        if(u .lt. .642) then
           chargepi = pj%subcode
           chargen = -pj%subcode
           subcode = pj%subcode
        else  
           chargepi = 0
           chargen = 0
           subcode = pj%subcode
        endif
        call cmkptc(kpion, regptcl, chargepi, a(1))
        call cmkptc(knuc, subcode, chargen, a(2))
        call c2bdcy(pj, a(1), a(2))
        np=2
      end
