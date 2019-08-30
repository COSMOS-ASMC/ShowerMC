!    ******************************************************************
!    *                                                                *
!    *   cbomegaDcy: Omega (baryon) decay.
!    *                                                                *
!    ******************************************************************
!
       subroutine cbomegaDcy(pj,  a,  np)
       implicit none
#include  "Zptcl.h"
#include  "Zcode.h"

       integer np               !output. no. of ptcls produced
       type(ptcl):: pj         ! input. kaon
       type(ptcl):: a(*)      ! output. produced ptcls
!
!        Omega- ---> Lambda + K- (67.8%)
!               ---> gzai0 + pi- (23.6%)
!               ---> gzai- + pi0 ( 8.6 %)
       integer subcode  ! 
       integer charge   !
       real*8 u

       charge = pj%subcode
       subcode = pj%subcode
       call rndc(u)
       if(u .lt. 0.678d0) then
!           lambda + K
          call cmkptc(klambda, subcode, 0, a(1))
          call cmkptc(kkaon, -subcode, subcode, a(2))
          call c2bdcy(pj, a(1), a(2))
          np=2
       elseif(u .lt. 0.914d0) then
!          gzai0 + pi
          call cmkptc(kgzai, subcode, 0, a(1))
          call cmkptc(kpion, -subcode, subcode, a(2))
          call c2bdcy(pj, a(1), a(2))
          np=2
       else
!          gzai- + pi0 ( 8.6 %)
          call cmkptc(kgzai, subcode, subcode, a(1))
          call cmkptc(kpion, 0, 0, a(2))
          call c2bdcy(pj, a(1), a(2))
          np=2
       endif
      end

