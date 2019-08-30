!    ******************************************************************
!    *                                                                *
!    *   clambdacDcy: Lamdac decay : mainly  for muon 
!    *                                                                *
!    ******************************************************************
!
       subroutine clambdacDcy(pj,  a,  np)
       implicit none
#include  "Zptcl.h"
#include  "Zcode.h"

       integer np               !output. no. of ptcls produced
       type(ptcl):: pj         ! input. kaon
       type(ptcl):: a(*)      ! output. produced ptcls
       real*8 u, w
       integer icon
!           lambdac --> lamuda +  mu+ +  neumu  (2 %)
!
       integer:: i

        call rndc(u)
        if(u .lt. .02) then
           if(pj%subcode .eq. regptcl) then
              call cmkptc(klambda, regptcl, 0, a(1))
              call cmkptc(kmuon,   antip,  1,  a(2))
              call cmkptc(kneumu,  regptcl, 0, a(3))
           else
              call cmkptc(klambda, regptcl, 0, a(1))
              call cmkptc(kmuon,   regptcl, -1, a(2))
              call cmkptc(kneumu,  antip, 0, a(3))
           endif
           call cnbdcy(3, pj%mass, a, 0, w, icon)
           np=3
        else
!         many modes with small branching ratios
!         and  no major decay mode but
!           p+anything or n+anything is about 50 %,50%
!        so we use Lambda0 + anything for simplicity
!        and for anything we take pi+,pi0
           call cmkptc(klambda, regptcl, 0, a(1))
           call cmkptc(kpion,   regptcl, 1, a(2))
           call cmkptc(kpion,   regptcl, 0, a(3))
           call cnbdcy(3, pj%mass, a, 0, w, icon)
           np=3
        endif
        do   i=1, np
           call cibst1(i, pj, a(i), a(i))
        enddo
      end
