!    ******************************************************************
!    *                                                                *
!    *   cDeltaDecay: Delta (3,3 reso) decay
!    *                                                                *
!    ******************************************************************
!
       subroutine cDeltaDecay(pj,  a,  np)
       implicit none
#include  "Zptcl.h"
#include  "Zcode.h"

       type(ptcl),intent(in):: pj       !  Delta
       type(ptcl),intent(inout):: a(*)    !  produced ptcls
       integer,intent(out):: np   ! no. of ptcls produced =2

       real(8):: u
       integer subcode, chargePi, chargeN
!     Delta+  --> p pi0    or n pi+
!     Delta-  --> n pi-    
!     Delta0  --> n pi0   or  p pi-

!        Ng mode neglected       
       call rndc(u)

       if( pj%charge > 0 ) then
!              Delta+  --> p pi0    or n pi+
          if( u < 0.5 ) then
                !   p pi0 
             chargeN = 1
             chargePi = 0
          else
             !   n pi+
             chargeN = 0
             chargePi = 1
          endif
          call cmkptc(kpion, regptcl, chargePi, a(1))
          call cmkptc(knuc,  regptcl, chargeN, a(2))
       elseif( pj%charge == 0 ) then
          !     Delta0  --> n pi0   or  p pi-
          if( u < 0.5 ) then
                !   n pi0 
             chargeN = 0
             chargePi = 0
             subcode = regptcl
          else
             !   p pi-
             chargeN = 1
             chargePi = -1
             subcode = antip
          endif
          call cmkptc(kpion, subcode, chargePi, a(1))
          call cmkptc(knuc,  regptcl, chargeN, a(2))
       else
!     Delta-  --> n pi-              
          call cmkptc(kpion, antip, -1, a(1))
          call cmkptc(knuc,  regptcl, 0, a(2))
       endif
       

       call c2bdcy(pj, a(1), a(2))
       np=2
       end
