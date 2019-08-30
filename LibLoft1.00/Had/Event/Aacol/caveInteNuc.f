!     *****************************************************************
!     *                                                               *
!     *  caveInteNuc:  get average no. of interacting nucleons
!     *  when a heavy ptcl interacts with air target   
!     *                                                               *
!     *****************************************************************
!  Note: this interacting nucleons do not include successive 
!        interactions; the number is smaller by <Nsuc> factor
!        than <N> = A sigma(pA)/sigma(A'A)
!
!
        subroutine caveInteNuc(pj, tgtMassN,  avn )
        implicit none

#include  "Zcode.h"
#include  "Zptcl.h"
#include  "Zheavyp.h"

        type(ptcl):: pj   !  projectile heavy
        real*8  avn  !  average no. of interacting nucleons 
        integer tgtMassN  ! target mass # ( # of nucleons)
        
        real*8  xspp, xspA, xsAA, tgtA, a3, sucave, pne, tgtZ
        integer ihg
!
        tgtA = tgtMassN
        tgtZ = tgtA*0.4   ! almost dummy
        a3 = tgtA**0.333333333
        ihg = Charge2heavyG(pj%charge)
        pne = pj%fm%p(4)/HeavyG2massN(ihg)

!        call cppXsec(pne, xspp)
!        call cxp2xAXsec(tgtA, xspp, xspA)
        call cinelx(pj, tgtA,tgtZ, xspA)
!                         this is 0.4d0*tgtA is rather dummy
        call cAAXsec2(pj, tgtA, tgtZ,  xsAA)
        avn = HeavyG2massN(ihg) * xspA /xsAA
!             this inlcude successive collision inside the target.
!            Cosmos needs the first collision inside the target
!            so divide this by the average number of successive
!            collisions. However, this theory has no firm basis, 
!            sot that we put a switch do use sucave =1
!              get <Nsuc>
        if(HowIntNuc .eq. 0) then
           call caveNoOfSucC(a3, xspp, sucave)
           avn = avn/sucave
        endif
      end
