      subroutine qgs01init
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
#include "Zmass.h"
!
!
      type (ptcl)::pj
      integer  ntp
      integer  ia, iz
      type (ptcl)::a(*)
      type (ptcl)::b(5)
      integer k, i
      integer MONIOU
      COMMON /AREA43/ MONIOU
            
      INTEGER          NSP
      COMMON /AREA12/  NSP
      INTEGER          ICH(95000)
      DOUBLE PRECISION ESP(4,95000)
      COMMON /AREA14/  ESP,ICH
      integer  icp, iap, iat, kicp
      integer  j
      real*8  E0
      MONIOU = 0

      call PSASETC 
!      call PSASET    !  withDiffcode version  
      call XXASET
      call PSAINI 
!      call QGSPSAINI !  withDiffcode version
      return

!     ****************
      entry qgs01event(pj, ia, iz, a, ntp)
!     ****************

      if( pj%code .ne. kgnuc ) then
         iap = 1               !set proj. mass number (1-for hadron)
      else
         iap = pj%subcode
      endif
      call ccoscode2QGS( pj, kicp )
      E0 = pj%fm%p(4)/iap    !to get energy per nucleon
      iat = ia
      icp=iabs(kicp)/2+1   !primary particle class
                            ! (1- pion, 2 - nucleon, 3 - kaon)

!        icp:  projectile class
!        iap:  projectile mass number
!        iat:  target mass number

      call XXAINI(E0, icp, iap, iat)
      call PSCONF

      ntp = 0
      DO  j = 1, NSP
!
!   0 pi0    1 pi+   -1 pi-
!   2 p  -2 pbar   3 n    -3 nbar  
!   4 K+   -4 K-  5  K0s   -5 K0L
!   6 eta    10 Lambda  -10 LambdaBar
!
         ntp = ntp + 1
         a(ntp)%fm%p(1)=ESP(4, j)
         a(ntp)%fm%p(2)=ESP(3, j)
         a(ntp)%fm%p(3)=ESP(2, j)
         a(ntp)%fm%p(4)=ESP(1, j)
         call cQGScode2cos(ICH(j), a(ntp))
      enddo
      call crot3mom(pj, a, ntp)
      end
      real*8 function  QSRAN(X)
      real*8  X                 !  not used
      real*8 u
      call rndc(u)
      QSRAN = u
      end
       

