      subroutine epsyncp(aPtcl,  B, up, mfp, path)
!      use modEMcontrol
      implicit none
#include "Zptcl.h"
#include "Zep3Vec.h"
!          samples path length for synchrotron emission.

       type(ptcl):: aPtcl  ! input. electron (E is in GeV)
       type(ep3Vec)::  B   ! input. 3 components of mag. field. in T.
      real*8 up    ! output. Upsilon value (Ee/Me * B/Bcr)
      real*8 mfp   ! output. mean free path  in m of emission of a photon
                   !           with arbitray energy
      real*8 path  ! output. sampled path in m. for emission 


      real*8    u,  cupsilon, cmBremMFP

      up = cupsilon(aPtcl, B)

       
!     mfp = cmBremMFP(aPtcl%fm%p(4), up, 0.d0) *1.d2
      mfp = cmBremMFP(aPtcl%fm%p(4), up, 0.d0)  ! m
      call rndc(u)
      path = -log(u) * mfp
      end
!    ***********
      subroutine epsynce(E, up, ephoton)
!    **************
      implicit none
      real*8 E  ! input. Electron energy in GeV
      real*8 up ! input. Upsilon value (=Ee/me *B/Bcr)
!
      real*8 ephoton ! output. sampled photon energy in GeV.

      real*8   x ! fractional energy.

      call cmBremE(up, x)
      ephoton = E* x
      end
!     ************
      subroutine epmpairp(aPtcl,  B, xaio, mfp, path)
      implicit none
#include "Zptcl.h"
#include "Zep3Vec.h"
!          samples path length for synchrotron emission.

       type(ptcl):: aPtcl  ! input. photon (E is in GeV)
       type(ep3Vec)::  B        ! input. 3 components of mag. field. in T.


       real*8 xaio               ! output. Xai value (Eg/Me * B/Bcr)/2
  !     eal*8 mfp   ! output. mean free path  in cm of magnetic pair creation
       real(8),intent(out)::mfp !  mean free path  in m of magnetic pair creation
 !     real*8 path  ! output. sampled path in cm. for pair creation
       real(8),intent(out)::path ! sampled path in m for pair creation

      real*8    u,  cmPairMFP, cxai

      xaio = cxai(aPtcl, B)

      
!     mfp = cmPairMFP(aPtcl%fm%p(4), xai) *1.d2   ! cm
      mfp = cmPairMFP(aPtcl%fm%p(4), xaio) ! m
      call rndc(u)
      path = -log(u) * mfp
      end
!    ***********
      subroutine epmpaire(E, xai, e1)
!    **************
      implicit none
      real*8 E  ! input. Photon energy in GeV
      real*8 xai ! input. xai value (=Ee/me *B/Bcr/2)
!
      real*8 e1 ! output. sampled electron energy  in GeV.

      integer  nc  ! no. of trial in cmPair
      real*8   x ! fractional energy.  >=0.5 

      call cmPairE(xai, x, nc)
      e1 = E* x
      end



