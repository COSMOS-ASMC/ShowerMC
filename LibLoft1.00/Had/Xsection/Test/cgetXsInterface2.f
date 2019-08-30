!    this must be modified depending on the int. model 
!    this is for dpmjet3
!    SEE also CosmosDat/UserHook/GetXsecDPM/Fit/toSpAbySp.f90
!   for getting true dpmjet3 x-sections > 50 GeV
!   for pp, pip, Kp subroutine cSppdpm,Spipdpm,SKpdpm in 
!   Cosmos/Particle/Event/Xsection/cSxpdpm.f90 can be used.

#include "BlockData/cblkGene.h"

      subroutine cgetXsIni
#include "Zevhnp.h"
      SxAbySxpOpt = 0  ! same as 2:  QGS2 like
      
      TotXSopt =  2    ! For pp   1: PDG  3: Totem 
                       !  2: between 1 and 3                
                       !  3 gives highest for pp.
                       !  pbarp uses PDG so pbap < pp.
      end subroutine cgetXsIni

!          This is somewhat totology; dpm cannot
!       get x-section on fly.   
!   For original dpmjet3 xsection, 
      subroutine cgetXsInterface2(pj, tg, xs)
      use modColInfo
      implicit none
#include  "Zglobalc.h"
#include  "Zcode.h"
#include  "Zptcl.h"
!#include  "Ztrack.h"
!#include  "Ztrackv.h"
!        
      type(ptcl):: pj    ! input projectile
      type(ptcl):: tg    ! input target
      real(8),intent(out):: xs(3) ! in mb . inela, total, ela

      real(8):: mfp

      real(8):: A, Z
      
!          upto next ------------, the stuf is not needed 
!      TrackBefMove%p = pj
!      if( tg%code == kgnuc ) then
!         TargetMassN = tg%subcode
!      else
!         TargetMassN =  1.
!      endif
!      TargetAtomicN = tg%charge
!!      call cmfpdpmjet3(xs, mfp)  ! this gives total Xs at E< 4.1 GeV!      call cmfpOther(xs, mfp)     ! only inelastic
!     ------------------------
      Z = tg%charge
      A = tg%subcode
      if(A < 0 ) A = 1
      call cinelx(pj, A, Z, xs(1))
      if( pj%code == kgnuc  .and. tg%code ==kgnuc )  then 
             ! For heavy pj and heavy target, only inela
         xs(2:3) = 0.
      elseif( pj%code /=  kphoton ) then
         call ctotx2(pj, A, Z, xs(2))
         xs(3)= xs(2)-xs(1)
      else
         xs(2:3)= 0.
      endif
      end      subroutine cgetXsInterface2

