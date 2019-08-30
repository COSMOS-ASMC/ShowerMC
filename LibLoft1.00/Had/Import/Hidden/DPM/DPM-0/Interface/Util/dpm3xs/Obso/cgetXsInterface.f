!    this must be modified depending on the int. model 
!      this is for dpmjet3
!    SEE also CosmosDat/UserHook/GetXsecDPM/Fit/toSpAbySp.f90
!   for getting true dpmjet3 x-sections > 50 GeV
!   for pp, pip, Kp subroutine cSppdpm,Spipdpm,SKpdpm in 
!   Cosmos/Particle/Event/Xsection/cSxpdpm.f90 can be used.
#include "BlockData/cblkGene.h"
      subroutine cgetXsIni
#include "Zevhnp.h"
!  #include "Zmanagerp.h"
      SxAbySxpOpt = 0
      end subroutine cgetXsIni
!          This is somewhat totology; dpm cannot
!       get x-section on fly.   
!   For original dpmjet3 xsection, 
      subroutine cgetXsInterface(pj, tg, xs)
      implicit none
#include  "Zglobalc.h"
#include  "Zcode.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
!        
      type(ptcl):: pj    ! input projectile
      type(ptcl):: tg    ! input target
      real(8),intent(out):: xs ! in mb

      real(8):: mfp

      real(8):: A, Z
      

      TrackBefMove.p = pj
      if( tg.code == kgnuc ) then
         TargetMassN = tg.subcode
      else
         TargetMassN =  1.
      endif
!      call cmfpdpmjet3(xs, mfp)  ! this gives total Xs at E< 4.1 GeV!      call cmfpOther(xs, mfp)     ! only inelastic
      Z = tg.charge
      if(tg.code == kgnuc ) then
         A = tg.subcode
      else
         A = 1
      endif
      call cinelx(pj, A, Z, xs)
      end      subroutine cgetXsInterface

