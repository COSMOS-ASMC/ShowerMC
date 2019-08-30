!          get xsection by Gheisha
      subroutine cxsecGheisha(pj, aa, zz, xs)
      implicit none
#include "Zptcl.h"

      type (ptcl):: pj  ! input. particle
      real*8 aa       ! input. average target A
      real*8 zz       ! input. average target Z
      real*8 xs       ! output. cross-section in mb

!            
!        m.f.p (kg/m**2) = abn /xsec(mb)
!
      real*8 abogn, toabn
      parameter (abogn=6.02d23, toabn=1.d28/abogn)

      integer  ghecode
      real*4 aa4, zz4, radl, rho, absl, mfp4, pin4
      real*8 pin8
      data rho/1./  ! to get MFP in g/cm2
      data radl/1./, absl/1./   ! dummy

      aa4 = aa
      zz4 = zz
      call gsmate(1, aa4, zz4, rho, radl, absl)
      call ccos2gheCode(pj, ghecode)
      call cpxyzp(pj%fm, pin8)
      pin4 = pin8
      call gpghei(1, ghecode, pin4, mfp4)
!           mfp4 is in g/cm2
!           convert it to mb.   kg/m2/toabn/aa = 1/mb 
      xs =toabn*aa/( mfp4 *10.)
      end

