      module Mulscat
#include "Zglobalc.h"  
      real(8)::massratio2       ! (m/me)^2
      real(8),parameter::finesc = 1./137.036 ! fine   struc. const
      real(8),parameter::hbarc = 197.327e-13 !  MeV cm
      real(8),parameter::cnst = 4*pi*Avogn *finesc**2 *(hbarc/0.511)**2
      real(8)::gamma, beta, gbeta, betasq !  usual ones gbet=gamma*beta
      real(8)::gamma1, beta1, gbeta1, betasq1 ! usual ones at start of segment 
      real(8)::gamma2, beta2, gbeta2, betasq2 ! usual ones at end of segment
      end module Mulscat

      subroutine epXc2(mediax, aPtcl,bPtcl, xc2)
      use Mulscat
      implicit none
#include "Zmedia.h"
#include "Zptcl.h"
#include "Zmass.h"
       type(epmedia)::  mediax
       type(ptcl)::  aPtcl  ! particle at the beggining of the segment
       type(ptcl)::  bPtcl  ! particle at the end of the segment
      
      massration2 = (aPtcl%mass/masel)**2
      
      gamma1 = aPtcle%fm%p(4)/aPtcl%mass
      gamma2 = bPtcle%fm%p(4)/bPtcl%mass
      beta1 = sqrt(1.0 - 1./gamma1**2)
      beta2 = sqrt(1.0 - 1./gamma2**2)
      gbetasq = gamm1*beta1 * gam

