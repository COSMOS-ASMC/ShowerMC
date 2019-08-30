      real*8 function epmuBrem(E, v)
      implicit none
#include  "Zmuint.h"
#include  "Zmass.h"
      real*8 E ! input. muon total  energy in GeV.
      real*8 v ! input. emitted fractional photon energy (Eg/E).
!
!      function values,( v dsigma/dv ) in mb.
!
!         compute v* dsigma(E,v)/dv in mb/target atom (single atom of
!         (Z,A)  ) of muon bremsstrahlung cross-section. (Z,A) must
!         have been given  by epmuSetCnst and  be usable thru
!         common /zmucom/
      real*8 epmuBremLogf
!
!       Differential muon brem x-section. 
!          

      epmuBrem = D* (4.*(1.-v)/3.d0 + v**2) *
     *     epmuBremLogf(E, v)  ! mb

      end
!     *******************
      real*8 function epmuBremLogf(E, v)
!     *******************
      implicit none
#include "Zmuint.h"
#include  "Zmass.h"
      real*8 E  ! input. muon total energy. in  GeV
      real*8 v  ! input fraction ; Eg/E

      real*8 MubyMe
      parameter (MubyMe = masmu/masele)

!
      real*8  delta, logf
!             
!         mim. momentum transfer in unit of Mmu
      delta = (masmu/E) * v /(1.d0-v)/2  


      logf =   Akm / (1.d0 + Akm2*delta)
      if(Z .gt.  10.) then
!            this correction not exist in Rosental
         logf = logf * 2./3.d0 /Z3
      endif
      epmuBremLogf = log(logf)
      end





