      real*8 function epmudsdvdr(E, v, r)
      implicit none
#include  "Zglobalc.h"
#include  "ZbasicCnst.h"
#include  "Zmuint.h"
#include  "Zmass.h"
      real*8 E ! input. muon total  energy in GeV.
      real*8 v ! input. emitted fractional virtual photon energy
!                        (Eg/E = (E+ + E-)/E)
      real*8 r ! input. r = (E+ - E-)/(E+ + E-)
!
!     v* dsigma/dvdr  is computed  in  mb.
!
      
      real*8  epmuphie,  epmuphim, rmax, epmuPairRmax, vmin, vmax
      real*8  epmuPairVmn, epmuvmax 



      vmin = epmuPairVmn(E)
      vmax = epmuvmax(E)


      if(v .le.  vmin .or. v .ge. vmax) then
         epmudsdvdr = 0.
      else
         rmax = epmuPairRmax(E, v)
         if(abs(r) .gt. rmax) then
            epmudsdvdr = 0.
         else
            epmudsdvdr =  2.d0/3.d0/pi*alpha*ar02*Zp2  *
     *          (1.-v) * (
     *            epmuphie(E,  v,  r) +
     *            (masele/masmu)**2 * epmuphim(E, v,  r)
     *                     )
         endif
      endif
      end
!
!     ****************************
      real*8 function epmuGzi(v, r)
      implicit none
#include "Zmass.h"
      real*8  v
      real*8  r
      epmuGzi  =(masmu/masele*v/2)**2 *(1-r**2)/(1-v)
      end
!     ****************************
      real*8 function epmuBeta(v)
      implicit none
      real*8 v

      epmuBeta = v**2/(1-v)/2
      end
!     ****************************
      real*8 function epmuYe(v, r)
      implicit none
      real*8 v, r
      
      real*8 beta, gzi, epmuBeta, epmuGzi,  r2
      

      beta = epmuBeta(v)
      gzi = epmuGzi(v, r)
      r2 = r*r
  
      epmuYe = (5.d0-r2 + 4.d0*beta*(1.+r2))/
     *  ( 2.*(1.+3.d0*beta)*log(3.d0+1.d0/gzi)
     *        -r2 -2*beta*(2.d0-r2) 
     *  )
      end
!     ****************************
      real*8 function epmuYm(v, r)
      implicit none
      real*8 v, r
      
      real*8 beta, gzi, epmuBeta, epmuGzi,  r2


      beta = epmuBeta(v)
      gzi = epmuGzi(v, r)
      r2 = r*r

      epmuYm = (4. + r2 + 3.d0*beta*(1 + r2) ) /
     * (  (1 + r2)*(1.5d0 + 2*beta)*log(3.d0 + gzi) +
     *     1.-1.5d0*r2
     *  )
      end
!     ****************************
      real*8 function epmuLe(E, v, r)
      implicit none
#include "Zmuint.h"
#include "Zmass.h"
      
      real*8 E, v, r
      
      real*8  gzi,  epmuGzi,  r2
      real*8 Ye, epmuYe, f1, f2, MebyMu
      parameter (MebyMu = masele/masmu)
  
      r2 = r*r
      Ye = epmuYe(v, r)
      gzi = epmuGzi(v, r)

      f1 = log( Ak/Z3*sqrt( (1+gzi)*(1+Ye) )/
     *           (  1 + 2*Akm2* (1+gzi)*(1+Ye) * MebyMu*(masele/E)/
     *                ( v*(1-r2) )
     *           )
     *        )


      f2 = log(1 + ( 1.5d0* MebyMu/Z3)**2 *(1+gzi)*(1+Ye) )/2

      epmuLe = f1 - f2
      end

!     ****************************
      real*8 function epmuLm(E, v, r)
      implicit none
#include "Zmuint.h"
#include "Zmass.h"
      
      real*8 E, v, r
      
      real*8  gzi,  epmuGzi,  r2
      real*8 Ym, epmuYm,  MebyMu
      parameter (MebyMu = masele/masmu)
  
      gzi = epmuGzi(v, r)
      Ym =  epmuYm(v, r)
      r2 = r*r
      epmuLm = log(  1.5d0 /MebyMu * Ak/(Z3*Z3) /
     *               (1+  2*Akm2* MebyMu*(1+gzi)*(1+Ym)*(masele/E)/
     *                 (v*(1-r2))
     *               ) 
     *            )
      end


!     ****************************
      real*8 function epmuphie(E, v, r)
      implicit none
#include "Zmuint.h"
#include "Zmass.h"
      
      real*8 E, v, r
      
      real*8 beta, gzi,  epmuGzi,  r2, epmuBeta,  epmuLe

      beta = epmuBeta(v)
      gzi = epmuGzi(v, r)
      r2 = r*r

      epmuphie  =  (
     *  ((2+r2)*(1+beta)+gzi*(3.d0+r2))*log(1+1.d0/gzi)
     *        + (1-r2-beta)/(1+gzi) -(3+r2)
     *             )  * epmuLe(E, v, r)

      end

!     ****************************
      real*8 function epmuphim(E, v, r)
      implicit none
#include "Zmuint.h"
#include "Zmass.h"
      
      real*8 E, v, r
      
      real*8 beta, gzi,  epmuGzi,  r2, epmuBeta,  epmuLm

      beta = epmuBeta(v)
      gzi = epmuGzi(v, r)
      r2 = r*r

      epmuphim = (
     *        ( (1+r2)*(1+1.5d0*beta) - (1+2*beta)/gzi *(1-r2)) *
     *            log(1+ gzi) 
     *          +
     *          gzi*(1- r2-beta)/(1+gzi) + (1+2*beta)*(1-r2)
     *           ) * epmuLm(E,  v,  r)
      end
!    ***********************************
      real*8 function epmuPairRmax(E, v)
      implicit  none
#include "Zmass.h"
      real*8 E  !  input. muon total energy
      real*8 v  !  input. (E+ + E-)/E
!      
!      compute  max of assymetry factor (E+ - E-)/(E- + E+)
!
      epmuPairRmax = (1. - 6.d0*(masmu/E)**2/(1-v)) *
     *       sqrt(1.  - 4.d0*masele/E/v)
      end
!     *************************************
      real*8 function epmuPairVmn(E)
      implicit none
#include "Zmass.h"
! #include "Zmuint.h"

      real*8 E  ! input muon total energy
      epmuPairVmn = 4*masele/E  ! why 4 ?
                          !  due to  epmuPairRmax which requires
                          !  this.  Probably,  correct one
                          !  says  2*masele/E
      end
            

