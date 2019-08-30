!     **********************************
      function epPrgene(media, Eg, x) result(ans)
      use BPLPM
      implicit none
!
!     Generic pair creation function for entire energy region.
!     The same note as for epPrgenex
!
!     Eg < PairNonSc ==> Original B.H is modified
!     Eg < PairScrE  ==> partial screening only 
!     PairWcrE< Ee < PairEgmaxL===>complete sc. 
!     Ee < Eg2H      ===> complete sc.+ LPM  table
!     Ee > Ee2H      ===> comp. + LPM  by rejection
!
!
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZepTrackp.h"

       type(epmedia)::  media  ! input media
      real(8),intent(in)::  Eg  !  photon energy in GeV.
      real(8),intent(in)::   x  !  x = Ee/Eg. (Ee one of pair
                                ! electron  energy)
      real(8)::ans  !  pair creation cross-section  ds/dx in mb


      real*8  epPairS
      real*8  epCompScrPrs
      real(8):: Egme
!
      integer i
      real*8 f, normf
!////////////////////
!      write(0,*) 'in epPrgene LPMeffct=',LPMeffect 
!////////////////
      Egme = Eg/masele
      if(Eg .lt. media%cnst%PairEgmin) then
         ans  = 0.
      elseif(Eg .le.  media%cnst%PrScrE) then
!       A)   normally  next one holds
!               ------|-----------|---------
!                     PrScrE     PairEgmaxL  
!                    =CompScr
!
!       B)   next situation may not happen
!               ------|----------|-------
!                    PrScrE     CompScrE
!                 =PairEgmaxL
!           partial screening or no screening
         ans = epPairS(media, Egme, x)
      elseif(Eg .le. media%cnst%PairEgmaxL .and.
     *       Eg .gt.   media%cnst%CompScrE) then
!          In the case of B, this dose not happen
!          complete screening without LPM
         ans = epCompScrPrs(media,  x)
      elseif( .not. LPMeffect ) then
         ans = epCompScrPrs(media,  x)
      else
!          complete screening +  LPM
         call  epNormLPMp(media,  normf)
         ans = epPairSH(media, Eg, x)*normf
      endif

      end
!     ************************************
      subroutine epPrgeneTX(xmin, xmax, tx)
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"


      real*8  xmin  ! input.  Ee/Eg min.
      real*8  xmax  ! input.  Ee/Eg max.
      real*8  tx    ! output.  Integral of pair function from
                    !         x= xmin to xmax
      external epPrgenex
      real*8   epPrgenex,  ans1, ans2
      real*8  Eg, d,  vt

      Eg = Egme *masele

      d=(xmax-xmin)/20.d0
      vt=xmax-d
      call k16pGaussLeg(epPrgenex, xmin, vt, 16,  ans1)
      call k16pGaussLeg(epPrgenex, vt, xmax, 16,  ans2)
      tx = ans1+ans2

      end
!     **************
      real*8 function epPrgeneSolv(v)
!     **************
      implicit none

!          used to solve total-cross-section * u = integral of
!          pair function from min to v.
!
!
!
      common/upsic/upsi,vmax
      real*8 upsi, vmax

      real*8 v

      real*8 ans

      call epPrgeneTX(v, vmax, ans)

      epPrgeneSolv = ans/upsi-1.d0

      end
