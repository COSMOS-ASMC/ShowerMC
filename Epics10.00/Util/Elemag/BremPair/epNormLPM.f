      subroutine epNormLPMc(media, norm)
      use BPLPM
      implicit none
!           This normalization const should be used for
!       complete screening + LPM function normalization.
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmass.h"      
#include "Zmedia.h"

       type(epmedia)::  media
      real*8  norm          ! output.  LPM function * norm 
                            !          becomes equal to complete
                            !          screening case.

      real*8 dsdxH, dsdxcs, dz

      real*8 epCompScrBrs

      real*8 xnorm/0.75d0/  ! normalization point
      integer i

!        normalize the cross-section at x=xnorm, Ee=EemaxL to the
!        complete screening case.
!//////////////
      write(0,*) ' this epNormLPMc should not be used now'
      stop
!////////////////

      dsdxH = 
     *   epBremSH(media,   media%cnst%CompScrE, xnorm) 
      dsdxcs =  epCompScrBrs(media,  xnorm) 
      norm = dsdxcs/dsdxH
!/////////////
      write(0,*) ' norm =',norm, 'dsdxH=',dsdxH, 
     *    ' dsdxcs=', dsdxcs
!/////////////
      end
!     ***************************
      subroutine epNormLPMp(media, norm)
      use BPLPM
      implicit none
!           This normalization const should be used for
!       complete screening + LPM function normalization
!       of pair creation.
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmass.h"      
#include "Zmedia.h"

       type(epmedia)::   media
      real(8),intent(out)::  norm     !   LPM function * norm 
                            !          becomes equal to complete
                            !          screening case.

      real*8 dsdxH, dsdxcs, dz  ! dz not used

!      real*8 epPairH, epCompScrPrs
      real*8  epCompScrPrs

      real*8 xnorm/0.75d0/  ! normalization point
      integer i

!        normalize the cross-section at x=xnorm, Eg=EgmaxL to the
!        complete screening case.

!      dz = media.Zeff
      dz = media%Z

      dsdxH =  epPairSH(media, media%cnst%PairEgmaxL, xnorm) 
      dsdxcs = epCompScrPrs(media,  xnorm) 

      norm = dsdxcs/dsdxH
      end
