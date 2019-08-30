      subroutine epBPnormXs( media, how )
      use BPLPM
!          normalize brems (and  pair;not yet) cross-sections in different
!        energy regions.
!     For a given medium, (media.cnst)

!          
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
       type(epmedia)::  media
      integer,intent(in)::how ! =1 assume Seltzer one is correct
                       !  and PS is normalized to Seltezer
                       ! CS(Tsai) is normazlied to the normalized PS
                       ! CS used for LPM is normalized to
                       ! normalized CS (Tsai)
                       ! =-1 assume CS is correct. and PS
                       ! is normalized to CS and Seltzer is
                       ! normalzied to the normalzied PS.
                       ! LMP is normalzied to CS
                       ! 0--> no normalization

      integer::i
      real(8),save:: xnorm=0.75d0  ! functions are normalized at this x=Eg/Ee
      real(8):: fs, fps1, fps2, fcs1, flpm 
      real(8):: epBrSfs, epBremS, epCompScrBrs
      real(8):: Eeme

!      write(0,*) ' bpsetEeme '
      Eeme = media%cnst%BrEemaxS2 / masele  ! upper E of Seltzer
      call epBPSetEeme(Eeme)
      fs = epBrSfs(media, Eeme, xnorm)
      fps1 = epBremS(media, Eeme,xnorm)

!      write(0,*) ' bpsetEeme -2 '
      Eeme = media%cnst%CompScrE / masele
      call epBPSetEeme(Eeme)
!      write(0,*) ' epBremS '
      fps2 = epBremS(media, Eeme, xnorm)
!      write(0,*) ' epCompScrBrs'
      fcs1 = epCompScrBrs(media, xnorm)

!      write(0,*) ' epBremSH'
      flpm = epBremSH(media, Eeme*masele, xnorm)

      media%cnst%how = how
      if( how == 1 ) then
         media%cnst%NormS = 1.
         media%cnst%NormPS = fs/fps1
         
         media%cnst%NormCS = media%cnst%NormPS * fps2/fcs1
         media%cnst%NormSH = media%cnst%NormCS * fcs1/flpm
      elseif( how == -1 ) then
         media%cnst%NormCS = 1.
         media%cnst%NormSH = fcs1/flpm
         media%cnst%NormPS = media%cnst%NormCS* fcs1/fps2
         media%cnst%NormS = media%cnst%NormPS * fps1/fs
      elseif(how == 0 ) then
         media%cnst%NormCS = 1.
         media%cnst%NormSH = 1.
         media%cnst%NormPS = 1.
         media%cnst%NormS = 1.
      else
         write(0,*) ' how =',how, ' to epNormBPxs is invalid'
         stop
      endif
!   We use how=1 cross-section. 
!    To convert how=-1 cross-section, 
!     we may divide cross-section by NormS
!   To convert how =0  cross-section,  
!     we may divide the cross-section by each  NormX

!///////////
      write(0,*) ' how =', how, ' NormS=', media%cnst%NormS,
     *  ' NormPS=', media%cnst%NormPS
      write(0,*) ' NormCS =', media%cnst%NormCS, ' NormSH=',
     *         media%cnst%NormSH
!////////////////
      end  subroutine epBPnormXs
  
      subroutine epBPSetEeme(Eemein)
      implicit none
#include "Zmedia.h"
#include "ZBPgene.h"      
      real(8),intent(in)::Eemein   ! Ee/me
      Eeme = Eemein
      end subroutine epBPSetEeme

      subroutine epBPSetEgme(Egmein)
      implicit none
#include "Zmedia.h"
#include "ZBPgene.h"      
      real(8),intent(in)::Egmein   ! Eg/me
      Egme = Egmein
      end subroutine epBPSetEgme

      subroutine epBPSetMedia(mediain)
      implicit none
#include "Zmedia.h"
#include "ZBPgene.h"
       type(epmedia)::  mediain
      media = mediain
      end subroutine epBPSetMedia
      
      subroutine epBPSetForce(forcein)
      implicit none
#include "Zmedia.h"
#include "ZBPgene.h"
      character(*),intent(in):: forcein   !  "?"-->(default) relevant routine is
                                          ! is selected depending on the Energy
                                          ! For others see inside epBrgene
      force = forcein
      end subroutine epBPSetForce
