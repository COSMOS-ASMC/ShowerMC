      subroutine epBrLSampP(media, Ee, prob)
      use modEMcontrol
!          partial screening region by ps1 or ps2
      implicit none
#include "Zmedia.h"
! #include "ZepTrackp.h"
       type(epmedia):: media
      real*8 Ee
      real*8 prob  ! output probability of Brems / X0

      real*8 ale
      

      if(Ee .le. media%cnst%BremEemin) then
         prob= 1.d-40
      else
         ale=log10(Ee)
         call kintp3(media%tbl%BrTXL(:,1),
     *   1, media%cnst%BremTXTL, media%cnst%BremLEemin,
     *   media%cnst%BremdETXL, ale, prob) 
      endif
      if( HowNormBrems == -1 ) then
         ! nothing to do
      else
         if( HowNormBrems == 1 ) then
            prob = prob/media%cnst%NormS
         elseif( HowNormBrems == 0 ) then
            prob = prob/media%cnst%NormPS
         else
            write(0,*) 'HowNormBrems =',HowNormBrems, 
     *              ' invalid in epBrLSampP'
            stop
         endif
      endif
      end
!     ************
      subroutine epBrLSampE(media, Ee, Eg)
!     ************
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
       type(epmedia):: media
      real*8 Ee,  Eg

      real*8 u, ale, us, ans, Ek, error
      real*8 dummy(1)
      data dummy/0.d0/
      save

      call rndc(u)
      Ek = Ee - masele
      ale = log10(Ee)
      if(u .gt. media%cnst%BremUminLA) then
!          region A
!&&&&&&&&&&&&&&&
!
!            this is a little bit faster  ( ~!0% )
!            than
!            using kpolintp2; accuracy seems OK
!          
         call k4ptdi(media%tbl%BrSTLA, 
     *        media%cnst%BremUszLA, 
     *        media%cnst%BremEsize,
     *        media%cnst%BremUszLA, 
     *        media%cnst%BremUminLA,
     *        media%cnst%BremLEemin,
     *        media%cnst%BremdULA,
     *        media%cnst%BremdEL, u,  ale,  ans)  
     
!         call kpolintp2(media.cnst.BremUminLA, 1,  media.cnst.BremdULA,
!     *        media.cnst.BremLEemin, 1,  media.cnst.BremdEL,
!     *        media.tbl.BrSTLA,   media.cnst.BremUszLA, 
!     *        media.cnst.BremUszLA, media.cnst.BremEsize,
!     *         3, 3, u, ale, ans, error) 
!         Eg= exp( ans*(1.-u))*media.cnst.BremEgmin
         Eg= exp( ans*(1.-u))*media%cnst%BremEgmin * Ek
      else
!         region B
         us = u**0.25d0
!&&&&&&&&&&&&&&&&&
         call k4ptdi(media%tbl%BrSTLB, 
     *        media%cnst%BremUszLB, 
     *        media%cnst%BremEsize,
     *        media%cnst%BremUszLB, 
     *        0.d0,
     *        media%cnst%BremLEemin,
     *        media%cnst%BremdULB,
     *        media%cnst%BremdEL, us,  ale,  ans)  
 
!        call kpolintp2(dummy, 1,  media.cnst.BremdULB,
!     *        media.cnst.BremLEemin, 1,  media.cnst.BremdEL,
!     *        media.tbl.BrSTLB,   media.cnst.BremUszLB, 
!     *        media.cnst.BremUszLB, media.cnst.BremEsize,
!     *         3, 3, us, ale, ans, error) 
!&&&&&&&&&&&&&&
!         Eg = exp(-ans*u)*(Ee - masele)
         Eg = exp(-ans*u)*Ek

      endif
      end

