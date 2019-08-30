      subroutine epBrSampP(media, Ee,  prob, path)
!     generic brems path sampling routine
      use modEMcontrol
      implicit none
#include "Zmedia.h"
!   #include "ZepTrackp.h"

      type(epmedia):: media  ! input.
      real*8 Ee                 ! input. Electron/positron energy in GeV
!     real(8),intent(in):: rrho ! rho(actual)/rho(norminal).
!     this is in modEMcontrol
      
! in Epics, this is media%rhoc.
! in Cosmos  rho(@ present loc)/media%rho
      real*8 prob          ! output. Brems prob. per r.l
      real*8 path          ! ouptut. sampled path in r.l

      real*8 u

!      logical::LPMworks   This is in modEMcontrol
!     rrho is to consider the density change in case of Air
!     So it is the same as media%rhoc for constant density media
!     in Epics, For Cosmos, (density at the current position)/(orignial density)
!         NOW:2019/04/24  rrho    and media%rhoc take the same value
      LPMworks = LPMeffect .and.
     * ( Ee > media%cnst%BremEeminLPM*Flpm/media%rhoc )
!     * ( Ee > media%cnst%BremEeminLPM*Flpm/min(media%rhoc, 1.d0) )
      
      if( LPMworks ) then
!              LPM region
         call epBrHSampP(media, Ee,  prob)
      else
         if( Ee <= EpartialSC ) then
               ! EpartialSC <= BrEemaxS2 so surely Seltzer region
            call epBrSSampP(media, Ee,  prob)
         elseif(Ee .le. media%cnst%BrScrE) then
!            screeinig region
            call epBrLSampP(media, Ee,  prob)
         elseif( Ee .gt. media%cnst%CompScrE ) then
!     complete screeing region;no LPM.  may not come here;
            call epBrCSampP(media, Ee, prob)
         else
            write(0,*) ' strange in epBrSamp '
            stop
         endif
      endif


      
      call rndc(u)
      path = - log(u)/prob
      end
!     *********************************** 
      subroutine epBrSampE(media, Ee,  Eg)
      use modEMcontrol
      implicit none
#include "Zmedia.h"
! #include "ZepTrackp.h"
       type(epmedia):: media  ! input.
       real*8 Ee                ! input. e-/e+ energy in GeV
!       real(8),intent(in):: rrho ! see epBrSamp  in modEMcontrol

       real*8 Eg                ! output. sampled Eg in GeV

!      logical LPMworks

!        no need to check here ; in modEMcontrol
!     LPMworks =LPMeffect .and.
!     * ( Ee > media%cnst%BremEeminLPM*Flpm/rrho )
!     * ( Ee > media%cnst%BremEeminLPM*Flpm/media%rhoc )
!     * ( Ee > media%cnst%BremEeminLPM*Flpm/min(media%rhoc, 1.d0))

      if( LPMworks ) then
         call epBrHSampE(media, Ee,  Eg)
      else
         if(Ee <= EpartialSC ) then
!            Seltzer region
            call epBrSSampE(media, Ee,  Eg)
         elseif(Ee <  media%cnst%BrScrE) then
!            partial screeinig region
            call epBrLSampE(media, Ee,  Eg)
         elseif(Ee  >= media%cnst%CompScrE ) then
!            complete screeing region
            call epBrCSampE(media, Ee, Eg)
         else
            write(0,*)  'logial error of Brems sampling'
            stop
         endif
      endif
      end
