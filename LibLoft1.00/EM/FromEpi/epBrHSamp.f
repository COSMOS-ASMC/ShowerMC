      subroutine epBrHSampP(media, Eein, prob)
      use modEMcontrol
      implicit none
#include "Zmedia.h"
! #include "ZepTrackp.h"
!        Brems sampling at LPM energies.
      type(epmedia):: media
      real(8),intent(in):: Eein
!      real(8),intent(in):: rrho  ! see epBrSamp  ;now in modEMcontrol
      real(8),intent(out):: prob  ! output probability of Brems / X0

      real*8 ale
      real(8):: Ee, cf
      
!     Ee = Eein*media%rhoc ; as of 2019/Apr/24, rrho and media%rhoc take the same value
      Ee = Eein*media%rhoc

      if(Ee  <  Flpm* media%cnst%BrEe1H ) then
         call cerrorMsg(
     *      'Energy is too low for LPM brems', 0)
      elseif(Ee <  media%cnst%BrEe2H) then
         ale=log10(Ee)
         call kintp3(media%tbl%BrTXH(:,1),
     *   1, media%cnst%BrneH, media%cnst%BrLEe1H,
     *   media%cnst%BrdEH, ale, prob) 
      else
!          use  1/sqrt(E) low
         prob =  media%tbl%BrTXH(media%cnst%BrneH,1)/
     *         sqrt((Ee/media%cnst%BrEe2H))
      endif
      if( HowNormBrems == -1 ) then
        ! nothing to do
      else
         if( HowNormBrems == 1 ) then
              ! normalzie to Seltzer
            prob = prob/ media%cnst%NormS
         elseif( HowNormBrems == 0 ) then
              ! no normalization
            if( Ee <= media%cnst%BrEemaxS2 ) then
               cf = media%cnst%NormS
            elseif( Ee <= media%cnst%CompScrE ) then ! same 
               cf = media%cnst%NormPS
            else
               cf = media%cnst%NormSH  ! LPM case not NormCs
            endif
            prob = prob/cf
         else
            write(0,*) 'HowNormBrems=',HowNormBrems, 'invalid'
            stop
         endif
      endif
      end
!     ************
      subroutine epBrHSampE(media, Eein,  Eg)
!     ************
      use modEMcontrol
      implicit none
#include "Zglobalc.h"
#include "Zmedia.h"
#include "Zmass.h"

       type(epmedia):: media
       real*8 Eein
!       real(8),intent(in):: rrho ! see previous sub.  ;now in modEMcontrol
       real(8),intent(out)::  Eg  ! GeV

      real*8 u, ale, us, ans, sqrtvm, Ek
      real*8 error
      real(8):: Ee

      Ee = Eein * media%rhoc
 !     Ee = Eein * rrho     ! same as above
      Ek = Eein - masele

      call rndc(u)
      ale = log10(Ee)
      if(Ee .lt. media%cnst%BrEe2H2) then
         if(u .gt. media%cnst%BrU1H) then
!          region A
!&&&&&&&&&&&&&
!            call k4ptdi(media.tbl.BrSTHA, 
!     *        media.cnst.Brnu1H, 
!     *        media.cnst.BrneH2,
!     *        media.cnst.Brnu1H, 
!     *        media.cnst.BrU1H,
!     *        media.cnst.BrLEe1H,
!     *        media.cnst.BrdU1H,
!     *        media.cnst.BrdEH2, u,  ale,  ans)  
            call kpolintp2(media%cnst%BrU1H, 1,  media%cnst%BrdU1H,
     *         media%cnst%BrLEe1H, 1,  media%cnst%BrdEH2, 
     *         media%tbl%BrSTHA, media%cnst%Brnu1H, 
     *         media%cnst%Brnu1H,  media%cnst%BrneH2,
     *           3, 3, u, ale, ans, error)
!&&&&&&&&&&&&
!             sqrtvm = sqrt(media.cnst.BrEgminH/media.rhoc/Ee)
!            sqrtvm = sqrt(media.cnst.BrEgminH*Ek/media.rhoc/Ee)
            sqrtvm = sqrt(media%cnst%BrEgminH*Ek/Eein)
            Eg =
     *      ( (1.d0-u) * ans + sqrtvm)**2 * Eein
         else
!          region B
            us= u**(1./media%cnst%BrPow)
!&&&&&&&&&&&&&&
!            call k4ptdi(media.tbl.BrSTHB, 
!     *        media.cnst.Brnu2H,
!     *        media.cnst.BrneH2,
!     *        media.cnst.Brnu2H,
!     *        media.cnst.BrU3H,
!     *        media.cnst.BrLEe1H,
!     *        media.cnst.BrdU2H,
!     *        media.cnst.BrdEH2, us,  ale,  ans)  
            call kpolintp2(media%cnst%BrU3H, 1,  media%cnst%BrdU2H,
     *         media%cnst%BrLEe1H, 1,  media%cnst%BrdEH2, 
     *         media%tbl%BrSTHB, media%cnst%Brnu2H,
     *         media%cnst%Brnu2H,  media%cnst%BrneH2,
     *           3, 3, us, ale, ans, error)
!&&&&&&&&&&&&&&
!            Eg= exp(- ans*u )*(Ee-masele)
            Eg= exp(- ans*u )*Ek
         endif
      else
!           use rejection method; employ cosmos function
!           neglect  Eg/Ee< 1.e-6
         call csetLPMCnst(media%s1, media%logs1, 
!     *     max(1.d-6, media.cnst.BrEgminH/media.rhoc/Ee),
     *     media%cnst%BrEgminH,
     *     media%X0g*Tokgpm2)     ! X0g in kg/m^2
!     *         media.cnst.BrEgminH*Ek/media.rhoc/Ee), ! old
!     *      media.X0g*Tokgpm2)     ! X0g in kg/m^2

!         call cbremErgLPM(Ee, media.rho*media.rhoc*Tokgpm3, Eg)  ! rho in kg/m^3
         call cbremErgLPM(Eein, media%rho*media%rhoc*Tokgpm3, Eg)  ! rho in kg/m^3
      endif
      end
