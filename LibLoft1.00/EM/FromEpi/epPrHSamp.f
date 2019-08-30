      subroutine epPrHSampP(media, Egin,  prob)
      use modEMcontrol
      implicit none
#include "Zmedia.h"
      type(epmedia):: media
      real(8),intent(in):: Egin ! GeV.   true Eg
!     real(8),intent(in):: rrho ! see epPrSampP  rho(true)/rho(nominal)
!      now in modEMcontrol      
      real*8 prob  ! output probability of Pair / X0

      real*8 ale

      real(8):: Eg

      Eg = Egin *media%rhoc 
!      Eg = Egin * rrho



      if(Eg .lt. media%cnst%PrEg1H) then
         write(0,*) ' Eg=',Eg,'<media%cnst%PrEg1H=',media%cnst%PrEg1H
         call cerrorMsg('Eg is too low for epPrHSamp',1)
!           next call results in undef for uiaev  etc.
!         call epfordebug('in epPrHSamp') 
         stop 1111

      elseif(Eg .lt. media%cnst%PrEg2H) then
         ale=log10(Eg)
         call kintp3(media%tbl%PrTXH,
     *   1, media%cnst%PrneH, media%cnst%PrLEg1H,
     *   media%cnst%PrdEH, ale, prob) 
      else
!        use 1/sqrt(E)
         prob =  media%tbl%PrTXH(media%cnst%PrneH)/
     *         sqrt((Eg/media%cnst%PrEg2H))
      endif
      end
!     ************
      subroutine epPrHSampE(media, Egin, Ee)
!     ************
      use modEMcontrol
!          samples higher energy pair electron
      implicit none
#include "Zglobalc.h"
#include "Zmedia.h"
#include  "Zmass.h"
      type(epmedia):: media
      real(8),intent(in):: Egin  ! GeV, gamma energy
!     real(8),intent(in):: rrho ! rho(true)/rho(norminal)
          ! now in modEMcontrol
      real(8),intent(out)::  Ee ! sampled pair electron energy (higher one) GeV

      real*8 u, ale,  ans
      real(8):: Eg
      real(8),parameter:: EminE2= masele*1.0001d0

      Eg = Egin * media%rhoc
!      Eg = Egin * rrho


      ale = log10(Eg)
 10   continue
      call rndc(u)
      if(Eg .lt. media%cnst%PrEg2H) then
         call k4ptdi(media%tbl%PrSTH, 
     *        media%cnst%Prnu1H,
     *        media%cnst%PrneH,
     *        media%cnst%Prnu1H,
     *        media%cnst%PrU1H,
     *        media%cnst%PrLEg1H,
     *        media%cnst%PrdU1H,
     *        media%cnst%PrdEH, u,  ale,  ans)  
         Ee = ans*Egin   ! here Egin 
      else
!           use rejection method; employ cosmos function
         call csetLPMCnst(media%s1, media%logs1,
     *       1.0d0, media%X0g*Tokgpm2) ! X0g in kg/m^2 1.0d0 is dummy
                                       ! not used for pair
         call cpairErgLPM(Egin, media%rho*media%rhoc*Tokgpm3, Ee)   ! rho in kg/m^3
!         call cpairErgLPM(Egin, media%rho*rrho*Tokgpm3, Ee) ! rho in kg/m^3
!            Ee is not nec. higher one.
!            make consistent with non LPM case.
         if(Egin - Ee .gt. Ee) then
            Ee = Egin- Ee
         endif
      endif
      if( Egin - Ee <= EminE2 ) goto 10  ! 2015.Jan.24

      end

