      subroutine epmudEdx(flagN, flagBr, flagPr, media, Emu, dEdx)
      implicit none
#include "Zmedia.h"
!          gives muon sum of the energy loss due to direct pair creation,
!          bremsstrahlung and nuclear interaction.
!          dE/dx  by ionization is not included here.
!
      integer flagN  ! input.  specifies how to treate nuclear interaction.
                     !         0--> n.i is neglected always.
                     !         1-->int(v0:vmax) of Emu* v*dsigma/dv is
                     !             put in dEdx.v0~0, vmax~1.0.
                     !             (v=E_loss/Emu)
                     !            I.e., n.i is treated as a continuous process.
                     !         2,3-->int(0:vmin) of Emu* v*dsigma/dv is
                     !             put in dEdx. vmin=cnst.muNVmin(=10^-3,-4)
                     !         (discrete process  by  v>vmin is sampled
                     !          for each case; 2--> no n.i is followed.
                     !          3--> n.i is followed)
                     !          
      integer flagBr  ! input.  The same meaning for bremsstrahlung as flagN.
      integer flagPr  ! input   //                   direct pair creation.
       type(epmedia):: media  ! input. media
      real*8  Emu     ! input. muon total energy
      real*8  dEdx    ! output.   dE/dx GeV/(g/cm^2). sum of loss due to
                      !           the three  process.
!
      real*8  dEdx1, dEdx2,  dEdx3
!
      call epmuNdEdx(flagN,  media, Emu, dEdx1)
      call epmuBrdEdx(flagBr, media, Emu, dEdx2)
      call epmuPrdEdx(flagPr, media, Emu, dEdx3)
      dEdx = dEdx1 + dEdx2 + dEdx3
      end
!     ********************
      subroutine epmuNdEdx(flag, media, Emu, dEdx)
!     ********************
      implicit none
#include "Zmedia.h"
!         compute dE/dx of muon by nuclear interaction
      integer flag         ! input. 0--> dEdx is always 0.
                           !        1--> dEdx is for all v
                           !        2--> dEdx is for v<vmin
                           !             nuclear interaction is 
                           !             not followed
                           !        3--> dEdx is for v<vmin
                           !             n.i is followed
       type(epmedia):: media  ! input. media
      real*8  Emu          ! input. muon total energy in GeV.
      real*8  dEdx         ! output. dE/dx GeV/(g/cm2)


      real*8 ale, pw

      if(flag .eq.0 .or.  (Emu .lt. media%cnst%muNEmin)) then
         dEdx = 0.
      else
         if(Emu .gt.  media%cnst%muNEmax1) then
            ale = log10(media%cnst%muNEmax1)
         else
            ale =  log10(Emu)
         endif
         if(flag .eq.  1) then
!              all is  continuous loss
            call kintp3(media%tbl%MuNdEdxt,
     *      1, media%cnst%muNTXT, media%cnst%muNLEmin,
     *     media%cnst%muNdETX, ale, dEdx)
            pw = media%cnst%muNpwdEdxt
         elseif(flag .eq. 2 .or. flag .eq. 3) then
!              v<vmin is regarded as continuous loss
            call kintp3(media%tbl%MuNdEdx0,
     *      1, media%cnst%muNTXT, media%cnst%muNLEmin,
     *     media%cnst%muNdETX, ale, dEdx)
            pw = media%cnst%muNpwdEdx0
         else
            call cerrorMsg('flag is invalid for epmuNdEdx',0)
         endif
         if(Emu .gt. media%cnst%muNEmax1) then
            dEdx = dEdx*(Emu/ media%cnst%muNEmax1)**pw
         endif
         dEdx = dEdx * Emu
      endif
      end
!     *******************
      subroutine epmuBrdEdx(flag, media, Emu, dEdx)
!     *******************
      implicit none
#include "Zmedia.h"
!         compute dE/dx of muon by brems
      integer flag         ! input. 0--> dEdx is always 0.
                           !        1--> dEdx is for all v
                           !        2/3--> dEdx is for v<vmin
       type(epmedia):: media  ! input. media
      real*8  Emu          ! input. muon total energy in GeV.
      real*8  dEdx         ! output. dE/dx GeV/(g/cm2)


      real*8 ale

      if(flag .eq.0 .or.  (Emu .lt. media%cnst%muBrEmin)) then
         dEdx = 0.
      else
         if(Emu .gt.  media%cnst%muBrEmax1) then
            ale = log10(media%cnst%muBrEmax1)
         else
            ale =  log10(Emu)
         endif
         if(flag .eq.  1) then
!              all is  continuous loss
            call kintp3(media%tbl%MuBrdEdxt,
     *      1, media%cnst%muBrTXT, media%cnst%muBrLEmin,
     *     media%cnst%muBrdETX, ale, dEdx)
         elseif(flag .eq. 2 .or. flag .eq. 3) then
!              v<vmin is regarded as continuous loss
            call kintp3(media%tbl%MuBrdEdx0,
     *      1, media%cnst%muBrTXT, media%cnst%muBrLEmin,
     *     media%cnst%muBrdETX, ale, dEdx)
         else
            call cerrorMsg('flag is invalid for epmuBrdEdx',0)
         endif
         dEdx = dEdx * Emu
      endif
      end

!     *******************
      subroutine epmuPrdEdx(flag, media, Emu, dEdx)
!     *******************
      implicit none
#include "Zmedia.h"
!         compute dE/dx of muon by pair creation
      integer flag         ! input. 0--> dEdx is always 0.
                           !        1--> dEdx is for all v
                           !        2,3--> dEdx is for v<vmin
       type(epmedia):: media  ! input. media
      real*8  Emu          ! input. muon total energy in GeV.
      real*8  dEdx         ! output. dE/dx GeV/(g/cm2)


      real*8 ale

      if(flag .eq.0 .or.  (Emu .lt. media%cnst%muPrEmin)) then
         dEdx = 0.
      else
         if(Emu .gt.  media%cnst%muPrEmax1) then
            ale = log10(media%cnst%muPrEmax1)
         else
            ale = log10(Emu)
         endif
         if(flag .eq.  1) then
!              all is  continuous loss
            call kintp3(media%tbl%MuPrdEdxt,
     *      1, media%cnst%muPrTXT, media%cnst%muPrLEmin,
     *     media%cnst%muPrdETX, ale, dEdx)
         elseif(flag .eq. 2 .or. flag .eq. 3) then
!              v<vmin is regarded as continuous loss
            call kintp3(media%tbl%MuPrdEdx0,
     *      1, media%cnst%muPrTXT, media%cnst%muPrLEmin,
     *     media%cnst%muPrdETX, ale, dEdx)
         else
           call cerrorMsg('flag is invalid for epmuPrdEdx',0)
         endif
         dEdx = dEdx * Emu
      endif
      end

