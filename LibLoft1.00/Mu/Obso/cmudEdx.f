      subroutine cmudEdx(flagN, flagBr, flagPr,  Emu, dEdx, dEdxF)
      implicit none

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
                     !         2/3-->int(0:vmin) of Emu* v*dsigma/dv is
                     !             put in dEdx. vmin=cnst.muNVmin(=10^-3)
                     !         (discrete process  by  v>vmin is sampled
                     !          for each case; 2--> no n.i is followed.
                     !          3--> n.i is followed)
                     !          
      integer flagBr  ! input.  The same meaning for bremsstrahlung as flagN.
      integer flagPr  ! input   //                   direct pair creation.
      real*8  Emu     ! input. muon total energy
      real*8  dEdx    ! output.   dE/dx GeV/(g/cm^2). sum of loss due to
                      !           the three  process.  restricted one 
                      !           if frag is 2/3
      real*8  dEdxF !  output  Full dE/dx. if flag =1, same as dEdx
!
      real*8  dEdx1, dEdx2,  dEdx3
      real*8  dEdx1F, dEdx2F,  dEdx3F
!
      call cmuNdEdx(flagN,  Emu, dEdx1)
      if(flagN .ne. 1) then
         call cmuNdEdx(1,  Emu, dEdx1F)
      else
         dEdx1F = dEdx1
      endif

      call cmuBrdEdx(flagBr, Emu, dEdx2)
      if(flagBr .ne. 1) then
         call cmuBrdEdx(1, Emu, dEdx2F)
      else
         dEdx2F = dEdx2
      endif
      call cmuPrdEdx(flagPr, Emu, dEdx3)
      if(flagPr .ne. 1) then
         call cmuPrdEdx(1, Emu, dEdx3F)
      else
         dEdx3F = dEdx3
      endif
      dEdx = dEdx1 + dEdx2 + dEdx3
      dEdxF = dEdx1F + dEdx2F + dEdx3F
      end
!     ********************
      subroutine cmuNdEdx(flag, Emu, dEdx)
!     ********************
      implicit none
#include "Zcmuint.h"
!         compute dE/dx of muon by nuclear interaction
      integer flag         ! input. 0--> dEdx is always 0.
                           !        1--> dEdx is for all v
                           !        2--> dEdx is for v<vmin
                           !             nuclear interaction is 
                           !             not followed
                           !        3--> dEdx is for v<vmin
                           !             n.i is followed
      real*8  Emu          ! input. muon total energy in GeV.
      real*8  dEdx         ! output. dE/dx GeV/(g/cm2)


      real*8 ale, pw

      if(flag .eq.0 .or.  (Emu .lt. muNEmin)) then
         dEdx = 0.
      else
         if(Emu .gt.  muNEmax1) then
            ale = log10(muNEmax1)
         else
            ale =  log10(Emu)
         endif
         if(flag .eq.  1) then
!              all is  continuous loss
            call kintp3(MuNdEdxt, 1, muNTXT, muNLEmin,
     *       muNdETX, ale, dEdx)
            pw = muNpwdEdxt
         elseif(flag .eq. 2 .or. flag .eq. 3) then
!              v<vmin is regarded as continuous loss
            call kintp3(MuNdEdx0, 1, muNTXT, muNLEmin,
     *      muNdETX, ale, dEdx)
            pw = muNpwdEdx0
         else
            call cerrorMsg('flag is invalid for cmuNdEdx',0)
         endif
         if(Emu .gt. muNEmax1) then
            dEdx = dEdx*(Emu/muNEmax1)**pw
         endif
         dEdx = dEdx * Emu
      endif
      end
!     *******************
      subroutine cmuBrdEdx(flag, Emu, dEdx)
!     *******************
      implicit none
#include "Zcmuint.h"
!         compute dE/dx of muon by brems
      integer flag         ! input. 0--> dEdx is always 0.
                           !        1--> dEdx is for all v
                           !        2/3--> dEdx is for v<vmin
      real*8  Emu          ! input. muon total energy in GeV.
      real*8  dEdx         ! output. dE/dx GeV/(g/cm2)


      real*8 ale

      if(flag .eq.0 .or.  (Emu .lt. muBrEmin)) then
         dEdx = 0.
      else
         if(Emu .gt.  muBrEmax1) then
            ale = log10(muBrEmax1)
         else
            ale =  log10(Emu)
         endif
         if(flag .eq.  1) then
!              all is  continuous loss
            call kintp3(MuBrdEdxt, 1, muBrTXT, muBrLEmin,
     *      muBrdETX, ale, dEdx)
         elseif(flag .eq. 2 .or. flag .eq. 3) then
!              v<vmin is regarded as continuous loss
            call kintp3(MuBrdEdx0, 1, muBrTXT, muBrLEmin,
     *      muBrdETX, ale, dEdx)
         else
            call cerrorMsg('flag is invalid for cmuBrdEdx',0)
         endif
         dEdx = dEdx * Emu
      endif
      end

!     *******************
      subroutine cmuPrdEdx(flag, Emu, dEdx)
!     *******************
      implicit none
#include "Zcmuint.h"
!         compute dE/dx of muon by pair creation
      integer flag         ! input. 0--> dEdx is always 0.
                           !        1--> dEdx is for all v
                           !        2/3--> dEdx is for v<vmin
      real*8  Emu          ! input. muon total energy in GeV.
      real*8  dEdx         ! output. dE/dx GeV/(g/cm2)


      real*8 ale

      if(flag .eq.0 .or.  (Emu .lt. muPrEmin)) then
         dEdx = 0.
      else
         if(Emu .gt.  muPrEmax1) then
            ale = log10(muPrEmax1)
         else
            ale = log10(Emu)
         endif
         if(flag .eq.  1) then
!              all is  continuous loss
            call kintp3(MuPrdEdxt, 1, muPrTXT, muPrLEmin,
     *      muPrdETX, ale, dEdx)
         elseif(flag .eq. 2 .or. flag .eq. 3) then
!              v<vmin is regarded as continuous loss
            call kintp3(MuPrdEdx0, 1, muPrTXT, muPrLEmin,
     *      muPrdETX, ale, dEdx)
         else
            write(0,*) ' flag = ', flag
            call cerrorMsg('flag is invalid for cmuPrdEdx',0)
         endif
         dEdx = dEdx * Emu
      endif
      end
