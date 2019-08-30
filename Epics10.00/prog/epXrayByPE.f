!      ************
      subroutine epPhotoElecGamma( Exray )
!     ************
       implicit none
#include  "Zcode.h"
#include  "Zcnfig.h"
#include  "Zmass.h"           

       real*8 eout, cost, cs, sn, sint, Exray
       logical kbtest
       type(epTrack)::  elec1, xray
!
!           get Photo-electron energy
!       call epphotoEe(Media(MediaNo).pe,   < v8.0
       call epphotoEe(Media(MediaNo),
     *      cTrack%p%fm%p(4), eout, cost)

       if(kbtest(Eabsorb, BitPhotoElec)) then
!            energy absorbed by atom is Eabs = Eshell= Eg-(Ee-Me)
          Move%dE = cTrack%p%fm%p(4) - (eout - masele)
          Move%dEeff= Move%dE
          Move%dEioni = Move%dE
          SumDe = SumDe + Move%dE
!                  regard it as deposited in the media
!          if(Det.cmp(Cn).CountDE .ge. 1) then >>>>>>>>>>>>>>light
             call epLightPreUserde(1, cTrack)
             if( Move%Abort /= 0 ) then
                if( Move%Abort /=3 ) then
                   call epempty ! empty the stack
                   call epSkipUpdateNo
                else
                   Move%Abort=0
                endif
                  ! no flag is needed. since called from epint
                return
             endif
!          endif                               <<<<<<<<<<<<<
          Exray = 0.
       else
!             bit 1 is not on; characteristic x-ray emmission;
!             This was neglected in v8.71 or earlier.
!             we  assume 
!                1)  p.e effect takes place for the largest possible
!                    shell energy (Say, if Eg> K-shell energy, L-shell
!                    p.e effect is neglected and all p.e effect  is assumed
!                    to take for K-shell. 
!                2)  For such p.e effect, vacancy of electron level is 
!                    filled by X-ray emission ; No Auger electron emmission
!                    is considered. 1)+2) are good approximation.
         Exray = max( cTrack%p%fm%p(4) - (eout - masele), 0.d0)
       endif
