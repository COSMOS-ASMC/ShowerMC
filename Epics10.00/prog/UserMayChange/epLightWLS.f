!  if you want to modify the wave length shift treatment, give wls(1) 
!  a value > 1 and implement a code fragment below by judging 
!  the value of wls(1).
      subroutine epLightWLS( wls ) 
! act as  a wave length  shifter. 
! given wave length is shited 
!   wls(1:4) control how to shift;  values of the WLS can be given in
!    Light10.dat etc 
!  (see UsterHook/TestLight/moreRealistic/LightData/Light10.dat )
! wls(1): 1.0--> default treatment.
!        if wave length given in cTrack is in between wls(2) ~ wls(3)
!        wave length is changed to wls(4) (nm).
!        else the light will disappear
! wls(1): > 1.0 : The user must make a procedure how to
!        shift the wave length. The coding must be given 
!        below (after %%%%%%  )
! wls(1): other. error.
!
! >>>>>>>>>>>>>>>>>>>>>>>light
      use modepLightPty
!c <<<<<<<<<<<<<<<<<<<<<light
      implicit none
#include "Zcode.h"
#include "ZepTrackp.h"
#include "ZepTrackv.h"
#include "Zcnfig.h"
      real(8),intent(in)::wls(4)

       type(epTrack)::  aTrack
      real(8)::wl0
      
      if( wls(1) == 1. ) then
         if( cTrack%wl >= wls(2) .and.
     *       cTrack%wl <= wls(3) ) then
            aTrack = cTrack  ! copy cuurent light
            aTrack%wl = wls(4)  ! change wl
                 ! convert wl into energy( wl0=wl refrac index=1.)
            call epLightwl2E(aTrack%wl, 1.d0,
     *       wl0,   aTrack%p%fm%p(4))
             ! set isotropic  angle
            call episoAngle(aTrack%w)
                 ! stack the light
            call eppush(aTrack)
         else
            !  no wave length shifted light appear
         endif
         !     %%%%%%%%%%%%%%%%%%
      elseif( wls(1) >= 1.d0) then
         write(0,*) ' you must add something in epLightWLS'
         write(0,*) ' current WLS values are :' , wls
         write(0,*) 'see .... part in UserHook/epLightWLS%f'
         stop
!         you can use these
!      current light is cTrack
!      Cn
!      MediaNo 
!      Media(MediaNo).rhoc
!      cLcompNo 
!      cPtyNo ( = Lcomp(cLcompNo)%comInfoNo )
!      Lcomp( cLcompNo )%refracN 
      else
         write(0,*) ' current WLS values are :' , wls
         write(0,*) ' and invalid'
         stop
      endif
      end subroutine epLightWLS
