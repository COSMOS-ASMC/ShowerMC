!     scintillation ligh is generated from energy deposit
!     when a charged partilce runs and deposit energy
!     This is different from ScintiFromCell
!         
      subroutine epLightScintiOnFly(info, aTrack)
      use modepLight
      use modepLightPty
      implicit none
#include "Zcode.h"
#include "ZepTrackp.h"
#include "ZepTrackv.h"
#include "Zcnfig.h"
#include "ZsepManager.h"

      integer,intent(in)::info
       type(epTrack)::  aTrack  ! charged ptcl at segment top



      real(8):: wl, wl0    ! wave length 
      real(8):: photons ! <no. of photons from scinti>
      integer:: nphotons ! with poisson fluctuatin
      integer:: npreal   ! really sampped
!      integer,pointer:: Nsmp   ! max of really sampled #
      integer:: Nsmp   ! max of really sampled #
       type(epTrack)::  scintilight
      
      integer:: i
      real(8):: u, weight

            !                        energy loss in in GeV
      photons = comInfo(cPtyNo)%NpPerMeV*Move%dEeff*1.e3 ! MeV-->N photons
      if( photons < usePoisson ) then
         call kpoisn(photons, nphotons)
      else
         nphotons = photons + 0.5
      endif
      Nsmp = comInfo( cPtyNo )%NpSample

      scintilight = aTrack  ! inherit some from charged ptcl
!             make light:       scinti's subcode is 1
      call cmkptc(klight, 1, 0, scintilight%p)
      if(nphotons > Nsmp) then
         ! use  Nsmp 
         npreal = Nsmp
         weight = 1. +  real(nphotons-Nsmp)/Nsmp
!///////////
!         write(0,*) ' w ', weight
!///////////
      else
         npreal = nphotons
         weight = 1.
      endif

      sumni = sumni + npreal
      sumniwi = sumniwi + npreal* weight

      scintilight%p%subcode = 1 
      scintilight%wgt = weight

      do i = 1, npreal
         if( info == 0 ) then
!               light generated uniformly along path
            call rndc(u)
            scintilight%pos%x = aTrack%pos%x + u*Move%dl*aTrack%w%x
            scintilight%pos%y = aTrack%pos%y + u*Move%dl*aTrack%w%y
            scintilight%pos%z = aTrack%pos%z + u*Move%dl*aTrack%w%z
!        else  
!           info =1 ==> use same pos as aTrack.pos
         endif
!              get w.l
!//////////
!         call Lcompchk("Za", cLcompNo)
!///////////
         if( comInfo( cPtyNo )%waveAF > 0  ) then
!                distribution exists; sample wl
            call csampAF( comInfo( cPtyNo )%waveAF, wl)
!/////////////
!            write(*,*) 'wl ', wl
!///////////// 
         else
!                use always  fixed wl
            wl = comInfo(cPtyNo)%peakWL
         endif
!             put energy for same use; wl and  wl0 is the same
!             we use always same w.l in any media
         call epLightwl2E(wl, 1.d0, wl0, scintilight%p%fm%p(4))  

!              set wl
         scintilight%wl = wl
!              set isotropic angle
         call episoAngle(scintilight%w)
!            call eptransVect(aTrack.w, xray.w, xray.w)  ! since iso.
!            call epe2p(xray)   ! p is not used for light
!              push light
         call eppush(scintilight)
      enddo
      end
