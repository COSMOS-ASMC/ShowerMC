!       interaction routines for light
!         light interaction:  Rayleigh
!                             attenuation (+ wls )
!                 fix which process and set path in  Move.
      subroutine epproLight

      use modepLightPty
      implicit none
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"
#include  "Zcnfig.h"


      real(8):: probRayl, tRayl, probAtten, tAtten
      real(8):: t    ! path length in r.l
      real(8):: qeff  ! Qeff of the  photo sensor

      if( cLcompNo <= 0 ) then
         ! non light material
         Move%dl = 0.
         Move%dt = 0.
         Move%dx = 0.
         Move%proc = "absorb"    ! instant absroption
      elseif( Lcomp( cLcompNo )%refracN == 0.) then
         !  this  component should be photo sensor
         Move%dl = 0.
         Move%dt = 0.
         Move%dx = 0.
         Move%proc = "pe"    ! p.e production
      else         
!///////////
!         call Lcompchk(' proc ', cptyNo)
!//////////////
         if( comInfo( cPtyNo )%Rayleigh == 0.) then
            probRayl = 0.
            tRayl = 10000.
         else
            call epLightRayleigh( 
     *      Lcomp(cLcompNo)%refracN,  probRayl, tRayl)
         endif
          ! attenuaiton
         call epLightAtten(Lcomp( cLcompNo ), 
     *        comInfo( cPtyNo ),   probAtten, tAtten)

         if(tRayl < tAtten)  then
            t = tRayl
            Move%proc = "rayl"
         else
            t = tAtten
            if( comInfo( cPtyNo )%WLS(1) == 0. ) then
               Move%proc = "absorb"
            elseif( comInfo( cPtyNo )%WLS(1) > 0. ) then
               Move%proc = "wls"
            else
               write(0,*) ' WLS=', comInfo( cPtyNo )%WLS, ' are invalid'
               write(0,*) ' cPtyNo=', cPtyNo
               stop
            endif
         endif
         Move%dt = t
!////////////////////
!         if(MediaNo <= 0) then
!            write(0,*) ' MediaNo=', MediaNo, ' in proc'
!            stop
!         endif
!//////////////
         Move%dx = Move%dt * Media(MediaNo)%X0g
         Move%dl = Move%dx * Media(MediaNo)%gtocm /
     *        Media(MediaNo)%rhoc
      endif
      end 

      subroutine epLightRayleigh( n, prob, path )
!      use modepLight
      use modepLightPty
      implicit none 
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"

      real(8),intent(in)::n  ! refracton index
      real(8),intent(out):: prob ! prob of Rayleigh scatt. /r.l
      real(8),intent(out):: path  ! sample path for //   in r.l


      real(8)::u, n2
      real(8),parameter::RaylC=24*3.141592**3*1.e27  ! to Rayleigh scat xsec in mb
      real(8):: xs   ! x-section in mb

!///////////////
!      if(MediaNo <= 0) then
!         write(0,*) 'MediaNo =', MediaNo, ' in raylei'
!         stop
!      endif
!///////////////
      if( Media(MediaNo)%gasF .eq. 1 ) then
!           media is gas
         n2 =n**2
         xs = Raylc * (1.0/(cTrack%wl*1.e-7))**4/   ! w.l in cm
     *        Media(MediaNo)%ndensity**2 *
     *        ( (n2-1)/(n2+2) )**2
      else 
         xs = 1.e-20  !  to be updated
      endif
      prob = xs*Media(MediaNo)%mbtoPX0

      call rndc(u)
      path  = -log(u)/prob
      end

      subroutine epLightAtten(each,  pty, prob, path )
!      use modepLight
      use modepLightPty
      implicit none 
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"

      type(eachInfo),intent(in)::each  ! variable part of property
      type(property),intent(in)::pty  ! common property 
      real(8),intent(out)::prob       ! prob /r.l
      real(8),intent(out)::path       ! sampled path in r.l


      real(8):: attenL, u

      if( pty%attenFile == "dummy") then
         if( pty%attenL == 0.) then
            attenL = 1.e12    ! no atten
         else     
            attenL= pty%attenL    !  cm
         endif
!      elseif( each%attenAF > 0 ) then
      elseif( pty%attenAF > 0 ) then
             ! fix attenL
!         call csampAF( pty%attenAF, attenL)
         call csampAFIntp(pty%attenAF, cTrack%wl, attenL)
      else
         write(0,*) ' strange: attenFile != "dummy" '
         write(0,*)
     *    ' but attenAF =0 in eachInfo; in pLightAtten '
         stop
      endif
!      attenL = attenL/ Media(MediaNo).X0  ! in r.l
      prob = Media(MediaNo)%X0/attenL     ! / r.l
      call rndc(u)
      path = - log(u)/prob
      end


!       psudo ptcl:      energy loss for scintillation
      subroutine epproEdepo
!      use modepLight
!      use modepLightPty
      implicit none
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"


      
      Move%dt = 0.
      Move%dx = 0.
      Move%dl = 0.
      end
!       psudo ptcl:   charged ptcl track for Cerenkov
      subroutine epprochgPath
!      use modepLight
!      use modepLightPty
      implicit none
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"

      
      Move%dt = 0.
      Move%dx = 0.
      Move%dl = 0.
      end
