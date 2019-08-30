!           this is moved from prog ; old one is in epbndry.f
!           name of the sub is not changed, though the file 
!           name was changed
!       We do some business for Light > 0 
      subroutine epnewComp(aTrack)
      use epModify
      use modXsecMedia, xmedia=>media, xelement=>element, 
     *    dontuse=>NoOfMedia
! >>>>>>>>>>>>>>>>>>>>>>>light
      use modepLightPty
! <<<<<<<<<<<<<<<<<<<<<light
      implicit none
#include "Zcode.h"
#include "ZepTrackp.h"
#include "ZepTrackv.h"
#include "Zcnfig.h"

       type(epTrack)::  aTrack  ! input.  a track in a new 
               ! component
      integer::modi ! modifier index
!         when a particle moves  to a new component
!         the following must be done
      Cn = aTrack%cn
!///////////////////
!      if(Cn > Det.nct) then
!         ! void
!         write(0,*) ' Cn=',Cn
!      endif
!/////////////////
      MediaNo = Det%Cn2media(Cn)
      mediumNo = MediaNo
!//////////////
!      if(MediaNo <= 0 ) then
!         write(0,*) ' MediaNo=', MediaNo, ' in newcomp'
!         stop
!      endif
!////////////////////
      Media(MediaNo)%rhoc = Det%cmp(Cn)%rhoc


!       r.l                         cm   / X0;
!            1/X0 = 1/(X0g/rho) = rho/X0g -->rho*rhoc/X0g
      if( aTrack%p%code >= -1)  then   ! <<<<<<<<<<>>>>>>light
         MaxPath =abs( Det%cmp(Cn)%MaxPathL/Media(MediaNo)%X0*
     *        Media(MediaNo)%rhoc)
      endif                           ! <<<<<<<<<<>>>>>>light  
!         Det.cmp(Cn).MaxPathL may be < 0; which means
!        that the value is set by Epics and abs should be used.
!         already dummy
!      call ep2cosCond(Media(MediaNo).Aeff, Media(MediaNo).Zeff)
      if( aTrack%p%code > 0 ) then  ! <<<<<<<<<<>>>>>>light
         call ep2cosCond           
      endif                         ! <<<<<<<<<<>>>>>>light   
!    >>>>>>>>>>>>>>>>>>>>>>>>>>>light
      if(Light > 0 ) then
!             constants different from default in
!             media file (e.g Birks coef.) should 
!             not be changed here by propterty data:
!           E.g  one BGO-a may be spcified as Light comp.
!             others BGO-x may not be so.  Then if Birks coef.
!             is changed for BGO-a, there is  no means to
!             restore it to default value.
! 
         cLcompNo = Det%cmp(Cn)%LightCompNo  
!///////////////
         if(Cn == 0 ) then
            write(0,*) ' Cn = 0 in newcomp cLcompNo =',cLcompNo
            write(0,*) ' MedaiNo =', MediaNo
         endif
!///////////////
         if( cLcompNo > 0 ) then
            cPtyNo = Lcomp(cLcompNo)%comInfoNo
            if( aTrack%p%code == klight) then
               call epLightwl2N( aTrack%wl, Lcomp(cLcompNo), 
     *                ComInfo( cPtyNo ), Media(MediaNo)%gasF,
     *                Lcomp( cLcompNo )%refracN )
            endif
         else
            cPtyNo = 0
         endif
      endif         
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      call epSetEmin(Cn)  ! next will be done

!      EminGamma =  Det.cmp(Cn).EminG
!      EminElec  =  Det.cmp(Cn).EminE
!      RecoilKEmin = Det.cmp(Cn).RecoilE
!      KEmin = KEminsave
!      EminH = Enminsave
!      if( Det.cmp(Cn).modifier > 0 ) then
!         modi =Det.cmp(Cn).modifier
!         if( IBITS( modify(modi)%kind, bitEmin, 1) > 0) then
!            KEmin = modify(modi)%Em%KEmin
!            EminH = modify(modi)%Em%Enmin
!         else
!            KEmin = KEminsave
!            EminH = Enminsave
!         endif
!      endif
!c////////////////////////////
!      write(0,'(a,i3, a, 1p,3g13.3)')
!     *   'Cn=',Cn,'Emn=', EminGamma, EminElec, RecoilKEmin
!//////////////////////

      end
      subroutine epLightGasF(cno, gasf)
!          returns gasF of component with #, cno
      implicit none
#include "Zcode.h"
#include "ZepTrackp.h"
#include "ZepTrackv.h"
#include "Zcnfig.h"
      integer,intent(in)::cno  ! component nubmer
      integer,intent(out):: gasf  ! obtined gas factor

      integer:: mno

      mno = Det%Cn2media(cno)
!/////////////
!      if( mno <= 0 ) then
!         write(0,*) ' mno=', mno, "in  gasf, cno=", cno
!         stop
!      endif
!///////////////
      gasf = Media(mno)%gasF
      end


