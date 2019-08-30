      subroutine epmuInte
      use modXsecMedia, xmedia=>media, xelement=>element,
     *  dontuse=>NoOfMedia

      implicit none
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"
#include  "Zcode.h"
#include  "Zevhnv.h"

      real*8 Et
      integer ia, iz
      real(8)::xs, mfp
!///////////
!      integer ns
!      record /epTrack/temp
!      integer,save:: nelec=0
!////////////
      if( Move%proc .eq. 'pair' ) then
         call epmuPrsmpE(Media(MediaNo), Move%Track%p%fm%p(4),
     *                   Et)
         cTrack%p%fm%p(4) = cTrack%p%fm%p(4) - Et  ! muon
         call epe2p(cTrack)         ! adjust momentum
         call eppush(cTrack)        ! push
         if(Media(MediaNo)%mu%MuPr .eq. 3) then
!              generate pair electrons; employ eppair
!              to do so, make cTrack a gamma of energy Et
            cTrack%p%fm%p(4) = Et
            call cmkptc(kphoton, 0, 0, cTrack%p)
            call epe2p(cTrack)       ! this may not be needed ; for safety
            call eppair
!/////////////////
!            nelec = nelec +1
!            call epqstn(ns)
!            call epqsTrack(ns, temp)
!            write(*,'(a, 2i3,i5, 1p,g15.4)') 'mpair',  nelec,
!     *       temp.p.code, temp.p.charge, temp.p.fm.p(4)
!            ns = ns -1
!            call epqsTrack(ns, temp)
!            nelec = nelec +1
!            write(*,'(a, 2i3,i5, 1p,g15.4)') 'mpair',  nelec,
!     *       temp.p.code, temp.p.charge, temp.p.fm.p(4)
!////////////////
!           &&&&&9.06
         elseif( Media(MediaNo)%mu%MuPr .eq. 2) then
            Move%dE=Et
            Move%dEeff=Et
            call epLightPreUserde(0, cTrack)
!           &&&&&&
         endif
      elseif( Move%proc .eq. 'brem' ) then   
         call epmuBrsmpE(Media(MediaNo), Move%Track%p%fm%p(4),
     *                   Et)
         cTrack%p%fm%p(4) = cTrack%p%fm%p(4) - Et  ! muon
         call epe2p(cTrack)         ! adjust momentum
         call eppush(cTrack)        ! push
         if(Media(MediaNo)%mu%MuBr .eq. 3) then
!             generate brems gamma; no deflection
            cTrack%p%fm%p(4) = Et
            call cmkptc(kphoton, kcasg, 0, cTrack%p)
            call epe2p(cTrack)  ! adjust momentum
            call eppush(cTrack) ! push
!           &&&&&9.06
         elseif( Media(MediaNo)%mu%MuBr .eq. 2) then
            Move%dE=Et
            Move%dEeff=Et
            call epLightPreUserde(0, cTrack)
!           &&&&&&
         endif
      elseif( Move%proc .eq. 'nuci' ) then   
         call epmuNsmpE(Media(MediaNo), Move%Track%p%fm%p(4),
     *                   Et)
         cTrack%p%fm%p(4) = cTrack%p%fm%p(4) - Et  ! muon
         call epe2p(cTrack)         ! adjust momentum
         call eppush(cTrack)        ! push
         if(Et .gt. 152.d-3) then
            if(Media(MediaNo)%mu%MuNI .eq. 3 ) then
!             generate gamma-N interaction; employ gamma interaction
!             routine 
               cTrack%p%fm%p(4) = Et
               call cmkptc( kphoton, 0, 0, cTrack%p)
               call epe2p(cTrack) ! adjust momentum
!
               call ep2cosPtcl( cTrack%p )
!                fix target A
!                  better treatment for photons;  v9.16 (next 2)
               call cfixModel( cTrack%p )  ! reset model since E<mu's
               call cgetPhotoPxs(ActiveMdl2, cTrack%p, 
     *              xmedia(mediumNo), xs, mfp)
               call cfixTarget(xmedia(mediumNo))
               call epcpTargetInfo
               call ep2cosCond2   

               call cphotop     ! Cosmos function
               call eppushPtcl(cTrack) ! use pos. info from this ptcl
!                 &&&&&9.06
            elseif( Media(MediaNo)%mu%MuNI .eq. 2) then
               Move%dE=Et
               Move%dEeff=Et
               call epLightPreUserde(0, cTrack)
!           &&&&&&
            endif
         endif
      endif   
      end
      
