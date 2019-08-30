      subroutine cepSampEIntL(media, pj)
!     before calling this
!     rrho =actual(density)/media%rho in mod  must
!     be fixed.
!     magnetic brems sampling is not done here
!     it must be done by fixing mag vector and call
!      xxx    cepMagBremPath(...)   xxx ==>    epsyncp
!     subroutine  epproe
!     *************
      use modMCScontrol
      use modEMcontrol  
      implicit none
!
!     electron:     brems, knock-on, and anihilation
!          synchrotron  radiation should be considered outside
!
#include  "Zmedia.h"
#include  "Zptcl.h"
      
! #include  "ZepTrackv.h"
! #include  "ZepTrackp.h"

!
      type(epmedia),intent(in):: media
      type(ptcl),intent(in):: pj
!     real(8),intent(in):: rrho !  rho(actual)/media%rho
!         in modEMcontrol      

      
      real*8 E, prob,  tbrem, tknock, tanihi, t, dt, dl, dx
!      ,
!     *   syncmfp
      real(8)::u
!      logical:: posok, numok
!      type(epPos):: posw
      integer:: i
       
      E = pj%fm%p(4)
!     sample path for brems;
      call epBrSampP(media,  E,  prob, tbrem)
      call csetIntInf(tbrem, .false., 'brem')

      if(Knockon) then
         if(pj%charge .eq. -1) then
            call epmollerp(media, E, RecoilKEmin, prob, tknock)
         else
            call epbhabhap(media, E, RecoilKEmin, prob, tknock)
         endif

         call csetIntInf(tknock, .false., 'knoc')
         
!         if(tbrem .le. tknock) then
!            t = tbrem
!            Move%proc='brem'
!         else
!           t = tknock
!            Move%proc='knoc'
!         endif
!      else
!         t = tbrem
!         Move%proc='brem'
      endif
!     if(cTrack%p%charge .eq. 1 .and. E .lt. Eanihi) then
      if( pj%charge .eq. 1 .and. E .lt. Eanihi ) then
!     call epanihip(Media(MediaNo), E, prob, tanihi)
         call epanihip( media, E, prob, tanihi)
         call csetIntInf( tanihi, .false., 'anih' )
!         if(tanihi .lt. t) then
!            t = tanihi
!            Move%proc='anih'
!         endif
      endif
!!!!  next should be done outside ; EPICS; epMCShard is made and called
!     after cepSampEIntL
      
!!!      if( doNewMCS ) then
!!!         if( (media%name /= 'sp' ) .and.
!!!     *       (media%name /= 'hollow')  ) then
!!!!            call cfixMCSmodel( cTrack%p ) ! see energetically OK ?
!!!            call cfixMCSmodel( pj ) ! see energetically OK ?
!!!
!!!            if( ActiveMCS /= 'Mol') then
!!!
!!!               if( MCSzCond > 0 ) then
!!!!                see if z is in the specifed region; else -> Mol
!!!                  posok = .false.
!!!                  call epl2w(cTrack%cn, cTrack%pos, posw)
!!!                  do i = 1, MCSzCond
!!!                     if(real(MCSzRange(i)) <= posw%z ) then
!!!                        if( posw%z <= imag(MCSzRange(i)) ) then
!!!                           posok = .true.
!!!                           exit
!!!                        endif
!!!                     endif
!!!                  enddo
!!!               else
!!!                  ! no check about Z
!!!                  posok=.true.
!!!               endif
!!!
!!!               if( MCSnumCond > 0 ) then
!!!                  if( MCSdebug ) then
!!!                     call epl2w(cTrack%cn, cTrack%pos, posw)
!!!                  endif
!!!                  ! see if cn is in the ranage
!!!                  numok = .false.
!!!                  do i = 1, MCSnumCond
!!!                     if(MCSnumRangeMin(i) <=  cTrack%cn ) then
!!!                        if( cTrack%cn <= MCSnumRangeMax(i) ) then
!!!                           numok = .true.
!!!                           exit
!!!                        endif
!!!                     endif
!!!                  enddo
!!!               else
!!!                  ! no check about cn
!!!                  numok = .true.
!!!               endif
!!!
!!!
!!!               if( MCSandor == 'and' ) then
!!!                  if( posok .and. numok ) then
!!!                     if( MCSrevert ) then
!!!                        ActiveMCS = 'Mol'
!!!                     else
!!!                        !  keep current one
!!!                     endif
!!!                  else
!!!                     if( .not. MCSrevert ) then
!!!                        ActiveMCS = 'Mol'
!!!                     else
!!!!                      keep current
!!!                     endif
!!!                  endif
!!!                  
!!!               else             ! 'or'
!!!                  if( posok .or.  numok ) then
!!!                     if(  MCSrevert ) then
!!!                        ActiveMCS = 'Mol'
!!!                     else
!!!!                         keep current
!!!                     endif
!!!                  else
!!!                     if( .not. MCSrevert ) then
!!!                        ActiveMCS = 'Mol'
!!!                     else
!!!!                          keep current
!!!                     endif
!!!                  endif
!!!               endif
!!!               if( MCSdebug ) then 
!!!                  if( ActiveMCS == 'Mol') then
!!!                     write(0,'(a, 1p,3g14.4)') 'mol: ',  posw
!!!                  else
!!!                     write(0,'(a, 1p,3g14.4)') 'hin: ',  posw
!!!                  endif
!!!               endif
!!!            else
!!!               ActiveMCS = 'Mol'
!!!            endif
!!!         else
!!!            ActiveMCS = 'Mol'
!!!         endif
!!!      endif
!!!      
!!!      if( ActiveMCS /= 'Mol' ) then
!!!         call cfixMixedConst(MediaNo, int(cTrack%p%charge))
!!!      endif
!!!
!!!
!!!      if( ActiveMCS == 'El_hin') then
!!!             ! sample hard scattering  mfp g/cm2.
!!!         call cgetLamh( KEeV, mfpHardgr)  ! mfp of hard cs
!!!!  in g/cm2
!!!         call rndc(u)
!!!         lHardgr = -log(u)* mfpHardgr  ! g/cm2
!!!         lHardrl =lHardgr / Media(MediaNo)%X0g    ! rl
!!!          !  r.l
!!!          ! next is  cm length. this might  be used later for soft
!!!          ! mcs treatment.
!!!         lHardcm = lHardgr* Media(MediaNo)%gtocm
!!!     *        / Media(MediaNo)%rhoc
!!!         if( lHardrl < t ) then
!!!            Move%proc="hcs"  ! hard Coulomb scattring
!!!            t = lHardrl
!!!         endif
!!!      endif
!!!
!!!      Move%dt = t
!!!      Move%dx = Move%dt * Media(MediaNo)%X0g
!!!      Move%dl = Move%dx * Media(MediaNo)%gtocm /
!!!     *    Media(MediaNo)%rhoc
!!!
!!!!                     only if X0 > 10 km; may be 30 km is o.k
!!!      if(Sync .eq.  2 .and.
!!!     *    Media(MediaNo)%X0/Media(MediaNo)%rhoc .gt. 10.d5 ) then
!!!!          sample synchrotron emission path
!!!         call epsyncp(cTrack%p, Bfield, Upsilon, syncmfp, dl)
!!!         dx = dl / Media(MediaNo)%gtocm *
!!!     *    Media(MediaNo)%rhoc
!!!         dt = dx / Media(MediaNo)%X0g
!!!         if(dt .lt.  t) then
!!!            Move%dt = dt
!!!            Move%dx = dx
!!!            Move%dl = dl
!!!            Move%proc = 'sync'
!!!         endif
!!!      else
!!!      endif
      end
