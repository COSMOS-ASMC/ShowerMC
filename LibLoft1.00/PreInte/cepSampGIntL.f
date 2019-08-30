      subroutine cepSampGIntL(media, pj)
      use modEMcontrol
!      subroutine epprog      
!     1)   Before this call   ciniIntInf must be called which is in
!          csetIntInf.f  (clear counters)
!     Then, this may be called. Gamma ray processes
!     except for magnetic pair creation  are  sampled
!       --if energy is appropiate for the process.
!     Five  processes, i.e, pair creation, compton, Photo-electric
!     effect, coherent scatt.  and Photo-production of hadronsns are
!     considered.
!           Each  process info is stored in modIntInf      
!     ***   All lengths for various processes in this sub. 
!           are given in r.l  ***
!                 
!     2)After this, if needed,the magnetic pair creation routine should
!     be called separetately after fixing B field. 
!
!            
!
!      
!   3)   The minimum of sampled interactin lengths must
!     be selected.
!         
!      

!       use modXsecMedia, xmedia=>media, xelement=>element,
!     *     dontuse=>NoOfMedia
       implicit none
!
!#include  "ZepTrackv.h"
!#include  "ZepTrackp.h"
!     #include  "Zevhnv.h"
#include "Zmedia.h"
#include "Zptcl.h"
#include "Zevhnv.h"
       
      type(epmedia),intent(in):: media  ! current media
      type(ptcl),intent(in)::  pj  ! projectile partilce
!      real(8),intent(in):: rrho ! rho(actual)/media%rho.  now in modEMcontrol 
!     if EPICS, it is media%rhoc
!     if Cosmos,  it is den(of the current partilc)/media%rho      
!
      
      
      real*8 tcomp, tphot, tpair, tgp, tcoh,  t
      real*8 E, prob, xprob(5), txray(5)
      real*8 xs, mfp
      real*8  pairmfp, dl, dx, tmpair, u
      real*4  xsec(5)  !  coh, incoh,  P.E  1/(g/cm2)
!      real*8  Excom1, Excom2  !now in ZepTrackp.h
      integer icon
!            where xcom data is used. 
!       Excom1: compton/photo abs/coherent scat
!       Excom2: pair; default is  not use xcom
!              both must be < 100 GeV
!      data Excom1/1.d-3/, Excom2/1.d-3/
!      save Excom1, Excom2

!     EupperBndCS in ZepTrackp.h
!     Excoms1,2     //
!     !             note: sampled EM process path is in r.l.
!     !  For the process which cannot occur, neglect.
      E = pj%fm%p(4)
      if(E .le. EupperBndCS) then   ! 
         xprob(3)= 0.
         if( media%xcom%size .gt. 0 .and.  E .lt. Excom1 ) then
!              
            call epXrayp(media, E, 1,  3,  xprob, txray)
!            tcomp=txray(2)
!            tcoh = txray(1)
!            tphot = txray(3)
            call csetIntInf(txray(1), .false., 'coh')    !coherent scat
            call csetIntInf(txray(2), .false., 'comp')   !compton 
            call csetintInf(txray(3), .false., 'phot')    ! photo electric
         else
!            tcoh = 1.e35
            call epcompp(media, E, prob, tcomp)
            call csetIntInf(tcomp, .false., 'comp')   
         endif
         if(Photo) then    ! p.e effect
            if(xprob(3)  .eq. 0.) then
               call epphotoEp(media,  E, prob, tphot) ! v8.0
               call csetIntInf(tphot, .false., 'phot')
            endif
!         else
!            tphot = 1.e35
         endif
!      else
!         tcomp =1.e35
!         tphot= 1.e35
!         tcoh = 1.e35
      endif
      if(E .gt. ElowerBndPair) then
         if( media%xcom%size .gt. 0 .and.
     *        E .lt. Excom2 ) then
            call epXrayp(media, E, 4, 5,  xprob, txray)
            prob= xprob(4)+xprob(5)
            call rndc(u)
            tpair = -log(u)/prob
         else
            call epPrSampP(media, E, prob, tpair)
         endif

         call csetIntInf(tpair, .false., 'pair')
!      else
!         tpair=1.e35
      endif
      if(IncGp > 0 .and. E .gt. 153.d-3) then ! > 0.152431798977028
!         call ep2cosCond
!     call cfixModel( cTrack%p )
         call cfixModel( pj )
!     call cgetPhotoPxs(ActiveMdl2, cTrack%p, xmedia(mediumNo),
         call cgetPhotoPxs(ActiveMdl2, pj, media,  xs, mfp)

!     prob = xs*Media(MediaNo)%mbtoPX0 ! prob/r.l
         prob = xs * media%mbtoPX0 ! prob/r.l         
         call rndc(u)
         tgp = -log(u)/prob     ! sampled path int r.l
         call csetIntInf(tgp, .false., 'photop') 
!      else
!         tgp=1.e35
      endif

!!!      if(MagPair .eq. 1) then
!!!         call epmpairp(cTrack%p, Bfield, Xai, pairmfp, dl)
!!!         dx = dl / Media(MediaNo)%gtocm *
!!!     *     Media(MediaNo)%rhoc
!!!         tmpair = dx / Media(MediaNo)%X0g
!!!      else
!!!         tmpair = 1.e35
!!!      endif
!!!      t=min(tpair, tcomp, tphot, tgp, tmpair, tcoh)
!!!      if(t .eq. tpair) then
!!!         Move%proc='pair'
!!!      elseif(t .eq. tcomp) then
!!!         Move%proc='comp'
!!!      elseif(t .eq. tphot) then
!!!         Move%proc='phot'
!!!      elseif(t .eq. tcoh) then
!!!         Move%proc='coh'
!!!      elseif(t .eq. tgp) then
!!!         Move%proc='photop'
!!!      else
!!!         Move%proc='mpair'
!!!      endif
!!!      Move%dt = t   ! in r.l
!!!      Move%dx = Move%dt * Media(MediaNo)%X0g
!!!      Move%dl = Move%dx * Media(MediaNo)%gtocm /
!!!     *     Media(MediaNo)%rhoc
      end
