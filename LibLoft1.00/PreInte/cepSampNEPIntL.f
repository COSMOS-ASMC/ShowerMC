!      from epcosIntF
!     subroutine epsmpNEPIntL(media)
      subroutine cepSampNEPIntL(media, pj)
!       now don't worry below
!             xmedia=>media is to avoid name                          
!       collision of media  in modXsecMedia and                       
!       media argument in the subroutine def. 
!     !      use modXsecMedia, xmedia=>media
      use modEMcontrol
      implicit none
#include "Zglobalc.h"      
#include "Zmedia.h"
#include "Zcode.h"
#include "Zptcl.h"
#include "Zevhnv.h"
#include  "Zevhnp.h"



!#include  "Ztrackp.h"
!#include  "Ztrack.h"
!#include  "Ztrackv.h"
!#include  "Zheavyp.h"
!     #include  "Zelemagp.h"

!#include  "Zcmuint.h"
!     **************************************************
!
      type(epmedia),intent(in):: media
    !     type(xsmedia),intent(inout):: media
      type(ptcl),intent(in):: pj


      

      real*8 mfp,  xs, mass

      real*8 collkgram, u, length, path, prob

!     call cdecayLeng(TrackBefMove, length)


      call cdecayLeng(pj, length)

      if(length .ne. Infty) then
         call csetIntInf(length, .true., 'decay')
      endif

!     if( TrackBefMove%p%code /= kmuon)  then
      if( length > 0 )   then  ! =0; stopping mu, pi,  etc
         if( pj%code /= kmuon)  then
!         call epgetxs(ActiveMdl2, TrackBefMove.p, media, xs, mfp)
!     call cGetXsec(ActiveMdl2, TrackBefMove%p, media, xs, mfp)
            call cGetXsec(ActiveMdl2, pj, media, xs, mfp)         
            if(xs == smallxs ) then
               collkgram = Infty
            elseif(xs == largexs) then
               collkgram = 0.
            else
               call rndc(u)
               collkgram = -mfp*log(u)
            endif
            call csetIntInf(collkgram, .false., 'coll')
!!!!!!!!!!
!            write(0,*) ' colllkgram ',collkgram, ' mfp=', mfp
!!!!!!!!!!!!            
         else
         ! muon
!     add muon pair, brems, n.i if requested.
                     ! pair
            if( media%mu%MuPr .ge. 2 .and.
     *           pj%fm%p(4) .gt.  media%cnst%muPrEmin ) then
               call epmuPrsmpP(media,  pj%fm%p(4), prob, path)
!     path is in r.l. convert  it to kg/m2
               path = path * media%X0kg
               call csetIntInf(path, .false., 'pair')
            endif
                  ! mu Brems
            if(media%mu%MuBr  .ge. 2 .and.
     *           pj%fm%p(4) .gt.  media%cnst%muBrEmin ) then
               call epmuBrsmpP(media,  pj%fm%p(4), prob, path)
!     path is in r.l. convert  it to kg/m2
               path = path*media%X0kg
               call csetIntInf(path, .false., 'brem')
            endif

            if( media%mu%MuNI  .ge. 2 .and.
     *           pj%fm%p(4) .gt.  media%cnst%muNEmin ) then
               call epmuNsmpP( media,  pj%fm%p(4), prob, path)
!               path is in r.l. convert  it to kg/m2
               path = path*media%X0kg
               call csetIntInf(path, .false., 'nuci')
            endif
         endif
               
         if( pj%charge /= 0 .and. Knockon) then
!                add knock on
!     call epKnockp(Media(MediaNo), cTrack%p, prob, path)
            call epKnockp(media, pj, prob, path)
            path = path * media%X0kg
            call csetIntInf(path, .false., 'knoc') ! cosmos
!!!!!!!
!            write(0,*) ' knock path kg/m2 ', path
!!!!!!!!!            
         endif
      endif
      end
!     !      Next 2 are originally inside epcosIntF.f
!      
!     The same save and querry are possible using epUI in Epics
!     or chookNEPInt
!      
      subroutine epsaveFirstCol
!        This is called only when the first nuclear interaction
!     or decay to save the collision info.
!      use modColInfo  -->Zpwork.h ( which needs Zptcl.h)
      implicit none
#include "ZepMaxdef.h"
#include "Zptcl.h"
#include "Zpwork.h"
!      
!       Zptcl.h is inside modColInfo
!
      type(ptcl)::  ptclSave(EPMAX_STACK)
      integer nptcls
      common /stackSavec/ ptclSave, nptcls
      
      integer i
      if( Nproduced >  EPMAX_STACK) then
         write(0,*) "Warning, first interacton generated: "
         write(0,*) "# of particls =",Nproduced, " > ",
     *               EPMAX_STACK, " so "
         write(0,*) " only first ", EPMAX_STACK, " ptcls are saved"
         write(0,*) 
     *     " for later retrieving (by calling epqFirstColPtcls)."
         write(0,*) " Simulation itself is not affected. "
         write(0,*)
     *     " If you need all the generated particle information "
         write(0,*) " please consider the use of epUI interface."
      endif
      nptcls = min(Nproduced, EPMAX_STACK)

      do i = 1,  nptcls
         ptclSave(i) = Pwork(i)
      enddo
      end
      subroutine epqFirstColPtcls(ptcls, n, m)
        implicit none
#include  "Zptcl.h"
#include "ZepMaxdef.h"
       type(ptcl)::  ptclSave(EPMAX_STACK)
      integer nptcls
      common /stackSavec/ ptclSave, nptcls
      integer n, m
       type(ptcl)::  ptcls(n)
      save
      if(n .lt. nptcls) then
         write(0,*)
     *    " n must be > ", nptcls, " in epqFirstColPtcls"
         stop  9999
      endif
      do m = 1, nptcls
         ptcls(m) = ptclSave(m)
      enddo
      m = nptcls
      end
