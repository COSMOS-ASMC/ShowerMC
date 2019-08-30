      subroutine cmuInte(pj, media)   
!        pair creation,  brems or hadron producion by muon  
      use modSetIntInf 
      use modEMcontrol
      use modColInfo
      implicit none
#include  "Zmedia.h"
!#include  "Ztrack.h"
!#include  "Ztrackv.h"
!     #include  "Ztrackp.h"
#include  "Zptcl.h"
#include  "Zpwork.h"
#include  "Zcode.h"
      type(ptcl),intent(inout):: pj  ! projectile muon
      type(epmedia),intent(in):: media
      
      real*8 Et, Ee, u
      character*80 msg
      integer:: chg
      
      type(ptcl):: pair(2), gamma
!      aTrack=MovedTrack

!     if( IntInfArray(ProcessNo)%process == 'mupair' ) then
      if( IntInfArray(ProcessNo)%process == 'pair' ) then
!           call cmuPrsmpE(MovedTrack%p%fm%p(4),  Et)
         call epmuPrsmpE( media, pj%fm%p(4),  Et)
         Nproduced = Nproduced + 1
         Pwork(Nproduced) = pj
         Pwork(Nproduced)%fm%p(4) = pj%fm%p(4) - Et ! muon; neglect angle
!          adjust momentum without chaning the directon
         call ce2pp( Pwork(Nproduced) )  ! not ce2p(..)
         if( media%mu%MuPr >=  2 ) then
!     if(MuPr .eq. 3 .or. Eabsorb(1) .ne. 0 ) then
!           Let's    almost always  generate pair.
!           make parir elecron from  a gamma of energy Et
            pair(1)= pj         ! copy parent info.
!             Is it better  to consider LPM ?  currently LPMworks=.false.
            LPMworks = .false.  
!            sammple  electron energy (higher one) ; use g
            call epPrSampE(media, Et, pair(1)%fm%p(4) ) 
            call rndc(u)
            if(u < 0.5d0 ) then
               chg=-1
            else
               chg = 1
            endif
            call cmkptc(kelec, -chg, chg, pair(1))
            pair(2)%fm%p(4) = Et -  pair(1)%fm%p(4)
            call cmkptc(kelec, chg, -chg, pair(2))

            Nproduced = Nproduced + 1
!              adjust momentum
            call ce2pp(pair(1))
            call ce2pp(pair(2))
            Pwork( Nproduced ) = pair(1)
            Nproduced = Nproduced + 1
            Pwork( Nproduced ) = pair(2)
         endif
!     elseif(IntInfArray(ProcessNo)%process ==  'mubrem' ) then
      elseif(IntInfArray(ProcessNo)%process ==  'brem' ) then   
!     call cmuBrsmpE(MovedTrack%p%fm%p(4),  Et)
!
         call epmuBrsmpE(media, pj%fm%p(4),  Et)
         Nproduced = Nproduced + 1
         Pwork(Nproduced) = pj
         Pwork(Nproduced)%fm%p(4) = pj%fm%p(4) - Et ! muon; neglect angle
         call ce2pp( Pwork(Nproduced) )
         if( media%mu%MuBr >=  2 ) then
!         if(MuBr .eq. 3 .or. Eabsorb(1) .ne. 0 ) then
            Nproduced = Nproduced + 1
            Pwork(Nproduced) = pj
!            make it gamma
            call cmkptc(kphoton, kcasg, 0, Pwork(Nproduced) )
            call ce2pp( Pwork(Nproduced) )
         endif
!      elseif(IntInfArray(ProcessNo)%process ==  'munuci' ) then
      elseif(IntInfArray(ProcessNo)%process ==  'nuci' ) then 
         call epmuNsmpE(media, pj%fm%p(4),  Et)
         pj%fm%p(4) = pj%fm%p(4) - Et ! muon
         call ce2pp( pj )
         Nproduced = Nproduced + 1
         Pwork(Nproduced) =  pj ! muon
         if(Et >  153.d-3 ) then
            if(  media%mu%MuNi >=  2 ) then
!            if(  MuNI .eq. 3  .or. Eabsorb(1) .ne. 0 )  then 
!     generate gamma-N interaction;
!     employ gamma interaction  routine
               gamma = pj 
               call cmkptc(kphoton, 0, 0, gamma)
               gamma%fm%p(4) = Et
               call ce2pp( gamma )
               call cfixTargetMuNI( media ) ! for mu n.i, we must fix
                           ! target nuceus
               call cfixModel(gamma) ! fix int. model
!     call cphotop( pj )  !
               call cphotop( gamma ) ! 
            endif
         endif
      else
         write(msg, *) ' in cinteMuon: process=',
     *        IntInfArray(ProcessNo)%process,
     *        ' undef. ProcessNo=',ProcessNo
         call cerrorMsg(msg, 0)
      endif   
      end
