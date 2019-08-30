!     in the last part, knock-on by a non electron charged particle (epNEPKnock)
!     eixts, tought it is not a hadronic iteraction.abos
      
      subroutine cinteNuc(pj)
!!!   use modXsecMedia  !! don't use now
      use modColInfo
      use modSetIntInf
      implicit none
#include  "Zptcl.h"
#include "Zpwork.h"       
!     #include  "Ztrack.h"
!     #include  "Ztrackv.h"
      type(ptcl),intent(in):: pj

!          cinteNuc.  Treat Nucelon interactions with air target.
      character*70 msg
      if(IntInfArray(ProcessNo)%process .eq. 'coll') then
         call chAcol(pj, TargetNucleonNo, TargetProtonNo, TargetXs, 
     *        Pwork, Nproduced)
!        Now  momenum in Pwork is defined in the same system as
!           pj         
      else
         write(msg, *) ' in cinteNuc: process=',
     *                IntInfArray(ProcessNo)%process,
     *               ' undef. ProcessNo=',ProcessNo
         call cerrorMsg(msg,0)
      endif
      end
!     ********************
      subroutine cintePion(pj)
      use modColInfo
      use modSetIntInf
      use modMuNuControl
!     ********************
!      use modXsecMedia
      implicit none
#include  "Zptcl.h"
#include  "Zpwork.h"
!     #include  "Ztrackv.h"
! for   MuonPolarization adn IncMuonPolari,  next is needed
#include  "Ztrackp.h" 
      type(ptcl),intent(in):: pj
!
      character*70 msg

      if(IntInfArray(ProcessNo)%process .eq. 'coll') then
         call chAcol(pj, TargetNucleonNo, TargetProtonNo, TargetXs,
     *     Pwork, Nproduced)

      elseif(IntInfArray(ProcessNo)%process .eq. 'decay') then
         if(pj%charge .eq. 0) then
            call cpi0Decay(pj, Pwork, Nproduced)
         else
            call cpiMuDecay(pj, IncMuonPolari, 
     *           Pwork, Nproduced, MuonPolarization)
         endif
      else
         write(msg, *) ' in cintePion: process=',
     *           IntInfArray(ProcessNo)%process,
     *               ' undef. ProcessNo=',ProcessNo
         write(0,*) ' Et=', pj%fm%p(4), ' chg=',
     *                pj%charge
         call cerrorMsg(msg,0)
      endif

      end
!     ********************
      subroutine cinteKaon(pj)
!     ********************
!!!!  use modXsecMedia
      use modColInfo
      use modSetIntInf
      use modMuNuControl
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"
      
!----      include 'Ztrack.h'
!#include  "Ztrack.h"
!----      include 'Ztrackv.h'
!#include  "Ztrackv.h"
!----      include 'Ztrackp.h'
#include  "Ztrackp.h"
!
      type(ptcl),intent(in):: pj      
      character*70 msg

      if(IntInfArray(ProcessNo)%process .eq. 'coll') then
         call chAcol(pj, TargetNucleonNo, TargetProtonNo, TargetXs,
     *     Pwork, Nproduced)
      elseif(IntInfArray(ProcessNo)%process .eq. 'decay') then
         call ckaonDecay(pj, IncMuonPolari, 
     *           Pwork, Nproduced, MuonPolarization)
      else
         write(msg, *) ' in cinteKaon: process=',
     *           IntInfArray(ProcessNo)%process,
     *               ' undef. ProcessNo=',ProcessNo
         call cerrorMsg(msg, 0)
      endif
      end
!     ********************
      subroutine cinteDmes(pj)
!     ********************
!!!   use modXsecMedia
      use modColInfo
      use modSetIntInf
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"
#include  "Zcode.h"
!----      include 'Ztrack.h'

!#include  "Ztrack.h"
!----      include 'Ztrackv.h'
!#include  "Ztrackv.h"

      type(ptcl),intent(inout):: pj
!
      character*70 msg
      integer icg
      icg =  pj%charge
      if(IntInfArray(ProcessNo)%process .eq. 'decay') then
         call cdDecay(pj,  Pwork, Nproduced)
      elseif(IntInfArray(ProcessNo)%process .eq. 'coll') then   
!            use kaon 
         if(icg .eq. 0) then
            call cmkptc(kkaon, k0l, 0, pj)
         else
            call cmkptc(kkaon, 0, icg, pj)
         endif
         call ce2pp(pj)         ! force to change incident to kaon
         
         call chAcol(pj, TargetNucleonNo, TargetProtonNo, TargetXs,
     *     Pwork, Nproduced)
      else
         write(msg, *) ' in cinteDmeson: process=',
     *           IntInfArray(ProcessNo)%process,
     *               ' undef. ProcessNo=',ProcessNo
         call cerrorMsg(msg, 0)
      endif
      end
!     ********************
      subroutine cinteHeavy(pj)
!!!   use modXsecMedia
      use modColInfo
      use modSetIntInf
!     ********************
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"
!#include  "Zcode.h"      
!#include  "Ztrack.h"
!#include  "Ztrackv.h"
#include  "Ztrackp.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"

      type(ptcl),intent(in):: pj
!
      character*70 msg

      if(IntInfArray(ProcessNo)%process .eq. 'coll') then
         call cheavyInt(pj, TargetNucleonNo, TargetProtonNo,
     *         TargetXs,   Pwork, Nproduced)
      else
         write(msg, *) ' in cinteHeavy: process=',
     *           IntInfArray(ProcessNo)%process,
     *               ' undef. ProcessNo=',ProcessNo
         call cerrorMsg(msg, 0)
      endif
      end
!     ********************
      subroutine cinteMuon(pj)
!!!   use modXsecMedia
      use modColInfo
      use modSetIntInf
      use modEMcontrol
      use modMuNuControl
!     ********************
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"
!#include  "Ztrack.h"
!#include  "Ztrackv.h"

#include  "Zlife.h"
#include  "Ztrackp.h"
#include  "ZmediaLoft.h"
      
      type(ptcl),intent(in):: pj
      
      real*8  capr, t0cap, u1, u2

!
      if(IntInfArray(ProcessNo)%process .eq. 'decay') then
!           see if negative muon is being captured
         if(pj%charge .eq. -1) then
            if(pj%fm%p(4) .le. pj%mass*1.001) then
!                   assume stopping muon
               call cgetCaprate( Media(MediaNo) )  ! v7.640
!               call cmucap(TargetNucleonNo, TargetProtonNo, capr)
               capr = Media(MediaNo)%colXs
               t0cap = 1./capr
               call rndc(u1)
               call rndc(u2)
               if( - log(u1)*t0mu .gt. - log(u2)*t0cap) then
!                       capture
                  IntInfArray(ProcessNo)%process="capt"
                  call cfixTarget( Media(MediaNo) )
                  call ccapnu(TargetNucleonNo, TargetProtonNo,
     *               Pwork, Nproduced)
               else
!!!!              MuonPolarization = 0.   ! stopping mu has no pol.
!                                          keep the value ;it's better.
                  call cmuNeuDcy(pj, MuonPolarization,
     *            Eabsorb,   Pwork, Nproduced)
!     *            Eabsorb(1),   Pwork, Nproduced)
                  
               endif
            else
               call cmuNeuDcy(pj, MuonPolarization,
     *            Eabsorb,   Pwork, Nproduced)
!
            endif
         else
            call cmuNeuDcy(pj, MuonPolarization,
     *            Eabsorb,   Pwork, Nproduced)
!     *            Eabsorb(1),   Pwork, Nproduced)
         endif
      else
          !  pair cre , brems, nuclear interaction
         call cmuInte(pj, Media(MediaNo)) 
      endif
      end
      subroutine csetMuonPol(val)
!         This may be used when you set muon as a primary for which
!        polarization is not fixed.
!
      use modMuNucontrol
      implicit none
!#include  "Ztrack.h"
!#include "Ztrackv.h"
      real*8 val
      MuonPolarization =min(1.d0, max(val, -1.d0 ))
      end
!     ********************
      subroutine cintennb(pj)
!     ********************may be used if 1ry is nn~
      use modColInfo
      use modSetIntInf
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"
      
!----      include 'Ztrack.h'
!#include  "Ztrack.h"
!----      include 'Ztrackv.h'
!#include  "Ztrackv.h"

      type(ptcl),intent(in):: pj !
      character*70 msg

      if(IntInfArray(ProcessNo)%process .eq. 'decay') then
         call cnnbdc(pj,  Pwork, Nproduced)
      else
         write(msg, *) ' in cintennb: process=',
     *           IntInfArray(ProcessNo)%process,
     *               ' undef. ProcessNo=',ProcessNo
         call  cerrorMsg(msg, 0)
      endif
      end
!     ********************
      subroutine cinteddb(pj)
!     ********************may be used if 1ry is dd~
      use modColInfo
      use modSetIntInf
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"
      
!----      include 'Ztrack.h'
!#include  "Ztrack.h"
!----      include 'Ztrackv.h'
!#include  "Ztrackv.h"
!
      type(ptcl),intent(in):: pj
      character*70  msg

      if(IntInfArray(ProcessNo)%process .eq. 'decay') then
         call cddbdc(pj,  Pwork, Nproduced)
      else
         write(msg, *) ' in cintennb: process=',
     *           IntInfArray(ProcessNo)%process,
     *               ' undef. ProcessNo=',ProcessNo
         call cerrorMsg(msg, 0)
      endif
      end
!     ********************
      subroutine cinteSigma(pj)
!     *****************
      use modColInfo
      use modSetIntInf
!!      use modXsecMedia
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"

!#include  "Ztrack.h"
!#include  "Ztrackv.h"
!
      type(ptcl),intent(in):: pj

      if(IntInfArray(ProcessNo)%process .eq. 'decay') then
         call csigmaDecay(pj,  Pwork, Nproduced)
      else
         call chAcol(pj, TargetNucleonNo, TargetProtonNo, TargetXs,
     *     Pwork, Nproduced)
      endif
      end
!     ********************
      subroutine cinteLambda(pj)
!     ****************
!     !      use modXsecMedia
      use modColInfo
      use modSetIntInf
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"

!#include  "Ztrack.h"
!#include  "Ztrackv.h"
!
!
      type(ptcl),intent(in):: pj
      if(IntInfArray(ProcessNo)%process .eq. 'decay') then
         call clambdaDcy(pj,  Pwork, Nproduced)
      else
         call chAcol(pj, TargetNucleonNo, TargetProtonNo, TargetXs,
     *     Pwork, Nproduced)
      endif
      end
!     ********************
      subroutine cinteGzai(pj)
!     *****************
!     use modXsecMedia
      use modColInfo
      use modSetIntInf
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"

!#include  "Ztrack.h"
!#include  "Ztrackv.h"
!
      type(ptcl),intent(in):: pj
      if(IntInfArray(ProcessNo)%process .eq. 'decay') then
         call cgzaiDecay(pj,  Pwork, Nproduced)
      else
         call chAcol(pj, TargetNucleonNo, TargetProtonNo, TargetXs,
     *     Pwork, Nproduced)
      endif
      end
!     ********************
      subroutine cinteBomega(pj)
!     *****************
!     !      use modXsecMedia
      use modColInfo
      use modSetIntInf
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"

!#include  "Ztrack.h"
!#include  "Ztrackv.h"
!
      type(ptcl),intent(in):: pj
      if(IntInfArray(ProcessNo)%process .eq. 'decay') then
         call cbomegaDcy(pj,  Pwork, Nproduced)
      else
         call chAcol(pj, TargetNucleonNo, TargetProtonNo, TargetXs,
     *     Pwork, Nproduced)
      endif
      end
!     *******************
      subroutine cinteEta(pj)
!     *****************
!     !      use modXsecMedia
      use modColInfo
      use modSetIntInf
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"

!#include  "Ztrack.h"
!#include  "Ztrackv.h"
!
      type(ptcl),intent(in):: pj
      if(IntInfArray(ProcessNo)%process .eq. 'decay') then
         call cetaDecay(pj,  Pwork, Nproduced)
      else
         call chAcol(pj, TargetNucleonNo, TargetProtonNo, TargetXs,
     *     Pwork, Nproduced)
      endif
      end
!     ********************
      subroutine cinteLambdac(pj)
!     ********************
!     !      use modXsecMedia
      use modColInfo
      use modSetIntInf
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"
      
!#include  "Ztrack.h"
!#include  "Ztrackv.h"
#include  "Zcode.h"
!
      character*80 msg

      type(ptcl),intent(inout):: pj
      if(IntInfArray(ProcessNo)%process .eq. 'decay') then
         call clambdacDcy(pj,  Pwork, Nproduced)
      elseif(IntInfArray(ProcessNo)%process .eq. 'coll') then   
!            use proton or p-bar
         if( pj%charge > 0 ) then
            call cmkptc(knuc, -1, 1, pj)
         else
            call cmkptc(knuc, 1, -1, pj)
         endif
         call ce2pp(pj)      ! not ce2p;   force incident to nucleon. 
         call chAcol(pj, TargetNucleonNo, TargetProtonNo, TargetXs,
     *     Pwork, Nproduced)
      else
         write(msg, *) ' in cinteLambdac: process=',
     *           IntInfArray(ProcessNo)%process,
     *               ' undef. ProcessNo=',ProcessNo
         call cerrorMsg(msg, 0)
      endif
      end
!     ********************
      subroutine cinterho(pj)
!     *******************
      use modColInfo
      use modSetIntInf
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"
      
!#include  "Ztrack.h"
!#include  "Ztrackv.h"
!
      type(ptcl),intent(in):: pj      
      character*80 msg

      if(IntInfArray(ProcessNo)%process .eq. 'decay') then
         call crhodc(pj,  Pwork, Nproduced)
      elseif(IntInfArray(ProcessNo)%process .eq. 'coll') then   
         Nproduced = 0     ! neglect coll.
      else
         write(msg, *) ' in cinterho: process=',
     *           IntInfArray(ProcessNo)%process,
     *               ' undef. ProcessNo=',ProcessNo
         call cerrorMsg(msg, 0)
      endif
      end
!     ********************
      subroutine cinteomega(pj)
!     ****************
      use modColInfo
      use modSetIntInf
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"
      
!#include  "Ztrack.h"
!#include  "Ztrackv.h"
!
      type(ptcl),intent(in):: pj
      character*80 msg

      if(IntInfArray(ProcessNo)%process .eq. 'decay') then
         call comgdc(pj,  Pwork, Nproduced)
      elseif(IntInfArray(ProcessNo)%process .eq. 'coll') then   
         Nproduced = 0     ! neglect coll.
      else
         write(msg, *) ' in cinteomega: process=',
     *           IntInfArray(ProcessNo)%process,
     *               ' undef. ProcessNo=',ProcessNo
         call cerrorMsg(msg, 0)
      endif
      end
!     ********************
      subroutine cintephi(pj)
!     *******************
      use modColInfo
      use modSetIntInf
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"
      
!#include  "Ztrack.h"
!#include  "Ztrackv.h"
!
      type(ptcl),intent(in):: pj      
      character*80 msg

      if(IntInfArray(ProcessNo)%process .eq. 'decay') then
         call cphidc(pj,  Pwork, Nproduced)
      elseif(IntInfArray(ProcessNo)%process .eq. 'coll') then   
         Nproduced = 0     ! neglect coll.
      else
         write(msg, *) ' in cintephi: process=',
     *           IntInfArray(ProcessNo)%process,
     *               ' undef. ProcessNo=',ProcessNo
         call cerrorMsg(msg, 0)
      endif
      end
!     ******************
      subroutine cinteEtap(pj)
!     *****************
      use modColInfo
      use modSetIntInf
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"

      type(ptcl),intent(in):: pj
      
      if(IntInfArray(ProcessNo)%process .eq. 'decay') then
         call cetapDecay(pj,  Pwork, Nproduced)
      else
!           negelct col.
         Nproduced = 0
      endif
      end
!     ********************
      subroutine cinteds(pj)
!     *******************
      use modColInfo
      use modSetIntInf
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"
      type(ptcl),intent(in):: pj      
      character*80 msg

      if(IntInfArray(ProcessNo)%process .eq. 'decay') then
         call cDsDecay(pj, Pwork, Nproduced)
      elseif(IntInfArray(ProcessNo)%process .eq. 'coll') then
         Nproduced = 0          ! neglect coll.
      else
         write(msg, *) ' in cinteds: process=',
     *           IntInfArray(ProcessNo)%process,
     *               ' undef. ProcessNo=',ProcessNo
         call cerrorMsg(msg, 0)
      endif
      end
!     ********************
      subroutine cinteDelta(pj)
!     ********************
      use modColInfo
      use modSetIntInf
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"
      type(ptcl),intent(in):: pj      

      if(IntInfArray(ProcessNo)%process .eq. 'decay') then
         call cDeltaDecay(pj, Pwork, Nproduced)
      else
         write(0,*)
     *     ' Delta resonance collision should not happen '
         stop
      endif
      end
      
!     ********************
      subroutine cinteTau(pj)
!     *******************
      use modColInfo
      use modSetIntInf
      implicit none
#include "Zptcl.h"
#include "Zpwork.h"
      type(ptcl),intent(in):: pj      
!
      character*80 msg

      if(IntInfArray(ProcessNo)%process .eq. 'decay') then
         call ctauNeuDcy(pj,  Pwork, Nproduced)
      elseif(IntInfArray(ProcessNo)%process .eq. 'coll') then
         Nproduced = 0          ! neglect coll.
      else
         write(msg, *) ' in cintetau: process=',
     *           IntInfArray(ProcessNo)%process,
     *               ' undef. ProcessNo=',ProcessNo
         call cerrorMsg(msg, 0)
      endif
      end
      

      subroutine epNEPKnock( pj )
      use modColInfo
      implicit none
#include  "Zcode.h"
#include  "Zmass.h"
!#include  "Ztrack.h"
!     #include  "Ztrackv.h"
#include "Zptcl.h"
#include "Zpwork.h"
!
      type(ptcl),intent(in):: pj  ! type(ptcl) is in modColInfo

      
      real*8 e1, er, tmp, cos1, cosr,  cs, sn, sinr, sine
      real(8):: dc(3)
      real(8):: p

!                                surv  elec  surv   elec
!      call  cKnockea(aTrack%p,  e1,   er,  cos1,   cosr)
      call  epKnockea(pj,  e1,   er,  cos1,   cosr)
!          survival paticle angle negligible always so
!          you may even put dc = (0, 0, 1)
       tmp=1.d0-cos1*cos1
       if(tmp .lt. 0.d0) then
          tmp=0.d0
          cos1=1.d0
       endif
       sine=sqrt(tmp)

       call kcossn(cs, sn)
       !    - is for  opposit side 
       dc(1) = cs*sine 
       dc(2) = sn*sine
       dc(3) = cos1
       Nproduced  = Nproduced + 1
       Pwork(Nproduced) = pj    ! copy the incident   
       Pwork(Nproduced)%fm%p(4) = e1  ! set energy
       p =  e1**2- pj%mass**2
       if(p < 0.) then
          p = 0.
          Pwork(Nproduced)%fm%p(4) = pj%mass
       else
          p= sqrt(p)
       endif
       Pwork(Nproduced)%fm%p(1:3) = dc(1:3)*p
      
!            knock on electron
       tmp=1.d0-cosr*cosr

       if(tmp .lt. 0.d0) then
          tmp=0.d0
          cosr=1.d0
       endif
       sinr=sqrt(tmp)
!          azimuthally opsit
       dc(1) = -cs*sinr
       dc(2) = -sn*sinr
       dc(3) = cosr
       Nproduced  = Nproduced + 1
       call cmkptc(kelec,  0,  -1,  Pwork(Nproduced))    ! make e-
       Pwork(Nproduced)%fm%p(4) = er  ! set E
       p = er**2 - masele**2
       if(p < 0.) then   ! would not hpappen since er > some min E> masele
          p = 0.
       else
          p = sqrt(p)
       endif
       Pwork(Nproduced)%fm%p(1:3) = dc(1:3)*p

!            we need rotation here. 
       call crot3mom(pj, Pwork, Nproduced)
       end
