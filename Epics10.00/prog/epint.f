       subroutine epint(icon)
       use modIntInfo
!       use modXsecMedia, xmedia=>media, xelement=>element,
!     *     dontuse=>NoOfMedia
       use modV1ry
       implicit none
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"
#include  "Zcode.h"
#include  "Zevhnv.h"
      integer icon  ! output.  always 1

      character*100 msg
      integer k,  ia, iz
      real*8 xs
      integer inela
       type(epPos):: temppos

!
!/////////////
!      logical show
!      common /showshow/show
!///////////////

      k = cTrack%p%code
!         almost dummy setting 
!      ia =Media(MediaNo).A
!      iz =Media(MediaNo).Z

#if defined (INTINFO)
      kintInfo = min(k, maxcodeForInt)
      if( codeAforInt(kintInfo) == 0 ) then
          !  we have to inform int info to epUI
          !  get current stack pos
         call epqstn(IntInfo1) ! product is put from posIntInfo1+1
         IntInfo1 = IntInfo1 + 1
      endif            
#endif


      if(k .eq. kphoton) then
         if(Move%proc .eq. 'comp') then
            call epcmpt
         elseif(Move%proc .eq. 'pair') then
            call eppair
         elseif(Move%proc .eq. 'phot') then
            call epphot
         elseif(Move%proc .eq. 'coh') then
            call epcoher
         elseif(Move%proc .eq. 'photop') then
!                     need not now
!            call ep2cosPtcl( cTrack%p )
!            call cfixTarget(xmedia(mediumNo)) ! 
!            call epcpTargetInfo  ! cp target info 
!     call ep2cosCond2  !
!                       
            call cphotop( cTrack%p )        ! Cosmos function
            call eppushPtcl(cTrack)  ! use pos. info from this ptcl
         elseif(Move%proc .eq. 'mpair') then
            call epmpair
         else
            write(msg,
     *       '('' proccess='',a4,'' for gamma undefined'')')
     *       Move%proc
            call cerrorMsg(msg,0)
         endif
      elseif(k .eq. kelec) then

         if(Move%proc .eq. 'brem') then
            call epbrem
         elseif(Move%proc .eq. 'knoc') then
            call epknoc
         elseif(Move%proc .eq. 'hcs' ) then
!           it has been done in  epdoMixedMCS2; 
!            so simply puth current track in stack
            call eppush(cTrack)
         elseif(Move%proc .eq. 'anih') then
            call epanih
         elseif(Move%proc .eq. 'sync') then
            call epsync
         else
            write(msg,
     *         '('' process='',a4, '' for e is undef.'')') Move%proc
            call cerrorMsg(msg, 0)
         endif
       !>>>>>>>>>>>>>>>>>>>>>>light
      elseif( k == klight) then
         if(Move%proc == "pe" ) then
                  ! photo electron generation at sensor
            call epLightAtSensor
         elseif( Move%proc == "rayl" ) then
                  ! Rayleigh scattering; use Xray region fomulat since
                  ! (1+cos^2)dcos
            call epcoher
         elseif( Move%proc == "absorb" ) then
                  ! absorbed. nothing to do;  not push any thing
         elseif( Move%proc == "wls" ) then
                   ! wave length shift
            call epLightPreWLS
         else
            write(0,*) ' light interacion=', Move%proc
            write(0,*) ' not defined '
            stop
         endif
        !<<<<<<<<<<<<<<<<<<<<<<
      else

         if(Move%proc .eq. 'knoc') then

            call epNEPknoc

         elseif(cTrack%p%code .eq. kmuon ) then
            if( Move%proc == 'decay') then
!               call ep2cosPtcl( cTrack%p ) ! set
!               call ep2coscond2            ! Move.Track
               call cinteMuon( cTrack%p )  ! capture is also treated
                           !   don't use nuc. interaction there
            elseif(Move%proc .eq.'pair') then
               call epmuInte
            elseif(Move%proc .eq.'brem') then
               call epmuInte
            elseif(Move%proc .eq.'nuci') then
               call epmuInte
            endif
         else  
            ! hadronic ptcls
!!?            call ep2cosPtcl( cTrack%p )
!            if(  k == kmuon .and. cTrack.p.charge == -1  .and.
!     *           cTrack.p.fm.p(4) <= cTrack.p.mass*1.001) then
!                        1.001 must be the same as in cinteMuon
!                 for decay of stopping mu-, 
!                 we must fix target for capture case; really captured
!                 or decay  is determined in cinteMuon; 

            call cfixTarget(xmedia(mediumNo))
!!?            call epcpTargetInfo ! cp target info 
!!?            call ep2cosCond2

            call cinteNEP( cTrack%p  )       ! cosmos eppp
            call eppushPtcl( cTrack )
         endif
      endif


      if( FirstC ) then
         call epSeeIf1stInt     ! see if really first interaction
!          if so, FirstC becomes F
         if(.not. FirstC) then  ! this is really 1st int.
            V1ry = 0
!            the interaction was recognized as the first one
!             convet vector into world coord.
            call epl2wTrack(cTrack, FirstIntTrack)
            FirstInt = FirstIntTrack%pos ! cp pos info. for compat
!            call epl2w(cTrack.cn, FirstInt, temppos)
!            FirstInt = temppos
            FirstMedia = Media(MediaNo)
!            FirstIntTrack.pos = temppos
            firstCn = cTrack%cn
            if(Light == 21 ) then
               call epLightIOwrite1stCol
            endif
         endif
      endif


#if defined (INTINFO)
      if( codeAforInt(kintInfo) == 0 ) then
         call epqstn(IntInfo2)
         codeAforInt(kintInfo) = kintInfo
         call epUI(codeAforInt(kintInfo), IntInfo1, IntInfo2)
      endif
#endif

      icon = 1
      end
      

      subroutine epcpTargetInfo      ! not needed now ?
!       copy target info in xmedia to Epics Media
!      use modXsecMedia, only: TargetNucleonNo, TargetProtonNo,
!     *  TargetXs, colElemNo
      use modColMediaInfo
      implicit none
#include  "ZmediaLoft.h"      
!#include  "ZepTrackv.h"
!#include  "ZepTrackp.h"
!#include  "Zcode.h"
!#include  "Zevhnv.h"
      Media(MediaNo)%colA = TargetNucleonNo
      Media(MediaNo)%colZ = TargetProtonNo
      Media(MediaNo)%colXs = TargetXs
      Media(MediaNo)%colElem = colElemNo
      end 
!     ******************
      subroutine epcoher
!        coherent scattering
!           since coherent scattering is effective at
!           low energies where angular distribution can be
!           approximated by (1+cos^2) dcos, we simply use this
       implicit none
#include  "ZepTrackv.h"
!
       type(epDirec)::  w

       real*8  eg, tmp, cosg
       real*8  cs, sn, sing

!             sample scattering angle from (1+cos^2)dcos
       call ksampRSA(cosg)
       tmp=1.d0-cosg*cosg
       sing = sqrt(tmp)
       call kcossn(cs,sn)

       w%x = cs*sing
       w%y = sn*sing
       w%z = cosg
       call eptransVect(cTrack%w,  w, cTrack%w)
!        energy unchaged;  
       call epe2p(cTrack)
       call eppush(cTrack)
       end
!     ******************
      subroutine epcmpt
       implicit none
#include  "ZepTrackv.h"
#include  "Zcode.h"

!
       type(epTrack)::  electron
       type(epDirec)::  w

       real*8 e1, eg, tmp, cosg, cose
       real*8 sine, cs, sn, sing

!             sample energies of compton elec. and gamma
       call epcompea(cTrack%p%fm%p(4), eg, e1, cosg, cose)

       tmp=1.d0-cose*cose
       if(tmp .lt. 0.d0) then
          tmp=0.d0
          cose=-1.d0
       endif
       sine = sqrt(tmp)
       call kcossn(cs,sn)

       electron = cTrack        ! copy everything from cTrack 

       electron%w%x = cs*sine
       electron%w%y = sn*sine
       electron%w%z = cose
!            w  get new direc-cos
       call eptransVect(cTrack%w,  electron%w,  electron%w)

       call cmkptc(kelec, regptcl, -1, electron%p)
       electron%p%fm%p(4) = e1
       call epe2p(electron)
!
!                treat gamma as counterpart of electron (negative d.c)
!
       tmp=1.d0-cosg*cosg
       if(tmp .lt. 0.) then
          cosg=-1.
          tmp=0.
       endif
       sing=sqrt(tmp)
       w%x = -cs*sing
       w%y = -sn*sing
       w%z = cosg
       call eptransVect(cTrack%w,  w, cTrack%w)
       cTrack%p%fm%p(4) = eg
       call epe2p(cTrack)
!            since gamma is likely to have large energy, save first

       call eppush(cTrack)
       call eppush(electron)

       end
!
!     ************
      subroutine eppair
!     ************
       implicit none
#include  "ZepTrackv.h"
#include  "Zcode.h"
#include  "Zmass.h"

!
       type(epTrack)::  elec1, elec2


       real*8  e1,  e2,  cos1, cos2
!       real*8 cs, sn,  u,  teta1, teta2
       real*8 cs, sn,  u,  teta2
       real*8 sin1, sin2
       
       integer ic
       
       real*8 Eg
!     
       Eg = cTrack%p%fm%p(4)
!           sample higher energy of pair
       call epPrSampE(Media(MediaNo),  Eg, e1)
!            assign charge
       call rndc(u)
       if(u .lt. .5) then
          ic=-1
       else
          ic=1
       endif
!            the other electron energy
       e2 = Eg - e1
!         sample angle; smaller enery electron must be put
!         last
       call epPairAng(e2, masele, teta2) ! teta2 < pi/2

       if(teta2 .lt. 0.03d0) then
          cos2 = 1. - teta2**2/2
          sin2 = teta2
       else
          cos2 = cos(teta2)
          sin2 = sin(teta2)
       endif
!
       sin1 = sin2 * sqrt(  (e2**2-masele**2)/(e1**2-masele**2) )
       if(sin1 .lt. 0.03d0) then
          cos1 = 1.- sin1**2/2
       else
          cos1 = sqrt(1.d0 - sin1**2)
       endif

!
!          the next simplified treatment is also no problem.
!       teta1=teta2 * e2/e1
!       if(teta1 .lt. 0.03d0) then
!          cos1 = 1. - teta1**2/2
!          sin1 = teta1
!       else
!          cos1 = cos(teta1)
!          sin1 = sin(teta1)
!       endif

       elec1 = cTrack           ! copy everything first
!               sample direction cos. of 1st

       call kcossn(cs,sn)
       elec1%w%x = cs*sin1
       elec1%w%y = sn*sin1
       elec1%w%z = cos1
!           
       call eptransVect(cTrack%w, elec1%w, elec1%w)

       call cmkptc(kelec, regptcl, ic, elec1%p)
       elec1%p%fm%p(4) = e1
       call epe2p(elec1)
!              push higher energy none
       call eppush(elec1)
!               lower energy electron

       elec2 = elec1
       elec2%p%fm%p(4) = e2
       call cmkptc(kelec, antip, -ic, elec2%p)

!           treat the other one as counter part    (negative d.c)
       elec2%w%x = -cs*sin2
       elec2%w%y = -sn*sin2
       elec2%w%z = cos2
       call eptransVect(cTrack%w, elec2%w, elec2%w)
       call epe2p(elec2)
       call eppush(elec2)
       end
!      ************
       subroutine epphot
!      ************
       implicit none
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"
#include  "Zcode.h"
#include  "Zcnfig.h"
#include  "Zmass.h"           

       real*8 eout, cost, cs, sn, sint, Exray
       logical kbtest
       type(epTrack)::  elec1, xray
!
!           get Photo-electron energy
!       call epphotoEe(Media(MediaNo).pe,   < v8.0
       call epphotoEe(Media(MediaNo),
     *      cTrack%p%fm%p(4), eout, cost)

       if(kbtest(Eabsorb, BitPhotoElec)) then
!            energy absorbed by atom is Eabs = Eshell= Eg-(Ee-Me)
          Move%dE = cTrack%p%fm%p(4) - (eout - masele)
          Move%dEeff= Move%dE
          Move%dEioni = Move%dE
          SumDe = SumDe + Move%dE
!                  regard it as deposited in the media
!          if(Det.cmp(Cn).CountDE .ge. 1) then >>>>>>>>>>>>>>light
             call epLightPreUserde(1, cTrack)
             if( Move%Abort /= 0 ) then
                if( Move%Abort /=3 ) then
                   call epempty ! empty the stack
                   call epSkipUpdateNo
                else
                   Move%Abort=0
                endif
                  ! no flag is needed. since called from epint
                return
             endif
!          endif                               <<<<<<<<<<<<<
          Exray = 0.
       else
!             bit 1 is not on; characteristic x-ray emmission;
!             This was neglected in v8.71 or earlier.
!             we  assume 
!                1)  p.e effect takes place for the largest possible
!                    shell energy (Say, if Eg> K-shell energy, L-shell
!                    p.e effect is neglected and all p.e effect  is assumed
!                    to take for K-shell. 
!                2)  For such p.e effect, vacancy of electron level is 
!                    filled by X-ray emission ; No Auger electron emmission
!                    is considered. 1)+2) are good approximation.
         Exray = max( cTrack%p%fm%p(4) - (eout - masele), 0.d0)
       endif
!         emitted electron
       elec1 = cTrack
       call cmkptc(kelec, regptcl, -1, elec1%p)
       call kcossn(cs,sn)
       sint = sqrt(1.d0-cost**2)
       elec1%w%x = cs*sint
       elec1%w%y = sn*sint
       elec1%w%z = cost
       elec1%p%fm%p(4) = eout 
       call eptransVect(cTrack%w, elec1%w, elec1%w)
       call epe2p(elec1)
       call eppush(elec1)
!         emitted xray; assume isotropic
       if( Exray  .gt.  0.) then
          xray = cTrack
          call cmkptc(kphoton, 0, 0, xray%p)
          call episoAngle( xray%w )

!          call rndc(cost)
!          cost = 2.0*cost-1.0
!          call kcossn(cs,sn)
!          sint = sqrt(1.-cost**2)
!          xray.w.x = cs*sint
!          xray.w.y = sn*sint
!          xray.w.z = cost
          xray%p%fm%p(4) = Exray
          call eptransVect(cTrack%w, xray%w, xray%w)
          call epe2p(xray)
          call eppush(xray)
       endif
       end
!      ************
       subroutine epbrem
!      ************
       implicit none
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"
#include  "Zcode.h"
#include  "Zmass.h"

!
       type(epDirec)::  w
       real*8 e1, eg, theta, cs, sn, cost, sint

       e1 = cTrack%p%fm%p(4)
!             sample brems gamma energy
       call epBrSampE(Media(MediaNo), e1, eg)

!             electron energy
       cTrack%p%fm%p(4) =  e1 - eg
       call epe2p(cTrack)
!          save electron. can assume electron dose not change angle
       call eppush(cTrack)
!         see if brem g angle needed
       if(AngleB) then
!          brems g angle relative to parent electron.
          call epBremAng(e1, masele, eg, Media(MediaNo)%Zeff, theta)

          if(theta .lt. 0.03d0) then
             sint = theta
             cost = 1.- theta**2 / 2
          else
             sint = sin(theta)
             cost = cos(theta)
          endif
          call kcossn(cs,sn)
          w%x = cs*sint
          w%y = sn*sint
          w%z = cost
          call eptransVect(cTrack%w,  w, cTrack%w)
       endif
       cTrack%p%fm%p(4) = eg
       call cmkptc(kphoton, 0, 0, cTrack%p)
       call epe2p(cTrack)       
       call eppush(cTrack)

       end

!     ************
      subroutine epanih
!     ************
      implicit none
#include  "ZepTrackv.h"

#include  "Zcode.h"
#include  "Zmass.h"

!
       type(epTrack)::   gamma
       type(epDirec)::   w
      real*8 Ee, eg1, eg2, cos1, cosr, tmp, sine, sinr
      real*8 cs, sn 
!        cpy parent info.
      gamma = cTrack
      Ee = cTrack%p%fm%p(4)
      call epanihiea(Ee, eg1, eg2, cos1, cosr)

      tmp=1.d0-cos1*cos1
      if(tmp .lt. 0.d0) then
         tmp=0.d0
         cos1=-1.d0
      endif
      sine=sqrt(tmp)
      call kcossn(cs,sn)
      w%x = cs*sine
      w%y = sn*sine
      w%z = cos1
      call eptransVect(cTrack%w, w, w)
!        save hi gamma
      call cmkptc(kphoton, 0, 0, gamma%p)
      gamma%p%fm%p(4) = eg1
#ifdef SUBSTREC
      gamma%w = w
#else
      call epsubvec( w, gamma%w)
#endif
      call epe2p(gamma)
      call eppush(gamma)
!       low gamma
      tmp=1.d0-cosr*cosr
      if(tmp .lt. 0.d0) then
         tmp=0.d0
         cosr=-1.d0
      endif
      sinr=sqrt(tmp)
      gamma%w%x = -cs*sinr
      gamma%w%y = -sn*sinr
      gamma%w%z = cosr
      call eptransVect(cTrack%w, gamma%w, gamma%w)
      gamma%p%fm%p(4) = eg2
      call epe2p(gamma)
      call eppush(gamma)
      end

!     ************
      subroutine  epknoc
      implicit none
#include  "ZepTrackp.h"
#include  "ZepTrackv.h"
#include  "Zcode.h"
#include  "Zmass.h"
           
!
      integer ic
       type(epDirec)::  w
       type(epTrack)::  survival
      real*8 Ee, e1, er, cos1, cosr, sine, cs, sn, sinr, tmp
      character*80 msg

      ic = cTrack%p%charge
      Ee = cTrack%p%fm%p(4)

      if(ic .eq. -1) then
         call epmollerea(Ee, RecoilKEmin, e1, er, cos1, cosr)
!         call epmollerea(Ee,  e1, er, cos1, cosr)  ! old
      elseif(ic .eq. 1) then
         call epbhabhae(Ee, RecoilKEmin, e1, er, cos1, cosr)
!         call epbhabhae(Ee, e1, er, cos1, cosr)  ! old
      else
         write(msg,*) ' charge =',ic,' for knocon'
         call cerrorMsg(msg, 0)
      endif

      tmp=1.d0-cos1*cos1
      if(tmp .lt. 0.d0) then
         tmp=0.d0
         cos1=1.d0
      endif
      sine=sqrt(tmp)
      call kcossn(cs,sn)
      w%x = cs*sine
      w%y = sn*sine
      w%z = cos1
      call eptransVect(cTrack%w, w, w)
      survival = cTrack
      survival%p%fm%p(4) = e1
#ifdef SUBSTREC
      survival%w = w
#else
      call epsubvec(w, survival%w)
#endif
      call cmkptc(kelec, -ic, ic, survival%p)
      call epe2p(survival)
      call eppush(survival)
!                knock on electron
      tmp=1.d0-cosr*cosr
      if(tmp .lt. 0.d0) then
         tmp=0.d0
         cosr=1.d0
      endif
      sinr = sqrt(tmp)
      survival%w%x = -cs*sinr
      survival%w%y = -sn*sinr
      survival%w%z = cosr
      survival%p%fm%p(4) = er
      call cmkptc(kelec, regptcl, -1, survival%p)
      call eptransVect(cTrack%w, survival%w, survival%w)
      call epe2p(survival)
      call eppush(survival)
      end
!     ********************************
      subroutine epsync
      use modEMcontrol
      implicit none
#include  "ZepTrackv.h"
#include  "Zcode.h"

       real*8 e1, eg

       e1 = cTrack%p%fm%p(4)
!             sample sync photon  energy  ; Upsilon is from modEMcontrol
       call epsynce(e1, Upsilon, eg)
!             electron energy
       cTrack%p%fm%p(4) =  e1 - eg
       call cadjm(cTrack%p, cTrack%p)  ! adjust momentum due to energy change
!          save electron. can assume electron dose not change angle
       call eppush(cTrack)
       cTrack%p%fm%p(4) = eg  ! no direction change
       call cmkptc(kphoton, 0, 0, cTrack%p)
       call epe2p(cTrack)       
       call eppush(cTrack)
       end
!     ********************************
      subroutine epmpair
!     magneic pair production
      use modEMcontrol
      implicit none
#include  "ZepTrackv.h"
#include  "Zcode.h"

       real*8 e1, eg, chg, u

       eg = cTrack%p%fm%p(4)
!             sample pair electron of higher energy; Xai is from modEMcontrol
       call epmpaire(eg, Xai, e1)
!            higher energy electron
       cTrack%p%fm%p(4) =  e1
!          save  higher energy electron.
!          can assume electron dose not change angle
       call rndc(u)
       if(u .lt. 0.5) then
          chg = -1
       else
          chg = 1
       endif
       call cmkptc(kelec, -chg, chg,p cTrack%p)
       call cadjm(cTrack%p, cTrack%p)
       call eppush(cTrack)
       cTrack%p%fm%p(4) = eg - e1  
       cTrack%p%charge = chg
       cTrack%p%subcode = -chg
       call cadjm(cTrack%p, cTrack%p) 
       call epe2p(cTrack)       
       call eppush(cTrack)
       end
!            to be removed and use same naem one in cinteNuc.f
!     ************
      subroutine epNEPknoc
!     ************
      implicit none
#include  "ZepTrackv.h"
#include  "Zcode.h"
!     
       type(epDirec)::  w 
       type(epTrack)::  aTrack
      real*8  e1,  er, cos1, cosr, tmp
      real*8 cs, sn, sinr


      call epKnockea(cTrack%p, e1, er, cos1, cosr)
!
!     We can neglect angle of survival particle completely
!      tmp=1.d0-cos1*cos1
!      if(tmp .lt. 0.d0) then
!         tmp=0.d0
!         cos1=1.d0
!      endif
!      sine=sqrt(tmp)

       call kcossn(cs,sn)

!      w.x = cs*sine
!      w.y = sn*sine
!      w.z = cos1
!           
      aTrack = cTrack
      aTrack%p%fm%p(4) = e1
      call epe2p(aTrack)
      call eppush(aTrack)
!                knock on electron
      tmp=1.d0-cosr*cosr
      if(tmp .lt. 0.d0) then
         tmp=0.d0
         cosr=1.d0
      endif
      sinr = sqrt(tmp)
      w%x = -cs*sinr
      w%y = -sn*sinr
      w%z = cosr
!           
      call eptransVect(cTrack%w, w, cTrack%w)
      call cmkptc(kelec, regptcl, -1, cTrack%p)
      cTrack%p%fm%p(4) = er
      call epe2p(cTrack)
      call eppush(cTrack)
      end
!    *************************************
      subroutine epmagDefR(aTrack, mag, r)
      implicit none
!       get magnetic deflecton radius.  This is
!       approximate one.

#include  "ZepTrack.h"

       type(epTrack):: aTrack  ! input. charged particle
       type(epPos)::  mag  !   innput. magnetic field vector in 
                          !           the local coordinate
                          ! field strength is in T.
      real*8  r   ! output. Radius of magnetic defletion.  cm 
                  !         rough value.

      real*8 maxb, temp

      maxb = max (abs(mag%x), abs(mag%y), abs(mag%z))
      if(maxb .ne. 0) then
         temp = aTrack%p%fm%p(4)**2-aTrack%p%mass**2
         if(temp .le. 0.) then
            r = 1.d-4
         else
            r = 333.d0* sqrt(temp)/maxb/
     *       abs(aTrack%p%charge)
            r= max(r, 1.d-4)
         endif
      else
         r = 1.d10
      endif
      end
