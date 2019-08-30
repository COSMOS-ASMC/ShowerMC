      subroutine epcoher(pj)
!        coherent scattering
!           since coherent scattering is effective at
!           low energies where angular distribution can be
!     approximated by (1+cos^2) dcos, we simply use this
       use modColInfo
       implicit none
#include "Zptcl.h"
#include "Zpwork.h"       
#include "Zcode.h"
       
       type(ptcl),intent(in):: pj

       type(ptcl):: gamma
       
       real(8):: w(3)

       real*8  eg, tmp, cosg
       real*8  cs, sn, sing
       real(8):: pabs
!             sample scattering angle from (1+cos^2)dcos
       call ksampRSA(cosg)
       call epcos2dir(cosg, cs, sn, w)
       pabs = sqrt( sum( pj%fm%p(1:3)**2 ))
       gamma = pj
       gamma%subcode = kcasg
       gamma%fm%p(1:3) = w(1:3)* pabs
       Nproduced = Nproduced + 1  ! assume clear before this call
!        energy unchaged;  
       Pwork( Nproduced ) = gamma
       call crot3mom(pj, Pwork, Nproduced)
       end
!     ******************
      subroutine epcmpt( pj )
      use modColInfo
      implicit none
#include  "Zcode.h"
#include "Zptcl.h"
#include "Zpwork.h"       

      type(ptcl):: pj
!

      type(ptcl):: electron
      type(ptcl):: gamma
      
       real*8 e1, eg, tmp, cosg, cose
       real*8 sine, cs, sn, sing, pabs
       real(8):: w(3)
       
!             sample energies of compton elec. and gamma
       call epcompea(pj%fm%p(4), eg, e1, cosg, cose)

       call epcos2dir(cose, cs, sn, w)


       call cmkptc(kelec, regptcl, -1, electron)
       electron%fm%p(4) = e1
       pabs = e1**2 - electron%mass**2
       if( pabs < 0.) then
          pabs = 0.
          e1 = electron%mass
       else
          pabs = sqrt(pabs)
       endif

       electron%fm%p(1:3) = w(1:3)* pabs
!
!     treat gamma as counterpart of electron (negative opposite
!       azimuth.       
!
       call cmkptc(kphoton, kcasg, 0, gamma)
       gamma%fm%p(4) = eg
       tmp=1.d0-cosg*cosg
       if(tmp .lt. 0.) then
          cosg=-1.
          tmp=0.
       endif
       sing=sqrt(tmp)

       gamma%fm%p(1) = -cs*sing * eg
       gamma%fm%p(2) = -sn*sing * eg
       gamma%fm%p(3) = cosg * eg
!     since gamma is likely to have larger energy, save first
!     (stacking will be made from leading particle)
       Nproduced =  Nproduced + 1
       Pwork( Nproduced ) = gamma
       Nproduced =  Nproduced + 1
       Pwork( Nproduced ) = electron
       call crot3mom(pj, Pwork, Nproduced)
       end
!     ************
      subroutine eppair(pj)
!     ************ pair creation by gamma
      implicit none
#include  "ZmediaLoft.h"      
#include  "Zptcl.h"
#include  "Zpwork.h"       
#include  "Zcode.h"
#include  "Zmass.h"

      type(ptcl),intent(in):: pj  !  gamma
!
      type(ptcl)::  elec


       real*8  e1,  e2,  cos1, cos2
!       real*8 cs, sn,  u,  teta1, teta2
       real*8 cs, sn,  u,  teta2
       real*8 sin1, sin2, pabs
       
       integer ic
       
       real*8 Eg
!     
       Eg = pj%fm%p(4)
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
       call epPairAng(e2, masele, teta2) ! teta2 < pi/2 ??

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

       elec = pj               ! copy everything first
       call cmkptc(kelec, ic, ic, elec)
!     sample direction cos. of 1st
       pabs = e1**2 - elec%mass**2
       if( pabs < 0.) then
          pabs = 0.
          e1 = elec%mass
       else
          pabs = sqrt(pabs)
       endif
       elec%fm%p(4) = e1          
       call kcossn(cs,sn)
       elec%fm%p(1)  = cs*sin1 * pabs
       elec%fm%p(2)  = sn*sin1 * pabs
       elec%fm%p(3) = cos1 * pabs
       ! save higher E first
       Nproduced = Nproduced  + 1
       Pwork( Nproduced ) = elec
!
       !               lower energy electron
       call cmkptc(kelec, -ic, -ic, elec)


!     treat the other one as counter part  (opposit azimuth)
       pabs = e2**2 - elec%mass**2
       if( pabs < 0.) then
          pabs = 0.
          e2 = elec%mass
       else
          pabs = sqrt(pabs)
       endif
       elec%fm%p(4) = e2       
       elec%fm%p(1) = -cs*sin2 * pabs
       elec%fm%p(2) = -sn*sin2 * pabs
       elec%fm%p(3) = cos2 * pabs

       Nproduced = Nproduced  + 1
       Pwork( Nproduced ) = elec

       call crot3mom(pj, Pwork, Nproduced)
       end
!      ************
       subroutine epphot(pj)
!      ************
       implicit none
#include  "Zcode.h"
#include  "Zptcl.h"
#include  "Zmass.h"       
#include  "ZmediaLoft.h"
#include  "Zpwork.h"
       
       type(ptcl),intent(in):: pj   ! gamma
       
       real*8 eout, cost, cs, sn, sint, Exray, pabs
       type(ptcl)::  elec, xray
!
!           get Photo-electron energy
!       call epphotoEe(Media(MediaNo).pe,   < v8.0
       call epphotoEe(Media(MediaNo), pj%fm%p(4), eout, cost)
       Exray = max( pj%fm%p(4) - (eout - masele), 0.d0)

!     emitted electron
       call cmkptc(kelec, regptcl, -1, elec)
       elec%fm%p(4) = eout
!     iso troic angle
       call episoAngle(elec%fm%p) ! direction cos
       pabs = eout**2 - elec%mass**2
       if(pabs < 0.) then
          pabs = 0.
       else
          pabs = sqrt(pabs)
       endif
       elec%fm%p(1:3) =  elec%fm%p(1:3) * pabs
!     emitted xray; assume isotropic
       Nproduced = Nproduced + 1
       Pwork(Nproduced)= elec
       if( Exray  .gt.  0.) then
          call cmkptc(kphoton, kChaX, 0, xray)
          xray%fm%p(4) = Exray
          call episoAngle( xray%fm%p )
          pabs = Exray
          xray%fm%p(1:3) = xray%fm%p(1:3) * pabs
          Nproduced = Nproduced + 1
          Pwork(Nproduced)= elec
       endif
       call crot3mom(pj, Pwork, Nproduced)
       end

!
************
      subroutine epbrem( pj )
      use modEMcontrol
!      ************
       implicit none
#include  "Zptcl.h"
#include  "ZmediaLoft.h"       
#include  "Zpwork.h"       
#include  "Zcode.h"
#include  "Zmass.h"       

       type(ptcl),intent(in):: pj

!
       type(ptcl):: elec, gamma
       real(8):: w(3)
       real*8 e1, eg, theta, cs, sn, cost, sint
       real(8):: pabs
       integer:: chg

       elec = pj
       call cmkptc(kphoton, 0, 0,  gamma)

       e1 = pj%fm%p(4)
!             sample brems gamma energy
       call epBrSampE(Media(MediaNo), e1, eg)
!             electron energy
       e1 =  e1 - eg
       elec%fm%p(4) = e1
       pabs = e1**2 - elec%mass**2
       if(pabs > 0.) then
          pabs = sqrt(pabs)
       else
          pabs = 0.
       endif
!        electron dose not change the direction       
       elec%fm%p(1) = 0.
       elec%fm%p(2) = 0.
       elec%fm%p(3) = pabs
!            save electron first
       Nproduced = Nproduced + 1
       Pwork( Nproduced ) = elec       
!         gamma       
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
       else
          sint = 0.
          cost = 1.d0
       endif
       
       call kcossn(cs,sn)
       w(1) = cs*sint
       w(2) = sn*sint
       w(3) = cost
       pabs = eg
       gamma%fm%p(1:3) = w(1:3)*pabs
       gamma%fm%p(4) = eg
!          save gamma 
       Nproduced = Nproduced + 1
       Pwork( Nproduced ) = gamma
       call crot3mom(pj, Pwork, Nproduced)

       end

!     ************
      subroutine epanih( pj )
!     ************
      implicit none
#include  "Zptcl.h"
#include  "Zpwork.h"
#include  "Zcode.h"
#include  "Zmass.h"

      type(ptcl),intent(in):: pj  ! e+
      
!
      type(ptcl)::   gamma
      real(8)::   w(3)

      real*8 Ee, eg1, eg2, cos1, cos2, tmp, sin1, sin2
      real*8 cs, sn 

      call cmkptc( kphoton, 0, 0, gamma)
      Ee = pj%fm%p(4)
      call epanihiea(Ee, eg1, eg2, cos1, cos2)

      call epcos2dir(cos1, cs, sn, w)
!        save hi gamma
      call cmkptc(kphoton, 0, 0, gamma)
      gamma%fm%p(4) = eg1
      gamma%fm%p(1:3)= w(1:3)* eg1
      Nproduced = Nproduced + 1
      Pwork( Nproduced ) = gamma
!       low gamma
      tmp=1.d0-cos2*cos2
      if(tmp .lt. 0.d0) then
         tmp=0.d0
         cos2=-1.d0
      endif
      sin2=sqrt(tmp)
      gamma%fm%p(1) = -cs*sin2 * eg2
      gamma%fm%p(2) = -sn*sin2 * eg2
      gamma%fm%p(3) = cos2 * eg2
      gamma%fm%p(4) = eg2
      Nproduced = Nproduced + 1
      Pwork( Nproduced ) = gamma

      call crot3mom(pj, Pwork, Nproduced)

      end

!     ************
      subroutine  epknoc( pj )
!     bhabha or moller scattering
      use modEMcontrol
      implicit none
#include  "Zptcl.h"
#include  "Zpwork.h"
#include  "Zcode.h"
#include  "Zmass.h"

      type(ptcl),intent(in):: pj
!
      integer ic
      real(8)::  w(3)
      type(ptcl)::  survival, elecr
      real*8 Ee, e1, er, cos1, cosr, sine, cs, sn, sinr, tmp
      character*80 msg
      real(8):: pabs
      
      ic = pj%charge
      Ee = pj%fm%p(4)

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

      call epcos2dir(cos1, cs, sn, w) 

      survival = pj

      pabs = e1**2 - masele**2
      if( pabs < 0.) then
         pabs = 0.
         e1=masele
      else
         pabs = sqrt(pabs)
      endif
      survival%fm%p(4) = e1   
      survival%fm%p(1:3) = w(1:3)* pabs

!     knock on electron
      pabs = er**2 - masele**2
      if( pabs < 0.) then
         pabs = 0.
         er=masele
      else
         pabs = sqrt(pabs)
      endif

      tmp=1.d0-cosr*cosr
      if(tmp .lt. 0.d0) then
         tmp=0.d0
         cosr=1.d0
      endif
      sinr = sqrt(tmp)
      call cmkptc(kelec, regptcl, -1, elecr)
      elecr%fm%p(1) = -cs*sinr * pabs
      elecr%fm%p(2) = -sn*sinr * pabs
      elecr%fm%p(3) = cosr * pabs
      elecr%fm%p(4) = er

      Nproduced = Nproduced  + 1
      PWork( Nproduced ) = survival
      Nproduced = Nproduced  + 1
      PWork( Nproduced ) = elecr

      call crot3mom(pj, Pwork, Nproduced)
      end
!     ********************************
      subroutine epsync(pj)
!          by electron
      use modEMcontrol
      implicit none
#include  "Zptcl.h"
#include  "Zpwork.h"      
#include  "Zcode.h"

      type(ptcl),intent(in):: pj
      
      
       real*8 e1, eg
       
       type(ptcl):: elec
       type(ptcl):: xray
       
       e1 = pj%fm%p(4)
!             sample sync photon  energy  ; Upsilon is from modEMcontrol
       call epsynce(e1, Upsilon, eg)
!       electron energy
       elec = pj
       elec%fm%p(4)=  e1 - eg
       call cadjm(elec, elec)  ! adjust momentum due to energy change
!          save electron. can assume electron dose not change angle

       call cmkptc(kphoton, ksync, 0, xray)
       xray%fm%p(1:3) = elec%fm%p(1:3)
       xray%fm%p(4) = eg        ! no direction change
!     at near erth  Eg (max) ~ 100(Ee/1e6)**2 (GeV)
!     i.e, Ee= 10^12 eV  Egmax~ 100keV
!           10^13 eV, Egmax ~0.01 GeV = 10 MeV.
!           10^15 eV, Egmax ~100 GeV.
!           10^19 eV,       ~10^19eV  over this Eg ~Ee (OK??)
       
       call cadjm(xray, xray) ! adjust momentum due to energy change

       Nproduced = Nproduced + 1
       Pwork( Nproduced ) = elec
       Nproduced = Nproduced + 1
       Pwork( Nproduced ) = xray
!!       call crot3mom(pj, Pwork, Nproduced)   !! not needed already in Exyz

       end
!     ********************************
      subroutine epmpair( pj )
!     magneic pair production
      use modEMcontrol
      implicit none
#include  "Zptcl.h"
#include  "Zpwork.h"
#include  "Zcode.h"

      type(ptcl),intent(in)::  pj    ! gamma


      type(ptcl):: elec1, elec2
      
      real*8 e1, eg,  u
      integer:: chg

      eg = pj%fm%p(4)
!             sample pair electron of higher energy; Xai is from modEMcontrol
      call epmpaire(eg, Xai, e1)
      call rndc(u)
      if(u .lt. 0.5) then
         chg = -1
      else
         chg = 1
      endif
      call cmkptc(kelec, chg, chg, elec1)
!     higher energy electron
      elec1%fm%p(1:3) = pj%fm%p(1:3)
!        save  higher energy electron.      
      elec1%fm%p(4) =  e1
!          can assume electron dose not change angle
      call cadjm(elec1, elec1)
      elec2 = elec1
      elec2%fm%p(4) = eg - e1
      elec2%charge = -chg
      elec2%subcode = -chg
      call cadjm(elec2, elec2)

      Nproduced = Nproduced + 1
      Pwork( Nproduced ) = elec1
      Nproduced = Nproduced + 1
      Pwork( Nproduced ) = elec2

!!      call crot3mom(pj,  Pwork, Nproduced)   !! no need.  already in Exyz
      end
