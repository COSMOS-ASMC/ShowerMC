!         These are obso. See epdoEMInte LibLoft/EM/..
!         Pwork assumes that higher energy particles are normally
!         in the last part, and stacked later from the last to
!         save the stack area.  Due to this fact, putting ptlcs
!         in Pwork must be reversed (realized from v6.10)
!        *****************
         subroutine cbrems(pj)
!        *****************
!     use modBremPairAng  for epPBS, not needed.
         use modEMcontrol
         implicit none
#include  "ZmediaLoft.h"
#include  "Zpwork.h"
#include  "Zcode.h"
! #include  "Ztrackp.h"
!#include  "Ztrack.h"
!#include  "Ztrackv.h"
         type(ptcl),intent(in):: pj

         real*8 beg, den, cvh2den
         real*8 theta, sint, cost, sn, cs
         type(coord)::dc1
         type(coord)::w

         type(ptcl):: gamma, elec
         
!         type(track)::aTrack

!!!!           get brem gamma energy
!!!!         if(LpmEffect .and. 
!!!!     *        TrackBefMove%p%fm%p(4) .gt. LpmBremEmin) then ! use BefMove%p
!!!          if(LpmEffect ) then  !  --> LPMeffect
!!!!            den = cthick2den(TrackBefMove.pos.depth)
!!!             den = cvh2den(MovedTrack%pos%height) ! this should be better
!!!                        !  kg/m3
!!!             rrho = den*1.0d-3/Media(MediaNo)%rho ! rho  in g/cm3
!!!          else
!!!             rrho = Media(MediaNo)%rhoc          ! (not used actually)
!!!          endif
!!!!  now we use Moved.. rather than TrackBefMove
         !  if  LPMworks=T, rrho  has been set 
!     call epBrSampE(Media(MediaNo), MovedTrack%fm%p(4),  beg)
         call epBrSampE(Media(MediaNo), pj%fm%p(4),  beg)
!              next one is called inside above sub. at very High E.
!         call cbremErgLPM(TrackBefMove%p%fm%p(4), den, beg)
!         else
!            call cbremsEnergy(MovedTrack%p%fm%p(4), beg)
!         endif

!     electron
         elec = pj
         elec%fm%p(4) =  pj%fm%p(4)  - beg
         call ce2pp( elec )   ! adjuext p(1:3) ; not cw2p.
!
!         save high energy one last
         Pwork(Nproduced+2) = pj

         gamma = pj
         call cmkptc(kphoton,  kcasg, 0,  gamma)
         gamma%fm%p(4) = beg
         
!           brems angle at low energy
         if(AngleB .and.
     *      pj%fm%p(4)  < EemaxBremAngle ) then
            call epBremAng(pj%fm%p(4),
     *      pj%mass,  beg,  Media(MediaNo)%Zeff, theta)
!             call cBremAng( TrackBefMove.p.fm.p(4), TrackBefMove.p.mass, 
!     *       beg, TargetAtomicN, theta)
             ! v7.646
!!!??            call cBremAng(1, TrackBefMove%p%fm%p(4), beg, theta)
            if(theta .lt. 0.03d0) then
               sint = theta
               cost = 1.- theta**2 / 2
            else
                sint = sin(theta)
                cost = cos(theta)
             endif
             call kcossn(cs,sn)
             w(1) = cs*sint
             w(2) = sn*sint
             w(3) = cost
             call ctransVectZ(MovedTrack%vec%w, w, dc1)
             call csetDirCos(dc1, aTrack)
          endif
          call ce2p(aTrack)
          Pwork(Nproduced+1) = aTrack%p
          Nproduced = Nproduced + 2
         end
!      *****************
      subroutine cknockon(cg)
      use modEMcontrol
!     *****************  moller and bhabha scattering
       implicit none

#include  "Zcode.h"
!#include  "Ztrack.h"
!     #include  "Ztrackv.h"
#include "Zpwork.h"
#include  "Zelemagp.h"
!     
       integer cg
       real*8 e1, er, tmp, cos1, cosr, sine, cs, sn, sinr
       type(coord)::dc
       type(coord)::dc1
       type(coord)::dcr
!
       type(track)::aTrack
!
       if(cg .eq. -1) then
!          call cmollerea(MovedTrack%p%fm%p(4), RecoilKineMinE,
!     *                   e1, er, cos1, cosr)
          call epmollerea(MovedTrack%p%fm%p(4), RecoilKineMinE,
     *                   e1, er, cos1, cosr)
       else
!          call cbhabhaea(MovedTrack%p%fm%p(4), RecoilKineMinE,
!     *    e1, er, cos1, cosr)
          call epbhabhaea(MovedTrack%p%fm%p(4), RecoilKineMinE,
     *    e1, er, cos1, cosr)
          
       endif

       tmp=1.d0-cos1*cos1
       if(tmp .lt. 0.d0) then
          tmp=0.d0
          cos1=1.d0
       endif
       sine=sqrt(tmp)
       call kcossn(cs,sn)
       dc%r(1) = cs*sine
       dc%r(2) = sn*sine
       dc%r(3) = cos1
!           convert angle
       call ctransVectZ(MovedTrack%vec%w, dc, dc1)
       aTrack = MovedTrack
       aTrack%p%fm%p(4) = e1
       call csetDirCos(dc1, aTrack)
       call ce2p(aTrack)

       Pwork(Nproduced + 2) = aTrack%p
!            knock on electron
       tmp=1.d0-cosr*cosr
       if(tmp .lt. 0.d0) then
          tmp=0.d0
          cosr=1.d0
       endif
       sinr=sqrt(tmp)
       dc%r(1) = -cs*sinr
       dc%r(2) = -sn*sinr
       dc%r(3) = cosr
       call ctransVectZ(MovedTrack%vec%w, dc, dcr)
       aTrack%p%fm%p(4) = er
       call csetDirCos(dcr, aTrack)
!
       if(cg .eq. 1) then
          call cmkptc(kelec,  0,  -1,  aTrack%p)
       endif

       call ce2p(aTrack)

       Pwork(Nproduced + 1) = aTrack%p

       Nproduced = Nproduced + 2
      end
!     *****************
      subroutine canihi
      use modEMcontrol
!     *****************
      implicit none

#include  "Zcode.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"

       real*8 e1, er, cos1, cosr, tmp, sine
       real*8 cs, sn, sinr
       type(track)::aTrack
       type(coord)::dc
       type(coord)::dc1
       type(coord)::dcr
!
!     call canihiea(MovedTrack%p%fm%p(4), e1, er, cos1, cosr)
       call epanihiea(MovedTrack%p%fm%p(4), e1, er, cos1, cosr)       
!               
       tmp=1.d0-cos1*cos1
       if(tmp .lt. 0.d0) then
          tmp=0.d0
          cos1=-1.d0
       endif
       sine=sqrt(tmp)
       call kcossn(cs,sn)
       dc%r(1) = cs*sine
       dc%r(2) = sn*sine
       dc%r(3) = cos1
!             
       call ctransVectZ(MovedTrack%vec%w, dc, dc1)
!            higher gamma
       aTrack = MovedTrack
       call csetDirCos(dc1, aTrack)
       call cmkptc(kphoton, kcasg, 0,  aTrack%p)
       aTrack%p%fm%p(4) = e1
       call ce2p(aTrack)


       Nproduced = Nproduced + 1
       Pwork(Nproduced) = aTrack%p
!       call cpush(aTrack)
!                low gamma
       tmp=1.d0-cosr*cosr
       if(tmp .lt. 0.d0) then
          tmp=0.d0
          cosr=-1.d0
       endif
       sinr=sqrt(tmp)
       dc%r(1) = -cs*sinr
       dc%r(2) = -sn*sinr
       dc%r(3) = cosr
       call ctransVectZ(MovedTrack%vec%w, dc, dcr)
       call csetDirCos(dcr, aTrack)
       aTrack%p%fm%p(4) = er
       call ce2p(aTrack)


       Nproduced = Nproduced + 1
       Pwork(Nproduced) = aTrack%p
!       call cpush(aTrack)
       end

!        *****************
         subroutine cmbrem
!          magnetic brems
!     *****************
         use modEMcontrol
         implicit none
#include  "Zcode.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"

         real*8 x, beg

         type(track)::aTrack

!           get brem gamma energy
         call cmBremE(Upsilon, x) 
         beg = x * MovedTrack%p%fm%p(4)  ! gamma energy

         aTrack = MovedTrack
         aTrack%p%fm%p(4) = MovedTrack%p%fm%p(4) -beg
         call ce2p(aTrack)
         Pwork( Nproduced + 2 ) = aTrack%p   
!           stack gamma (probably lower energy than e)
         aTrack%p%fm%p(4) = beg
         call cmkptc(kphoton,  kcasg, 0,  aTrack%p)
         call ce2p(aTrack)
         Pwork(Nproduced + 1) = aTrack%p

         Nproduced = Nproduced + 2         
         end


