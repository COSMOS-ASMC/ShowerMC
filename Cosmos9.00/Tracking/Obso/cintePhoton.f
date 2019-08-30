      subroutine cintePhoton
      use modSetIntInf
      implicit none
! #include  "Ztrack.h"
! #include  "Ztrackv.h"
      integer i
      character*70 msg

!///////////////////
!              to see possibility of photon projectile 
!              with Primary is 'gamma' ;    see last
!              part of this file.
!      call cpredpmjet( MovedTrack.p,  Pwork, Nproduced)
!//////////

      if(IntInfArray(ProcessNo)%process .eq. 'pair') then
         call cpair
      elseif(IntInfArray(ProcessNo)%process .eq. 'compt') then
         call ccompt
      elseif(IntInfArray(ProcessNo)%process .eq. 'photoe') then
         call cphotoEE  
      elseif(IntInfArray(ProcessNo)%process .eq. 'photop') then
         call cphotop( MovedTrack%p )
      elseif(IntInfArray(ProcessNo)%process .eq. 'cohs') then
         call ccohs
      elseif(IntInfArray(ProcessNo)%process .eq. 'mpair') then
         call cmpair
      else
         write(msg, *) ' process for photon; Proc#=',
     *        ProcessNo,
     *        IntInfArray(ProcessNo)%process, ' undef'
         call cerrorMsg(msg,  1)
         write(0, '("energy =",g12.4, "where=",i4, " w=",g12.3)') 
     *    MovedTrack%p%fm%p(4), MovedTrack%where, MovedTrack%wgt
         write(0, * ) " coszenith=", MovedTrack%vec%coszenith
         write(0, *) ' MoveStat=',MoveStat, 'No Of inte=',
     *       NumberOfInte
         do i = 1, NumberOfInte
            write(0,*) IntInfArray(i)%process, ' dt=',
     *         IntInfArray(i)%thickness,IntInfArray(i)%length
         enddo
         write(0,*) 
         write(0, * ) " (dep,h)B//A==", TrackBefMove%pos%depth,
     *            TrackBefMove%pos%height,
     *            MovedTrack%pos%depth,  MovedTrack%pos%height
         stop
      endif
      end
!     *******************
      subroutine cpair
!     *******************
      implicit none

#include  "ZmediaLoft"      
#include  "Zcode.h"
#include  "Zmass.h"
#include  "Ztrackp.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"

!
      real*8 e1, e2, u, eg, cs, sn
      real*8 teta, teta1, teta2, cos1, sin1, cos2, sin2
      integer ica
      real*8 den, cvh2den
      type(track)::aTrack
      type(coord)::dc
      type(coord)::dce
      real*8 temp
!     
      eg = MovedTrack%p%fm%p(4)
      call epPrSampE( Meida(MediaNo), eg, e1)
!      if(LpmEffect .and. eg .gt. LpmPairEmin) then
!         den = cthick2den(TrackBefMove.pos.depth)
!         den = cvh2den(TrackBefMove%pos%height)  ! better
!         call cpairErgLPM(eg, den, e1)
!      else
!         call cpairEnergy(eg, e1)
!      endif

      e2=eg - e1
      if( e1 .gt. e2) then
!          store higher energy ptcl later
         temp = e1
         e1= e2
         e2 = temp
      endif
!         now e1 < e2
      
!            assign charge for e1
      call rndc(u)
      if(u .gt. .5) then
         ica=1
      else
         ica=-1
      endif
!     
      aTrack = MovedTrack
!      if(eg .lt. 100.e-3) then
!          take  pair angle if < 100 MeV
      if(eg .lt. 10) then
!          take  pair angle if < 10GeV
!     call cPairAng(1, eg, e1, teta1) ! teta1 < pi/2.  this is for Tsai
!     call cPairAng(e1, masele, teta1)
         call epPairAng(e1, masele, teta1)         
          if(teta1 .lt. 0.03d0) then
             cos1 = 1. - teta1**2/2
             sin1 = teta1
          else
             cos1 = cos(teta1)
             sin1 = sin(teta1)
          endif
!
          sin2 = sin1 * sqrt(  (e1**2-masele**2)/(e2**2-masele**2) )
          if(sin2 .lt. 0.03d0) then
             cos2 = 1.- sin2**2/2
          else
             cos2 = sqrt(1.d0 - sin2**2)
          endif
          call kcossn(cs, sn)

!         teta = masele/eg
!         teta1=teta* e2/eg
!         teta2=teta* e1/eg
!         cos1=1. - teta1**2/2
!         cos2=1. - teta2**2/2
!               sample direction cos. of 1st
!         sin1=teta1
!         sin2=teta2
!
         dc%r(1) = cs * sin1
         dc%r(2) = sn * sin1
         dc%r(3)=  cos1
         call ctransVectZ(MovedTrack%vec%w, dc, dce)
         call cmkptc(kelec, 0, ica, aTrack%p)
         aTrack%p%fm%p(4) = e1
         call csetDirCos(dce, aTrack)
         call ce2p(aTrack)
         Nproduced = Nproduced + 1
         Pwork(Nproduced) = aTrack%p
!            another electron (higher E)

         dc%r(1) = -cs*sin2
         dc%r(2) = -sn*sin2
         dc%r(3) = cos2
         call ctransVectZ(MovedTrack%vec%w, dc, dce)
         aTrack%p%fm%p(4) = e2
         call cmkptc(kelec, 0, -ica,  aTrack%p)
         call csetDirCos(dce, aTrack)
         call ce2p(aTrack)
         Nproduced = Nproduced + 1
         Pwork(Nproduced) = aTrack%p
      else
!          neglect pair angle
         aTrack%p%fm%p(4) = e1
         call ce2p(aTrack)
         call cmkptc(kelec, 0, ica, aTrack%p)
         Nproduced = Nproduced + 1
         Pwork(Nproduced) = aTrack%p
!     
         aTrack%p%fm%p(4) = e2
         call ce2p(aTrack)
         call cmkptc(kelec, 0, -ica, aTrack%p)
         Nproduced = Nproduced + 1
         Pwork(Nproduced) = aTrack%p
      endif
      end
!     ***********
      subroutine ccompt
!     ***********
      implicit none

#include  "Zcode.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
!
      type(track)::aTrack
      real*8 eg, e1, cs, sn, cosg, cose
      real*8 sine, tmp, sing
      type(coord)::dc
      type(coord)::dce
      type(coord)::dcg
!  
!     call ccomptea(MovedTrack%p%fm%p(4), eg, e1, cosg, cose)
      call epcomptea(MovedTrack%p%fm%p(4), eg, e1, cosg, cose)

      call kcossn(cs,sn)
!           electron direction
      tmp=max(1.d0-cose*cose, 0.d0)
      sine=sqrt(tmp)
      dc%r(1)=cs*sine
      dc%r(2)=sn*sine
      dc%r(3)=cose
      call ctransVectZ(MovedTrack%vec%w, dc, dce)
      aTrack = MovedTrack
      call cmkptc(kelec, 0, -1, aTrack%p)
      aTrack%p%fm%p(4) = e1
      call csetDirCos(dce, aTrack)
      call ce2p(aTrack)
      Nproduced = Nproduced + 1
      Pwork(Nproduced) = aTrack%p
!            gamma dicrection
      tmp=max(1.d0-cosg*cosg, 0.d0)
      sing=sqrt(tmp)
      dc%r(1) = -cs*sing
      dc%r(2) = -sn*sing
      dc%r(3) = cosg
      call ctransVectZ(MovedTrack%vec%w, dc, dcg)
      aTrack%p%fm%p(4) = eg
      call cmkptc(kphoton, kcasg, 0, aTrack%p)
      call csetDirCos(dcg, aTrack)
      call ce2p(aTrack)
      Nproduced = Nproduced + 1
      Pwork(Nproduced) = aTrack%p


      end
!     ****************
      subroutine ccohs
!     ****************
!        coherent scattering
      implicit none
#include  "Zcode.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
 !     ******************
!        since coherent scattering is effective at
!        low energies where angular distribution can be
!        approximated by (1+cos^2) dcos, we simply use
   
      type(track)::aTrack
      type(coord)::w
      real*8 cosg, tmp, cs, sn, sing
!             sample scattering angle from (1+cos^2)dcos
      call ksampRSA(cosg)
      tmp=1.d0-cosg*cosg
      sing = sqrt(tmp)
      call kcossn(cs,sn)
      
      w%r(1) = cs*sing
      w%r(2) = sn*sing
      w%r(2) = cosg
      call ctransVectZ(MovedTrack%vec%w,  w, w)
      aTrack = MovedTrack
!        energy unchaged;
      call csetDirCos(w, aTrack)
      call ce2p(aTrack)
      Nproduced = Nproduced + 1
      Pwork(Nproduced) = aTrack%p
      end
!     ****************
      subroutine cphotoEE
!     ****************
!        photo electric effect
      implicit none
#include  "Zcode.h"
#include  "Zmass.h"
#include  "ZmediaLoft.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"

!     ******************
      type(track)::aTrack
      type(coord)::w
      real*8  cs, sn, sing

      real*8 eout, eg, rEg, rEe, cost, a, tmp
!          essentially no shell down to 1keV.
      eg = MovedTrack%p%fm%p(4)

      call   epphotoEe(Media(MediaNo)  eg, eout, cost)
      aTrack = MovedTrack

      tmp=1.d0-cost*cost
      sing = sqrt(tmp)
      call kcossn(cs,sn)
      
      w%r(1) = cs*sing
      w%r(2) = sn*sing
      w%r(2) = cost
      call ctransVectZ(MovedTrack%vec%w,  w, w)
      aTrack%p%fm%p(4) = eout
      call cmkptc(kelec, 0, -1, aTrack%p)
      call csetDirCos(w, aTrack)
      call ce2p(aTrack)
      Nproduced = Nproduced + 1
      Pwork(Nproduced) = aTrack%p
      end
!     *******************
      subroutine cmpair
!         magnetic pair creation
!     *******************
      implicit none
#include  "Zcode.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
!
      real*8 e1, e2, u
      integer ica, nc
      type(track)::aTrack
!       
      call cmPairE(Xai, e2, nc)
!           e2 is higher energy fraction; change to real energy
!         store later in  working array, then higher one is
!        stored first in the stack to save the memory.
      e2 = MovedTrack%p%fm%p(4) * e2
      e1=MovedTrack%p%fm%p(4) - e2
!            assign charge for e1
      call rndc(u)
      if(u .gt. .5) then
         ica=1
      else
         ica=-1
      endif
!     
      aTrack = MovedTrack
      aTrack%p%fm%p(4) = e1
      call ce2p(aTrack)
      call cmkptc(kelec, 0, ica, aTrack%p)
      Nproduced = Nproduced + 1
      Pwork(Nproduced) = aTrack%p

      aTrack%p%fm%p(4) = e2
      call ce2p(aTrack)
      call cmkptc(kelec, 0, -ica, aTrack%p)
      Nproduced = Nproduced + 1
      Pwork(Nproduced) = aTrack%p
      end
      
!/////////////
      subroutine cpredpmjet(pj,  a, ntp)
!#include  "Zair.h"
!#include "Zptcl.h"
      type(ptcl)::pj
      integer ntp
      type(ptcl)::a(*)
!      TargetNucleonNo=14
!      TargetProtonNo= 7 
!      call cdpmjet( pj, TargetNucleonNo, TargetProtonNo,
!     *   a, ntp)
!     write(0,*) ' ntp =',ntp
      write(0,*) ' cpredpmjet should not be called'
      stop
      end
!/////////////
