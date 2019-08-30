!            get max movable lenghth of a ptcl.
      subroutine cmaxMovLen(leng, thick)
!       leng:  real*8. output.  max movable length in m.
!      thick:  real*8. output.  thickness corresponding to leng in kg/m2.
!                         however, note;
!                             AlmostVacT, if Reverse=0 and height >AlmostVacH
!                             0. if Reverse = 1.
!                             0. if Reverse = 2 and  height > AlmostVacH
!                                  
!   A)   if a charged ptcl
!    A-1) compute (radius of gyro circle)/LamorDiv
!      where LamorDiv is 10 in default. Also compute length where
!      dB < 1%. Take minimum of both. 
!    A-2)   if not Reverse mode,  compute maximum gramage where cascade
!           scatteing remains very small; For a high energy electron,
!           density change must be kept small so get minimum of 
!           of the both gramage.  Convert the gramage into length.
!           if A-1) is shorter than A-2, take A-1 length and compute
!           corresponding gramage.  If A-2 is shorter,  leng and thick
!           are already obtained.
!    A-3)   if Reverse mode =1,  take A-1) and make thick = 0
!           if Reverse mode =other, take A-1) and compute thick 
!           corresponding to A-1)
!  elsse if
!   B)   a neutral particle, 
!     B-1) assume a  large length; rmg
!          If not Revesr mode, 
!            for neutrinos use such B-1) and corresponding thickness
!            (thickness is not used at all)
!            for photons, 
!                 if E> mag pair region, get length wherer dB<1 %
!                 and take min of this and rmg. (gramage not yet computed)
!            
!                 if E > LPM region, get gramage where LPM xsection
!                     remains const
!                 else  take gramage= X0*5
!                 compute  corresponding length
!                 if its rmg< length, use rmg and compute corresponding
!                 thickness else use already computed length and thick
!         if Reverse mode = 1, use leng=rmg and thick=0
!         if Reverse mode =other, use leng =rmg and corresponding thickness
!
!
      use modAtmosDef
      use modEfield
      use modEMcontrol
      implicit none
#include  "Zcode.h"
#include  "Ztrack.h"
! next  is incuded by the above  one ; Mag is defined in Ztrakv.h
!       as type(magfield):: Mag      
! #include  "Zmagfield.h"      
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zelemagp.h"
#include  "ZmediaLoft.h" 
      real*8 leng, thick
!
      real*8  ztrunc, rmg, rmgmax, L 
      real*8  clen2thick, erg
      integer jcut
      integer ka
      
      ka =  TrackBefMove%p%code
      erg = TrackBefMove%p%fm%p(4) - TrackBefMove%p%mass 
      if( erg .le. 0.) then
         leng = 0.
         thick = 0.
         return  !  *****
      endif
!             fix energy dependent truncation path
      if(TrackBefMove%p%charge .ne. 0) then
!            magnetic deflection
         call cmagDefR(TrackBefMove, Mag,  rmg)  ! get radius approx.
         rmg = rmg/LamorDiv     !  this is almost streight movable
!            get max length within which B is almost const (dB< 1 %)
         call clengSmallBC(TrackBefMove, rmgmax)
         
         rmg = min(rmg, rmgmax)

!                E field effect.  Length: momentum change < 1%.
         if( HowEfield > 0 ) then
            call cmaxEfEffLen(TrackBefMove, L)  ! L in m
            rmg = max(min(L, rmg), 0.1d0)
         endif
!            mul. scatt and lpm
         if(Reverse .eq. 0) then
!               scattering effect; streight and scattered line must be
!               not so much different ; path < ztrunc (kg/m2)
!             1 kg/m2 = 1000g/10^4 cm2 = 0.1g/cm2
            call cmaxCasLen(TrackBefMove, ztrunc)
            
            if(TrackBefMove%p%code .eq. kelec .and. 
     *         TrackBefMove%p%fm%p(4) .gt. LpmBremEmin
     *         .and. LpmEffect ) then
               ztrunc =
     *         max( min(TrackBefMove%pos%depth/10., ztrunc), 
     *              30.d0 )
            endif

            if(ka .eq. kmuon .or. ka .eq. kpion .or. 
     *         ka .eq. kkaon ) then
!                    energy loss < erg's 1%-->5*erg g/cm2 = 50kg /m2
               ztrunc =min(ztrunc, 50.d0*erg)
            endif

            ztrunc =   min(ztrunc, maxstep(TrackBefMove%where))
            call cthick2len(TrackBefMove,
     *           ztrunc, leng, thick, jcut)
            if(rmg .lt. leng) then
               thick = clen2thick(TrackBefMove%pos%height,
     *              TrackBefMove%vec%coszenith, rmg)

               leng = rmg
            endif
         elseif(Reverse .eq. 1) then
            leng = rmg
            thick = 0.  !  thick is not used at all
         else
            leng = rmg
            thick = clen2thick( TrackBefMove%pos%height,
     *         TrackBefMove%vec%coszenith, rmg)    !  used for +dE/dx 
         endif
      else
!               neutral
         rmg = 1.d5
         if(Reverse .eq. 0) then
            if(TrackBefMove%p%code .eq. kneumu .or.
     *          TrackBefMove%p%code .eq. kneue) then
               leng = rmg           ! means very large
               thick = AlmostVacT   ! not used
            else
               if(TrackBefMove%p%code .eq. kphoton) then
                  if(erg .gt. MagPairEmin .and. MagPair .ne. 0) then
                     call clengSmallBC(TrackBefMove, rmgmax)  
                     rmg = min(rmg, rmgmax)
                  endif
               endif
!                  next one cannot be compliled by ifort at
!                  opteron.  reason unknown so it is modifed.
!               if(TrackBefMove.p.code .eq. kphoton .and.
!     *              TrackBefMove.p.fm.p(4) .gt. LpmPairEmin .and. 
!     *              LpmEffect) then

               if(TrackBefMove%p%code .eq. kphoton .and. LpmEffect
     *          .and. TrackBefMove%p%fm%p(4) .gt. LpmPairEmin ) then
                  if(TrackBefMove%pos%height .lt. AlmostVacH) then
                     ztrunc = TrackBefMove%pos%depth/10.
                  else
                     ztrunc = AlmostVacT  
                  endif
               else
!     ztrunc = X0*5
                  ztrunc =Media(MediaNo)%X0g*50.d0  ! 5X0 kg/m2
               endif
               call cthick2len(TrackBefMove,
     *            ztrunc, leng, thick, jcut)
!                 thick may have been changed to shorter one.
               if(rmg .lt. leng) then
                  thick = clen2thick(TrackBefMove%pos%height,
     *            TrackBefMove%vec%coszenith, rmg)
                  leng = rmg
               else
!                 leng and thick are given
!                  thick = AlmostVacT   ! not used.
!                  leng = rmg           ! strange
               endif
            endif
         elseif(Reverse .eq. 1) then
            leng = rmg
            thick = 0.
         else
            leng = rmg
!c            if(TrackBefMove.pos.height .gt. AlmostVacH) then
!c               thick = 0.
!c            else
               thick = clen2thick(TrackBefMove%pos%height,
     *         TrackBefMove%vec%coszenith, rmg)
!c            endif
         endif
      endif
      end
!     **********************
      subroutine cmaxCasLen(aTrack, kgpm2)
      implicit none
!       get max. movable length for cascade so
!       that the scattering deflection can be
!        neglected
#include "Ztrack.h"

#include  "Ztrackp.h"
!c         #include  "Ztrackv.h"
#include  "Zelemagp.h"
#include  "ZmediaLoft.h"

      type(track)::aTrack ! input.
      real*8 kgpm2  ! output. length kg/m2

!
      real*8  ek, ttrunc


      ek = aTrack%p%fm%p(4) - aTrack%p%mass
      if(ek .gt. 1.d-3) then
!             1MeV -->5e-3 r.l   --> 360*5e-3 kg/m2~  1.8 kg/m2
         ttrunc=min( ek*5.0d0, 1.d0)
      else
!                1e-3 r.l  --> 365e-3 ~  0.3 kg/m2
         ttrunc = max(1.d-3, ek*2.d0)
      endif
      kgpm2= ttrunc*Media(MediaNo)%X0g*10.d0    !   kg/m2
      end
      subroutine cmaxEfEffLen(aTrack, L)
!           max movable length L in m
!          electric field exist.  length for
!          momentum loss/gain is < 1%.
!  dp/dt = ZeE =Z*eval*Ef  (GeV/c/s)
!         dp = Z* eval*Ef*dt = Z*eval*Ef*dt = Z*eval*Ef*L/c (GeV/c)
! 
!          dp/p ~  Z*eval*Ef*L/c/p 
!      L = dp/p /(Z*eval*Ef)*c*p
!             if Ef=1000 V/m, p=0.1GeV/c, dp/p~0.01, Z=1
!        L = 0.01 *3e8*0.1 /(0.3*1000) = 3e5/300 = 1000 m
!            p=1MeV==> 10 m
!      Ef=1.e5    50 MeV   
!            L  = 10 m
      use modBEfield 
      implicit none
#include "Ztrack.h"
      
      type(track)::aTrack  ! input current particle track
      real(8),intent(out):: L  ! length for dp/p < 1%

      real(8):: p, Ef
      real(8),parameter:: dpbyp = 0.01d0
      p = sqrt(dot_product(aTrack%p%fm%p(1:3),aTrack%p%fm%p(1:3)))

      call cgetEfield(aTrack)  ! ans is in Efld

      Ef = sqrt(dot_product( Efld(1:3), Efld(1:3))  )
!      dpbyp =  aTrack.p.chare*eval*Ef*L/3e8/p
      if(Ef == 0. ) then
         L = 1.e10
      else
         L =abs( dpbyp*p*3.0e8/(aTrack%p%charge*eval*Ef) )
      endif
      end
!    *************************************
      subroutine cmagDefR(aTrack, mag, r)
!       get magnetic deflecton radius.  This is
!       approximate one.
      implicit none

#include  "Ztrack.h"
! #include  "Zmagfield.h"
      
      type(track)::aTrack  ! input. charged particle
      type(magfield)::mag  ! innput. magnetic field
      real*8  r   ! output. Radius of magnetic defletion.  m

      real*8 maxb, EK

      if(aTrack%p%charge .eq. 0) then
         r = 1.e30
      else
               !  2013.Aug
         maxb = aTrack%vec%w%r(1)*mag%x +
     *        aTrack%vec%w%r(2)*mag%y +
     *        aTrack%vec%w%r(3)*mag%z
         maxb = abs(maxb)

!         maxb = max (abs(mag.x), abs(mag.y), abs(mag.z)) !  old one

         if(maxb .ne. 0) then
            EK = aTrack%p%fm%p(4)-aTrack%p%mass
            if( EK < 3*aTrack%p%mass ) then
!               r is smaller than true Lamor radius which
!               would be obtained with momentum
!               since K.E < P; so get p here
               EK= sqrt( sum(aTrack%p%fm%p(1:3)**2) ) ! momentum p
            endif
            r = 3.33d0*EK/maxb/abs(aTrack%p%charge)
!     for 1keV/c Z=1, maxb=0.3G,  r=10 cm
!           r= 3.33e-6 /0.3e-4 =  11e-2 m = 0.11 m ~ 10cm
!            r= max(r, 1.d-2)  
!     r= max(r, 0.1d0)  ! 2013.Aug.
!     if B~10^13 G = 10^9 T and p = 1keV/c
!     r= 3.33 e-6/1e9= 3.33e-15 m ~ proton radius. It should have
!     happened other effect before moving such distance;
!     so put r > 1micron m 
            r = max(r, 1.d-6)  !     2019 Jul. 
         else
            r = 1.e30
         endif
      endif
      end
!      ***********************
      subroutine clengSmallBC(aTrack, r)
!       get length where the change of magnetic
!     field can be regarded as small < 1 %
      use modAtmosDef
      use modEMcontrol
      implicit none

#include  "Ztrack.h"
! #include  "Zmagfield.h"                  
!  #include  "Ztrackp.h"
      ! MagChgDist
! #include  "Zearth.h"

      type(track)::aTrack ! input. r is obtaiend at this ptcl is
!                         located.
      real*8  r  ! output. within this length (m), geomag can be
!                      regarged as constant.

!     at the surface of Earth, it is about 20 km = MagChgDist
!     at larger radial distance, it becomes larger
!
      r =   aTrack%pos%radiallen/Eradius * MagChgDist

      end
