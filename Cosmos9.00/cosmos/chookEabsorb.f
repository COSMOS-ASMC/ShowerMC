      subroutine chookEabsorbi(info)
!           init for each event
      implicit none
#include "Zmaxdef.h"
#include "Zcode.h"
#include "Ztrack.h"
! #include "Zmagfield.h"
#include "Ztrackv.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"
#include "Zabsorb.h"
      integer info ! not used now
      integer i
      do i = 0, NoOfSites + 1
!            i = 0 and NoOfSite+1 will not be used now. 
         dEbydEdx(i) = 0.
         dEbyDeath(i) = 0.
         dEbyDeathG(i) = 0.
         dEbyDeathE(i) = 0.
         dEbyDeathMuPiK(i) = 0.
         dEbyDeathNeu(i) = 0.
         dEbyDeathP(i) = 0.
         dEbyDeathNut(i) = 0.
         dEbyDeathO(i) = 0.
      enddo
      do i = 1, 7
!          i for g,e,mu,pi,K, n, (neu and others)
!             Energy crashing to the given layer
!             default is the NoOfSites. Can be changed by Eabosrb(2)
         Ecrash(i) = 0.
!             Energy escaping to the upper bound( default is injection
!             height + 1m)  
         Espace(i) = 0.
      enddo
!          counter used for checking energy conservation 
!          at the multiple producion.   
!           negative value is generated energy is < ininial energy
!           positive value is generated energy is > initial energy 
      do i = 1, 2
!         i=1 for max Ebreak, i=2 Relative Break value for that event
         MaxEbreak(i) = 0.
!         i=1 for max RelEgreak, i=2 E Break value for that event
         MaxRelEbreak(i) = 0.
      enddo

      SumEdiff = 0.
      SumAbsEdiff = 0.

      end

!>   This is called when Eabsorb != 0 and
!!   when a charged particle runs from a 
!!   to b and deposits energy dE (GeV) to the  Air.
!! @param a
!! @param b
!! @param dE
!! @param info
!! @todo test 2016/2/16
      subroutine chookEabsorb( a, b, dE, info )
      implicit none
#include "Zmaxdef.h"
#include "Zobs.h"
#include "Ztrack.h"
#include "Zabsorb.h"
!   This is called when Eabsorb != 0 and
!   when a charged particle runs from a 
!   to b and deposits energy dE (GeV) to the  Air.
!             
      type(track)::  a  ! input. charged particle track info. at a.
      type(track)::  b  ! input. charged particle track info. at b.
      real*8  dE  ! input.  energy deposit in GeV; no weight is applied yet
      integer info ! in/out. not used now

      if(a%where .ge. 1)  then
         dEbydEdx(a%where)= dEbydEdx(a%where) + dE*a%wgt
      endif
      end
      subroutine chookEabsorbD( a, dE, info )
      implicit none
#include "Zmaxdef.h"
#include "Zcode.h"
#include "Zobs.h"
#include "Ztrack.h"
#include "Zabsorb.h"
!   This is called when Eabsorb != 0 and
!   when a  particle energy becomes < Emin (info=0) for 
!   its traveling time becomes > limit or (info=2)
!   its angle relative to the 1ry becomes > limit.(info=4) 
!   Whether this is called or not depends on the   particle
!   and bit in Eabsorb.  dE is energy that can be regarded
!   as absorbed in the Air. (GeV) eventually.
!   bit 1 is the LSB of Eabsorb.
!     
!   bit   particle
!    1     photon: used to absorb shell energy at 
!          photoelectric effect. This bit is not used in the Air.  
!    2    photon.  
!    3    e+/e-
!    4    proton
!    5    neutron
!    6    anti-N
!    7    decaying prtcl 
!    8    others
!             
!***** Normally Eabsorb=6 (110 in bit pattarn) is enough.****
!
!
      type(track)::  a  ! input.  a particle that is < Emin 
                          ! at birth
      real*8  dE  ! input.  energy which is supposed to  be emitted by 
                  ! the dying particle
      integer info ! in/out. not used now.
      if(a%where .ge. 1) then
         if(a%p%code .eq. kneue .or. a%p%code .eq. kneumu) then
           !  neutrino
            dEbyDeathNeu(a%where) = dEbyDeathNeu(a%where) + dE*a%wgt
         elseif(a%p%code .eq. knuc .and. a%p%charge .eq. 0 .and.
     *          a%p%subcode .eq. regptcl ) then
            !  low E nutron
            dEbyDeathNut(a%where) = dEbyDeathNut(a%where) + dE*a%wgt
         else
!                  next one and above are kept same as older versions for 
!                  compativilty 
            dEbyDeath(a%where) = dEbyDeath(a%where) + dE*a%wgt

!                  we further put details 
            if(a%p%code .eq. kphoton ) then
               dEbyDeathG(a%where) = dEbyDeathG(a%where) + dE*a%wgt
            elseif(a%p%code .eq. kelec ) then
               dEbyDeathE(a%where) = dEbyDeathE(a%where) + dE*a%wgt
            elseif( a%p%code .le. kkaon )  then
               dEbyDeathMuPiK(a%where) = dEbyDeathMuPiK(a%where)
     *             + dE*a%wgt
            elseif(a%p%code .eq. knuc .and. a%p%charge .eq. 1 ) then
!                 p
               dEbyDeathP(a%where) = dEbyDeathP(a%where) + dE*a%wgt
            else
!                pbar, nbar, heavy,  others 
               dEbyDeathO(a%where) = dEbyDeathO(a%where) + dE*a%wgt
            endif
         endif
      endif
      end

      subroutine chookEabsorbB(a,info)
      use modEMcontrol
      implicit none
#include "Zmaxdef.h"
#include "Zcode.h"
#include "Ztrack.h"
! #include "Zmagfield.h"
#include "Ztrackp.h"
#include "Ztrackv.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"
#include "Zabsorb.h"
!  This is called when Eabsorb(1) !=0 and
!    a particle crosses an observation level
!  info=2:  normal observation level
!      =1:  upper boundary (whcih is equal to or higher than the injection
!           height)
!      =3:  lower boundary (which is equal to or lower than the last
!           observation depth).  
! 
      type(track)::  a  ! input.  a particle that croses the observation level
      integer info  ! input.  see above

      integer code
      integer lv
      logical count

      count =.false. 
      if(info .eq. 2 .and. a%where .eq. EabsorbL )  then
!     Eabsorb(2) is now  EabsorbL;
!          the layer number where the user want to
!           take sum of particle energy reaching there from above.
!                ptcl comes to the specified level OR
!                upper boundary.  This condition neglect
!                ptcles reaching the real lower boundary
!               (if Eabosrb(2) is not NoOfSites+1)
         if( a%vec%coszenith .gt. 0.)   then
!          if info=2,  this is to count energy  reaching
!          the level from above, we discard going up ptcls
!          (but it may come down, so some over count may happen) 
            count = .true.
            lv = 1
         endif
      elseif( info .eq. 1 ) then
         lv = 2      ! escape to space
         count = .true.
      endif
! 
!            
      if(count) then
         code = a%p%code
!              for other ptcls than g,e,mu, pi, K, N, use 7
         if(code .gt. 7) code=7
!             at the last layer we see sum of the particle energy
!             what energy we should use here is somewhat annoying 
!             point.  We use total energy here.  
         if(lv .eq. 1) then
!            Ecrash(code) = Ecrash(code) +  a.p.fm.p(4)*a.wgt
            Ecrash(code) = Ecrash(code) + (a%p%fm%p(4)-a%p%mass)*a%wgt
         else
!           Espace(code) = Espace(code) + a.p.fm.p(4)*a.wgt
            Espace(code) = Espace(code) +  (a%p%fm%p(4)-a%p%mass)*a%wgt
         endif
      endif
      end
      subroutine chookEabsorbC(a, n, p, info)
!     use modXsecMedia
      use modColInfo
      implicit none
#include "Zmaxdef.h"
#include "Zcode.h"
#include "Zmass.h"
#include "Zobs.h"
#include "Ztrack.h"
!  #include "Zmagfield.h"      
#include "Ztrackv.h"
#include "Zabsorb.h"

      type(track)::  a    ! input.  incident particle track
                           !  that made the collision
      integer n     ! input. number of procuded ptcls in the collision
      type(ptcl)::  p(n) ! input. produced ptpcls.
      integer info  ! input. not used now.
      
      real*8 Eout, diff, reldiff, Ein
      real*8 diff1, diff2, Ein1, Ein2
      integer i
!          We check conservation above  this energy (GeV).
      real*8 Ebig/5.d3/
      save
!          since target is not well known and  nucleon there
!          has Fermi momentum, we simply assume rest mass.
!          For E>5TeV, Ar mass ~40GeV/5000 < 4/500 < 1 % 
!          So we can see the conservation neglecting mass term
!            
      if( a%p%fm%p(4) .gt. Ebig ) then
         Ein1 = a%p%fm%p(4) + masn*(TargetNucleonNo-TargetProtonNo) +
     *       masp*TargetProtonNo
         Ein2 = a%p%fm%p(4) + masp
         Eout = 0.
         do i = 1, n
            Eout = Eout + p(i)%fm%p(4)
         enddo

         diff1 = Eout - Ein1
         diff2 = Eout - Ein2
         Ein = Ein1
         diff = diff1
!              take smaller diff one for Ein
         if( abs(diff2) .lt. abs(diff1)) then
            diff = diff2
            Ein = Ein2
         endif
         reldiff = (Eout/Ein -1.)

         if( abs(reldiff) .gt. 0.02 ) then
            call chookEabsorbW(a, n, p, Eout, diff)
         endif

         if( abs(diff) .gt. abs(MaxEbreak(1)) ) then
            MaxEbreak(1) = diff
            MaxEbreak(2) = reldiff
         endif

         if( abs(reldiff) .gt. abs(MaxRelEbreak(1)) ) then
            MaxRelEbreak(1) = reldiff
            MaxRelEbreak(2) = diff
         endif

         SumEdiff = SumEdiff + diff
         SumAbsEdiff = SumAbsEdiff + abs(diff)
      endif
      end
      subroutine chookEabsorbW(a, n, p, Eout, diff)
!     use modXsecMedia
      use modColInfo
      implicit none
#include "Zmaxdef.h"
#include "Zcode.h"
#include "Zmass.h"
#include "Zobs.h"
#include "Ztrack.h"
! #include "Zmagfield.h"      
#include "Ztrackv.h"
#include "Zabsorb.h"

      type(track)::  a    ! input.  incident particle track
                           !  that made the collision
      integer n     ! input. number of procuded ptcls in the collision
      type(ptcl)::  p(n) ! input. produced ptpcls.
      real*8 Eout, diff
      integer nevent, ntevent


      call cqEventNo(nevent, ntevent)

      write(0,*) "=================  Warning: Ebreak ============="
      write(0,*) " event #=", nevent, " 1ry Energy=", inci%p%fm%p(4)
      write(0,*) "incident(not 1ry) code=", a%p%code,
     *     " subcode=", a%p%subcode, " Ein~", a%p%fm%p(4)
      write(0,*) "Target(A,Z)=", TargetNucleonNo, TargetProtonNo
      write(0,*) "No. of generted particles =", n
      write(0,*) "Sum of outgoing ptcl Energy:Eout=", Eout
      write(0,*) "dE=(Eout-Ein)=", diff, " dE/Ein~",  
     *          diff/a%p%fm%p(4)
      write(0,*) "(Eout-Ein)/1ryE=",  diff/inci%p%fm%p(4)
      write(0,*) "Height=", a%pos%height," depth=",a%pos%depth/10.
      write(0,*) " where=", a%where, " weight=", a%wgt
      write(0,*) "================================================"

      end
