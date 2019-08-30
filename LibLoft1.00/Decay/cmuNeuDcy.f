!       mu---> e + neue + neumu; inclusive treatment 
!      some kind of Exclusive treatment
!     
      subroutine cmuNeuDcy(pj, polari, inclusive, a, np)
      implicit none

#include  "Zptcl.h"
#include  "Zmass.h"
#include  "Zcode.h"
      
      integer np     ! output.  no. of produced particles
      real*8 polari  ! input.   muon polarization
      integer inclusive !input.  0--> inclusive treatment
                        !       >0--> Exclusive treatment
                        !   Only Electron is correctly sampled
                        !   4 momenta of two neutrinos are sampled
                        !   so that the total 4 momentum be conserved.
                        !   Since the formula assumes Me=0, we 
                        !   add electron mass to the total energy
      type(ptcl)::pj  ! input. muon
      type(ptcl):: a(*)  ! output.  produced ptcls
!
      type(ptcl):: twoNeu
      real*8 po, f, ecm, cosa, lam
      integer  i, charge, subcode
      logical  try
   
      real*8 pab

      try = .true.
!           make 3 ptcls; neu_e, neu_mu, e
      subcode = -pj%charge
      call cmkptc(kneue, subcode, 0, a(1))
      subcode = pj%charge
      call cmkptc(kneumu, subcode, 0, a(2))
      charge = pj%charge
      call cmkptc(kelec, 0, charge,  a(3))
      if(pj%charge .eq. 1) then
!                for mu+, polarization should be made to inverse sign
!                because e+ goes same as mu+ polarizaion while
!                neue trino opposit.
         po=-polari
      else
         po=polari
      endif
      do while ( try )
!           we fix electron first which is detectabl
!               electron ; paralell to neu_mu
!         x^2(3-2x + (1-2x)Pcos) dxdcos
!            first integrate by cos and get energy distribution
!         x^2(3-2x) dx =2x^2(1-x)+x^2: this is the same as  neu_mu energy
!            (f=x=2E/Mmu)

         call csNeumuEMu(f)
         ecm=f*pj%mass/2
         if(ecm .le. a(3)%mass) cycle
            
!         angluar distribution is
!        ( 1 + lam* cos ) dcos type
!
         lam = (1.-2*f)/(3.-2*f) *po

         call ksampLin(lam, 1.d0, -1.d0, 1.d0, cosa)
!             set px,py,pz
         pab =sqrt( ecm**2 - a(3)%mass**2 )
         call cpCos2pxyz(cosa, pab, a(3)%fm)
         a(3)%fm%p(4) = ecm
!         -------------
         if( inclusive .eq. 0 ) then
!            inclusive treatment
!           sample energy in f=2e*/mmu: e* is at muon rest sytem
!           of neue
            call csampNeueEMu(f)
            ecm=f*pj%mass/2 
!              sample decay angle of neue at muon rest system
!            (1+Pcos)dcos
            call ksampLin(po, 1.d0, -1.d0, 1.d0, cosa)
!              set random momentum about azimuth (px,py,pz)
            call cpCos2pxyz(cosa, ecm, a(1)%fm)
            a(1)%fm%p(4) = ecm
!                 since inclusive, no conservation tried
!          neumu; sample energy in f=2e*/mmu: e* is at muon rest sytem
            call csNeumuEMu(f)
            ecm=f*pj%mass/2
!              sample decay angle of neumu at muon rest system
!             (1+lambda P cos) dcos with
!             lambda = (1-2f)/(3-2f)
            lam=(1.-2*f)/(3.-2*f)
            call ksampLin(lam*po, 1.d0, -1.d0, 1.d0, cosa)
!             set px,py,pz
            call cpCos2pxyz(cosa, ecm, a(2)%fm)
            a(2)%fm%p(4) = ecm
         else
!             exclusive; but not possible
!             we are interested in only electron energy so that
!             other two neus are sampled from remaining  4 momenta
            twoNeu%fm%p(4) = pj%mass - a(3)%fm%p(4)
            if(twoNeu%fm%p(4) .le.  0.) cycle
            twoNeu%fm%p(1) = - a(3)%fm%p(1)
            twoNeu%fm%p(2) = - a(3)%fm%p(2)
            twoNeu%fm%p(3) = - a(3)%fm%p(3)
            twoNeu%mass =
     *        sqrt( pj%mass**2 -2*pj%mass*a(3)%fm%p(4) + a(3)%mass**2)

            call c2bdcy(twoNeu, a(1), a(2))

!              another choice may be to sample all 3 by inclusive way
!           and forced conservation by next conformal transformation
!            call cnbdc2(3, pj.mass, a)
!             since we assume massless ; no further transformation
!             is needed
         endif
! **************
!         pab =sqrt(a(3).fm.p(1)**2 + a(3).fm.p(2)**2
!     *   + a(3).fm.p(3)**2 )
!         write(*,'("cm ",5G14.5)' ) 
!     *   a(1).fm.p(4)+a(2).fm.p(4)+a(3).fm.p(4),
!     *   2* a(3).fm.p(4)/pj.mass, a(3).fm.p(3)/pab,
!     *   po, lam
!     *   a(1).fm.p(1)+a(2).fm.p(1)+a(3).fm.p(1),
!     *   a(1).fm.p(2)+a(2).fm.p(2)+a(3).fm.p(2),
!     *   a(1).fm.p(3)+a(2).fm.p(3)+a(3).fm.p(3)
! *************
         np = 3
!               boost to lab.
         do i = 1, np
            call cibstPol(i, pj, a(i), a(i) )
         enddo
!  ***************
!         write(*,'("lb ",4G14.5)' ) 
!     *   a(1).fm.p(4), a(2).fm.p(4), a(3).fm.p(4), 
!     *   a(1).fm.p(4)+a(2).fm.p(4)+a(3).fm.p(4)
! ***************
!         since approximation here is to assume Me=0
!         we get Ee < Me sometimes (normally 5/10^5 or less )
!         so we reject such an event
         try =  a(3)%fm%p(4) .le. a(3)%mass
      enddo
      end
