       subroutine ccapnu(atms, z, a, np)
       implicit none

#include  "Zptcl.h"
#include  "Zmass.h"
#include  "Zcode.h"
       integer  atms, z
       integer np       ! output.
       type(ptcl):: a  ! output. 
!
!         atms: input.  atomic mass
!            z: input.   charge of the capturer
!         give neu(mu) which comes from captured
!         negative muon.
!        er:   recoil nutron average kinetic energy.
       real*8 er, xm, e, cs, sn, u, cost, x, sint

         parameter (er=15.e-3, xm=er/masmu/(masmu-er) )
!
         if(atms .eq. 1) then
             e=masmu
         else
!                sample x from x*exp(-x)dx with mean xm.
             call ksgmrm(1.00d0, xm, x)
!               energy of neutrino; peaked at around 100 MeV
             e= 1./(1./masmu + x )
         endif
!          make neutrino of muon type
         call cmkptc(kneumu, regptcl, 0, a)
         call kcossn(cs, sn)
         call rndc(u)
         cost=2*u-1.
         sint = sqrt(1.d0 - cost**2)
         a%fm%p(1) = cs*sint*e
         a%fm%p(2) = sn*sint*e
         a%fm%p(3) = cost*e
         a%fm%p(4) = e
         np=1
       end
