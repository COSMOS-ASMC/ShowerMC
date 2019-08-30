!    ******************************************************************
!    *                                                                *
!    *   cpimu:  pi-->mu decay with polarization
!    *                                                                *
!    ******************************************************************
!
!
!
      subroutine cpiMuDecay(pj,  mupol, a, np, polari)
      implicit none
!----      include '../../Zptcl.h'
#include  "Zptcl.h"
#include  "Zcode.h"
      
      type(ptcl):: pj    !  input. pion
      logical mupol    !  if T, polarization is taken into account
      type(ptcl):: a(*)   ! to store decay product
      integer np     !    no. of ptcls stored in a 
      real*8   polari
      integer subcode, charge
!
!           pi --> mu+neumu (muon should be set last)
      subcode =  -pj%charge
      call cmkptc(kneumu,  subcode, 0,  a(1))
      charge =  pj%charge
      call cmkptc(kmuon,  0,  charge,   a(2))
      if(pj%fm%p(4) <= pj%mass) then
!                stopping pion.
!          In  some case, Ppi /=0 seems to result in.
!          reset Ppi
         pj%fm%p(1:3) = 0.
         pj%fm%p(4) = pj%mass
      endif
      call c2bdcy(pj, a(1), a(2))
!               set polarization of muon
      if(mupol) then
         call cpimuPolari(pj, a(2), polari)
      else
         polari=0.
      endif
      np=2
      end
