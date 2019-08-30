!          test cpiMuDecay      
      implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
      
      type(ptcl):: pj    !  input. pion
      logical mupol    !  if T, polarization is taken into account
      type(ptcl):: a(2)   ! to store decay product
      integer np     !    no. of ptcls stored in a 
      real*8   polari
      integer subcode, charge
      
      integer i

      pj.fm.p(1) =0.
      pj.fm.p(2) =0.
      pj.fm.p(3) =0.
      call cmkptc(kpion, -1, 1,  pj)
      mupol = .true.
      do  i = 1, 1000000
         call cpiMuDecay(pj,  mupol, a, np, polari)
         write(*,*) a(2).fm.p(4)-a(2).mass, polari
      enddo
      end
