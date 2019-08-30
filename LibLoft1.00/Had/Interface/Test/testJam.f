      implicit none
#include "Zcode.h"      
#include "Zptcl.h"      
      type(ptcl):: pj  !  input. projectile
      type(ptcl):: a(1000)
      integer ia, iz, ntp, i, now(2)
      real*8 sig, tgA, tgZ

      call cmkSeed(0, now)   ! make seed using timer and hostname            
      call rnd1r(now)

      call cmkptc(9, 56, 26, pj)
      pj.fm.p(1) =  0.
      pj.fm.p(2) =  0.
      pj.fm.p(3) =  1000.
      pj.fm.p(4) = sqrt( pj.fm.p(3)**2 + pj.mass**2)
      ia = 207
      iz =  82
      tgA = ia
      tgZ = tgA*0.4  ! if pj is not heavy, can be arbitrary.
      call cinelx(pj,  tgA, tgZ, sig)
      call cjamEvent(pj, ia, iz, sig, a, ntp)
      do i = 1, ntp
         write(*,*) i, a(i).code, a(i).charge, a(i).fm.p(4)
      enddo
      end
