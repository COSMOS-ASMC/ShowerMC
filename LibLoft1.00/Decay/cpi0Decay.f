!    ******************************************************************
!    *                                                                *
!    *   cpi0Decay:  one pi0 is made to decay into two gammas         *
!    *                                         or g e e~              *
!    ******************************************************************
!
!
      subroutine cpi0Decay(pj, a, np)
      implicit none
!----      include '../../Zptcl.h'
#include  "Zptcl.h"
!----      include '../../Zcode.h'
#include  "Zcode.h"

      integer np         ! output. no. of ptcles produced
      type(ptcl):: pj   ! input. pi 0 
      type(ptcl):: a(*) ! output. to contain produced ptcls


      real*8 u,  w
      integer icon, i

      call rndc(u)
      if(u .lt. .98798) then
         call cmkptc(kphoton, kdirectg, 0, a(1))
         a(2) = a(1)
!           pi0--> 2 gamma (special routine for massless case)
         call c2bdc0(pj, a(1), a(2))
         np=2
      else
!            pi0-->g+ e + e~
!        because the 3 body decay prob. is small, rapid k3bdcy needs
!        not be used.
!           3  body pure phase space
         call cmkptc(kphoton, kdirectg, 0, a(1))
         call cmkptc(kelec, regptcl, -1, a(2))
         call cmkptc(kelec, antip,  1, a(3))
         call cnbdcy(3, pj%mass, a, 0, w, icon)
         np=3
!             boost
         do i = 1, np
            call cibst1(i, pj, a(i), a(i))
         enddo
      endif
      end

