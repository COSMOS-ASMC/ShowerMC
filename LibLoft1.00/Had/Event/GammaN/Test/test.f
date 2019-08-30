!          particles are produced in cms.
       implicit none
#include  "Zptcl.h"
#include  "Zevhnv.h"
       type(ptcl):: a(10)

       integer i, ntp, icon, j
       real*8 sum

       Cmsp.mass = 1.67
       read(*, *) Cmsp.mass
       write(*, '(a,g15.5)' ) '# cmsp=', Cmsp.mass
       do i = 1, 1000
          call cg2pi(1, a, ntp, icon)
          sum = 0.
          do j= 1, 4
             sum = sum + a(j).fm.p(4)
          enddo
          write(*, *) (sngl( a(j).fm.p(4)), j=1, ntp), sngl( sum), icon
       enddo

       end

 
!       **************
        subroutine cg2pi(ic, a, ntp, icon)
!       **************
!          particles are produced in cms.
       implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zevhnv.h"

       type(ptcl):: a(*)

       integer ic, ntp

        real*8 w
        integer icon
!                   gp-->p pi+ pi- or gn --> n pi+ pi-
        character*70 msg
       
!            in cms.
       call cmkptc(knuc, regptcl, ic, a(1))
       call cmkptc(kpion, 0, -1, a(2))
       call cmkptc(kpion, 0, 0, a(3))
       call cmkptc(kpion, 0, 1, a(4))
       call cnbdcy(4, Cmsp.mass, a,  0, w, icon)
       ntp = 4
       end

