!     *****************************************************************
!     *                                                               *
!     * csptxy:  set ptx, pty 
!     *                                                               *
!     *****************************************************************
!
!
      subroutine csptxy(a,  nt)
!        a(nt):  type ptcl. Input.  At this moment, pt is assumed
!                to be in a(i).fm.p(3) (=ptz) position. 
!
       implicit none

#include  "Zptcl.h"
!
       integer nt
       type(ptcl):: a(nt)
!
       integer i
       real*8 cs, sn, pt
!
       do  i=1, nt
!          sample cos and sin of random azimuthal angle
         call kcossn(cs,sn)
         pt = a(i)%fm%p(3)
!          store ptx,pty
         a(i)%fm%p(1) = pt * cs
         a(i)%fm%p(2) = pt * sn
       enddo
      end
