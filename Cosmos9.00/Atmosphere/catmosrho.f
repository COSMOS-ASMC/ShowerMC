! -----------------------------------------
!         density of air at x from starting point rs and angle coss
!        
      real*8 function catmosrho(x)
      use modAtmosDef
      implicit none
      real*8 x

!----      include 'Zearth.h'
! #include  "Zearth.h"
      real*8 cvh2den, zp
      common /ccatmosrho/ coss, rs
      real*8 coss, rs
!
      real*8 cnewh
!
!         radius at x
      
      zp = cnewh(rs, coss, x) - Eradius
      catmosrho = cvh2den(zp)
      end
