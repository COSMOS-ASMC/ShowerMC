#include "ZsaveStruc.h"
!cc            test ctransVectZx,  ciTransVectZx
!       implicit none
!       include 'Zcoord.h'
!       record /coord/ wa, wb, dir1, dir2
!       real*8 norm
! 
! 
! 
!       wa.r(1) = .50
!       wa.r(2) = 0.2
!       wa.r(3) = sqrt(1. - (wa.r(1)**2+wa.r(2)**2))
!
!       wb.r(1) = wa.r(2)
!       wb.r(2) = -wa.r(1)
!       wb.r(3) = 0.
!       norm = sqrt( wb.r(1)**2 + wb.r(2)**2 + wb.r(3)**2)
!       wb.r(1) = wb.r(1)/norm
!       wb.r(2) = wb.r(2)/norm
!       wb.r(3) = wb.r(3)/norm
!
!
!       dir1.r(1) = sqrt(2.)/3.
!       dir1.r(2) = dir1.r(1)*1.1
!       dir1.r(3) = sqrt(1. - dir1.r(1)**2 - dir1.r(2)**2)
!       write(*,*) dir1.r(1), dir1.r(2), dir1.r(3)
!       call ctransVectZx(1, wa, wb, dir1, dir2)
!       write(*,*) dir2.x, dir2.y, dir2.z
!       call ciTransVectZx(1, wa, wb, dir2, dir2)
!       write(*,*) dir2.x, dir2.y, dir2.z
!       call ctransVectZx(2, wa, wb, dir1, dir1)
!       write(*,*) dir1.r(1), dir1.r(2), dir1.r(3)
!       end
!-------------------------------------------------------------
       subroutine ctransVectZx(init, zax, xax, dir1, dir2)
!
!
!     /usage/ call ctransVectZx(init, zax, xax, dir1, dir2)
!         Directions cosines(dir1) are given in a system
!         (=R).
!         The z axis of  R has dircetion cos zax in B.
!         The x axis of  R has direction cos xax in B.
!
!         This subroutine transform the angles so that (dir1)
!         be the direction cosines in the B-system,
!         and put the result into dir2.
!         dir2 can be the same one as dir1, or zax.
!    init:  if zax and xax  are different from those from the 
!           previous call, give 1  else  give a diff. value.
!    zax:  /coord/
!    xax:  /coord/
!    dir1: /coord/
!    dir2; /coord/
!  *** Note *** If you are ok with an arbitrary x axis, use ctransVectZ.
!
       implicit none

#include  "Zcoord.h"
       type(coord)::xax
       type(coord)::zax
       type(coord)::dir1
       type(coord)::dir2
       integer init
!
       type(coord)::yvec
       type(coord)::temp
       real*8 norm

#ifdef USESAVE
       save yvec
#endif

       character*70 msg
!            y -axis unit vector
       if(init .eq. 1) then
          yvec%r(1) = zax%r(2) * xax%r(3) - zax%r(3) * xax%r(2)
          yvec%r(2) = zax%r(3) * xax%r(1) - zax%r(1) * xax%r(3)
          yvec%r(3) = zax%r(1) * xax%r(2) - zax%r(2) * xax%r(1)
!           need not  normalize; but check it
          norm = yvec%r(1)**2 + yvec%r(2)**2 + yvec%r(3)**2
          if(abs(norm-1.0) .gt. 1.e-4) then
             write(msg, *)
     *             'ctransVectZx: input dir. cos. are not orthogonal'
             call cerrorMsg(msg, 0)
          endif
       endif
!
       temp%r(1) = dir1%r(1) * xax%r(1) + dir1%r(2) *yvec%r(1) +
     *     dir1%r(3) *zax%r(1)
       temp%r(2) = dir1%r(1) * xax%r(2) + dir1%r(2) *yvec%r(2) + 
     *     dir1%r(3) *zax%r(2)
       temp%r(3) = dir1%r(1) * xax%r(3) + dir1%r(2) *yvec%r(3) + 
     *     dir1%r(3) *zax%r(3)
!
                ! @jaxa, temp undefined warning comes out
                ! since temp.sys is not used.
!       dir2 = temp
       dir2%r(:) =temp%r(:)
       end
!      ===================================================
       subroutine ciTransVectZx(init, zax, xax, dir1, dir2)
!      ===================================================
!   This is the inverse of ctransVectZx
!
!     /usage/ call ciTransVectZx(init, zax, xax, dir1, dir2)
!   Suppose 2 system B and R.  
!         Directions cosines(dir1) are given in B
!
!         The z axis of  R has dircetion cos zax in B.
!         The x axis of  R has direction cos xax in B.
!
!         This subroutine transform the angles so that (dir1)
!         be the direction cosines in the R-system,
!         and put the result into dir2.
!         dir2 can be the same one as dir1, or zax.
!    init:  if zax and xax  are different from those from the preious call, give 1
!           else  give a diff. value.
!    zax:  /coord/
!    xax:  /coord/
!    dir1: /coord/
!    dir2; /coord/
!
       implicit none
#include  "Zcoord.h"
       type(coord)::xax
       type(coord)::zax
       type(coord)::dir1
       type(coord)::dir2
       integer init
!
       type(coord)::yvec
       type(coord)::temp
       real*8 norm
#ifdef USESAVE
       save yvec
#endif
       character*70 msg
!            y -axis unit vector
       if(init .eq. 1) then
          yvec%r(1) = zax%r(2) * xax%r(3) - zax%r(3) * xax%r(2)
          yvec%r(2) = zax%r(3) * xax%r(1) - zax%r(1) * xax%r(3)
          yvec%r(3) = zax%r(1) * xax%r(2) - zax%r(2) * xax%r(1)
!           need not  normalize; but check it
          norm = yvec%r(1)**2 + yvec%r(2)**2 + yvec%r(3)**2
          if(abs(norm-1.0) .gt. 1.e-4) then
             write(msg, *)
     *             'ctransVectZx: input dir. cos. are not orthogonal'
             call cerrorMsg(msg, 0)
          endif
       endif
!
       temp%r(1) = dir1%r(1) * xax%r(1) + dir1%r(2) *xax%r(2) +
     *  dir1%r(3)  *xax%r(3)
       temp%r(2) = dir1%r(1) * yvec%r(1)+ dir1%r(2) *yvec%r(2) + 
     *  dir1%r(3) *yvec%r(3)
       temp%r(3) = dir1%r(1) * zax%r(1) + dir1%r(2) *zax%r(2) +
     *    dir1%r(3)  *zax%r(3)
       dir2 = temp
       end

             
