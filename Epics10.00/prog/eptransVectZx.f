#include "ZsaveStruc.h"
!cc            test eptransVectZx,  epiTransVectZx
!       implicit none
!       include 'ZepDirec.h'
!       record /epDirec/ wa, wb, dir1, dir2
!       real*8 norm
! 
! 
! 
!       wa.x = .50
!       wa.y = 0.2
!       wa.z = sqrt(1. - (wa.x**2+wa.y**2))
!
!       wb.x = wa.y
!       wb.y = -wa.x
!       wb.z = 0.
!       norm = sqrt( wb.x**2 + wb.y**2 + wb.z**2)
!       wb.x = wb.x/norm
!       wb.y = wb.y/norm
!       wb.z = wb.z/norm
!
!
!       dir1.x = sqrt(2.)/3.
!       dir1.y = dir1.x*1.1
!       dir1.z = sqrt(1. - dir1.x**2 - dir1.y**2)
!       write(*,*) dir1.x, dir1.y, dir1.z
!       call eptransVectZx(1, wa, wb, dir1, dir2)
!       write(*,*) dir2.x, dir2.y, dir2.z
!       call ciTransVectZx(1, wa, wb, dir2, dir2)
!       write(*,*) dir2.x, dir2.y, dir2.z
!       call eptransVectZx(2, wa, wb, dir1, dir1)
!       write(*,*) dir1.x, dir1.y, dir1.z
!       end
!-------------------------------------------------------------
       subroutine eptransVectZx(init, zax, xax, dir1, dir2)
!
!
!     /usage/ call eptransVectZx(init, zax, xax, dir1, dir2)
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
!  *** Note *** If you are ok with an arbitrary x axis, use eptransVect.
!
       implicit none

#include  "ZepDirec.h"
       type(epDirec)::  xax, zax, dir1, dir2
       integer init
!
       type(epDirec)::  yvec, temp
       real*8 norm

#ifdef USESAVE
       save yvec
#endif

       character*70 msg
!            y -axis unit vector
       if(init .eq. 1) then
          yvec%x = zax%y * xax%z - zax%z * xax%y
          yvec%y = zax%z * xax%x - zax%x * xax%z
          yvec%z = zax%x * xax%y - zax%y * xax%x
!           need not  normalize; but check it
          norm = yvec%x**2 + yvec%y**2 + yvec%z**2
          if(abs(norm-1.0) .gt. 1.e-4) then
             write(0, *)
     *             'eptransVectZx: input dir. cos. are not orthogonal'
             write(0,*) ' norm =',norm
             write(0,*) 'xax%x y,z =', xax%x, xax%y, xax%z
             write(0,*) 'zax%x y,z =', zax%x, zax%y, zax%z
             call epfordebug('eptransVecZx')
             stop
          endif
       endif
!
       temp%x = dir1%x * xax%x + dir1%y *yvec%x +
     *     dir1%z *zax%x
       temp%y = dir1%x * xax%y + dir1%y *yvec%y + 
     *     dir1%z *zax%y
       temp%z = dir1%x * xax%z + dir1%y *yvec%z + 
     *     dir1%z *zax%z
!
       dir2 = temp
       end
!      ===================================================
       subroutine epiTransVectZx(init, zax, xax, dir1, dir2)
!      ===================================================
!   This is the inverse of eptransVectZx
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
!    zax:  /epDirec/
!    xax:  /epDirec/
!    dir1: /epDirec/
!    dir2; /epDirec/
!
       implicit none
#include  "ZepDirec.h"
       type(epDirec)::  xax, zax, dir1, dir2
       integer init
!
       type(epDirec)::  yvec, temp
       real*8 norm
#ifdef USESAVE
       save yvec
#endif
       character*70 msg
!            y -axis unit vector
       if(init .eq. 1) then
          yvec%x = zax%y * xax%z - zax%z * xax%y
          yvec%y = zax%z * xax%x - zax%x * xax%z
          yvec%z = zax%x * xax%y - zax%y * xax%x
!           need not  normalize; but check it
          norm = yvec%x**2 + yvec%y**2 + yvec%z**2
          if(abs(norm-1.0) .gt. 1.e-4) then
             write(0, *)
     *             'epitransVectZx: input dir. cos. are not orthogonal'
             write(0,*) ' norm =',norm
             write(0,*) 'xax%x y,z =', xax%x, xax%y, xax%z
             write(0,*) 'zax%x y,z =', zax%x, zax%y, zax%z
             call epfordebug('epitransVecZx')
             stop
          endif
       endif
!
       temp%x = dir1%x * xax%x + dir1%y *xax%y +
     *  dir1%z  *xax%z
       temp%y = dir1%x * yvec%x+ dir1%y *yvec%y + 
     *  dir1%z *yvec%z
       temp%z = dir1%x * zax%x + dir1%y *zax%y +
     *    dir1%z  *zax%z
       dir2 = temp
       end
!      implicit none
!      real(8)::Normal(3), Xnormal(3), xx
!c      Normal=(/0., 0, 1./)
!c      Normal=(/0., 1., 0./)
!c      Normal=(/1.,0., 0./)
!      Normal(1)=0.4d0
!      Normal(2) = sqrt(1.d0 - 0.4d0**2 -0.1d0**2)
!      Normal(3) = 0.1d0
!      call epgetNormal2vec(Normal, Xnormal)
!      write(0,*) Xnormal
!      call cscalerProd(Xnormal, Normal, xx)
!      write(0,*) xx
!      end
      subroutine epgetNormal2vec(vec, z)
      implicit none
!           get a  vector normal to a given vector
      real(8),intent(in)::vec(3)  ! give vector
      real(8),intent(out)::z(3)   ! one vector normal to vec
      
      real(8)::cosf,sinf, cost, sint
      real(8)::fai
      real(8)::yax(3,3), zax(3,3), rot(3,3), z0(3)
      
      fai = atan2(vec(2), vec(1))
      cosf = cos(fai)
      sinf = sin(fai)

      cost = vec(3)
      sint = sqrt(1.d0 - cost**2)
!         apply rot to (0,0,1)
      z0=(/0.d0, 0.d0, 1.d0/)
!      rotattion matrix around y by a=(pi/2 - theta)
!       cos(a) = sint   sin(a) = cost 
      call cgetRotMat3(2,  sint, cost, yax)
!        around z by -fai  
      call cgetRotMat3(3, cosf, -sinf, zax)
!          fai x  a      rot =   apply yax then zax
      call cmultRotMat3(zax, yax, rot)  
      call capplyRot3(rot, z0, z)
      end

      
