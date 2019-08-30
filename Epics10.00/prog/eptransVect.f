!c            test eptransVect
!       implicit none
!#include "ZepDirec.h"
!
!       record /epDirec/ wa, dc, ans1, ans2
!
!       real*8 cst, snt, fai, dfai
!       real(8),parameter::pi=asin(1.0d0)*2.0d0
!       integer,save:: init
!!      change next 3 for test 
!       wa.x= 0.
!       wa.y= -1.d-7
!       wa.z= -sqrt(1.d0-wa.x**2-wa.y**2)
!       
!       write(*,'(a, 1p,3g18.12)') 'wa=',wa.x, wa.y, wa.z
!       cst=wa.z
!       snt=sqrt(1.d0-cst**2)
!!               Z        
!!               |   wa     dc        
!!               |    *      %
!!               |  *   % 
!!               *% -------------->  Y
!!              /
!!             /
!!            X
!!     wa in (X,Y,Z)  system
!!     dc in wa system where wa is z-direction
!!     we get dc in (X,Y,Z) system.
!!       
!      init = 1 
!      do dfai=0.d0, 360.d0, 10.d0
!          fai = dfai/180.d0*pi
!          dc.x=snt*cos(fai)
!          dc.y=snt*sin(fai)
!          dc.z=cst
!          write(*,*) 'dfai=', dfai
!          write(*,'(a, 1p, 3g15.6)') "local  ",  dc.x, dc.y, dc.z
!          call eptransVect(wa, dc, ans1)
!          write(*,'(a, 1p,3g15.6)') "world  ",
!     *         ans1.x, ans1.y, ans1.z
!          call epitransVect(init, wa, ans1, ans2)
!          write(*,'(a, 1p,3g15.6)') "w->loc ",  ans2.x, ans2.y, ans2.z
!          call epitransVect(0, wa, ans1, ans2)
!          write(*,'(a, 1p,3g15.6)') " again ",  ans2.x, ans2.y, ans2.z
!          init = 2
!       enddo
!       end
       subroutine eptransVect(zax, dir1, dir2)
       implicit none
#include "ZepDirec.h"
       type(epDirec)::   zax   ! input
       type(epDirec)::   dir1  ! input
       type(epDirec)::   dir2  ! output. can dir1 / zax

       real*8   epsx

       parameter ( epsx=1.d-10 )

!
!         Directions cosines(dir1) are given in a system
!         (=R) whose z-axis has direction cosines (zax) in
!         a certain system(=B).
!         This subroutine transform the angles so that (dir1)
!         be the direction cosines in the B-system,
!         and put the result into dir2.
!         The x and y
!         axes of the R-sysetm are chosen so that the transformation
!         becomes simplest. This does not guarantee that the dir2
!         have the same sing as the original one when zax.z is
!         1.0 or close to 1.0.   
!      dir2 can be the same one as dir1, or zax.
!      dir1 need not be  the direction cosine, but can be momentum 
!          or arbitrary  vector.  zax must be direction cos.
!
       real*8 w1a, w2a, w3a, dc1, dc2, dc3, el2, em2, d, a, b, c
       real*8 tmpa, tmpb, tmpc, temp

      w1a=zax%x
      w2a=zax%y
      w3a=zax%z
      dc1=dir1%x
      dc2=dir1%y
      dc3=dir1%z
!
      el2=w1a**2
      em2=w2a**2
      d=1.+w3a
      if(abs(d) .gt. epsx) then
         a=el2/d - 1.
         b=w1a*w2a/d
         c=em2/d - 1.
         tmpa=a*dc1 + b*dc2 + w1a*dc3
         tmpb=b*dc1 + c*dc2 + w2a*dc3
      else
         tmpa= dc2
         tmpb= dc1
      endif
      tmpc=w1a*dc1 + w2a*dc2 + w3a*dc3
      dir2%x=tmpa
      dir2%y=tmpb
      dir2%z=tmpc
      temp = dir2%x**2 + dir2%y**2 + dir2%z**2
      if(abs(dir2%x) .gt. 1.d0 .or. 
     *         abs(dir2%y) .gt. 1.d0 .or. 
     *         abs(dir2%z) .gt. 1.d0 .or. 
     *         abs(temp -1.d0) .gt. epsx) then
!           renormalize
         temp = sqrt(temp)
         if(temp .eq. 0.) then
            dir2%x = 0.d0
            dir2%y = 0.d0
            dir2%z = 1.d0
         else
            dir2%x = dir2%x/temp
            dir2%y = dir2%y/temp
            dir2%z = dir2%z/temp
         endif
      endif
      end
      subroutine epitransVect(init, zax, dir1, dir2)
       implicit none
#include "ZepDirec.h"
       integer,intent(in):: init ! give 1 if zax is diff. from prev. call
       type(epDirec)::   zax    ! input  in B
       type(epDirec)::   dir1  ! input   in B
       type(epDirec)::   dir2  ! output. can dir1 / zax  in R

       real*8   epsx

       parameter ( epsx=1.d-10 )

!
       real(8),save::  el2, em2, d
       real(8),save:: abc(3)

       if(init == 1 ) then
          el2 = zax%x**2
          em2 = zax%y**2
          d = 1.d0 + zax%z
          if(abs(d) > epsx) then
             abc(1)=el2/d - 1.d0
             abc(2)=zax%x*zax%y/d
             abc(3)=zax%x
          else
             abc(1) =0.
             abc(2) =1.d0
             abc(3) =0.
          endif
          call epitransVectZx(init, zax, abc, dir1, dir2)
       else
          call epitransVectZx(init, zax, abc, dir1, dir2)
       endif
       end
