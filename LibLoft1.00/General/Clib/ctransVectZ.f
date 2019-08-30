!            test ctransVectZ
!      implicit none
!      include 'Zcoord.h'
!      record /coord/ wa, dc, ans
!      real*8 cst, snt, fai
!     wa.r(1)=.50
!      wa.r(2)=0.5
!      wa.r(3)=sqrt(1. - (wa.r(1)**2+wa.r(2)**2))
!      cst=0.8
!      snt=sqrt(1.-cst**2)
!      do fai=0., 2*3.1415, .1
!         dc.r(1)=snt*cos(fai)
!         dc.r(2)=snt*sin(fai)
!         dc.r(3)=cst
!         call ctransVextZ(wa, dc, ans)
!         write(*,*) sngl(ans.r(1)), sngl(ans.r(2)),
!    *    sngl(ans.r(3))
!      enddo
!      end
       subroutine ctransVectZ(zax, dir1, dir2)
       implicit none

#include  "Zcoord.h"
       type(coord)::zax
       type(coord)::dir1
       type(coord)::dir2
       real*8  sml, epsx, av

       parameter (sml=0.001d0,  epsx=1.d-8, av=.985d0)

!
!   /usage/ call ctransVectZ(zax, dir1, dir2)
!         Directions cosines(dir1) are given in a system
!         (=R) whose z-axis has direction cosines (zax) in
!         a certain system(=B).
!         This subroutine transform the angles so that (dir1)
!         be the direction cosines in the B-system,
!         and put the result into dir2.
!         The x and y
!         axes of the R-system are chosen so that the transformation
!         becomes simplest. This does not guarantee that the dir2
!         have the same sing as the original one when xax(1) is
!         1.0 or close to 1.0.   If you have to avoid such, 
!         use, ctransVectZ2. (For magnetic deflection, you need this).
!      dir2 can be the same one as dir1, or zax.
!      dir1 need not be  the direction cosine, but can be momentum 
!          or arbitrary  vector.  zax must be direction cos.
!
       real*8 w1a, w2a, w3a, dc1, dc2, dc3, el2, em2, d, a, b, c
       real*8 tmpa, tmpb, tmpc
!       real*8 norm
          w1a=zax%r(1)
          w2a=zax%r(2)
          w3a=zax%r(3)
          dc1=dir1%r(1)
          dc2=dir1%r(2)
          dc3=dir1%r(3)
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
!                  check result
!              norm = tmpa**2 + tmpb**2 + tmpc**2 
!              if(abs(norm-1.d0) .gt. 1.e-5) then
!                     renormalize
!                   norm = sqrt(norm)
!                   tmpa=tmpa/norm
!                   tmpb=tmpb/norm
!                   tmpc=tmpc/norm
!              endif
          dir2%r(1)=tmpa
          dir2%r(2)=tmpb
          dir2%r(3)=tmpc
        end
