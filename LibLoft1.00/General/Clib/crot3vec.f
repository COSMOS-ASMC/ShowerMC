!c     test crot3vec, cirot3vec
!      implicit none
!#include "Zptcl.h"
!        
!      integer i
!      type(fmom):: wa, dc, ans
!      real*8 cst, snt, fai
!      wa.p(1)=0.1
!      wa.p(2)=-0.80
!      wa.p(3)=-sqrt(1. - (wa.p(1)**2+wa.p(2)**2))
!      
!      wa.p(1:3) = wa.p(1:3)*10.
!      cst=0.35d0
!      snt=sqrt(1.-cst**2)
!      i = 1
!      do fai=0., 2*3.1415, .1
!
!         dc.p(1)=snt*cos(fai)*10.
!         dc.p(2)=snt*sin(fai)*10.
!         dc.p(3)=cst*10.
!         write(*,'(a,1p, 4g15.5)')
!     *     'bef  crot',  dc.p(1:3),  dc.p(1)**2 + dc.p(2)**2
!         call crot3vec(wa, dc, ans)
!         write(*,'(a,1p, 4g15.5)')
!     *     'aft crot',  ans.p(1:3), ans.p(1)**2 + ans.p(2)**2
!         call cirot3vec(i, wa, ans, ans)
!         write(*,'(a, 1p, 4g15.5)') 
!     *     'aft cirot', ans.p(1:3), ans.p(1)**2 + ans.p(2)**2
!       enddo
!       end
!      ************************************
       subroutine crot3vec(zax, vec1, vec2)
!      ************************************
       implicit none
#include "Zptcl.h"
       type(fmom):: zax, vec1, vec2
!   dangerous.  vec2 dose not recieve 4th coponent !!!
       real*8  sml, epsx, av

       parameter (sml=0.001,  epsx=1.e-4, av=.985)

!
!   /usage/ call crot3vec(zax, vec1, vec2)
!         4 momentum  vec1 are given in a system
!         (=R) whose z-axis has pallarel to 4 momentum (zax) in
!         a certain system(=B).
!         This subroutine transform the angles so that (vec1)
!         be the 4 momeuntum in the B-system,
!         and put the result into vec2.
!   (Note 4th component is acutally not needed and not changed)
!         The x and y
!         axes of the R-system are chosen so that the transformation
!         becomes simplest. This does not guarantee that the vec2
!         have the same sing as the original one when xax(1) is
!         1.0 or close to 1.0. 
!      vec2 can be the same one as vec1, or zax.
!   This is essentially the same as ctransVectZ.
!
       real*8  el2, em2, d, a, b, c
       real*8 tmpa, tmpb, tmpc, pabs
       type(fmom):: zaxdir
!
!           get direction cos from zax
       pabs = sqrt(zax%p(1)**2 + zax%p(2)**2 + zax%p(3)**2)
       if(pabs .eq. 0.) then
          zaxdir%p(1) = 0.
          zaxdir%p(2) = 0.
          zaxdir%p(3) = 1.
       else
          zaxdir%p(1) = zax%p(1)/pabs
          zaxdir%p(2) = zax%p(2)/pabs
          zaxdir%p(3) = zax%p(3)/pabs
       endif
!
       el2=zaxdir%p(1)**2
       em2=zaxdir%p(2)**2
       d=1.+zaxdir%p(3)
       if(abs(d) .gt. epsx) then
          a=el2/d - 1.
          b=zaxdir%p(1)*zaxdir%p(2)/d
          c=em2/d - 1.
          tmpa=a*vec1%p(1) + b*vec1%p(2) + zaxdir%p(1)*vec1%p(3)
          tmpb=b*vec1%p(1) + c*vec1%p(2) + zaxdir%p(2)*vec1%p(3)
       else
          tmpa= vec1%p(2)
          tmpb= vec1%p(1)
       endif
       tmpc=zaxdir%p(1)*vec1%p(1) + zaxdir%p(2)*vec1%p(2) +
     *     zaxdir%p(3)*vec1%p(3)
       vec2%p(1) = tmpa
       vec2%p(2) = tmpb
       vec2%p(3) = tmpc
       end

!    ***************************************
      subroutine cirot3vec(init, p1, p2, po)
      implicit none
#include "Zptcl.h"
!    ***************************************
!     
!         suppose p1 and p2 are given in the same system.
!      This program rotate the  z-axis so that
!      it coinsides with p1 and transform p2 as seen from p1.
!      The result is put in po.
!      Note that the x and y are arbitrary chosen.
!      the 4-th  comp. is not used and unchanged.
!
      integer init  ! input. give 1 if p1 is diff. from prev. call.
      type(fmom)::p1  ! input.
      type(fmom)::p2  ! input. 
      type(fmom)::po  ! output. po may be p2
!
      type(fmom)::x, y, px
      real*8 xnorm, ynorm, pnorm, sum, ax, ay, az

!c #ifdef   USESAVE
      save x, y, xnorm, ynorm, pnorm
!c #endif
      
!
!        form x and y with  p1 being z.
!
      if(init .eq. 1) then
         pnorm = sqrt(p1%p(1)**2 + p1%p(2)**2 + p1%p(3)**2)
      endif
      if( pnorm .eq. 0.0 ) then
         po = p2
      else
         if(init .eq. 1) then 
            ax = abs(p1%p(1))
            ay = abs(p1%p(2))
            az = abs(p1%p(3))

            if(az .gt. ay) then
               if(az .gt. ax) then
                  sum = -(p1%p(1) + p1%p(2)) 
                  x%p(1) = 1.
                  x%p(2) = 1.
                  x%p(3) = sum/p1%p(3)
               else
                  sum = -(p1%p(2) + p1%p(3)) 
                  x%p(2) = 1.
                  x%p(3) = 1.
                  x%p(1) = sum/p1%p(1)
               endif
            elseif(ay .gt. ax ) then
               sum = -(p1%p(1) + p1%p(3)) 
               x%p(1) = 1.
               x%p(3) = 1.
               x%p(2) = sum/p1%p(2)
            else
               sum = -(p1%p(2) + p1%p(3)) 
               x%p(2) = 1.
               x%p(3) = 1.
               x%p(1) = sum/p1%p(1)
            endif
!             form y
            y%p(1) = p1%p(2)* x%p(3) - p1%p(3) * x%p(2)
            y%p(2) = p1%p(3)* x%p(1) - p1%p(1) * x%p(3)
            y%p(3) = p1%p(1)* x%p(2) - p1%p(2) * x%p(1)
!
!            since p1, x, y are not unit vector, we have to normalize
!
            xnorm = sqrt(x%p(1)**2 + x%p(2)**2 + x%p(3)**2)
            ynorm = sqrt(y%p(1)**2 + y%p(2)**2 + y%p(3)**2)
         endif
!
!          if init != 1, we come here directly
!
!          get projection of p2 on to x, y, p1
!
         px%p(1) = (p2%p(1)* x%p(1) + p2%p(2)* x%p(2) +
     *        p2%p(3) * x%p(3))/xnorm
         px%p(2) = (p2%p(1)* y%p(1) + p2%p(2)* y%p(2) +
     *        p2%p(3) * y%p(3))/ynorm
         px%p(3) = (p2%p(1)* p1%p(1) + p2%p(2)* p1%p(2) +
     *        p2%p(3) * p1%p(3))/pnorm
!
         po%p(1:3) = px%p(1:3)
         po%p(4) = p2%p(4)
      endif
      end






      

