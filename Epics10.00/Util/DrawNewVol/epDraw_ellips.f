!     *****************************
      subroutine epDraw_ellips(comp, p, n)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                              !   an ellipse in local coord.
                              ! (x,y,z)= gpsep is a separator
                              ! to be converted to a blank line
                              ! dimension of p must be >+ (nvccl+2)*2
      integer  n              ! output.  number of (x,y,z) data
                              ! put in p.


      integer ia, ib, ig, ih1, ih2
      parameter( ia = 1,  ib = 2,  ig=3, ih1=4, ih2= 5)

      real*8 a, b, g, h1,  h2
      integer n1, n2, nsv1, nsv2
      logical kdgtest

      real*8  h, dh,  aa, bb, sq, teta, dteta,  hpi, hmin, hmax

!
      hpi = pi/2.d0
      a = Volat( comp%vol + ia)
      b = Volat( comp%vol + ib)
      g = Volat( comp%vol + ig)
      h1 = Volat( comp%vol + ih1)
      h2 = Volat( comp%vol + ih2)

      
      hmin = max(h1, g*cos(pamax*Torad))
      hmax = min(h2, g*cos(pamin*Torad))
      h = hmin
      n = 0
      sq = sqrt(1.d0 - (h/g)**2)

      aa = a *  sq
      bb = b *  sq

      call epdrawElps(aa, bb, h, thetamax, thetamin, p(n+1), n1)
      nsv1 = n + 1
      n = n + n1
      teta =asin(h/g)
      dteta = pi * 15./180.
      do while (.true.)
         dh =  g * (sin(min(teta+dteta, hpi) )- sin(teta) )
         if(dh .eq. 0.) goto 10
         if( h .lt. 0. .and. h + dh + dh/4. .gt. 0.) then
            h = min(0.d0, hmax)
         else
            h = min(h+dh, hmax)
         endif
         sq = sqrt(1.d0 - (h/g)**2)
         aa = a *  sq
         bb = b *  sq
         call epdrawElps(aa, bb, h, thetamax, thetamin, p(n+1), n1)
         nsv2 = n+1
         n = n + n1
         if(h .eq. hmax .or. teta .eq. hpi) goto 10
         teta= min(teta + dteta, hpi)
      enddo
 10   continue
      n = n + 1
      p(n)%x = gpsep
!   
      if(kdgtest(howcyl, 1) .and. h1 .gt. g*cos(pamax*Torad)) then
         call epdrawCylEdg(p(nsv1), n1, h1, p(n+1), n2)
         n = n + n2
         n = n + 1
         p(n)%x = gpsep
      endif

      if(kdgtest(howcyl, 2) .and. h2 .lt. g*cos(pamin*Torad) ) then
         call epdrawCylEdg(p(nsv2), n1, h2, p(n+1), n2)
         n = n + n2
         n = n + 1
         p(n)%x = gpsep
      endif
      end

