!     *******************
      subroutine epDraw_hypb(comp, p, n)
      implicit none

#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                               !   hypb in local coordnate.

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  

       integer ia, ib, ig, ih1, ih2
       parameter( ia = 1,  ib = 2,  ig=3, ih1=4, ih2= 5)

       real*8 a, b, g, h1,  h2
       integer n1, n2, nt
       logical kdgtest

       real*8  h, dh, hg, aa, bb, sq, xdp, hp

!
       a = Volat( comp%vol + ia)
       b = Volat( comp%vol + ib)
       g = Volat( comp%vol + ig)
       h1 = Volat( comp%vol + ih1)
       h2 = Volat( comp%vol + ih2)

      
       h = h1
       n = 0
       nt = 0
       hg = (h/g)**2
       sq = sqrt(1.+ hg)
       aa = a *  sq
       bb = b *  sq
       call epdrawElps(aa, bb, h, thetamax, thetamin, p(n+1), n1)
       n = n + n1
       do while (.true.)
!          xdp = x''(z=h)
!             0 < xdp < a/g**2; for larger xdp, smaller step of dh
          xdp = a/g**2 / sq * (1.0-hg/(1.0+hg))
          dh =max(0.5/xdp, (h2-h1)/20.)
          if( h .lt. 0. .and. h + dh + dh/4. .gt. 0.) then
             hp = min(0.d0, h2)
          else
             hp = min(h+dh, h2)
          endif
          hg = (hp/g)**2
          xdp = (xdp +  a/g**2 / sq * (1.0-hg/(1.0+hg)))/2.
          dh =max( 0.5/xdp, (h2-h1)/20.)
          if( h .lt. 0. .and. h + dh + dh/4. .gt. 0.) then
             h = min(0.d0, h2)
          else
             h = min(h+dh, h2)
          endif
          if(h .ne. 0. .and. h + dh/4 .gt. h2) then
             h = h2
          endif

          hg = (h/g)**2
          sq = sqrt(1.+ hg)
          aa = a *  sq
          bb = b *  sq
          call epdrawElps(aa, bb, h, thetamax, thetamin, p(n+1), n1)
          n = n + n1
          if(h .eq. h2) goto 10
      enddo
 10   continue
      nt = n - n1
      n = n + 1
      p(n)%x = gpsep
!   
      if(kdgtest(howcyl, 1)) then
         call epdrawCylEdg(p, n1, h1, p(n+1), n2)
         n = n + n2
!         n = n + 1     ! extra 3rd line; not needed
!         p(n).x = gpsep
      endif

      if(kdgtest(howcyl, 2)) then
         call epdrawCylEdg(p(nt+1), n1, h2, p(n+1), n2)
         n = n + n2
!         n = n + 1      ! extra 3rd line; not needed
!         p(n).x = gpsep
      endif
      end
