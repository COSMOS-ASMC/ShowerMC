!     *******************
      subroutine epDraw_eparab(comp, p, n)
      implicit none

#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                               !   eparab in local coordnate.

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  

       integer ia, ib,  ih1, ih2
       parameter( ia = 1,  ib = 2,  ih1=3, ih2= 4)

       real*8 a, b,  h1,  h2


       logical kdgtest

       real*8  h, dh,  aa, bb, sq, hp, xa, ae, error
       integer nt, n1, n2
       data error/0.03d0/
!
       a = Volat( comp%vol + ia)
       b = Volat( comp%vol + ib)
       h1 = Volat( comp%vol + ih1)
       h2 = Volat( comp%vol + ih2)
       ae = min(a, b)
       h = h1
       if(h .eq. 0.) then
          h = min(h2*0.01d0, 0.01d0)
       endif
       n = 0
       nt = 0
       sq = sqrt(h)
       aa = a *  sq
       bb = b *  sq
       xa = sqrt(error)*2*ae
       call epdrawElps(aa, bb, h, thetamax, thetamin, p(n+1), n1)
       n = n + n1
       do while (.true.)
          dh =max(xa, (h2-h1)*0.1)
          if( h .lt. 0. .and. h + dh + dh/4. .gt. 0.) then
             hp = min(0.d0, h2)
          else
             hp = min(h+dh, h2)
          endif
          sq = sqrt(h)
          dh =max(xa, (h2-h1)*0.1)
          if( h .lt. 0. .and. h + dh + dh/4. .gt. 0.) then
             h = min(0.d0, h2)
          else
             h = min(h+dh, h2)
          endif
          if(h .ne. 0. .and. h + dh/4 .gt. h2) then
             h = h2
          endif
          sq = sqrt(h)
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
         n = n + 1
         p(n)%x = gpsep
      endif

      if(kdgtest(howcyl, 2)) then
         call epdrawCylEdg(p(nt+1), n1, h2, p(n+1), n2)
         n = n + n2
         n = n + 1
         p(n)%x = gpsep
      endif
      end

