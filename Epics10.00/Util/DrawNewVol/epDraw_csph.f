      subroutine epDraw_csph(comp, p, n)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                               !   a csph in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line
                               ! dimension of p must be >+ (nvccl+2)*2
      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  


      real*8  zmx, zmn, h

      integer ir, ih, irh
      parameter (ir = 1,  ih = 2, irh=3)
      integer nvtxx, n2
      logical kdgtest

      h  =  Volat( comp%vol+ih) 
      zmx = min( Volat( comp%vol+ir )*cos(pamin*Torad),
     *            Volat( comp%vol + ir))
      zmn = max(Volat( comp%vol+ir )*cos(pamax*Torad), h)

      call epdrawSphere(Volat( comp%vol + ir),
     *  zmn, zmx, p, n)
      if( zmn .eq. h ) then
         if(kdgtest(howcyl, 2) ) then
!               inquire # of vertexes.
            call epqsphere(nvtxx)
!                  draw edge
            call epdrawCylEdg(p(1), nvtxx, h, p(n+1), n2)
            n = n + n2
            n = n + 1
            p(n)%x = gpsep
         endif
      endif
      end
      subroutine epDraw_dcsph(comp, p, n)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                               !   a csph in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line
                               ! dimension of p must be >+ (nvccl+2)*2
      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  



       integer ir, ih1, ih2, irh1, irh2
       parameter (ir = 1,  ih1 = 2, irh1=3, ih2=4, irh2=5 )

       real*8 r, h1, h2
       real*8  zmx, zmn
       integer nvtxx, n2, nsv
       logical kdgtest

       h1 = Volat( comp%vol + ih1)
       h2 = Volat( comp%vol + ih2)
       zmx = min( Volat( comp%vol+ir )*cos(pamin*Torad),  h2)
       zmn = max( Volat( comp%vol+ir )*cos(pamax*Torad),  h1)

       r = Volat( comp%vol + ir)
       call epdrawSphere(r,  zmn, zmx,   p, n)
       nsv = n

       if( zmn .eq. h1 ) then
          if(kdgtest(howcyl, 2) ) then
!               inquire # of vertexes.
             call epqsphere(nvtxx)

!                  draw edge
             call epdrawCylEdg(p(1), nvtxx, h1, p(n+1), n2)
             n = n + n2
             n = n + 1
             p(n)%x = gpsep
          endif
       endif
       if( zmx .eq. h2 ) then
          if(kdgtest(howcyl, 1) ) then
!               inquire # of vertexes.
             call epqsphere(nvtxx)
!                  draw edge
             call epdrawCylEdg(p(nsv-nvtxx), nvtxx, h2, p(n+1), n2)
             n = n + n2
             n = n + 1
             p(n)%x = gpsep
          endif          
       endif
       end
