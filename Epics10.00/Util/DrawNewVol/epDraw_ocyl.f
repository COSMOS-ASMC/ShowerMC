      subroutine epDraw_ocyl(comp, p, n)
      implicit none
!          this version neglect thetamin, max
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                               !   a ocyl in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line
                               ! dimension of p must be >+ (nvccl+2)*2
      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  

       integer ir, ih, isa, iea,
     *         ircossa, irsinsa, ircosea, irsinea
       parameter( ir = 1,  ih = 2,  isa=3, iea=4, 
     *            ircossa=5, irsinsa=6, ircosea=7, irsinea=8 )

       real*8 r, h, sa, ea
       real*8 t1, t2
      integer n1, n2, nsv1, nsv2
      logical kdgtest
      real*8 x

       r = Volat( comp%vol + ir)
       h = Volat( comp%vol + ih)
       sa= Volat( comp%vol + isa)
       ea= Volat( comp%vol + iea)

       n = 0
       if(ea < sa) then
          ea = ea + 360.
       endif
       t1 = sa
       t2 = ea
       call epdrawCcl(r,  0.d0, t1, t2, p(n+1), n1)
       nsv1 = n+1
       n = n + n1
       call epdrawCcl(r,  h, t1, t2, p(n+1), n2)
       nsv2 = n + 1
       n =n + n2
       n = n + 1
       p(n)%x = gpsep
       if(kdgtest(howcyl, 1)) then
          call epdrawCylEdg(p(nsv1), n1, 0.d0, p(n+1), n2)
          n = n + n2 
       endif
       if(kdgtest(howcyl, 2)) then
          call epdrawCylEdg( p(nsv2), n1, h, p(n+1), n2)
          n = n + n2
       endif
!       if(scut) then
!           draw startng cut  square
         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = Volat( comp%vol + ircossa)
         p(n)%y = Volat( comp%vol + irsinsa)
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = h

         n = n + 1
         p(n)%x = Volat( comp%vol + ircossa)
         p(n)%y = Volat( comp%vol + irsinsa)
         p(n)%z = h

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
!      endif
!      if(ecut) then
!           draw endign cut  square
         n = n + 1
         p(n)%x = Volat( comp%vol + ircosea)
         p(n)%y = Volat( comp%vol + irsinea)
         p(n)%z = 0.

         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = Volat( comp%vol + ircosea)
         p(n)%y = Volat( comp%vol + irsinea)
         p(n)%z = h

         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = h

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!     endif

      end
