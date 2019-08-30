!     *****************************
      subroutine epDraw_cecyl(comp, p, n)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                              !   a cecyl in local coordnate.
                              ! (x,y,z)= gpsep is a separator
                              ! to be converted to a blank line
                              ! dimension of p must be >+ (nvccl+2)*2
      integer  n              ! output.  number of (x,y,z) data
                              ! put in p.  

       integer ira, irb, ih, isa, iea,
     * ix0, iy0, ix1, iy1
       parameter( ira = 1,  irb=2, ih = 3,  isa=4, iea=5,
     *   ix0=6, iy0=7, ix1=8, iy1=9)

       real*8 ra, rb, h, sa, ea


       real*8 t1, t2
       type(epPos)::  mp

      integer n1, n2, nsv
      logical kdgtest
      logical isinside
      logical isinside2, scut, ecut
      real*8 x
      isinside(x) = mod(thetamin-thetamax+360.d0, 360.d0) .ge.
     *               mod(x-thetamax+360.d0, 360.d0)
      isinside2(x) = mod(ea-sa+360.d0, 360.d0) .ge.
     *               mod(x-sa+360.d0, 360.d0)

       ra = Volat( comp%vol + ira)
       rb = Volat( comp%vol + irb)
       h = Volat( comp%vol + ih)
       sa= Volat( comp%vol + isa)
       ea= Volat( comp%vol + iea)
       scut = .false.
       ecut = .false.
       mp%x = ( Volat( comp%vol + ix0) + Volat( comp%vol + ix1) )/2.
       mp%y = ( Volat( comp%vol + iy0) + Volat( comp%vol + iy1) )/2.
       

       n = 0
       if(.not. isinside(sa) .and. isinside(ea) ) then
          t1 = thetamax
          t2 = ea
          ecut = .true.
       elseif(isinside(sa) .and. isinside(ea)) then
          if(isinside2(thetamin)) then
!              there is two regions
             scut = .true.
             t1 =  sa
             t2 = mod(thetamin,360.d0)
             if(t2 .lt. t1) t2 = t2+360.d0
             call epdrawElps(ra, rb, 0.d0, t1, t2, p(n+1), n1)
             nsv = n+1
             n = n + n1
             call epdrawElps(ra, rb, h, t1, t2, p(n+1), n2)
             n =n + n2
             n = n + 1
             p(n)%x = gpsep
             if(kdgtest(howcyl, 1)) then
                call epdrawCcylEdg(p(nsv), n1, 0.d0,  mp, p(n+1), n2)
                n = n + n2 
             endif
             if(kdgtest(howcyl, 2)) then
                call epdrawCcylEdg(p(nsv+n1), n1, h,  mp, p(n+1), n2)
                n = n + n2
                n = n + 1
                p(n)%x = gpsep
             endif


             t1 =mod(thetamax, 360.d0)
             t2 = ea
             if(t2 .lt. t1) t2 = t2+360.d0
             call epdrawElps(ra, rb, 0.d0, t1, t2, p(n+1), n1)
             nsv = n+1
             n = n +  n1
             call epdrawElps(ra, rb, h, t1, t2, p(n + 1 ), n2)
             n =n+ n2
             n = n + 1
             p(n)%x = gpsep
             ecut = .true.
             if(kdgtest(howcyl, 1)) then
                call epdrawCcylEdg(p(nsv),  n1, 0.d0, mp, p(n+1), n2)
                n = n + n2 
             endif
             if(kdgtest(howcyl, 2)) then
                call epdrawCcylEdg(p(nsv+n1),  n1, h, mp, p(n+1), n2)
                n = n + n2
                n = n + 1
                p(n)%x = gpsep
             endif
             goto 100
          else
             t1 = sa
             t2 = ea
             scut = .true.
             ecut = .true.
          endif
       elseif(isinside2(thetamax) .and. isinside2(thetamin)) then
          t1 =mod(thetamax, 360.d0)
          t2 =mod(thetamin, 360.d0)
       elseif(isinside2(thetamin) .and. isinside(sa) ) then
          t1 = sa
          t2 = mod(thetamin, 360.d0)
          scut = .true.
       else
          return  !  no part to be drawn
       endif
       if(t2 .lt. t1) t2=t2+360.d0
       call epdrawElps(ra, rb, 0.d0, t1, t2, p(n+1), n1)
       nsv = n+1
       n = n + n1
       call epdrawElps(ra, rb, h, t1, t2, p(n+1), n2)
       n = n +  n2
       n = n + 1
       p(n)%x = gpsep
 100   continue
       if(kdgtest(howcyl, 1)) then
          call epdrawCcylEdg(p(nsv),  n1, 0.d0, mp, p(n+1), n2)
          n = n + n2 
       endif
       if(kdgtest(howcyl, 2)) then
          call epdrawCcylEdg( p(nsv+n1), n1, h,  mp, p(n+1), n2)
          n = n + n2
          n = n + 1
          p(n)%x = gpsep
       endif


       if(scut) then
!           draw startng cut  square
         n = n + 1
         p(n)%x = mp%x
         p(n)%y = mp%y
         p(n)%z = 0.

         n = n + 1
         p(n)%x = Volat( comp%vol + ix0)
         p(n)%y = Volat( comp%vol + iy0)
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = mp%x 
         p(n)%y = mp%y
         p(n)%z = h

         n = n + 1
         p(n)%x = p(n-3)%x
         p(n)%y = p(n-3)%y
         p(n)%z = h

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif
      if(ecut) then
!           draw endign cut  square
         n = n + 1
         p(n)%x = Volat( comp%vol + ix1)
         p(n)%y = Volat( comp%vol + iy1)
         p(n)%z = 0.

         n = n + 1
         p(n)%x = mp%x
         p(n)%y = mp%y
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = p(n-3)%x
         p(n)%y = p(n-3)%y
         p(n)%z = h

         n = n + 1
         p(n)%x = mp%x
         p(n)%y = mp%y
         p(n)%z = h

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif
      end
