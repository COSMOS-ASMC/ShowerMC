!     *****************************
      subroutine epDraw_ccone(comp, p, n)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                              !   a ccone in local coordnate.
                              ! (x,y,z)= gpsep is a separator
                              ! to be converted to a blank line
                              ! dimension of p must be >+ (nvccl+2)*2
      integer  n              ! output.  number of (x,y,z) data
                              ! put in p.  




       integer ira, irb, ih, irap, ik, isa, iea, ix0, iy0, ix1, iy1
       integer inx, iny, inz, ikp
       parameter (inx = 12, iny=13, inz=14, ikp=15)
       parameter( ira = 1,  irb = 2,  ih = 3, irap=4, ik=5,
     *      isa=6,  iea=7, ix0=8, iy0=9, ix1=10, iy1=11)


       real*8 ra, rb, h,  rap, sa, ea
       type(epPos)::  p1, p2, p3, p4 
       real*8   rbp, x1, y1, x2, y2


       real*8 t1, t2
       type(epPos)::  mp1, mp2


      integer n1, n2, nsv1
      logical kdgtest
      logical isinside
      logical isinside2, scut, ecut

      real*8 x, k
      isinside(x) = mod(thetamin-thetamax+360.d0, 360.d0) .ge.
     *               mod(x-thetamax+360.d0, 360.d0)
      isinside2(x) = mod(ea-sa+360.d0, 360.d0) .ge.
     *               mod(x-sa+360.d0, 360.d0)

!           check some values
       ra = Volat( comp%vol +  ira)
       rb = Volat( comp%vol +  irb)
       h = Volat( comp%vol +  ih)
       rap = Volat( comp%vol +  irap)
       sa= Volat( comp%vol +  isa)
       ea= Volat( comp%vol +  iea)
       k = Volat( comp%vol +  ik)
       rbp = rb*k
       x1 =  Volat( comp%vol +  ix0)
       y1 =  Volat( comp%vol +  iy0) 
       x2 =  Volat( comp%vol +  ix1)
       y2 =  Volat( comp%vol +  iy1)
!          define the cut plane
       p1%x =  x1
       p1%y =  y1
       p1%z =  0.
       p2%x =  p1%x*k
       p2%y =  p1%y*k
       p2%z =  h
       p3%x =  x2*k
       p3%y =  y2*k
       p3%z = h

       p4%x = x2
       p4%y = y2
       p4%z = 0.


       scut = .false.
       ecut = .false.
       mp1%x = ( x1 + x2)/2.
       mp1%y = ( y1 + y2)/2.
       mp1%z = 0.

       mp2%x = mp1%x*k
       mp2%y = mp1%y*k
       mp2%z = h

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
             nsv1 = n+1
             n = n + n1
             call epdrawElps(rap, rbp, h, t1, t2, p(n+1), n2)
             n =n + n2
             n = n + 1
             p(n)%x = gpsep
             if(kdgtest(howcyl, 1)) then
                call epdrawCcylEdg(p(nsv1), n1, 0.d0, mp1, p(n+1), n2)
                n = n + n2 
             endif
             if(kdgtest(howcyl, 2)) then
                call epdrawCcylEdg(p(nsv1+n1), n1, h, mp2, p(n+1), n2)
                n = n + n2
                n = n + 1
                p(n)%x = gpsep
             endif

             t1 =mod(thetamax, 360.d0)
             t2 = ea
             if(t2 .lt. t1) t2 = t2+360.d0
             call epdrawElps(ra, rb, 0.d0, t1, t2, p(n+1), n1)
             nsv1 = n+1
             n = n +  n1
             call epdrawElps(rap, rbp, h, t1, t2, p(n + 1 ), n2)
             n =n+ n2
             n = n + 1
             p(n)%x = gpsep
             ecut = .true.
             if(kdgtest(howcyl, 1)) then
                call epdrawCcylEdg(p(nsv1),  n1, 0.d0, mp1, p(n+1), n2)
                n = n + n2 
             endif
             if(kdgtest(howcyl, 2)) then
                call epdrawCcylEdg(p(nsv1+n1),  n1, h, mp2, p(n+1), n2)
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
       nsv1 = n+1
       n = n + n1
       call epdrawElps(rap, rbp, h, t1, t2, p(n+1), n2)
       n = n +  n2
       n = n + 1
       p(n)%x = gpsep
 100   continue
       if(kdgtest(howcyl, 1)) then
          call epdrawCcylEdg(p(nsv1),  n1, 0.d0, mp1, p(n+1), n2)
          n = n + n2 
       endif
       if(kdgtest(howcyl, 2)) then
          call epdrawCcylEdg( p(nsv1+n1), n1, h, mp2, p(n+1), n2)
          n = n + n2
          n = n + 1
          p(n)%x = gpsep
       endif


       if(scut) then
!           draw startng cut  square
         n = n + 1
         p(n) = mp1

         n = n + 1
         p(n) =  p1

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n) = mp2

         n = n + 1
         p(n) = p2

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif
      if(ecut) then
!           draw endign cut  square
         n = n + 1
         p(n) = p4

         n = n + 1
         p(n) = mp1

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n) = p3

         n = n + 1
         p(n) = mp2

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif
      end






