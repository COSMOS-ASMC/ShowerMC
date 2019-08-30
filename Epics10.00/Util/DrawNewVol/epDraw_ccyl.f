!     *****************************
      subroutine epDraw_ccyl(comp, p, n)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                              !   a ccyl in local coordnate.
                              ! (x,y,z)= gpsep is a separator
                              ! to be converted to a blank line
                              ! dimension of p must be >+ (nvccl+2)*2
      integer  n              ! output.  number of (x,y,z) data
                              ! put in p.  

       integer ir, ih, isa, iea,
     *         ircossa, irsinsa, ircosea, irsinea
       parameter( ir = 1,  ih = 2,  isa=3, iea=4, 
     *            ircossa=5, irsinsa=6, ircosea=7, irsinea=8 )

       real*8 r, h, sa, ea
       real*8 t1, t2
       type(epPos)::  mp
       real*8 rm, da

      integer n1, n2, nsv
      logical kdgtest
      logical isinside
      logical isinside2, scut, ecut
      real*8 x
      isinside(x) = mod(thetamin-thetamax+360.d0, 360.d0) .ge.
     *               mod(x-thetamax+360.d0, 360.d0)
      isinside2(x) = mod(ea-sa+360.d0, 360.d0) .ge.
     *               mod(x-sa+360.d0, 360.d0)

       r = Volat( comp%vol + ir)
       h = Volat( comp%vol + ih)
       sa= Volat( comp%vol + isa)
       ea= Volat( comp%vol + iea)
       scut = .false.
       ecut = .false.
       da = mod(ea-sa+360.d0, 360.d0)/2.d0
       rm = r*cos(da*Torad)
       da = (sa+da)*Torad
       mp%x = rm * cos(da)
       mp%y = rm * sin(da)
       

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
             call epdrawCcl(r,  0.d0, t1, t2, p(n+1), n1)
             nsv = n+1
             call epdrawCcl(r,  h, t1, t2, p(n+n1+1), n2)
             n =n + n1 + n2
             n = n + 1
             p(n)%x = gpsep
             if(kdgtest(howcyl, 1)) then
                call epdrawCcylEdg(p(nsv), n1, 0.d0,  mp, p(n+1), n2)
                n = n + n2 
             endif
             if(kdgtest(howcyl, 2)) then
                call epdrawCcylEdg( p(nsv+n1), n1, h,  mp, p(n+1), n2)
                n = n + n2
                n = n + 1
                p(n)%x = gpsep
             endif


             t1 =mod(thetamax, 360.d0)
             t2 = ea
             if(t2 .lt. t1) t2 = t2+360.d0
             call epdrawCcl(r,  0.d0, t1, t2, p(n+1), n1)
             nsv = n+1
             call epdrawCcl(r,  h, t1, t2, p(n+n1+1), n2)
             n =n+ n1 + n2
             n = n + 1
             p(n)%x = gpsep
             ecut = .true.
             if(kdgtest(howcyl, 1)) then
                call epdrawCcylEdg(p(nsv), n1, 0.d0, mp, p(n+1), n2)
                n = n + n2 
             endif
             if(kdgtest(howcyl, 2)) then
                call epdrawCcylEdg( p(nsv+n1), n1, h, mp, p(n+1), n2)
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
       call epdrawCcl(r,  0.d0, t1, t2, p(n+1), n1)
       nsv = n+1
       call epdrawCcl(r,  h, t1, t2, p(n+n1+1), n2)
       n = n + n1 + n2
       n = n + 1
       p(n)%x = gpsep
 100   continue
       if(kdgtest(howcyl, 1)) then
          call epdrawCcylEdg(p(nsv), n1, 0.d0, mp, p(n+1), n2)
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
         p(n)%x = r*cos(sa*Torad)
         p(n)%y = r*sin(sa*Torad)
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
         p(n)%x = r*cos(ea*Torad)
         p(n)%y = r*sin(ea*Torad)
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
      subroutine epdrawCcylEdg(peri, n1, h, centin,  p, n)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
      integer n1    ! input. number of points in peri
       type(epPos)::  peri(n1)  ! input. the vertex of cyl-like object
                               !  (circle or part of it, etc      )
      real*8 h    !i nuput. place of the surface height.
                  !   h=0 for floor, h> 0 for ceil.
       type(epPos)::  centin     ! center of the cut pos. (x,y) are used. 
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
               ! a surface which cuts the cyl-like object.  The 
               ! surface makes the right angle with the axis
               ! of the cyl-like object at height h. (local coordinate).
               ! The size will be 2*n1+2 
      integer n ! output. number of vertex put in p.

      integer i
       type(epPos)::  cent

      cent=centin
      cent%z = h

      do i = 1, n1-1
         p(i) = cent
      enddo
      n = n1 
      p(n)%x = gpsep

      do i = 1, n1
         p( n + i ) = peri(i)
      enddo
      n = n + n1
      n = n + 1
      p(n)%x = gpsep
      end
