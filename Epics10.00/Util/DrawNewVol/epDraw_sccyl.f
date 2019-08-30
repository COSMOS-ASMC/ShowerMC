!     *****************************
      subroutine epDraw_sccyl(comp, p, n)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
#include "ZepDirec.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                              !   a sccyl in local coordnate.
                              ! (x,y,z)= gpsep is a separator
                              ! to be converted to a blank line
                              ! dimension of p must be >+ (nvccl+2)*2
      integer  n              ! output.  number of (x,y,z) data
                              ! put in p.  

      integer ir, ih, in1x, in1y, in1z, in2x, in2y, in2z, isa, iea
      integer maxz, minz, isapx, isapy, isapz1, isapz2,
     *        ieapx, ieapy, ieapz1, ieapz2, imaxz, iminz

      parameter (ir = 1,  ih = 2,  in1x= 3, in1y=4, in1z=5)
      parameter (in2x= 6, in2y=7, in2z=8, maxz=9, minz=10  )
      parameter (isa = 11, iea=12, isapx=13, isapy=14, isapz1=15,
     *      isapz2=16, ieapx=17, ieapy=18, ieapz1=19, ieapz2=20,
     *      imaxz=21, iminz=22)

      real*8 r, h, sa, ea
       type(epDirec)::  n1, n2





      real*8 t1, t2
       type(epPos)::  mp
      real*8 rm, da, k2, cost, sint, h1, h2

      integer j1, j2, nsv1, nsv2, i
      logical kdgtest
      logical isinside
      logical isinside2, scut, ecut
      real*8 x
      isinside(x) = mod(thetamin-thetamax+360.d0, 360.d0) .ge.
     *               mod(x-thetamax+360.d0, 360.d0)
      isinside2(x) = mod(ea-sa+360.d0, 360.d0) .ge.
     *               mod(x-sa+360.d0, 360.d0)


!     
      r = Volat( comp%vol + ir)
      h = Volat( comp%vol + ih)
      n1%x = Volat( comp%vol + in1x)
      n1%y = Volat( comp%vol + in1y)
      n1%z = Volat( comp%vol + in1z)

      n2%x = Volat( comp%vol + in2x)
      n2%y = Volat( comp%vol + in2y)
      n2%z = Volat( comp%vol + in2z)

      k2 = h* n2%z

      sa =  Volat( comp%vol + isa )
      ea =  Volat( comp%vol + iea ) 

      scut = .false.
      ecut = .false.
      da = mod(ea-sa+360.d0, 360.d0)/2.d0
      rm = r*cos(da*Torad)
      da = (sa+da)*Torad
      mp%x = rm * cos(da)
      mp%y = rm * sin(da)
      cost = mp%x/r
      sint = mp%y/r
      h1 = -r*(n1%x*cost + n1%y*sint )/n1%z
      h2 = (k2-r*(n2%x*cost + n2%y*sint ))/n2%z

      n = 0
      if(.not. isinside(sa) .and. isinside(ea) ) then
          t1 = thetamax
          t2 = ea
          ecut = .true.
       elseif(isinside(sa) .and. isinside(ea)) then
          if(isinside2(thetamin)) then
!              there are two regions
             scut = .true.
             t1 =  sa
             t2 = mod(thetamin,360.d0)
             if(t2 .lt. t1) t2 = t2+360.d0
             call epdrawCcl(r, 0.d0, t1, t2, p(n+1), j1)
             nsv1 = n+1
!               adjust z
             do i = nsv1, nsv1+ j1 - 1
                cost = p(i)%x/r
                sint = p(i)%y/r
                p(i)%z = -r*( n1%x*cost + n1%y*sint )/n1%z
             enddo

             n = n + j1

             call epdrawCcl(r,  h, t1, t2, p(n+1), j2)
             nsv2 = n + 1
!               adjust z
             do i = nsv2, nsv2+ j2 - 1
                cost = p(i)%x/r
                sint = p(i)%y/r
                p(i)%z = (k2-r*(n2%x*cost + n2%y*sint ))/n2%z
             enddo
             n = n +  j2
             n = n + 1
             p(n)%x = gpsep
             if(kdgtest(howcyl, 1)) then
                call epdrawCcylEdg(p(nsv1), j1, h1,  mp, p(n+1), j2)
                n = n + j2 
             endif
             if(kdgtest(howcyl, 2)) then
                call epdrawCcylEdg( p(nsv2), j1, h2,  mp, p(n+1), j2)
                n = n + j2
                n = n + 1
                p(n)%x = gpsep
             endif


             t1 =mod(thetamax, 360.d0)
             t2 = ea
             if(t2 .lt. t1) t2 = t2+360.d0
             call epdrawCcl(r,  0.d0, t1, t2, p(n+1), j1)
             nsv1 = n+1
!               adjust z
             do i = nsv1, nsv1+ j1 - 1
                cost = p(i)%x/r
                sint = p(i)%y/r
                p(i)%z = -r*( n1%x*cost + n1%y*sint )/n1%z
             enddo
             n = n+ j1 

             call epdrawCcl(r,  h, t1, t2, p(n+1), j2)
             nsv2 = n+1
!               adjust z
             do i = nsv2, nsv2+ j2 - 1
                cost = p(i)%x/r
                sint = p(i)%y/r
                p(i)%z = (k2-r*(n2%x*cost + n2%y*sint ))/n2%z
             enddo

             n =n+ j2
             n = n + 1
             p(n)%x = gpsep
             ecut = .true.
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
       call epdrawCcl(r,  0.d0, t1, t2, p(n+1), j1)
       nsv1 = n+1
!               adjust z
       do i = nsv1, nsv1+ j1 - 1
          cost = p(i)%x/r
          sint = p(i)%y/r
          p(i)%z = -r*( n1%x*cost + n1%y*sint )/n1%z
       enddo

       n = n + j1
       call epdrawCcl(r,  h, t1, t2, p(n+1), j2)
       nsv2 = n+1
!               adjust z
       do i = nsv2, nsv2+ j2 - 1
          cost = p(i)%x/r
          sint = p(i)%y/r
          p(i)%z = (k2-r*(n2%x*cost + n2%y*sint ))/n2%z
       enddo
       n = n + j2
       n = n + 1
       p(n)%x = gpsep
 100   continue
       if(kdgtest(howcyl, 1)) then
          call epdrawCcylEdg(p(nsv1), j1, h1, mp, p(n+1), j2)
          n = n + j2 
       endif
       if(kdgtest(howcyl, 2)) then
          call epdrawCcylEdg( p(nsv2), j1, h2,  mp, p(n+1), j2)
          n = n + j2
          n = n + 1
          p(n)%x = gpsep
       endif


       if(scut) then
!           draw starting cut  square
         n = n + 1
         p(n)%x = mp%x
         p(n)%y = mp%y
         p(n)%z =  h1

         n = n + 1
         p(n)%x = r*cos(sa*Torad)
         p(n)%y = r*sin(sa*Torad)
         cost = p(n)%x/r
         sint = p(n)%y/r
         p(n)%z =-r*(n1%x*cost + n1%y*sint )/n1%z

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = mp%x 
         p(n)%y = mp%y
         p(n)%z = h2

         n = n + 1
         p(n)%x = p(n-3)%x
         p(n)%y = p(n-3)%y
         cost = p(n)%x/r
         sint = p(n)%y/r
         p(n)%z = (k2-r*(n2%x*cost + n2%y*sint ))/n2%z


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
         cost = p(n)%x/r
         sint = p(n)%y/r
         p(n)%z =-r*(n1%x*cost + n1%y*sint )/n1%z


         n = n + 1
         p(n)%x = mp%x
         p(n)%y = mp%y
         cost = p(n)%x/r
         sint = p(n)%y/r
         p(n)%z = h1

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = p(n-3)%x
         p(n)%y = p(n-3)%y
         cost = p(n)%x/r
         sint = p(n)%y/r
         p(n)%z = (k2-r*(n2%x*cost + n2%y*sint ))/n2%z

         n = n + 1
         p(n)%x = mp%x
         p(n)%y = mp%y
         cost = p(n)%x/r
         sint = p(n)%y/r
         p(n)%z = h2

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif
      end

