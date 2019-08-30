!     *****************************
      subroutine epDraw_chute(comp, p, n)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                              !   an chutee in local coord.
                              ! (x,y,z)= gpsep is a separator
                              ! to be converted to a blank line
                              ! dimension of p must be >+ (nvccl+2)*2
      integer  n              ! output.  number of (x,y,z) data
                              ! put in p.


      integer ia, ib, ih
      parameter( ia = 1,  ib = 2,  ih=3)

      real*8 a, b, h, dx
      integer n1, i
      logical kdgtest

      integer nvqellips
      parameter (nvqellips = nvccl/4+4)
       type(epPos)::  ptemp(nvqellips)

!

      a = Volat( comp%vol + ia)
      b = Volat( comp%vol + ib)
      h = Volat( comp%vol + ih)
      n =0

      if( kdgtest(how, 1) )then
         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = a
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = 0.
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = a
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif

      if( kdgtest(how, 3) )then
         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = 0.
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = h

         n = n + 1
         p(n)%x = 0.
         p(n)%y = b
         p(n)%z = h

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif
      if( kdgtest(how, 2) .or. kdgtest(how, 5)  .or.
     *  kdgtest(how, 6 ) ) then
         call epdrawElps(a, h, 0.d0, 180.d0, 270.d0, ptemp, n1)
         if(kdgtest(how, 2) ) then
!           since center is at (a,0,h), we have to 
!           exchange y<->z and shift
            do i = 1, n1-1
               p(n+i)%x = ptemp(i)%x + a
               p(n+i)%y = 0.d0
               p(n+i)%z = ptemp(i)%y + h
            enddo
            p(n+n1) = ptemp(n1)
            n = n + n1

            do i = 1, n1-1
               p(i+n)%x = ptemp(i)%x + a
               p(i+n)%y = 0.
               p(i+n)%z = 0.
            enddo
            n = n + n1
            p(n)%x = gpsep
            n = n + 1
            p(n)%x = gpsep
         endif
         if(kdgtest(how, 5) ) then
            do i = 1, n1-1
               p(n+i)%x = ptemp(i)%x + a
               p(n+i)%y = b
               p(n+i)%z = ptemp(i)%y + h
            enddo
            p(n+n1) = ptemp(n1)
            n = n + n1

            do i = 1, n1-1
               p(i+n)%x = ptemp(i)%x + a
               p(i+n)%y = b
               p(i+n)%z = 0.
            enddo
            n = n + n1
            p(n)%x = gpsep
            n = n + 1
            p(n)%x = gpsep
         endif

         if( kdgtest(how, 6) ) then
            do i = 1, n1-1
               p(i+n)%x = ptemp(i)%x + a
               p(i+n)%y = 0.
               p(i+n)%z = ptemp(i)%y + h
            enddo
            n = n + n1
            p(n)%x = gpsep
            do i = 1, n1-1
               p(i+n)%x = ptemp(i)%x + a
               p(i+n)%y = b
               p(i+n)%z = ptemp(i)%y + h
            enddo
            n = n + n1
            p(n)%x = gpsep
            n = n + 1
            p(n)%x = gpsep
         endif
      endif
      end

