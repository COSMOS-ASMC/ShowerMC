      subroutine epDraw_torus(comp, p, n)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"


       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                               !   a elightg in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line
                               ! dimension of p must be >+ (nvccl+2)*2
      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  

      integer i
       type(epPos)::  pp

      call epDraw_torus0(comp, p, n)
      if( comp%struc(1:7) == "torus_y" .or. 
     *    comp%struc(1:7) == "torus_x"  ) then
         do i = 1, n
            if( p(i)%x /= gpsep )   then
               call epc2v_torus(comp, p(i), pp)
               p(i) = pp
            endif
         enddo
      endif
      end

      subroutine epDraw_torus0(comp, p, np)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"

       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                               !   a elightg in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line
                               ! dimension of p must be >+ (nvccl+2)*2
      integer  np               ! output.  number of (x,y,z) data
                               ! put in p.  

      

      real(8):: R, sr
      integer,parameter:: n=nvccl  ! circle is approximated by n edges
      real(8),parameter :: da=2*pi/n
 
      real(8):: cosd, sind, temp

      integer i,  nn, j, j1, j2

      R = Volat( comp%vol +  1) 
      sr =Volat( comp%vol +  2) 

         ! make the circle perpendicular to the (x-y) plane
         !  next is ccl on (xy); center is (0,0)
         !                         below is not opposit !
      call epDrawCcl(sr, 0.d0,   thetamax, thetamin,
     *      p, nn)
!           modify origin and plane
         !  x-->x y-->z z-->y  =0
      do j = 1, nn
         p(j)%z = p(j)%y 
         if( p(j)%x /= gpsep ) then
            p(j)%x = p(j)%x + R
         endif
         p(j)%y = 0.
      enddo         
      np = nn

      cosd = cos(da)
      sind = sin(da)
      j2  = 0 
      do i = 2, n + 1
            !     rotate by exp(i*da) (only,xy)
            !    (x+iy)*(cos(da) + isin(da))
            !    (xcos -ysin, ycos +xsin)
         j1 = j2 + 1
         j2 = j1 + nn -1
         do j = j1, j2
            if(p(j)%x /= gpsep) then
               p(j+nn)%x = p(j)%x*cosd - p(j)%y*sind
               p(j+nn)%y = p(j)%y*cosd + p(j)%x*sind
               p(j+nn)%z = p(j)%z
            else
               p(j+nn)%x = p(j)%x
            endif   
         enddo
         np = np + nn
      enddo
      np = np + 1
      p(np)%x = gpsep
      np = np + 1
      p(np)%x = gpsep
      end
