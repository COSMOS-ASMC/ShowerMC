      subroutine epDraw_fpolygon(comp, p, n)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"


       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                               !  (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  

      integer i
       type(epPos)::  pp
      call epfpolygonCnst(comp)
      call epDraw_fpolygon0(comp, p, n)
      if( comp%struc(1:11) == "fpolygon_yz" .or. 
     *    comp%struc(1:11) == "fpolygon_zx"  ) then
         do i = 1, n
            if( p(i)%x /= gpsep )   then
               call epc2v_fpolygon(comp, p(i), pp)
               p(i) = pp
            endif
         enddo
      endif
      end

      subroutine epDraw_fpolygon0(comp, p, np)
      use fpolygon
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

      

 
      integer i, j
      np =0 
      j = comp%vol+2
      do i = 1, npoly*2+1, 2
         np = np + 1
         if( i == npoly*2+1 ) then
            p(np)%x = Volat(j+1)
            p(np)%y = Volat(j+2)
            p(np)%z = 0.
         else
            p(np)%x = Volat(j+i)
            p(np)%y = Volat(j+i+1)
            p(np)%z = 0.
         endif
      enddo
      np = np + 1
      p(np)%x = gpsep

      do i = 1, npoly*2+1, 2
         np = np + 1
         if( i == npoly*2+1 ) then
            p(np)%x = Volat(j+1)
            p(np)%y = Volat(j+2)
            p(np)%z = height
         else
            p(np)%x = Volat(j+i)
            p(np)%y = Volat(j+i+1)
            p(np)%z = height
         endif
      enddo
      np = np + 1
      p(np)%x = gpsep
      np = np + 1
      p(np)%x = gpsep

      end
