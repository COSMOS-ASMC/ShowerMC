      subroutine epDraw_polyhed0(comp, p, n)
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
      call eppolyhed0Cnst(comp)
      call epDraw_polyhed00(comp, p, n)
      if( comp%struc(1:11) == "polyhed0_yz" .or. 
     *    comp%struc(1:11) == "polyhed0_zx"  ) then
         do i = 1, n
            if( p(i)%x /= gpsep )   then
               call epc2v_polyhed0(comp, p(i), pp)
               p(i) = pp
            endif
         enddo
      endif
      end

      subroutine epDraw_polyhed00(comp, p, np)
      use modpolyhed0
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"

       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
           !      volume in local coordnate.
           ! (x,y,z)= gpsep is a separator
          ! to be converted to a blank line
        ! dimension of p must be >+ (nvccl+2)*2
      integer  np    ! output.  number of (x,y,z) data
                     ! put in p.  

      integer i, j


      np =0 
      j = comp%vol+2
      do i = 1, npoly*3+1, 3
         np = np + 1
         if( i == npoly*3+1 ) then
            p(np)%x = Volat(j+1)
            p(np)%y = Volat(j+2)
            p(np)%z = Volat(j+3)
         else
            p(np)%x = Volat(j+i)
            p(np)%y = Volat(j+i+1)
            p(np)%z = Volat(j+i+2)
         endif
      enddo
      np = np + 1
      p(np)%x = gpsep

      j = comp%vol+2 + npoly*3
      do i = 1, npoly*3+1, 3
         np = np + 1
         if( i == npoly*3+1 ) then
            p(np)%x = Volat(j+1)
            p(np)%y = Volat(j+2)
            p(np)%z = Volat(j+3)
         else
            p(np)%x = Volat(j+i)
            p(np)%y = Volat(j+i+1)
            p(np)%z = Volat(j+i+2)
         endif
      enddo
      np = np + 1
      p(np)%x = gpsep
      np = np + 1
      p(np)%x = gpsep
      end
