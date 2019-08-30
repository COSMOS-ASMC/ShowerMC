      subroutine epLightGetPipeSurfN(comp, attr, aPos, surfn)
!         get surface number of a given point (x,y,z)
!         which  is assumed to be very close to one of the
!         3 surfaces of a cyl
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
       type(Component)::  comp   !input component
       type(epPos)::  aPos       ! input. track of which pos is examined
      real(8),intent(in):: attr(3)  ! r1, r2, h
      integer,intent(out)::surfn  ! output. surface number.  (1~4)
                     !   1--> bottom, 2 top, 3 inner side  4 outer side

       type(epPos)::  cpos

      call epv2c_pipe(comp, aPos, cpos)
      call epLightGetPipeSurfN0(attr(1), attr(2), attr(3), 
     *  cpos%x, cpos%y, cpos%z,  surfn)

      end subroutine epLightGetPipeSurfN


      subroutine epLightGetPipeSurfN0(r1, r2, h, x, y, z, surfn)
      use modepLightMaxDef
      implicit none
#include "ZepTrackv.h"

!         get surface number of a given point (x,y,z)
!         which  is assumed to be very close to one of the
!         3 surfaces of a cyl
!
      real(8),intent(in)::r1, r2, h  !inner and outer radius and hight (cm)
      real(8),intent(in):: x, y, z  !  given point
      integer,intent(out)::surfn  ! output. surface number.  (1~4)
                     !   1--> bottom, 2 top, 3 inner side  4 outer side

      real*8 o(4) 
      real*8  mind
      integer i

      if( z == 0.) then
         surfn = 1
      elseif( z == h ) then
         surfn = 2
      else
         o(1) = z**2
         o(2) = (h-z)**2
         o(3) = x**2 + y**2 - r1**2
         o(4) = x**2 + y**2 - r2**2
         mind = min( abs(o(1)), abs(o(2)), abs(o(3)), abs(o(4) ))
         do i = 1, 4
            if(mind .eq. abs(o(i)) ) exit
         enddo
         surfn = i
         if(abs(o(i)) .gt. EpsLength2**2 ) then
            write(0,*) ' warning: surface check  for pipe CrossMode=',
     *      CrossMode
            write(0,*) ' r1, r2, h =', r1,r2,  h
            write(0,*) ' x,y,z=', x, y, z
            write(0,*) ' boundary=',Move%boundary
            write(0,*) ' closest surf. # is =', i, mind
         endif   
      endif
      end   subroutine epLightGetPipeSurfN0
