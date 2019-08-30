      subroutine epLightGetCylSurfN(comp, aPos, r, h, surfn)
      use modepLightMaxDef
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"

!         get surface number of a given point (x,y,z)
!         which  is assumed to be very close to one of the
!         3 surfaces of a cyl
!
       type(Component)::  comp  ! input componennt
      real(8),intent(in)::r, h  ! radius and hight (cm)
       type(epPos)::   aPos  ! input. given point
      integer,intent(out)::surfn  ! output. surface number.  (1~3)
                     !   1--> bottom, 2 top, 3 side
       type(epPos)::  cpos

      call epv2c_cyl(comp, aPos, cpos)  ! convert to canonical                                   
      call  epLightGetCylSurfN0( r, h, 
     *       cpos%x, cpos%y, cpos%z, surfn)

      end     subroutine  epLightGetCylSurfN


      subroutine epLightGetCylSurfN0(r, h,  x, y, z, surfn)
      use modepLightMaxDef
      implicit none
#include "ZepTrackv.h"

!         get surface number of a given point (x,y,z)
!         which  is assumed to be very close to one of the
!         3 surfaces of a cyl
!
      real(8),intent(in)::r, h  ! radius and hight (cm)
      real(8),intent(in):: x, y, z  !  given point
      integer,intent(out)::surfn  ! output. surface number.  (1~3)
                     !   1--> bottom, 2 top, 3 side

      real*8 o(3) 
      real*8  mind
      integer i

      if( z == 0.) then
         surfn = 1
      elseif( z == h ) then
         surfn = 2
      else
         o(1) = z**2
         o(2) = (h-z)**2
         o(3) = x**2 + y**2 - r**2
         mind = min( abs(o(1)), abs(o(2)), abs(o(3)))
         do i = 1, 3
            if(mind .eq. abs(o(i)) ) exit
         enddo
         surfn = i
         if(abs(o(i)) .gt. EpsLength2**2 ) then
            write(0,*) ' warning: surface check  for cyl; CrossMode=',
     *      CrossMode
            write(0,*) ' r, h =', r, h
            write(0,*) ' x,y,z=', x, y, z
            write(0,*) ' boundary=',Move%boundary
            write(0,*) ' closest surf. # is =', i, mind
         endif   
      endif
      end subroutine  epLightGetCylSurfN0
