      subroutine epLightGetPrismSurfN(comp, aPos, surfn)
      use prism
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
       type(Component)::  comp   !input component
       type(epPos)::  aPos    ! input. track of which pos is examined 
      integer,intent(out):: surfn  ! surface # obtained

       type(epPos)::  cpos
      call epv2c_prism(comp, aPos, cpos)  ! convert to canonical
      call  epLightGetPrismSurfN0(a, b, c, h,
     *       cpos%x, cpos%y, cpos%z, surfn)
      end       subroutine epLightGetPrismSurfN

      subroutine epLightGetPrismSurfN0(a, b, c, h, x, y, z, surfn)
      implicit none
!         get surface number of a given point (x,y,z)
!         which  is assumed to be very close to one of the
!         5 surfaces of a prism
!
      real(8),intent(in):: a, b, c, h  !input. prism parameters
!               see ../prog/NewVol/Fig/NewVol2.pdf
      real(8),intent(in):: x, y, z  !input    given point
                          ! in canoical space
      integer,intent(out):: surfn  ! output. surface number. 
!       surf #:   bottom 1:  \ 2:  / 5: x-z @y=0: 3; x-z@y=b; 4
!       
!               *         
!        5    *    *   2
!           *   3     *     (back 4)
!          --------------*
!               1
      real(8),parameter::eps=1.d-5   ! 1.d-6 is dangerous
      integer i

      if( abs(z) <= eps ) then
         surfn = 1
      elseif( abs(y) <= eps ) then
         surfn = 3
      elseif( abs(y-b) <=  eps ) then
         surfn = 4
      elseif(abs( (x-c)*h/c+h -z ) <= eps ) then
         surfn = 5
      else
         surfn = 2 
      endif
      end      subroutine epLightGetPrismSurfN0
