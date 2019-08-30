      subroutine epLightGetOctagonSurfN(comp, aPos, surfn)
      use octagon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
       type(Component)::  comp   !input component
       type(epPos)::  aPos    ! input. track of which pos is examined 
      integer,intent(out):: surfn  ! surface # obtained

       type(epPos)::  cpos
!      call epoctagonCnst(comp)
      call epv2c_octagon(comp, aPos, cpos)  ! convert to canonical
      call  epLightGetOctaSurfN0(a, b, c, d,
     *       cpos%x, cpos%y,c pos%z, surfn)
      end       subroutine epLightGetOctagonSurfN

      subroutine epLightGetOctaSurfN0(a, b, c, d, x, y, z, surfn)
      use modepLightMaxDef
      implicit none
!         get surface number of a given point (x,y,z)
!         which  is assumed to be very close to one of the
!         10 surfaces of a octagon
!
      real(8),intent(in):: a, b, c, d  !input. octagon parameters
!               see ../prog/NewVol/Fig/NewVol2.pdf
      real(8),intent(in):: x, y, z  !input    given point
      integer,intent(out):: surfn  ! output. surface number. 
!           z  6 
!           |_____  
!         7 /      \10
!           |      |   5 is oposit of 2
!         3 |  2   | 4
!           |      |
!         8  \_____/9 _____ y
!               1
!         if 7~10 shrink, the numbering is the same as for box.
!
      real*8 o(10) 
      real*8  mind
      real(8),parameter::eps=0.d0
      integer i

      if( abs(z) <= eps ) then
         surfn = 1
      elseif( abs(z-c) <= eps ) then
         surfn = 6
      elseif( abs(x) <=  eps ) then
         surfn = 5
      elseif( abs(x-a) <= eps ) then
         surfn = 2
      elseif( abs(y)  <= eps ) then
         surfn = 3
      elseif( abs( y-b) <= eps ) then
         surfn = 4
      else
  !             surface # vs pos.
         o(1) = z  
         o(2) = a-x
         o(3) = y
         o(4) = b-y
         o(5) = x
         o(6) = c-z
         o(7) = z - y - c + d  !     z =  y+c-d
         o(8) = z + y - d        !   z= -y + d
         o(9) = z - y +(b-d)     !   z= y - (b-d)
         o(10) = z+y-(b-d)-c     !   z= -y + (b-d)+c
         mind = 1.d20
         do i = 1, 10
            if( mind > abs(o(i))) then
               mind = abs( o(i) )
               surfn = i
            endif
         enddo
         if(abs(o(surfn)) .gt. EpsLength2 ) then
            write(0,*) ' warning: surface check #=',surfn, o(surfn)
            write(0,*) ' a, b,c, d=', a,b,c,d
            write(0,*) ' x, y,z=', x,y,z
            write(0,*) ' for octagnon'
            write(0,*) ' '
         endif   
      endif
      end      subroutine epLightGetOctaSurfN0
