!       get a xrossin point of a half line with a
!       curved part of a cylinder.
!
      subroutine epxpcylwall(r, h, sa, ea, pos, dir,
     *      length, icon)
      implicit none
#include "Zglobalc.h"
#include "ZepPos.h"
#include "ZepDirec.h"

      real*8 r  !  input. radius of the cylinder
      real*8 h  !  input. heigth of //
      real*8 sa  ! input. starting angle of the cylinder (deg)
      real*8 ea  ! input. ending angle //
       type(epPos)::  pos  ! input. starting pos. of the line
       type(epDirec)::  dir  ! input. lines direction cos.
      real*8 length ! output. length from pos to  x.p 
                    !                                 
      integer icon  !  output. 0 x.p is obtained. length  can be used
                    !            to get the x.p as pos+length*dir
                    !            pos is inside r **( i.e., z is
                    !            and  angle region is not judged.
                    !          1 x.p is obtained. length can be used
                    !            to get the x.p as pos+length*dir
                    !            pos is outside r
                    !         -1 no x.p
!
!                x.p is obtained is on the curved part only.
!
       type(epPos)::  xp
      real*8  aa, bb, cc, dd, ang 

      logical isinside
      real*8 x, leng2
      isinside(x) = mod(ea-sa+360.d0, 360.d0) .ge.
     *               mod(x-sa+360.d0, 360.d0)

      
      aa = dir%x**2 + dir%y**2
      if(aa .ne. 0.) then
         bb = (pos%x*dir%x + pos%y*dir%y)
         cc = (pos%x**2+pos%y**2) - r**2
         dd =  bb**2 - aa*cc
         if(dd .ge. 0.d0) then
            dd = sqrt(dd)

            if(cc .ge. 0.) then
!                 pos is outside of the circle
               leng2 = (-bb + dd)/aa
               if(leng2 .lt. 0.) then
                  icon = -1
                  return ! *********
               endif
               length = (-bb - dd)/aa

               xp%x = pos%x + length*dir%x
               xp%y = pos%y + length*dir%y
               xp%z = pos%z + length*dir%z
               ang= atan2(xp%y, xp%x)*Todeg
               if(xp%z .lt. 0. .or. xp%z .gt. h .or.
     *          .not. isinside(ang)) then
                  length = leng2
                  xp%x = pos%x + length*dir%x
                  xp%y = pos%y + length*dir%y
                  xp%z = pos%z + length*dir%z
                  if(xp%z .lt. 0. .or. xp%z .gt. h) then
                     icon = -1
                  else
                     ang= atan2(xp%y, xp%x)*Todeg
                     if(.not. isinside(ang) )then
                        icon = -1
                     else
                        icon = 1
                     endif
                  endif
               else
                  icon = 1
               endif
            else
!                 pos. is inside of the circle
               length = (-bb + dd)/aa
               xp%x = pos%x + dir%x*length
               xp%y = pos%y + dir%y*length
               xp%z = pos%z + dir%z*length
               if(xp%z .lt. 0. or. xp%z .gt. h) then
                  icon = -1
               else
                  ang = atan2(xp%y, xp%x)*Todeg
                  if( .not. isinside(ang))  then
                     icon = -1
                  else
                     icon = 0
                  endif
               endif
            endif
         else
            icon = -1
         endif
      else
         icon = -1
      endif
      end
