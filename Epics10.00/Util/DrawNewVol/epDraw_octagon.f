      subroutine epDraw_octagon(comp, p, n)
      use octagon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)      ! output. (x,y,z) to describe
                               ! octagon in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  
      integer i
      real(8):: temp, temp2

      call epDraw_octagn0(comp, p, n)
      if( comp%struc == 'octagon' .or. 
     *    comp%struc == 'octagon_w' ) then
      else
!      elseif(comp.struc == 'octagon_x') then
         do i = 1, n
            if( p(i)%x /= gpsep ) then
               call epc2v_octagon(comp, p(i), p(i))
!               temp = b-p(i).y 
!               p(i) = epPos(temp, p(i).x, p(i).z)
            endif
         enddo
!      elseif(comp.struc == 'octagon_z') then
!         do i = 1, n
!            if( p(i).x /= gpsep ) then
!               temp = p(i).x
!               temp2 = b- p(i).y
!               p(i) = epPos(p(i).z, temp2, temp) 
!            endif
!         enddo
!      else
!         write(0,*) ' comp.struc=',comp.struc, ' error'
!         stop
      endif
      end

      subroutine epDraw_octagn0(comp, p, n)
      use octagon
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)      ! output. (x,y,z) to describe
                               ! octagon in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  

      logical kdgtest
 

      integer base, i, i1, i2, j1, j2

      call epoctagonCnst(comp)

      n = 0

 !             x =0
      n = n + 1
!      p(n).x = 0.           
!      p(n).y = d     
!      p(n).z = 0.           
!       assignment like  above sometimes fails (say p(n).y = 0
!          may result while d !=0; compiler bug ???)
!         so we use following notation
      p(n) = epPos(0.d0, d, 0.d0)

      n = n + 1
!      p(n).x =  0.
!      p(n).y = b-d       !---
!      p(n).z =  0.
      p(n) = epPos(0.d0, b-d, 0.d0)

      n = n + 1
!      p(n).x= 0.
!      p(n).y= b
!      p(n).z= d         ! /
      p(n) = epPos(0.d0, b, d)

      n = n + 1
!      p(n).x= 0.         ! |
!      p(n).y= b
!      p(n).z= c - d  
      p(n) = epPos(0.d0, b, c-d)

      n = n + 1
      p(n)= epPos( 0.d0, b-d, c)
!      p(n).x= 0.        ! \
!/////////
!      write(0,*) ' p(n).x =',p(n).x, b-d
!//////////
!      p(n).z= c  
!      p(n).y= b-d
!/////////
!      write(0,*) ' p(n).y =',p(n).y, b-d
!//////////
      n = n + 1
!      p(n).x= 0.         ! --
!      p(n).y= d
!      p(n).z= c  
      p(n) = epPos(0.d0, d, c)

      n = n + 1
!      p(n).x= 0.         ! /
!      p(n).y= 0.
!      p(n).z= c - d 
      p(n) = epPos(0.d0, 0.d0, c-d)

      n = n + 1
!      p(n).x= 0.         ! |
!      p(n).y= 0.
!      p(n).z= d 
      p(n) = epPos(0.d0, 0.d0, d)

      n = n + 1
!      p(n).x= 0.        ! \
!      p(n).y= d
!      p(n).z= 0. 
      p(n) = epPos(0.d0, d, 0.d0)
!////////
!         write(0,*)  p(n).x, p(n).y, p(n).z, d
!//////////
      n = n + 1
!      p(n).x = gpsep
!      p(n).y = 0.   ! dummy
!      p(n).z = 0.   !//
      p(n) = epPos( gpsep, 0.d0, 0.d0)

!        -----------  x = a
      i2 = n
      do i = 1, i2
         n = n + 1
         p(n) = p(i)
         if(p(n)%x /= gpsep ) then
            p(n)%x = a
         endif
      enddo
      n = n + 1
!      p(n).x = gpsep
!      p(n).y =  0.  !dummy
!      p(n).z =  0.  ! //
      p(n) = epPos( gpsep, 0.d0, 0.d0)
!        =============

!         side 
       n = n + 1
       i1 = n
!////////  next 2 lines behave strange; should be compiler problem
!        (  --> d is not set )
!          if we insert some write statement, it disappears
!
 !!      p(n) = o      !\
 !!      p(n).z = d    
!          next two are also NG
 !!       p(n) = o      !\
 !!      p(n)%z = d    

!          so above 2 are replaced by next
!       p(n).x = 0.
!       p(n).y = 0.
!       p(n).z = d
       p(n) = epPos(0.d0, 0.d0, d)

       n = n + 1   
!       p(n).x = 0.
!       p(n).y = d
!       p(n).z = 0.
       p(n) = epPos(0.d0, d, 0.d0)

       n = n + 1    ! --
!       p(n).x = 0.
!       p(n).y = b - d
!       p(n).z = 0. 
       p(n) = epPos(0.d0, b-d, 0.d0)

       n = n + 1    
!       p(n).x = 0.   ! /
!       p(n).y = b
!       p(n).z = d
       p(n) = epPos(0.d0, b, d)

       n = n + 1
!       p(n).x = gpsep
       p(n) = epPos(gpsep, 0.d0, 0.d0)

       i2 = n
!/////////////
!       write(0,*) ' lower side n=', i2-i1
!       do i = i1,i2
!          write(0,*) p(i).x, p(i).y, p(i).z
!       enddo
!//////////////
!        -----------   upper half
!       j1 = n +1
       do i = i1, i2
          n = n + 1
          p(n) = p(i)
          p(n)%z = (c -p(i)%z) 
       enddo
!       j2 = n
!/////////////
!       write(0,*) ' upper side n=', j2-j1 
!       do i = j1,j2
!          write(0,*) p(i).x, p(i).y, p(i).z
!       enddo
!//////////////
       n = n +1
!       p(n).x = gpsep
       p(n) = epPos(gpsep, 0.d0, 0.d0)
 !        =============   x=a
       i2 = n
       do  i = i1, i2
          n = n + 1
          p(n) = p(i)
          if(  p(i)%x /= gpsep ) then
             p(n)%x = a
          endif
       enddo
      end
