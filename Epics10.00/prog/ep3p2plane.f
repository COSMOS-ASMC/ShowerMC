      subroutine epxpLand4vp(p1, p2, p3, p4, a, b, l, icon, jcon)
      implicit none
#include "ZepPos.h"
!        get a crossing point of a half line with a plane formed
!        by the 4 points on a plane. and judge if the x point 
!        is inside of the "square" formed by p1,p2,p3,p4.
       type(epPos)::  p1,p2,p3,p4 ! input.  
       type(epPos)::  a  ! input.  starting point of the half line
       type(epPos)::  b  ! input.  direction cos of the line |b| =1
      real*8  l         ! output. a+lb is the crossing point.
                        !        l may be negative

      integer icon  ! output. 0: p is inside of the "square"
                    !         1: p is on  p1-p2
                    !         2: p is on p2-p3
                    !         3: p is on p3-p4
                    !         4: p is on p4-p1 
                    !         5: p is outside of the "square"
                    !         6: line is parallel to the plane
      integer jcon  ! output. has meaning if icon = 1~4
                    !         0: p is not on the vertex
                    !         i: p=pi (i=1~4)

       type(epPos)::  n, p
      real*8 k
      call ep3p2plane(p1, p2, p3, n, k)
      call epxpLandP(a, b, n, k, l, icon)
      if(icon .le. 2) then
         p%x = a%x + b%x*l
         p%y = a%y + b%y*l
         p%z = a%z + b%z*l
         call episPin4(p1, p2, p3, p4, n, p, icon, jcon)
      else
         icon = 6
      endif
      end
      subroutine epxpLand4vp2(p1, p2, p3, p4, a, b, n, k, l, icon, jcon)
      implicit none
#include "ZepPos.h"
!        get a crossing point of a half line with a plane 
!        expressed by normal verctor and distance form origin to
!        the plane. 
       type(epPos)::  p1, p2, p3, p4 ! input. 4 points forming a plane
       type(epPos)::  a  ! input.  starting point of the half line
       type(epPos)::  b  ! input.  direction cos of the line |b| =1
       type(epPos)::  n  ! input.  normal vector
      real*8 k          ! input. scaler product of n and any
                        ! point on the plane.  (normally distance
                        ! from the origin to the plane
      real*8  l         ! output. a+lb is the crossing point.
                        !   l may be negative.

      integer icon  ! output. 0: p is inside of the "square"
                    !         1: p is on  p1-p2
                    !         2: p is on p2-p3
                    !         3: p is on p3-p4
                    !         4: p is on p4-p1 
                    !         5: p is outside of the "square"
                    !         6: line is parallel to the plane
      integer jcon  ! output. has meaning if icon = 1~4
                    !         0: p is not on the vertex
                    !         i: p=pi (i=1~4)

       type(epPos)::  p

      call epxpLandP(a, b, n, k, l, icon)
      if(icon .le. 2) then
         p%x = a%x + b%x*l
         p%y = a%y + b%y*l
         p%z = a%z + b%z*l
         call episPin4(p1, p2, p3, p4, n, p, icon, jcon)
      else
         icon = 6
      endif
      end
      subroutine epxpLand3vp(p1, p2, p3, a, b, l, icon, jcon)
      implicit none
#include "ZepPos.h"
!        get a crossing point of a half line with a plane formed
!        by the 3 points. and judge if the x point 
!        is inside of the "triangle" formed by p1,p2,p3
       type(epPos)::  p1,p2,p3 ! input.  
       type(epPos)::  a  ! input.  starting point of the half line
       type(epPos)::  b  ! input.  direction cos of the line |b| =1
      real*8  l         ! output. a+lb is the crossing point.  l  >=0
                        !        

      integer icon  ! output. 0: p is inside of the "triangle"
                    !         1: p is on  p1-p2
                    !         2: p is on p2-p3
                    !         3: p is on p3-p1 
                    !         4: p is outside of the "trianble"
                    !         5: line is parallel to the plane
      integer jcon  ! output. has meaning if icon = 1~3
                    !         0: p is not on the vertex
                    !         i: p=pi (i=1~3)

       type(epPos)::  n, p
      real*8 k

      call ep3p2plane(p1, p2, p3, n, k)
      call epxpLandP(a, b, n, k, l, icon)
      if(icon .le. 2) then
         p%x = a%x + b%x*l
         p%y = a%y + b%y*l
         p%z = a%z + b%z*l
         call episPin3(p1, p2, p3, n, p, icon, jcon)
!
!      integer icon  ! output. 0: p is inside of the triangle formed by p1,2,3.
!                    !         1: p is on the line, p1-p2.
!                    !         2: p is on the p2-p3  line
!                    !         3: p is on the p3-p1 line. 
!                    !         4: p is outside of the triangle.
!     integer jcon  ! output. has meaning only if icon =1~3.
!                    !         0--> the point is not on the vertex
!                    !         1--> p = p1
!                    !         2--> p = p2
!                    !         3--> p = p3 
!
      else
         icon = 5
      endif
      end
      subroutine epxpLand3vpB(p1, p2, p3, a, b, l, towhich, icon)
!/////////////
!      use sqTccl,only:check 
!///////////
      implicit none
#include "ZepPos.h"
!        get a crossing point of a half line with a plane formed
!        by the 3 points. and judge if the x point 
!        is inside of the "triangle" formed by p1,p2,p3
       type(epPos)::  p1,p2,p3 ! input.  
       type(epPos)::  a  ! input.  starting point of the half line
       type(epPos)::  b  ! input.  direction cos of the line |b| =1
      real(8),intent(out):: l         ! output. a+lb is the crossing point.  l  >=0
                        !        
      real(8),intent(out):: towhich  ! scaler proc. of plane's  normal vecter  and the b
                                   ! obtained if icon <=3. 
                                   ! can be used  the line direction is going inward
                                   ! or outward of an object.
      integer icon  ! output. 0: p is inside of the "triangle"
                    !         1: p is on  p1-p2
                    !         2: p is on p2-p3
                    !         3: p is on p3-p1 
                    !         4: p is outside of the "trianble"
                    !         5: line is parallel to the plane

       type(epPos)::  n, p
      real*8 k

      integer::jcon 

      call ep3p2plane(p1, p2, p3, n, k)
      call epxpLandP(a, b, n, k, l, icon)
!///////////////
!      if(check == 1 ) then
!         write(0,*) ' n, k=', n, k
!         write(0,*) ' landP l,icon =',l, icon
!      endif
!/////////////
      if(icon <= 2) then
         p%x = a%x + b%x*l
         p%y = a%y + b%y*l
         p%z = a%z + b%z*l
         call episPin3(p1, p2, p3, n, p, icon, jcon)
!///////////
!         if(check == 1 ) then
!            write(0,*) 'xp =',p
!            write(0,*) 'Pin3 icon, jcon =', icon, jcon
!         endif
!////////////
!
!      integer icon  ! output. 0: p is inside of the triangle formed by p1,2,3.
!                    !         1: p is on the line, p1-p2.
!                    !         2: p is on the p2-p3  line
!                    !         3: p is on the p3-p1 line. 
!                    !         4: p is outside of the triangle.
!     integer jcon  ! output. has meaning only if icon =1~3.
!                    !         0--> the point is not on the vertex
!                    !         1--> p = p1
!                    !         2--> p = p2
!                    !         3--> p = p3 
         if( icon <= 3 ) then
            call epscalerProd(n, b, towhich)
         else
            towhich = 0.
         endif
      else
         icon = 5
      endif
      end
      
      subroutine epxpLandP(a, b, n, k, l, icon)
      implicit none
#include "ZepPos.h"
!        get a crossing point of a half line with a given plane.
       type(epPos)::  a  ! input.  starting point of the half line
       type(epPos)::  b  ! input.  direction cos of the line |b| =1
       type(epPos)::  n  ! input   normal vector of the plane 
      real*8  k         ! input.  n.x=k where x is any point on the plane
      real*8  l         ! output. a+lb is the crossing point.  
      integer icon      ! output. 0: crossing point obtained. (l>=0)
                        !         1: crossing point is at the backside (l<0).
                        !         2: line seems to be on the plane (l=0)
                        !         3: line seems to be parallel to the plane


      real*8 na, nb
      real*8 temp, eps
      data eps/1.d-8/

      call epscalerProd(n, a, na)
      call epscalerProd(n, b, nb)
      temp = k - na
      if(nb /= 0.) then
         l = temp/nb
         if(l .lt. 0.) then
            icon = 1
         else
            icon = 0
         endif
      else
         if( abs(temp) < eps) then
            icon = 2
            l  = 0.
         else
            icon = 3
         endif
      endif
      end

      subroutine episPin3(p1, p2, p3, n, p, icon, jcon)
      implicit none
#include "ZepPos.h"
       type(epPos)::  p1 ! input.  p1,p2,p3 are  non 
       type(epPos)::  p2 ! input. colinear points
       type(epPos)::  p3 ! input. to form a plane
       type(epPos)::  n  ! input. normal vector of the plane
       type(epPos)::  p  ! input given point which is on the plane fomred by
                        !  three points p1, p2, p3

      integer icon  ! output. 0: p is inside of the triangle formed by p1,2,3.
                    !         1: p is on the line, p1-p2.
                    !         2: p is on the p2-p3  line
                    !         3: p is on the p3-p1 line. 
                    !         4: p is outside of the triangle.
      integer jcon  ! output. has meaning only if icon =1~3.
                    !         0--> the point is not on the vertex
                    !         1--> p = p1
                    !         2--> p = p2
                    !         3--> p = p3 
      real*8 x1, y1, x2, y2, x3, y3, x, y
      real*8 nmax, temp, a1, b1, a2, b2, u, v, p12, p23, p31
      real*8 fxy

      fxy(a1, b1, a2, b2, u, v) = (v-b1)*(a2-a1) - (u -a1)*(b2-b1)
      

      nmax = max(abs(n%x), abs(n%y), abs(n%z))
      if(nmax .eq. abs(n%x)) then
!             use projection on to y-z plane
         x1 = p1%y
         y1 = p1%z
         x2 = p2%y
         y2 = p2%z
         x3 = p3%y
         y3 = p3%z
         x = p%y
         y = p%z
      elseif(nmax .eq. abs(n%y)) then
!            use projection on to x-z plane
         x1 = p1%x
         y1 = p1%z
         x2 = p2%x
         y2 = p2%z
         x3 = p3%x
         y3 = p3%z
         x = p%x
         y = p%z
      else
!             use projection on to x-y plane
         x1  =p1%x
         y1 = p1%y
         x2  =p2%x
         y2 = p2%y
         x3  =p3%x
         y3 = p3%y
         x  =p%x
         y = p%y
      endif

      p12 = fxy( x1, y1, x2, y2, x, y)
      p23 = fxy( x2, y2, x3, y3, x, y) 
      temp = p12*p23
      if(temp .lt. 0.) then
         icon = 4
         goto 100
      elseif(temp .eq. 0.) then
         if(p12 .eq. 0.) then
            icon = 1
         else
            icon = 2
         endif
         if(x1 .eq. x .and. y1 .eq. y) then
            jcon = 1
         elseif(x2 .eq. x .and. y2 .eq. y) then
            jcon = 2
         elseif(x3 .eq. x .and. y3 .eq. y) then
            jcon = 3
         else
            if(p12 .eq. 0.) then
               if( (x .ge. min(x1,x2) .and. x .le. max(x1,x2))
     *             .and.
     *          (y .ge. min(y1, y2) .and. y .le. max(y1, y2))) then
                  jcon =0
               else
                  icon = 4
               endif
            else
!                 p23 = 0
               if( (x .ge. min(x2,x3) .and. x .le. max(x2,x3))
     *             .and.
     *          (y .ge. min(y2, y3) .and. y .le. max(y2, y3))) then
                  jcon =0
               else
                  icon =  4
               endif
            endif
         endif
         goto 100
      endif
      p31 = fxy(x3,y3, x1, y1, x, y) 
      temp = p23*p31
      if(temp .lt. 0.) then
         icon = 4
         goto 100
      elseif(temp .eq. 0.) then
!            p31 should be 0         
         icon =3
         if( (x .ge. min(x1,x3) .and. x .le. max(x1,x3))
     *             .and.
     *      (y .ge. min(y1, y3) .and. y .le. max(y1, y3)) ) then
            jcon =0
         else
            icon =  4
         endif
         goto 100
      endif
      temp = p31* p12
      if(temp .lt. 0.) then
         icon = 4
      else
         icon = 0
      endif
 100  continue
      end

      subroutine episPin4(p1, p2, p3, p4, n, p, icon, jcon)
      implicit none
#include "ZepPos.h"
       type(epPos)::  p1 ! input.  p1,p2,p3,p4 form a "square" plane.
       type(epPos)::  p2 ! input.  any 3 of which are not colinear
       type(epPos)::  p3 ! input.  "square" may be skewed
       type(epPos)::  p4 ! input.  2 of them may be the same.
       type(epPos)::  n  ! input. normal vector of the plane
       type(epPos)::  p  ! input given point which is on the plane fomred by
                        !  three points p1, p2, p3, p4

      integer icon  ! output. 0: p is inside of the "square"
                    !         1: p is on  p1-p2
                    !         2: p is on p2-p3
                    !         3: p is on p3-p4
                    !         4: p is on p4-p1 
                    !         5: p is outside of the "square"
      integer jcon  ! output. has meaning if icon = 1~4
                    !         0: p is not on the vertex
                    !         i: p=pi (i=1~4)
      call episPin3(p1, p2, p3, n, p, icon, jcon)
      if(icon .eq. 4) then
         call episPin3(p1, p3, p4, n, p, icon, jcon)
         if(icon .eq. 4) then
            icon = 5
         elseif(icon .eq. 2) then
            icon = 3
            if(jcon .ne. 0) then
               jcon = 4
            endif
         elseif(icon .eq. 3) then
            icon = 4  ! jcon=4 or 1 will not happen
         endif   
      elseif(icon .eq. 3 .and. jcon .eq. 0) then
         icon = 0
      endif
      end

      subroutine ep3p2plane(p1, p2, p3, n, k)
      implicit none
!        3 points to palne;
!        get normal vector and const to represent a plane
!        formed by a  three non-colinear points
!
#include "ZepPos.h"
       type(epPos)::  p1  ! input. first point
       type(epPos)::  p2  ! input. second point
       type(epPos)::  p3  ! input. 3rd point.  They should not be on a line.
       type(epPos)::  n   ! output.a vector normal to the plane  formed by
                         !      p1,p2,p3.. it is a unit vector (from v9.136) 
      real*8 k           ! output. const to represent the plane.
                         ! vec(n).vec(x) = k where x is any point on the plane.

       type(epPos)::  p2p1, p3p2
       type(epPos)::  p3p1
      real(8):: norm
      call epdiff3vec(p2, p1, p2p1)
!      call epdiff3vec(p3, p2, p3p2)
      call epdiff3vec(p3, p1, p3p1)
!      call epvectorProd(p2p1, p3p2, n)
      call epvectorProd(p2p1, p3p1, n)
      call epscalerProd(n, n, norm)
      norm = sqrt(norm)
      n%x = n%x/norm
      n%y = n%y/norm
      n%z = n%z/norm
      call epscalerProd(n, p1, k)
      end
      
      subroutine epdiff3vec(p1, p2, p3)
      implicit none
#include "ZepPos.h"
       type(epPos)::  p1  ! input.  
       type(epPos)::  p2  ! input.  
       type(epPos)::  p3  ! output p3 = p1 - p2

      p3%x = p1%x - p2%x
      p3%y = p1%y - p2%y
      p3%z = p1%z - p2%z
      end
      subroutine epscalerProd(p1, p2, ans)
      implicit none
#include "ZepPos.h"
       type(epPos)::  p1  ! input.  
       type(epPos)::  p2  ! input.  
      real*8 ans  ! output ans = p1.p2

      ans = p1%x * p2%x + p1%y * p2%y +  p1%z * p2%z
      end

      subroutine epvectorProd(p1, p2, p3)
      implicit none
#include "ZepPos.h"
       type(epPos)::  p1  ! input.  
       type(epPos)::  p2  ! input.  
       type(epPos)::  p3  ! output. p1 x p2

      p3%x = p1%y*p2%z - p1%z* p2%y
      p3%y = p1%z*p2%x - p1%x* p2%z
      p3%z = p1%x*p2%y - p1%y* p2%x

      end
      subroutine epwhichside(nm, k, p, diff)
!        suppose a plane expressed by a normal vector nm and 
!        scaler value k (k = x*nm where x is any point on the plane)
!        this routine see on which side a given point p is located
!        with respect to the plane.  diff > 0 means p is on the
!        same side as the normal vector is directed.
!        < 0 means opposit side.  
!        0 means p is on the plane. 
!        If 3 points are given to form the plane and the order
!        they were given to ep3p2plane is anti-clock wise if seen from a one
!        side of the plane, diff > 0 means that the point is in the same side.
!        
      implicit none
      real(8),intent(in):: nm(3)
      real(8),intent(in):: k
      real(8),intent(in):: p(3)
      real(8),intent(out)::diff 

      real(8)::sp
      call epscalerProd(p, nm, sp)
      diff = sp - k
      end
      subroutine epwhichside0(p1, p2, p3, p, nm, k,  diff)
!        a plane is  expressed by p1, p2, p3.
!        This routine sees on which side a given point p is located
!        with respect to the plane.  diff > 0 means p is on the
!        same side as the normal vector (nm) is directed.
!        < 0 means opposit side.  
!        0 means p is on the plane. 
!        If the order of  p1, p2, p3 is anti-clock wise when
!        seen from a one  side of the plane, diff > 0 means that the 
!        point is in the same side.
!       

      implicit none
      real(8),intent(in):: p1(3), p2(3), p3(3)  ! input non colinear  3 points
      real(8),intent(in):: p(3)   ! a point to be examined
      real(8),intent(out):: nm(3) ! normal vector of the plane
      real(8),intent(out):: k !  scaler value k (k = x*nm where x is any point on the plane)

      real(8),intent(out)::diff 

      real(8)::sp
      call ep3p2plane(p1, p2, p3, nm, k)
      call epscalerProd(p, nm, sp)
      diff = sp - k

      end



