 subroutine kisInsidePolVol0(p, q, n, x, inside)
   implicit none
   integer,intent(in)::n  ! number of points
!   type point3D
!      real(8)::xyz(3)
!   end type point3D
!
!   type(point3D),intent(in)::p(n) ! n points on a plane
!   type(point3D),intent(in)::q(n) ! n points on a plane
!          !        p(1), p(2) are anti clockwise
!          !        q(1),..         //
!          !
   real(8),intent(in)::p(3,n)  ! p(n) is on a circle. 
   real(8),intent(in)::q(3,n)  ! q(n) is //
                   !  each circle is approximated by p,q to
      !form surfaces. The side wall is formed by connecting
      !  pi and qi
      !  Let's take the normal vector of a plane as directed 
      !  outwards of the volume. Then,  pi must round left-rotating
      !  screw while qi righ-rotating one. 
      !  
      !  Triangles formed by connecting (p1 q2 q1 p1),
      ! (p1 p2 q2 p1) ... (pi qi+1, qi,pi)  (pi,pi+1,qi+1,pi) 
      ! forms side wall. The points must round counter-clock wise. 
      ! 
   real(8),intent(in)::x(3)    ! given point. to be judged 
                               ! inside of the volume or not
   integer,intent(out):: inside ! 0 inside
                       ! 1 on the surface.
                       ! 2 outside
   
   real(8):: nv(3)   ! normal vector of plane (outward going)
   real(8):: k       ! nv*p (p is on the plane)
   real(8):: diff
   
   integer:: n1, n2, n
   integer:: l1, l2
!     form a triangle reprsenting a plane  by pi
!    take 3 points 
   if(n < 3 ) then
      write(0,*) ' # of points =', n , ' too small for kisInsidePolVol0'
      stop
   endif
                     ! n 3     4    5   6  7  8  9
   n1 = 1            !   1     1    1   1  1  1  1
   n2 = 1 + n/3      !   2     2    2   3  3  3  4
   n3 = 1 + 2*n/3    !   3     3    4   5  5  6  7
!       use n1, n2, n3-th points to form triangle 
!          get normal vector
   call ep3p2plane(p(1,n3), p(1,n2), p(1,n1), nvp, k)
   call epwhichside(nvp, k, x, diff)
   if( diff > 0.d0 ) then
      inside = -1
      return    ! *************
   elseif( diff == 0.d0 ) then
      ! should be on the plane; see if x's  projection is
      ! inside p polygon
      call epselPerpenPlane(nvp, l1, l2) 
!         see if x is inside a polygon projected on l1-l2 plane
      call kinout( p(l1, 1), 3, p(l2,1), 3, n, x(l1), x(l2), eps, 'conv', icon)    
      if(icon == 1 ) then
         inside = 2
      elseif( icon == -1 .or. icon == 0 ) then
         inside = 1
      else
         write(0,*) 'kisInsidePolVol0(p, q, n, x, inside)'
         stop
      endif
   endif
!     same for qi
!     get normal vector
   call ep3p2plane(q(1,n1), q(1,n2), q(1,n1), nvq, k)
   call epwhichside(nvq, k, x, diff)
   if(diff > 0.d0 ) then
      inside = -1
      return    ! *************
   elseif( diff == 0.d0 ) then
      ! should be on the plane; see if x's  projection is
      ! inside p polygon
      call epselPerpenPlane(nvq, l1, l2) 
!         see if x is inside a polygon projected on l1-l2 plane
      call kinout( q(l1, 1), 3, q(l2,1), 3, n, x(l1), x(l2), eps, 'conv', icon)    
      if(icon == 1 ) then
         inside = 2
      elseif( icon == -1 .or. icon == 0 ) then
         inside = 1
      else
         write(0,*) 'kisInsidePolVol0(p, q, n, x, inside)'
         stop
      endif
   endif
!         x is inbetween two circle planes
!    see if inside  the wall 
   do i = 1, n

 end subroutine kisInsidePolygonVol

 subroutine epselPerpenPlane(nvp, l1,l2) 
! select a plane (x-y, y-z, z-x) which is 
! most close to the perpendicular plane with
! respect to a given  vector nvp
   implicit none
   real(8),intent(in):: nvp(3) ! given  vector
   integer,intent(out):: l1 
   integer,intent(out):: l2   ! l1=1,l2=2 for x-y
                          !   l1=2 l2=3 for y-z
                          !   l1=3 l2=1 for z-x  
   integer:: lmax, lmin, l
   lmax = maxloc(nvp)   ! which has max value
   lmin = minloc(nvp)    !which has min value
   if( abs( nvp(lmax)) > abs(nvp(lmin) ) ) then
         ! get really |max|
      l = lmax
   else
      l = lmin
   endif
   l1 = l + 1               ! l     l1   l2 
   if(l1 > 3 ) l1 = 1       ! 3     1     2  x-y
   l2 = l - 1               ! 1     2     3  y-z
   if(l2 < 1 ) l2 = 3       ! 2     3     1  z-x
 end subroutine epselPerpenPlane
