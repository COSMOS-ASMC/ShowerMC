subroutine kxplPolyhed0(p, q, n, convex, &
     pos, dir, el, insidevol, crosspos)
!    find the  crossing point (xp) of a line with a given polyhedra
!    (not general polyhedra). xp is the nearest one.
!
  implicit none
  integer,intent(in):: n   ! # of vertexes in p, q
  real(8),intent(in):: p(3,n)  ! vertex list.  
  real(8),intent(in):: q(3,n)  ! vertex list. 
  integer,intent(in):: convex  ! 1--> planes by p and q are voth convex
                               ! 0--> p and/or q concave. 
     ! for p,q see kisInPolyhed0.f90.  p(1)->p(2)->.. p(n)->p(1)
     ! must be coplaner. (q too).
 ! 3D volume is formed these two surfaces and side wall
 ! fomred by connecting p(i) and q(i). Further, the 4 vertexes
 ! formed at i,i+1 of p and q make a square (not exactly plane) 
 ! which is divided into two triagles; these trianles make the
 ! wall of the polyhedron. 

  real(8),intent(in):: pos(3) ! the line passes this point
  real(8),intent(in):: dir(3) ! with this direction cos
  real(8),intent(out)::el ! length from pos to the nearest
                           ! xp.

  integer,intent(out):: insidevol !  0: pos is inside
                               !  1: pos is outside
  integer,intent(out):: crosspos  !  1: xp is on the p-plane
                               !  6: xp is on the q-plane
                               !  2: xp is side wall.
                               !  -1: no xp.
  real(8):: normalv(3), a, b, c, d  
  integer,parameter::req = 1  ! request so that 
                            !  normal vector is outgoing
  integer::icon
  real(8),parameter::pi=3.14159265358979323846d0 ! asin(1.d0)*2
  real(8):: xp(3), ela(3), eltemp, sumteta
  character(len=4):: cond
  integer::nx, icon2a(3), ll(1)
!  real(8), parameter:: elmin = -1.d-5
  real(8), parameter:: elmin = 0.
  integer::im(1)

  nx = 0  ! # of x.p
  crosspos=-1
!!!!!!!!!
!  write(0,*) '-------------------'
!  write(0,*) ' top of kxplPoly '
!!!!!!!!!!!!!
!      judge if a given 'pos' is inside
  call kisInPolyhed0(p, q, n, pos, insidevol)
!!!!!!!!!!!!
!  if( insidevol == 0 ) then
!     write(0,*) 'pos', pos, ' is inside since insidevol=',insidevol
!  endif
!!!!!!!!!!
  ! 0->inside.  1-->on the surface,  2--> outside.
  !     we chage it as 1-->0,  2-->1
  if(insidevol == 1) then
     insidevol = 0
  elseif(insidevol == 2) then
     insidevol = 1
  endif
!   get plane formed by p's
  call  kgetNormalVec(p, n, req, normalv, d, icon)
  if(icon /= 0 ) then
     write(0,*) ' error input to kgetNormalVec in kxplPolyhed0'
     write(0,*) ' n=', n
     write(0,* )' px', p(1,1:n)
     write(0,* )' py', p(2,1:n)
     write(0,* )' pz', p(3,1:n)
     stop
  endif
!   ax +bx + cx = d is the p plane
  a = normalv(1)
  b = normalv(2)
  c = normalv(3)
!     seee xp with the plane
  call kxplp(pos(1), pos(2), pos(3), dir(1), dir(2), dir(3), &
      a, b, c, d,  el, icon)
!!!!!!!!
!  write(0,*) ' aft kxplp for p plane, icon=',icon, ' el=',el
!!!!!!!!!
  if( icon == 0 .and. el > 0. ) then
     ! xp found.   see if it is inside the pi's 
     xp(:) = pos(:) + el*dir(:)
     call kinout3(p, n, xp, normalv, sumteta, cond)
     if( abs(sumteta) > pi )  then
        ! inside (if outside ~0) )
        if( insidevol == 0 .and. convex == 1 ) then
           crosspos = 6
           return
        else
           nx = nx + 1
           icon2a(nx) = 6
           ela(nx) = el
        endif
     endif
  endif
!!!!!!!!
!  write(0,*) ' aft kxplp for p plane, crosspos=',crosspos
!!!!!!!!!

!   get plane formed by q's
  call  kgetNormalVec(q, n, req, normalv, d, icon)
  if(icon /= 0 ) then
     write(0,*) ' error input to kgetNormalVec in kxplPolyhedra0'
     write(0,*) ' n=', n
     write(0,*) ' qx', q(1,1:n)
     write(0,*) ' qy', q(2,1:n)
     write(0,*) ' qz', q(3,1:n)
     stop
  endif
!   ax +bx + cx = d is the q plane
  a = normalv(1)
  b = normalv(2)
  c = normalv(3)
!     see xp with the plane
  call kxplp(pos(1), pos(2), pos(3), dir(1), dir(2), dir(3), &
       a, b, c, d, eltemp, icon)
!!!!!!!!
!  write(0,*) ' aft kxplp for q plane, icon=',icon, ' eltemp=',eltemp
!!!!!!!!!
  if( icon == 0 .and. eltemp > 0. ) then
     xp(:) = pos(:) + eltemp*dir(:)
     ! xp found.   see if it is inside q
     call kinout3(q, n, xp, normalv, sumteta, cond)
     if( abs(sumteta) > pi )  then
       ! inside (if outside ~0) )
!!!!!!
!           write(0,*) ' for q pl: inside '
!           write(0,*) ' insidevold =',insidevol, 'eltemp=',eltemp
!!!!!!!!!!
         if( insidevol == 0 .and. convex==1 ) then
            crosspos = 1
            el = eltemp
            return
         endif
         nx = nx + 1
         icon2a(nx) = 1
         ela(nx) = eltemp
     endif
  endif
!!!!!!!!
!  write(0,*) ' bef  side chk, crosspos=',crosspos,' el=',el, ' nx=',nx
!!!!!!!!!

  call kxplPolyhed0Side( p, q, n, convex, pos, dir, eltemp, &
       insidevol, icon)
!!!!!!!!
!  write(0,*) ' aft  side chk, insidevol=',insidevol, &
!       ' icon =',icon, ' eltemp=',eltemp
!!!!!!!!!

  if(icon /= -1 .and. eltemp > 0. ) then
     nx = nx + 1
     icon2a(nx) = icon
     ela(nx) = eltemp
  endif
  if( nx > 0 ) then
     ll = minloc( ela(1:nx) )
     el = ela(ll(1))
     crosspos = icon2a(ll(1))
  endif
!  write(0,*) &
!   ' at exit of kxpl, insicon1=',insidevol, ' icon2=',crosspos
!  if(crosspos > 0 ) then
!     write(0,*) ' el=',el, ' isidevol=',insidevol, ' crosspos=',crosspos
!  endif
!  write(0,*) '//////////////////////'
!!!!!!!!!!!!
end subroutine kxplPolyhed0

 subroutine kxplPolyhed0Side( p, q, n, convex, pos, dir, el, icon1, &
 icon2)
!         get xp on the side wall of polyhedra0
      implicit none
  integer,intent(in):: n   ! # of vertexes in p, q
  real(8),intent(in):: p(3,n)  ! vertex list.  see kxplPolyhedra
  real(8),intent(in):: q(3,n)  ! vertex list. p
  integer,intent(in):: convex  ! 1-->p and q are convex
                               ! 0-->p and/or  q are concave
  real(8),intent(in):: pos(3) ! the line passes this point
  real(8),intent(in):: dir(3) ! with this direction cos
  real(8),intent(out)::el ! length from pos to the nearest
                           ! xp.
  integer,intent(in):: icon1   ! 0 if pos is inside  polyhedra else 1
  integer,intent(out)::icon2   ! -1  if no x.p.  else 2  

  integer i, j, icon, nx, icons, l1, l2
  real(8):: eltemps, towhichs
  real(8):: eltemp, towhich
!         form a triangle using p, q
!        p(j)---- p(j+1)
!        /   \     /
!       /     \   /
!     q(j)     \ /
!          -- -  *   q(j+1)
  el=1.0d20
  icon2 = -1
  nx = 0
  do j = 1, n
     if( j < n) then
        l1=j
        l2=j+1
     else
        l1=n
        l2=1
     endif
!   for the 10^6  completely the smae random poins and directions
!     by expLand3vpB: 
!        real 1m34.135s
!        user 1m31.840
!        sys  0m2.189s
!     by kxpland3vpB
!        real 4m14.346s
!        user 4m11.589s
!        sys  0m2.446s
     call  epxpLand3vpB(p(1,l1), q(1,l1), q(1,l2),  &
             pos, dir, eltemp,  towhich, icon)
     if(icon < 4 .and. eltemp >0.) then
!       call  kxpLand3vpB(q(1,j), p(1,j), q(1,j+1),  &
!             pos, dir, eltemp,  towhich, icon)
!       if(icon < 4 .and. icon >=0 .and.   eltemp >0. ) then
        nx = nx + 1
        el=min(eltemp, el)
        icon2 = 2
        if(convex == 1) then
           if( icon1 == 0 ) exit
!              if( towhich < 0. ) exit  ! from outside
   !           if( nx == 2 ) exit
        endif
     endif
     call  epxpLand3vpB(q(1,l2), p(1,l2), p(1,l1),  &
             pos, dir, eltemp,  towhich, icon)
     if(icon < 4 .and. eltemp >0.) then
!        call  kxpLand3vpB(q(1,j+1), p(1,j), p(1,j+1),  &
!             pos, dir, eltemp,  towhich, icon)
!        if(icon < 4 .and. icon >=0 .and.   eltemp >0. ) then
        nx = nx + 1
        el =min( el, eltemp)
        icon2 = 2
        if(convex == 1) then
           if( icon1 == 0 ) exit
 !             if( towhich < 0. ) exit  ! from outside
        endif
     endif
  enddo
end subroutine kxplPolyhed0Side
