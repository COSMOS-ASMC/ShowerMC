!  kIsInPolyhedra0 is a subroutine to judge that a given point is inside
!  of a given polyhedra or not.
!  The test program here is to generate a number of points and judge
!  if they are inside or not. 
!  The program write for a one given point
!     inside  x y z
!  where inside=0 if x y z is inside of the polyhedra
!               1  //         on the surface //
!               2  //         outside of the //
!  To use the test program,
!  1) Uncomment main and readvertex program located at the last part
!  2) Fix the 'digitize' and 'rotate' in the main program.
!     If 'digitize' is to be used, its typical value will be 10.
!     Then, the each test point is digitized by int( digitize*x) / digitize.
!     This is usefull when the user wants to know what happens when
!     the point is exactly on the surface (including just on a vertex).
!     If such one is not needed, 0 may be given.
!     'rotate' is used to incline p and q polygon. 
!  3) Prepare a data file to represent  p and q.
!  4) Compile the program, run.  Output should be redirected to a file.
!  5) The error ouput consisting of 3 columns may by used to draw
!     the volume shape (say by gnuplot).


! Data preparation: 
!  Pprepare vertex list p and q  (each coplanar n verteces:
!  the last vertex and the 1st one will be connected to make the polygon
!  closed.)  The two polygon's normal vectors need not be directed to
!  the same direction.  Let the verteces of one polygon be p(3,n)
!  and the other q(3,n). (3 is for x,y,z).  The data file would contain
!  
!  p1x p1y p1z
!  p2x p2y p2z
!  ... 
!  pnx pny pnz
!                  one blank line
!  q1x q1y q1z
!  q2x q2y q2z
!  ...
!  qnx qny qnz
!
!Here,  n <=100.  The user may give a file name arbitrary, but it must
!  be copied to temp.dat.

! The volume is fomed by surfaces of polygon p and q and side walls
! formed by  connecting p(:,i) and q(:,i). One square like wall is
! generally not a plane, so two triagnle planes are formed by dividing
! the square by the diagonal line. Therefore, the lines formed by
! connecting p's and q's must not cross. 
!  If we follow the vertexes like,  p1->p2->...pn->p1,
!  and make the right-handed screw rotation globally 
!  the same as this rotation, the screw should go outwards
!  of the volume and this direction is defined as the normal
!  vector of the plane (p1->p2..->p1). Then, the normal
!  vector of the plane (q1->q2->..q1) would go inwards
!  of the volume.  If it is ok if this relation of two
!  normal vectors is inverse.
! 
! To test a special point for debugging purpose,  max index of the
!  do loop be 0.  And activate last do part of the main program.
!  After starting the run, give the (x,y,z) value.
! (Activate bugbug, bugbug2 and uncomment lines inbetween !////  lines
!  will be useful).
module modQisInpolyhed0
! A  simple polyhedron
  implicit none
  real(16),save::minmax(3,2)  ! (:,1) xmin, ymin, zmin
                             ! (:,2) xmax, ymax, zmax
  real(16),parameter::pi=3.14159265358979323846264338327950q0 ! asin(1.d0)*2
  real(16),parameter::twopi = 2*pi
!  logical,save::bugbug=.false.
!  logical,save::bugbug2=.false.
end module modQisInpolyhed0

 subroutine kQisInPolyhed0(p, q, n, x, inside)
   use modQisInpolyhed0
   implicit none
   integer,intent(in)::n  ! number of points
   real(16),intent(in)::p(3,n)  ! p(n) are vertex points on a plane
   !  No check is performed here for coplanarity; for that 
   !  use kCheckCoPlan.f90.
   real(16),intent(in)::q(3,n)  ! q(n) are another such one
      !  They form two surfaces of the volume of which other 
      !  surfaces are formed  as follows:
      !  p(:,i) is connected to q(:,i); then they will form  
      !  a number of ~  squares. (i.e., p(:,i) and q(:,i) 
      !  must be ~ nearest points among p and q.  Each square
      !  might not be a real square (some skew may exist).
      !  So the real surfaces are formed by a number of triangles.
      !  
      ! One p is connected to two q's to form a triangle 
      ! One q is //               p's  //
      !***  When the volume is seen from outside, p1,p2,...pn should
      ! globally rotate  counter-clock wise and the normal vector
      ! is directed to the outwords of the volume. (right-hand screw
      ! moving direction is the normal vector, rotation of the screw
      ! is the same as the  global order of  p1, p2...
      ! Hence, the normal vector of plane formed by qi is directed
      ! to the inwards of the volume. The opposit case is alos
      ! OK. 
      !  
      ! The side wall trianle points are so taken that their normal 
      ! vector is always directed to outwards. 
      !** p(n) --> p(1) (and q(n)-->q(1)) forms a edge of the volume.
   real(16),intent(in)::x(3)    ! given point. to be judged 
                               ! inside/outside of the volume
   integer,intent(out):: inside ! 0 inside but sometimes on the surface
                       ! 1 on the surface. should be regarded as inside 
                       ! 2 outside.
                       ! -1.  strange input.
   !   method. As described in S.Nagashima, Vol.27 No.7 1986
   !           Jyoho Shori Gakkai Ronbunn-shi p.744
   !   Suppose a sphere of radius 1 with the center being the given
   !   point X.
   !   For every plane forming the surface of the volume, 
   !   connect each of  verteces with X and get crossing-point with
   !   the sphere to form a figure. Get the area of the figure formed
   !   on the sphere. (Area is signed. If X is the same side of the
   !   normal vector of the surface, + else -. Normal vectors should
   !   be given in the same manner for all surfaces). If the sum of
   !   the area is -4pi, X is inside, if it is 0, outside.
   !   The area is simply calulated
   !   by computing the angle between traiangls formed by the surface
   !   verteces and X.
   !    
   real(16):: nvs(3)  ! normal vector of a surface
           ! formed by 3 points on the polyhedron.
           ! two triagnle areas are formed with X
   real(16):: ks ! nvs(:) *r(:) (r is arbitrayy  point on the plane)
   real(16):: nv(3)   ! normal vector of 3 planes formed by a
           !  triangle area and X:
   integer  ns, n1, n2, n3
   real(16):: k(3)  
   real(16):: area, suba
   real(16):: diff
   real(16):: v1(3), v2(3)
   integer:: i, i2, i3
   real(16):: ans, alpha, sumteta   
   integer:: icon, jcon
   if(n < 3 ) then
      write(0,*) &
       ' # of points =', n , ' too small for kisInPolyhed0'
      stop
   endif
   

   minmax(1,1)=min ( minval( p(1,:) ), minval( q(1,:)  ) )  ! xmin
   minmax(2,1)=min ( minval( p(2,:) ), minval( q(2,:)  ) )  ! ymin
   minmax(3,1)=min ( minval( p(3,:) ), minval( q(3,:)  ) )  ! zmin
   minmax(1,2)=max ( maxval( p(1,:) ), maxval( q(1,:)  ) )  ! xmax
   minmax(2,2)=max ( maxval( p(2,:) ), maxval( q(2,:)  ) )  ! ymax
   minmax(3,2)=max ( maxval( p(3,:) ), maxval( q(3,:)  ) )  ! zmax
!////////
!   write(0,*) 'min =',minmax(:,1)   
!   write(0,*) 'max=',minmax(:,2)
!//////////
   inside = 2   
   do i = 1, 3
      if( x(i) < minmax(i,1) ) then
         return
      elseif( x(i) > minmax(i,2) ) then
         return
      endif
   enddo

   area = 0.

   call kQisInPolyhed0Side(p, q, n, x, suba, inside)
!!!!!!!!!!!
!   write(0,*) ' aft kisInPoly .Side, suba=',suba, ' inside=',inside
!!!!!!!!!!!
   if(inside == 1 ) then
      return   ! *************
   endif
!////////////
!   write(0,*) ' side wall area=',suba
!   write(0,*) '-----------------'
!////////////
   area = area + suba



! for top cap
!         get outwards goint unit normal vect nv:
!         1 below request proper direction 
   call kQisInPolyhed0Cap(p, n, x, nv, ks,  suba, inside)
   if( inside /= 0 ) then
      return ! **********
   endif

   call epQwhichside(nv, ks, x, diff)
   suba =sign( suba -(n-2)*pi, -diff)
   area = area + suba
!//////////
!   write(0,*) ' suba by top =', suba
!//////////

!////////
!   write(0,*) '-------------'
!   write(0,*) " bottom "
!!!//////
   call kQisInPolyhed0Cap(q, n, x, nv, ks,  suba, inside)
   if( inside /= 0 ) then
      return ! **********
   endif

   call epQwhichside(nv, ks, x, diff)   ! diff is without -  for bottom
   suba =sign( suba -(n-2)*pi, diff)

   area = area + suba
   if( abs(area) < twopi ) then
!       area ~ +/- 0
      inside = 2
   else
!   area ~ 4pi or ~ -4pi
      inside = 0
   endif
!///////////
!   write(0,*) ' suba by bot =', suba
!      write(0,*) ' area 4 ', area, ' inside=',inside
!////////////

 end subroutine kQisInPolyhed0


 subroutine kQisInPolyhed0Side(p, q, n, x, suba, inside)
   use modisInpolyhed0
   implicit none
   integer,intent(in):: n
   real(16),intent(in):: p(3,n), q(3,n)
   real(16),intent(in):: x(3)
   real(16),intent(out):: suba
   integer,intent(out):: inside  ! inside  1  --> on the side wall
                                 ! inside  0 -->  furthe check needed

   real(16):: nv(3)
   integer:: i, i2
   real(16):: teta(3), ans
   character(5)::condi
!      side triangles formed by p and q
   suba = 0.
   do i = 1, n
       if( i < n ) then
         i2 = i + 1
      else
         i2 = 1
      endif
      nv(1) = -2. ! nv should be computed inside
      call kQProjAReaOfTriang(p(1,i),  q(1,i), q(1,i2), x, nv, &
           teta, ans, condi)
!////////
!      if(bugbug) then
!         write(0,*) ' i',i, ' condi =',condi, ' ans=',ans
!      endif
!///////////
      if(condi /= 'ok') then
         if( condi == 'onin' .or. condi == 'edge') then
            inside = 1
            return
         elseif(condi /= 'ng') then
            ans = 0
         else
            write(0,'(a, 1p,3g15.5)') 'strange data for kisInPolyhed0'//  &
                 ' p(1,i), q(1,i), q(1,i2)=',  p(1,i), q(1,i), q(1,i2)
            write(0,*) 'i,i2=', i, i2
            stop
         endif
      endif   
      suba = suba + ans
      nv(1) = -2. ! nv should be computed inside
      call kQProjAReaOfTriang(p(1,i),  q(1,i2), p(1,i2), x, nv,  &
           teta, ans, condi)
!////////
!      if(bugbug) then
!         write(0,*) ' i s',i, ' icon =',condi, ' ans=',ans
!     endif
!///////////
      if(condi /= 'ok') then
         if( condi == 'onin' .or. condi == 'edge') then
            inside = 1
            return
         elseif(condi /= 'ng') then
            ans = 0
         else
            write(0,'(a, 1p,3g15.5)') 'strange data for kisInPolyhedra0'//  &
                 ' p(1,i), q(1,i), q(1,i2)=',  p(1,i), q(1,i), q(1,i2)
            write(0,*) 'i,i2=', i, i2
            stop
         endif
      endif   
      suba = suba + ans
   enddo
   inside = 0
 end subroutine kQisInPolyhed0Side

 subroutine kQisInPolyhed0Cap(p, n, x, nv, ks, suba, inside)
   use modisInpolyhed0
   implicit none
   integer,intent(in):: n
   real(16),intent(in):: p(3, n)
   real(16),intent(in):: x(3)

   real(16),intent(out):: nv(3), ks
   real(16),intent(out):: suba
   integer,intent(out):: inside ! if inside=1, x is on the plane
                             !   if inside =2, x is outside
                             !   if 0, further check needed 
   

   integer:: jcon, i, i2, i3
   real(16):: ans, sumteta
   character(5):: condi


   call kQgetNormalVec(p, n, 1, nv, ks, jcon)
   if( jcon /= 0 ) then
      write(0,'(a)') &
      ' Stranage input to kisInPolyhed0Cap: p is not coplanar ?:'
      write(0,'(1p, 3g16.7)' ) p(:,:)
      stop
   endif
!/////////
!   write(0,*) ' nv=',nv, ' ks = ',ks
!//////////////
   suba = 0.
   call kQgetAngleBtw2Planes(1, x, p(1,1), p(1,2), p(1,n), nv, ans, condi)
   if( condi /= 'ok' ) then
      select case(condi)
      case('x1-2',  'x2-3') 
         inside = 1
      case('ng')   
         write(0,*) ' strange in kisInPolyhed0'
         stop
      case default
         ! on the p plane.
         call kQinout3(p, n, x, nv, sumteta, condi)
         select case(condi)
         case('in', 'edge','vtx')   
            inside = 1
         case('out') 
            inside = 2
         case default
            write(0,*) 'stragne 111, condi from kinout3=',condi
            stop
         end select
      end select
      return
   endif
!/////////
!   write(0,*) 'T 1 angle  =',ans, 'condi=',condi
!////////////
   suba = suba + ans
   do i = 2, n
      i3 = i-1
      if(i /= n) then
         i2 = i+ 1
         call kQgetAngleBtw2Planes(2, x, p(1,i), p(1,i2), p(1,i3),  &
         nv, ans,condi)
!/////////
!         write(0,*) 'T 2 angle for i =', i,ans, 'condii =',condi
!////////////
      else
         i2 = 1
         call kQgetAngleBtw2Planes(3, x, p(1,i), p(1,i2), p(1,i3), & 
              nv, ans, condi)
! /////////
!   write(0,*) 'T 3 angle for i =', i,ans, ' condi=',condi
!////////////
      endif
      if( condi /= 'ok' ) then
         select case(condi)
         case('x1-2',  'x2-3') 
            inside = 1
         case('ng')   
            write(0,*) ' strange in kisInPolyhed0'
            stop
         case default
         ! on the p plane.   ! will not happen
            call kQinout3(p, n, x, nv, sumteta, condi)
            select case(condi)
            case('in', 'edge','vtx')   
               inside = 1
            case('out') 
               inside = 2
            case default
               write(0,*) 'stragne 222, condi from kinout3=',condi
               stop
            end select
         end select
         return
      endif
      suba = suba + ans
   enddo
 end subroutine kQisInPolyhed0Cap

! program main
!   use modisInpolyhed0
!   implicit none
!   real(8),parameter:: digitize=50.  ! give 0 if don't want so
!   logical,parameter:: rotate =.false.
!
!
!   real(8):: ux, uy, uz
!   integer:: inside 
!   integer:: i
!   integer,parameter::m = 100
!   integer:: n
!   real(8):: p(3,m), q(3,m)
!   real(8):: pd(3,m), qd(3,m)
!
!   real(8):: x(3), xd(3)
!
!   real(8):: cost, sint
!   real(8):: rm(3,3), rm1(3,3), rm2(3,3)
!   real(8):: xmin, ymin, xmax, ymax, zmin, zmax
!
!   call readvertex(p, q, m,  n)
!
!   cost=cos(pi/4.)
!   sint=sin(pi/4.)
!   call cgetRotMat3(1, cost, sint, rm)
!   cost = cos(pi/3.)
!   sint = sin(pi/3.)
!   call cgetRotMat3(3, cost, sint, rm1)
!   call cmultRotMat3(rm, rm1, rm2)
!   cost = cos(-pi/5.) 
!   sint = sin(-pi/5.)
!   call cgetRotMat3(2, cost, sint, rm1)
!   call cmultRotMat3(rm2, rm1, rm) 
!
!   do i = 1, n
!      call capplyRot3(rm, p(1,i), pd(1,i))
!   enddo
!
!   cost = cos(-pi/4.) 
!   sint = sin(-pi/4.)
!   call cgetRotMat3(2, cost, sint, rm1)
!   call cmultRotMat3(rm2, rm1, rm) 
!   do i = 1, n
!      call capplyRot3(rm, q(1,i), qd(1,i))
!   enddo
!   if(.not. rotate) then
!      pd(:,:) = p(:,:)
!      qd(:,:) = q(:,:)
!   endif
!
!   do i = 1, n
!      write(0,'(1p,3g15.5)') pd(:,i)
!   enddo
!   write(0,'(1p,3g15.5)') pd(:,1)
!   write(0,*)
!   do i = 1, n 
!      write(0,'(1p,3g15.5)') qd(:,i)
!   enddo
!   write(0,'(1p,3g15.5)') qd(:,1)
!
!   xmin = min( minval(pd(1,1:n)), minval(qd(1,1:n)) ) -0.25
!   xmax = max( maxval(pd(1,1:n)), maxval(qd(1,1:n)) ) + 0.25
!   ymin = min( minval(pd(2,1:n)), minval(qd(2,1:n)) ) -0.25
!   ymax = max( maxval(pd(2,1:n)), maxval(qd(2,1:n)) ) + .25
!   zmin = min( minval(pd(3,1:n)), minval(qd(3,1:n)) ) -0.25
!   zmax = max( maxval(pd(3,1:n)), maxval(qd(3,1:n)) ) + 0.25   
!   write(0,*) xmin, xmax
!   write(0,*) ymin, ymax
!   write(0,*) zmin, zmax
!
!   do i = 1,  100000
!      call rndc(ux)
!      call rndc(uy)
!      call rndc(uz)
!      x(1) = (xmax-xmin)*ux + xmin
!      x(2) = (ymax-ymin)*uy + ymin
!      x(3) = (zmax-zmin)*uz + zmin
!      if( digitize > 0. ) then
!         x(1) = int( digitize*x(1) ) / digitize
!         x(2) = int( digitize*x(2) ) / digitize
!         x(3) = int( digitize*x(3) ) / digitize
!      endif
!      call kisInPolyhed0(pd, qd, n, x, inside)
!      write(*,'(i3,1p,3g14.5)') inside, x(:)
!   enddo
!
!
!!  do 
!!      write(0,*) ' Enter x '
!!      read(*,*) x   !7.3000        6.8000        3.9000
!!      if(x(1) == -100.) stop
!!      write(0,*) x
!!      call kisInPolyhed0(p, q, n, x, inside)
!!      write(0,*) inside
!!  enddo
! end program main
! subroutine readvertex(p, q, m, n)
!  implicit none
!  integer,intent(in)::m
!  integer,intent(out)::n
!  real(8),intent(out)::p(3,m), q(3,m)
!
!  real(8):: x(3)
!  character(80):: input
!  integer::i
!  integer,parameter::fn=11
!
!  open(fn,file="temp.dat")
!  n = 0
!  do
!     read(fn,'(a)') input
!     if( input == "#") cycle 
!     if(input == '' ) exit
!     n = n + 1
!     if(n > m) then
!        write(0,*) ' too many vertex'
!        stop
!     endif
!     read(input, *) p(:,n)
!  enddo
!  do i = 1, n
!     read(fn,'(a)') input
!     if( input == "#") cycle 
!     read(input,*) q(:, i)
!  enddo
!
! end subroutine readvertex
