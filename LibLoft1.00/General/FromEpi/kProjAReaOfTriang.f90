 subroutine kProjAReaOfTriang(p1, p2, p3, x, nv, teta, area, cond)
!     Points p1,p2,p3 form a triangle area. Suppose a unit sphere
!     of which the center is at x, and projected area of the 
!     triagle formed by connecting x and p_i.  This
!     computes the sureface area of the projected triagle
!     which is the Sum(teta_i) -(n-2)pi = teta1+teta2+teta3-pi
!     (n=3). teta_i is the angle between two triangle surfaces
!     formed by x and p_i's 
   ! If a right-handed screw is rotated as p1->p2->p3->p1
   ! the screw moving direction should be the normal 
   ! vector direction which goes outwards from the surface.
   implicit none
   real(8),intent(in):: p1(3)  !  vertex of a triangle
   real(8),intent(in):: p2(3)  !  vertex of the triangle
   real(8),intent(in):: p3(3)  !  vertex of the triangle
   real(8),intent(in):: x(3)   ! point to be judeged
   real(8),intent(inout):: nv(3)   ! outwards going normal
             ! unit vector of the surface formed by p_i. 
             !  (or surface containing all p_i's)
             ! If nv(1) < -1, nv is to be computed newly
             ! here else  assumed that  the value is
             ! alredy computed there.
   real(8),intent(out):: teta(3)  !  angle between two surfaces 
                  !   having the common line x-pi.
   real(8),intent(out):: area  ! signed  area
   character(*),intent(out):: cond  ! length >= 5
               ! 'ok'   area obtained.
               ! 'ng'   p1 p2 p3 do not form a triagle
               ! 'edge'  x is on the triangle plane  and
               !       on one of the p1-p2, p2-p3 or p3-p1 segment 
               ! 'onin'    x is on the triagle plane and inside of 
               !       of the triangle area
               ! 'onout'  x is on the triagle plane but outside of
               !       the triagle 
   real(8),parameter::pi=3.14159265358979323846d0 ! asin(1.d0)*2
   real(8):: ks, diff, alpha
   integer:: jcon
   character(5):: condi

   if( nv(1) < -1.0 ) then
!            outwards going normal unit vector of the surface
!          formed by p1,p2,p3
      call ep3p2plane0(p1, p2, p3, nv, ks, alpha, jcon)
      if( jcon /= 0 ) then
         cond = 'ng'
         area = 0.
         return
      endif
   endif
   call kgetAngleBtw2Planes(1, x, p1, p2, p3, nv, teta(1), condi)
!///////
!   if(bugbug2) then
!      write(0,*) ' teta(1)=', teta(1), ' condi=', condi
!   endif
!////////

   if( condi /= 'ok' ) then
      select case(condi)
      case('x-12-', 'x-13-') 
         cond = 'onout'   ! in this case, x is finally outside
      case('x1-2', 'x1-3', 'x2-3') 
         cond = 'edge'
      case default
         cond = condi
      end select
      area = 0.
      return
   endif

   call kgetAngleBtw2Planes(2, x, p2, p3, p1, nv, teta(2), condi)
!///////
!   if(bugbug2) then
!      write(0,*) ' teta(2)=', teta(2), ' condi=', condi
!   endif
!////////
   if( condi /= 'ok' ) then
      select case(condi)
      case('x-12-', 'x-13-') 
         cond = 'onout'   ! in this case, x is finally outside
      case('x1-2', 'x1-3', 'x2-3') 
         cond = 'edge'
      case default
         cond = condi
      end select
      area = 0.
      return
   endif

   call kgetAngleBtw2Planes(3, x, p3, p1, p2, nv, teta(3), condi)
!///////
!   if(bugbug2) then
!      write(0,*) ' teta(3)=', teta(3), ' condi=',condi
!   endif
!////////
   if( condi /= 'ok' ) then
      select case(condi)
      case('x-12-', 'x-13-') 
         cond = 'onout'   ! in this case, x is finally outside
      case('x1-2', 'x1-3', 'x2-3') 
         cond = 'edge'
      case default
         cond = condi
      end select
      area = 0.
      return
   endif

   call epwhichside(nv, ks, x, diff)
   area = sign(sum(teta(:)) -pi, -diff)
   cond = 'ok'
 end subroutine kProjAReaOfTriang
