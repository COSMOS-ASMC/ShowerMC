 subroutine kgetAngleBtw2Planes(i, x, p1, p2, p3, nv, teta, cond)
   implicit none

   integer,intent(in):: i  ! see the calling seq. above
           ! let n1 be the normal vector formed by x->p1->p3-x
           ! and n2 bet the //                     x->p1->p2-x
           ! i=1-> n1 n2 are newly computed and n1 is saved, n2 is saved
           ! i=2-> use -n2(saved) as n1 and compute n2 newly
           ! i=3-> use -n2(saved) as n1 and use -n1(saved) as n2 

   real(8),intent(in):: p1(3)  !  vertex of a triangle
   real(8),intent(in):: p2(3)  !  vertex of the triangle
   real(8),intent(in):: p3(3)  !  vertex of the triangle
   real(8),intent(in):: x(3)   !  point to see p_i's
   real(8),intent(in):: nv(3)   !unit outwards going normal of
                         ! triangle plane formed by p.
       ! in some case, p may be collinear as a part of a polygon.
       ! In such a case, nv should be the normal vec. of such a
       ! polygon.
   real(8),intent(out):: teta   !   rad. angle obtained.
   character(*),intent(out):: cond ! must be len is at least 5.
             !  'ok'.     teta obtained. x and p1-p2-p3 are not
             !            coplanar 
             !  'ng'      strange input. (p1=p2 etc)

             !         Other cases (below): x and p_i are coplanar. 
             !  'x1-2'. x is on the line p1--p2.  x is inbetween p1,p2.
             !         0<= teta <=1.
             !  'x-12-'. x is on the line p1--p2 but outside of p1-p2.
             !            teta<0 or teta > 1. 
             !  'x1-3'   x is on p1--p3
             !           teta has same meaning as above.
             !  'x-13-'  see x-12-.
             !  'x2-3'   x is on p2--p3    teta=pi
             !  'onin'   x is on the plane and inside of triagnle p1-p2-p3; 
             !           teta=pi.
             !  'onout'  x is on the plane  but outside of the 
             !           triangle.  teta = 0.


!                  p2
!                  
!         p3         p1
!
!           x
   real(8),parameter::pi=3.14159265358979323846d0 ! asin(1.d0)*2
   real(8),parameter::twopi = 2*pi

   real(8):: n1(3) 
   real(8),save:: n2(3)
   real(8),save::n1save(3)
   real(8):: n1xn2(3), r(3, 3), sumteta
   real(8):: ks, alpha, n1n2
   character(4)::condi
   integer:: jcon

   real(8),parameter::eps= 1.d-10

   cond = 'ok'
   if( i == 1 ) then
! normal vector of the  triangle surface formed by  x p1 p3 x   
      call ep3p2plane0(x,   p1, p3, n1, ks, alpha, jcon)  !
      if( jcon == 1 .or. jcon== 2 ) then
           !  p1/= p3 but x p1 p3 are collinear and  x is on the
           !  p1-p3 segment. jcon=1==>  0<= alpha <= 1
           !  jcon=2==>   alpha<0 or alpha > 1
         teta = alpha
         if( alpha >=0 .and. alpha <= 1.0) then
            cond = 'x1-3'
         else
            cond = 'x-13-'
         endif
         return
      elseif( jcon == -1) then
         write(0,*) ' strange in kgetAngleBtw2Planes:::'
         write(0,*) 'x=',x(:)
         write(0,*) 'p1=',p1(:)
         write(0,*) 'p3=',p3(:)
         teta =0.
         cond='ng'
         return
      endif
      n1save = n1
   else
      n1 = -n2
   endif
! normal vector of the  triangle surface formed by  x p1 p2 x   
   if( i /= 3 ) then
      call ep3p2plane0(x,   p1, p2, n2, ks, alpha, jcon)  ! now outwards
      if( jcon == 1 .or. jcon == 2 ) then
           !  p1/= p2 but x p1 p2 are collinear and  x is on the
           !  p1-p2 segment. jcon=1==>  0<= alpha <= 1
           !  jcon=2==>   alpha<0 or alpha > 1
         teta = alpha
         if( alpha >= 0. .and. alpha<=1.0) then
            cond = 'x1-2'
         else
            cond = 'x-12-'
         endif
         return
      elseif( jcon == -1) then
         write(0,*) ' strange in kgetAngleBtw2Planes:::'
         write(0,*) 'x=',x(:)
         write(0,*) 'p1=',p1(:)
         write(0,*) 'p2=',p2(:)
         cond = 'ng'
         return
      endif
   else
      n2 = -n1save
   endif
   n1n2 =max(min( dot_product(n1,n2), 1.d0), -1.d0)
   teta = acos(n1n2)
!/////////
!   if(bugbug2) then
!      write(0,*) ' n1=',n1(:)
!      write(0,*) ' n2=',n2(:)
!      write(0,*) ' n1n2=',n1n2
!      write(0,*) ' teta=',teta, ' acos=',acos(n1n2)
!   endif
!////////////
!   we may proceed without following job if we don't
!   pay attention to that
!   the point is  exactly on the cap surface
   if( teta <1.d-7  .or. (abs(teta-pi) < eps )  ) then
      ! this is not enough to judge that x and p's are
      ! coplanar, (think the case of 3 p's are collinear;
      ! otherwise, this is coplanar condtion).
      if(abs(sum(n1(:)+nv(:))) < eps .or. &
         abs(sum(n2(:)+nv(:))) < eps ) then
         ! x is on the plane
         ! see if inside of the triangle
         r(:,1) =  p1(:)
         r(:,2) =  p2(:)
         r(:,3) =  p3(:)
         call kinout3(r, 3, x, n1, sumteta, condi)
!///////////
!         if(bugbug2) then
!            write(0,*) ' kinout3 condi=',condi, ' sumteta=',sumteta
!         endif
!////////////
         if( condi == 'in' ) then
            cond = 'onin'
            teta = 0.
         elseif( condi == 'vtx' .or. condi =='edge' ) then
            cond = 'x2-3'
            teta = pi
         elseif(condi  == 'out' ) then
            !  x is outside of the triangle
            cond = 'onout'
            teta =  pi
         else
            cond='ng'  ! strange input
         endif
         return
      endif
   endif
!        vector product of n1, n2
   call epvectorProd(n1, n2,  n1xn2)
!    if n1xn2  and nv have same direction, 
   if( dot_product( n1xn2, nv) > 0. ) then
!     then,  teta should be  2pi - teta
      teta = twopi - teta
   endif
   cond = 'ok'
 end subroutine kgetAngleBtw2Planes
