subroutine ep3p2plane0(x, p1, p2, n, k, alpha, icon)
  implicit none
!        3 points to palne; 
!        get normal vector and const to represent the plane
!        formed by  three non-colinear points
!  diff from ep3p2plane is that this checks 3 points are colinear
!  or not. If normal vector cannot be obtained, non zero icon is given.
!
      real(8),intent(in):: x(3), p1(3), p2(3)
      real(8),intent(out):: n(3) ! unit normal vector.
      real(8),intent(out):: k ! output. const to represent the plane.
           ! vec(n).vec(r) = k where r is any point on the plane.
           ! k is the distance from the origin to the plane.
      real(8),intent(out):: alpha ! if icon =1,  alpha is such that
                      !  (x-p1) = alpha*(p2-p1)
                      ! alpha< 0: x is far from p2 than p1
                      !      0  : x = p1
                      !   0~ 1  : x is between p1 and p2
                      !      1  : x is p2
                      !     >1  : x is far from p1 than p2    
      integer,intent(out):: icon  ! if 0, n, k obtained.
                  ! if  1,  x, p1, p2 are colinear. but  p1 /= p2 
                  !    -1,  p1 == p2 alpha has no meaning

      real(8):: p1x(3), p2p1(3)
!      real(8):: a(3)
      real(8):: temp 
      integer i

      p1x(:) = p1(:) - x(:)
      p2p1(:)= p2(:) - p1(:)
      call epvectorProd(p1x, p2p1, n)
      temp =  dot_product( n(:), n(:) )
      if( temp == 0. ) then
        icon = 1
!        where( p2p1(:) /= 0. ) 
!           a(:) = -p1x(:)/p2p1(:)
!        endwhere
        do i = 1, 3
           if( p2p1(i) /= 0 )  then
              alpha = -p1x(i)/p2p1(i)
              if(alpha < 0.  .or. alpha > 1.) then
                 icon = 2
              endif
              return
           endif
        enddo
        icon = -1
        alpha = 0.
      else
         n(:) = n(:)/sqrt(temp)
         k = dot_product(n(:), p1(:)) 
         icon = 0 
      endif
    end subroutine ep3p2plane0
