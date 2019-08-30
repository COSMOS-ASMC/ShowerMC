      subroutine kseeinbox(x0,x1, y0,y1, z0, z1, x,y,z, icon)
      implicit none
!          see if x,y,z is inside a box
      real(8),intent(in)::x0, x1  ! box edge in  x direction . from x0 to x1
      real(8),intent(in)::y0, y1   ! box edge in  y direction . from y0 to y1
      real(8),intent(in)::z0, z1   ! box edge in  z direction . from z0 to z1
      
      real(8),intent(in)::x, y, z  ! given point
      integer,intent(out)::icon   ! 0.  point is inside (including on the boundary)
                                  ! 1.  outside

      if( x < x0  .or. x > x1) then
         icon = 1   
      elseif( y < y0  .or. y > y1) then
         icon = 1
      elseif( z < z0 .or.  z > z1 ) then 
         icon = 1
      else
         icon = 0
      endif
      end subroutine kseeinbox
