      subroutine kseeincyl(x0, y0, z0, r, h, x,y,z, icon)
      implicit none
!          see if x,y,z is inside a cylinder
      real(8),intent(in)::x0, y0, z0  ! cylinder bottom center 
      real(8),intent(in)::r  ! cylinder radius
      real(8),intent(in)::h  ! cylinder height
      
      real(8),intent(in)::x, y, z  ! given point
      integer,intent(out)::icon   ! 0.  point is inside (including on the boundary)
                                  ! 1.  outside

      if( z < z0 .or. z > h ) then
         icon = 1
      else
        if( (x-x0)**2 + (y-y0)**2 > r**2 ) then
           icon = 1
        else
           icon = 0
        endif
      endif
      end
