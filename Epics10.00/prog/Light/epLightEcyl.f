      subroutine epLightGetEcylSurfN(ra, rb,  h,  x, y, z, surfn)
      use modepLightMaxDef
      implicit none
!         get surface number of a given point (x,y,z)
!         which  is assumed to be very close to one of the
!         3 surfaces of a cyl
!
      real(8),intent(in)::ra, rb, h  ! radius and hight (cm)
      real(8),intent(in):: x, y, z  !  given point
      integer,intent(out)::surfn  ! output. surface number.  (1~3)
                     !   1--> bottom, 2 top, 3 side

      real*8 o(3) 
      real*8  mind
      integer i

      if( z == 0.) then
         surfn = 1
      elseif( z == h ) then
         surfn = 2
      else
         o(1) = z**2
         o(2) = (h-z)**2
         o(3) = (x/ra)**2 + (y/rb)**2 - 1.
         mind = min( abs(o(1)), abs(o(2)), abs(o(3)))
         do i = 1, 3
            if(mind .eq. abs(o(i)) ) exit
         enddo
         surfn = i
         if(abs(o(i)) .gt. EpsLength2**2 ) then
            write(0,*) ' warning: surface check i=',i, o(i)
         endif   
      endif
      end
