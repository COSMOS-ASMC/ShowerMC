      subroutine epLightGetBoxSurfN(a, b, c, x, y, z, surfn)
      use modepLightMaxDef
      implicit none
!         get surface number of a given point (x,y,z)
!         which  is assumed to be very close to one of the
!         6 surface of a box 
!
      real(8),intent(in):: a, b, c  !input.  3 edge lengths of a canonical box.
      real(8),intent(in):: x, y, z  !input    given point
      integer,intent(out):: surfn  ! output. surface number.  (1~6)
                     !       0 means (x,y,z) is "far" from the surface

      real*8 o(6) 
      real*8  mind
      integer i

      if( z .eq. 0.) then
         surfn = 1
      elseif( z .eq. c ) then
         surfn = 6
      elseif( x .eq. 0.) then
         surfn = 5
      elseif( x .eq. a ) then
         surfn = 2
      elseif( y .eq. 0. ) then
         surfn = 3
      elseif( y .eq. b ) then
         surfn = 4
      else
         o(1) = z
         o(2) = a-x
         o(3) = y
         o(4) = b-y
         o(5) = x
         o(6) = c-z
         mind = min( abs(o(1)), abs(o(2)), abs(o(3)), 
     *            abs(o(4)), abs(o(5)),abs(o(6)) )
         do i = 1, 6
            if(mind .eq. abs(o(i)) ) exit
         enddo
         surfn = i
         if(abs(o(i)) .gt. EpsLength2 ) then
            write(0,*) ' warning: surface check i=',i, o(i)
            write(0,*) ' a, b,c=', a,b,c
            write(0,*) ' x, y,z=', x,y,z
            write(0,*) ' for box'
            write(0,*) '  '
         endif   
      endif
      end
