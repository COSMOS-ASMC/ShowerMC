      subroutine kxplsq(x0, y0, wxi, wyi, a, b,  el1, el2, icon)
      implicit none
!         length to the xing point of a directeced line
!         with a squre
      real(8),intent(in)::x0, y0  !  the line starts from here
      real(8),intent(in)::wxi, wyi !  direction cos of  the line (may not be normalized)
      real(8),intent(in)::a, b    !  the x and y direction edge length of the square
                                  !  the left bottom corner is at (0,0)
      
      real(8),intent(out):: el1, el2  !  if x-point exist, the length to the 
                                 !  nearer one and further one (>0)
     
      integer,intent(out):: icon  !  0 if x-sing points xists; x0,y0 in square
                                  !  if x0,y0 is on the edge, it is not xing point
                                  !  1 two x-sing pints;       x0,y0 outside
                                  !  -1  no xing
      real(8):: x, y, el, n, wx, wy
      real(8):: ela(2)
      integer:: xpc, in
      
      xpc = 0
 

      n = wxi**2 + wyi**2
      if( n > 0.d0 ) then
!       suppose  true direction cos: Wx Wy
!       and   wxi = sin Wx  wyi=sin Wx
!       wxi**2 + wyi**2 = sin**2 ; so true 
!        Wx = wxi/sin  Wy = wyi/sin
         n = sqrt(n)   ! this is sin
         wx = wxi/n
         wy = wyi/n
      else
         icon = -1
         return
      endif

      if(x0 >= 0.0 .and. x0<= a .and. y0>=0. and. y0<= b) then
         in = 0
      else
         in = 1
      endif

      if(wy /= 0.0) then
      ! see y = 0 line
         el = -y0/wy
         if(el > 0.) then
            x = x0 + el*wx
            if(x >= 0.0 .and. x <= a ) then
               xpc = xpc + 1
               ela(xpc) = el
               if(in == 0 ) goto 100
            endif
         endif
      !  see if y=b line
         el = (b-y0)/wy
         if(el > 0.) then
            x = x0 + el*wx
            if(x >= 0.0 .and. x <= a ) then
               xpc = xpc + 1
               ela(xpc) = el
               if( in == 0 ) goto 100
            endif
         endif
         if(xpc == 2) goto 100
      endif
      if(wx /= 0.) then
      ! see x = 0 line
         el = -x0/wx
         if(el > 0.) then
            y = y0 + el*wy
            if(y >= 0.0 .and. y <= b ) then
               xpc = xpc + 1
               ela(xpc) = el
               if(in == 0 )  goto 100
            endif
         endif
      !  see if x=a line
         el = (a-x0)/wx
         if(el > 0.) then
            y = y0 + el*wy
            if(y >= 0.0 .and. y <= b ) then
               xpc = xpc + 1
               ela(xpc) = el
               if( in == 0 ) goto 100
            endif
         endif
      endif
 100  continue
      if(xpc == 0) then
         icon = -1
      elseif(xpc == 1 ) then
         el1 = ela(1)/n
         icon = 0
      else
         el1 = min(ela(1), ela(2))/n
         el2 = max(ela(1), ela(2))/n
         icon = 1
      endif

      end

