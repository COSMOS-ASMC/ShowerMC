      subroutine kxplsph(x0, y0, z0, l, m, n, r, el, icon)
      implicit none
      real*8  x0, y0, z0 ! input. the line passes this point
      real*8  l, m, n  !  input.  direc cos.  of  the line
      real*8  r        !  input.  radius of the sphere
      real*8  el       !  output. el>=0 distance to the
                       !          sphere  from  x0,y0,z0
      integer icon    !  output. icon =0.  x-point exists 
                      !                  x0,.. is inside
                      !          icon = 1  x-point exists
                      !                  x0.. is outside
                      !                =-1.  no x-point

      real*8  rsqr, r0l, d
      integer icon1, icon2 
      
      rsqr = x0**2 + y0**2 + z0**2 -r**2
      if(rsqr .le. 0.) then
!          inside
         icon2 = 0
      else
         icon2 = 1
      endif
      r0l = x0*l + y0*m + z0*n
      d = r0l**2 - rsqr
      if(d .ge. 0.) then
         d = sqrt(d)
         el = -r0l - d
         if(el .ge. 0.) then
            icon1 = 0
         else
            el = -r0l + d
            if(el .ge. 0.) then
               icon1 = 0
            else
               icon1 = 1
            endif
         endif
      else
         icon1 = 1
      endif
!
      if(icon2 .eq. 0) then
         icon = 0
      elseif(icon1 .eq. 0) then
         icon = 1
      else
         icon = -1
      endif
      end

