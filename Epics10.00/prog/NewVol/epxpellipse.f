      subroutine epxpellipse(x0, y0, wx, wy,  a, b, ox, oy,
     *     length, icon)
!          crossing point of a line with an ellipse.
      implicit none
      real*8 x0, y0  ! input given line passes this point
      real*8 wx, wy  ! input given line has this direction cos.
      real*8 a, b        ! input ellipse lies on x-y with 'a' radius on x
                         !       and 'b' radius on y.
      real*8 ox, oy      ! input  center of ellips

      real*8 length      ! output. if cross, the length to the ellipse >= 0
      integer icon       ! output. =0, cross. point is inside  ellipse
                         !          1  //     // outside ellipse
                         !         -1  no cross.

      real*8 xx, yy
      real*8 aa, bb, cc, dd
      real*8 leng1, leng2

      xx = x0 - ox
      yy = y0 - oy

      
      aa = (wx/a)**2 + (wy/b)**2
      if(aa .eq. 0.) then
         icon = -1
         return !  *******
      endif
      bb = ( xx*wx/a**2 + yy*wy/b**2)
      cc = (xx/a)**2 + (yy/b)**2 - 1.d0
      
      dd = bb**2 - aa*cc

      if(dd .lt. 0.) then
         icon = -1
         return !  ********
      endif
      dd = sqrt(dd)
      leng1 =( -bb - dd ) / aa
      leng2 =( -bb + dd ) / aa 

      if(leng1 .lt. 0. .and. leng2 .ge. 0.) then
         length = leng2
      elseif(leng1 .ge. 0.) then
         length = leng1
      else
         icon = -1
         return                 ! ********
      endif
      if( cc .gt. 0.) then
         icon = 1
      else
         icon = 0
      endif
      end
