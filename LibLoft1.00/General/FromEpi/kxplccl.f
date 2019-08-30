      subroutine kxplccl(x0, y0, wxi, wyi, r,  el1, el2, icon)
      implicit none
!         length to the xing point of a directeced line
!         with a circ l
      real(8),intent(in)::x0, y0  !  the line starts from here
      real(8),intent(in)::wxi, wyi  !  direction cos of  the line
                                 !  may not be normalized
                                 ! such as (sinT *cosF, sinT*sinF)  
      real(8),intent(in)::r       !  circle radius (cente is  at 0,0)
      
      real(8),intent(out):: el1, el2  !  if x-point exist, the length to the 
                                 !  nearer one and further one (>0)
                                 ! if icon =1, el1< el2
     
      integer,intent(out):: icon  !  0 if x-sing points xists; x0,y0 in ccle
                                  !  if x0,y0 is on the ccl, it is not xing point
                                  !  1 two x-sing pints;       x0,y0 outside
                                  !  -1  no xing

      real(8):: a, b, c, d
      real(8):: wx, wy, n
      
      n = wxi**2 + wyi**2   !( sin cosf , sin sinf)  n=sin^2
      if(n == 0.) then
         icon = -1
         return
      endif
      n = sqrt(n)           ! n = sin
      wx = wxi/n
      wy = wyi/n

      a = wx**2 + wy**2
      b = wx*x0 + wy*y0
      c = x0**2 + y0**2 -r**2
!        a l^2 + 2 bl + c =  0
      d = b**2 - a*c
      if(d <= 0.) then
         icon = -1
         return
      endif
      d = sqrt(d)
      el1 = (-b-d)/a/n
      el2 = (-b+d)/a/n   ! el2 > el1 since d>0
      if( el1 > 0. ) then
         icon = 1
      elseif(el2 > 0.) then
         el1 = el2
         icon = 0.
      else
         icon = -1
      endif

      end

         

      
         
      
