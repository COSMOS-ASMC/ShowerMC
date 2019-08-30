      subroutine kScsplInteg(x, y, n,  coef, nc, a, b, s)
      implicit none
      integer n   ! input. size of x,y
      real(4):: x(n), y(n)  ! input.
      integer nc   ! input. size of coef
      real(4):: coef(nc, 3)  ! input. which is the output from 
                          !   kcsplCoeff.
      real(4):: a    ! input. lower limt of the integral
      real(4):: b    ! input. upper //
      real(4):: s    ! output. obtained integral value.
      

      real(4) xx, x2, x1, xa, xb, sa, sb, sab, h
      integer ind, ia, ib, ifla, iflb, ip, ier, i

      ind = 1
      ia = 1
      ifla = 0
      iflb = 0
      xx = a
      if (a .gt. b) xx = b
      do while( ind .le. 2 )   !  until ind=3
         x1 = xx - x(ia)
         do i=ia,n-1
            ip = i
            x2 = xx - x(i+1)
            
            if (x2 .lt. 0.0) go to 30
            if (i .lt. n-1) x1 = x2
         enddo
         ip = n - 1
         if (x2 .gt. 0.0) ier = 1
 30      continue
         if (x1 .lt. 0.0) then
            ier = 1
            x1 = - x1
            if (ind .eq. 1) ifla = 1
            if (ind .eq. 2) iflb = 1
         endif
         if(ind .ne. 2) then
            ia = ip
            xa = x1
            xx = b
            if (a .gt. b) xx = a
         endif
         ind = ind + 1
      enddo

      ib = ip
      xb = x1
!          integral from a to b.
      sa = y(ia)+
     *  xa*(coef(ia,1)/2.0 + 
     *       xa*(coef(ia,2)/3.0+xa*coef(ia,3)/4.0))

      sa = sa * xa
      if (ifla .eq. 1) sa = - sa
      sab = 0.0

      if (ib-1 .ge. ia) then
         do  i=ia,ib-1
            h = x(i+1) - x(i)
            sab = 
     *      sab+h*(y(i+1)+y(i)-
     *      (coef(i+1,2)+coef(i,2))*h*h/6.0)/2.0
         enddo
      endif

      sb = y(ib)+
     *   xb*(coef(ib,1)/2.0 +
     *        xb*(coef(ib,2)/3.0+xb*coef(ib,3)/4.0))
      sb = sb * xb
      if (iflb .eq. 1) sb = - sb
      s = sb + sab - sa
      if (a .gt. b) s =  - s

      end


