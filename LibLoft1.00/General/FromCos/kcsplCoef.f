      subroutine kcsplCoef(x, y, n, coef, nc)
!      compute the coefficients of the cubic spline
      implicit none
      integer n  ! input. size of x, y
      real*8 x(n)  ! input. 
      real*8 y(n)  ! input.
      integer nc   ! input. size of  1st dim. of coef. >= n-1
      real*8 coef(nc, 3)  ! output. obtained coeff.

      character*100 msg

      real*8  ec(4), d(2)


      integer i, ii, ks, ke, ip, i1, i2, k, j
      real*8 a2, a1, x1, x2, x3, h, hy, dy1, dy2, y1, y2, hh, h1, h2
      real*8 pinv
!
      if (n .lt. 2 .or. nc .lt. n-1 ) then
         write(0,*) ' n=', n, ' nc=',nc
         write(0,*) " if (n .lt. 2 .or. nc .lt. n-1 ) wrong"
         call cerrorMsg('input is wrong for kcslpCoef',0)
      endif

      do  i = 1, n-1
         if ( x(i) .eq. x(i+1)) then
            call cerrorMsg(
     *     'x must be ascending: kcsplCoef',1)
            write(msg,*) ' i=',i, ' x(i)=',x(i), ' x(i+1)=',x(i+1)
            call cerrorMsg(msg, 0)
         endif
      enddo


      ii = 2
      ks = 1
      ke = min(4, n)
      ip = 1
      do  i = 1, 2
         i1 = 2 * i - 1
         i2 = 2 * i
         if (i .ne. 1) then
            ks = max(1, n-3)
            ke = n
            ip = n
         endif
         a2 = 0.
         d(i) = 0.
         do k = ks, ke
            if (ip .ne. k) then
               a1 = 1.0d0
               do  j = ks, ke
                  if (j .ne. ip .and. j .ne. k) then
                     x1 = x(ip) - x(j)
                     x2 = x(k) - x(j)
                     a1 = a1 * x1 / x2
                  endif
               enddo
               x3 = x(k) - x(ip)
               d(i) = d(i) + a1 * y(k) / x3
               a2 = a2 - 1.0d0 / x3
            endif
         enddo
         d(i) = d(i) + y(ip) * a2
         if (i .eq. 2) ii = n
         h = x(ii) - x(ii-1)
         ec(i1) = 1.0d0
         hy = y(ii) - y(ii-1)
         ec(i2) = 6.0d0 * (hy / h - d(i)) / h
         if (i .eq. 2) ec(i2) = - ec(i2)
      enddo       
      
      if (n .ne. 2) then
         h1 = x(2) - x(1)
         y1 = y(2) - y(1)
         do i = 2, n-1
            h2 = x(i+1) - x(i)
            y2 = y(i+1) - y(i)
            hh = h1 + h2
            coef(i,1) = h2 / hh
            coef(i,2) = 1.0d0 - coef(i,1)
            coef(i,3) = 6.0d0 * (y2 / h2 - y1 / h1) / hh
            h1 = h2
            y1 = y2
         enddo
      endif

      coef(1,1) = - ec(1) * 0.5d0
      coef(1,2) =   ec(2) * 0.5d0
      if (n .ne. 2) then
         do k = 2, n-1
            pinv = 2.0d0 + coef(k,2) * coef(k-1,1)
            coef(k,1) = - coef(k,1) / pinv
            coef(k,2) =
     *        (coef(k,3) - coef(k,2) * coef(k-1,2)) / pinv
         enddo
      endif
      dy1 = (ec(4) - ec(3) * coef(n-1,2)) /
     *               (2.0d0 + ec(3) * coef(n-1,1))
      do  i = 1, n-1
         k = n - i
         dy2 = coef(k,1) * dy1 + coef(k,2)
         h = x(k+1) - x(k)
         coef(k,3) = (dy1 - dy2) / (6.0d0 * h)
         coef(k,2) = 0.5d0 * dy2
         coef(k,1) = (y(k+1) - y(k)) / h -
     *         (coef(k,2) + coef(k,3) * h) * h
         dy1 = dy2
      enddo
      end
