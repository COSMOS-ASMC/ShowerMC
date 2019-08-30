!             test kpolintp2test kpolintp2
!
!      implicit none
!      integer m, n, adj
!      parameter (m = 20, n=20, adj=m+5)
!      real*8 xa(m),  ya(n), za(adj, n), x, y, z, error
!      integer i, j, ma, na
!      ma = 4
!      na = 4
!      do i =1, m
!         xa(i) = i/5.d0
!      enddo
!c
!      do j=1, n
!         ya(j) = j/3.d0
!      enddo
!      do i = 1, m
!         do j = 1,  n
!            za(i, j) = exp(xa(i))* sin(ya(j))
!         enddo
!      enddo
!c
!      do i=1, m
!         x= xa(i)+.23
!         do  j=1, n
!            y = ya(j) +.2
!!            call kpolintp2(xa, 1, 1/5.d0, ya, 1, 1/3.d0,
!            call kpolintp2(xa, 1, 0.d0, ya, 1, 0.0,
!     *       za, adj, m, n, ma, na, x, y, z, error)
!            write(*, *) x, y, exp(x)*sin(y), z, error
!!            write(*, *) sngl(x), sngl(y),
!!     *       sngl(abs(z-exp(x)*sin(y))/exp(x)*sin(y))
!         enddo
!         write(*,*)
!      enddo
!      end
!
!       two dimensional interpolation, using kpolintp
!
      subroutine kpolintp2(xa, xstep, dx, ya, ystep, dy,
     *  za, adj, m, n,  ma, na,  x, y, z, error)
!     ********************************************************
!      real*8 xa(xstep, m). input.
!      real*8 ya(ystep, n). input.
!      real*8 dx.           input. if non zero, xa is assumed to be
!                                  equistepd with dx.
!      real*8 dy.           input. the same as dx for y.
!      integer adj.         input. adjustable dim. fo z >= m
!      real*8 za(adj, n). input.  function values at (xa, ya)
!      real*8 x.            input.  
!      real*8 y.            input.  
!      real*8 z.            output. estimated func. value at (x, y).
!
!     integer ma.          input.  a max of ma points around x is used
!     integer na.          input.  a max of na points around y is used.
!                    ma, na <= maxp(=10)
!
      implicit none
      integer xstep, ystep,  n, m, na, ma, adj
      real*8 xa(xstep, m), ya(ystep, n), za(adj, n)
      real*8 x, y, z, error, dx, dy

      integer lx, ly, kx, ky, i, maxp
      parameter (maxp = 10)  ! max # of points. 
      real*8  ztmp(maxp)
      if(ma .gt. maxp .or. na .gt. maxp) then
         write(*, *) "kpolintp2: too large number of point"
         stop
      endif

      if(dx .eq. 0.) then
         call kdwhereis(x, m, xa, xstep, lx)
      else
         lx = (x - xa(1,1))/dx+1
      endif

      if(dy .eq. 0.) then
         call kdwhereis(y, n, ya, ystep, ly)
      else
         ly = (y-ya(1,1))/dy +1
      endif
      kx = min(max(lx - (ma -1)/2,1), m+1-ma)   ! max of ma points from kx
      ky = min(max(ly - (na -1)/2,1), n+1-na)   ! max of na points from ky

!         make ztmp: scan along x direction
      do i = kx, min(m, kx + ma-1)
         if(dy .eq. 0.) then
            call kpolintp(ya(1, ky), ystep, za(i, ky), adj, 
     *      min(na, n- ky +1), y, ztmp(i-kx+1), error)
         else
            call kpolintpeqs( (ya(1, 1)+(ky-1)*dy), dy, za(i, ky), adj, 
     *      min(na, n- ky +1), y, ztmp(i-kx+1), error)
         endif
      enddo
      if(dx .eq. 0.) then
         call kpolintp(xa(1, kx), xstep, ztmp,  1, 
     *      min(ma, m-kx +1),  x,  z, error)
      else
         call kpolintpeqs( (xa(1, 1)+(kx-1)*dx),  dx, ztmp,  1, 
     *      min(ma, m-kx +1),  x,  z, error)
      endif         
      end



      
