!             test kpolintp2Stest kpolintp2S
!  single precision version of kpolintp2
!      implicit none
!      integer m, n, adj
!      parameter (m = 20, n=20, adj=m+5)
!      real*4 xa(m),  ya(n), za(adj, n), x, y, z, error
!      integer i, j, ma, na
!      ma = 4
!      na = 4
!      do i =1, m
!         xa(i) = i/5.
!      enddo
!
!      do j=1, n
!         ya(j) = j/3.
!      enddo
!      do i = 1, m
!         do j = 1,  n
!            za(i, j) = exp(xa(i))* sin(ya(j))
!         enddo
!      enddo
!
!      do i=1, m
!         x= xa(i)+.23
!         do  j=1, n
!            y = ya(j) +.2
!            call kpolintp2S(xa, 1, 1/5., ya, 1, 1/3.,
!     *       za, adj, m, n, ma, na, x, y, z, error)
!c            write(*, *) x, y, exp(x)*sin(y), z, error
!            write(*, *) x, y,
!     *       abs(z-exp(x)*sin(y))/exp(x)*sin(y)
!         enddo
!         write(*,*)
!      enddo
!      end

!       two dimensional interpolation, using kpolintpS
!
      subroutine kpolintp2S(xa, xstep, dx, ya, ystep, dy,
     *  za, adj, m, n,  ma, na,  x, y, z, error)
!     ********************************************************
!      real*4 xa(xstep, m). input.
!      real*4 ya(ystep, n). input.
!      real*4 dx.           input. if non zero, xa is assumed to be
!                                  equistepd with dx.
!      real*4 dy.           input. the same as dx for y.
!      integer adj.         input. adjustable dim. fo z >= m
!      real*4 za(adj, n). input.  function values at (xa, ya)
!      real*4 x.            input.  
!      real*4 y.            input.  
!      real*4 z.            output. estimated func. value at (x, y).
!
!     integer ma.          input.  a max of ma points around x is used
!     integer na.          input.  a max of na points around y is used.
!                    ma, na <= maxp(=10)
!
      implicit none
      integer xstep, ystep,  n, m, na, ma, adj
      real*4 xa(xstep, m), ya(ystep, n), za(adj, n)
      real*4 x, y, z, error, dx, dy

      integer lx, ly, kx, ky, i, maxp
      parameter (maxp = 10)  ! max # of points. 
      real*4  ztmp(maxp)
!
      if(ma .gt. maxp .or. na .gt. maxp) then
         write(*, *) "kpolintp2S: too large number of points"
         stop
      endif
      if(dx .eq. 0.) then
         call kwhereis(x, m, xa, xstep, lx)
      else
         lx = (x - xa(1,1))/dx+1
      endif

      if(dy .eq. 0.) then
         call kwhereis(y, n, ya, ystep, ly)
      else
         ly = (y-ya(1,1))/dy +1
      endif
      kx = min(max(lx - (ma -1)/2,1), m+1-ma)   ! max of ma points from kx
      ky = min(max(ly - (na -1)/2,1), n+1-na)   ! max of na points from ky

!         make ztmp: scan along x direction
      do i = kx, min(m, kx + ma-1)
         if(dy .eq. 0.) then
            call kpolintpS(ya(1, ky), ystep, za(i, ky), adj, 
     *      min(na, n- ky +1), y, ztmp(i-kx+1), error)
         else
            call kpolintpSeqs( (ya(1, 1)+(ky-1)*dy), dy, za(i, ky), adj, 
     *      min(na, n- ky +1), y, ztmp(i-kx+1), error)
         endif
      enddo
      if(dx .eq. 0.) then
         call kpolintpS(xa(1, kx), xstep, ztmp,  1, 
     *      min(ma, m-kx +1),  x,  z, error)
      else
         call kpolintpSeqs( (xa(1, 1)+(kx-1)*dx),  dx, ztmp,  1, 
     *      min(ma, m-kx +1),  x,  z, error)
      endif         
      end

