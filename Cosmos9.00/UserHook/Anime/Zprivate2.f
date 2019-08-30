      integer n, maxp
      parameter ( maxp=30000, n = 1000 )
      real*4  x(maxp, n), y(maxp, n), z(maxp,n)
      real*4  pixel, mulalpha
      integer*2  codex(maxp, n), chgx(maxp, n)
      integer idx(maxp), si(maxp)
      integer jmin, jmax
      character*100 dir 
      common /Zprivate/  x, y, z, pixel, mulalpha, idx, si,  
     *   jmin, jmax,   codex, chgx
      common /Zprivate2/ dir

      
