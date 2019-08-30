      integer n, maxp
!      parameter ( maxp=25000, n =2000 )
      parameter ( maxp=30000, n =500 )
      real*4  x(maxp, n), y(maxp, n), z(maxp,n)
      real*4  pixel0, pixel, mulalpha, rmax2, pixelinc
      real*4  aax, aay, aaz, bbx, bby, bbz, ccx, ccy, ccz
      integer*2  codex(maxp, n), chgx(maxp, n)
      integer idx(maxp), si(maxp)
      integer jmin, jmax, deltaN, maxppt, offset
      character*100 dir 
      common /Zprivate/  x, y, z, pixel0, pixel, mulalpha, rmax2,
     *    pixelinc,  aax, aay, aaz, bbx, bby, bbz, ccx, ccy, ccz,
     *    idx, si,  jmin, jmax, deltaN,  maxppt, codex, chgx, offset
      common /Zprivate2/ dir

      
