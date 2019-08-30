      integer n, maxp
#if defined  G95
      parameter ( maxp = 60000, n = 1300 )
#elif defined G77
      parameter ( maxp=25000, n =2000 )
!     parameter ( maxp=100000, n =500 )
#else
      parameter ( maxp=10000, n= 1000)
#endif
      real*4  x(maxp, n), y(maxp, n), z(maxp,n)
      real*4  aax, aay, aaz, bbx, bby, bbz, ccx, ccy, ccz
      integer*2  codex(maxp, n), chgx(maxp, n)
      integer*2  thinc(n)
      integer*2  maxthin
      real*4  pixel,  rmax2
      integer idx(n), si(maxp)
      integer codesel(9), chgsel, select, ncodes
      integer jmin, jmax,  maxppt, offset, outtype
      logical split
      character*100 dir 
      common /ZprivateA/  x
      common /ZprivateB/  y
      common /ZprivateC/  z
      common /ZprivateE/  codex, chgx
      common /ZprivateD/  pixel,  rmax2, split, outtype,
     *    aax, aay, aaz, bbx, bby, bbz, ccx, ccy, ccz,
     *    idx, si,  jmin, jmax,  maxppt, offset, 
     *    chgsel, select, codesel, ncodes,
     *    thinc,  maxthin
      common /Zprivate2/ dir


      

      
