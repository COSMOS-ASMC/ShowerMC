      integer n, npitbl
      real*8  nx
      parameter(n=101, npitbl= 46,  nx = n-1)
      real*8 uconst, dx
      real*8 ndndx(n), intendndx(n)
      real*8 ndndxn(n), intendndxn(n)
      real*8 intendndx2(n)
      real*8 ppsx(n), pipsx(n), xval(npitbl) 
      real*8 ppsxn(n)
      common /cinippxc/intendndx, ndndx, uconst, dx, ppsx,
     *        intendndx2,  pipsx, xval,
     *   intendndxn, ndndxn, ppsxn




