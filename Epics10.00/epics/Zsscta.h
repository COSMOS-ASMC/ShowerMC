!           /Zssctana/
      real*8 pxstp, pystp, e1, x1, y1, z1, wx1, wy1, wz1

      real*8 pzstpx(nstpl), pzstpy(nstpl),  boxz0(maxcmp), boxc(maxcmp)
      integer nct, nev, k1, ic1,  nevmax, lx, ly,
      integer jxstp, jystp,  lusex, lusey

        common /Zsscta/ pxstp(ntstpx), pystp(ntstpy),
     *  pzstpx(nstpl), pzstpy(nstpl),  boxz0(maxcmp), boxc(maxcmp), 
     *  e1, x1, y1, z1, wx1, wy1, wz1,
     *  jxstp(maxcmp), jystp(maxcmp),  lusex(nstpl), lusey(nstpl),
     *  nct,   nev,  k1, ic1, nevmax, lx, ly



!
