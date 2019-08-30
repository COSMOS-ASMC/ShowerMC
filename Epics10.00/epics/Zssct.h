       real*8 stpw, unitx, unity, gap, 
       integer nx, ny, maxcmp, nstrpx
           parameter(stpw=.5, unitx=8.5, unity=8.5, gap=1., nx=6, ny=6,
     * maxcmp=200,   nstrpx=unitx/stpw+.5, nstrpy=unity/stpw+.5,
     * ntstpx=nstrpx*nx, ntstpy=nstrpy*ny, unitxt=unitx+gap,
     * unityt=unity+gap, nstpl=10)
           common /Zssct/ jstpx(maxcmp), jstpy(maxcmp),
     1     nstpla,
     2     decomp(maxcmp), tmcomp(maxcmp),side(maxcmp),
     3     dexstp(ntstpx,nstpl), deystp(ntstpy,nstpl)
           common /Zssctc/ mat(maxcmp)
           character*8 mat
