           parameter(stpw=.1, unitx=28., unity=28.,
     * maxcmp=200,   nstrpx=unitx/stpw+.5, nstrpy=unity/stpw+.5,
     * nstpl=10)
           common /Zessct/jstpx(maxcmp), jstpy(maxcmp),
     1     nstpla,
     2     decomp(maxcmp), tmcomp(maxcmp),
     3     dexstp(nstrpx,nstpl), deystp(nstrpy,nstpl)
           common /Zesstc/ mat(maxcmp)
           character*8 mat
