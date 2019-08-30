c     ****************************************************************
c     *                                                              *
c     * mfgmf0: initialize routine for computing geomagnetic field   *
c     *         effect (give gmf constants and 1ry direction)        *
c     * mfgeom: compute deflection and deflection angle of the       *
c     *         particle at given path travelled                     *
c     *                                                              *
c     ************************ tested 83.11.17 ***********************
c
c  /usage/   call mfgmf0(bh, bv, w01,  w02, w03, bt)
c                 this must be called when the 1ry direction is fixed
c                 this must be called when starting the ptcl following
c                 in mf2, may be called at mfcas0 with w1=w1a and w2=w2a
c            call mfgeom(e, iz, dst, dx, dy, tx, ty)
c
c coordinate system:
c
c          let y* be geomagnetic north
c              x*                west
c              z* vertical axis going downwards
c
c      then  gmf is // to y*-z* plane
c
c          let the unit vector of 1ry paticle be z#=(w01, w02, w03) in
c          * system.  ( # is for vector)
c          form x axis by  b# x z#  where b# is the vector of gmf so
c          that  x# = b# x z# / ! b# x z# ! (if b#//z#, x should be
c          x* ).    y axis is on b#-z# plane, i.e y# = z# x x#.
c          the angle between b# and z# is alfa.   gmf effect is
c          computed easily in the frame where x axis is rotated
c          so that z axis coincide with b#.    then the deflection etc
c          is transformed to (x,y,z) system where coordinate of ptcls
c          are measured.
c
c
c  *** mfgmf0 ***
c
c -- input --
c     bh:  horizontal component of gmf             (gauss)
c     bv:  vertical   //               (may be < 0)  //
c    w01:  1ry ptcle direction cosine to x* axis
c    w02:  1ry ptcle direction cosine to y* axis
c    w03:  1ry ptcle direction cosine to z* axis
c
c -- output --
c     bt:  b*sin(alfa) where b is sqrt(bv**2 + bh**2)
c
c
c
c      if e-m cascade is treated 1 dimensionally during ptcl following
c      and 3 dim. sampling is done only at observation plane, the
c      treatmnet of gmf here is not perfect because ptcl direction is
c      afected by scattering in the path.
c
c      think 2 segments of path. at 1st segment, ptcl dirction is given
c      by w1,w2,w3. using this value, gmf deflection is computed at the
c      end of 1st segment.  defelction at the end of 2nd segment must be
c      computed using new w1,w2,w3 which must be corrected by deflection
c      angle at the end of 1st segment.  the defelction angle by gmf is
c      correctly included in this routine, but scattering angle cannot
c      be included because it is not sampled until observation point.
c
c      that is,  both effects are so small that they are treated
c      independently.  the net effect is obtained by adding the both
c      effects.
c
c
c
c  *** mfgmf  ***
c
c -- input --
c     e: energy of ptcl in tev
c    iz: charge of ptcl
c   dst: distance travelled in cm
c
c -- output --
c    dx: deflection of x in (x,y,z) system  (change as compared to
c                                            the streight movement)
c    dy:  //           y
c    tx: deflection angle (x component)
c    ty:  //               y
c
c
c
c
c
      subroutine mfgmf0(bh, bv, w01, w02, w03, bt)
c
      dimension tm(3,3)
      data small /1.e-3/
      save b, tm, cosa, sina
      real*8 dst
c
c          mag. of gmf
      b=sqrt(bh**2 + bv**2)
c          cos(alfa)
      cosa=(bh * w02  +  bv * w03)/b
      cosa2=cosa**2
      sina=sqrt(1. - cosa2)
      bt=b * sina
      sb2= w01**2 + (bh/b-w02)**2 + (bv/b-w03)**2
      if(sb2 .gt. small) then
          cota=cosa/sina
          tm(1,1)= ( bh*w03 - bv*w02 ) /bt
          tm(1,2)= -cota*w01
          tm(1,3)= w01
          tm(2,1)= bv*w01/bt
          tm(2,2)= bh/bt - cota*w02
          tm(2,3)= w02
          tm(3,1)= -bh*w01/bt
          tm(3,2)= bv/bt - cota*w03
          tm(3,3)= w03
      else
          tm(1,1)= 1.
          tm(1,2)= 0.
          tm(1,3)= 0.
          tm(2,1)= 0.
          tm(2,2)=bv/b
          tm(2,3)=bh/b
          tm(3,1)=0.
          tm(3,2)= -bh/b
          tm(3,3)= bv/b
      endif
      return
c
c
c     ************
      entry mfgeom(p, iz, dst, w1i, w2i, w3i, dx, dy, tx, ty)
c     ************
c
      if(iz .eq. 0) then
         dx=0
         dy=0
         tx=0
         ty=0
      else
         z=iz
         wt=z*b*dst*3.e-10/p
         wt2=wt**2/2
         ak=wt*dst/2
c
         dx=  (cosa*w2i - sina*w3i) *ak
         dy= - ak*cosa*w1i
c
c        tx=  (cosa*w2i - sina*w3i) *wt - wt2 *w1i
         tx=  (cosa*w2i - sina*w3i) *wt
c        ty= -wt*cosa*w1i - (cosa2*w2i - sina2*w3i)*wt2
         ty= -wt*cosa*w1i
c        tz= sina*wt*w1i + (csa *w2i - sina2*w3i)*wt2
         tz= sina*wt*w1i
      endif
      return
c
c       convert (x,y,z) in 1ry system into * system
c     ************
      entry mfptos(x, y, z, xs, ys, zs)
c     ************
c
      tmp1= tm(1,1)*x + tm(1,2)*y + tm(1,3)*z
      tmp2= tm(2,1)*x + tm(2,2)*y + tm(2,3)*z
      zs  = tm(3,1)*x + tm(3,2)*y + tm(3,3)*z
      xs  = tmp1
      ys  = tmp2
      return
c
c       convert (x,y,z) in * system into  1ry system
c     ************
      entry mfstop(xs, ys, zs, x, y, z)
c     ************
c
      tmp1= tm(1,1)*xs + tm(2,1)*ys + tm(3,1)*zs
      tmp2= tm(1,2)*xs + tm(2,2)*ys + tm(3,2)*zs
      z  = tm(1,3)*xs + tm(2,3)*ys + tm(3,3)*zs
      x  = tmp1
      y  = tmp2
      end
