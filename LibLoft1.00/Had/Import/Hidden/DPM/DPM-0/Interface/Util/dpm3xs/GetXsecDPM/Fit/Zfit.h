      integer maxleft, maxright
      integer maxbin
      parameter (maxbin=2000)
      integer nparam, nregion
!    
      parameter (nparam = 5)
      integer minout, minsave
      parameter (minout=56, minsave=7)
      real*8 oparam(nparam)	

      real*8 drx1, drx2  ! drawing region
      real*8 param(nparam)
      real*8 x(maxbin), y(maxbin)
      real*8 maxval, maxdep
      integer maxpos 
      integer npoint, from, to, op
      real*8 chisq, x1h
      common /Zfitc/  x, y,  drx1, drx2,
     * maxval, maxdep,
     * param, 
     * oparam,   chisq, x1h,  npoint,  maxpos, 
     * from, to, op, maxleft, maxright 
