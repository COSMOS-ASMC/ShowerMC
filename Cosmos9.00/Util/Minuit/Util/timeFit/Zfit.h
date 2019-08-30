      integer maxbin
      parameter (maxbin=2000)
      integer nparam, nregion
c    
      parameter (nparam = 3, nregion=2)
      integer minout, minsave
c      parameter (minout=56, minsave=57)
      parameter (minout=56, minsave=7)
      real*8 oparam(nparam)	
      real*8 x1(nregion), x2(nregion)      ! fitting region

      real*8 drx1(nregion), drx2(nregion)  ! drawing region
      real*8 param(nparam, nregion) 
	   real*8 low(nparam, 2)  ! 2 is (g,e,m) and hadron
	   real*8 up(nparam, 2)   ! //
      real*8 x(maxbin), y(maxbin)
      integer npoint
      real*8 chisq, x1h
      common /Zfitc/  x, y, x1, x2, drx1, drx2, param, low, up, 
     *   oparam,   chisq, x1h,  npoint
