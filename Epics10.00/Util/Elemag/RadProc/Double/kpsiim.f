!     ****************************************************************
!     *                                                              *
!     *  kpsiim:  sum of 1/( (n+x)**2 + y**2 ) from n=0 or 1 to inf  *
!     *          related to the imaginary part of Psi function       *
!     *                                                              *
!     ****************************************************************
!
!    /usage/
!              f = kpsiim(x, y, n1, eps)
!
!     x:    input argument
!     y:    //
!    n1:    0 or 1.  n=n1 to infinity is summed.
!   eps:    specifies relative error for f
!
! ***  note ***
!       for n1=0, (x,y) must not be (0,0)
!
!       kpsiim is related to the imaginary part of psi-function psi(z=x+
!     iy), that is, im( psi(z=x+iy) ) is y*kpsiim (n1=0;  (x,y)^=(0,0))
!
!
!
      real*8 function kpsiim( x, y, n1, eps )
      implicit none
      real*8 x
      real*8 y
      integer n1
      real*8 eps
!
!
      real*8 as, bs, asmbs, asbs4, s1, s2
      real*8 s3, s4, se1, se2, s, fn, tmpo, fns, tmp
!
!     c0,c1,..c4 are zeta function at n=2,3,..,6...
!
      real*8 c0, c1, c2, c3, c4
      data c0,c1,c2,c3,c4/1.6449340668d0, 1.2020569032d0,
     *    1.0823232337d0, 1.0369177551d0, 1.0173430620d0/
!
      as = x*x
      bs = y*y
      asmbs = as-bs
      asbs4 = as*bs*4.
      s1 = x+x
      s2 = 3.*as-bs
      s3 = (s1+s1)*asmbs
      s4 = (5.*as-bs)*asmbs-asbs4
      se1 = -s1*(3.*asmbs*asmbs-asbs4)
      se2 = -(as+bs)*s4
      s = c0-s1*c1+s2*c2-s3*c3+s4*c4
      fn = 1.
      tmpo = 1.

      do while(.true.)
         fns=fn*fn
         tmp=fns*fn
         tmp=tmp*tmp
         tmp=(se1*fn+se2)/tmp/((fn+x)**2+bs)
         s=s+tmp
         if(abs(tmp*fn/s)+abs(tmpo*fn/s) .lt. eps) goto 10
         tmpo=tmp
         fn=fn+1.
      enddo
!
   10 continue

      if(n1 .eq. 0) s = s+1./(as+bs)
      kpsiim=s
      end
