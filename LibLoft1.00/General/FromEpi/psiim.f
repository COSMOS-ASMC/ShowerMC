!     ****************************************************************
!     *                                                              *
!     *  psiim:  sum of 1/( (n+x)**2 + y**2 ) from n=0 or 1 to inf   *
!     *          related to the imaginary part of psi function       *
!     *                                                              *
!     ****************************************************************
!
!    /usage/
!              f = psiim(x, y, n1, eps)
!
!     x:    input argument
!     y:    //
!    n1:    0 or 1.  n=n1 to infinity is summed.
!   eps:    specifies relative error for f
!
! ***  note ***
!       for n1=0, (x,y) must not be (0,0)
!
!       psiim is related to the inaginary part of psi-function psi(z=x+
!     *y), that is, im( psi(z=x+i*y) ) is y*psiim (n1=0;  (x,y)^=(0,0)
!
!
!
!
!
      function psiim(x,y,n1,eps)
!
!
!
!     c0,c1,..c4 are zeta function at n=2,3,..,6...
!
      data c0,c1,c2,c3,c4/1.6449340668,1.2020569032,1.0823232337,1.0369
     177551,1.0173430620/
!
      as=x*x
      bs=y*y
      asmbs=as-bs
      asbs4=as*bs*4.
      s1=x+x
      s2=3.*as-bs
      s3=(s1+s1)*asmbs
      s4=(5.*as-bs)*asmbs-asbs4
      se1=-s1*(3.*asmbs*asmbs-asbs4)
      se2=-(as+bs)*s4
      s=c0-s1*c1+s2*c2-s3*c3+s4*c4
      fn=1.
      tmpo=1.
    5 continue
      fns=fn*fn
      tmp=fns*fn
      tmp=tmp*tmp
      tmp=(se1*fn+se2)/tmp/((fn+x)**2+bs)
      s=s+tmp
      if(abs(tmp*fn/s)+abs(tmpo*fn/s) .lt. eps) goto 10
      tmpo=tmp
      fn=fn+1.
      go to 5
!
   10 continue
      if(n1 .eq. 0) s=s+1./(as+bs)
      psiim=s
      return
      end
