!
!           Takahashi &  Mori  double exponential formulas
!     

       subroutine ktminte(func, a,b, epsa, epsr, nmin, nmax,
     *                    ans, error, n,   icon)
       implicit none
       real*8 func    ! input. integrand
       external func
       real*8 a       ! input.  lower limit of integration.
       real*8 b       ! input.  upper limit 
       real*8 epsa    ! input.  absolute error you may permit
       real*8 epsr    ! input.  relative error //
       integer nmin   ! input.  min  of the function evaluation points
       integer nmax   ! input.  max  of  //
       real*8  ans    ! output. approximate integral value
       real*8  error  ! output. estimated error
       integer n      ! output  number of function evaluations made
       integer icon   ! output. condition code.
                      !        0  o.k 
                      !        1  ans may have some error > your input.
                      !        2  some input is n.g
!             Be carefule if the integrand has sharp peaks
!             between (a, b). If the peaks are at a or b, 
!             no problem.
       integer nwk
       parameter (nwk = 350)
       real*8 x(nwk), w(nwk)
       real*8 xa, xb
       data itab,ims/0,5/

!     -----------------------------------------------------------------
!     set the constants and tabulate the sample points and weight
!     functions  once for all when the subroutine is called for the
!     first time.
!     -----------------------------------------------------------------
      if(itab .ne. 0) go to 300
      in=2**ims*10
      maxm=in*2+1
      em4=1.d-4
      epsi=dmach(epsi)
      azero=1.0d-30*epsi
      w0=1.5d0
      fm=dfmin(fm)
      p=-dlog(fm)/w0
      tmax=dlog(p)
      h=tmax/dfloat(in)
      do 1010 i=1,in
      t=h*dfloat(i)
      s=w0*dsinh(t)
      c=dcosh(s)
      x(i)=dexp(-s)/c
      if(x(i) .le. 0.0d0) go to 1015
      w(i)=dcosh(t)/c/c*w0
 1010 continue
      go to 1016
 1015 in=i-1
 1016 continue
      itab=1
  300 continue
!     -----------------------------------------------------------------
!     check parameter errors and initialize variables.
!     -----------------------------------------------------------------
      icon=30000
      if(epsa .lt. 0.00d0 .or. (epsr .lt. 0.00d0) .or. (nmin .gt. nmax)
     *     .or.(nmin .lt. 0)) go to 8000
      n=0
      err=epsa
      s=0.0d0
      icon=0
      aa=dmin1(a,b)
      bb=dmax1(a,b)
      rr=(b-a)*0.5d0
      r=dabs(rr)
!     -----------------------------------------------------------------
!     determine the finite interval for applying trapezoidal rule and
!     compute rough estimate of the integral.
!     -----------------------------------------------------------------
      if(r .lt. azero) go to 8000
      xa = (aa+bb)*0.5d0
      xb = r
      v1 = func(xa)*w0
      n=1
      fmax=dabs(v1)
      is=2**ims
      inc=is*10
      im=is
      n1=0
      n2=0
      k1=0
      k2=0
      ea=epsa
      er=epsr
      ep=epsi
      h=tmax/10.0d0
      do 5 i=is,in,im
      xi=x(i)
      if(n1 .gt. 0) go to 2
      xb = -r*xi
      xa = aa-xb
      v = func(xa)*w(i)
      n=n+1
      v1=v1+v
      w1=dabs(v)
      fmax=dmax1(fmax,w1)
      v2=dmin1(dmax1(ea,er*fmax,ep*fmax),em4*fmax)
      if(w1 .le. v2)       go to 1
      k1=0
      go to 3
    1 k1=k1+1
      if(k1 .ge. 2) n1=i-im
      go to 3
    2 if(n2 .gt. 0) go to 9
    3 if(n2 .gt.0) go to 5
      xb = r*xi
      xa = bb- xb
      v = func(xa)*w(i)
      n=n+1
      v1=v1+v
      w2=dabs(v)
      fmax=dmax1(fmax,w2)
      v2=dmin1(dmax1(ea,er*fmax,ep*fmax),em4*fmax)
      if(w2 .le. v2) go to 4
      k2=0
      go to 5
    4 k2=k2+1
      if(k2 .ge. 2) n2=i-im
    5 continue
      if(n1.eq.0) n1=in
      if(n2.eq.0) n2=in
      if(k1 .ne. 0) go to 6
      icon=icon+1000
      ea=w1
      er=ea
      if(dabs(v1) .gt. azero) er=er/dabs(v1)
    6 if(k2 .ne. 0) go to 9
      icon=icon+2000
      ea=dmax1(ea,w2)
      er=ea
      if(dabs(v1) .gt. azero) er=er/dabs(v1)
    9 continue
      fact=1.0d0
      v1=h*v1*rr
!     -----------------------------------------------------------------
!     set the values of test criteria.
!     -----------------------------------------------------------------
      ep=ea
      if(dabs(v1) .gt. azero) ep=ea/dabs(v1)
      ep=dmax1(ep,er)
      if(ep .lt. 1.d-4) go to 101
      fact=1.d0
      go to 200
  101 if(ep .lt. 1.d-5) go to 102
      fact=0.9d0
      go to 200
  102 if(ep .lt. 1.d-10) go to 103
      fact=0.8d0
      go to 200
  103 if(ep .le. 0.0d0) go to 104
      fact=0.75d0
      go to 200
  104 fact=0.75d0
  200 continue
      facti=1.d0/fact
      ep=epsi**0.75d0
      ea=ea**fact
      er=er**fact
!     -----------------------------------------------------------------
!     compute the integral by trapezoidal rule.
!     -----------------------------------------------------------------
      k1=0
      er1=v1
      k2=max0(n1,n2)
      do 20 j=1,ims
      v2=0.0d0
      im=is
      is=is/2
      do 30 i=is,k2,im
      xi=x(i)
      if( i .gt. n1) go to 17
      xb = -r*xi
      xa = aa-xb
      v =func(xa)*w(i)
      n=n+1
      v2=v2+v
   17 if(i .gt. n2) go to 30
      xb =  r*xi
      xa =  bb- xb
      v = func(xa)*w(i)
      n=n+1
      v2=v2+v
   30 continue
      v2=0.5d0*(v1+ rr*h*v2)
!     -----------------------------------------------------------------
!     test the convergency.
!     -----------------------------------------------------------------
      err=dabs(v1-v2)
      if(n .lt. nmin) go to 32
      if(err .le. dmax1(ea,er*dabs(v2))) go to 100
      if(err .lt. er1) go to 34
      k1=k1+1
      if(k1 .ge. 2 .and. j .ge. 3) go to 2024
      go to 36
   34 k1=0
   36 continue
      if(err .lt. ep*dabs(v2)) go to 2022
   32 continue
      if(n+n-1 .gt. nmax) go to 2000
      h=h*0.5d0
      v1=v2
      er1=err
   20 continue
!     -----------------------------------------------------------------
!     set the resultants.
!     -----------------------------------------------------------------
      if(n+n-1 .le. maxm) icon=icon+5000
 2000 icon=icon+20000
      go to 120
 2022 icon=icon+10000
      err=(ep*dabs(v2))**facti
      go to 121
 2024 icon=icon+10000
      err=er1
      go to 120
  100 if(icon .ne. 0) icon=icon+10000
      if(icon .eq.0 .and.(epsa+epsr) .eq. 0) icon=10000
  120 if(err .gt. azero) err=err**facti
  121 s=v2
      err=dmax1(err,epsi*dabs(v2))
!     -----------------------------------------------------------------
!     exit.
!     -----------------------------------------------------------------
 8000 call mgssl(icon,mcode)
      return
      end
