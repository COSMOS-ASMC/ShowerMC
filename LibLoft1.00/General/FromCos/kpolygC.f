!  this is copy of kpolyg.f from Epics.  C is added in the last part of
! eahc name.
!c               to test kpolygC
!      implicit none
!      real*8 x,y, kpolygC
!      integer i, n
!      do n=0, 0
!         do i=1, 2500
!            x =-0.10 + (i-1)* 0.01
!            y = kpolygC(n, x)
!            write(*,*) x, y
!         enddo
!         write(*, *) 
!      enddo
!      end
!
!     *******************************************
!     *                                         *
!     *  kpolygC:  n-th logarithmic derivative    *
!     *          of gamma function.             *
!     *                                         *
!     *******************************************
!
! usage:
!         y= kpolygC(n,x)
!
!    n:  0,1,2... specifies n-th derivative
!    x:  argument. may be negative.
!        kpolygC(0,x) = psi(x) = d ln(Gamma(x))/dx
!         d  (G'/G ) /dx
!        kpolygC(1, x) =d psi(x)/dx etc
!      Note  psi(x) !=   dln(Gamma(x+1))/dx
!
!
       real*8 function kpolygC(n,x)
       implicit none
       real*8 x
       integer n
!
!
!
!
!     Computes n-th derivative of psi function with 7-8 significant dig
!     ts. psi(x) is logarithmic derivative of gamma function. psi(1)=-g
!     where g=0.57721566---Euler const. in Rossi notation, psi(s+1)=
!     dlog(s!)/ds=polyg(0,s+1.), where s! means factorial of s.  x must
!     not be 0,-1,-2,---, and for practical reason must be > -20, n<
!     nmax, (of course must be >= 0). if so message printed, and polyg=
!     1.0e35.
!
!
      integer,parameter::nmax=50 
      real*8 f(nmax)
      real*8 b(6)/8.333333e-2, -1.388889e-3, 3.3068783e-5,
     1         -8.267196e-7,  2.087676e-8, -52.841901e-11/
      logical first/.true./
      integer i, nt, j, m
      real*8 z, s, znt, z2, sn, sum, zn
      save first 

      if(first) then
!         compute factorial up to nmax
         first=.false.
         f(1)=1.
         do i=2,nmax
            f(i)=float(i-1)*f(i-1)
         enddo   
      endif   
!
      z=x
      if(z .gt. 0. .or. z - aint(z) .ne. 0. .and. 
     *           (n .ge. 0  .and.  n .lt. nmax)  ) then
         s=0.
         nt=n+1
         do while (z .lt. 4.)
            znt=1.
            do i=1,nt
               znt=znt*z
            enddo   
            s=s+1./znt
            z=z+1.
         enddo
         s=f(nt)*s
         z2=z*z
         zn=1.
         do while (nt .gt. 1)
            nt=nt-1
            zn=zn*z
         enddo   
         sn=-1.
         if(mod(n,2).ne.0) sn=1.
         sum=0.
         do i=1,6
            j=7-i
            m=j+j+n
            sum=(sum+f(m)*b(j))/z2
         enddo 
         sum=((sum+0.5*f(n+1)/z)/zn+s)*sn
         if(n .eq. 0) then
            kpolygC=log(z)+sum
         else   
            kpolygC=sum+f(n)/zn*sn
         endif   
      else   
         write(*,
     *   '( ''***error input to kpolygC: (z,n)='',e18.8,i10)') z,n
         kpolygC=1.e36
      endif   
      end
