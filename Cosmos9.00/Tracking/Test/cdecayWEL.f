!  decay with costant rate of energy loss.
!
!       decay probability function can be expressed as
!   
!   p(x)dx=   dx 1/(x-sqrt(x**2-1))**p /sqrt(x**2-1)
!
!     the range of x is 1 to g=E0/m.
!   p is almost independent of energy and for muons
!     0.7 to 10.  For larger p's, we can use usual 
!    exp probability.  
!
      subroutine cdecayWEL(pin, g, x)
      implicit none
      real*8 pin ! input.  see above. should be 0.1<pin<10.
      real*8 g  ! input.  E0/m.  1<= g 
      real*8 x ! output. sampled x.  decay length 'l' is related to this
               !           by x=g(1-al) where a is dE/dl/E0 (/m).
!
!     Method:  
!        If x > 5, we use p(x)=(2x)**p/x
!           x < 5, we use (2x)**p/x + 1/sqrt(2(x-1)) and rejection
!   To decide which side, we compare  int(x=1 to  5) of p(x)dx and
!                                     int(x=5 to g) of   (2x)**p/x 
!
!  if  g< 5,  we use rejection method only.
!
!    log( int(x=1 to 5)) is approximated by 4-th order polynomial:
!         sum c_i p**i  (i=0 to 4)
!
!
!    c0        .77099
!    c1        1.3470  
!    c2        .12049 
!    c3       -.57001E-02 
!
      real*4 int1, int2, ans, xm, tf, rf, p
      real*8 u
      
      p = pin
      if(p .gt. 10.) then
         call cerrorMsg('p> 10 for cdecayWEL', 0)
      endif
      if(g .gt. 5.) then
         if(p .le. 0.1d0) then
            ans = 0.771
         else
            ans =((-0.57001E-02*p + 0.12049  )*p+1.3470)*p
     *        +0.77099
         endif
         xm= 5.
!              int(1 to 5) of p(x)dx
         int1 = exp(ans)
!              int(5 to g) of p(x)dx~ (2x)**p/xdx
         int2 = 2.0**p/p * (g**p - xm**p)
      else
         int1 = 1.  ! dummy value
         int2 = 0.  ! //  so that int1 > int2
         xm = g
      endif
      call rndc(u)
      if(u  .lt. int1/(int1+int2)) then
!            use rejection.
!           integral of 1/sqrt(2(x-1)) from 1 go xm
       int1 = sqrt(2.0*(xm-1.0))
!           integral of (2x)**p/x from (1 to xm)
       int2 = 2.0**p/p *(xm**p-1.0)
       do while (.true.)
!
          call rndc(u)
          if(u .lt. int1 /(int1+int2)) then
!           use dx/sqrt(2(x-1))
             call rndc(u)
             x = (u*int1)**2/2.0+1.0
          else
!            use dx(2x)**p/x 
             call rndc(u)
             x =( p*int2*u/2.0**p+ 1.)**(1./p)
          endif
          call rndc(u)
          tf = 1./(x-sqrt(x*x-1.0))**p/sqrt(x*x-1.0)
          rf = (2.0*x)**p/x + 1.0/sqrt(2*(x-1.0))
          if(u .lt. tf/rf) goto 10
       enddo
 10    continue
      else
!         use (2x)**p/x dx
         call rndc(u)
         x =( p*int2*u/2.0**p+ xm**p)**(1./p)
      endif
      end



