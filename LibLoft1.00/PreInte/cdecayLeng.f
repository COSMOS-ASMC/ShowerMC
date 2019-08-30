!     ******************************************************************
!     *                                                                *
!     * cdecayLeng: samples decay length of a given particle
!     *                                                                *
!        samples decay length of a particle.
!        
!        For a charged particle, we consider the change of the life
!     time due to the energy change by ionization loss.-->
!     This is now obso.

!     subroutine  cdecayLeng(aTrack, length)
      subroutine  cdecayLeng(pj, length)
      implicit none
#include  "Zglobalc.h"
#include  "Zptcl.h"      
! #include  "Ztrack.h"
! #include  "Ztrackv.h"

      

!     type(track)::aTrack ! input%a  track of a decaying particle
      type(ptcl),intent(inout):: pj     ! a  ptclt 
      real(8),intent(out):: length ! sampled decay length in m

      real*8 g, u,  gbeta, dedt, dedtf, dedl, ctau, rho
      real*8 a, x, p, gmin, ctaumax, pmin, cvh2den
      type(fmom)::gb

      data gmin/200.d0/, ctaumax/50.d0/, pmin/10.d0/

      if(pj%fm%p(4) .lt. pj%mass) then
         pj%fm%p(4) = pj%mass
      endif

      call cgetctau(pj, ctau)
      if(ctau .eq. Infty) then
         length = Infty
      else
         g = pj%fm%p(4)/pj%mass
         if(g .gt. gmin) then
            call rndc(u)
            length = - ctau* g *log(u)
         else
            gbeta =  sqrt(g**2-1.d0)
!            call cgetlf(aTrack.p, gb)
!            gbeta = sqrt(gb.p(1)**2 + gb.p(2)**2 + gb.p(3)**2)
!!            if(aTrack%p%charge .eq. 0 .or. ctau .lt. ctaumax) then
            call rndc(u)
            length = - ctau * gbeta * log(u)
!!            elseif(FromEpics) then
!                  no need to consider energy loss.
!                  this may not be good for rock of large length
!                 In that case we must consider decay length
!                 in Epics.  (don't call cdecayLeng and manage
!                  decay in Epics).
!!               call rndc(u)
!!               length = - ctau * gbeta * log(u)
!!            else
!!               call rndc(u)
!!               length = - ctau * gbeta * log(u)
!                   eloss consideration somethat wrong.
!               rho = cvh2den(aTrack.pos.height)
!               call cdedxInAir(aTrack.p, rho, dedt, dedtf) ! dedt in GeV/(kg/m^2)
!               dedl= dedt*rho   !    GeV/m
!               a = dedl/aTrack.p.fm.p(4) ! 1/m
!               p =1.0d0/( a * ctau * g )
!               if( p .gt. pmin) then
!                  call rndc(u)
!                  length = - ctau* gbeta *log(u)
!               else   
!                  call cdecayWEL(p, g, x)
!                  length =(1.0d0 -  x/g)/a
!                  if(length .lt. 0.) then
!                       may happpen when a <<< 1.
!                     length = 0.d0
!                  endif
!               endif
!!            endif
         endif
      endif
      end
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
      subroutine cdecayWEL(pin, g, x)  ! obso
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
!           integral of 1/sqrt(2(x-1)) from 1 to xm
       int1 = sqrt(2.0*(xm-1.0))
!           integral of (2x)**p/x from (1 to xm)
       int2 = 2.0**p/p *(xm**p-1.0)
!              xm**p = 1.0
       if(int1 .eq. 0. .and. int2 .eq. 0.) then
          x = 1.0
          return  ! **********
       endif

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
          if(x .eq. 1.0) goto 10
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
