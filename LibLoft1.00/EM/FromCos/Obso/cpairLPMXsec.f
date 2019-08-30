!           test cpairLPMXsec
!      implicit none
!      integer j, i
!      real*8  rho, e, xs
!      do  j=1, 11
!          rho=1. * 10.**(-(j-1)/2.)
!          do i=1, 14
!             e=1.e8 * 10.**( (i-1)/3.0)
!             call cpairLPMXsec(e, rho, xs)
!             write(*, *) sngl(e), sngl(xs)
!          enddo
!          write(*, *)
!      enddo 
!      end
       subroutine cpairLPMXsec(eg, rhoin, xs)
       implicit none
!          compute pair creation probability of very high energy
!          gamma rays in air.   Landau effect is taken into account
!    eg.  input.   gamma energy in gev.   5.e8 < eg < 1.e13
!   rhoin.input.   air density in kg/m**3   
!    xs.  output.  probability  /r.l
!        
!         input outside of this range will be accepted 
!         properly.
!  
!
        real*8  eg, rhoin, xs
!
!       Since the probability of pair creation f(v, Eg, rho)dv
!      (v = Ee/Eg) scales as f(v, Eg*rho)dv,
!      we made an approximate formula for f(v, x)dv for x
!      =1 to 10^10 GeV*gm/cm3
!
        real*8 x, xlog
        integer i
        real*8 xxx(6)
        data ( xxx(i), i=  1,   6)/
     1       59.363582    ,  -19.567178    ,   2.4813214    ,
     2  -0.15128851    ,  0.44302046E-02, -0.50465366E-04               
     * /   
    
       x = eg* rhoin * 1.e-3   !  GeV g/cm3
       if(x .lt. 1.e5) then
           xs = 0.777
       elseif(x .lt. 1.e10) then
           xlog = log(x)
           xs = 0.
           do i =6, 2, -1
              xs =(xs + xxx(i))* xlog
           enddo
           xs = xs + xxx(1)
           xs = exp(xs)
       else
           xs = 2.12e-2 * (x/1.e10)**(-0.5)
       endif

       end
