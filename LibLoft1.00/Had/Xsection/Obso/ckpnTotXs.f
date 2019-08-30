!        K+ n total xsection , elastic xs
      subroutine ckpnTotXs(p, xs)
      implicit none
#include "Zmass.h"
      real*8 p  ! input.  momentum of K+. in GeV
      real*8 xs     ! output. total Kn cross section in mb
      integer np, i, m
      real*8  error
      parameter (np=26, m=5)
      real*8  px(np), mb(np)
      real*8  Z, Y1, Y2, B, s, rts, s0, eta1, eta2

      parameter(Z=17.91, Y1=7.1, Y2=13.45, B=0.308, 
     *   s0=5.38**2, eta1=0.458, eta2=0.545)

           data ( px(i), i=           1 ,          np )/
     1  0.499, 0.596, 0.769, 0.919, 0.967,  
     2  1.071, 1.233, 1.492, 2.001, 2.423,  
     3  3.087, 3.883, 5.140, 7.929, 13.895,  
     4  20.630, 32.647, 51.663, 80.719, 124.520,  
     5  207.359, 392.263, 811.319, 2005.980, 4959.750,  
     6  7848.650  
     * /   

           data ( mb(i), i=           1 ,          np )/
     1  10.185, 11.888, 14.742, 17.207, 18.529,  
     2  19.290, 19.418, 18.897, 18.389, 17.896,  
     3  17.652, 17.647, 17.642, 17.752, 17.860,  
     4  17.852, 18.084, 18.319, 18.683, 19.055,  
     5  19.828, 20.910, 21.901, 23.088, 24.670,  
     6  25.501  
     * /   

      save    
      if(p .gt. 250.) then
!         s =(M+E)^2 - p^2  = M^2 + m^2 +2ME
!           = 
         s = masp**2 + maskc**2 + 2*masp*sqrt(p**2 + maskc**2)
         rts = sqrt(s)
         xs =Z + B*log(s/s0)**2 + Y1*(1./s)**eta1 +Y2*(1./s)**eta2 
      elseif( p .gt. 0.6) then   
!                  take log x only
         call kpolintplogxyFE(px, 1, mb, 1, np, m, 1,  p, xs, error) 
      else
         call ckpnElaXs(p, xs)
      endif
      end
!         
      subroutine ckpnElaXs(p, xs)
!           pi+ p elastic cross section in mb
      implicit none
      real*8 p ! input.  momentum of n in GeV
      real*8  xs   ! output np elastic xs. mb.

       integer np, m, i
       parameter (np=30, m=5)
       real*8 px(np), mb(np)
       real*8 error
       real*8 xssave/-1./

           data ( px(i), i=           1 ,           np )/
     1  0.499, 0.619, 0.759, 0.908, 0.992,  
     2  1.142, 1.348, 1.695, 2.053, 2.717,  
     3  4.086, 5.990, 8.451, 11.331, 15.387,  
     4  21.710, 33.066, 49.724, 75.734, 112.445, 
     5  146.967, 197.049, 271.020, 445.603, 661.604,  
     6  982.309, 1678.060, 2903.360, 5023.390, 7749.220  
     * /   
           data ( mb(i), i=           1 ,           np )/
     1 5.160,5.942,6.661,7.318,7.567,
     2 7.316,6.884,6.263,5.510,4.718,
     3 4.150,3.505,3.146,2.901,2.711,
     4 2.550,2.432,2.366,2.365,2.461,
     5 2.477,2.561,2.648,2.793,2.926,
     6 3.025,3.169,3.342,3.408,3.668
     * /   


       save
! 
       if( p .gt. 197.) then
!           assume prop.to total
          if( xssave .lt. 0.) then
             call ckpnTotXs(px(22), xssave)
          endif
          call ckpnTotXs(p, xs)
          xs = xs * mb(22)/xssave
       elseif(p .gt. px(1)) then
!                 take log x only
          call kpolintplogxyFE(px, 1, mb, 1, np, m, 1, p, xs, error) 
       else
!            get value at 0.1
          xs = mb(1)
       endif
       end
      subroutine ckpnInelaXs(p, xs)
      implicit none
      real(8),intent(in)::p
      real(8),intent(out)::xs

      real(8)::txs, exs
      call ckpnTotXs(p, txs)
      call ckpnElaXs(p, exs)
      xs =max( txs - exs, 0.d0)
      end
      
