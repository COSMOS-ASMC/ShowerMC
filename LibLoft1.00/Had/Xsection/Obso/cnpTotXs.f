!        np total xsection
      subroutine cnpTotXs(p, xs)
      implicit none
#include "Zmass.h"
      real*8 p  ! input.  momentum of n. in GeV
      real*8 xs     ! output. total np cross section in mb
      integer np, i, m
      real*8  error
      parameter (np=54, m=5)
      real*8  px(np), mb(np)
      real*8  Z, Y1, Y2, B, s, rts, s0, eta1, eta2
      parameter(Z= 35.45, Y1=42.53, Y2=33.34, B=0.308, 
     *   s0=5.38**2, eta1=0.458, eta2=0.545)

           data ( px(i), i=           1 ,           np )/
     1   0.1001    ,  0.1126    ,  0.1224    ,  0.1337    ,
     2   0.1489    ,  0.1716    ,  0.1912    ,  0.2119    ,
     3   0.2395    ,  0.2590    ,  0.2856    ,  0.3135    ,
     4   0.3407    ,  0.3739    ,  0.4207    ,  0.4595    ,
     5   0.4996    ,  0.5512    ,  0.6082    ,  0.6646    ,
     6   0.7227    ,  0.7781    ,  0.8378    ,  0.9246    ,
     7    1.051    ,   1.225    ,   1.400    ,   1.584    ,
     8    1.827    ,   2.161    ,   2.506    ,   2.949    ,
     9    3.470    ,   4.063    ,   4.950    ,   5.767    ,
     a    6.854    ,   8.432    ,   10.07    ,   12.09    ,
     b    15.02    ,   18.75    ,   23.07    ,   28.66    ,
     c    36.68    ,   47.41    ,   61.27    ,   78.80    ,
     d    99.85    ,   130.3    ,   174.4    ,   236.8    ,
     e    318.3    ,   395.5                                            
     * /   
           data ( mb(i), i=           1 ,           np )/
     1    1465.    ,   1264.    ,   1097.    ,   946.4    ,
     2    782.0    ,   619.0    ,   489.9    ,   395.0    ,
     3    301.3    ,   242.9    ,   197.0    ,   152.1    ,
     4    122.7    ,   100.7    ,   78.76    ,   65.88    ,
     5    56.14    ,   48.73    ,   44.17    ,   42.57    ,
     6    41.28    ,   38.34    ,   35.40    ,   34.12    ,
     7    34.12    ,   36.07    ,   38.13    ,   40.30    ,
     8    41.06    ,   42.61    ,   43.41    ,   43.42    ,
     9    43.42    ,   42.64    ,   42.12    ,   42.13    ,
     a    41.11    ,   40.86    ,   40.62    ,   39.64    ,
     b    39.16    ,   39.90    ,   39.66    ,   39.67    ,
     c    39.43    ,   39.44    ,   38.73    ,   39.46    ,
     d    39.23    ,   39.48    ,   40.47    ,   39.99    ,
     e    40.00    ,   41.00                                            
     * /   

      if(p .gt. 390.) then
!         s =(m+E)^2 - p^2  = m^2 + p^2+m^2+2mE-p^2
!           =  2m(m+E)
         s = (masp + sqrt(p**2+masp**2))*masp*2.
         rts = sqrt(s)
         xs =Z + B*log(s/s0)**2 + Y1*(1./s)**eta1 -Y2*(1./s)**eta2 
      elseif( p .gt. 0.4) then   
         call kpolintpFE(px, 1, mb, 1, np, m, p, xs, error) 
      else
         call cnpElaXs(p, xs)
      endif
      end
! 
      subroutine cnpElaXs(p, xs)
!           np elastic cross section in mb
      implicit none
      real*8 p ! input.  momentum of n in GeV
      real*8  xs   ! output np elastic xs. mb.

       integer np, m, i
       parameter (np=37, m=5)
       real*8 px(np), mb(np)
       real*8 error
       real*8 xssave/-1./
       data ( px(i), i=           1 ,           np )/
     1   0.9860E-02,  0.1271E-01,  0.1942E-01,  0.3230E-01,
     2   0.4935E-01,  0.7025E-01,  0.1000    ,  0.1187    ,
     3   0.1294    ,  0.1472    ,  0.1692    ,  0.1988    ,
     4   0.2386    ,  0.2833    ,  0.3257    ,  0.3785    ,
     5   0.4803    ,  0.5145    ,  0.5818    ,  0.6548    ,
     6   0.7628    ,  0.8543    ,  0.9906    ,   1.132    ,
     7    1.243    ,   1.413    ,   1.575    ,   1.712    ,
     8    2.180    ,   2.960    ,   3.788    ,   5.068    ,
     9    8.427    ,   16.91    ,   37.07    ,   66.05    ,
     *    100.0                                                         
     * /   
       data ( mb(i), i=           1 ,           np )/
     1    8452.    ,   7556.    ,   6039.    ,   4437.    ,
     2    3306.    ,   2362.    ,   1549.    ,   1169.    ,
     3    1022.    ,   813.3    ,   621.7    ,   438.5    ,
     4    301.1    ,   198.5    ,   136.3    ,   100.1    ,
     5    60.07    ,   50.88    ,   42.83    ,   40.02    ,
     6    34.96    ,   32.67    ,   31.88    ,   30.92    ,
     7    31.31    ,   29.80    ,   26.85    ,   24.78    ,
     8    21.65    ,   18.23    ,   15.92    ,   12.92    ,
     9    11.43    ,   10.17    ,   9.629    ,   9.575    ,
     *    9.578                                                         
     * /   
       save
       if( p .gt. 100.) then
!            assume prop. to Total
          if(xssave .lt. 0.) then
             call cnpTotXs(px(np), xssave)          
          endif
          call  cnpTotXs(p, xs)
          xs = xs * mb(np)/xssave
       elseif(p .gt. 0.01) then
          call kpolintplogxyFE(px, 1, mb, 1, np, m, 3, p, xs, error) 
       else
!            get value at 0.01
          xs = mb(1)
       endif
       end
      subroutine cnpInelaXs(p, xs)
      implicit none
      real(8),intent(in)::p
      real(8),intent(out)::xs

      real(8)::txs, exs
      call cnpTotXs(p, txs)
      call cnpElaXs(p, exs)
      xs =max( txs - exs, 0.d0)
      end
