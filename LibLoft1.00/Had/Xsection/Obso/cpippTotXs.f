!        pi+p total xsection , elastic xs
      subroutine cpippTotXs(p, xs)
      implicit none
#include "Zmass.h"
      real*8 p  ! input.  momentum of n. in GeV
      real*8 xs     ! output. total np cross section in mb
      integer np, i, m
      real*8  error
      parameter (np=53, m=5)
      real*8  px(np), mb(np)
      real*8  Z, Y1, Y2, B, s, rts, s0, eta1, eta2
 
 
      parameter(Z=20.86, Y1=19.24, Y2=6.03, B=0.308, 
     *   s0=5.38**2, eta1=0.458, eta2=0.545)
           data ( px(i), i=           1 ,           np )/
     1   0.1190    ,  0.1275    ,  0.1465    ,  0.1649    ,
     2   0.1811    ,  0.2050    ,  0.2197    ,  0.2389    ,
     3   0.2573    ,  0.2703    ,  0.2854    ,  0.3058    ,
     4   0.3180    ,  0.3356    ,  0.3595    ,  0.3850    ,
     5   0.4185    ,  0.4594    ,  0.5068    ,  0.5645    ,
     6   0.6198    ,  0.6509    ,  0.6737    ,  0.7146    ,
     7   0.7848    ,  0.8327    ,  0.8922    ,  0.9944    ,
     8    1.092    ,   1.193    ,   1.285    ,   1.364    ,
     9    1.439    ,   1.542    ,   1.636    ,   1.718    ,
     a    1.805    ,   2.011    ,   2.241    ,   2.401    ,
     b    2.573    ,   2.825    ,   3.242    ,   4.066    ,
     c    5.201    ,   7.377    ,   9.961    ,   13.85    ,
     d    19.08    ,   27.20    ,   39.73    ,   59.20    ,
     e    76.10                                                         
     * /   
           data ( mb(i), i=           1 ,           np )/
     1    9.845    ,   13.81    ,   19.99    ,   30.48    ,
     2    44.11    ,   75.42    ,   100.4    ,   137.3    ,
     3    174.4    ,   193.5    ,   200.7    ,   190.5    ,
     4    174.4    ,   143.8    ,   113.1    ,   88.09    ,
     5    62.45    ,   45.45    ,   33.24    ,   24.57    ,
     6    19.84    ,   17.15    ,   15.77    ,   14.36    ,
     7    15.21    ,   17.41    ,   20.14    ,   23.30    ,
     8    26.13    ,   30.54    ,   34.97    ,   38.40    ,
     9    40.67    ,   38.20    ,   35.14    ,   32.84    ,
     a    30.84    ,   29.27    ,   29.89    ,   30.83    ,
     b    30.99    ,   30.19    ,   29.10    ,   27.61    ,
     c    26.89    ,   25.38    ,   24.59    ,   24.32    ,
     d    23.81    ,   23.43    ,   23.17    ,   23.40    ,
     e    23.39                                                         
     * /   


      save    
      if(p .gt. 30.) then
!         s =(M+E)^2 - p^2  = M^2 + m^2 +2ME
!           = 
         s = masp**2 + maspic**2 + 2*masp*sqrt(p**2 + maspic**2)
         rts = sqrt(s)
         xs =Z + B*log(s/s0)**2 + Y1*(1./s)**eta1 -Y2*(1./s)**eta2 
      elseif( p .gt. 0.6) then   
         call kpolintplogxyFE(px, 1, mb, 1, np, m, 3,  p, xs, error) 
      else
         call cpippElaXs(p, xs)
      endif
      end
!         
      subroutine cpippElaXs(p, xs)
!           pi+ p elastic cross section in mb
      implicit none
      real*8 p ! input.  momentum of n in GeV
      real*8  xs   ! output np elastic xs. mb.

       integer np, m, i
       parameter (np=54, m=5)
       real*8 px(np), mb(np)
       real*8 error
       real*8 xssave/-1./

           data ( px(i), i=           1 ,           np )/
     1   0.1190    ,  0.1275    ,  0.1465    ,  0.1649    ,
     2   0.1811    ,  0.2050    ,  0.2197    ,  0.2389    ,
     3   0.2573    ,  0.2703    ,  0.2854    ,  0.3058    ,
     4   0.3180    ,  0.3356    ,  0.3595    ,  0.3850    ,
     5   0.4185    ,  0.4594    ,  0.5068    ,  0.5645    ,
     6   0.5900    ,  0.6137    ,  0.6414    ,  0.6638    ,
     7   0.6903    ,  0.7395    ,  0.7805    ,  0.8039    ,
     8   0.8444    ,  0.8914    ,  0.9273    ,  0.9888    ,
     9    1.075    ,   1.146    ,   1.241    ,   1.349    ,
     a    1.438    ,   1.556    ,   1.659    ,   1.931    ,
     b    2.283    ,   2.821    ,   3.716    ,   5.167    ,
     c    7.511    ,   11.03    ,   16.03    ,  19.8157 ,
     d    29.373   ,   42.482   ,   59.9484 ,    84.1809,    
     e    134.34   ,   240.072
     * /   


           data ( mb(i), i=           1 ,           np )/
     1    9.845    ,   13.81    ,   19.99    ,   30.48    ,
     2    44.11    ,   75.42    ,   100.4    ,   137.3    ,
     3    174.4    ,   193.5    ,   200.7    ,   190.5    ,
     4    174.4    ,   143.8    ,   113.1    ,   88.09    ,
     5    62.45    ,   45.45    ,   33.24    ,   24.57    ,
     6    21.01    ,   18.74    ,   16.11    ,   14.29    ,
     7    12.74    ,   10.67    ,   9.417    ,   8.846    ,
     8    8.396    ,   8.752    ,   9.612    ,   11.30    ,
     9    12.73    ,   14.35    ,   16.01    ,   17.31    ,
     a    17.40    ,   16.34    ,   13.83    ,   10.99    ,
     b    8.740    ,   7.710    ,   6.627    ,   5.695    ,
     c    5.075    ,   4.740    ,   4.359    ,   3.94511 ,
     d    3.46975  ,   3.24744  ,   3.20098  ,   3.27183 ,
     e    3.294    ,    3.30042
                            
     * /   

       save
! 
       if( p .gt. 240.) then
!           assume prop.to total
          if( xssave .lt. 0.) then
             call cpippTotXs(px(np), xssave)
          endif
          call cpippTotXs(p, xs)
          xs = xs * mb(np)/xssave
       elseif(p .gt. 0.1) then
          call kpolintplogxyFE(px, 1, mb, 1, np, m, 3, p, xs, error) 
       else
!            get value at 0.1
          xs = mb(1)
       endif
       end
      subroutine cpippInelaXs(p, xs)
      implicit none
      real(8),intent(in)::p
      real(8),intent(out)::xs

      real(8)::txs, exs
      call cpippTotXs(p, txs)
      call cpippElaXs(p, exs)
      xs =max( txs - exs, 0.d0)
      end
