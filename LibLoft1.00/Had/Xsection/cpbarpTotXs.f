!        pbarp total xsection, pbap ela

!      At low energies, xs is diff.from pp natually.
!      As energy goes higher xs coicides with the one for
!      pp. 
!      However, at E>10^16 eV, xs becomes little bit
!      larger than the one for pp. (due to PDG parameterization) 
      subroutine cpbarpTotXs(p, xs)
      use modpdgXs, only : cSigmaT
      implicit none
#include "Zmass.h"
      real*8 p  ! input.  momentum of pbar. in GeV
      real*8 xs     ! output. total pbarp cross section in mb
      integer np, i, m
      real*8  error
      real(8),parameter::Pnorm=300.d0
      parameter (np=34, m=5)
      real*8  px(np), mb(np)

           data ( px(i), i=           1 ,           np )/
     1   9.7289E-02,  0.1485    ,  0.1925    ,  0.2330    ,
     2   0.2821    ,  0.3558    ,  0.4426    ,  0.5736    ,
     3   0.6850    ,  0.8180    ,   1.004    ,   1.319    ,
     4    1.686    ,   2.098    ,   2.911    ,   4.207    ,
     5    5.918    ,   9.286    ,   13.79    ,   21.64    ,
     6    32.15    ,   49.08    ,   70.95    ,   109.8    ,
     7    144.3    ,   192.2    ,   274.0    ,   418.3    ,
     8    702.7    ,   1133.    ,   1877.    ,   2986.    ,
     9    5015.    ,  1.5568E+04                                        
     * /   
           data ( mb(i), i=           1 ,           np )/
     1    445.7    ,   380.5    ,   321.4    ,   283.2    ,
     2    249.5    ,   214.1    ,   179.9    ,   153.6    ,
     3    134.6    ,   123.0    ,   113.7    ,   103.9    ,
     4    95.54    ,   86.43    ,   75.75    ,   66.40    ,
     5    59.76    ,   54.64    ,   51.30    ,   47.16    ,
     6    45.22    ,   43.82    ,   42.92    ,   42.03    ,
     7    41.82    ,   41.39    ,   41.84    ,   42.07    ,
     8    42.31    ,   43.22    ,   43.70    ,   45.12    ,
     9    46.59    ,   50.73                                            
     * /   
      if(p .gt. Pnorm) then
         xs = cSigmaT('pb', 'p', p)
      elseif( p .gt. 0.1) then   
         call kpolintplogxyFE(px, 1, mb, 1, np, m, 3, p, xs, error) 
      else
         xs =mb(1)
      endif
      end      subroutine cpbarpTotXs
      subroutine cpbarpElaXs(p, xs)
!           pbarp elastic cross section in mb
      implicit none
      real*8 p ! input.  momentum of pbar in GeV
      real*8  xs   ! output pbarp elastic xs. mb.

       integer np, m, i
       parameter (np=35, m=5)
       real*8 px(np), mb(np)
       real*8 error
       real(8),parameter::Pnorm=300.d0
       real(8),save:: xsn, Norm, ppp
       real(8):: xspel, xspt, xspit 
       real(8),save::xsnorm

       logical,save:: first=.true. 

           data ( px(i), i=           1 ,           np )/
     1   9.8633E-02,  0.1485    ,  0.2299    ,  0.3062    ,
     2   0.4078    ,  0.5358    ,  0.6850    ,  0.8758    ,
     3    1.166    ,   1.575    ,   1.986    ,   2.505    ,
     4    3.572    ,   4.957    ,   6.784    ,   9.286    ,
     5    11.87    ,   15.18    ,   21.64    ,   32.60    ,
     6    56.27    ,   83.59    ,   127.6    ,   289.4    ,
     7    486.1    ,   936.0    ,   1551.    ,   3241.    ,
     8    6680.    ,  1.4542E+04,  3.3434E+04,  6.4372E+04,
     9   1.5845E+05,  6.2031E+05,  1.7030E+06                           
     * /   
           data ( mb(i), i=           1 ,           np )/
     1    106.0    ,   96.90    ,   86.76    ,   78.90    ,
     2    68.79    ,   60.62    ,   53.69    ,   47.81    ,
     3    43.25    ,   37.51    ,   32.19    ,   26.90    ,
     4    20.89    ,   16.91    ,   13.91    ,   11.81    ,
     5    10.80    ,   9.772    ,   9.078    ,   8.568    ,
     6    7.960    ,   7.513    ,   7.358    ,   7.170    ,
     7    7.135    ,   7.252    ,   7.488    ,   8.239    ,
     8    9.064    ,   10.13    ,   11.21    ,   12.07    ,
     9    13.43    ,   15.42    ,   17.42                               
     * /   
       save
       if( p >= Pnorm ) then
          if( first ) then
             ppp = Pnorm
             call cppElaXs(PPP, xspel)
             call cppTotXs(ppp, xspt)
             call cpbarpTotXs(ppp, xspit)
             call kpolintpFE(px, 1, mb, 1, np, m,
     *           ppp, xsnorm, error) 
             xs = xspit*xspel/xspt
             Norm = xs - xsnorm
             first = .false.
          endif
          ppp = p
          call cppElaXs(PPP, xspel)
          call cppTotXs(ppp, xspt)
          call cpbarpTotXs(ppp, xspit)
          xs = xspit*xspel/xspt
          xs = xs - Norm
       elseif(p .gt. 0.1) then
          call kpolintplogxyFE(px, 1, mb, 1, np, m, 3, p, xs, error) 
       else
!            get value at ~0.1
          xs = 106.
       endif
       end      subroutine cpbarpElaXs
      subroutine cpbarpInelaXs(p, xs)
      implicit none
      real(8),intent(in)::p
      real(8),intent(out)::xs

      real(8)::txs, exs
      call cpbarpTotXs(p, txs)
      call cpbarpElaXs(p, exs)
      xs =max( txs - exs, 0.d0)
      end      subroutine cpbarpInelaXs
