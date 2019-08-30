!        pi-p total xsection , elastic xs
      subroutine cpimpTotXs(p, xs)
      use modpdgXs, only: cSigmaT
      implicit none
#include "Zmass.h"
      real*8 p  ! input.  momentum of n. in GeV
      real*8 xs     ! output. total np cross section in mb
      integer np, i, m
      real*8  error
      parameter (np=108, m=5)
      real*8  px(np), mb(np)
      real(8),parameter:: Pnorm = 100.

          data ( px(i), i=           1 ,           np )/
     1 0.158, 0.168, 0.180, 0.193, 0.206, 0.216, 0.225, 0.237,
     2 0.247, 0.259, 0.264, 0.271, 0.282, 0.293, 0.301, 0.307,
     3 0.316, 0.328, 0.339, 0.348, 0.360, 0.373, 0.385, 0.396,
     4 0.411, 0.432, 0.446, 0.461, 0.485, 0.511, 0.539, 0.564,
     5 0.594, 0.615, 0.625, 0.642, 0.651, 0.672, 0.685, 0.693,
     6 0.703, 0.734, 0.744, 0.760, 0.766, 0.775, 0.793, 0.808,
     7 0.817, 0.836, 0.857, 0.873, 0.895, 0.917, 0.951, 0.977,
     8 0.983, 1.010, 1.028, 1.050, 1.070, 1.095, 1.107, 1.133,
     9 1.152, 1.194, 1.260, 1.332, 1.386, 1.453, 1.492, 1.557,
     a 1.640, 1.737, 1.834, 1.942, 2.042, 2.200, 2.302, 2.418,
     b 2.540, 2.689, 2.836, 3.002, 3.173, 3.409, 3.726, 4.081,
     c 4.446, 4.754, 5.210, 5.721, 6.372, 6.883,
     d    8.814    ,
     e    11.61    ,   16.06    ,   23.46    ,   32.47    ,
     f    46.73    ,   70.64    ,   99.20    ,   130.0    ,
     g    175.5    ,   218.0    ,   284.3    ,   401.2    ,
     h    518.1                                                         
     * /   
           data ( mb(i), i=           1 ,           np )/
     1 13.60, 16.34, 19.99, 24.29, 30.30, 36.96, 44.13, 50.01,
     2 57.84, 63.74, 67.18, 70.42, 71.46, 68.90, 64.93, 62.13,
     3 57.36, 53.36, 48.12, 44.52, 40.45, 36.29, 33.56, 31.29,
     4 29.20, 27.18, 26.32, 26.16, 26.63, 27.25, 28.06, 29.45,
     5 30.84, 32.67, 34.37, 36.91, 39.13, 41.94, 44.22, 46.08,
     6 47.57, 47.05, 45.75, 44.43, 41.80, 40.05, 38.42, 36.88,
     7 35.96, 36.28, 37.39, 39.32, 42.16, 45.92, 50.66, 54.82,
     8 57.46, 59.66, 56.91, 53.99, 50.22, 45.91, 42.98, 40.43,
     9 38.62, 37.15, 36.16, 36.52, 36.79, 36.16, 35.43, 34.51,
     a 34.89, 35.23, 35.74, 36.24, 36.11, 35.86, 34.92, 34.33,
     b 33.70, 32.85, 32.64, 32.87, 32.36, 31.63, 30.92, 30.40,
     c 29.99, 29.51, 28.83, 28.52, 27.94, 27.77, 
     d 26.89    ,
     e 26.02    ,   25.29    ,   24.98    ,   24.77    ,
     f 24.36    ,   24.16    ,   24.15    ,   24.25    ,
     g  24.34    ,   24.64    ,   24.74    ,   25.14    ,
     h  25.55                                                         
     * /   
      save    
      if(p .gt. Pnorm) then
         xs = cSigmaT('pi-', 'p', p)
      elseif( p .gt. 0.2) then   
!         call kpolintplogxyFE(px, 1, mb, 1, np, m, 0,  p, xs, error) 
         call kpolintpFE(px, 1, mb, 1, np, m,   p, xs, error) 
      else
         call cpimpElaXs(p, xs)
      endif
      end      subroutine cpimpTotXs
!         
      subroutine cpimpElaXs(p, xs)
      use modpdgXs, Mpdg=>M
!           pi- p elastic cross section in mb
      implicit none
#include "Zmass.h"
      real*8 p ! input.  momentum of n in GeV
      real*8  xs   ! output np elastic xs. mb.

       integer np, m, i
       parameter (np=52, m=5)
       real*8 px(np), mb(np)
       real*8 error
       real*8 xssave/-1./
       real(8),parameter:: Pnorm=100.0
       real(8),parameter:: sm=(maspic + masp + Mpdg)**2
       real(8),parameter:: smpp=(masp + masp + Mpdg)**2
       real(8),save::  Norm
       real(8):: PPP             
       real(8):: xspel, xspt, xspit 
       real(8),save::xsnorm

       logical,save:: first = .true.

           data ( px(i), i=           1 ,           np )/
     1   0.1320    ,  0.1553    ,  0.1731    ,  0.1845    ,
     2   0.2017    ,  0.2171    ,  0.2280    ,  0.2418    ,
     3   0.2527    ,  0.2680    ,  0.2871    ,  0.3076    ,
     4   0.3232    ,  0.3480    ,  0.3785    ,  0.4056    ,
     5   0.4454    ,  0.5187    ,  0.5611    ,  0.6159    ,
     6   0.6374    ,  0.6760    ,  0.7242    ,  0.7608    ,
     7   0.7915    ,  0.8480    ,  0.9085    ,  0.9354    ,
     8   0.9823    ,   1.007    ,   1.052    ,   1.089    ,
     9    1.128    ,   1.221    ,   1.354    ,   1.531    ,
     a    1.988    ,   2.361    ,   2.904    ,   3.901    ,
     b    5.191    ,   7.546    ,   10.29    ,   14.03    ,
     c    19.51    ,   27.40    ,   38.66    ,   63.24    ,
     d    100.9    ,   170.9    ,   263.4    ,   366.3                  
     * /   
           data ( mb(i), i=           1 ,           np )/
     1    21.81    ,   17.87    ,   15.43    ,   13.60    ,
     2    11.46    ,   11.37    ,   12.90    ,   16.87    ,
     3    19.70    ,   22.07    ,   24.94    ,   23.46    ,
     4    19.14    ,   15.05    ,   11.98    ,   10.78    ,
     5    9.501    ,   10.06    ,   11.74    ,   13.93    ,
     6    16.12    ,   18.75    ,   20.09    ,   18.44    ,
     7    16.12    ,   14.80    ,   14.98    ,   19.05    ,
     8    23.44    ,   25.33    ,   23.63    ,   20.50    ,
     9    16.38    ,   13.70    ,   12.27    ,   10.47    ,
     a    8.894    ,   8.098    ,   7.434    ,   6.578    ,
     b    5.566    ,   5.006    ,   4.689    ,   4.304    ,
     c    3.871    ,   3.626    ,   3.467    ,   3.301    ,
     d    3.220    ,   3.180    ,   3.325    ,   3.577                  
     * /   

       save
       if( p .gt. Pnorm) then
          if( first ) then
             ppp = Pnorm
             call cppElaXs(PPP, xspel)
             call cppTotXs(ppp, xspt)
             call cpimpTotXs(ppp, xspit)
             call kpolintpFE(px, 1, mb, 1, np, m,
     *           ppp, xsnorm, error) 
             xs = xspit*xspel/xspt
             Norm = xs - xsnorm
             first=.false.
          endif
          ppp = p
          call cppElaXs(PPP, xspel)
          call cppTotXs(ppp, xspt)
          call cpimpTotXs(ppp, xspit)
          xs = xspit*xspel/xspt
          xs = xs - Norm
       elseif(p .gt. px(1)) then
!          call kpolintplogxyFE(px, 1, mb, 1, np, m, 0, p, xs, error) 
          call kpolintpFE(px, 1, mb, 1, np, m,  p, xs, error) 
       else
!            get value at 0.1
          xs = mb(1)
       endif
       end    subroutine cpimpElaXs
      subroutine cpimpInelaXs(p, xs)
      implicit none
      real(8),intent(in)::p
      real(8),intent(out)::xs

      real(8)::txs, exs
      call cpimpTotXs(p, txs)
      call cpimpElaXs(p, exs)
      xs =max( txs - exs, 0.d0)
      end   subroutine cpimpInelaXs

