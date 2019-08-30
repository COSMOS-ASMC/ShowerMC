!        K+ p total xsection , elastic xs
      subroutine ckppTotXs(p, xs)
      use modpdgXs, only: cSigmaT
      implicit none
#include "Zmass.h"
      real*8 p  ! input.  momentum of n. in GeV
      real*8 xs     ! output. total np cross section in mb
      integer np, i, m
      real*8  error
      parameter (np=40, m=5)
      real*8  px(np), mb(np)
      real(8),parameter:: Pnorm=10.

           data ( px(i), i=           1 ,          np )/
     1   0.1005    ,  0.1706    ,  0.2619    ,  0.3586    ,
     2   0.4511    ,  0.5370    ,  0.6587    ,  0.7537    ,
     3   0.8003    ,  0.8846    ,  0.9779    ,   1.034    ,
     4    1.065    ,   1.131    ,   1.166    ,   1.263    ,
     5    1.388    ,   1.466    ,   1.549    ,   1.728    ,
     6    1.938    ,   2.173    ,   2.377    ,   2.653    ,
     7    3.319    ,   3.913    ,   4.824    ,   6.409    ,
     8    8.305    ,   10.19    ,   12.50    ,   15.33    ,
     9    20.07    ,   27.34    ,   41.97    ,   59.50    ,
     a    85.63    ,   126.3    ,   203.9    ,   302.5                  
     * /   

           data ( mb(i), i=           1 ,          np )/
     1    10.98    ,   11.55    ,   11.91    ,   12.33    ,
     2    12.54    ,   12.56    ,   12.46    ,   12.59    ,
     3    13.13    ,   14.05    ,   15.20    ,   16.17    ,
     4    16.95    ,   17.81    ,   18.31    ,   18.54    ,
     5    18.33    ,   18.02    ,   17.81    ,   17.89    ,
     6    17.73    ,   17.55    ,   17.49    ,   17.36    ,
     7    17.28    ,   17.20    ,   17.25    ,   17.25    ,
     8    17.33    ,   17.38    ,   17.43    ,   17.45    ,
     9    17.50    ,   17.71    ,   17.97    ,   18.30    ,
     a    18.67    ,   19.19    ,   19.92    ,   20.62                  
     * /   
      save    
      if(p .gt. Pnorm ) then
         xs = cSigmaT('K+', 'p', p)
      elseif( p .gt. 0.6) then   
!                  take log x only
         call kpolintplogxyFE(px, 1, mb, 1, np, m, 1,  p, xs, error) 
      else
         call ckppElaXs(p, xs)
      endif
      end      subroutine ckppTotXs
!         
      subroutine ckppElaXs(p, xs)
      use modpdgXs, Mpdg => M
      implicit none
!           K+ p elastic cross section in mb
#include "Zmass.h"
      real*8 p ! input.  momentum of n in GeV
      real*8  xs   ! output np elastic xs. mb.

       integer np, m, i
       parameter (np=27, m=5)
       real*8 px(np), mb(np)
       real*8 error
       real(8),parameter:: Pnorm=100.0
       real(8),save:: Norm
       real(8)::spip, spp, Epp, PnormPP, PPP             
       real(8):: xspel, xspt, xspit 
       real(8),save::xsnorm
       logical,save:: first = .true.
           data ( px(i), i=           1 ,           np )/
     1   0.1005    ,  0.1706    ,  0.2619    ,  0.3586    ,
     2   0.4511    ,  0.5370    ,  0.6330    ,  0.7725    ,
     3   0.9106    ,   1.011    ,   1.203    ,   1.410    ,
     4    1.669    ,   1.996    ,   2.497    ,   3.107    ,
     5    4.066    ,   5.427    ,   7.428    ,   9.819    ,
     6    14.92    ,   22.57    ,   34.81    ,   52.64    ,
     7    83.26    ,   133.7    ,   250.5                               
     * /   
           data ( mb(i), i=           1 ,           np )/
     1    10.98    ,   11.55    ,   11.91    ,   12.33    ,
     2    12.54    ,   12.56    ,   12.46    ,   12.17    ,
     3    12.06    ,   11.56    ,   10.62    ,   9.627    ,
     4    8.319    ,   6.958    ,   5.701    ,   4.811    ,
     5    4.129    ,   3.525    ,   3.366    ,   3.312    ,
     6    3.387    ,   2.992    ,   2.492    ,   2.358    ,
     7    2.433    ,   2.429    ,   2.686                               
     * /   

       save
! 
       if( p .gt. Pnorm) then
          if( first ) then
             ppp = Pnorm
             call cppElaXs(PPP, xspel)
             call cppTotXs(ppp, xspt)
             call ckppTotXs(ppp, xspit)
             call kpolintpFE(px, 1, mb, 1, np, m,
     *           ppp, xsnorm, error) 
             xs = xspit*xspel/xspt
             Norm = xs - xsnorm
             first=.false.
          endif
          ppp = p
          call cppElaXs(PPP, xspel)
          call cppTotXs(ppp, xspt)
          call ckppTotXs(ppp, xspit)
          xs = xspit*xspel/xspt
          xs = xs - Norm
       elseif(p .gt. px(1)) then
!                 take log x only
          call kpolintplogxyFE(px, 1, mb, 1, np, m, 1, p, xs, error) 
       else
!            get value at 0.1
          xs = mb(1)
       endif
       end      subroutine ckppElaXs
      subroutine ckppInelaXs(p, xs)
      implicit none
      real(8),intent(in)::p
      real(8),intent(out)::xs

      real(8)::txs, exs
      call ckppTotXs(p, txs)
      call ckppElaXs(p, exs)
      xs =max( txs - exs, 0.d0)
      end      subroutine ckppInelaXs

