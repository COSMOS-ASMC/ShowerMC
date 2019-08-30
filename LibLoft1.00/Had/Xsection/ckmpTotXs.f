!        K- -p total xsection , elastic xs
      subroutine ckmpTotXs(p, xs)
      use modpdgXs, only: cSigmaT
      implicit none
#include "Zmass.h"
      real*8 p  ! input.  momentum of n. in GeV
      real*8 xs     ! output. total np cross section in mb
      integer np, i, m
      real*8  error
      parameter (np=45, m=3)
      real*8  px(np), mb(np)
      real(8),parameter:: Pnorm=10.

           data ( px(i), i=           1 ,           np )/
     1   0.9980E-01,  0.1155    ,  0.2860    ,  0.3614    ,
     2   0.4023    ,  0.4349    ,  0.4500    ,  0.48      ,
     *   0.5158    ,
     3   0.5913    ,  0.6455    ,  0.6810    ,  0.7291    ,
     4   0.7618    ,  0.7883    ,  0.8358    ,  0.9258    ,
     5   0.9816    ,   1.036    ,   1.082    ,   1.164    ,
     6    1.228    ,   1.321    ,   1.429    ,   1.522    ,
     7    1.622    ,   1.788    ,   1.933    ,   2.162    ,
     8    2.454    ,   3.026    ,   3.769    ,   4.927    ,
     9    7.207    ,   9.422    ,   12.68    ,   16.99    ,
     a    21.37    ,   29.33    ,   43.11    ,   68.50    ,
     b    99.22    ,   153.1    ,   236.3    ,   369.9 
     * /   
           data ( mb(i), i=           1 ,           np )/
     1    107.3    ,   98.89    ,   92.45    ,   74.10    ,
     2    68.94    ,   63.22    ,   49.47    ,   42.0,
     *    40.81    ,
     3    36.01    ,   34.48    ,   34.15    ,   35.16    ,
     4    39.84    ,   41.20    ,   39.84    ,   44.71    ,
     5    48.76    ,   51.41    ,   47.83    ,   40.22    ,
     6    33.67    ,   30.28    ,   31.78    ,   33.02    ,
     7    33.83    ,   32.39    ,   30.14    ,   29.71    ,
     8    28.31    ,   26.60    ,   26.21    ,   24.86    ,
     9    23.58    ,   22.69    ,   21.94    ,   21.21    ,
     a    20.91    ,   20.61    ,   20.22    ,   20.22    ,
     b    20.32    ,   20.71    ,   21.01    ,   21.73    
     * /   

      save    
      if(p .gt. Pnorm) then
         xs = cSigmaT('K-', 'p', p)
      elseif( p .gt. 0.2) then   
!         call kpolintplogxyFE(px, 1, mb, 1, np, m, 3,  p, xs, error) 
         call kpolintpFE(px, 1, mb, 1, np, m,  p, xs, error) 
      else
         call ckmpElaXs(p, xs)
      endif
      end      subroutine ckmpTotXs
!         
      subroutine ckmpElaXs(p, xs)
      use modpdgXs, Mpdg => M
!           pi+ p elastic cross section in mb
      implicit none
#include "Zmass.h"
      real*8 p ! input.  momentum of n in GeV
      real*8  xs   ! output np elastic xs. mb.

       integer np, m, i
       parameter (np=40, m=5)
       real*8 px(np), mb(np)
       real*8 error

       real(8),parameter:: Pnorm=50.0
       real(8)::spip, spp, Epp, PnormPP, PPP             
       real(8):: xspel, xspt, xspit, Norm
       real(8),save::xsnorm

       logical,save:: first = .true.


           data ( px(i), i=           1 ,           np )/
     1   0.9980E-01,  0.1155    ,  0.1460    ,  0.1918    ,
     2   0.2658    ,  0.3278    ,  0.3794    ,  0.4286    ,
     3   0.4478    ,  0.4611    ,  0.5010    ,  0.5714    ,
     4   0.6330    ,  0.6744    ,  0.7081    ,  0.7544    ,
     5   0.7844    ,  0.8439    ,  0.9395    ,   1.041    ,
     6    1.114    ,   1.205    ,   1.271    ,   1.321    ,
     7    1.478    ,   1.736    ,   1.886    ,   2.183    ,
     8    2.589    ,   3.419    ,   4.094    ,   5.000    ,
     9    6.380    ,   8.974    ,   13.85    ,   23.79    ,
     a    40.07    ,   61.84    ,   96.36    ,   165.5                  
     * /   
           data ( mb(i), i=           1 ,           np )/
     1    107.3    ,   98.89    ,   80.41    ,   60.54    ,
     2    45.36    ,   37.60    ,   32.71    ,   32.08    ,
     3    28.72    ,   24.39    ,   21.01    ,   18.19    ,
     4    16.05    ,   14.93    ,   14.03    ,   15.82    ,
     5    18.54    ,   20.32    ,   21.32    ,   21.21    ,
     6    18.36    ,   15.37    ,   12.74    ,   10.46    ,
     7    8.671    ,   8.629    ,   8.028    ,   6.983    ,
     8    5.844    ,   4.845    ,   4.337    ,   3.977    ,
     9    3.683    ,   3.394    ,   3.009    ,   2.605    ,
     a    2.435    ,   2.423    ,   2.459    ,   2.543                  
     * /   

       save
! 
       if( p .gt. Pnorm) then
          if( first ) then
             ppp = Pnorm
             call cppElaXs(PPP, xspel)
             call cppTotXs(ppp, xspt)
             call ckmpTotXs(ppp, xspit)
             call kpolintpFE(px, 1, mb, 1, np, m,
     *           ppp, xsnorm, error) 
             xs = xspit*xspel/xspt
             Norm = xs - xsnorm

             first=.false.
          endif
          ppp = p
          call cppElaXs(PPP, xspel)
          call cppTotXs(ppp, xspt)
          call ckmpTotXs(ppp, xspit)
          xs = xspit*xspel/xspt
          xs = xs - Norm
       elseif(p .gt. px(1)) then
!          call kpolintplogxyFE(px, 1, mb, 1, np, m, 3, p, xs, error) 
          call kpolintpFE(px, 1, mb, 1, np, m, p, xs, error) 
       else
!            get value at 0.1
          xs = mb(1)
       endif
       end      subroutine ckmpElaXs
      subroutine ckmpInelaXs(p, xs)
      implicit none
      real(8),intent(in)::p
      real(8),intent(out)::xs

      real(8)::txs, exs
      call ckmpTotXs(p, txs)
      call ckmpElaXs(p, exs)
      xs =max( txs - exs, 0.d0)
      end      subroutine ckmpInelaXs

