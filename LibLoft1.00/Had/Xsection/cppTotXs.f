!        pp total xsection, pp ela.xs
      subroutine cppTotXs(p, xs)
      use modpdgXs, only : cSigmaT, csOfstu
      implicit none
#include "Zmass.h"
#include "Zevhnp.h"
      real*8 p  ! input.  momentum of p. in GeV
      real*8 xs     ! output. total pp cross section in mb
      integer np, i, m
      real*8  error
      parameter (np=39, m=5)
      real*8  px(np), mb(np)
           data ( px(i), i=           1 ,           np )/
     1   0.1015    ,  0.1370    ,  0.1900    ,  0.2363    ,
     2   0.2746    ,  0.3147    ,  0.3706    ,  0.4306    ,
     3   0.5003    ,  0.5892    ,  0.7329    ,  0.8871    ,
     4   0.9894    ,   1.089    ,   1.150    ,   1.248    ,
     5    1.336    ,   1.552    ,   1.904    ,   2.305    ,
     6    2.752    ,   3.331    ,   4.376    ,   5.908    ,
     7    7.760    ,   10.62    ,   15.56    ,   21.88    ,
     8    29.54    ,   38.80    ,   68.81    ,   100.8    ,
     9    151.8    ,   213.4    ,   292.1    ,   433.8    ,
     a    671.2    ,   1067.    ,   2082.                               
     * /   
           data ( mb(i), i=           1 ,           np )/
     1    486.6    ,   316.8    ,   155.7    ,   91.67    ,
     2    69.23    ,   52.00    ,   35.88    ,   28.12    ,
     3    24.76    ,   23.61    ,   23.11    ,   24.89    ,
     4    27.97    ,   31.60    ,   35.69    ,   39.47    ,
     5    44.59    ,   47.77    ,   48.02    ,   46.77    ,
     6    44.83    ,   43.19    ,   42.07    ,   40.75    ,
     7    40.53    ,   40.11    ,   39.26    ,   39.06    ,
     8    39.06    ,   38.65    ,   38.65    ,   38.65    ,
     9    38.85    ,   39.26    ,   39.47    ,   40.53    ,
     a    40.97    ,   42.07    ,   43.42                               
     * /   

      save    


      real(8):: s

!      if(p .gt. 2000.) then
      if(p .gt. 50. ) then
         if( TotXSopt == 1 ) then  ! pdg parameterization (lowest value)
            xs = cSigmaT('p','p', p) 
         else
            s = csOfstu(p, masp, masp)
            if( TotXSopt == 2 ) then  ! midlle of 1 and 3
               xs = 42.6*(s)**(-0.46) - 33.4*(s)**(-0.545) + 35.5 +
     *           0.295*log(s/29.1)**2
            elseif( TotXSopt == 3) then  ! TOTEM (COMPETE) fitting      
               xs = 42.6*(s)**(-0.46) - 33.4*(s)**(-0.545) + 35.5 +
     *         0.307*log(s/29.1)**2
            else
               write(0,*) ' TotXSopt=', TotXSopt, ' invalid' 
               stop
            endif
         endif
      elseif( p .gt. 0.46) then   
         call kpolintplogxyFE(px, 1, mb, 1, np, m, 3,  p, xs, error) !         call kpolintpFE(px, 1, mb, 1, np, m,  p, xs, error) 
      else
         call cppElaXs(p, xs)
      endif
      end
!         
      subroutine cppElaXs(p, xs)
      use modpdgXs, only : csOfstu
!           pp elastic cross section in mb
      implicit none
#include "Zmass.h"
      real*8 p ! input.  momentum of n in GeV
      real*8  xs   ! output np elastic xs. mb.

       integer np, m, i
       parameter (np=38, m=5)
       real*8 px(np), mb(np)
       real*8 error
       real*8 xssave/-1./
      real(8):: s

           data ( px(i), i=           1 ,           np )/
     1   9.9947E-02,  0.1382    ,  0.1875    ,  0.2380    ,
     2   0.2746    ,  0.3108    ,  0.3452    ,  0.3834    ,
     3   0.4258    ,  0.4638    ,  0.5349    ,  0.6110    ,
     4   0.7181    ,  0.8682    ,  0.9637    ,   1.132    ,
     5    1.396    ,   1.656    ,   1.984    ,   2.399    ,
     6    2.929    ,   3.823    ,   4.802    ,   6.510    ,
     7    9.982    ,   13.79    ,   20.35    ,   30.91    ,
     8    45.20    ,   68.65    ,   102.3    ,   145.4    ,
     9    204.6    ,   296.3    ,   490.0    ,   1067.    ,
     a    1621.    ,   2055.                                            
     * /   
           data ( mb(i), i=           1 ,           np )/
     1    499.9    ,   316.6    ,   155.5    ,   92.89    ,
     2    70.49    ,   53.77    ,   42.57    ,   33.88    ,
     3    28.44    ,   25.84    ,   23.86    ,   23.11    ,
     4    22.74    ,   23.35    ,   23.97    ,   24.74    ,
     5    24.09    ,   23.21    ,   21.77    ,   19.58    ,
     6    16.96    ,   14.61    ,   12.39    ,   10.74    ,
     7    9.910    ,   9.394    ,   8.535    ,   7.631    ,
     8    7.311    ,   6.967    ,   6.817    ,   6.851    ,
     9    6.776    ,   6.992    ,   7.062    ,   7.284    ,
     a    7.516    ,   7.716                                            
     * /   
           
       save


 
!       if( p .gt. 2000.) then
       if( p .gt. 50. ) then
          s = csOfstu(p, masp, masp)
          xs = 11.4 - 1.52*log(s)+0.130*log(s)**2  ! TOTEM
       elseif(p .gt. 0.1) then
          call kpolintplogxyFE(px, 1, mb, 1, np, m, 3, p, xs, error) 
!          call kpolintpFE(px, 1, mb, 1, np, m,  p, xs, error) 
       else
!            get value at 0.1
          xs = mb(1)
       endif
       end
      subroutine cppInelaXs(p, xs)
      implicit none
      real(8),intent(in)::p
      real(8),intent(out)::xs

      real(8)::txs, exs
      call cppTotXs(p, txs)
      call cppElaXs(p, exs)
      xs =max( txs - exs, 0.d0)
      end
