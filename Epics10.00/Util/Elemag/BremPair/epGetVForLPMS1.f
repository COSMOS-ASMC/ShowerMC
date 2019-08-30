!   usage: echo  Fe | ./epGetV.out  etc
#include "ZepicsBD.h"
      implicit none
!
!
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"
#include "ZepTrackp.h"

      character*50  file
      character*16  basemedia
      integer nerg
      integer io, i, norm
      real*8  Ee, x, f, epBrgenex, xmax
      real*8  xmin, tprob, xx, fxx, Nc
      real*8  E1
      common /landuc/  s1, alogs1, sconst, x0g, conv2mb, effz
      real s1, alogs1, sconst, x0g, conv2mb, effz

      real*8 Tomb   ! to mb conversion
      parameter (Tomb = 1.d27/N0)
      real er, eps
      real v, E, s,  smigb, gzai, gmigdl, psimig

      data er/1.e-4/, eps/1.e-4/


      io = 10
      read(*,*)   basemedia

      file ="../../../Data/BaseM/"//basemedia

      open(io, file=file, action ='read')
      call epBPZpartH(media)
      call epBPgeneini(io, media)
  

      E = 3.
      do while (E < 1000.)
!      get v for ss=s =1 
!         ss=sqrt(sbrem2(v,E,s)) ! so sbrem2= 1
!      tmp=sconst*v
!      sbrem2=tmp/(1.-v)/e/gzai(s)   ! gzai(s)=1 for s=1
!        hence
!       1 = sconst*v/(1-v)/E      
!      then v= 
         v = 1./(  1. + sconst/E)
         write(0,*) E,v
!             verify
         s = smigb(v, E, s, er)
         write(*,*) E, v, s
         E = E*10**0.1
      enddo
      end
