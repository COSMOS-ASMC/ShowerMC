      subroutine epSetMu(media, mu)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
       type(epmedia)::  media  ! input. media
       type(mubpn)::  mu     ! input. must be media.mu
      

      real*8 z, z3
      
!      z = media.Zeff
      z = media%Z
      z3 = z**(1.d0/3.d0)
!          const for muon sampling function.
!         
!       dv /v**(2-b) /(1+av**b)**2
!      where       
!         a= pa*Emu**(-0.40) + qa
!         b= pb*Emu**(qb) + 0.92
!    pa, pb, qb, ra are given as follows
!
      mu%pa = 49.232*z**(-0.05075)+0.7494-0.0188*z
      mu%ra = 124.9*z**(-0.0197)+1.2969
      mu%pb =  0.0871*z**(-0.0923)
      mu%qb = -0.3459*z**0.05615
!        used for rejection at brems sampling
      mu%Ak = 189.d0
      mu%Akm = mu%Ak*masmu/masele/z3
      mu%Akm2 = mu%Akm*sqrt(exp(1.d0))
      mu%PointLike = 0.22
!      mu.Shadow = media.Aeff**(0.1d0)
      mu%Shadow = media%A**(0.1d0)
      if(media%Z .gt.  10.) then
         mu%logf0 = log(mu%Akm*2./3.d0/z3)
      else
         mu%logf0 = log(mu%Akm)
      endif
      end


