      subroutine cSetMu(Zeffin, Aeffin)
      implicit none
#include "Zcmuint.h"
#include "Zmass.h"

      real*8 Zeffin ! input.
      real*8 Aeffin  ! input.

      real*8 z, z3

      muNLEmin = log10(muNEmin)
      muPrLEmin = log10(muPrEmin)
      muBrLEmin = log10(muBrEmin)
      
      z = Zeffin
      Zeff = z
      z3 = z**(1.d0/3.d0)
      Zeff3 = z3
!          const for muon sampling function.
!         
!       dv /v**(2-b) /(1+av**b)**2
!      where       
!         a= pa*Emu**(-0.40) + qa
!         b= pb*Emu**(qb) + 0.92
!    pa, pb, qb, ra are given as follows
!
      mupa = 49.232*z**(-0.05075)+0.7494-0.0188*z
      mura = 124.9*z**(-0.0197)+1.2969
      mupb =  0.0871*z**(-0.0923)
      muqb = -0.3459*z**0.05615
!        used for rejection at brems sampling
      muAk = 189.d0
      muAkm = muAk*masmu/masele/z3
      muAkm2 = muAkm*sqrt(exp(1.d0))
      muPointLike = 0.22
      muShadow = Aeffin**(0.1d0)
      if(Zeff .gt.  10.) then
         mulogf0 = log(muAkm*2./3.d0/z3)
      else
         mulogf0 = log(muAkm)
      endif
      end



