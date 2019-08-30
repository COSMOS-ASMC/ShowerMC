#include "ZcosmosBD.h"
      program testxs
      use modpdgXs, only : csOfstu
      implicit none      
!  #include "ZcosmosExt.h"
#include "Zmass.h"
!
!       sigma tot, inel, ela test for pp
!
      real*8  p,  Ek, xsT, XsI, XsE, roots, s


      Ek = 0.01
      write(*,*) '# Ek(GeV) p roots  St Sin Sel Sela'
      do while (Ek .lt. 1.e9)
         p = sqrt( (Ek+ maspic)**2 - maspic**2)
         s = csOfstu(p, maspic, masp)
         roots =sqrt(s)
         call cpimpTotXs(p, xsT)
         call cpimpInelaXs(p,  xsI)
         call cpimpElaXs(p, xsE)
         write(*,'(1p,7g13.4)') p, Ek, roots, xsT, xsI, xsE
         Ek = Ek*10.**0.02
      enddo

      end





