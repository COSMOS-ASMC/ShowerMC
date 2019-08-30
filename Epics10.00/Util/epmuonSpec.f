
!         unfinished


      real*8 function epmuonSpec00(p)
!           muon (+/-) flux * p**2.7 in /cm^2 s sr (GeV/c)**1.7
      implicit none
      real*8 p  ! input momentum in GeV/c
      
      real*8 lp, f
      lp = log10(p)
      if(lp .lt. -0.1938) then
         f = -1./(-0.1938 + 0.455938) *(lp + 0.455938)  - 3.0
         epmuonSpec00 = 10.**f
      elseif(lp .lt. 3.30107) then
         epmuonSpec00 = 0.052802*exp(- ((lp -1.197406)/0.444815)**2/2) +
     *       0.106497*exp(- ((lp -1.925)/0.7)**2/2 )
      else
         epmuonSpec00 = 0.01* (p/10.**3.30107)**(-1.)
      endif
      end
      real*8 function epmuonSpec0(p)
      implicit none
      real*8 p ! input. muon momentum in GeV/c
      real*8 epmuonSpec00
      epmuonSpec0 = epmuonSpec00(p)/p**2.7
      end
      real*8 function epmuonSpec(p, cosz)
      implicit none
      real*8 p !input
      real*8 cosz ! input

      real*8  ps
      
