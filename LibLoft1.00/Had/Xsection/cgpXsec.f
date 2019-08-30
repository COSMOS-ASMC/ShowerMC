!      real*8 ein, xs, a
!      ein =0.15
!      a = 14.6
!      do i = 1, 500
!         call cgpXsec(a, ein, xs)
!         ein = ein*10.0**0.025
!      enddo
!      end
!      xsection in mb for Photo production of hadrons at Neucleaus.
      subroutine  cgpXsec(a, energy, xs)
      implicit none
!       a: input. real*8  Target Mass No.
!     energy: input real*8.  gamma energy in GeV
!     xs:  output. real*8.   cross-section in mb
!
      real*8 a, energy, xs

      real*8 loge, pw

!
!           gp x-section,   xs in mb
      call cgpxs1(energy, xs)

      if( a .eq. 1.) then
!          nothing to do
      elseif(energy .lt. 2.) then    
         xs = xs*a*(1.02-.135*energy)
      elseif(energy .lt. 1.e6 ) then
!          before v7.50
!         xs = xs*a*0.75    ! this is close to A**0.91 dependence
!          pw = around 0.91 and decreasing to 0.82 
!             (by Bezrukov and Bugaev, Sov.J. of Nuc. phys.
!              33(5), may, 1981 
         loge = log (energy)
         pw = (-0.00062749*loge + 0.004126)*loge +0.24338+0.6666 
         xs = xs* a**pw
      else
!            extraplation to reach A**2/3 dependence
         pw = 0.18*(energy/1.e6)**(-0.07) +0.6666
         xs = xs* a**pw
      endif
      end
