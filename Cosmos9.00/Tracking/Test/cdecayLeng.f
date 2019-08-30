!     ******************************************************************
!     *                                                                *
!     * cdecayLeng: samples decay length of a given particle
!     *                                                                *
!        samples decay length of a particle.
!        
!        For a charged particle, we consider the change of the life
!     time due to the energy change by ionization loss.

      subroutine  cdecayLeng(aTrack, length)
      implicit none
#include  "Zglobalc.h"
#include  "Ztrack.h"


      type(track)::aTrack ! input%a  track of a decaying particle
      real*8  length  ! output. length to wthich the decay takes place (m).

      real*8 g, u,  gbeta, dedt, dedl, ctau, rho
      real*8 a, x, p, gmin, ctaumax, pmin, cvh2den
      type(fmom)::gb

      data gmin/200.d0/, ctaumax/50.d0/, pmin/10.d0/

      if(aTrack%p%fm%p(4) .le. aTrack%p%mass) then
         length = 0.
         return   ! *******
      endif

      call cgetctau(aTrack%p, ctau)
      if(ctau .eq. Infty) then
         length = Infty
      else
         g = aTrack%p%fm%p(4)/aTrack%p%mass
         if(g .gt. gmin) then
            call rndc(u)
            length = - ctau* g *log(u)
         else
            call cgetlf(aTrack%p, gb)
            gbeta = sqrt(gb%p(1)**2 + gb%p(2)**2 + gb%p(3)**2)
            if(aTrack%p%charge .eq. 0 .or. ctau .lt. ctaumax) then
               call rndc(u)
               length = - ctau * gbeta * log(u)
            else
               rho = cvh2den(aTrack%pos%height)
               call cdedxInAir(aTrack%p, rho, dedt) ! dedt in GeV/(kg/m^2)
               dedl= dedt*rho   !    GeV/m
               a = dedl/aTrack%p%fm%p(4) ! 1/m
               p =1.0d0/( a * ctau * g )
               if( p .gt. pmin) then
                  call rndc(u)
                  length = - ctau* gbeta *log(u)
               else   
                  call cdecayWEL(p, g, x)
                  length =(1.0d0 -  x/g)/a
               endif
            endif
         endif
      endif
      end
