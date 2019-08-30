!     *******************************************************
!     *                                                              *
!     *   compton scattering in  matter. E in GeV. Length in r.l    *
!     *                                                              *
!     *  ccomptPath:  get prob. for compton scattering   /r.l
!     *  ccomptea  samples energy of compton electron  and angles
!     *                                                   
!       This does not include decrease of the total cross-section
!       at low energies where atomic electron binding becomes
!       important
!     ****************************************************************
!
!  As compared to the old table method, this is ~10 % slower than.
!  However, no limitation at 20 keV.

!
      subroutine ccomptPath( Eg, p, path)
      implicit none
#include "Zmass.h"
      real*8 Eg             ! input. Gamma energy in GeV
      real*8 p              ! output. probability ( number of 
                 ! occurence ) of compton scattering  per r.l
      real*8 path    ! output. sampled path in r%l 


      real*8 u, g
!             tp* 3/8 /cconst= total cross section normalized
!             by thomson cross section.
!             cconst=3/8  * thomson * n0 * z/a * x0ing
!

      real*8 cconst/2.7429630/
      save cconst
!
      g  =Eg/masele
      if(g .lt. 0.1d0) then
!          p=(  (5.2*g-2.)*g +1. )*2.66666*media.basearea
          p=(  (5.2*g-2.)*g +1. )*2.66666*cconst
      else
          p=(  (1. - (g+1.)*2/g**2)*log(g*2+1.) + .5 + 4./g -
     *        1./(g*2+1.)**2/2 ) /g * cconst
      endif
      call rndc(u)
      path = - log(u)/p
      end
!
!
!     ***********
      subroutine ccomptea(Eg, Egout, Ee, cosg, cose)
      implicit none
#include  "Zmass.h"
!     ***********
      real*8 Eg  ! input. gamma energy in GeV.
      real*8 Egout  ! output. scattered gamma energy in GeV
      real*8 Ee     ! output. scattered electron energy in GeV
      real*8 cosg   ! output. cos of scattered gamma angle 
      real*8 cose   ! output. cos of scattered electron energy

      real*8 xmin, a1, a2, x, temp, u, sin2g, cos2e, g
!             x = Egout/Eg
      g  =Eg/masele
      xmin =  1.d0/( 1.d0 + 2*g)
      a1 = - log(xmin)
      a2 = (1.d0 -xmin*xmin) /2
      do while (.true.)
         call rndc(u)
         if( u .lt. a1/(a1+a2) ) then
!            sample form 1/x dx
            call rndc(u)
            x = xmin * exp( a1 * u )
         else
!              sample from x dx in (xmin~1)
            call ksampLin(1.d0, 0.d0, xmin, 1.d0, x)
         endif
!             rejection by (1- xsin^2/(1+x^2)
         temp = (1 - x)/x/g
         sin2g =  temp*(2.d0-temp)
         call rndc(u)
         if(u .lt. (1. - x*sin2g/(1+x*x))) goto 10
      enddo
 10   continue
      Egout = Eg*x
      Ee = Eg- Egout + masele
!        give cos of phton and electron
      cosg= 1. - temp
!           tan(el)=cot(gm/2)/(1+g);
!           cot(t/2)=  +-sqrt( (1+cos(t))/(1-cos(t)) ) so that
      cos2e=(1.d0-cosg) / (  1.d0-cosg +(1.d0+cosg)/(1.d0+g)**2 )
!           electron angle is always 0 to 90 deg.
      cose=min( max(sqrt(cos2e), 1.d-10), 0.9999999999d0)

      end
