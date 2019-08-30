!        test  cnkg and cmnkg
!
!      real*8 cnkg, s, r 
!      s = 0.4
!      do while (s .lt. 1.8)
!         r = 0.01
!         do while (r .lt. 50.)
!            write(*, *) sngl(r), sngl(cnkg(s, r))
!            r = r * 10.**0.1
!         enddo
!         write(*, *)
!         s = s + 0.2
!      enddo
!      end
!      
!      *****************************************************************
!      *                                                               *
!      * cnkg: lateral distribution of electrons by Greisen's formula  *
!      *                                                               *
!      *****************************************************************
!
!  It is not recomended to use cnkg as it is, because it give too wide
!  lateral distribution. You must use about 1/2 Moliere Unit to get
!  rather good distribution.  
!         

!
!    (time: if don't use polynomial approximation, it will be 1.9 times)
!
       real*8 function cnkg(s, r)
       implicit none
       real*8 s ! input.  age of shower (one dimensional age)
       real*8 r ! input.  distance from the center of shower. in moliere unit
!      function value.  normalized number of electrons /(m.u)**2
                  !  (integral (nkg)*2pir*dr from 0 to infinity = 1 )
!                     if s <=0 or  s > 2.25, rho=0 is put

!
       real*8 cnkgcs
!
       if(s .gt. 0.  .and. s .lt. 2.25) then
           cnkg = r**(s-2.) * (1.+r)**(s-4.5) * cnkgcs(s)
       else
           cnkg = 0.
       endif
       end
       real*8 function cnkgcs(s)
       implicit none
!      parameter (pi=3.1415, tpii=1./pi/2.)
!          csnkg=gma(4.5 - s)
!    *     /gma(s)/gma(4.5-s*2)  * tpii
!          polynom approx.
!             better than .05 % except s<.15 or s>2.0
       real*8 s
       real*8 temp
           temp = (((((((-.4315466e-01*s+0.2827905 )*s-.5464292 )*s+
     *     .1768150 )*s+0.1253462 )*s+0.1909647 )*s+0.2234282 )*s-
     *     .1260477e-01)*s+0.7873029e-03
           cnkgcs = temp/s
       end
!      *****************************************************************
!      *                                                               *
!      * cmnkg: lateral distribution of electrons by pair electrons    *
!      *       created by m.c                                          *
!      *                                                               *
!      *****************************************************************
!
!

!
! -- output --
!
!  *** note ***
!   if s < 0.6 or  s > 2, nkg is used.
!   for 10000 call's in the mnkg region, 181 msec is needed on m380.
!   if polynomial approx for r0,rx,cs are not used, it becomes 400
!   msec.
!
      real*8 function cmnkg(s, r)
      implicit none
      real*8 s  ! input. 1 dim age.
      real*8 r  ! input.  distance from the center of shower. in moliere unit
!
!     function value:  normalized number of electrons /(m.u)**2
!           (integral (rho)*2pir*dr from 0 to infinity = 1 )
      
       real*8 a, b, r0, cs
!
       call cmnkgcs(s, a, b, r0, cs)
       cmnkg = (r/r0)**(s-a) * (1.+r/r0)**(s-b) *cs
       end
!      ************************
       subroutine cmnkgcs(s, a, b, r0, cs)
       implicit none
       real*8 s, a, b, r0, cs
!
       real*8 cnkgcs, rx

       if(s .lt. .6  .or. s .gt. 2.) then
!            use nkg
           cs = cnkgcs(s)
           a = 2.
           b = 4.5
           r0 = 1.
       else
!          tmp=(-3.309*log10(s) - .285)*log10(s) - .249
!          r0=10.**tmp
!             .1 %  approx of r0 at s=(.4-2.0).
           r0=(((((((0.2272372 *s-2.547489 )*s+ 12.33055 )*s-
     *       33.50755)*s+ 55.33528  )*s-55.46729  )*s+ 30.93687 )*s-
     *       7.405076  )*s+0.6610533
!              this give too large lateral for small s
!          rx=.29* s**1.415/r0
!          rx for m.c lateral (s=.4 to 2.) good within .2 %)
!         rx=((((((-.3179712    *s+ 2.959941    )*s-11.47338    )*s+
!    *    24.40712)*s-30.39885    )*s+ 22.81418    )*s-9.357027    )*s+
!    *    1.880529
!
!             corrected one. rx=.29*s**(1.40+.1/s**4)/r0 (s=.5 to 2.)
           rx=(((((-.3310270 *s+ 2.953910)*s-10.35767)*s+
     *     19.40102 )*s-19.71657 )*s+ 10.79047 )*s-2.225853
           a=-.21*s + 2.146
           b= (2*s+2.-a) + (s+2.-a)/rx
!
!          cs=gma(b-s)/gma(s+2.-a)/gma(a+b-2.-s*2)*tpii / r0**2
!           cs*s for m.c lateral (good at s=.4 to 2 within .3 %)
          cs=((((((-.9450091    *s+ 8.993343    )*s-35.02486    )*s+
     *    71.67270)*s-82.55823    )*s+ 56.46587    )*s-21.05954    )*s
     *    + 3.595105
          cs = cs/s
       endif
       end
