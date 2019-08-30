!     ****************************************************************
!     *
!     * cmollerPath:  moller scattering prob. / r.l
!     * cmollerea: energy of survival and recoil electrons and angles
!     *
!
!
       subroutine cmollerPath(ein, w, prob, path)
       implicit none
#include  "Zmass.h"

        real*8 ein  ! input.  electron energy in GeV
        real*8 w    ! input.  minimum kinetic energy of recoil electron
                    !   to be treated in GeV.  (around 200 kev)
        real*8 prob ! output.  prob. per r%l
        real*8 path ! output.  sampled  path in r%l
!       --------------------------------


        real*8  em, beta2,  t0,  g, u
!        constm=.03*z/a*x0inkgpm
        real*8 constm/5.4859d0/
        real*8 cmollertx
!
        save constm


        g=ein/masele
        beta2=1.0d0 - 1./g**2
        t0=ein-masele
        if(t0 .gt. 0.) then
           em= w/t0
        else
           em=1000.
        endif
        if(em .ge. 0.5d0) then
           prob= 1.d-35
        else
!           ( .3*z/a) *x0ing = media.basearea*2; prob /r.l
          prob =
!     *      epmollertx(g, em)*masele/t0/beta2 * media.basearea*2.0
     *      cmollertx(g, em)*masele/t0/beta2 * constm
       endif
       call rndc(u)
       path = - log(u)/prob
       end
!      ************
       subroutine cmollerea(ein, w,  es, er, coss, cosr)
       implicit none
#include  "Zmass.h"
        real*8 ein  ! input.  electron energy in GeV
        real*8 w    ! input.  minimum kinetic energy of recoil electron
                    !   to be treated( in GeV).  (around 200 kev)
        real*8 es   ! output.  survival(higher) electron energy in GeV
        real*8 er   ! output.  recoiled electron energy in GeV
        real*8 coss ! output.  cos angle of the survival electron 
        real*8 cosr ! output.  cos angle of the recoiled //
!       --------------------------------

        real*8 g, em, t0, u, ep, ge, tr,  gr, gs 
        real*8 cmollerrf
        g = ein/masele
        t0=ein-masele
        if(t0 .gt. 0.) then
           em= w/t0
        else
           em=1000.
        endif
        if(em .ge. 0.50d0) then
           er=masele
           es=ein
           tr=0.
        else
!                   rejection method
!                *** until loop*** 
          do while (.true.)
             call rndc(u)
             ep=1.d0/ (  (1.d0-em*2.d0)*u/em + 2.d0 )
             ge = cmollerrf(g, ep)
             call rndc(u)
             if    ( u .lt. ge)
     *            goto 100
          enddo
 100      continue
          tr=ep*t0
          er= tr + masele
       endif
       es = ein - er + masele
       if(es .lt. masele) then
           es=masele
           er=max(ein-es+masele, masele)
           tr=er-masele
       endif
!                    angle part
       gr = er/masele
       gs = es/masele
       if(g .gt. 1.) then
          cosr =sqrt( (gr-1.0)*(g+1.0)/(gr+1.)/(g-1.0))
          coss = sqrt( (gs-1.0)*(g+1.)/(gs+1.)/(g-1.0))
       else
          cosr = 1.0
          coss=  1.0
       endif

       end

      real*8  function cmollerG(g,x)
      implicit none
      real*8  g  ! input gamma factor of electron
      real*8  x  ! input. w/T0.  w is cut-off kinetic energy.
!      ds/dx=  cmollerG / x**2 * 2pir^2mZ/beta^2/T0
      real*8 temp
      real*8 g1, g2, gsave
      data gsave/0.d0/
      save g1, g2, gsave
      if(g .ne. gsave) then
         gsave = g
         g1= (g-1.)/g 
         g2 = (2.0*g-1.0)/g**2
      endif
      temp= x/(1.0-x)
      cmollerG =  (g1*x)**2 +
     *      temp*(temp-g2) + 1.0
      end

      real*8 function cmollerrf(g,x)
      implicit none
      real*8 g  ! input.
      real*8 x  ! input.

      real*8 cmollerG
!        rejection function; G(g, x)/G(g,0.5d0)
      cmollerrf= cmollerG(g,x)/cmollerG(g,0.5d0)
      end

      real*8 function cmollertx(g, xm)
      implicit none
!      a part of the moller scattering total x-section 
!         if 2pi r^2 m/beta^2Z/T0  is
!      multiplied, you will get the total x-section for 
!      a charge Z atom.
!      
!         if 2pir^2/beta^2/T0 *Z Na/A = 0.300/beta^2* Z/A*m/T0 is 
!      multiplied, you will  get the total cross
!      section in 1/(g/cm^2), and its inverse is the  m.f.p in g/cm^2.
!
      real*8 g  ! input gamma factor
      real*8 xm  ! input.  w/T0.

      cmollertx = ((g-1.0)/g)**2 *(0.5-xm) + (1./xm-2.0) -
     *  (1./(1.0-xm)-2.0) + (2.0*g-1.0)/g**2*log(4.0*xm*(1.0-xm))
      end
