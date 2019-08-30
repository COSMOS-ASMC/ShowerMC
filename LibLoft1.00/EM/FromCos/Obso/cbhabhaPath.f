!     ****************************************************************
!     *
!     * cbhabhaPath:  Bhabha scattering prob. / r.l
!     * cbhabhaea: energy of survival posititon and recoil 
!                   electrons and angles
!
       subroutine cbhabhaPath(ein, w, prob, path)
       implicit none
#include  "Zmass.h"
        real*8 ein  ! input.  positron energy in GeV
        real*8 w    ! input.  minimum kinetic energy of recoil electron
                    !   to be treated (in GeV).  (around 200 kev)
        real*8 prob ! output.  prob. per r%l
        real*8 path ! output.  sampled  path in r%l
!       --------------------------------
        real*8  em,  t0,  g, u

        real*8 cbhabhatx
!            constm=.03*z/a*x0inkgpm2
        real*8 constm
        save constm
        data  constm/5.4859d0/

!

        g=ein/masele
        t0=ein-masele
        if(t0 .gt. 0.) then
           em= w/t0
        else
           em=1000.
        endif
        if(em .ge. 1.0d0) then
           prob= 1.d-35
        else
!           ( .3*z/a) *x0ing = media.basearea*2; prob /r.l
          prob =
!     *      epbhabhatx(g, em)*masele/t0 * media.basearea*2.0
     *      cbhabhatx(g, em)*masele/t0 * constm
       endif
       call rndc(u)
       path = - log(u)/prob
       end
!      ************
       subroutine cbhabhaea(ein, w,  es, er, coss, cosr)
       implicit none
#include  "Zmass.h"
        real*8 ein  ! input.  positron energy in GeV
        real*8 w    ! input.  minimum kinetic energy of recoil electron
                    !   to be treated( in GeV).  (around 200 kev)
        real*8 es   ! output.  survival positron energy in GeV
        real*8 er   ! output.  recoiled electron energy in GeV
        real*8 coss ! output.  cos angle of the survival positron
        real*8 cosr ! output.  cos angle of the recoiled //
!
!            this is the same as Moller scat. case. except fof
!           rejection function. and range of em, ep
!       --------------------------------

        real*8 g, em, t0, u, ep, ge, tr,  gr, gs 
        real*8 cbhabharf

        g = ein/masele
        t0=ein-masele
        if(t0 .gt. 0.) then
           em= w/t0
        else
           em=1000.
        endif
        if(em .ge. 1.d0) then
           er=masele
           es=ein
           tr=0.
        else
!                   rejection method
!                *** until loop*** 
          do while (.true.)
             call rndc(u)
             ep=1.d0/ (  (1.0d0-em)*u/em + 1.0d0 )
             ge=cbhabharf(w, em, g, ep)
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

      real*8  function cbhabhaG(g, x)
      implicit none
      real*8  g  ! input gamma factor of electron
      real*8  x  ! input. w/T0.  w is cut-off kinetic energy.
!      ds/dx=  cbhabhaG / x**2 * 2pir^2mZ/T0
      real*8 temp
      real*8 c1, c2, c3, c4, y, gsave, beta2
      data gsave/0.d0/

      save  c1, c2, c3, c4, y, gsave, beta2

      if(g .ne. gsave) then
         gsave = g
         beta2 = 1.d0 - 1.d0/g**2
         y=1.d0/(1.0d0 + g)
         temp = 1.d0-2.0d0*y
         c4= temp**3
         c3 = c4 + temp**2
         c2 = temp*(3.0d0 + y**2)
         c1 = 2.d0 - y**2
      endif
      cbhabhaG = (( (c4*x -c3)*x + c2)*x -c1) *x +1./beta2
      end

      real*8 function cbhabharf(w, em, g, x)
      implicit none
#include "Zmass.h"
      real*8 w  ! input.  cut-off energy in GeV.
      real*8 em ! input.  w/T0
      real*8 g  ! input.
      real*8 x  ! input.

      real*8 cbhabhaG
!        rejection function; G(g, x)/G(g,x)_max 

      real*8 g1,  gm
      real*8 brp
      parameter (brp=masele*0.999)

      gm = cbhabhaG(g, em)

      if(w .gt. brp)  then
!          check if g1 > gm         
         g1  = cbhabhaG(g, 1.d0)
         gm = max(g1, gm) 
      endif
      cbhabharf= cbhabhaG(g,x)/gm
      end

      real*8 function cbhabhatx(g, em)
      implicit none
!      a part of the Bhabha scattering total x-section 
!         if Z* 2pi r^2 m/T0  is
!      multiplied, you will get the total x-section for 
!      a charge Z atom.
!      
!         if  2pir^2 m/T0 *Z Na/A = 0.300*Z/A*m/T0 is 
!      multiplied, you will  get the total cross
!      section in 1/(g/cm^2), and its inverse is the  m.f.p in g/cm^2.
!
      real*8 g  ! input gamma factor
      real*8 em  ! input.  w/T0.

      real*8 beta2, y, c1, c2, c3, c4, temp

      beta2 = 1.d0 - 1.d0/g**2
      y=1.d0/(1.0d0 + g)
      temp = 1.d0-2.0d0*y
      c4= temp**3
      c3 = c4 + temp**2
      c2 = temp*(3.0d0 + y**2)
      c1 = 2.d0 - y**2
      
      cbhabhatx = c1*log(em) +c2*(1.d0-em) - c3*(1.0d0-em**2)/2.0d0
     *        + c4*(1.0d0 -em**3)/3.0d0 
     *        + (1.d0-em)/beta2/em
      end
