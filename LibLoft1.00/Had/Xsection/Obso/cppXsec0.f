!   collection of inelastic cross-section at low energy.
!       you should not use this one  e > 2.e6 GeV for pbar-p
!       and for all others e>  1 TeV. p.d formula is too bad
!       in such a region.
!
!        cppXsec0   p   p      n p
!        cpbarpXsec0  p_b p      n_b  n    n_b p   p_b n
!        cpiMinuspXsec0  pi- p
!        cpiPluspXsec0  pi+ p
!        ckMinuspXsec0  k-  p      k0_b n
!        ckPluspXsec0  k+  p k+ n  k- n   k0_b  p   k0 n
!        cpiMAirReacX  pi- air reaction xs
!        cpiMAirAbsX   pi- air absorption xs
!        cpiPAirReacX   pi+ air reaction xs
!        cpiPAirAbsX    pi+ air absorption xs
!
!       ***********************************************************
!       *
!       * p_p inelastic cross-section
!       *
!       ****************** tested 88.08.11 *****************k.k ***
!
!   E: input. real*8.  E kinetic energy in GeV
!              for E<300 MeV XS=1.e-5
!              for E>100 GeV Xs by Particle Data Book Formula
!
!   XS:output. Inelastic cross-section in mb
!
        subroutine cppXsec0(e, xs)
        implicit none
        real*8 e, xs
!
        real*8 ee
        integer i, icon, l
        real*8 lp
        real*8 ctotppX, celappX

        real*8 ea(21), xsa(21)
!           ea in GeV
      data (ea    (i),i=   1,  21)/
     1 0.301, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700,
     2 0.800, 0.900, 1.000, 1.500, 2.000, 5.000, 7.000,10.000,20.000,
     3 40.000,50.000, 100.001/
      data (xsa   (i),i=   1,  21)/
     1 0.000, 0.500, 1.000, 2.500, 5.000, 7.350,11.000,15.000,18.500,
     222.000,23.000,23.500,25.000,25.500,28.000,29.000,29.500,30.000,
     330.000,29.500,31.000/
      if(e .le. 0.)  then
         xs = 0.

      elseif(e .lt. .30101d0) then
         xs=1.d-5
      elseif(e .lt. 100.) then
         ee = e
         call kfrge(ea, 1, 21, ee, l, icon)
         xs=(xsa(l)-xsa(l-1))/(ea(l)-ea(l-1))*(ee-ea(l-1)) +
     *        xsa(l-1)
      else
         lp = log(e)            ! e = p      
         xs = ctotppX(e, lp) - celappX(e, lp)
      endif
      end
      real*8 function ctotppX(p, logp)
      implicit none
        real*8 p, logp
        real*8 a, b, n,  c, d
        real*8 cXsec
        data a/48./, b/0./, n/0./, c/.522/, d/-4.51/
        ctotppX = cXsec(p, logp, a,  b, n, c, d)
      end
      real*8 function celappX(p, logp)
      implicit none
        real*8 p, logp
        real*8 a, b, n,  c, d
        real*8 cXsec
        data a/11.9/, b/26.9/, n/-1.21/, c/0.169/, d/-1.85/
        celappX = cXsec(p, logp, a, b, n,  c, d)
      end
!     **********************************
!     This particle data book formulation 
!
      real*8 function cXsec(p, logp, a, b, n, c, d)
!        elastic/inelastic  high energy  cross-secion formula.
!        a + bp**n + c*log(p)**2 + d*log(p) 
!        p is momentum in GeV/c
      implicit none
      real*8 p, logp, a, b, n, c, d
!
      cXsec = (c*logp + d)*logp + b*p**n + a
      end
!       
!      
!       ***********************************************************
!       *
!       * p_b p inelastic cross-section
!       *
!       ******************                 *****************k.k ***
!
!   e: input.  e kinetic energy in GeV
!
!   xs:output. inelastic cross-section in mb
!
        subroutine cpbarpXsec0(e, xs)
        implicit none
        real*8 e, xs
        real*8 lp
        real*8 ctotpbpX, celapbpX


!               due to the formula given in ptcl data 88. p121

        if(e .lt. 100.) then
           xs =  (e/100)**(-.55)*3 +32
        else
           lp = log(e)   !  e = p here.
           xs =  ctotpbpX(e, lp) - celapbpX(e, lp)
        endif
        end
      real*8 function ctotpbpX(p, logp)
      implicit none
        real*8 p, logp
        real*8 a, b, n,  c, d
        real*8 cXsec
        data a/38.4/, b/77.6/, n/-0.64/, c/0.26/,  d/-1.2/
        ctotpbpX = cXsec(p, logp, a,  b, n, c, d)
      end
      real*8 function celapbpX(p, logp)
      implicit none
        real*8 p, logp
        real*8 a, b, n,  c, d
        real*8 cXsec
        data a/10.2/, b/52.7/, n/-1.16/, c/0.125/, d/-1.28/
        celapbpX = cXsec(p, logp, a,  b, n, c, d)
      end

!       ***********************************************************
!       *
!       * pi- p inelastic cross-section
!       *
!       ****************** tested 88.08.11 *****************k.k ***
!
!   e: input.  pi- k.e     in GeV
!              for e<10mev  xs=1.e-5
!
!   xs:output. inelastic cross-section in mb
!
        subroutine cpiMinuspXsec0(e, xs)
        implicit none
        real*8 e, xs
        integer i, l, icon
        real*8 lp
        real*8 ctotpimpX, celapimpX

        real*8 ea(25), xsa(25), ee

      data (ea    (i),i=   1,  25)/
     1 0.010, 0.020, 0.030, 0.060, 0.100, 0.200, 0.300, 0.400, 0.450,
     2 0.500, 0.550, 0.670, 0.870, 1.000, 1.200, 1.350, 1.500, 2.000,
     3 2.800, 3.000, 5.000, 7.000,10.000,50.000, 100.001 /
      data (xsa   (i),i=   1,  25)/
     1 0.000, 4.000, 6.000,10.000,15.000,15.300,16.500,18.000,19.500,
     222.500,30.000,26.000,33.500,32.000,25.000,23.500,25.000,27.000,
     326.000,24.500,23.500,22.500,22.500,21.500,21.000/
        if(e .le. 10.d-3) then
           xs=1.d-5
        elseif(e .lt. 100.) then
           ee= e
           call kfrge(ea, 1, 25, ee, l, icon)
           xs=(xsa(l)-xsa(l-1))/(ea(l)-ea(l-1))*(ee-ea(l-1)) +
     *         xsa(l-1)
        else
           lp = log(e)   !  e = p here.
           xs = ctotpimpX(e, lp) - celapimpX(e, lp)
        endif
      end
      real*8 function ctotpimpX(p, logp)
      implicit none

        real*8 p, logp
        real*8 a, b, n,  c, d
        real*8 cXsec
        data a/33./, b/14./, n/-1.36/, c/0.456/, d/-4.03/
        ctotpimpX = cXsec(p, logp, a,  b, n, c, d)
      end
      real*8 function celapimpX(p, logp)
      implicit none
        real*8 p, logp
        real*8 a, b, n,  c, d
        real*8 cXsec
        data a/1.76/, b/11.2/, n/-0.64/, c/0.043/, d/0./
        celapimpX = cXsec(p, logp, a,  b, n, c, d)
      end

!       ***********************************************************
!       *
!       * pi+ p inelastic cross-section
!       *
!       ****************** tested 88.08.11 *****************k.k ***
!
!   e: input.  pi+ k.e     in GeV
!              for e<100mev xs=1.e-5
!              for e>100gev  p.d
!
!   xs:output. inelastic cross-section in mb
!
        subroutine cpiPluspXsec0(e, xs)
        implicit none 
        real*8 e, xs
        
        real*8 ea(18), xsa(18), ee
        integer i, icon, l
        real*8 lp, ctotpippX, celapippX
!
      data (ea    (i),i=   1,  18)/
     1 0.100, 0.300, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100,
     2 1.200, 1.300, 1.500, 5.000,10.000,30.000,50.000,70.000, 100.001/
      data (xsa   (i),i=   1,  18)/
     1 1.d-5, 4.500, 7.500, 8.000,10.500,12.000,13.500,16.000,18.000,
     2 19.500,19.800,20.000,20.200,19.500,19.200,19.200,19.500,20.000/
!
        if(e .lt. 100.d-3) then
           xs=1.d-5
        elseif(e .lt. 100.) then
           ee= e
           call kfrge(ea, 1, 18, ee, l, icon)
           xs=(xsa(l)-xsa(l-1))/(ea(l)-ea(l-1))*(ee-ea(l-1)) +
     *         xsa(l-1)
        else
           lp = log(e)   !  e = p here.
           xs = ctotpippX(e, lp) - celapippX(e, lp)
        endif
      end
      real*8 function ctotpippX(p, logp)
      implicit none
        real*8 p, logp
        real*8 a, b, n,  c, d
        real*8 cXsec
        data a/16.4/, b/19.3/, n/-0.42/, c/0.19/, d/0./
        ctotpippX = cXsec(p, logp, a,  b, n, c, d)
      end
      real*8 function celapippX(p, logp)
      implicit none
        real*8 p, logp
        real*8 a, b, n,  c, d
        real*8 cXsec
        data a/0./, b/11.4/, n/-0.4/, c/0.079/, d/0./
        celapippX = cXsec(p, logp, a,  b, n, c, d)
      end

!       ***********************************************************
!       *
!       * k-  p inelastic cross-section
!       *
!       ****************** tested 88.08.11 *****************k.k ***
!
!   e: input.  k- k.e   GeV
!              for e<40.e-3 MeV xs=1.e-5
!              for e>100gev xs  P.D formula
!
!   xs:output. inelastic cross-section in mb
!
        subroutine ckMinuspXsec0(e, xs)
        implicit none
        real*8 e, xs
        real*8 ee
        integer i, icon, l

        real*8 lp
        real*8 ctotkmpX, celakmpX
        real*8 ea(23), xsa(23)

!
      data (ea    (i),i=   1,  23)/
     1 0.040, 0.070, 0.080, 0.085, 0.090, 0.100, 0.120, 0.130, 0.150,
     2 0.170, 0.200, 0.250, 0.300, 0.400, 0.500, 0.600, 0.640, 1.500,
     3 5.000,10.000,20.000,50.000, 100.0001/
      data (xsa   (i),i=   1,  23)/
     1 1.d-5,30.000,38.000,39.800,40.000,40.200,40.000,39.000,35.000,
     230.000,22.500,17.500,18.000,20.000,22.500,26.000,27.000,25.000,
     321.000,20.000,19.000,18.500,18.0000/
!
       if(e .le. 40.e-3)then
           xs=1.e-5
       elseif(e .lt. 100.) then
           ee = e
           call kfrge(ea, 1, 23, ee, l, icon)
           xs=(xsa(l)-xsa(l-1))/(ea(l)-ea(l-1))*(ee-ea(l-1)) +
     *         xsa(l-1)
       else
           lp = log(e)   !  e = p here.
           xs = ctotkmpX(e, lp) - celakmpX(e, lp)
       endif
      end
      real*8 function ctotkmpX(p, logp)
      implicit none
        real*8 p, logp
        real*8 a, b, n,  c, d
        real*8 cXsec
        data a/32.1/, b/0./, n/0./, c/0.66/, d/-5.6/
        ctotkmpX = cXsec(p, logp, a,  b, n, c, d)
      end
      real*8 function celakmpX(p, logp)
      implicit none
        real*8 p, logp
        real*8 a, b, n,  c, d
        real*8 cXsec
        data a/7.3/, b/0./, n/0./, c/0.29/, d/-2.4/
        celakmpX = cXsec(p, logp, a,  b, n, c, d)
      end

!       ***********************************************************
!       *
!       * k+  p inelastic cross-section
!       *
!       ****************** tested 88.08.11 *****************k.k ***
!
!   e: input.  k+ k.e    GeV
!              for e<100.e-3  xs=1.e-5
!              for e>100.gev xs by p.d
!
!   xs:output. inelastic cross-section in mb
!
        subroutine ckPluspXsec0(e, xs)
        implicit none
        real*8 e, xs
        integer i, icon, l
        real*8 ea(13), xsa(13), ee
        real*8 lp
        real*8 ctotkppX, celakppX


      data (ea    (i),i=   1,  13)/
     1 0.100, 0.200, 0.300, 0.500, 1.000, 2.000, 3.000, 4.000, 5.000,
     210.000,20.000,50.000, 100.0001 /
      data (xsa   (i),i=   1,  13)/
     1 1.d-5, 2.500, 5.000, 7.500,10.000,12.000,13.000,13.500,14.000,
     2 15.000,15.300,15.800,16.200/

        if(e .le. 100.e-3)then
           xs=1.d-5
        elseif(e .lt. 100.) then
           ee = e
           call kfrge(ea, 1, 13, ee, l, icon)
           xs=(xsa(l)-xsa(l-1))/(ea(l)-ea(l-1))*(ee-ea(l-1)) +
     *         xsa(l-1)
        else
           lp = log(e)   !  e = p here.
           xs = ctotkppX(e, lp) - celakppX(e, lp)
        endif
      end
      real*8 function ctotkppX(p, logp)
      implicit none
        real*8 p, logp
        real*8 a, b, n,  c, d
        real*8 cXsec
        data a/18.1/, b/0./, n/0./, c/0.26/, d/-1.0/
        ctotkppX = cXsec(p, logp, a,  b, n, c, d)
      end
      real*8 function celakppX(p, logp)
      implicit none
        real*8 p, logp
        real*8 a, b, n,  c, d
        real*8 cXsec
        data a/5.0/, b/8.1/, n/-1.8/, c/0.16/, d/-1.3/
        celakppX = cXsec(p, logp, a,  b, n, c, d)
      end

!       ***********************************************************
!       *
!       * pi- air reaction cross-section
!       * (total - elastic)
!       *
!       ****************** tested 88.08.11 *****************k.k ***
!
!   e: input.  pi- k.e   GeV
!              for e>500.mev xs at e=500 mev is given
!
!   xs:output. reaction  cross-section in mb
!
        subroutine cpiMAirReacX(e, xs)
        implicit none
        real*8 e, xs
        integer i, icon, l
        real*8 ea(23), xsa(23), ee


!            MeV
      data (ea    (i),i=   1,  23)/
     1   0.,  10.,  20.,  30.,  40.,  50.,  70.,  80., 100., 120., 140.,
     2 160., 180., 200., 240., 260., 280., 300., 340., 360., 400., 450.,
     3 500.0001/
      data (xsa   (i),i=   1,  23)/
     1   0.,  92., 138., 172., 205., 230., 356., 414., 452., 460., 462.,
     2 463., 461., 460., 416., 391., 359., 328., 276., 241., 224., 218.,
     3 207./
           ee=min(500.d0, e*1.d3)
           call kfrge(ea, 1, 23, ee, l, icon)
           xs=(xsa(l)-xsa(l-1))/(ea(l)-ea(l-1))*(ee-ea(l-1)) +
     *         xsa(l-1)
        end

!       ***********************************************************
!       *
!       * pi+ air reaction cross-section
!       * (total - elastic)
!       *
!       ****************** tested 88.08.11 *****************k.k ***
!
!   e: input.  pi+ k.e  GeV
!              for e>500.mev xs at e=500 mev is given
!
!   xs:output. reaction  cross-section in mb
!
        subroutine cpiPAirReacX(e, xs)
        implicit none
        real*8 e, xs

        integer i, icon, l
        real*8 ea(24), xsa(24), ee
!           mev
      data (ea    (i),i=   1,  24)/
     1   0.,  10.,  20.,  30.,  40.,  50.,  70.,  80.,  90., 100., 120.,
     2 140., 160., 180., 200., 240., 260., 280., 300., 340., 360., 400.,
     3 450., 500.0001/
      data (xsa   (i),i=   1,  24)/
     1   0.,  57.,  94., 126., 155., 195., 316., 379., 408., 416., 437.,
     2 443., 448., 443., 437., 393., 370., 345., 328., 253., 230., 207.,
     3 190., 190./
           ee=min(500.d0, e*1.d3)
           call kfrge(ea, 1, 24, ee, l, icon)
           xs=(xsa(l)-xsa(l-1))/(ea(l)-ea(l-1))*(ee-ea(l-1)) +
     *         xsa(l-1)
        end
!       ***********************************************************
!       *
!       * pi- air absorption x-section
!       *
!       ****************** tested 88.08.11 *****************k.k ***
!
!   e: input.  pi- k.e   GeV
!              for e>500.mev xs at e=500 mev is given
!
!   xs:output. reaction  cross-section in mb
!
        subroutine cpiMAirAbsX(e, xs)
        implicit none
        real*8 e, xs
        integer i, icon, l
        real*8 ea(13), xsa(13), ee
!           mev
      data (ea    (i),i=   1,  13)/
     1   0.,  10.,  20.,  30.,  40.,  50.,  70., 100., 140., 180., 200.,
     2 300., 500.001/
      data (xsa   (i),i=   1,  13)/
     1   0.,  91., 136., 159., 163., 172., 184., 193., 195., 193., 184.,
     2 115.,   0./
           ee=min(500.d0, e*1.d3)
           call kfrge(ea, 1, 13, ee, l, icon)
           xs=(xsa(l)-xsa(l-1))/(ea(l)-ea(l-1))*(ee-ea(l-1)) +
     *         xsa(l-1)
        end
!       ***********************************************************
!       *
!       * pi+ air absorption x-section
!       *
!       ****************** tested 88.08.11 *****************k.k ***
!
!   e: input.  pi+ k.e   tev
!              for e>500.mev xs at e=500 mev is given
!
!   xs:output. reaction  cross-section in mb
!
        subroutine cpiPAirAbsX(e, xs)
        implicit none
        real *8 e, xs
!
        integer i, icon, l
        real * 8 ea(13), xsa(13), ee
!           mev
      data (ea    (i),i=   1,  13)/
     1   0.,  10.,  20.,  30.,  40.,  50.,  70., 100., 140., 180., 200.,
     2 300., 500.001/
      data (xsa   (i),i=   1,  13)/
     1   0.,  47.,  92., 115., 129., 140., 161., 178., 184., 184., 182.,
     2 115.,   0./
           ee=min(500.d0, e*1.d3)
           call kfrge(ea, 1, 13, ee, l, icon)
           xs=(xsa(l)-xsa(l-1))/(ea(l)-ea(l-1))*(ee-ea(l-1)) +
     *         xsa(l-1)
        end
