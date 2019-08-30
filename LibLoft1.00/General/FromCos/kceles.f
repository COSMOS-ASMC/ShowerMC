!   ********************** kceles *********************************
!             this file containes programs for celestial
!             coordinate transformation routines.
!  1) for subroutines which are related to the vernal equinox
!     position (implicitly or explicitly), the accuracy is kept within
!     the order of 1/100 degrees which is enough for our purpose.
!     for other subroutines, double precision accuracy is kept.
!  3) the input/output parameters and function values are all
!     in double precision.
!  4) many of subroutines use a block common $kcele where some
!     fundamental constants and common variables are defined.
!     $kcele is defined as a separate module and included by
!     include statement to each subroutine.
!  5) many subroutines require that kcelei should have been called
!     before hand.  kcelei should be called once for a particular
!     observation place.(*)
!  6) for subroutines which are related to the a.s array
!     coordinate,  kadthi should have been called.(*)
!  7) the vernal equinox employed is a new one recommended by
!     iau.
!  8) ut1 is approximated by utc=lst-dtgmt (local standard time
!     - difference from the greenwich time).  the difference
!     ut1 - utc (nearly =) dut1 is not considered. it is at most
!     1 sec.
!  9) to improve the accuracy of the calculation for the apparent
!     position of stars,  we have to introduce the equation of
!     equinoxes, dut1, nutaion, proper motion of stars, parallax, polar
!     motion etc. (@)
! 10) one exceptional case for the correction of parallax is
!     the case of the apparent position of the moon.  it is
!     considered in kmoont.
!
!  @) in conversion routines from b1950.0 to j2000 (kb50j2) or
!     j2000 to b1950.0 (kj2b50), neglecting the proper motion and
!     e-term results in an error of 5/100 degrees.
!  *) two initializations are managed by jreadi in tips,
!     so that the tips user need not worry about that.
!
!        call kcelei(lat, lon, dtgmt, height)  inform place once.
!                                 latitude, longitude (deg), gmt
!                                 difference(hour), and height (m)
!        call kadthi(astox)      inform a.s array x-axis angle
!                                 measured anticlockwise from south.
!
! -----------------------------------------------------------------
!  kcosd  : function to compute cos with argument in degree.
!  kadbp : to get arrival direction of a.s by plane front approx.
!  kadth : convert arrival direction in array system to
!          direction cos. in horizontal system.
!  kadthi: to initialize 'arrival direction to horizontal coord.
!          transformation' routine
!  kb50j2: convert (dec,ra)b1950.0 into those for j2000
!  kcelei: to initialize celestial coordinate conversion routines
!  kceleq: inquire current latitude, longitude, dtgmt, and height
!  kctoq : ecliptic to equatorial coordinate transformation
!  kdcmjd: decompose mjd into y,m,d,h,m,s
!  kdhtoh: for given declination and time from meridian,
!          get d.c in horizontal system at observation place.
!  kdifva: get vector angle difference
!  kdtoa : convert direction cosines into polar angles
!  kdtoe : convert equatorial coord. given in direction cos.
!          into equatorial declination and r.a.
!          (this can be used for latitude, longitude relation
!           in other coord. system such as ecliptic one).
!  kdtoh : convert horizontal coord. given in direction cos.
!          into height and azimuth.
!  kdzth2: from dec. and zenith, give time from meridian.
!  kdztoh: from dec. and zenith. angle, give time  from
!          meridian. (see for diff. from kdzth2, program)
!  kedtgd: equatorial angle in direction cos to galactic angle
!          in direction cos
!  keqtog:  equatorial angle to galactic angle conversion
!  ketod : inverse of kdtoe
!  ketoh : equatorial to horizontal coordinate transformation
!  kgcrc : compute coordinate of observation place in geocentric
!          rectangular system.
!  kgcttc: conversion from the geocentric to topocentric
!          coordinate (topocenter --> origin is at the
!                                observation place)
!  kgdted:  inverse of kedtgd.
!  kgtoeq:  inverse of keqtog.
!  khtad : inverse of kadth.
!  khtod : inverse of kdtoh
!  khtoe : horizontal to equatorial coordinate transformation
!  kjtmjd: julian day to modified julian day conversion
!  kjxjy : position at j'x' is converted to that at j'y'
!  kj2b50: convert (dec,ra)j2000.0 into those for b1950.0
!  kj2j90: inverse of kj90j2.
!  kj2tox: convert equatorial vector in j2000 into those at
!          time x which has been specified by pij matrix
!          obtained in kpmtrx
!  kj90j2: convert dec, ra in 1990 japanese ephemeris into j2000
!  kmjd  : get modified julian day
!  kmjdst: convert mjd into the mean local siderial time
!  kmjdtj: mjd to julian day conversion
!  kmjdym: functions the same as kdcmdj
!  kmobec: get mean obliquity of the ecliptic plane
!         (inclination angle of ecliptic plane)
!  kmoon : get ecliptic angle of the moon in geocentric coord.
!  kmoont: get topocentric equatorial coord. of the moon.
!  kmtoj2: mjd to time measured from j2000.0 conversion
!  knormv: normalize vector and get vector length
!  koangl: get opening angle between two (dec, ra)
!  komega: get solid angle around a source given by komeg0
!          for a given point in celestial  coord.
!  komeg0: initialization for komega and komeg1
!  komeg1: get theta between a source given by komeg0 and
!          a given point in celestial coord.
!  kpmtrx: compute precession matrix from a given mjd.
!  kqtoc : equatorial to ecliptic  coordinate transformation
!  ksided: to compute local mean siderial time
!  ksidet: to compute local mean siderial time.
!  kside0: to compute siderial time of greenwich at 0hut of a
!          given day
!  ksun  : compute the sun position in ecliptic coord.
!  ksuneq: compute the sun position in equatorial coord.
!  ktu   : to compute elapsed days form 1900/1/0  0h ut
!          to a given day
!  kxtoj2: from observed equatorial vector to j2000 vector
!  kvtoa : get polar angles from unnormalized 3-d vector
!  ksind  : function to get sin with an argument in degree
!        ------------------------------------------------------------
!
         subroutine kcelei(tlat, tlon, dtgmt, height)
!
!   tlat: input. latitude of the place in deg. (+for northern hemisph..
!   tlon: input. longitude of the place in deg. (+ for east, - for west
!   dtgmt: input. time difference from gmt of the place. in hour.
!          + for earlier time
!  height: input. height of the place above the sea level (m).
!
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
!
             tlons=tlon
             tlats=tlat
             coslat=cos(tlats*Torad)
             sinlat=sin(tlats*Torad)
!
             dtgmts=dtgmt
             heighs=height
!                get geocentric rectangular coordinate of the
!                place. (ug, vg, wg)
             call kgcrc(tlats, tlons, heighs, ug, vg, wg)
!                correction to geocenter (in japan)
             ug=ug-136.
             vg=vg+521.
             wg=wg+681.
             return
       entry kceleq(tla, tlo, dt, h)
!            inquire current const for kcelei
           tlo=tlons
           tla=tlats
           dt=dtgmts
           h=heighs
           return
        end
!           compute coordinate of a given place
!           in the geocentric rectangular
!           system.
        subroutine kgcrc(fai, al, h, u, v, w)
!           fai: input. latitude  in deg.
!           al: input. longitude in deg.
!           h: input. height in m
!           u, v, w: output. coordinate of the place in the geocentric
!                  rectangular system.
!                   (in m).
           implicit real*8 (a-h, o-z)
           real*8 ksind, kcosd
           real*8 n
!                besselian ellipsoid
          data ae/6377.397155d03/
          data e2/0.006674372230614d0/
!
          sinf=ksind(fai)
          n=ae/sqrt(1.d0 - e2*sinf**2)
          cosf=kcosd(fai)
          u=(n+h)*cosf * kcosd(al)
          v=(n+h)*cosf * ksind(al)
          w=(n*(1.d0-e2)+h) * sinf
         end
!        *************************************************************
!        *
!        *  ksidet: compute mean siderial time
!        *          at a given longitude and time
!        *
!        *************************************************************
!
!   /usage/
!              call ksidet(year, month, day, time,  st)
!
!  year:   input.  integer*4.  say 87
! month:   input.  //          say 12
!   day:   input.  //          say 23
!  time:   input.        local time of the place. in hour
!                 12h.25m.36s.21==>12.+25./60.+36.21/3600.
!    st:   output. siderial time in degree.
!                 at longitude tlong (common), at time time (hour)
!
       subroutine ksidet(year, month, day, time, st)
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
           integer*4  year, month, day
!
           call ktu(year, month, day, ed)
           call kside0(ed, st0)
!              elapsed time from  0ut  in deg.
           et= (time-dtgmts)*15.d0
           std= st0 + tlons + sidcor* et
!               mod(360)
           st=  mod(std, 360.d0)
       end
!         ---------------- get modified julian day ----------
!      year:input  like 87
!      month:input like 11
!      day:  input like 23     these are integers
!      time: input like (hh + mm/60.d0+ sec/3600.d0)
!      mjd:  output.  modified julian days
!
!  ****   kcelei should have been called before hand.
!
       subroutine kmjd(year, month, day, time, mjd)
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
           integer*4  year, month, day
           real*8 mjd
           call ktu(year, month, day, ed)
!                  get modified julian day
           mjd=ed + (time-dtgmts)/24.d0 + 15019.5d0
       end
!        *************************************************************
!        *
!        *  ksided: compute mean siderial time
!        *          at a given longitude and time
!        *
!        *************************************************************
!
!   /usage/
!              call ksided(time, st0, st)
!
!  time:   input.  local time of the place. in hour
!                 12h.25m.36s.21==>12.+25./60.d0+36.21/3600.d0
!   st0:   input.  siderial time of greenwich in degree.
!                 at 0 hour.
!    st:   output. siderial time in degree.
!                 at longitude tlong (common), at time time (hour)
!
       subroutine ksided(time, st0, st)
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
!              elapsed time from 0ut  in deg.
           et= (time-dtgmts)*15.d0
           std= st0 + tlons + sidcor* et
!               mod(360)
           st=  mod(std, 360.d0)
       end
!      -------------------------------------------------------------
!              compute siderial time of greenwich at ut 0h of the
!              day.
!   ed: input.  elapsed day since 1900.1.0 noon of ut.
!  st0: output.
       subroutine kside0(ed,  st0)
         implicit real*8 (a-h, o-z)
!                consts.. to get siderial time in degree.
           parameter  (
     1      c1= (45.836d0/3600.d0+38.d0/60.d0 + 6.0d0)*15.d0,
     2      c2= 8640184.542d0/3600.d0 *15.d0,
     3      c3= 0.0929d0/3600.d0*15.d0 )
!
!
           tu=ed/36525.d0
           sd=( c3*tu + c2)* tu + c1
           st0=mod(sd, 360.d0)
       end
!          get elapsed day from 1900/1/0 noon ut at greenwich
!          if 15019.5 is added.  this is equivalent to the
!          modified julian day of the day at 0:0:0.
!       year, month, day:  input. integer*4.  year must be such as
!          87
!        ed:  output.
       subroutine ktu(iyear, month, day, ed)
         implicit real*8 (a-h, o-z)
           integer*4  year, month, day
           year=mod(iyear, 1900)
           if(month .le. 2) then
               iy=year-1
               im=month+12
           else
               iy=year
               im=month
           endif
           ed= 365*iy + 30*im + day + int(3*(im+1)/5) + int(iy/4)
     1         - 33.5d0
       end
!      ----------------------------------------------------------
!             convert horizontal to equatorial coord.
!      /usage/ call khtoe(st, hx, hy, hz, ex, ey, ez)
!          st: input.  siderial time in deg.
!  (hx,hy,hz): input.  direction cos of horizontal angle
!  (ex,ey,ez): output. //               equatorial angle
!        ** note **
!            to convert these angle into theta and fai
!            use call kdtoe(ex, ey, ez, teta, fai)
!
!           convert to equatorial to horizontal
!      /usage/ call ketoh(st, ex, ey, ez, hx, hy, hz)
!    (ex,  .) is input and (hx,..) is output
!            to convert these angle into theta and fai
!            use call khtoa(hx, hy, hz, teta, fai)
!
       subroutine khtoe(st, hx, hy, hz, ex, ey, ez)
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
           coss=cos(st*Torad)
           sins=sin(st*Torad)
           ex=hx*coss*sinlat - hy * sins + hz*coss*coslat
           ey=hx*sins*sinlat + hy*coss + hz * sins*coslat
           ez=-hx*coslat  + hz*sinlat
       end
       subroutine ketoh(st, ex, ey, ez, hx, hy, hz)
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
           coss=cos(st*Torad)
           sins=sin(st*Torad)
           hx=ex*sinlat*coss + ey*sinlat*sins - ez*coslat
           hy=-ex*sins + ey*coss
           hz= ex*coslat*coss  + ey*coslat*sins + ez * sinlat
       end
       subroutine ketod( delta, alfa,  ex, ey, ez)
         implicit real*8 (a-h, o-z)
!           convert r.a and declination into direction cos
!       this can be used to get other direction cosines
!       from latitude and longitude (say ecliptic coordinate)
#include  "Zkcele.h"
          cosa=cos(alfa*Torad)
          sina=sin(alfa*Torad)
          cosd=cos(delta*Torad)
          sind=sin(delta*Torad)
!
          ex=cosd*cosa
          ey=cosd*sina
          ez=sind
          return
       entry khtod(h, a, hx, hy, hz)
          cosh=cos(h*Torad)
          sinh=sin(h*Torad)
          cosa=cos(a*Torad)
          sina=sin(a*Torad)
          hx=cosh*cosa
          hy=-cosh*sina
          hz=sinh
       end
       subroutine kdtoe(ex, ey, ez, delta, alfa)
         implicit real*8 (a-h, o-z)
!       this can be used to get other latitude and longitude
!       from direction cosines.     (say ecliptic coordinate)
            call kdtoa(ex, ey, ez, t, f)
            delta=90.d0-t
            alfa=mod(f+360.d0, 360.d0)
       end
       subroutine kdtoh(hx, hy, hz, teta, fai)
         implicit real*8 (a-h, o-z)
            call kdtoa(hx, hy, hz, t, f)
            teta=90.d0-t
            fai=mod(360.d0-f, 360.d0)
       end
!              compute direction cos of arrival direction
!          of a shower in horizontal coordinate.
!       /usage/   call kadthi(astox) init.
!                 call kadth(ax, ay, az, hx,hy,hz)
!     astox: input. angle (in deg) between south to the x-axis of the
!            array
!    (ax,ay,ax): input. direction cos of arrival direction in array
!                coordinate  (x-y plane is assumed to be horizontal)
!    (hx, hy, hz):output.  direction cos of arrival direction in
!                 horizontal coord.
       subroutine kadthi(astox)
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
          cossx=cos(astox*Torad)
          sinsx=sin(astox*Torad)
       end
       subroutine kadth(ax, ay, az, hx, hy, hz)
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
          hx=ax*cossx - ay*sinsx
          hy=ax*sinsx + ay*cossx
          hz=az
       end
!         inverse of kadth
       subroutine khtad(hx, hy, hz, ax, ay, az)
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
          ax=hx*cossx + hy*sinsx
          ay=-hx*sinsx + hy*cossx
          az=hz
       end
!      *************************************************************
!      *  kdhtoh: for given declination and time from culmination,
!      *          get d.c in horizontal system at observation
!      *          given by kcelei
!      *************************************************************
!      del: input. declination in degree of source
!        h: input. time (in hour) from the meridian (+ -)
!  w1,2,3:  output.  direction cos in horizontal system
!
       subroutine kdhtoh(del, h, w1, w2, w3)
!
       implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
           data delsv/-100.d0/
           save delsv, cosd, sind
!
           if(del .ne. delsv) then
               cosd=cos(del*Torad)
               sind=sin(del*Torad)
               delsv=del
           endif
           fai=h*tofai *Torad
           sinfai=sin(fai)
           cosfai=cos(fai)
           w1=sinlat*cosd*cosfai - coslat*sind
           w2=cosd*sinfai
           w3=coslat*cosd*cosfai + sinlat*sind
       end
!      ************
!        from declination and zenith angle give time from
!        culmination
!      del:input. declination in degree
!       w3:input. 3rd component of the direction cos of zenith
!        h: output. time in hour from the culmination
!     icon: output =0  ---> h obtained
!                   1  ---> such zenith angle cannot happen
       subroutine kdztoh(del, w3, h, icon)
!      ************
!
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
           data delsv/-100.d0/
           save delsv, cosd, sind

           if(del .ne. delsv) then
               cosd=cos(del*Torad)
               sind=sin(del*Torad)
               delsv=del
           endif
           if(cosd .eq. 0. .or. coslat .eq. 0.) then
              if(abs(w3-sinlat*sind) .le. 1.d-5) then
                 h=12.d0
                 icon=0
              else
                 icon=1
              endif
           else
              cosx=( w3 - sind*sinlat )/cosd/coslat
              if(abs(cosx) .le. 1.d0) then
                 icon=0
                 fai=acos(cosx)*Todeg
                 h=fai*toh
              elseif(abs(w3- sinlat*sind) .le. 1.d-5) then
                 icon=0
                 h=12.d0
              else
                 icon=1
              endif
           endif
       end
!        from declination and zenith angle give time from
!        culmination
!      del:input. declination in degree
!       w3:input. 3rd component of the direction cos of zenith
!        h: output. time in hour from the culmination
!     icon: output =0  ---> h obtained
!                   1  ---> max zenith angle should be smaller
!                           than the w3 and w3 is changed to
!                           the max value. h corresponds to
!                           such value .
!                   2-----> min zenith angle realized is larger
!                           than w3 so that no h is obtained
!       w3: output.   see above (icon=1)
!      1------------------------0     cos (zenith)
!           !        !
!           !        !
!           !min &   max zenith observable
!        !      !        !      <---w3 position
!        2      0        1      <---icon
       subroutine kdzth2(del, w3, h, icon)
!      ************
!
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
           data delsv/-100.d0/
           save delsv, cosd, sind

           if(del .ne. delsv) then
               cosd=cos(del*Torad)
               sind=sin(del*Torad)
               delsv=del
           endif
           if(cosd .eq. 0.d0 .or. coslat .eq. 0.d0) then
              if(abs(w3-sinlat*sind) .le. 1.d-5) then
                 h=12.d0
                 icon=0
              elseif(w3 .lt. sinlat*sind) then
                 w3=sinlat*sind
                 h=12.d0
                 icon=1
              else
                 icon=2
              endif
           else
              cosx=( w3 - sind*sinlat )/cosd/coslat
              if(abs(cosx) .le. 1.d0) then
                 icon=0
                 fai=acos(cosx)*Todeg
                 h=fai*toh
              elseif(abs(w3- sinlat*sind) .le. 1.d-5) then
                 icon=0
                 h=12.d0
              elseif(w3 .lt. sinlat*sind) then
                 h=12.d0
                 w3=sinlat*sind
                 icon=1
              else
                 icon=2
              endif
           endif
       end
!      ---------- equatorial angle to galactic angle -----------
!  /usage/  call keqtog(dec, ra, glat, glon)
!
!   dec:  input.  declination in degree
!    ra:  input.  right ascension in degree
!  glat:  output. galactic latitude in degree
!  glon:  output. galactic longitude in degree
!
       subroutine keqtog(dec, ra, glat, glon)
         implicit real*8 (a-h, o-z)
!
             call ketod(dec, ra, ex, ey, ez)
             call kedtgd(ex, ey, ez, gx, gy, gz)
             call kdtoa(gx, gy, gz, glat, glon)
             glat=90.d0-glat
             glon=mod(glon+360.d0, 360.d0)
        end
!  --------------  inverse of above -------------
        subroutine kgtoeq(glat, glon, dec, ra)
         implicit real*8 (a-h, o-z)
             call ketod(glat, glon, gx, gy, gz)
             call kgdted(gx, gy, gz, ex, ey, ez)
             call kdtoe(ex, ey, ez, dec, ra)
        end
!         ------- galactic direction cos to equatorial direction cos
!    /usage/  call kgdted(gx, gy, gz, ex, ey, ez)
!
!   gx, gy, gz: input.  galactic direction cos
!   ex, ey, ez: output. equatorial //
!
       subroutine kgdted(gx, gy, gz, ex, ey, ez)
!
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
           logical first/.true./
           parameter (
     +     t=(12.d0+49d0/60.d0)*15.d0*Torad, f=62.6d0*Torad)
!
           if(first) then
               first=.false.
               cos33=cos(33.d0*Torad)
               sin33=sin(33.d0*Torad)
               sint=sin(t)
               cost=cos(t)
               sinf=sin(f)
               cosf=cos(f)
               a11=-sint
               a12=-cost*cosf
               a13=cost*sinf
               a21=cost
               a22=-sint*cosf
               a23=sint*sinf
               a32=sinf
               a33=cosf
               b11=-sint
               b12=cost
               b21=-cosf*cost
               b22=-cosf*sint
               b23=sinf
               b31=sinf*cost
               b32=sinf*sint
               b33=cosf
           endif
           gxp=cos33*gx + sin33*gy
           gyp=-sin33*gx + cos33*gy
           gzp=gz
           ex=a11*gxp +a12*gyp + a13*gzp
           ey=a21*gxp +a22*gyp + a23*gzp
           ez=         a32*gyp + a33*gzp
           return
! ---------------------  e to g ----------------------
      entry kedtgd(ex, ey, ez, gx, gy, gz)
           if(first) then
               first=.false.
               cos33=cos(33.d0*Torad)
               sin33=sin(33.d0*Torad)
               sint=sin(t)
               cost=cos(t)
               sinf=sin(f)
               cosf=cos(f)
               a11=-sint
               a12=-cost*cosf
               a13=cost*sinf
               a21=cost
               a22=-sint*cosf
               a23=sint*sinf
               a32=sinf
               a33=cosf
               b11=-sint
               b12=cost
               b21=-cosf*cost
               b22=-cosf*sint
               b23=sinf
               b31=sinf*cost
               b32=sinf*sint
               b33=cosf
           endif
           gxp=b11*ex + b12*ey
           gyp=b21*ex + b22*ey + b23*ez
           gzp=b31*ex + b32*ey + b33*ez
           gx= cos33*gxp -sin33*gyp
           gy= sin33*gxp + cos33*gyp
           gz=gzp
           return
       end
!     real*8  mjd, time
!     integer*4  y, m, d
!     mjd=47096.3114004629d0
!     call kmjdym(mjd, y, m, d, time)
!     write(*,*) ' y,m,d=',y, m, d, ' time=',time
!     end
      subroutine kmjdym(mjd, y, m, d, time)
!           not used in janzos.   see kdcmjd
         implicit real*8 (a-h, o-z)
          real*8  mjd,  jd
          integer*4 y, m, d, a, b, c, e, f, g, h
          time=(mjd- int(mjd))*24.d0
          jd=mjd+ 2400000.5d0
          a=int(jd)+68569
          b=int( a/36524.25)
          c=a - int(36524.25d0*b+0.75d0)
          e=int ( (c+1)/365.2425d0)
          f= c - int(365.25d0*e)+31
          g=int (f/30.59d0)
          d=f -int (30.59d0*g)+ 0.5d0 + (jd -int(jd))
          h=int(g/11.d0)
          m=g-12*h+2
          y=100*(b-49)+ e + h
          if(d .eq. 32) then
             d=1
             y=y+1
             m=1
          endif
       end
!         this include kdcmjd
!                 and  kmjdst
!
!     real*8  mjd, time, st
!     integer*4  y, m, d
!       mjd=47161.0000004629d0
!       call kdcmjd(mjd, y, m, d, time)
!       write(*,*) ' y,m,d=',y, m, d, ' time=',time
!       call kcelei(36., 135., 9.)
!       call kmjdst(mjd, st)
!       write(*,*) ' st=',st
!     end
!     ****************************************************************
!     *
!     *  kdcmjd: decompose mjd into year, month, day, and hour (ut)
!     *
!     ****** tested 89.03.13******************************************
!      this one is based on honda.
!      same result is obtained  as kmjdym made by k.k
!
!  /usage/ call kdcmjd(mjd, iy, im, id, ihr, imn, sec
!
!       mjd: input.  real*8.   modified julian day
!        iy: output.  integer*4.  year.  say, 1989
!        im: //                   month.  say, 11
!        id: //                   day.    say, 23
!       ihr: //       //          hour    say, 12
!       imn: //       //           min.   say  15
!      sec:  //                   sec.   say  12.4567
!
!
      subroutine kdcmjd(mjd,iy,im,id,ihr,imn,sec)
         implicit real*8 (a-h, o-z)
!     input  mjd :  modified julian day (fractional) real*8
!     output  iy : year  im month id day of month
        real*8  mjd
!
        time =(mjd- int(mjd))*24.0d0
        ihr  = int(time)
        time =(time - ihr)*60.d0
        imn  = int(time)
        sec  =(time - imn)*60.d0
!
        jd= int(mjd + 2400001.0d0)
        l = jd + 68569
        n = 4*l/146097
        l = l - (146097*n+3)/4
        iy = 4000*(l+1)/1461001
        l = l - 1461*iy/4+31
        im = 80*l/2447
        id = l - 2447*im/80
        l = im/11
        im = im + 2 - 12*l
        iy = 100*(n-49) + iy + l
      end
!     ****************************************************************
!     *
!     *  kmjdst: mjd to siderial time conversion. (local)
!     *          this is based old equinox. for j2000 ephemeris
!     *          use the new one with the same name
!     ****** tested 89.03.13******************************************
!
!  /usage/ call kmjdst(mjd, st)
!
!        kcelei must have been called beforehand.
!
!       mjd: input.  real*8.   modified julian day
!        st: output.           siderial time at specified longitude
!                             (normalized to 1.)
!     subroutine kmjdst(mjd, st)
!        implicit real*8 (a-h, o-z)
!          real*8 mjd
!          data siddy0 /0.671262d0/
!
!          du = (mjd - 40000.0d0)
!          tu = du/36525.0d0
!
!          sidsol= sidcor +( sidcr2 + sidcr3*tu )*tu
!          sidday= siddy0 + sidsol*du + tlons/360.0d0
!
!          sidt = sidday-int(sidday)
!          if(sidday.lt.0.0d0) then
!              sidt =1.0d0 + sidt
!          endif
!          st= sidt
!      end
       subroutine kmjdst(mjd, st)
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
          real*8 mjd
          parameter (c1=(50.54841d0/3600.d0 + 41.d0/60.+18.d0),
     *    c2=8640184.812866d0/3600.d0,
     *    c3=0.0931047d0/3600.d0, c4=-0.0000062d0/3600.d0)
!             get j2000 time
          call kmtoj2(mjd, t)
          am=   ((c4*t+ c3)*t + c2)*t + c1
!            approximate ut1 (dut1  is < 1 sec) if mjd is correct
          ut1=( mjd - int(mjd) )*24.d0
          st=( ut1+ 12.d0 + am )*15.d0 + tlons
          st=mod(st, 360.d0)
          if(st .lt. 0.d0) then
             st=st+360.d0
          endif
       end
       subroutine kmjdtj(mjd, jd)
           real*8 mjd, jd
           jd=mjd+2400000.5d0
       end
       subroutine kjtmjd(jd, mjd)
           real*8 jd, mjd
           mjd=jd-2400000.5d0
       end
!                 test of j2000 conversion
!          declination  is changed about 0.05 deg
!          right ascension is changed about 0.12 deg
!
!       real*8 pij(3,3), mjd, ex, ey, ez, ex2, ey2, ez2,
!     * d, r, do, ao
!c          mjd of about 90.04.20.
!       mjd=48000.d0
!       call kpmtrx(mjd, pij)
!       write(*,*) (pij(i,1), i=1,  3 )
!       write(*,*) (pij(i,2), i=1,  3 )
!       write(*,*) (pij(i,3), i=1,  3 )
!       call ketod(35.d0, 100.d0, ex, ey, ez)
!       write(*,*)  ' ex , ey , ez =', ex , ey , ez
!       call kxtoj2(pij, ex, ey, ez, ex2, ey2, ez2)
!       write(*,*)  ' ex2, ey2, ez2=', ex2, ey2, ez2
!       call kj2tox(pij, ex2, ey2, ez2, ex, ey, ez)
!       write(*,*)  ' ex , ey , ez =', ex , ey , ez
!       call kdtoe(ex2, ey2, ez2, d, r)
!       write(*,*) ' d, r=',d, r
!       do=0.
!       ao=0.
!       call ketod(do, ao, ex, ey, ez)
!       call kxtoj2(pij, ex, ey, ez, ex2, ey2, ez2)
!       call kdtoe(ex2, ey2, ez2, d, r)
!       write(*,*) ' do, ao=',do, ao
!       write(*,*) ' d, r=',d, r
!       do=50.
!       ao=0.
!       call ketod(do, ao, ex, ey, ez)
!       call kxtoj2(pij, ex, ey, ez, ex2, ey2, ez2)
!       call kdtoe(ex2, ey2, ez2, d, r)
!       write(*,*) ' do, ao=',do, ao
!       write(*,*) ' d, r=',d, r
!       do=25.
!       ao=0.
!       call ketod(do, ao, ex, ey, ez)
!       call kxtoj2(pij, ex, ey, ez, ex2, ey2, ez2)
!       call kdtoe(ex2, ey2, ez2, d, r)
!       write(*,*) ' do, ao=',do, ao
!       write(*,*) ' d, r=',d, r
!       end
!      ***********************************************************
!      *
!      * kpmtrx:  compute precession matrix defined in p.423-425
!      *          of 'japanese ephemeris 1990'
!      *
!      ***********************************************************
!     /usage/ call kpmtrx(mjd, pij)
!   mjd:  input.   modified julian day (real*8)
!   pij:  output.  pij(3,3)   precession matrix
!                  pij(i,j)=nij
!       suppose  directional cosines are
!       given in the rectangular equatorial
!       coordinate based on the  equinox at the instance  of
!       the observation. (the observed directional cosines are
!       such ones provided that the siderial time contains the
!       effect of the precession).
!       the matrix computed here can be used to obtain the
!       directional cosines at j2000 by
!
!              d2  =  pji dobsv
! **** this is computed only at first event of each run *****
!      because change is very small.
!
      subroutine kpmtrx(mjd, pij)
       implicit real*8 (a-h,o-z)
       real*8 mjd
       dimension pij(3,3)
       parameter (pi=3.14159265358979324d0,
     *            Torad=pi/3600.d0/180.d0)
!       days from 2000y01m1.5d (julian day,2451545.0)
!        in unit of 100 julian years
       call kmtoj2(mjd, t)
       t2 = t*t
       t3 = t*t2
!
       zeta = 2306.2181d0*t + 0.30188d0*t2 + 0.017998d0*t3
       za =   2306.2181d0*t + 1.09468d0*t2 + 0.018203d0*t3
       teta = 2004.3109d0*t - 0.42665d0*t2 - 0.041833d0*t3
!
       zeta = zeta* Torad
       za = za* Torad
       teta = teta*Torad
!
       pij(1,1) = cos(zeta)*cos(teta)*cos(za) -sin(zeta)*sin(za)
       pij(1,2) =-sin(zeta)*cos(teta)*cos(za) -cos(zeta)*sin(za)
       pij(1,3) =                           -sin(teta)*cos(za)
       pij(2,1) = cos(zeta)*cos(teta)*sin(za) +sin(zeta)*cos(za)
       pij(2,2) =-sin(zeta)*cos(teta)*sin(za) +cos(zeta)*cos(za)
       pij(2,3) =                           -sin(teta)*sin(za)
       pij(3,1) = cos(zeta)*sin(teta)
       pij(3,2) =-sin(zeta)*sin(teta)
       pij(3,3) = cos(teta)
      end
      subroutine kmtoj2(mjd, t)
!       days from 2000y01m1.5d (julian day,2451545.0)
!        in unit of 100 julian years
         real*8  t, mjd
         t  = (mjd - 51544.5d0)/36525d0
      end
!     *****************************************************************
!     *
!     *  kj2tox: convert equatorial vector in j2000 into those at
!     *          time x which has been specified by pij matrix
!     *          obtained in kpmtrx
!     *
!     *****************************************************************
!
       subroutine kj2tox(pij, ex2, ey2, ez2, ex, ey, ez)
!
!    pij(3,3): input.   precession matrix obtained in kpmtrx
!    ex2,...ez2: input. directional cosines in j2000 ephemeris
!    ex,...ez: output.  directional cosines at x
!
       implicit real*8 (a-h,o-z)
       dimension pij(3,3)
!
       ex= pij(1,1)*ex2 + pij(1,2)*ey2 + pij(1,3)*ez2
       ey= pij(2,1)*ex2 + pij(2,2)*ey2 + pij(2,3)*ez2
       ez= pij(3,1)*ex2 + pij(3,2)*ey2 + pij(3,3)*ez2
      end
!     *****************************************************************
!     *
!     *  kxtoj2: convert equatorial vector at x which has been specifie
!     *          by pij matrix computed in kpmtrx into the one
!     *          in j2000 ephemeris
!     *
!     *****************************************************************
!
       subroutine kxtoj2(pij, ex, ey, ez, ex2, ey2, ez2)
!
!    pij(3,3): input.   precession matrix obtained in kpmtrx
!    ex,...ez: input.   directional cosines at x
!    ex2,...ez2: output.directional cosines in j2000 ephemeris
!
       implicit real*8 (a-h,o-z)
       dimension pij(3,3)
!
       ex2= pij(1,1)*ex + pij(2,1)*ey + pij(3,1)*ez
       ey2= pij(1,2)*ex + pij(2,2)*ey + pij(3,2)*ez
       ez2= pij(1,3)*ex + pij(2,3)*ey + pij(3,3)*ez
      end
!                    test kmoon
!         implicit real*8 (a-h, o-z)
!       real*8 mjd/48234.00068287d0/
!       call kmoon(mjd, elat, elon, rm)
!       write(*,*) ' elat=',elat, ' elon=',elon, ' rm=',rm
!       mjd=33250.0d0
!       call kctoq(mjd, 0.2535301d0, 1.5274972d0,
!    *  0.0260904d0,  ex, ey, ez)
!       write(*,*) ex, ey,ez
!       end
!       ************************************************************
!       *
!       * kmoon: compute the moon position in the ecliptic coord.
!       *        (apparent geocentric ecliptic coordinate)
!       *   (accuracy is better than 1/100 degree)
!       ************************************************************
! /usage/  call kmoon(mjd,  elat, elon, rmoon)
!
!   mjd: input.  real*8.  modified julian day of the place, time
!   elat: output.         apparent ecliptic latitude of  the moon
!                         (degree)
!   elon: output.         apparent ecliptic longitude of them moon
!                         (degree)
! rmoon:  output.         distance from the center of the earth
!                         to that of the moon. (gravitation center)
!                         (in m).
!                         the constant of the sin parallax is
!                         sinpi=ae/r where ae= the equatorial
!                         radius of the earth, r the distance
!                         to the moon. (between the gravitational
!                         center). ae=6378.140 km
!
!     **note** this coordinate may be converted into that in the
!          apparent geocentric  equatorial coord, and then
!          converted to the apparent topocentric equatorial coord.
!
        subroutine kmoon(mjd, elat, elon, rmoon)
          implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
          real*8 mjd, ksind, kcosd
!
          real*8 c(62), d(62), e(62)
          real*8 f(46), g(46), h(46)
!
          data  c/
     1               1.2740d0,  0.6583d0,  0.2136d0,  0.1856d0,
     2    0.1143d0,  0.0588d0,  0.0572d0,  0.0533d0,  0.0459d0,
     3    0.0410d0,  0.0348d0,  0.0305d0,  0.0153d0,  0.0125d0,
     4    0.0110d0,  0.0107d0,  0.0100d0,  0.0085d0,  0.0079d0,
     5    0.0068d0,  0.0052d0,  0.0050d0,  0.0048d0,  0.0040d0,
     6    0.0040d0,  0.0040d0,  0.0039d0,  0.0037d0,  0.0027d0,
     7    0.0026d0,  0.0024d0,  0.0024d0,  0.0022d0,  0.0021d0,
     8    0.0021d0,  0.0021d0,  0.0020d0,  0.0018d0,  0.0016d0,
     9    0.0012d0,  0.0011d0,  0.0009d0,  0.0008d0,  0.0008d0,
     a    0.0007d0,  0.0007d0,  0.0007d0,  0.0006d0,  0.0005d0,
     b    0.0005d0,  0.0005d0,  0.0005d0,  0.0004d0,  0.0004d0,
     c    0.0004d0,  0.0004d0,  0.0003d0,  0.0003d0,  0.0003d0,
     d    0.0003d0,  0.0003d0,  0.0003d0/
!
          data  d/
     1               107.248d0,  51.668d0, 317.831d0, 176.531d0,
     2    292.463d0,  86.16 d0, 103.78 d0,  30.58 d0, 124.86 d0,
     3    342.38 d0,  25.83 d0, 155.45 d0, 240.79 d0, 271.38 d0,
     4    226.45 d0,  55.58 d0, 296.75 d0,  34.5  d0, 290.7  d0,
     5    228.2  d0, 133.1  d0, 202.4  d0,  68.6  d0,  34.1  d0,
     6      9.5  d0,  93.8  d0, 103.3  d0,  65.1  d0, 321.3  d0,
     7    174.8  d0,  82.7  d0,   4.7  d0, 121.4  d0, 134.4  d0,
     8    173.1  d0, 100.3  d0, 248.6  d0,  98.1  d0, 344.1  d0,
     9     52.1  d0, 250.3  d0,  81.0  d0, 207.0  d0,  31.0  d0,
     a    346.0  d0, 294.0  d0,  90.0  d0, 237.0  d0,  82.0  d0,
     b    276.0  d0,  73.0  d0, 112.0  d0, 116.0  d0,  25.0  d0,
     c    181.0  d0,  18.0  d0,  60.0  d0,  13.0  d0,  13.0  d0,
     d    152.0  d0, 317.0  d0, 348.0  d0/
!
           data  e/
     1           -4133.3536d0,  8905.3422d0,  9543.9773d0,  359.9905d0,
     2 9664.0404d0,   638.635d0, -3773.363d0,13677.331d0, -8545.352d0,
     3 4411.998d0,  4452.671d0,  5131.979d0,  758.698d0, 14436.029d0,
     4 -4892.052d0,-13038.696d0,14315.966d0,-8266.71d0, -4493.34d0,
     5 9265.33d0,   319.32d0,   4812.66d0,    -19.34 d0,13317.34d0,
     6 18449.32d0,  -1.33 d0,  17810.68d0,   5410.62d0, 9183.99d0,
     7 -13797.39d0,  998.63d0, 9224.66d0,  -8185.36d0,  9903.97d0,
     8  719.98 d0, -3413.37d0, -19.34d0, 4013.29d0, 18569.38d0,
     9 -12678.71d0, 19208.02d0, - 8586.0d0, 14037.3d0,-7906.7d0,
     a 4052.0 d0,   -4853.3d0,  278.6 d0,   1118.7d0, 22582.7d0,
     b 19088.0d0,  -17450.7d0, 5091.3d0,   -398.7d0, -120.1d0,
     c 9584.7  d0, 720.d0,     -3814.0d0, -3494.7d0,18089.3d0,
     d 5492.0d0,   -40.7d0,    23221.3d0/
!
         data f/
     1            0.2806d0, 0.2777d0, 0.1732d0, 0.0554d0,
     2  0.0463d0, 0.0326d0, 0.0172d0, 0.0093d0, 0.0088d0,
     3  0.0082d0, 0.0043d0, 0.0042d0, 0.0034d0, 0.0025d0,
     4  0.0022d0, 0.0021d0, 0.0019d0, 0.0018d0, 0.0018d0,
     5  0.0017d0, 0.0016d0, 0.0015d0, 0.0015d0, 0.0014d0,
     6  0.0013d0, 0.0013d0, 0.0012d0, 0.0012d0, 0.0011d0,
     7  0.0010d0, 0.0008d0, 0.0008d0, 0.0007d0, 0.0006d0,
     8  0.0006d0, 0.0005d0, 0.0005d0, 0.0004d0, 0.0004d0,
     9  0.0004d0, 0.0004d0, 0.0004d0, 0.0003d0, 0.0003d0,
     a  0.0003d0, 0.0003d0/
!
         data g/
     1           215.147d0,  77.316d0,   4.563d0, 308.98 d0,
     2 343.48d0, 287.90d0,  194.06 d0,  25.6  d0,  98.4  d0,
     3   1.1 d0, 322.4 d0,  266.8  d0, 188.0  d0, 312.5  d0,
     4 291.4 d0, 340.0 d0,  218.6  d0, 291.8  d0,  52.8  d0,
     5 168.7 d0,  73.8 d0,  262.1  d0,  31.7  d0, 260.8  d0,
     6 239.7 d0,  30.4 d0,  304.9  d0,  12.4  d0, 173.0  d0,
     7 312.9 d0,   1.0 d0,  190.0  d0,  22.0  d0, 117.0  d0,
     8  47.0 d0,  22.0 d0,  150.0  d0, 119.0  d0, 246.0  d0,
     9 301.0 d0, 126.0 d0,  104.0  d0, 340.0  d0, 270.0  d0,
     a 358.0 d0, 148.0 d0/
!
         data h/
     1            9604.0088d0,  60.0316d0, -4073.3220d0, 8965.374d0,
     2   698.667d0, 13737.362d0,14375.997d0, -8845.31d0,-4711.96d0,
     3 -3713.33d0,  5470.66d0, 18509.35d0,  -4433.31d0, 8605.38d0,
     4  13377.37d0,  1058.66d0, 9244.02d0, -8206.68d0, 5192.01d0,
     5  14496.06d0,   420.02d0, 9284.69d0, 9964.00d0, - 299.96d0,
     6  4472.03d0,    379.35d0, 4812.68d0, -4851.36d0,19147.99d0,
     7 -12978.66d0, 17870.7d0, 9724.1d0, 13098.7d0, 5590.7d0,
     8 -13617.3d0,  -8485.3d0, 4193.4d0, -9483.9d0, 23282.3d0,
     9  10242.6d0,   9325.4d0,  14097.4d0, 22642.7d0,18149.4d0,
     a  -3353.3d0,  19268.0d0/
!
          t=(mjd-42412.d0)/365.25d0
          t=t + (0.0317d0*t+1.43d0)*1.d-6
!
          a = 0.0040d0*ksind(93.8d0  - 1.33d0*t)
     *       +0.0020d0*ksind(248.6d0 - 19.34d0*t)
     *       +0.0006d0*ksind(66.d0   + 0.2d0*t)
     *       +0.0006d0*ksind(249.d0  -19.3d0*t)
!
          b= 0.0267d0*ksind(68.64d0  - 19.341d0*t)
     *      +0.0043d0*ksind(342.d0   - 19.36d0*t)
     *      +0.0040d0*ksind( 93.8d0  -  1.33 d0*t)
     *      +0.0020d0*ksind(248.6d0  - 19.34d0*t)
     *      +0.0005d0*ksind(358.d0 -   19.4d0*t)
!               longitude
          tmp=124.8754d0+4812.67881d0*t +
     *        6.2887d0*ksind(338.915d0+ 4771.9886d0*t+a)
!
           do   i=1, 62
              tmp=tmp+ c(i)*ksind( d(i) + e(i)*t )
           enddo
          elon=mod(tmp,360.d0)
!               latitude
          tmp=5.1282d0*ksind(236.231d0 + 4832.0202d0*t+b)
!
           do   i=1, 46
              tmp=tmp+ f(i)*ksind( g(i) + h(i)*t )
           enddo
          elat=tmp
!
          sinpi= 0.9507d0
     *          + 0.0518d0*kcosd(338.92d0 + 4771.989d0*t)
     *          + 0.0095d0*kcosd(287.2 d0 - 4133.35 d0*t)
     *          + 0.0078d0*kcosd( 51.7 d0 + 8905.34 d0*t)
     *          + 0.0028d0*kcosd(317.8 d0 + 9543.98 d0*t)
     *          + 0.0009d0*kcosd( 31.0 d0 +13677.3  d0*t)
     *          + 0.0005d0*kcosd(305.0 d0 - 8545.4  d0*t)
     *          + 0.0004d0*kcosd(284.0 d0 - 3773.4  d0*t)
     *          + 0.0003d0*kcosd(342.0 d0 + 4412.0  d0*t)
!
          rmoon=ae/(sinpi*Torad)
       end
       real*8 function ksind(x)
!                sin x with x in degree
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
         ksind=sin(x*Torad)
       end
       real*8 function kcosd(x)
!                cos x with x in degree
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
         kcosd=cos(x*Torad)
       end
!       ************************************************************
!       *
!       * kctoq: ecliptic to equatorial coordinate transformation
!       * kqtoc: eqatorial to ecliptic  coordinate transformation
!       *
!       ************************************************************
! /usage/  call kctoq(mjd, cx, cy, cz, ex, ey, ez)
!
!   mjd: input.  real*8.  modified julian day of the place, time
!  cx,cy,cz:  input.        directional cosines in ecliptic coordinate.
!  ex,ey,ez:  output.      //        equatorial coordinate.
       subroutine kctoq(mjd, cx, cy, cz, ex, ey, ez)
!
          implicit real*8 (a-h, o-z)
          real*8  mjd
!           get mean obliquity of the ecliptic
          call kmobec(mjd, cose, sine)
          ex=cx
          ey=cy*cose - cz*sine
          ez=cy*sine + cz*cose
       end
!          eqatorial to ecliptic  coordinate transformation
!          call kqtoc(mjd, ex, ey, ez, cx, cy, cz)
       subroutine kqtoc(mjd, ex, ey, ez, cx, cy, cz)
          implicit real*8 (a-h, o-z)
          real*8  mjd
!           get mean obliquity of the ecliptic
         call kmobec(mjd, cose, sine)
         cx=ex
         cy=ey*cose + ez*sine
         cz=-ey*sine + ez*cose
       end
!      ********************************************************
!      *
!      *  kmobec: get mean obliquity of teh ecliptic plane
!      *
!      *         (inclination angle of ecliptic plane)
!      ********************************************************
! /usage/ call kmobec(mjd, cose, sine)
!  mjd: input. real*8  modified jd.
!  cose: output. cosine of obliquity
!  sine: output. sine of obliquity
!
       subroutine kmobec(mjd, cose, sine)
          implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
          real*8  mjd
!              get cos and sin of inclination angle of the
!              mean ecliptic plane.
!               get time from j2000.
           call kmtoj2(mjd, t)
!           eps=(((0.0000005036d0*t -0.000000164d0)*t - 0.01300417d0)*t
           eps=(((0.0000005036d0*t -0.00000164d0)*t - 0.01300417d0)*t
     *         +23.439291d0)*Torad
           cose=cos(eps)
           sine=sin(eps)
       end
!                   test ksun
!        real*8 mjd
!        mjd=48234.00068287d0
!        call ksun(mjd, slon, rsun)
!        write(*,*) ' slon=',slon, ' rsun=',rsun
!       end
!       ************************************************************
!       *
!       *  ksuneq: compute sun's position in the geocentric
!       *          equatorial coordinate
!       *
!       *        accuracy is about 1/100. degree
!       *  no need to convert this one into topocentric one in
!       *  this accuracy.
!       ************************************************************
! /usage/ call ksuneq(mjd, ex, ey, ez)
!
!     mjd: input. real*8.   modified jd.
!    ex,ey,ez: output. directional cosines in the geocentric
!                      equatorial coordinate
!
        subroutine ksuneq(mjd, ex, ey, ez)
          implicit real*8 (a-h, o-z)
          real*8 mjd
!             get ecliptic longitude in deg
          call ksun(mjd, slon, rsun)
!             convert to directional cos. in ecliptic
          call ketod(0.d0,slon, cx, cy, cz)
!             convert to equatorial coordinate
          call kctoq(mjd, cx, cy, cz, ex, ey, ez)
        end
!       ************************************************************
!       *
!       *  ksun: compute sun's position in the geocentric
!       *        ecliptic coordinate
!       *
!       *        accuracy is about 1/100. degree
!       *  no need to convert this one into topocentric one in
!       *  this accuracy.
!       ************************************************************
! /usage/ call ksun(mjd, slon, rsun)
!
!   mjd:  input. real*8.  modified julian day when the position is
!                wanted.
!  slon:  output.         ecliptic longitude in degree.
!  rsun:  output.         distance between the centers of earth and
!                  sun in m.
!  *** note ***
!         if the earth position in the heriocentric ecliptic coordinate
!       is wanted, it is given by slon+180.
!
!       the ecliptic latitude is always very small (< 1''), and is to
!       be set to zero in our approximation.
!       to convert the ecliptic coordinate into the equatorial
!       coordinate, use, ketod(0.d0, slon,  cx, cy, cz),
!       kctoq(mjd,cx, cy, xz, ex, ey, ez), kdtoe(ex, ey,ez, dec, ra)
!
!
       subroutine ksun(mjd, slon, rsun)
          implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
          real*8 mjd, ksind, kcosd
!
          t=(mjd-42412.d0)/365.25d0
          t=t + (0.0317d0*t+1.43d0)*1.d-6
          tmp = 279.0358d0 +360.00769d0*t
     1         +(1.9159d0-0.00005d0*t)*ksind(356.531d0+359.991d0*t)
     2         + 0.02d0             *ksind(353.06d0 + 719.981d0*t)
     3         - 0.0048d0           *ksind(248.64d0 - 19.341d0*t)
     4         + 0.0020d0           *ksind(285.0d0 + 329.64d0*t)
     5         + 0.0018d0           *ksind(334.2d0 -4452.67d0*t)
     6         + 0.0018d0           *ksind(293.7d0 - 0.20d0*t)
     7         + 0.0015d0           *ksind(242.4d0 + 450.37d0*t)
     8         + 0.0013d0           *ksind(211.1d0 + 225.18d0*t)
     9         + 0.0008d0           *ksind(208.0d0 + 659.29d0*t)
     a         + 0.0007d0           *ksind( 53.5d0 +  90.38d0*t)
     b         + 0.0007d0           *ksind( 12.1d0 -  30.35d0*t)
     c         + 0.0006d0           *ksind(239.1d0 + 337.18d0*t)
     d         + 0.0005d0           *ksind( 10.1d0 -   1.50d0*t)
     e         + 0.0005d0           *ksind( 99.1d0 -  22.81d0*t)
     f         + 0.0004d0           *ksind(264.8d0 + 315.56d0*t)
     g         + 0.0004d0           *ksind(233.8d0 + 299.30d0*t)
     h         + 0.0004d0           *ksind(198.1d0 + 720.02d0*t)
     i         + 0.0003d0           *ksind(349.6d0 + 1079.97d0*t)
     k         + 0.0003d0           *ksind(241.2d0 -44.43d0*t)
!
       slon=mod(tmp, 360.d0)
!
       q=(-0.007261d0+0.0000002d0*t)*kcosd(356.53d0 + 359.991d0*t)
     *     + 0.000030d0
     1    - 0.000091d0 * kcosd(353.1d0 + 719.98d0*t)
     2    + 0.000013d0 * kcosd(205.8d0 + 4452.67d0*t)
     3    + 0.000007d0 * kcosd( 62.d0  + 450.4d0*t)
     4    + 0.000007d0 * kcosd(105.d0  + 329.6d0*t)
!
       rsun=10.**q * aunit
      end
      subroutine kadbp(nftch,dx,dy,dz,dt,wt,u,v,w,tz,icon)
      implicit real*8  (a-h, o-z)
      dimension dt(nftch),dx(nftch),dy(nftch),dz(nftch),wt(nftch)
!----------------------------------------------------------------------
      sww=0.d0
      swx=0.d0
      swy=0.d0
      swz=0.d0
      swt=0.d0
      sxy=0.d0
      syz=0.d0
      szx=0.d0
      sxt=0.d0
      syt=0.d0
      sx2=0.d0
      sy2=0.d0
       do   i=1,nftch
        sww=sww+wt(i)
        swx=swx+dx(i)*wt(i)
        swy=swy+dy(i)*wt(i)
        swz=swz+dz(i)*wt(i)
        swt=swt+dt(i)*wt(i)
        sxy=sxy+dx(i)*dy(i)*wt(i)
        syz=syz+dy(i)*dz(i)*wt(i)
        szx=szx+dz(i)*dx(i)*wt(i)
        sxt=sxt+dx(i)*dt(i)*wt(i)
        syt=syt+dy(i)*dt(i)*wt(i)
        sx2=sx2+dx(i)*dx(i)*wt(i)
        sy2=sy2+dy(i)*dy(i)*wt(i)
       enddo
!
      a1=sww*sx2-swx*swx
      a2=sww*sxy-swx*swy
      a3=sww*szx-swx*swz
      a4=sww*sxt-swx*swt
      b1=a2
      b2=sww*sy2-swy*swy
      b3=sww*syz-swy*swz
      b4=sww*syt-swy*swt
      ab=a1*b2-a2*b1
      if(abs(ab).gt.1.d-6) then
        p=(a2*b3-a3*b2)/ab
        q=(a2*b4-a4*b2)/ab
        r=(a3*b1-a1*b3)/ab
        s=(a4*b1-a1*b4)/ab
        aa=p*p+r*r+1.0d0
        bb=0.6*(p*q+r*s)
        cc=0.09*(q*q+s*s)-1.0d0
        t1=-0.5d0*bb/aa
        t2=bb*bb-4.d0*aa*cc
        if(t2.lt.0.d0) then
!             % solution is imaginary %
          icon=2
        else
!             % solvable ] %
          t2=0.5d0*sqrt(t2)/aa
          w=t1+t2
          if(w.lt.0.d0) w=t1-t2
          u=p*w+0.3d0*q
          v=r*w+0.3d0*s
          tz=(u*swx+v*swy+w*swz+0.3d0*swt)/(0.3d0*sww)
          if(abs(u).le.1.0d0 .and. abs(v).le.1.0d0) then
            icon=0
          else
!                 % solution is not normalized %
            icon=1
          endif
        endif
!            direction cosines cannot be determined
      else
        icon=3
      endif
      return
      end
!  -------------------------------------------
!          normalize 3-d vector.
       subroutine knormv(a1, b1, c1, fn1)
            real*8 a1, b1, c1, fn1
            fn1=sqrt( a1**2+b1**2+c1**2)
            a1=a1/fn1
            b1=b1/fn1
            c1=c1/fn1
       end
!           get theta and fai of unnormalized 3-d vector
       subroutine kvtoa(vx, vy, vz, teta, fai)
         implicit real*8 (a-h, o-z)
!
           d=sqrt( vx**2 + vy**2 + vz**2)
           call kdtoa(vx/d, vy/d, vz/d, teta, fai)
       end
!             get theta and fai of direction cos.
       subroutine kdtoa(vx, vy, vz, teta, fai)
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
           if(vz .gt. 1.d0) then
              teta=0.
           else
              teta=acos(vz)
           endif
           if(abs(teta) .lt. 1.d-4) then
               fai=0.
           else
               fai=atan2(vy, vx)
           endif
!            to degree
           teta=Todeg*teta
           fai=Todeg*fai
       end
! ------------------------------------------------------------------
!           get difference of 2 vector angles.
!    (a1, b1, c1)  1 vector (direction cos).
!    (a2, b2, c2)  another //   these are input
!     difax, difay: projected angle difference. in deg. output
!    difa:  absolute difference. in deg. output
!           if direction cos's are invalid, -1. is put
       subroutine  kdifva(a1, a2, b1, b2, c1, c2, difax,
     * difay, difa)
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
!
            tetap(ww1, ww2)= asin( ww1/ sqrt(abs(1.d0-ww2**2)))* Todeg
!              component of angle difference
            tmp= a1*a2+b1*b2+c1*c2
            if(tmp .ge. 1.00d0) then
                difa=-1.d0
            else
                difax= tetap(a1,b1) - tetap(a2, b2)
                difay= tetap(b1,a1) - tetap(b2, a2)
                difa=acos( tmp  )*Todeg
            endif
       end
!      **************************************************************
!      *
!      *  komeg0:  init for komega and/or komeg1
!      *  komega:  get solid angle around a source given by komeg0
!      *           for a given point in celestial coord.
!      *  komeg1:  get theta between a source given by komeg0
!      *           and a given point in celestial coord.
!      *
!      **************************************************************
!      odec: input.   declination in deg
!       ora: input.   r.a in deg
!       dec: input.   declination id deg
!        ra: input.   r.a in deg
!     omega: output.  solid angle spaned by (dec, ra) around
!                     (odec, ora)
!      teta: output.  opening angle between (dec, ra) and
!                     (odec, ora).  in deg
       subroutine komeg0(odec, ora)
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
           parameter (coeff=pi*2.)
           save ex, ey, ez
           call ketod(odec, ora, ex, ey, ez)
           return
       entry komega(dec, ra, omega)
           call ketod(dec, ra, rx, ry, rz)
           cost=ex*rx + ey*ry + ez*rz
           omega= coeff * (1.d0- cost)
           return
       entry komeg1(dec, ra, teta)
           call ketod(dec, ra, rx, ry, rz)
           cost=ex*rx + ey*ry + ez*rz
           teta=acos(cost)* Todeg
       end
!      ********************** get opening angle between two
!      ********************** celestial coordinates
!   odec: input. declination in degree
!    ora: input. corresponding r.a in deg
!    dec: input. another declination in deg
!     ra: input. corresponding r.a in deg
!   teta: output. opening angle between two. in deg.
       subroutine koangl(odec, ora, dec, ra, teta)
         implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
           call ketod(odec, ora, ex, ey, ez)
           call ketod(dec, ra, rx, ry, rz)
           cost=ex*rx + ey*ry + ez*rz
           teta=acos(cost)* Todeg
       end
!         implicit real*8 (a-h, o-z)
!          real*8 mjd
!          call kcelei(36.005845953d0, 139.190535564d0, 9.d0,
!     *       875.8d0)
!
!          write(*,*) ' ug,vg,wg=', ug, vg, wg
!          time=21.d0 + 15.d0/60.d0
!          call kmjd(79, 9, 8, time, mjd)
!          call kmjdst(mjd, st)
!          write(*,*) ' mjd=',mjd
!          write(*,*) ' st=',st-tlons
!          ex=3942.937d3
!          ey=-4398.962d3
!          ez=4754.583d3
!          call knormv(ex,ey, ez, rs)
!          write(*,*) ' ex*rs=',ex*rs, ' ey*rs=',ey*rs,' ez*rs=',ez*rs
!          write(*,*) ' rs=',rs
!          call kdtoe(ex, ey, ez, de, ra)
!          write(*,*) ' ra=',ra, ' de=',de
!          call kgcttc(mjd, ex, ey, ez, rs, ex2, ey2, ez2)
!          call kdtoe(ex2, ey2, ez2, de, ra)
!          write(*,*) ' ra=',ra, ' de=',de
!         end
!         *************************************************************
!         *
!         *
!         * kgcttc: conversion from the geocentric to topocentric
!         *         coordinate (topocenter --> origin is at the
!         *                      observation place)
!         *
!         *  for our purpose, this may be needed only for the moon.
!         *************************************************************
!
!  /usage/ call kgcttc(mjd, ex, ey, ez, rs, tex, tey, tez)
!
!   mjd :       input. real*8 modified jd.
!   ex, ey, ez: input. direction cosines in the geocentric
!                      coordinate
!   rs.         input. distance between the center of the earth and
!                      the object (say, moon) (in m)
!   tex, tey, tez. output. directions cosines in the topocentric
!                      coordinate
!
!   ***note*** the observation place should have been specified by
!              call kcelei.
!              small correction due to the difference of the
!              geocenter and the center of the earth for the
!              local place is adjusted by using very approximate
!              values.
!
          subroutine kgcttc(mjd, ex, ey, ez, rs, tex, tey, tez)
           implicit real*8 (a-h, o-z)
#include  "Zkcele.h"
             real*8 mjd, ksind, kcosd
!             get mean greenwich siderial time in degree
             call kmjdst(mjd, st)
             st=st-tlons
             csg=kcosd(st)
             sng=ksind(st)
!              get geocentric equatorial coord. of the observation
!              place
             a=ug*csg - vg*sng
             b=ug*sng + vg*csg
             c=wg
!              convert the object position into topocentric coord.
             da=ex*rs - a
             db=ey*rs - b
             dc=ez*rs - c
!                normalize the vector, da,db,dc
             call knormv(da, db, dc, dummy)
             tex=da
             tey=db
             tez=dc
          end
!                    test kmoont
!         implicit real*8 (a-h, o-z)
!       real*8 mjd/48234.00068287d0/
!       common /$$$/ ex1, ey1, ez1
!       dmax=-1.
!       dmin=100.
!       do 100 i=1, 24
!         call kcelei(30.d0,  90.d0, 8.d0, 4300.d0)
!         call kmoont(mjd, ex, ey, ez)
!         call kdtoe(ex, ey, ez, dec, ra)
!         write(*,*) ' dec=', dec, ' ra=',ra
!         mjd=mjd+ 1./24.
!         call   kdifva(ex, ex1,ey,ey1, ez, ez1,  difax,
!    *     difay, difa)
!         dmax=max(dmax, difa)
!         dmin=min(dmin, difa)
! 100   continue
!        write(*,*) ' dmax=', dmax, ' dmin=', dmin
!       end
!       ************************************************************
!       *
!       * kmoont: compute the moon position in the topocentric
!       *         equatorial coordinate
!       *   (accuracy is better than 1/100 degree)
!       * kcelei must have been called.
!       ************************************************************
! /usage/  call kmoont(mjd, ex,  ey, ez)
!
!   mjd: input.  real*8.  modified julian day
!   ex, ey, ez.  apparent topocentric directional cosines in the
!                         topocentric equatorial coordinate.
!
!    ** geocentric and topocentric diff. is order of .5 deg
!
!
        subroutine kmoont(mjd, ex, ey, ez)
          implicit real*8 (a-h, o-z)
!         common /$$$/ ex1, ey1, ez1
          real*8 mjd
!             get ecliptic latitude  and longitude of the moon
          call kmoon(mjd,  elat, elon, rmoon)
!             convert to direction cosine in geocentric ecliptic
          call ketod( elat, elon, cx, cy, cz)
!             convert to geocentric equatorial coordinate
          call kctoq(mjd, cx, cy, cz, ex, ey, ez)
!$$$$$$$$$$$$$$$
!          call kdtoe(ex, ey, ez, dec, ra)
!          write(*,*) ' dec=', dec, ' ra=', ra
!          ex1=ex
!          ey1=ey
!          ez1=ez
!$$$$$$$$$$$$$$
!             convert to topocentric coord.
          call kgcttc(mjd, ex, ey, ez, rmoon, ex, ey, ez)
        end
!     ***********************************************************
!     *
!     *  kb50j2: convert dec, ra in b1950.0 into j2000
!     *   proper motion and e-term neglected
!     *  5/100 degree error may be included
!     *
!     *  kj2b50: inverse of the above
!     *
!     *
!     ***********************************************************
!       dec, ra: input. (dec, ra)(b1950.0)
!       dec2, ra2: output (dec, ra) (j2000)
      subroutine kb50j2(dec, ra, dec2, ra2)
        implicit real*8 (a-h, o-z)
             call ketod(dec, ra, ex, ey, ez)
!            write(*,*) ' ex, ,, ', ex, ey, ez
             ex2=.9999256782d0*ex-0.011182061d0*ey-0.0048579477d0*ez
             ey2=0.0111820609d0*ex+.9999374784d0*ey-0.0000271765d0*ez
             ez2=0.0048579479d0*ex-0.0000271474d0*ey+.9999881997d0*ez
!            write(*,*) ' ex2 ,, ', ex2, ey2, ez2
             call kdtoe(ex2, ey2, ez2, dec2, ra2)
      end
!       dec2, ra2: input  (dec, ra) (j2000)
!       dec, ra: output (dec, ra)(b1950.0)
      subroutine kj2b50(dec2, ra2, dec, ra)
        implicit real*8 (a-h, o-z)
             call ketod(dec2, ra2, ex2, ey2, ez2)
!            write(*,*) ' ex2, ,, ', ex2, ey2, ez2
           ex=.9999257080d0*ex2+0.0111789382d0*ey2+0.0048590039d0*ez2
           ey=-0.0111789382d0*ex2+.9999375133d0*ey2-0.0000271579d0*ez2
           ez=-0.0048590038d0*ex2-0.0000271626d0*ey2+.9999881946d0*ez2
!            write(*,*) ' ex ,, ', ex, ey, ez
           call kdtoe(ex, ey, ez, dec, ra)
      end
!     ***********************************************************
!     *
!     *  kj90j2: convert dec, ra in 1990 japanese ephemeris into j2000
!     *
!     *  kj2j90: inverse of the above
!     *
!     *
!     ***********************************************************
!       dec, ra: input. (dec, ra)(j1990.5)
!       dec2, ra2: output (dec, ra) (j2000)
      subroutine kj90j2(dec, ra, dec2, ra2)
        implicit real*8 (a-h, o-z)
             call ketod(dec, ra, ex, ey, ez)
!            write(*,*) ' ex, ,, ', ex, ey, ez
             ex2=.99999732d0*ex-0.00212430d0*ey-0.00092315d0*ez
             ey2=0.00212430d0*ex+.99999774d0*ey-0.00000098d0*ez
             ez2=0.00092315d0*ex-0.00000098d0*ey+.99999957d0*ez
!            write(*,*) ' ex2 ,, ', ex2, ey2, ez2
             call kdtoe(ex2, ey2, ez2, dec2, ra2)
      end
!       dec2, ra2: input  (dec, ra) (j2000)
!       dec, ra: output (dec, ra)(b1950.0)
      subroutine kj2j90(dec2, ra2, dec, ra)
        implicit real*8 (a-h, o-z)
             call ketod(dec2, ra2, ex2, ey2, ez2)
!            write(*,*) ' ex2, ,, ', ex2, ey2, ez2
           ex=.99999732d0*ex2+0.00212430d0*ey2+0.00092315d0*ez2
           ey=-0.00212430d0*ex2+.99999774d0*ey2-0.00000098d0*ez2
           ez=-0.00092315d0*ex2-0.00000098d0*ey2+.99999957d0*ez2
!            write(*,*) ' ex ,, ', ex, ey, ez
           call kdtoe(ex, ey, ez, dec, ra)
      end
!        ***********************************************************
!        *
!        * kjxjy: convert steller positions for date jx into jy
!        *
!        *
!        ***********************************************************
!
!  /usage/  call kjxjy(mjd1, mjd2, dec1, ra1, dec2, ra2)
!
!    mjd1: input. real*8  modified jd of j'x'
!    mjd2: input. real*8  modified jd of j'y'
!    dec1: input.         declination of a source in deg. at j'x'
!    ra1 : input.         raight ascensiion of the source in deg.
!    dec2: ouput.         declination in deg at j'y'
!    ra2 : ouput.         right ascension in deg at j'y'
!
         subroutine kjxjy(mjd1, mjd2, dec1, ra1, dec2, ra2)
            implicit real*8 (a-h, o-z)
            real*8 mjd1, mjd2
            dimension pij1(3,3), pij2(3,3)
!
!              compute precession matrix
            call kpmtrx(mjd1, pij1)
!
!           do 100 j=1, 3
!              write(*,*) (pij1(i, j), i=1, 3)
! 100       continue
            call kpmtrx(mjd2, pij2)
!           do 200 i=1, 3
!              write(*,*) (pij2(i, j), j=1, 3)
! 200       continue
!                 get vectors
            call ketod(dec1, ra1, ex1, ey1, ez1)
!                 convert     x-->j2000
            call kxtoj2(pij1, ex1, ey1, ez1, ex2, ey2, ez2)
!                 convert  j2000--->y
            call kj2tox(pij2, ex2, ey2, ez2,  ex, ey, ez)
!               convert to dec, ra
            call kdtoe(ex, ey, ez, dec2, ra2)
         end

