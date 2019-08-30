!        Zkcele.h
!     unit here is SI. (m,s,kg)
       real*8 sidcor, sidcr2, sidcr3, pi, Torad, Todeg,
     *            oneday, tofai, toh, vlight, taua, ae,
     *            ctaua, aunit
      parameter (
     * pi=3.141592653589793238d0, Torad=pi/180.d0, Todeg=1.d0/Torad,
     * oneday=24.d0, tofai=360.d0/oneday, toh=oneday/360.d0,
     * sidcor=1.002737909350795d0, sidcr2=5.9006d-11, sidcr3=5.9d-15,
     * vlight=2.99792458d8,
     * taua=499.004782d0,
     * ae=6378140.d0,
     * ctaua=vlight*taua, aunit=ctaua)

       real*8 tlons, tlats, dtgmts, sinlat, coslat, sinsx, cossx,
     *        heighs, ug, vg, wg
!
       common /Zkcele/  tlons, tlats, dtgmts, sinlat, coslat,
     * sinsx, cossx, heighs, ug, vg, wg
!
!  ae:    equatorial radius for earth
!  taua:  light-time for unit distance
!  ctaua=aunit. unit distance
!  tlons: longitude of observation place in deg
!  tlats: latitude of  // (- for south)
!  dtgmts: difference of the local time from greenwich time in hour
!  sinlat: sin of tlats
!  coslat: cos of tlats
!  sinsx: sin of south to x-axis of the detector
!  cossx: cos  //
!  heighs: height of the observation place in m
!  ug, vg, wg: coordinate of the observation place in the
!              geocentric rectangular coordinate. (m)
