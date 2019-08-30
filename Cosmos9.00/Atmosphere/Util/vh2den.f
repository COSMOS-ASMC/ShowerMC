      implicit none
      real*8  vh, den, vhmin, vhmax, dh
      real*8  cthick2h
      integer i
      
      vhmin =  20000.
      vhmax =  1.
      dh = 10.**0.025
      vh = vhmin
      do i = 1, 1000000
         if(vh .lt. vhmax) goto 100
         den = cthick2h(vh)
         write(*, *) vh, den
         vh = vh/dh
      enddo
 100  continue
      end
c       Units are in SI.

c   below cthick2h1 and cthick2h2 are so made that they coincide at
c   t=230

c
c      -------------------------------------
       real*8 function cvh2den(z)
c      --------------------------vertical height to density
c        at z > 50 km, very bad. but ok for our purpose
       implicit none
c----       include 'Zstdatmos.h'
#include  "Zstdatmos.h"
       real*8 z
       real*8 zsave, ans
       data zsave/-1.d30/
       save zsave, ans


       if(z .ne. zsave) then
          if(z .ge. hc) then
             ans = 10.d0*exp( (hb-z)/hn)/hn 
          else
             ans = 10.d0*hlhmi*( (ha-z)/hl)**(hmi-1.d0)
          endif
          zsave = z
       endif
       cvh2den = ans
       end
c
c     -------------------------------------------------------------
      real*8 function cthick2h(t)
      implicit none
      real*8 t
#include  "Zstdatmos.h"
      external cblkStdAtmos


      if(t .lt. t0) then
         cthick2h =  hb - log(t/10.d0) * hn    !    for t<t0 = 2300
      else
         cthick2h = ha - hl * (t/10.d0)**hm      !    for t > t0 
      endif
      end
c      -------------------------------------
       real*8 function cthick2den(t)
c      -------------------------------------
       implicit none
c----       include 'Zstdatmos.h'
#include  "Zstdatmos.h"
       real*8 t

       if(t .lt. t0) then
          cthick2den = t/hn
       else
          cthick2den =10.d0*hlhmi*
     *         (t/10.d0)**(1.d0 - hm)
       endif
       end
c      -------------------------------------
       real*8 function cvh2denp(z)
c      -------------------------------------

c          d rho/dz
       implicit none
c----       include 'Zstdatmos.h'
#include  "Zstdatmos.h"
       real*8 z
       real*8 zsave, ans, cvh2den
       data zsave /-1.d30/
       save   ans, zsave

       if(z .ne. zsave) then
          if(z .ge. hc) then
c             cvh2denp = - 10.d0*exp( (hb-z)/hn)/hn/hn
c             cvh2denp = cnst1*exp( (hb-z)/hn) 
             ans = - cvh2den(z)/hn
          else
c              cvh2denp =- 10.d0*hlhmi*(hmi-1.d0)/hl*
c     *        ( (ha-z)/hl)**(hmi-2.d0)
c              cvh2denp =cnst2 *
c     *          ( (ha-z)/hl)**cnst3
             ans = -(hmi-1.0d0)/(ha-z)*cvh2den(z)
          endif
          zsave= z
       endif
       cvh2denp = ans
       end
c      ----------------------------------
       real*8 function cvh2scaleh(z)
c      ----------------------------------
c          get scale height defined by 
c         -   rho/ (d rho/d z).  This has discontinuity at
c       z = 11.1 km
       implicit none
       real*8 z
       real*8 cvh2den, cvh2denp
c
       cvh2scaleh = - cvh2den(z)/ cvh2denp(z)
       end
c      -------------------------------------
       real*8 function cvh2den2p(z)
c      -------------------------------------
c          d(d rho/dz)/dz
       implicit none
c----       include 'Zstdatmos.h'
#include  "Zstdatmos.h"
       real*8 z
c       logical firsttime/.true./
c       save firsttime
c       real*8 cnst1, cnst2, cnst3
c       save cnst1, cnst2, cnst3
       real*8 zsave, ans, cvh2denp, cvh2den
       data zsave/-1.d30/
       save zsave
c
c       if(firsttime) then
c          cnst1 = 10.d0/hn/hn/hn
c          cnst2 = 10.d0*hlhmi*(hmi-1.d0)*(hmi-2.d0)
c     *          /hl/hl 
c          cnst3 = hmi -3.d0
c          firsttime = .false.
c       endif
       if(z .ne. zsave) then
          if(z .ge. hc) then
c            cvh2den2p =  10.d0*exp( (hb-z)/hn)/hn/hn/hn
c            cvh2den2p =  cnst1*exp( (hb-z)/hn)
             ans =  cvh2den(z)/hn/hn
          else
c              cvh2den2p = 10.d0*hlhmi*(hmi-1.d0)*(hmi-2.d0)
c     *          /hl/hl *
c     *        ( (ha-z)/hl)**(hmi-3.d0)
c             cvh2den2p = cnst2 * ( (ha-z)/hl )** cnst3
             ans = -(hmi-2.d0)/(ha -z) * cvh2denp(z)
          endif
          zsave =z
       endif
       cvh2den2p = ans
       end
c      -------------------------------------
       real*8 function cvh2den3p(z)
c      -------------------------------------
c          rho'''(z)
       implicit none
c----       include 'Zstdatmos.h'
#include  "Zstdatmos.h"
       real*8 z
c       logical firsttime/.true./
c       save firsttime
c       real*8 cnst1, cnst2, cnst3
c       save cnst1, cnst2, cnst3
       real*8 zsave, ans, cvh2den2p, cvh2den
       data zsave/-1.d30/
       save zsave, ans
c
c       if(firsttime) then
c          cnst1 = - 10.d0/hn/hn/hn/hn
c          cnst2 =- 10.d0*hlhmi*(hmi-1.d0)*(hmi-2.d0)*(hmi-3.d0)
c     *          /hl/hl/hl
c          cnst3 = hmi -4.d0
c          firsttime = .false.
c       endif
       if(z .ne. zsave) then
          if(z .ge. hc) then
c             cvh2den3p = - 10.d0*exp( (hb-z)/hn)/hn/hn/hn/hn
c             cvh2den3p =  cnst1*exp( (hb-z)/hn)
             ans = - cvh2den(z)/hn/hn/hn
          else
c               cvh2den3p =- 10.d0*hlhmi*(hmi-1.d0)*(hmi-2.d0)*(hmi-3)
c     *          /hl/hl/hl *
c     *        ( (ha-z)/h1)**(hmi-4.d0)
c             cvh2den3p = cnst2 * ( (ha-z)/hl )** cnst3
             ans = -(hmi-3.0d0)/(ha-z) * cvh2den2p(z)
          endif
          zsave =z
       endif
       cvh2den3p = ans
       end
c      ---------------------------------------
       real*8 function cvh2thick(z)
c      ---------------------------------------
       implicit none
c----       include 'Zstdatmos.h'
#include  "Zstdatmos.h"
       real*8 z

       if(z .gt. hc) then
          cvh2thick = 10.d0* exp( (hb-z)/hn) 
       else
          cvh2thick = 10.d0* ((ha-z)/hl)**hmi 
c
c            the next factor is -1.5d-11 and is neglected.
c     *              -  10.d0* ((ha-hc)/hl)**hmi +
c     *                10.d0*exp( (hb-hc)/hn)
c
       endif
       end
       block data cblkStdAtmos
#include  "Zstdatmos.h"
c            main frame cosmos data. 
c            d rho/dz is discontinuous at hc, but others are
c            continuous.  Thickness values are almost the same
c            as the standard obained by cstdatmos0-multi-seg.f
c
      data ha/4512224.7572830657d-2/, hb/4541933.9782793734d-2/,
     *     hl/1244541.6061892177d-2/, hm/.1854893358365053d0/,
     *     hn/6.3300000000746224d3/, t0/230.000458235099d1/,
     *     hc /10.996296495761545d3/
      data hmi/5.3911455097420/, hlhmi/4.3318322850207d-04/

c      coef. to have continuous density values for h=h(t) and its derivative
c      about t=t0.  However, the absolute thickness is little bit larger
c      than the standard, so that we don't use this one.
c
c      data ha/25289.093750000d0/, hb/48490.145953733/,
c     *     hl/2712.61d0/, hm/0.32397778012688d0/,
c     *     hn/6800 .0/, hmi/3.0866314338235/, hc/11.1d3/,
c     *     hlhmi/1.1378815278872D-03/, t0/2443.3759999999/
      end
