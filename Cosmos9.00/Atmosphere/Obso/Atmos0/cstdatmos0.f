!       Units are in SI.

!   below cthick2h1 and cthick2h2 are so made that they coincide at
!   t=230

!
!      -------------------------------------
       real*8 function cvh2den(z)
!      --------------------------vertical height to density
!        at z > 50 km, very bad. but ok for our purpose
       implicit none
!----       include 'Zstdatmos.h'
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
!
!     -------------------------------------------------------------
!      real*8 z, cvh2temp
!      do z = 0., 1250.d3, 1d3
!         write(*, *) sngl(z/1.d3), sngl(cvh2temp(z))
!      enddo
!      end

      real*8 function cvh2temp(z)
!          vettical height to temperatur (Kelvin)

      real*8 z   ! input.  vertical height in m
!        output is temperature of the atmospher in Kelvin

      real*8 ans, error, hx
      integer  j, m,  loc, nh

      parameter(nh=17)

      real*8 h(nh), t(nh)  ! h (km) vs t 
      
      data h/91, 100, 110, 120, 130, 140, 160, 180, 200, 250,
     *      300, 350, 400, 450, 500, 550, 600/	
      data t/186.87, 	195.08,	240.0,	360.0,	469.27,	559.63,
     *	696.29,	790.07,	854.56,	941.33,	976.01,	990.06,	995.83,
     *	998.22, 999.24,	999.67,	999.85/
      data m/2/

      if(z .lt. 11.1d3) then
         cvh2temp = (216.65d0 - 288.15d0)/11.1d3 * z + 288.15d0
      elseif(z .lt. 20.0d3) then
         cvh2temp = 216.65d0
      elseif(z .lt. 32.2d3) then
         cvh2temp = (228.756d0 - 216.65d0)/(32.2d3-20.d3) * (z-20.d3)
     *            + 216.65d0
      elseif(z .lt. 47.4d3) then
         cvh2temp = (270.65d0-228.756d0)/(47.4d3-32.2d3) * (z-32.2d3)
     *            + 228.756d0
      elseif(z .lt. 51.d3) then
         cvh2temp = 270.650
      elseif(z .lt. 72.d3) then
         cvh2temp = (214.263d0 - 270.65d0)/(72.d3- 51.0d3) * (z-51.d3)
     *            + 270.65d0
      elseif(z .lt. 86.d3) then
         cvh2temp = (186.87d0-214.263d0)/(86.d3-72.d3)* (z -72.d3)
     *            + 214.263d0
      elseif( z .lt. 91.d3) then
         cvh2temp = 186.87d0
      elseif( z .lt. 600d3) then
         hx = z/1.d3  ! in km
         call kdwhereis(hx, nh, h, 1,  loc)
         j = min(max(loc - (m-1)/2,1), nh+1-m)  ! max of m points from j
         call kpolintp(h(j), 1, t(j), 1, m, hx, ans, error)
         cvh2temp = ans
      else
         cvh2temp = 1000.
      endif
      end

!---------------------------------------------
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
!      -------------------------------------
       real*8 function cthick2den(t)
!      -------------------------------------
       implicit none
!----       include 'Zstdatmos.h'
#include  "Zstdatmos.h"
       real*8 t

       if(t .lt. t0) then
          cthick2den = t/hn
       else
          cthick2den =10.d0*hlhmi*
     *         (t/10.d0)**(1.d0 - hm)
       endif
       end
!      -------------------------------------
       real*8 function cvh2denp(z)
!      -------------------------------------

!          d rho/dz
       implicit none
!----       include 'Zstdatmos.h'
#include  "Zstdatmos.h"
       real*8 z
       real*8 zsave, ans, cvh2den
       data zsave /-1.d30/
       save   ans, zsave

       if(z .ne. zsave) then
          if(z .ge. hc) then
!             cvh2denp = - 10.d0*exp( (hb-z)/hn)/hn/hn
!             cvh2denp = cnst1*exp( (hb-z)/hn) 
             ans = - cvh2den(z)/hn
          else
!              cvh2denp =- 10.d0*hlhmi*(hmi-1.d0)/hl*
!     *        ( (ha-z)/hl)**(hmi-2.d0)
!              cvh2denp =cnst2 *
!     *          ( (ha-z)/hl)**cnst3
             ans = -(hmi-1.0d0)/(ha-z)*cvh2den(z)
          endif
          zsave= z
       endif
       cvh2denp = ans
       end
!      ----------------------------------
       real*8 function cvh2scaleh(z)
!      ----------------------------------
!          get scale height defined by 
!         -   rho/ (d rho/d z).  This has discontinuity at
!       z = 11.1 km
       implicit none
       real*8 z
       real*8 cvh2den, cvh2denp
!
       cvh2scaleh = - cvh2den(z)/ cvh2denp(z)
       end
!      -------------------------------------
       real*8 function cvh2den2p(z)
!      -------------------------------------
!          d(d rho/dz)/dz
       implicit none
!----       include 'Zstdatmos.h'
#include  "Zstdatmos.h"
       real*8 z
!       logical firsttime/.true./
!       save firsttime
!       real*8 cnst1, cnst2, cnst3
!       save cnst1, cnst2, cnst3
       real*8 zsave, ans, cvh2denp, cvh2den
       data zsave/-1.d30/
       save zsave
!
!       if(firsttime) then
!          cnst1 = 10.d0/hn/hn/hn
!          cnst2 = 10.d0*hlhmi*(hmi-1.d0)*(hmi-2.d0)
!     *          /hl/hl 
!          cnst3 = hmi -3.d0
!          firsttime = .false.
!       endif
       if(z .ne. zsave) then
          if(z .ge. hc) then
!            cvh2den2p =  10.d0*exp( (hb-z)/hn)/hn/hn/hn
!            cvh2den2p =  cnst1*exp( (hb-z)/hn)
             ans =  cvh2den(z)/hn/hn
          else
!              cvh2den2p = 10.d0*hlhmi*(hmi-1.d0)*(hmi-2.d0)
!     *          /hl/hl *
!     *        ( (ha-z)/hl)**(hmi-3.d0)
!             cvh2den2p = cnst2 * ( (ha-z)/hl )** cnst3
             ans = -(hmi-2.d0)/(ha -z) * cvh2denp(z)
          endif
          zsave =z
       endif
       cvh2den2p = ans
       end
!      -------------------------------------
       real*8 function cvh2den3p(z)
!      -------------------------------------
!          rho'''(z)
       implicit none
!----       include 'Zstdatmos.h'
#include  "Zstdatmos.h"
       real*8 z
!       logical firsttime/.true./
!       save firsttime
!       real*8 cnst1, cnst2, cnst3
!       save cnst1, cnst2, cnst3
       real*8 zsave, ans, cvh2den2p, cvh2den
       data zsave/-1.d30/
       save zsave, ans
!
!       if(firsttime) then
!          cnst1 = - 10.d0/hn/hn/hn/hn
!          cnst2 =- 10.d0*hlhmi*(hmi-1.d0)*(hmi-2.d0)*(hmi-3.d0)
!     *          /hl/hl/hl
!          cnst3 = hmi -4.d0
!          firsttime = .false.
!       endif
       if(z .ne. zsave) then
          if(z .ge. hc) then
!             cvh2den3p = - 10.d0*exp( (hb-z)/hn)/hn/hn/hn/hn
!             cvh2den3p =  cnst1*exp( (hb-z)/hn)
             ans = - cvh2den(z)/hn/hn/hn
          else
!               cvh2den3p =- 10.d0*hlhmi*(hmi-1.d0)*(hmi-2.d0)*(hmi-3)
!     *          /hl/hl/hl *
!     *        ( (ha-z)/h1)**(hmi-4.d0)
!             cvh2den3p = cnst2 * ( (ha-z)/hl )** cnst3
             ans = -(hmi-3.0d0)/(ha-z) * cvh2den2p(z)
          endif
          zsave =z
       endif
       cvh2den3p = ans
       end
!      ---------------------------------------
       real*8 function cvh2thick(z)
!      ---------------------------------------
       implicit none
!----       include 'Zstdatmos.h'
#include  "Zstdatmos.h"
       real*8 z

       if(z .gt. hc) then
          cvh2thick = 10.d0* exp( (hb-z)/hn) 
       else
          cvh2thick = 10.d0* ((ha-z)/hl)**hmi 
!
!            the next factor is -1.5d-11 and is neglected.
!     *              -  10.d0* ((ha-hc)/hl)**hmi +
!     *                10.d0*exp( (hb-hc)/hn)
!
       endif
       end
       block data cblkStdAtmos
#include  "Zstdatmos.h"
!            main frame cosmos data. 
!            d rho/dz is discontinuous at hc, but others are
!            continuous.  Thickness values are almost the same
!            as the standard obained by cstdatmos0-multi-seg.f
!
      data ha/4512224.7572830657d-2/, hb/4541933.9782793734d-2/,
     *     hl/1244541.6061892177d-2/, hm/.1854893358365053d0/,
     *     hn/6.3300000000746224d3/, t0/230.000458235099d1/,
     *     hc /10.996296495761545d3/
      data hmi/5.3911455097420/, hlhmi/4.3318322850207d-04/

!      coef. to have continuous density values for h=h(t) and its derivative
!      about t=t0.  However, the absolute thickness is little bit larger
!      than the standard, so that we don't use this one.
!
!      data ha/25289.093750000d0/, hb/48490.145953733/,
!     *     hl/2712.61d0/, hm/0.32397778012688d0/,
!     *     hn/6800 .0/, hmi/3.0866314338235/, hc/11.1d3/,
!     *     hlhmi/1.1378815278872D-03/, t0/2443.3759999999/
      end
      
  
