!          lat, long, horizontal, vertical, def. angle
!
         program Geomag
         implicit none

#include  "Zcoord.h"
#include  "Zmagfield.h"
         real*8 year, lat, long, h
         type (coord):: llh
         type (magfield):: b, hva
         integer icon, igrf
         real*8  tand
!
         write(*,*)
     *  'Enter year(1996.5) lat(-90~90) long(-180~180)  height(0)'
         read(*, *) year, lat, long, h
         write(*,*)
     *   "Enter 1 or 2 for 1) IGRF data 2) Dod WMM data"
         read(*,*) igrf
         call csetCoord('llh', lat, long, h, llh)
         if(igrf .eq. 1) then
            call crdGeomag('../Data/Geomag/igrf', year)
         else
            call crdGeomag('../Data/Geomag/wmm', year)
         endif

         call cprGeomag

         call cgeomag(year, llh, b, icon)
         write(*,*) ' icon=',icon
         write(*,*) ' ned=', b%n, b%e, b%d
         call cned2hva(b, hva)
         tand = atan( hva%v/hva%h) * 180./3.141592
         call ctransMagTo('xyz', llh, b, b)
         write(*,*) ' x=',b%x, ' y=', b%y, ' z=', b%z
         write(*, *)
#ifdef  UNIONMAP
     *         sngl(hva%h), sngl(hva%v), sngl(hva%a)
#else
     *         sngl(hva%x), sngl(hva%y), sngl(hva%z)
#endif
         write(*,*) ' fukaku =', tand, ' deg' 
         call cprintMagF(b)
         call cprintMagF(hva)
         end
