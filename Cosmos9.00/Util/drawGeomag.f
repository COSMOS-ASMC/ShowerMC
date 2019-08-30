!             makes a table of
!          lat, long, horizontal, vertical, def. angle
!
!          for lat = -90 to 90 step 5 deg
!          long =-180 to 180 step 10 deg.
!  the resultant data can be drawn by, say, gnuplot. 
         program drawGeomag
         implicit none

#include  "Zcoord.h"
#include  "Zmagfield.h"
         real*8 year, lat, long, h
         type (coord):: llh
         type (magfield):: b, hva
         integer icon, ilat, ilong, igrf
         character*64 file
!
         year=2000.
         h =0.
         call cerrorMsg(
     *   'Enter year of geomagnetism (=2000) and height(=0 m)', 
     *    1)
         read(*,*) year, h

         call cerrorMsg(
     *   'Enter 1, 2 or 3 for 1) igrf, 2) wmm or other data', 1)
         read(*,*) igrf

         if(igrf .eq. 1) then
            file = '../Data/Geomag/igrf'
         elseif(igrf .eq. 2) then
            file = '../Data/Geomag/wmm'
         else
            call cerrorMsg('Enter file path for geomag data',1)
            read(*,*) file
         endif

         call crdGeomag(file, year)
         write(*,'("#     lat   long   H   V (Tesla) Ang(deg) ")')
         do ilat = -90, 90,  5
            lat = ilat
            do ilong = -180, 180, 10
               long = ilong
               call csetCoord('llh', lat, long, h, llh)
               call cgeomag(year, llh, b, icon)
               call cned2hva(b, hva)
!              call ctransMagTo('xyz', llh, b, b)
!              write(*,*) ' x=',b.x, ' y=', b.y, ' z=', b.z
               write(*, *) float(ilat), float(ilong), 
#ifdef  UNIONMAP
     *         sngl(hva%h), sngl(hva%v), sngl(hva%a)
#else
     *         sngl(hva%x), sngl(hva%y), sngl(hva%z)
#endif
!c             call cprintMagF(b)
!c             call cprintMagF(hva)
            enddo
            write(*,*)
         enddo   
         end
