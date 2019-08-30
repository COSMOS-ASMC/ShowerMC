      implicit none
#include  "Zglobalc.h"
#include  "Zcoord.h"
#include  "Zmagfield.h"
      real*8 yearin, lat, long, height
      type(coord):: llh
      type(magfield):: h
      integer icon
                                  
      character*100 filepath
      filepath='../../../Data/Geomag/igrf'


      yearin= 2005.0
      write(*,*)  'y=',yearin
      write(*,*)  ' Change above values, if want or enter /'
      read(*,*)  yearin

      call  crdGeomag(filepath, yearin)

      lat = 30.d0
      long = 140.d0
      height = 0.d0

      do while (.true.)
         write(*,*) ' lat=',lat, ' long=',long, ' h=', height
         write(*,*)  ' Change above values, if want  or "99 /" to stop'
         read(*,*)  lat, long, height
         if( lat > 90. ) stop

         call csetCoord('llh', lat, long, height, llh)
         call cgeomag(yearin, llh, h,  icon)

         write(*,*)  ' icon=',icon, ' h%x,y,z=',h%x, h%y, h%z
         write(*,*)  ' |B|=',sqrt(h%x**2 + h%y**2 + h%z**2)
      enddo

      end
