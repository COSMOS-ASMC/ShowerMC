      implicit none
#include  "Zglobalc.h"
#include  "Zcoord.h"
#include  "Zmagfield.h"
! #include  "Zearth.h"
      real*8 yearin, lat, long, height
      type(coord):: llh
      type(magfield):: h
      integer icon
                                  
      character*100 filepath
      filepath='../../../Data/Geomag/igrf'


      yearin= 2010.0
      write(0,*)  'y=',yearin
      write(0,*)  ' Change above values, if want or enter /'
      read(*,*)  yearin

      call  crdGeomag(filepath, yearin)

      lat = 30.d0
      long = 140.d0
      height = 0.d0

      do long = -180., 180., 2.49999999999d0
         do lat = -89.9, 89.9, 2.49999999999d0
            call csetCoord('llh', lat, long, height, llh)
            call cgeomag(yearin, llh, h,  icon)
            if(icon /= 0 ) then
               write(0,*) ' icon = ', icon, 'from cgeoma'
               stop
            endif
!            write(*,*)  ' icon=',icon, ' h.x,y,z=',h.x, h.y, h.z
!            write(*,*) ' |B|=',sqrt(h.x**2 + h.y**2 + h.z**2)
            write(*,'(1p, 3g13.4)') 
     *         long, lat, sqrt(h%x**2 + h%y**2 + h%z**2)
         enddo
         write(*,*) " "
      enddo
      end
