c
c     test cthik2len.f
c
      implicit none
#include "Ztrack.h"
      record /track/ aTrack
      real*8  len, z, cosz, thick, cutt, thick2
      real*8  thick3
      real*8  clen2thickEx,   clen2thickAp
      integer ios, jcut, i, nc
c          read segmented atmosphere data
      call creadAtmosD
c          manipulate data
      call catmosCnst1
                  
      aTrack.pos.height =20.d3
      aTrack.vec.coszenith = 1.0d0
      aTrack.p.code = 1
      thick = 0.1d0
      nc = 0
      do i= 1, 10000
         call cthick2len(aTrack, thick, len, cutt, jcut)
         thick2= clen2thickEx(aTrack.pos.height,
     *   aTrack.vec.coszenith,  len, 10)
         thick3= clen2thickAp(aTrack.pos.height,
     *   aTrack.vec.coszenith,  len)
         
         if( (cutt-thick2) / cutt .gt. 1.d-4 .or.
     *       (cutt-thick3) / cutt .gt. 1.d-4 ) then
            write(0,*) thick, len, cutt, jcut, thick2, thick3
         endif
         if( jcut .eq. 1 .and. nc .lt. 10) then
            write(0,*) thick, len, cutt, jcut, thick2, thick3
            nc = nc +1
         elseif( jcut .eq. 1) then
            stop
         endif
         write(*,*)  thick, len
         thick = thick + 1.0d0
      enddo
      end
