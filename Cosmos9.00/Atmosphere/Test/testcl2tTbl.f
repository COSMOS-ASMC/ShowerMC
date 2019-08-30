c     cl2tTbl: make table for length to thickness conversion.
c
c     Thickness of air corresponding to a given length along
c     a given direction from a given height.
c
c        test program
c
      implicit none
#include "Zearth.h"
#include "Zatmos.h"

      real*8  cosz,  z2, leng, cosx,  clenbetween2h
      real*8 cnewcos, clen2thickT, ct2lT, thick, r1, r2
      real*8 clen2thickEx
      integer loc

      cosz = 0.0
      Hbase = -1.d0
      write(*, *) ' enter cosz at base'
      read(*, *) cosz
      Htop = 30.d3

      r1 = Htop + Eradius
      r2 = Hbase + Eradius

      leng = clenbetween2h(r2, r1, cosz)
      cosx = cnewcos(r2, cosz, leng)
      write(*, *)' cos at top=', cosx,' length to top=', leng

      call cl2tTbl(Htop, Hbase, cosx, LenStep,
     *     LenTbl, HeightTbl, CosTbl,  ThickTbl, maxl2t, NumStep)
      write(*,*) 'LenStep=',LenStep, ' NumStep=', NumStep


      r1 = Hbase + Eradius
      z2 = 10.d3
      r2 = z2 + Eradius
      do while(z2 .gt. -0.3d0)
         leng = clenbetween2h(r1, r2, cosz)
         cosx = cnewcos(r1, cosz, leng)
         write(*, *) ' enter length from h=', z2, ' cos=', cosx
         read(*,*) leng
         thick = clen2thickT(z2, cosx, leng)
         write(*, *) ' thick=', thick,  ' leng=', leng
         write(*,*) ' which should be =', ct2lT(z2, cosx, thick)
         write(*, *) ' exact thick=', clen2thickEx(z2, cosx, leng, 6)
         write(*,*) ' enter  height'
         read(*, *) z2
         r2 = z2 + Eradius
      enddo
      end
#include "../cl2tTbl.f"
