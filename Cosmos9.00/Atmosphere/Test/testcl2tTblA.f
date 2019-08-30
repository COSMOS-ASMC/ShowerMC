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
      real*8 clen2thickEx, clen2thickTA, ct2lTA, thicke
      integer loc



      cosz = 0.0
      Hbase = 1400.d0
      Htop = 100.d3

      write(*, *) ' enter cosz at base'
      read(*, *) cosz

      call creadAtmosD
      r1 = Htop + Eradius
      r2 = Hbase + Eradius
      write(0,*) ' Re=',Eradius
      write(0,*) ' r1=',r1, ' r2=',r2
      leng = clenbetween2h(r2, r1, cosz)
      cosx = cnewcos(r2, cosz, leng)
      write(0, *)' cos at top=', cosx,' length to top=', leng

      write(0,*) 'LenStep=',LenStep
      write(0,*) 'maxl2t=', maxl2t, ' NumStep=',NumStep

      call cl2tTbl(Htop, Hbase, cosx, cosz, LenStep,
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
         thick = clen2thickTA(z2,  leng)
         thicke = clen2thickT(z2, cosx, leng)
         write(*, *) ' thick=', thick, thicke, ' leng=', leng
         write(*,*) ' which should be =', ct2lTA(z2,  thick), 
     *      ct2lT(z2, cosx, thicke)
         write(*, *) ' exact thick=', clen2thickEx(z2, cosx, leng, 6)
         write(*,*) ' enter  height'
         read(*, *) z2
         r2 = z2 + Eradius
      enddo
      end
#include "../cl2tTbl.f"
#include "../cl2tTblA.f"
