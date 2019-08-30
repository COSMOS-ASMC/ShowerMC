      subroutine toasciiH(cnum, num, ir, zf)
      implicit none
#include "Ztrack.h"
      integer cnum, num, ir(2)
      type(track)::zf

      write(*, '("#")')
      write(*,
     *    '("#                    ********************** ")')
      write(*, '("#")')
      write(*,
     * '("#    cum event #  current  #       seeds")')
      write(*,'(i10,i10,8x,2i11)') cnum, num, ir
      write(*,
     *'("# 1ry code subc chg       E(GeV)     cos_zen  1st z(g/cm2)")')
      write(*,'(i8,2i5,4x,g12.3,3x,f8.5,g12.3)')
     *    zf.p.code, zf.p.subcode, zf.p.charge,
     *    zf.p.fm.p(4), zf.vec.coszenith,
     *    zf.pos.depth/10.
      end
      subroutine toasciiN(nlow, pp)
      implicit none
#include "Zprivate.h"
#include "Ztrack.h"
      integer nlow
      type(parent)::pp

      if(nlow .gt. 0) then
         write(*,'("#  nodal point info:")')
         write(*,'(a)')
     *  '#  # of child   code         E          H(m)'//
     *  '     Depth(g/cm2)      cos_zen'
         write(*,'(i10,i10,3x,g12.3,2x,g12.4,g12.3,2x,f9.5)')
     *   nlow, pp.code,  pp.erg, pp.height,
     *   pp.depth/10., pp.coszenith
         write(*,'("  code subc chg   Energy(GeV)")')
      endif
      end

      subroutine toasciiC(cc)
!      print child at nodal point
      implicit none
#include "Zprivate.h"
      type(child)::cc
      
      write(*,'(3i5, g12.3)')
     *      cc.code, cc.subcode, cc.charge, cc.fm(4)
      end

      subroutine toasciiNp(n)
      implicit none
      integer n
      write(*, '("#  # of high energy particles")')
      write(*, *)   n
      write(*, '(a)')
     * '# where code subc chg  Eenergy(GeV)'//
     * '      x       y(cm)            T(n)'
      end


      subroutine toasciiHE(he)
!      print high energy particles
      implicit none
#include "Zprivate.h"
      type(ob)::he

      write(*, '(i5,3i5, g14.4, 2g12.3, f12.3)')
     * he.where, he.code, he.subcode, he.charge,
     * he.erg, he.x/100., he.y/100., he.atime
      end


