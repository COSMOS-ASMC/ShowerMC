!          cxyz2prim:  xyz to primary system coord. conversion
!          cxyz2primD: xyz to primary system for Direction cos. 
!
      subroutine cxyz2prim(base, a, b)
      implicit none
!
#include  "Zcoord.h"
#include  "Zobs.h"
#include  "Zpos.h"
#include  "Zmagfield.h"
#include  "Zobsv.h"

!                                  origin of primary system is the
!                                  origin of this deepest detector
      type(coord)::base ! base detector coord. whose origin is the
                          ! origin of 1ry system.
      type(coord)::a   ! input.  coord. in 'xyz'
      type(coord)::b   ! output. transformed coord. in 'prim'

      real(8):: temp(3)
      temp(:) = a%r(:) - base%r(:)
      b%r(:) = matmul(Txyz2prim(:,:), temp)
      b%sys = 'prim'
!  old one (<= v7.637)
!      record /coord/ tempa
!       coord in "xyz"   from the origin of det.
!       tempa.r(1) = a.r(1) - base.r(1)
!       tempa.r(2) = a.r(2) - base.r(2)
!       tempa.r(3) = a.r(3) - base.r(3)
!       call cscalerProd(tempa, Xprimary, b.r(1))
!       call cscalerProd(tempa, Yprimary, b.r(2))
!       call cscalerProd(tempa, Zprimary, b.r(3))
!       b.sys = 'prim'
       end
!      -----------------------------------
       subroutine cxyz2primD(a, b)
      implicit none
!
#include  "Zcoord.h"
#include  "Zobs.h"
#include  "Zpos.h"
#include  "Zmagfield.h"
#include  "Zobsv.h"

       type(coord)::a
       type(coord)::b
       real(8):: temp(3)
       temp(:) = matmul(Txyz2prim(:,:), a%r(:))
       b%r(:) = temp(:)
       b%sys = 'prim'
!        old version (<= 7.637)
!       record /coord/temp
       
!          b can be a

!       call cscalerProd(a, Xprimary, temp.r(1))
!       call cscalerProd(a, Yprimary, temp.r(2))
!       call cscalerProd(a, Zprimary, temp.r(3))
!       b = temp
!       b.sys = 'prim'
       end

      subroutine cprim2xyz(base, a, b)
      implicit none
!           inverse of xyz2prim
#include  "Zcoord.h"
#include  "Zobs.h"
#include  "Zpos.h"
#include  "Zmagfield.h"
#include  "Zobsv.h"


!                                  origin of base detector
      type(coord)::base ! base detector coord. whose origin is the
                          ! origin of 1ry system.
      type(coord)::a   ! input.  coord. in 'prim'
      type(coord)::b   ! output. transformed coord. in 'xyz'

      real(8):: temp(3)
      temp(:) = matmul(Tprim2xyz(:,:), a%r(:))
      b%r(:) = temp(:) + base%r(:)
      b%sys = 'xyz'
      end
!      -----------------------------------
       subroutine cprim2xyzD(a, b)
      implicit none
!
#include  "Zcoord.h"
#include  "Zobs.h"
#include  "Zpos.h"
#include  "Zmagfield.h"
#include  "Zobsv.h"

       type(coord)::a
       type(coord)::b
       real(8):: temp(3)
       temp(:) = matmul(Tprim2xyz(:,:), a%r(:))
       b%r(:) = temp(:)
       b%sys = 'xyz'
       end
