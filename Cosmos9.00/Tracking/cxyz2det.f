!       cxyz2det:  xyz to det
!       cdet2xyz:  det to xyz
!       cxyz2detD: xyz to det for Direction cos.
!       cdet2xyzD: det to xyz for Direction cos.
!
      subroutine cxyz2det(det, a, b)
!          convert coord value in the "xyz" system into "det" system.
!        det: /coord/  input. detector coord in "xyz"
!          a: /coord/  input. coord in 'xyz'
!          b: /coord/  outupt. coord in 'det'
!
      implicit none

#include  "Zcoord.h"
#include  "Zobs.h"
#include  "Zpos.h"
#include  "Zmagfield.h"
#include  "Zobsv.h"

      type(coord)::a
      type(coord)::b 
      type(coord)::det

!      record /coord/ tempa
      real(8):: temp(3)
      temp(:) = a%r(:) - det%r(:)
      b%r(:) = matmul( Txyz2det(:,:), temp(:))
      b%sys = 'det'
!      real*8 leng
!       coord in "xyz"   from the origin of det.

!       tempa.r(1) = a.r(1) -  det.r(1)
!       tempa.r(2) = a.r(2) -  det.r(2)
!       tempa.r(3) = a.r(3) -  det.r(3)
!          make it direction cos.
!       call c3DV2DDCos(tempa, tempa, leng)
!          maket it direction cos in "det" system.
!       call ciTransVectZx(1, DetZaxis, DetXaxis, tempa, b)
!          to 3D vectors
!      
!       b.r(1) = b.r(1) *leng
!       b.r(2) = b.r(2) *leng
!       b.r(3) = b.r(3) *leng
!       b.sys = 'det'
       end
!      *******************************
       subroutine cdet2xyz(det, a, b)
!      *******************************
       implicit none

#include  "Zcoord.h"
#include  "Zobs.h"
#include  "Zpos.h"
#include  "Zmagfield.h"
#include  "Zobsv.h"

!
! why              type(coord)::  
        type(coord)::a
        type(coord)::b
        type(coord)::det
        real(8)::temp(3)
        temp(:) = matmul( Tdet2xyz(:,:), a%r(:))
        b%r(:) = temp(:) + det%r(:)
        b%sys = 'xyz'
!        record /coord/ tempa
!        real*8 leng
!           to direction cos in 'det'
!        call c3DV2DDCos(a, tempa, leng)
!           to direction cos in 'zyz'
!        call ctransVectZx(1, DetZaxis, DetXaxis,tempa, b)
!           to xyz sys
!        b.r(1) = b.r(1) *leng + det.r(1)
!        b.r(2) = b.r(2) *leng + det.r(2)
!        b.r(3) = b.r(3) *leng + det.r(3)
!        b.sys = 'xyz'
        end
!      -----------------------------------
        subroutine  cxyz2detD(a, b)
       implicit none

#include  "Zcoord.h"
#include  "Zobs.h"
#include  "Zpos.h"
#include  "Zmagfield.h"
#include  "Zobsv.h"


!              
       type(coord)::a
       type(coord)::b
       real::temp(3)
       temp(:) =matmul(Txyz2det(:,:), a%r(:)) 
       b%r(:) = temp(:)
       b%sys = 'det'
!       call ciTransVectZx(1, DetZaxis, DetXaxis, a, b)
!       b.sys = 'det'
       end
!      -----------------------------------------------------
       subroutine cdet2xyzD(a, b)
       implicit none
#include  "Zcoord.h"
#include  "Zobs.h"
#include  "Zpos.h"
#include  "Zmagfield.h"
#include  "Zobsv.h"

        type(coord)::a
        type(coord)::b
        real(8):: temp(3)
        temp(:) = matmul( Tdet2xyz(:,:), a%r(:))
        b%r(:) = temp(:)
!        call ctransVectZx(1, DetZaxis, DetXaxis, a, b)       
        end
