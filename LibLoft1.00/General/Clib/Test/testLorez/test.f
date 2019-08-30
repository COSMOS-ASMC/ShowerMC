      implicit none
#include "Zptcl.h"
#include "Zcode.h"
!
      type(ptcl):: p, po
      type(fmom)::  gb, gbi
      real(8):: pA(4), pB(4)
      real*8  g, error
      integer i, j, k
      g=1.
      call cmkptc(knuc, 0, 1, p)

!      p%fm%p(1)= 1.5d0
!      p%fm%p(2)= -8.5d1
!      p%fm%p(3)= -5.d2
!      po = p
      write(0,*) ' enter px, py, pz'
      read(*,*)  p%fm%p(1), p%fm%p(2), p%fm%p(3)
      write(0,*) ' input=', p%fm%p(1), p%fm%p(2), p%fm%p(3)
      call cpm2e(p, p)
      gb%p(1)=0.     
      gb%p(2)=0.
      do  i=1, 50
         if(g  .lt. 1.d5) then
            gb%p(3)= g * sqrt(1.d0-1.d0/g/g)
         else
            gb%p(3) = g - 0.5d0/g -1.d0/8.d0/g/g/g
         endif
         gb%p(4)= g
!           inverse
         gbi%p(1)=0.
         gbi%p(2)=0.
         gbi%p(3)=-gb%p(3)
         gbi%p(4)=g
         pA(:)= p%fm%p(:)
         do k = 1, 5
            call clorez( gb, p,  po)
!         write(*, *) ' converted po', 
!     *          po%fm%p(1), po%fm%p(2), po%fm%p(3), po%fm%p(4)
!     call clorez2( gbi, po,  po)
            call clorez( gbi, po,  po)
            pB(:)= po%fm%p(:)
            do j = 3, 4
               if( abs( pA(j) ) > 1. ) then
                  error =(pA(j) - pB(j))/pA(j)
               else
                  error =(pA(j) - pB(j))
               endif
               write(*,'(2i3, 1p, g13.6, g14.6, g14.6)')
     *          k, j, g, 
     *          pA(j),  error
            enddo
            p = po
         enddo
         g=g*10.d0**0.25d0
      enddo
      end program
