!
!      
      subroutine procLat(h,  idxf, code, layer)
!         This treats one 1-D histogram with web index (idxr, idxf)
!         
      implicit none

#include "../../../Hist/Z90histc.h"
#include "../../../Hist/Z90histo.h"
#include "../../../Hist/Z90hist1.h"
      
      type(histogram1):: h  ! input 1 D histogram
      integer idxf  ! web fai bin index (1~nf) 
      integer code  ! ptcl code
      integer layer ! at which layer (among max of 3 layers)
      integer maxsize  ! max histogram size
      integer n        ! actual histogram size
      
      integer i, j, nbin
      parameter (maxsize=2000)
      real*8 x(maxsize), y(maxsize)
      integer icon
      integer  kwhistxy

      save


      nbin= kwhistxy(h, x, y, maxsize)
      write(*,'(3i4)')  idxf, code, layer
      do i=1,nbin
         write(*,'(1p2E11.3)') x(i), y(i)
      enddo
      write(*, '("0 0")')
      end
