      subroutine ksortinv(idx, n)
      implicit none
!        invert sorted idx order
      integer n
      integer idx(n)
      integer i, j

      do i = 1, n/2
         j =idx(i)
         idx(i) = idx(n-i+1)
         idx(n-i+1) = j
      enddo
      end
