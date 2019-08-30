      real*8 cvh2den, h
      h = 10.d3
      do while (h .lt. 5000.d3) 
         write(*,*) sngl(h),sngl( cvh2den(h) )
         h = h *10.**.1d0
      enddo
      end
