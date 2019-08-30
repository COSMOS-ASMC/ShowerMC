!      real*8 xmin, xmax
!      call kdmachmnmx(xmin, xmax)
!      write(*,*) xmin, xmax
!      end

!     ********************************
      subroutine kdmachmnmx(xmin, xmax)
!     ********************************
      implicit none
!           returns machine min and max of double 
!           precision floating numbers: 
!       This is approx one, but will be usable.
!
      real*8 xmin   ! machine min. true one will be xmin/1.1 or so.
      real*8 xmax   ! machine max. true one will be xmax*2 or so.
      real*8 d1,  d2
      integer i
      character*25 aaa

      d2 = 1.d-70
      do i = 1, 10000
         d1 = d2/5.0
         write(aaa,'(g25.20)') d1
         if(index(aaa, '.0') .ne. 0 .or.
     *      index(aaa, '0. ') .ne. 0 )  goto 10
         d2 = d1
      enddo
 10   continue
      xmin = d2
      d1 = -log10(xmin)
      if(d1 .gt. 310.d0) then
!          sun, hp. 
         xmax = 1.d0/(xmin*10.5d16) ! for safety against overflow.
      else
!          next, alpha
         xmax = 1.d0/(xmin*4.d0)
      endif
      end
      
