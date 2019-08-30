      subroutine test1(cc)
      integer cc(0:2)
      real abc
      integer xx
      integer*2  yy
      common /fff/  abc(0:10), xx(-1:2), yy(2)
      
      integer i

      do i = 0,2
         write(*,*) cc(i) 
      enddo

      do i = 0,2
         abc(i) = i*2
         xx(i-1) = -i
      enddo
      yy(1) = 123
      call test2
      end
      subroutine test2
      real abc
      integer xx
      integer*2  yy
      common /fff/  abc(0:10), xx(-1:2), yy(2)
      integer i 
      do i = 0, 2
         write(*,*)  abc(i)
         write(*,*)  xx(i-1)
      enddo
      end



