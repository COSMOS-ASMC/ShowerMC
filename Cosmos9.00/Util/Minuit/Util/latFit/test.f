      real*8  x(4,3,4)
      data x(1,1,:)/-1.,2.,-3.,4/
      do i =1, 4
         write(0,*) x(1,1,i)
      enddo
      end
