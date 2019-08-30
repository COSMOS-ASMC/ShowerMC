      subroutine sub2(a,b)
      implicit none
      real(8)::a, b

      real(8)::x, y
      common /abc/ x, y(10)

      x = a
      y(10)=b
      end subroutine sub2

      subroutine sub3
      real(8)::x, y
      common /abc/ x, y(10)
      write(0,*) x, y(10)
      end subroutine sub3

