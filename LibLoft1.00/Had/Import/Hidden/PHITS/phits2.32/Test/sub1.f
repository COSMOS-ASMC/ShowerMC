      subroutine sub1
      implicit none
      real(8)::a, b

      real(8)::x, y
      common /abc/ x, y(10)

      a = y(1)**2
      b = x
      call sub2(a, b)
      end subroutine sub1

