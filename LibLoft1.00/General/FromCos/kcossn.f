!        testing kcossn
!      include  'rnd.f'
!          ---------------
!      implicit none
!      integer i
!      real*8 cs, sn
!      do i=1, 50000
!          call kcossn(cs, sn)
!         write(*, *) sngl(asin(sn))
!      enddo
!      end
      subroutine kcossn(cs,sn)
!     to generate cos(phy),sin(phy), where phy is uniform in (0,2pi).
!      reuired subprogram.  rndc 
      implicit none
      real*8 cs, sn
!
      real*8 u, v, a, b, c
      c = 5.d0

      do  while (c .gt. 1.d0)
          call rndc(u)
          call rndc(v)
          v=v+v-1.
          a=u*u
          b=v*v
          c=a+b
      enddo    

      cs=(a-b)/c
      sn=2.*u*v/c
      end
