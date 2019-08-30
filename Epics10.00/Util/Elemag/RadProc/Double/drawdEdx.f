!
      implicit none
      real*8 rho, z, a, x0, h
      real*8 e, fe, ans, x, cvh2den
      external fe
      common /landu1/ e
      x0 =364.
      z = 7.25
      a = 14.5
      h = 1.d3
      read(*,*) e
      e = e/1.e9
      do while (h .lt. 200.d3)
	 rho= cvh2den(h)
         call zpart( z, a,  rho*1.d-3)  ! must be in g/cm^3
         call k16pGaussLeg(fe, 0.d0, 1.d0, 16, ans)
         write(*, *) e*1d9,  h/1.d3, e*ans/x0 *rho*1.d9*1.d3   ! ev/km
         h = h *10.**0.01
      enddo
      end
      real*8 function fe(x)
      real*8 x, fbrem
      fe =  fbrem(x)* x
      end
fbrem(x)* x
      end
