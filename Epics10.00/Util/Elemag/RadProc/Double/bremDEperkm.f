!
      implicit none
      real*8 rho, z, a
      real*8 e, fe, ans, x
      external fe
      common /landu1/ e

      z = 7.25
      a = 14.5
      rho = 1.e-5
      e = 1.e15/1.e9
      x= rho * e	
      e = 1.e19/1.e9      ! standard
      do while (x  .lt. 1.e-3 * 1.e22/1.e9)
	 rho = x/e
         call zpart( z, a,  rho)
         call k16pGaussLeg(fe, 1.d-4, 1.d0, 16, ans)
         write(*, *) sngl(x), sngl(ans)
         x = x* 10.**.1
      enddo
      end
      real*8 function fe(x)
      real*8 x, fbrem
      fe =  fbrem(x)* x
      end
 fbrem(x)* x
      end
