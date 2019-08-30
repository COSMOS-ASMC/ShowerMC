!        get s
      real*8 rho, s, z, a, e, v, smigb, x, smigp
      integer i
      real*8 er/1.e-3/

      a = 14.5
      z = 7.25
      e = 1.e19/1.e9
      x = 1.e4
      do while(x .lt. 1.5e22/1.e9 * 1.e-3)
         rho = x/e
         call zpart(z, a, rho)
         v  = 1.e-4
         do while(v .lt. 1.0) 
!            s = smigb(v, e, s, er)
            s = smigp(v, e, s, er)
     
            write(*, *) sngl(x*1.e3), sngl(v), sngl(s)
            v = v * 10.**0.02
         enddo
         write(*,*)
         x = x * sqrt(10.d0)
      enddo
      end
