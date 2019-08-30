!
      implicit none
      real*8 rho, v, fbrem, z, a
      real*8 e
      common /landu1/ e
      read(*, *) rho
      z = 7.25
      a = 14.5

      call zpart( z, a,  rho)
      e = 1.e15/1.e9
      do  while(e .lt. 1.1e22/1.e9)
         v = 0.9999
         do while(v .gt. 1.e-6)
            write(*, *) sngl(v), sngl(fbrem(v)*v)
            v = v/10.**0.05
         enddo
         write(*,*)
         e = e* 10.
      enddo
      end
