!
      implicit none
      real*8 rho, v, fpair, z, a
      real*8 e
      common /landu1/ e
      read(*, *) rho
      z = 7.25
      a = 14.5

      call zpart( z, a,  rho)
      e = 1.e17/1.e9
      do  while(e .lt. 1.1e22/1.e9)
         v = 0.00
         do while(v .lt. .50001)
            write(*, *) sngl(v), sngl(fpair(v))
            v = v + 0.005
         enddo
         write(*,*)
         e = e* 10.
      enddo
      end
   enddo
      end
