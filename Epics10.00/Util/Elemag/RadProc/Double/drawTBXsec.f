!
      implicit none
      real*8 rho,  z, a, vmin, ans
      real*8 e
      common /landu1/ e
      real*8 v(100)
      integer ne, nrho

      z = 7.25
      a = 14.5
      rho = 1.d-3

      nrho = 0
      do while (rho .gt. 1.e-9)
         nrho = nrho + 1
         call zpart( z, a,  rho)
         e = 1.e15/1.e9
         ne = 0
         do  while(e .lt. 1.1e22/1.e9)
            vmin=max(1.d-3/e, 1.d-4)
            call totcb(vmin, 1.d0, ans)
!            write(*,*) sngl(e), sngl(ans)
            ne = ne + 1
            v(ne) = log10(ans)
            e = e* 10.**.1
         enddo
         write(*, *)
         rho = rho/10.**0.5
         call kmkDataStm2b(v, ne, nrho, 'bremxs', 'f9.4',9)
      enddo

      end

 
)
      enddo

      end

 
