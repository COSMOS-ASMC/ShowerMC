      implicit none
      real al183z, e, ccz, emass, bcoef, fz, z333
      common / bpcom /  al183z,e,ccz,emass,bcoef,fz,z333
      real  v, fbrem, z
      real  vmaxv, zpart

      read(*, *) z


      v = zpart(z)
      e = .711e-3
      do  while(e .lt. 1.e5)
         v = vmaxv(e)*0.99999d0

         do while(v .gt. 1.d-6)
            write(*, *) v, fbrem(v)*v, e
            v = v/10.d0**0.05d0
         enddo
         write(*,*)
         e = e* 10.d0
      enddo
      end






ddo
      end













