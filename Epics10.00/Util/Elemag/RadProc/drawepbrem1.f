      implicit none
      real*8  v, epbrem1, z
      real*8 e, epmaxv1

      read(*, *) z


      call  epbrem1ini(z)
      e = .711d-3
      do  while(e .lt. 1.d5)
         v = epmaxv1(e)*0.99999d0

         do while(v .gt. 1.d-6)
            write(*, *) sngl(v), sngl(epbrem1(e, v)*v), sngl(e)
            v = v/10.d0**0.05d0
         enddo
         write(*,*)
         e = e* 10.d0
      enddo
      end
  enddo
      end

