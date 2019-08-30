      implicit none
      real*8  v, eppair1, z
      real*8 e, epmaxv1, vmax

      read(*, *) z


      call  epbrem1ini(z)
      e = 1.2d-3
      do  while(e .lt. 1.d5)
         v = 0.9999d0
         vmax = epmaxv1(e)
         do while(v .gt. 0.5d0)
            if(v .lt . vmax) then
               write(*, *) sngl(v), sngl(eppair1(e, v))
            endif
            v = v - 0.005d0
         enddo
         write(*,*)
         e = e* 10.d0
      enddo
      end
    enddo
      end

