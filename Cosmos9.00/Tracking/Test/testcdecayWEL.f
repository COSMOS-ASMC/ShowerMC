      implicit none
      real*8 p, x, g
      integer i, imax

      write(0,*) ' Enter p, g, imax'

      read(*,*)  p, g, imax

      do i=1, imax
         call cdecayWEL(p, g, x)
         write(*,*) sngl(x)
      enddo
      end
