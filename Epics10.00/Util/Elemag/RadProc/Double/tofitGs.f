!          to make approximate formula for G(s).
!
!     since    G(s) --> 12pis^2,
!        we output  G(s)/ (12pis^2/(1+ 12pis^2)) and approximate the
!        results with a polynomial.
!
      implicit none

      real*8 s, gmigdl
      real*8 const, pi, temp
      
      pi =  asin(1.0d0)*2 
      write(*, *) pi
      const = pi * 12.0d0
      do s = 0.d0, 1.2d0, 0.01d0
         if(s .ne. 0.) then
            temp = gmigdl(s, 1.d-5)/( const*s**2/(1.d0 + const*s**2))
         else
            temp = 1.
         endif
         write(*, *) sngl(s), sngl(temp)
      enddo
      end

      
