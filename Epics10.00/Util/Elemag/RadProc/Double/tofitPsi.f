!          to make approximate formula for Psi(s).
!
!     since    Psi(s) --> 6s
!        we output  Psi(s)/( 6s/(1+6s) )
!
      implicit none

      real*8 s, psimig
      real*8 temp
      
      do s = 0.d0, 1.1d0, 0.01d0
         if(s .ne. 0.) then
            temp = psimig(s, 1.d-5)/( 6.d0*s/(1.d0 + 6.d0*s))
         else
            temp = 1.
         endif
         write(*, *) sngl(s), sngl(temp)
      enddo
      end

      
