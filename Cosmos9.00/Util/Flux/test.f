!        generate isotropic point on a sphere and extract
!        vertical ones. What dose it look like ?
!        should be uniform. in  the circle
      implicit none
      real*8  cs, sn, sint, cost, x, y, z, cosx, u

      integer i

      do i = 1, 1000000

         call kcossn(cs, sn)
         call rndc(cost)
         sint = sqrt(1.d0- cost**2)

         x = sint* cs
         y = sint* sn
         z = cost 

         call rndc(cosx)
         if(cosx .gt. 0.99d0) then
            call rndc(u)
            if(cost .gt. u) then
               write(*,*)  sngl(x), sngl(y)
            endif
         endif
      enddo
      end
