!      real*8 x
!      do i = 1, 1000000
!         call ksampRSA(x)
!         write(*,*) x
!      enddo
!      end
      subroutine ksampRSA(costheta)
      implicit none
!      sample cos from  (1+cos^2)dcos
      real*8  costheta  ! output  cos ofsampled angle  
      real*8  u1, u2, u3

      call rndc(u1)
      if(u1 .lt. 0.75d0) then
         call rndc(u2)
         costheta = 2.0*u2 -1.
      else
!         second term. take max of 3 u's.
         call rndc(u1)
         call rndc(u2)
         call rndc(u3)
         costheta = max(u1, u2, u3)
         call rndc(u1)
         if(u1 .lt. 0.5) then
            costheta = -costheta
         endif
      endif
      end
   
                                                                              
