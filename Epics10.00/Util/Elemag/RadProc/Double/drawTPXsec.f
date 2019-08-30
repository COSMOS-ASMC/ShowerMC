!
      implicit none
      real*8 rho,  z, a,  ans
      real*8 e, cvh2den
      common /landu1/ e
      real*8 v(100), h, x0
      integer ne, nrho
      read(*,*) e
      e = e/1.e9
      x0 = 364.
      z = 7.25
      a = 14.5
      h = 10.d3
      do while (h .lt. 1001.d3)
         write(*,*) '#  h=', h/1.d3, ' km'
         rho = cvh2den(h)
!         write(*,*)'# rho=', rho*1000./1.d6 , 'g/cm3'
         call zpart( z, a,  rho*1.d-3)  ! rho shouud be g/cm3
!         do  while(e .lt. 1.1e21/1.e9)
            call totcp(0.d0, 1.d0, ans)
            ans = ans / x0 * rho *1.d3  !  prob/km
            if(ans .gt. 1.d-37) then
               write(*,*) sngl(e), sngl(h), sngl(1.d0/ans)
            endif
!            e = e* 10.**.1
!         enddo
!         write(*, *)
          h = h * 10.d0**0.01d0
      enddo

      end

 

      enddo

      end

 
