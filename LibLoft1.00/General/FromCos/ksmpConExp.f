!      implicit none
!      integer n,  i, j
!
!      real*8 sumx,  sum, xa(500)
!
!      read(*,*) n, sumx
!
!      do i = 1, 10000
!         call ksmpConExp(n, sumx, xa)
!         sum = 0
!         do j= 1, n
!            write(*,*) sngl(xa(j))
!            sum = sum + xa(j)
!         enddo
!         if( abs(sum-sumx) .gt. 1.d-5 )  then
!            write(0,*)  sum, sumx
!            stop
!         endif
!      enddo
!      end
       subroutine ksmpConExp(n, sumx, xa)
!
!         This samples  x of n particles so that sum of x
!         becomes sumx.  The pre-assumed one particle distribution of
!         each particle is propotional to exp(-b xi/xc) dxi
!         but there is a constraint of delta(sum(xi) - sumx).
!       sampling is can be done by using (sum-xi)^(m-2) dxi
!       for givne sum = sum(xi) (i=1, m)
!      n should be b* sumx/xc;   For xc=0.1 and b=10, n=sumx*100.
!       

      implicit none
      integer n     ! input  number of particles to be generated
      real*8 sumx   ! input  sum of x
      real*8 xa(n)  ! output.  x of n particles generated 
      
      real*8 x, u, xk
      integer j
      
      x = sumx
      do j = 1, n-1
         call rndc(u)
         xk = x*(1.0- u**(1.0/(n-j)))
         xa(j) = xk
         x = x - xk
      enddo
      xa(n) = x
      end


