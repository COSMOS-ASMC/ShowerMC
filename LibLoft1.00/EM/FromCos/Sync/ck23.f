!c      test ck23: K2/3
!      implicit none
!      real*8 x, y, v, kv, ck23
!      integer icon
!
!      x = 1.d-4
!      v = 2.d0/3.d0
!
!      do while (x .lt. 100.) 
!         call dbkr(x, v, kv, icon)
!         y = ck23(x)
!         write(*, *) sngl(x), sngl(kv/y), icon
!         x = x * 10.**0.1
!      enddo
!      end
!c
      real*8 function ck23(x)
!          compute K2/3(x) within 5 digit accuracy.
!         this uses the result of makepolfork23
      implicit none
      real*8 x
      real*8 pi/3.14159265/, xl
      
      if(x .lt. 1.d-1) then
         ck23 = 1.074764 * x**(-2./3)
         if(x .gt. 0.001) then
            ck23 =ck23 * ( (( 24.87215 *x -4.845528)*x -0.2195235)*x          
     *          +  1.000107 )
         endif    
      else
         ck23 = exp(-x)* sqrt(pi/2/x)
         if(x .lt. 30.) then
            xl = log(x)
            ck23 =  ck23* ( (( -0.2016376E-02*xl + 0.1788832E-01)*xl
     *          -0.5794207E-01)*xl +  1.072769 )
         endif
      endif
      end
