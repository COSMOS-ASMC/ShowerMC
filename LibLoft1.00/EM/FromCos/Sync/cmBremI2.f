!      implicit none
!      real*8 u, s
!      real*8 cmBremI2
!      u = 1.e-5
!      do while (u .lt. 300.) 
!         s = cmBremI2(u)
!         write(*,*) u, s
!         u = u * 10.**.1
!      enddo
!      end
	real*8 function cmBremI2(u)
!          This program give approximate integral value of
!          v^2 2zeta K(2zeta) dv from v = 0 to 1.
!          Original integaral is done by inteSpec1.f and
!          result was approximated by polynomials.
!
	implicit none
	real*8  u  ! input. Upsilon
!       function value  ! output. integral value 

        real*8 lu, ans
	if(u .lt. 3.16d-5) then
!                ans = ans * u**2 *  3.9257
                ans = u**2 *  3.9257     ! 2012.Feb.4
	elseif(u .lt.  1.584d-2) then
		lu = log(u)
		ans = ((-0.5219002E-02 *lu - 0.1067265)* lu     
     *          -0.7421499 )*lu  -0.7647935     
                ans = ans * u**2 *  3.9257
        elseif(u .lt. 1.) then
		lu = log(u)
		ans = ((0.1696074E-01*lu + 0.1197285)*lu
     *         +  0.1439367E-01)*lu  +   0.5595450E-01 
                ans = ans * u**2 *  3.9257
	elseif(u .lt. 300.) then
		lu = log(u)

		ans = (( 0.2566103E-02*lu -0.3215703E-01 )*lu
     *          +  0.9286942E-01)*lu + 0.1500733     

        else
	   write(0,*) ' error input to function cmBremI2(u)'
	   write(0,*) ' u=',u
	   call cerrorMsg(' oh my god', 0)
        endif	
        cmBremI2 = ans
	end
