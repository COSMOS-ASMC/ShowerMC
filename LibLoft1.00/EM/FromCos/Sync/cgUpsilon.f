!	real*8 x, y
!	real*8 gupsilon
!	x = 0.001
!	do while (x .lt. 500.) 
!	    y = gupsilon(x)
!	    write(*, *) sngl(x), sngl(y)
!	    x = x * 10.**0.1
!        enddo	
!	end
!    ******************************************************************
!         g(upsilon) for total energy loss function of synchrotron
!    ******************************************************************
!
	real*8 function cgUpsilon(u)
	implicit none
	real*8 u  ! input. upsilon

	real*8 z

	if(u .lt. 0.005) then
		cgUpsilon = u ** 2 * (1. -5.953*u)
	elseif(u .lt. 100.) then
!                 polynomial approx. within 5 % error. roughly
		z = log(u)
                cgUpsilon = 
     *          exp(((((-0.6305280E-04 * z -0.1159820E-03)*z 
     *           +0.9467500E-02 )*z-0.6861282E-01)*z  +1.080512)*z
     *           -1.770289)       
	else
	      cgUpsilon = 0.5563*u**0.66666666
	endif	
	end



