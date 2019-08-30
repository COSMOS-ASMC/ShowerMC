!         test ctXai
!	program testcsynkTXai
!	implicit none
!	real*8 x, ctXai, y
!
!	x = 1.e-2
!	do while (x .lt. 4000.) 
!		y = ctXai(x)
!		write(*, *) sngl(x),sngl(y)
!	        x = x * 10.**0.1
!	enddo	
!	end
!       *************************************************************
!
!           ctXai: T(X)
!
!         magnetic pair creation function.
!         This is a polynomial fit to the T(x)=0.16/x *K1/3(2/3x/)**2 
!          (Eq.3.3d) of Erber.
!        Note that Table VI seems half of this value.
!
!       *************************************************************
	real*8 function ctXai(x)
	implicit none
	real*8 x      ! input.  x = Xai = hv/m * H/Hc/2
	
	real*8 z

	if(x .lt. 0.25) then
		ctXai = 0.46 * exp(-4.d0/3.d0/x)
        elseif(x .lt. 1000.) then	
           z = log(x)
           ctXai =exp( (((((((((-0.1055848E-06*z + 0.3326052E-05)*z 
     *            -0.4298832E-04)*z +  0.3164778E-03 )*z
     *            -0.1791347E-02)*z +  0.1028788E-01 )*z 
     *            -0.5477092E-01)*z +  0.2251815 )*z    
     *            -0.6962186)*z     + 1.227734)*z -2.443217)       
	
	else
	   ctXai = 0.60 * x**(-0.3333333333)
	endif	
	end


