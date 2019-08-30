!      program testBremI1
!      implicit none
!      real*8 upsilon, xc, y, cmBremI1
!      xc = 0.
!c      do while (xc .lt. 4.e-3)
!         upsilon = 100.
!         do while (upsilon .gt. 1.e-6)
!            y = cmBremI1(upsilon, xc)
!            write(*, *) sngl(upsilon), sngl(y)
!            upsilon = upsilon / 10.**0.1
!            if(y .lt. 1.e-5)  upsilon = 1.e-9
!         enddo
!         write(*, *)
!c      
!c      enddo
!      end
!     **********************************************************
	real*8 function cmBremI1(upsilon, xc)
	implicit none
        real*8 upsilon  ! input. Upsilon value
        real*8 xc   !   input.  fix the cut off, below which the synchroron 
                    !           is treated as energy loss only.
                    ! one of   3.1627 x 10-3, 10-3, ... 10-6, or 0.0
	            ! if xc is not one of these, neareset one
                    ! is  chosen.

!  some values;
!     U     31.6  10      1      .1            .01  
!    I1           1.99   3.6      4.85          5.19
!    I2      .19   .227   .15     1.76e-2    3.53e-4 
!      
!        Psi0:  the integral  x= Xc to 1 of 4 k(z)/z * 1/(2+ 3Uz)^2 dz
!        where      z = 2/3/U*  x/(1-x).
!
!         Xc (=Eg/Ee) must be one
!        3.162x10-3, 10-3, 3.162x10-4, 10-4,... 10-6, 0.
!      
!        The 'cross-section' for emitting fractional gamma ray energy x ~ x+dx,
!        in the unit distance  is given by 
!        p(Ee, x, U)dx = root(3)/2/pi (SyncConvR)/Ee* U * (1-x)/x K(2zeta)dx
!        (/meter, if Ee is in GeV). ( U is   E/m * B/Bc ).
!        This is integrated from x = Xc to 1 to give the total 'crosssection'
!        The following tbl containes  11 coefficients for the polynomial 
!        approximation of the part of the integral as a function 
!        of Upsilon. (in log-log) 
!
!        polynomial value times root(3)/2/pi SyncConvR/Ee *U
!        is the total  cross-section (= probability/m).
!        
!        tbl(*, i) i=1, 8 contains table of 11 coefficients of polynaomials
!        Integral (x =Xc to 1 ) of   (1-x)/x K( 2/3/U  x/(1-x))/ 4 dx
!        where U= upsilon value.   Xc= 3.1627 x 10-3 to 10-6 step /3.1627
!        (i=1 is for 3.162x10-3, i=2 for 10-3, ...).  
!
	real*8  psi, x, xc0coef(8)
	integer ixc, i
!
	real*8  tbl(11, 8)
	real*8  psi3_3(11), psi_3(11), psi3_4(11), psi_4(11)
	real*8  psi3_5(11), psi_5(11), psi3_6(11), psi_6(11)
	equivalence (psi3_3(1), tbl(1,1)) 
	equivalence (psi_3(1),  tbl(1,2)) 
	equivalence (psi3_4(1), tbl(1,3)) 
	equivalence (psi_4(1),  tbl(1,4)) 
	equivalence (psi3_5(1), tbl(1,5)) 
	equivalence (psi_5(1),  tbl(1,6)) 
	equivalence (psi3_6(1), tbl(1,7)) 
	equivalence (psi_6(1),  tbl(1,8)) 


	

          data (psi3_3(i), i=  1,   11)/
     1  -0.36519777    , -0.16437488    , -0.42191193E-01,
     2   0.45068651E-02,  0.37406044E-03, -0.11694408E-03,
     3  -0.21981366E-04,  0.50330259E-05,  0.56416748E-06,
     4  -0.59483304E-07, -0.10075934E-07                                
     * /   
          data (psi_3(i), i=  1,   11)/
     1  -0.27621962    , -0.17806948    , -0.36641237E-01,
     2   0.34820468E-02,  0.23404694E-03, -0.11645243E-03,
     3  -0.11166093E-05,  0.51644288E-05,  0.66054175E-07,
     4  -0.94357089E-07, -0.73799132E-08                                
     * /   
          data (psi3_4(i), i=  1,   11)/
     1  -0.21981296    , -0.18463515    , -0.33139695E-01,
     2   0.22214970E-02,  0.96012771E-04, -0.47834498E-04,
     3   0.15951965E-04,  0.32866129E-05, -0.42901554E-06,
     4  -0.10473240E-06, -0.53301833E-08                                
     * /   
          data (psi_4(i), i=  1,   11)/
     1  -0.18243710    , -0.18712415    , -0.31745835E-01,
     2   0.90498800E-03,  0.12296262E-03,  0.46807678E-04,
     3   0.22577025E-04,  0.42917421E-06, -0.75768289E-06,
     4  -0.97835715E-07, -0.37327138E-08                                
     * /   

          data (psi3_5(i), i=  1,   11)/
     1  -0.15645188    , -0.18758306    , -0.32142658E-01,
     2  -0.27076403E-03,  0.30784765E-03,  0.14090033E-03,
     3   0.21276677E-04, -0.25452860E-05, -0.96816286E-06,
     4  -0.87676147E-07, -0.26691205E-08                                
     * /   
          data (psi_5(i), i=  1,   11)/
     1  -0.13939747    , -0.19218187    , -0.32177166E-01,
     2   0.43072437E-03,  0.50884199E-03,  0.78017665E-04,
     3   0.33119508E-05, -0.20977899E-05, -0.47432646E-06,
     4  -0.37755362E-07, -0.10680872E-08                                
     * /   
          data (psi3_6(i), i=  1,   11)/
     1  -0.12864838    , -0.19609994    , -0.31528914E-01,
     2   0.11811036E-02,  0.56385559E-03,  0.84414004E-05,
     3  -0.73982020E-05, -0.91893847E-06, -0.85709610E-07,
     4  -0.61186589E-08, -0.19199417E-09                                
     * /   
          data (psi_6(i), i=  1,   11)/
     1  -0.12129256    , -0.19786006    , -0.31152534E-01,
     2   0.13900376E-02,  0.58057791E-03, -0.13423486E-04,
     3  -0.10636349E-04, -0.54026892E-06,  0.34550944E-07,
     4   0.35276755E-08,  0.72269160E-10                                
     * /   
!               for Xc=0.  for U>3.e-3.
         data ( xc0coef(i), i=  1,   8)/
     1  -0.10535316    , -0.20019514    , -0.30586898E-01,
     2  0.13258432E-02,  0.58454445E-03, -0.13505923E-04,
     3  -0.97296726E-05, -0.53504562E-06                                
     * /   
   
	if(xc .eq. 0.) then
	    if(upsilon .lt. 3.e-3) then
	       psi =  0.269
            else
	       psi = xc0coef(8)
	       x = log(upsilon)
	       do i = 7, 1, -1
		  psi = x * psi + xc0coef(i)
	       enddo	
	    endif
	 else
!              
	    ixc =  log10(3.1627e-3/ xc)*2 + 1.1
	    if(abs (10.**(-2 - ixc/2.0) / xc  -1.) .gt. 0.05) then
	       call cerrorMsg('input xc is invalid: cmBremI1', 0)
	    endif	
	    if(ixc .lt. 1 .or. ixc .gt. 8) then
	       call cerrorMsg
     *        ('input xc is too large or too small: cmBremI1', 0)
	    endif
!	
	    
	    psi = tbl(11, ixc)
	    x = log(upsilon)
	    do i = 10, 1, -1
	       psi = x * psi + tbl(i, ixc)
	    enddo	
	 endif

	cmBremI1 =4* exp(psi)
	end

