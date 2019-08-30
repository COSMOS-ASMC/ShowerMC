#
#
#        This gives standard atmosphere data. The scale height is approximated
#   by a  number of stright lines as a function of height. The data gives 
#   height,  temperatur, etc at  each nodal point.
#   At height lesss than 11km, special formula is used, however.
#
#         The scale height, H, is expressed by H = H0 + a(z-z0)
#                                                = kT/mg        
#   in each region.
#   We neglect height dependence of gravitational accelleration g,  
#   and the average mass of  air molecules, m.   
#   Since the table below gives T(z)= T0 + b(z-z0) at the nodal points,
#   we can first get b,  and then a by a = dH/dz = k/mg * b. 
#   H at a nodal point, z,   is obtained as H(z) =kT(z)/mg.
#   The density is given by
#              rho = rho0 * (1+ a(z-z0)/H(z0))**(-1-1/a)    (a != 0)
#                  = rho0 * exp(- (z-z0)/H)          (a =0; hence H is const)
#
#   The gramage between  given heights, z1 and z2  is by
#
#              d = d0 *(fd(z2) - fd(z1))  where
#         
#            fd(z) = (1+ a(z-z0)/H(z0))**(-1/a)                       (a != 0)
#
#                  =  exp(-(z-z0)/H )                                 (a = 0)
#
#  and d0 = rho0*H(z0).
#
#
#
# height(m)   	Temp(K)	Press(hP) 	density(kg/m3)  
#-----------------------------------------------------------------
	-400	290.75	1062.2		1.2790
	0	288.15	1013.25		1.2250
	3.e3	268.659	7.0121e2	0.9889
	6.e3	249.187	4.7217e2	0.7528
	11.1e3	216.65	223.46		0.35932
	20.0e3	216.65	55.293		8.891e-2
	32.2e3	228.756	8.6314		1.3145e-2
	47.4e3	270.65	1.1022		1.4187e-3
	51.0e3	270.65	7.0458e-1	9.0696e-4
	72.0e3	214.263	3.8362e-2	6.2374e-5
	86.0e3	186.87	3.7388e-3	6.958e-6
	91.0e3	186.87	1.5381e-3	2.860e-6
	110.0e3	240.0	7.1042e-5	9.708e-8
	130.0e3	469.27	1.2505e-5	8.152e-9
	160.0e3	696.29	3.0359e-6	1.233e-9
	250.0e3	941.33	2.4767e-7	6.073e-11
	300.0e3	976.01	8.7704e-8	1.916e-11
	400.0e3 995.83	1.4518e-8	2.803e-12
	500.0e3	999.24	3.0236e-9	5.215e-13
        600.0e3	999.85	8.2130e-10	1.137e-13
	700.0e3	999.97	3.1908e-10	3.070e-14
	800.0e3	999.99	1.7036e-10	1.136e-14
	900.0e3	1000.0	1.0873e-10	5.759e-15	
	1.0e6	1000.0	7.5138e-11	3.561e-15
