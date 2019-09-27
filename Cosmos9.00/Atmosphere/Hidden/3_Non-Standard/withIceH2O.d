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
# height(m)   	Temp(K)	 density(kg/m3)    media
#-----------------------------------------------------------------
	-400	290.75	1.2790     ice
	0	288.15	1.2250       
	1950    275.4	1.0715
	2000    275.15	1.068        
	3.e3	268.7   0.9889     Air
	5.e3    255.6   0.8315     H2O*0.01
	5.1e3   255.0   0.8236     Air        
	6.e3	249.2	0.7528
	11.1e3	216.65	0.35932
	20.0e3	216.65	8.891e-2
	32.2e3	228.756	1.3145e-2
	47.4e3	270.65	1.4187e-3
	51.0e3	270.65	9.0696e-4
	72.0e3	214.263	6.2374e-5
	86.0e3	186.87	6.958e-6
	91.0e3	186.87	2.860e-6
	110.0e3	240.0	9.708e-8
	130.0e3	469.27	8.152e-9
	160.0e3	696.29	1.233e-9
	250.0e3	941.33	6.073e-11
	300.0e3	976.01	1.916e-11
	400.0e3 995.83	2.803e-12
	500.0e3	999.24	5.215e-13
        600.0e3	999.85	1.137e-13
	700.0e3	999.97	3.070e-14
	800.0e3	999.99	1.136e-14
	900.0e3	1000.0	5.759e-15
	1.0e6	1000.0	3.561e-15
