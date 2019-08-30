         block data
         include "Zfit.h"
         
	 data x1(1)/0.03/, x2(1)/2.0/
c          for LDD fitting
c	 data x1(2)/2.0/, x2(2)/30.0/
c          for FDD fitting
	 data x1(2)/1.0/, x2(2)/50.0/

         data drx1(1)/0.01/, drx2(1)/1.5/
         data drx1(2)/1.5/, drx2(2)/100./

         data param(1,1)/150.d0/
	 data param(2,1)/1.05d0/  
         data param(3,1)/0.001d0/
         data param(1,2)/150.d0/
	 data param(2,2)/1.05d0/  
         data param(3,2)/0.001d0/
c              g,e,m  (*,1) * is for param 
         data low(1,1)/1.d0/, low(2,1)/0.3d0/, low(3,1)/-0.2d0/
         data up(1,1)/1000.d0/, up(2,1)/3.2d0/,  up(3,1)/0.4d0/   
c                h  (*,2)  * is for param 
         data low(1,2)/10.d0/, low(2,2)/0.3d0/, low(3,2)/-0.2d0/
	 data up(1,2)/5000.d0/,  up(2,2)/3.2d0/, up(3,2)/0.4d0/   

	 end

