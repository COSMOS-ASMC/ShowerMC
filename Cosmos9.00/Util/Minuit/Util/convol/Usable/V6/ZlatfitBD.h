         block data
         include "Zlatfit.h"
c             fitting region         
	 data x1(1)/0.03/, x2(1)/0.5/
	 data x1(2)/0.3/, x2(2)/2.0/
	 data x1(3)/1.0/, x2(3)/15.0/
	 data x1(4)/7.0/, x2(4)/60.0/
c	 data x1(4)/7.0/, x2(4)/30.0/  N.G
c             drawing region
         data drx1(1)/0.01/, drx2(1)/0.3/
         data drx1(2)/0.3/, drx2(2)/1.0/
         data drx1(3)/1.0/, drx2(3)/10./
         data drx1(4)/10./, drx2(4)/100./
c          a/r^(b+cr**pw)     
         data param(1,1)/0.1d0/
	 data param(2,1)/1.0d0/  
         data param(3,1)/0.01d0/
         data param(4,1)/0.5d0/

         data param(1,2)/0.1d0/
	 data param(2,2)/1.0d0/  
         data param(3,2)/0.01d0/
         data param(4,2)/0.5d0/

         data param(1,3)/1.d0/
	 data param(2,3)/2.0d0/  
         data param(3,3)/0.1d0/
         data param(4,3)/0.5d0/

         data param(1,4)/100.d0/
	 data param(2,4)/1.0d0/  
         data param(3,4)/1.0d0/
         data param(4,4)/0.50d0/

c              g,e,m,h  (*,g,e,m,u1) * is for param 
c              g
         data low(1,1)/1.d-5/, low(2,1)/-8.0d0/, low(3,1)/0.0d0/
         data up(1,1)/1.0d4/, up(2,1)/7.0d0/,  up(3,1)/3.0d0/   
         data low(4,1)/0.0d0/
         data up(4,1)/1.0d0/
c              e 
         data low(1,2)/1.d-5/, low(2,2)/-8.0d0/, low(3,2)/0.0d0/
         data up(1,2)/1.0d4/, up(2,2)/7.0d0/,  up(3,2)/3.0d0/   
         data low(4,2)/0.0d0/
         data up(4,2)/1.0d0/
c              mu 
         data low(1,3)/1.d-5/, low(2,3)/-8.0d0/, low(3,3)/0.0d0/
         data up(1,3)/1.0d4/, up(2,3)/7.0d0/,  up(3,3)/3.0d0/   
         data low(4,3)/0.0d0/
         data up(4,3)/1.0d0/

c              h
         data low(1,4)/1.d-5/, low(2,4)/-8.0d0/, low(3,4)/0.0d0/
         data up(1,4)/1.0d4/, up(2,4)/7.0d0/,  up(3,4)/3.0d0/   
         data low(4,4)/0.0d0/
         data up(4,4)/1.0d0/

	 end

