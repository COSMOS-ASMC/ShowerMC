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
c          a/r^(b+cr**pw)   (abcpw,code  ,region)
c           gamma rgn1
         data param(1,1,1)/0.1d0/
	 data param(2,1,1)/0.0d0/  
         data param(3,1,1)/1.0d0/
         data param(4,1,1)/0.5d0/
c           g  rgn 2
         data param(1,1,2)/0.1d0/
	 data param(2,1,2)/0.0d0/  
         data param(3,1,2)/1.0d0/
         data param(4,1,2)/0.5d0/
c            g rgn 3
         data param(1,1,3)/0.1d0/
	 data param(2,1,3)/0.0d0/  
         data param(3,1,3)/1.0d0/
         data param(4,1,3)/0.5d0/
c              g rgn 4
         data param(1,1,4)/0.1d0/
	 data param(2,1,4)/2.0d0/  
         data param(3,1,4)/1.0d0/
         data param(4,1,4)/0.25d0/
c            elec rgn 1
         data param(1,2,1)/0.1d0/
	 data param(2,2,1)/0.0d0/  
         data param(3,2,1)/2.5d0/
         data param(4,2,1)/0.5d0/
c            e rgn 2
         data param(1,2,2)/0.1d0/
	 data param(2,2,2)/0.0d0/  
         data param(3,2,2)/2.5d0/
         data param(4,2,2)/0.1d0/
c            e rgn 3         
         data param(1,2,3)/0.1d0/
	 data param(2,2,3)/0.0d0/  
         data param(3,2,3)/1.0d0/
         data param(4,2,3)/0.4d0/
c            rgn 4
         data param(1,2,4)/0.1d0/
	 data param(2,2,4)/0.0d0/  
         data param(3,2,4)/1.0d0/
         data param(4,2,4)/0.35d0/
c             mu rgn1
         data param(1,3,1)/1.d0/
	 data param(2,3,1)/0.0d0/  
         data param(3,3,1)/0.7d0/
         data param(4,3,1)/0.5d0/
c            mu rgn 2
         data param(1,3,2)/1.d0/
	 data param(2,3,2)/0.0d0/  
         data param(3,3,2)/1.0d0/
         data param(4,3,2)/0.5d0/
c           mu  rgn 3
         data param(1,3,3)/1.d0/
	 data param(2,3,3)/0.0d0/  
         data param(3,3,3)/1.0d0/
         data param(4,3,3)/0.5d0/
c            mu  rgn 4
         data param(1,3,4)/1.d0/
	 data param(2,3,4)/0.0d0/  
         data param(3,3,4)/1.0d0/
         data param(4,3,4)/0.5d0/
c            hadron rgn 1
         data param(1,4,1)/100.d0/
	 data param(2,4,1)/0.0d0/  
         data param(3,4,1)/1.0d0/
         data param(4,4,1)/0.20d0/
c               rgn 2
         data param(1,4,2)/100.d0/
	 data param(2,4,2)/1.0d0/  
         data param(3,4,2)/1.0d0/
         data param(4,4,2)/0.50d0/
c              rgn 3
         data param(1,4,3)/100.d0/
	 data param(2,4,3)/0.50d0/  
         data param(3,4,3)/0.0d0/
         data param(4,4,3)/0.8d0/
c             rgn 4
         data param(1,4,4)/100.d0/
	 data param(2,4,4)/0.0d0/  
         data param(3,4,4)/1.0d0/
         data param(4,4,4)/0.50d0/

c              g,e,m,h  (*,g,e,m,u1) * is for param 
c              g
         data  up(1,1,:)/4*1.0d4/
         data low(1,1,:)/4*1.d-5/
         data  up(2,1,:)/1.d0,  1.d0,    1.3d0,3.5d0/
         data low(2,1,:)/-1.d0, -1.25d0, -1.d0, 1.4d0/
         data  up(3,1,:)/2.7d0, 2.5d0, 2.d0,1.5d0/
         data low(3,1,:)/0.7d0, 0.2d0, 0.18d0, 0.2d0/
         data  up(4,1,:)/0.8d0, 0.7d0, 0.8d0, 0.38d0/
         data low(4,1,:)/0.03d0, 0.07d0, 0.2d0, 0.02d0/
c                e
         data  up(1,2,:)/4*1.0d4/
         data low(1,2,:)/4*1.d-5/
         data  up(2,2,:)/0.8d0, 1.2d0, 2.0d0, 3.5d0/
         data low(2,2,:)/-0.8d0, -0.9d0, -1.d0, -3.d0/
         data  up(3,2,:)/3.d0, 2.92d0, 2.8d0, 2.7d0/
         data low(3,2,:)/2.d0, 2.d0, 0.18d0, 0.1d0/
         data  up(4,2,:)/0.6d0, 0.17d0, 0.6d0, 0.6d0/
         data low(4,2,:)/0.12d0, 0.03d0, 0.05d0, 0.05d0/
c               mu
         data  up(1,3,:)/4*1.0d4/
         data low(1,3,:)/4*1.d-5/
         data  up(2,3,:)/0.25d0, 0.35d0, 0.7d0, 1.d0/
         data low(2,3,:)/-1.5d0, -2.d0, -1.5d0, -2.d0/
         data  up(3,3,:)/1.d0, 2.5d0, 2.d0, 1.5d0/
         data low(3,3,:)/0.5d0, 0.1d0, 0.05d0, 0.1d0/
         data  up(4,3,:)/0.9d0, 0.95d0, 0.9d0, 0.8d0/
         data low(4,3,:)/0.08d0, 0.06d0, 0.13d0, 0.2d0/
c                h
         data  up(1,4,:)/4*1.0d4/
         data low(1,4,:)/4*1.d-5/
         data  up(2,4,:)/0.9d0, 0.7d0, 0.5d0, 1.d0/
         data low(2,4,:)/-1.d0, 0.2d0,-1.d0,-1.0d0/
         data  up(3,4,:)/2.4d0, 0.1d0, 2.8d0, 1.5d0/
         data low(3,4,:)/0.5d0, 0.d0, 0.1d0, 0.05d0/
         data  up(4,4,:)/0.5d0, 0.99d0, 0.7d0, 0.6d0/
         data low(4,4,:)/0.d0, 0.7d0,  0.01d0,0.05d0/

	 end

