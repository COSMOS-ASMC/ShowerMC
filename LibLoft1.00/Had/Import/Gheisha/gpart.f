
      subroutine gpart
      implicit none
c.
c.
c
      real *4 emass
      real *4 emmu
      real *4 pmass
c      data emass/0.51099906e-03/
      data emass/0.511e-03/
      data emmu/0.105658389/
      data pmass/0.93827231/
c
c
c.
c.    -------------------------------------------------------------------
c.
      call gspart( 1,'gamma$       ',1,0.      , 0.,1.000000e+15)
      call gspart( 2,'positron$    ',2,emass   , 1.,1.000000e+15)
      call gspart( 3,'electron$    ',2,emass   ,-1.,1.000000e+15)
      call gspart( 4,'neutrino$    ',3,0.      , 0.,1.000000e+15)
      call gspart( 5,'muon +$      ',5,emmu    , 1.,2.197030e-06)
      call gspart( 6,'muon -$      ',5,emmu    ,-1.,2.197030e-06)
      call gspart( 7,'pion 0$      ',3,0.134973, 0.,0.840000e-16)
      call gspart( 8,'pion +$      ',4,0.139567, 1.,2.603000e-08)
      call gspart( 9,'pion -$      ',4,0.139567,-1.,2.603000e-08)
      call gspart(10,'kaon 0 long$ ',3,0.49767 , 0.,5.183000e-08)
      call gspart(11,'kaon +$      ',4,0.493646, 1.,1.237100e-08)
      call gspart(12,'kaon -$      ',4,0.493646,-1.,1.237100e-08)
      call gspart(13,'neutron$     ',3,0.939566, 0.,8.960000e+02)
      call gspart(14,'proton$      ',4,pmass   , 1.,1.000000e+15)
      call gspart(15,'antiproton$  ',4,pmass   ,-1.,1.000000e+15)
      call gspart(16,'kaon 0 short$',3,0.49767 , 0.,8.922000e-11)
      call gspart(17,'eta$         ',3,0.5488  , 0.,7.479742e-19)
      call gspart(18,'lambda$      ',3,1.11563 , 0.,2.631000e-10)
      call gspart(19,'sigma +$     ',4,1.18937 , 1.,0.800000e-10)
      call gspart(20,'sigma 0$     ',3,1.19255 , 0.,7.400000e-20)
      call gspart(21,'sigma -$     ',4,1.19743 ,-1.,1.479000e-10)
      call gspart(22,'xi 0$        ',3,1.3149  , 0.,2.900000e-10)
      call gspart(23,'xi -$        ',4,1.32132 ,-1.,1.639000e-10)
      call gspart(24,'omega -$     ',4,1.67243 ,-1.,0.822000e-10)
      call gspart(25,'antineutron$ ',3,0.939566, 0.,8.960000e+02)
      call gspart(26,'antilambda$  ',3,1.11563 , 0.,2.631000e-10)
      call gspart(27,'antisigma -$ ',4,1.18937 ,-1.,0.800000e-10)
      call gspart(28,'antisigma 0$ ',3,1.19255 , 0.,7.400000e-20)
      call gspart(29,'antisigma +$ ',4,1.19743 , 1.,1.479000e-10)
      call gspart(30,'antixi 0$    ',3,1.3149  , 0.,2.900000e-10)
      call gspart(31,'antixi +$    ',4,1.32132 , 1.,1.639000e-10)
      call gspart(32,'antiomega +$ ',4,1.67243 , 1.,0.822000e-10)
      call gspart(33,'tau +$       ',4,1.7841  , 1.,3.040000e-13)
      call gspart(34,'tau -$       ',4,1.7841  ,-1.,3.040000e-13)
      call gspart(35,'d +$         ',4,1.8693  , 1.,1.062000e-12)
      call gspart(36,'d -$         ',4,1.8693  ,-1.,1.062000e-12)
      call gspart(37,'d 0$         ',3,1.8645  , 0.,4.280000e-13)
      call gspart(38,'anti d 0$    ',3,1.8645  , 0.,4.280000e-13)
      call gspart(39,'ds+$         ',4,1.9693  , 1.,4.360000e-13)
      call gspart(40,'ds-$         ',4,1.9693  ,-1.,4.360000e-13)
      call gspart(41,'lambda c +$  ',4,2.2849  , 1.,1.790000e-13)
      call gspart(42,'w +$         ',4,81.000  , 1.,9.400000e-26)
      call gspart(43,'w -$         ',4,81.000  ,-1.,9.400000e-26)
      call gspart(44,'z 0$         ',3,92.400  , 0.,7.740000e-26)
      call gspart(45,'deuteron$    ',4,1.875613,+1.,1.000000e+15)
      call gspart(46,'triton$      ',4,2.81448 ,+1.,1.000000e+15)
      call gspart(47,'alpha$       ',4,3.727417,+2.,1.000000e+15)
      call gspart(48,'geantino$    ',6,0.      , 0.,1.000000e+15)
      call gspart(200,'pseudo$     ',6,0.      , 0.,1.000000e+15)
c
  99  return
      end
