*cmz :  3.14/13 22/06/90  08.37.41  by  rene brun
*-- author :    rene brun   23/05/90
      subroutine hadini
c
      common /hadrun/ runtes,eftes
      common /redver/ irii(17),ikii(17),ieii(17),
     +        thresh(268)
*      integer*2 ich,ibar,k1,k2
      common/abltis/am(110),ga(110),tau(110),ich(110),
     +              ibar(110),k1(110),k2(110)
      common/reacn/siin(296),wkn(5184)
      common /reach/umo(296),plabf(296),siinh(296),wkh(5184),
     +              nrk(2,268),nure(30,2)
*      integer*2 nzk
      common/split/ nzk(460,3),wt(460)
*
      dimension iriic(17),ikiic(17),ieiic(17),nktemp(536)
c
      dimension umopi( 92),umokc( 68),umop( 39),umon( 63),umok0( 34),
     +plapi( 92),plakc( 68),plap( 39),plan( 63),plak0( 34), spikp1(
     +315),spikpu(278),spikpv(372), spikpw(278),spikpx(372),spikp4(
     +315), spikp5(187),spikp6(306),skmpel(102),spikp7(289),skmnel(
     +68),spikp8(187), spikp9(143),spikp0(169),spkpv(143), sappel(105)
     +,spikpe(399),sapnel( 84),spikpz(273), sanpel( 84),spikpf(273),
     +spkp15(187),spkp16(255), nrkpi( 164),nrkkc( 134),nrkp( 70),nrkn(
     +116),nrkk0( 52), nurec(60)
     
     
     
     
     
c
      equivalence (nktemp(  1),nrkpi(1)),(nktemp(165),nrkkc(1))
      equivalence (nktemp(299),nrkp (1)),(nktemp(369),nrkn (1))
      equivalence (nktemp(485),nrkk0(1))
      dimension umoc(296),plabfc(296)
      equivalence (umoc(  1),umopi(1)),(umoc( 93),umokc(1))
      equivalence (umoc(161),umop (1)),(umoc(200),umon (1))
      equivalence (umoc(263),umok0(1))
      equivalence (plabfc(  1),plapi(1)),(plabfc( 93),plakc(1))
      equivalence (plabfc(161),plap (1)),(plabfc(200),plan (1))
      equivalence (plabfc(263),plak0(1))
c
      dimension amz(16),gaz(16),tauz(16),ichz(16),ibarz(16),
     +  k1z(16),k2z(16),nzk1(153),nzk2(153),nzk3(153),wtz(153)
c
c*** reaction channel cross section data
c****pi- p data
c**** pi+ n data
      data plapi/ 0.,.3,. 5,.6,.7,.8,.9,.95,1.,1.15,1.3,1.5,1.6,1.8,2.,
     +2.3,2.5,2.8,3.,3.5,4., 0.,.285,.4,.45,.5,.6,.7,.75,.8,.85,.9,1.,
     +1.15,1.3,1.5,1.6, 1.8,2.,2.3,2.5,2.8,3.,3.5,4.,4.5, 0.,.285,.4,.
     +45,.5,.6,.7,.75,.8,.85,.9,1.,1.15,1.3,1.5,1.6, 1.8,2.,2.3,2.5,2.
     +8,3.,3.5,4.,4.5, 0.,.3,. 5,.6,.7,.8,.9,.95,1.,1.15,1.3,1.5,1.6,1.
     +8,2., 2.3,2.5,2.8,3.,3.5,4./
     
     
     
      data plakc/ 0.,0.58,0.80,1.01,1.23,1.45,1.68,1.94,2.18,2.42,2.68,
     +2.96,3.24, 3.51,3.84,4.16,4.49, 0.,0.58,0.80,1.01,1.23,1.45,1.68,
     +1.94,2.18,2.42,2.68,2.96,3.24, 3.51,3.84,4.16,4.49, 0.,0.58,0.80,
     +1.01,1.23,1.45,1.68,1.94,2.18,2.42,2.68,2.96,3.24, 3.51,3.84,4.
     +16,4.49, 0.,0.58,0.80,1.01,1.23,1.45,1.68,1.94,2.18,2.42,2.68,2.
     +96,3.24, 3.51,3.84,4.16,4.49/
     
     
     
      data plak0/ 0.,0.58,0.80,1.01,1.23,1.45,1.68,1.94,2.18,2.42,2.68,
     +2.96,3.24, 3.51,3.84,4.16,4.49, 0.,0.58,0.80,1.01,1.23,1.45,1.68,
     +1.94,2.18,2.42,2.68,2.96,3.24, 3.51,3.84,4.16,4.49/
     
     
c                 pp   pn   np   nn
      data plap/ 0.,1.06,1.34,1.63,1.92,2.2,2.5,2.8,3.1,3.43,3.75,4.07,
     +4.43, 0.,1.06,1.34,1.63,1.92,2.2,2.5,2.8,3.1,3.43,3.75,4.07,4.43,
     +0.,1.06,1.34,1.63,1.92,2.2,2.5,2.8,3.1,3.43,3.75,4.07,4.43/
     
c    app   apn   anp   ann
      data plan / 0.0,0.001,0.1,.2,.3,.4,.5,.6, .74,1.06,1.34,1.63,1.
     +92,2.2,2.5,2.8,3.1,3.43,3.75,4.07,4.43 , 0.0,0.001,0.1,.2,.3,.4,.
     +5,.6, .74,1.06,1.34,1.63,1.92,2.2,2.5,2.8,3.1,3.43,3.75,4.07,4.
     +43 , 0.0,0.001,0.1,.2,.3,.4,.5,.6, .74,1.06,1.34,1.63,1.92,2.2,2.
     +5,2.8,3.1,3.43,3.75,4.07,4.43 /
     
     
      data umopi/ 1.08,1.233,1.302,1.369,1.496,1.557,1.615,1.6435 ,1.
     +672,1.753,1.831,1.930,1.978,2.071,2.159,2.286,2.366,2.482,2.56,2.
     +735,2.90, 1.08,1.222,1.302,1.3365,1.369,1.434,1.496,1.527,1.557 ,
     +1.586,1.615,1.672,1.753,1.831,1.930,1.978,2.071,2.159,2.286,2.
     +366,2.482,2.560,2.735,2.90,3.06, 1.08,1.222,1.302,1.3365,1.369,1.
     +434,1.496,1.527,1.557 ,1.586,1.615,1.672,1.753,1.831,1.930,1.978,
     +2.071,2.159,2.286,2.366,2.482,2.560,2.735,2.90,3.06, 1.08,1.233,
     +1.302,1.369,1.496,1.557,1.615,1.6435 ,1.672,1.753,1.831,1.930,1.
     +978,2.071,2.159,2.286,2.366,2.482,2.56,2.735,2.90/
     
     
     
     
      data umokc/ 1.44, 1.598,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,
     +2.7,2.8,2.9,3.0, 3.1,1.44, 1.598,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,
     +2.5,2.6,2.7,2.8,2.9,3.0, 3.1,1.44, 1.598,1.7,1.8,1.9,2.0,2.1,2.2,
     +2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0, 3.1,1.44, 1.598,1.7,1.8,1.9,2.0,
     +2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0, 3.1/
     
     
     
     
      data umok0/ 1.44, 1.598,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,
     +2.7,2.8,2.9,3.0, 3.1,1.44, 1.598,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,
     +2.5,2.6,2.7,2.8,2.9,3.0, 3.1/
     
     
c                 pp   pn   np   nn
      data umop/ 1.88,2.102,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.,3.1,3.2,
     +1.88,2.102,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.,3.1,3.2, 1.88,2.
     +102,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.,3.1,3.2/
     
c    app   apn   anp   ann
      data umon / 1.877,1.87701,1.879,1.887,1.9,1.917,1.938,1.962,
     +2.,2.102,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.,3.1,3.2, 1.877,1.
     +87701,1.879,1.887,1.9,1.917,1.938,1.962, 2.,2.102,2.2,2.3,2.4,2.
     +5,2.6,2.7,2.8,2.9,3.,3.1,3.2, 1.877,1.87701,1.879,1.887,1.9,1.
     +917,1.938,1.962, 2.,2.102,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.,3.1,
     +3.2/
     
c*** reaction cahnnel state particles
      data nrkpi/ 13,1,15,21,81,0, 13,54,23,53,13,63,13,58,23,57,13,65,
     +1,32,53,31,54,32,53,33,53,35, 63,32, 13,8,23,1,17,15,21,24,22,15,
     +82,0, 61,0,13,55,23,54,14,53,13,64, 23,63,13,59,23,58,14,57,13,
     +66,23,65,1,31,8,32,1,33,1,35,54,31,55, 32,54,33,53,34,54,35,
     +14,1,23,8,17,24,20,15,22,24,83,0, 62,0,14,54,23,55,13,56,14,63,
     +23,64,14,58,23,59,13,60,14,65,23,66,8,31,1,34,8,33,8,35,55,31,54,
     +34,55,33,56,32,55,35, 14,8,24,20,84,0, 14,55,23,56,14,64,14,59,
     +23,60,14,66,8,34,56,31,55,34,56,33,56,35, 64,34 /
     
     
     
     
     
     
     
     
      data nrkkc/ 15,1,89,0, 24,53,15,54,1,36,1,40,1,44,36,63,15,63,45,
     +53,44,54, 15,8,24,1,91,0, 24,54,15,55,8,36,1,37,8,40,1,41,8,44,1,
     +45,36,64,37,63,15,64,24,63,45,54,44,55,93,0, 16,1,25,8, 17,23,21,
     +14, 20,13,22,23, 90,0, 38,1,39,8,16,54,25,55,1,42,8,43,16,63,25,
     +64,39,64,38,63,46,54, 47,55,8,47,1,46,52,0,51,0, 16,8,17,14,20,
     +23,22,14,92,0, 8,38,16,55,25,56,8,42,16,64,38,64,46,55,47,56,8,
     +46,94,0 /
     
     
     
     
     
c
c   k0 p   k0 n   ak0 p   ak/ n
c
      data nrkk0/ 24,8, 106,0,15,56,24,55,37,8,41,8,45,8,37,64,24,64,
     +44,56,45,55, 25,1,17,13, 22,13,21,23, 107,0,39,1,25,54,16,53,43,
     +1,25,63,39,63,47,54,46,53,47,1,103,0/
     
     
c   pp  pn   np   nn
      data nrkp/1,1,85,0, 8,53,1,54,1,63,8,57,1,58,54,54,53,55,63,54,
     +64,53, 1,8,86,0, 8,54,1,55,8,63,1,64,8,58,1,59,64,54,63,55,54,55,
     +53,56,77,0, 8,8, 95,0,8,55,1,56,8,64,8,59,1,60,55,55,54,56,64,55,
     +63,56/
     
     
c     app   apn   anp   ann
      data nrkn/ 1,2,17,18,15,16,8,9,13,14,99,0,87,0, 1,68,8,69,2,54,9,
     +55,102,0, 2,63,9,64,1,75,8,76,53,67,54,68, 55,69,56,70,63,68,64,
     +69,75,54,76,55, 2,8,18,20, 16,24,14,23, 101,0,88,0, 2,55,9,56,1,
     +67,8,68,2,64,8,75,2,59,8,72,68,55,67,54,69,56, 1,9,18,21,15,25,
     +13,23,100,0, 96,0,2,53,9,54,1,69,8,70,1,76,9,63,1,73,9,58,55,70,
     +53,68,54,69/
     
     
     
c*** channel cross section
      data spikp1/ 0.,300.,40.,20.,13.,8.5,8.,9.5,12.,14.,15.5,20.,17.,
     +13.,10., 9.,8.5,8.,7.8,7.3,6.7, 9*0.,.23,.35,.7,.52,.4,.3,.2,.15,
     +.13,.11, .09,.07, 0.,0.033,.8,1.35,1.35,.5,15*0., 3*0.,.00,0.80,
     +2.2,3.6,4.6,4.7,3.5,2.4,1.8,1.4,.75,.47,.25,.13,.08, 6*0.,
     +0.,1.2,3.3,5.4,6.9,7.3,5.3,3.6,2.7,2.2,1.1,.73,.4,.22,.12, 9*0.,
     +.0,0.,2.0,4.4,6.8,9.9,7.9,6.0,3.8,2.5,2.,1.4,1.,.6,.35,10*0.,
     +.25,.55,.75,1.25,1.9,2.,1.8,1.5,1.25,1.,.8,6*0., 4*0.,.4,.85,1.1,
     +1.85,2.8,3.,2.7,2.2,1.85,1.5,1.2,6*0., 6*0.,.5,1.2,1.7,3.4,5.2,6.
     +4,6.1,5.6,5.2,6*0., 2*0.,.0,1.0,3.3,5.2,4.45,3.6,2.75,1.9,1.65,1.
     +3,.95,.6,.45,6*0., 3*0.,.0,.45,1.4,1.5,1.1,.85,.5,.3,.2,.15,8*0.,
     +5*0.,.0,.0,.6,.8,.95,.8,.7,.6,.5,.4,6*0., 5*0.,.0,.00,.85,1.2,1.
     +4,1.2,1.05,.9,.7,.55,6*0., 5*0.,.0,.00,1.,1.5,3.5,4.15,3.7,2.7,2.
     +3,1.75,6*0., 10*0.,.5,2.0,3.3,5.4,7. /
     
     
     
     
     
c**** pi+ n data
      data spikpu/ 0.,25.,13.,11.,10.5,14.,20.,20.,16.,14.,19.,28.,17.
     +5,13.5,12., 10.5,10.,10.,9.5,9.,8.,7.5,7.,6.5,6., 0.,48.,19.,15.,
     +11.5,10.,8.,6.5,5.5,4.8,4.2,7.5,3.4,2.5,2.5,2.1, 1.4,1.,.8,.6,.
     +46,.3,.2,.15,.13, 11*0.,.95,.65,.48,.35,.2,.18,.17,.16,.15,.1,.
     +09,.065,.05,.04, 12*0.,.2,.25,.25,.2,.1,.08,.06,.045,.03,.02,.01,
     +.005,.003, 12*0.,.3,.24,.18,.15,.13,.12,.11,.1,.09,.08,.05,.04,.
     +03, 0.,0.16,.7,1.3,3.1,4.5,2.,18*0., 3* .0,0.0,0.0,4.0,11.0,11.4,
     +10.3,7.5,6.8,4.75,2.5,1.5,.9,.55,.35,13*0.,.1,.34,.5,.8,1.1,2.25,
     +3.3,2.3,1.6,.95,.45,.28,.15,10*0., 2*0.,.17,.64,1.,1.5,2.1,4.25,
     +6.2,4.4,3.,1.8,.9,.53,.28,10*0., 2*0.,.25,.82,1.3,1.9,2.8,5.5,8.,
     +5.7,3.9,2.35,1.15,.69,.37,10*0., 7*0.,.00,.34,1.5,3.47,5.87,6.23,
     +4.27,2.6,1.,.6,.3,.15,6*0./
     
     
      data spikpv/ 7*0.,.00,.16,.75,1.73,2.93,3.12,2.13,1.3,.5,.3,.15,.
     +08,6*0., 10*0.,.2,.6,.92,2.4,4.9,6.25,5.25,3.5,2.15,1.4,1.,.7,13*
     +0.,.13,.4,.62,1.6,3.27,4.17,3.5,2.33,1.43,.93,.66,.47,13*0.,.07,.
     +2,.31,.8,1.63,2.08,1.75,1.17,.72,.47,.34,.23,17*0.,.33,1.,1.8,2.
     +67,5.33,6., 5.53,5.,17*0.,.17,.5,.9,1.83,2.67,3.0,2.77,2.5,3*0.,
     +0.,.0,0.,1.0,3.3,2.8,2.5,2.3,1.8,1.5,1.1,.8,.7,.55,.3,10*0.,
     +9*0.,.1,.4,1.,1.4,2.2,2.5,2.2,1.65,1.35,1.1,.8,.6,.4,3*0.,
     +9*0.,.15,.6,1.5,2.1,3.3,3.8,3.3,2.45,2.05,1.65,1.2,.9,.6,3*0.,
     +9*0.,.10,.2,.5,.7,1.3,1.55,1.9,1.8,1.55,1.35,1.15,.95,.7,3*0.,
     +9*0.,.00,.2,.5,.7,1.3,1.55,1.9,1.8,1.55,1.35,1.15,.95,.7,3*0.,
     +14*0.,.2,.5,.85,2.,2.15,2.05,1.75,1.,17*0.,.13,.33,.57,1.33,1.43,
     +1.36,1.17,.67,17*0.,.07,.17,.28,.67,.72,.69,.58,.33,3*0., 14*0.,.
     +4,.7,1.,1.6,1.8,2.3,1.9,1.7 /
     
c**** pi- p data
      data spikpw/ 0.,25.,13.,11.,10.5,14.,20.,20.,16.,14.,19.,28.,17.
     +5,13.5,12., 10.5,10.,10.,9.5,9.,8.,7.5,7.,6.5,6., 0.,48.,19.,15.,
     +11.5,10.,8.,6.5,5.5,4.8,4.2,7.5,3.4,2.5,2.5,2.1, 1.4,1.,.8,.6,.
     +46,.3,.2,.15,.13, 11*0.,.95,.65,.48,.35,.2,.18,.17,.16,.15,.1,.
     +09,.065,.05,.04, 12*0.,.2,.25,.25,.2,.1,.08,.06,.045,.03,.02,.01,
     +.005,.003, 12*0.,.3,.24,.18,.15,.13,.12,.11,.1,.09,.08,.05,.04,.
     +03, 0.,0.16,.7,1.3,3.1,4.5,2.,18*0., 3* .0,0.0,0.0,4.0,11.0,11.4,
     +10.3,7.5,6.8,4.75,2.5,1.5,.9,.55,.35,13*0.,.1,.34,.5,.8,1.1,2.25,
     +3.3,2.3,1.6,.95,.45,.28,.15,10*0., 2*0.,.17,.64,1.,1.5,2.1,4.25,
     +6.2,4.4,3.,1.8,.9,.53,.28,10*0., 2*0.,.25,.82,1.3,1.9,2.8,5.5,8.,
     +5.7,3.9,2.35,1.15,.69,.37,10*0., 7*0.,.00,.34,1.5,3.47,5.87,6.23,
     +4.27,2.6,1.,.6,.3,.15,6*0./
     
     
      data spikpx/ 7*0.,.00,.16,.75,1.73,2.93,3.12,2.13,1.3,.5,.3,.15,.
     +08,6*0., 10*0.,.2,.6,.92,2.4,4.9,6.25,5.25,3.5,2.15,1.4,1.,.7,13*
     +0.,.13,.4,.62,1.6,3.27,4.17,3.5,2.33,1.43,.93,.66,.47,13*0.,.07,.
     +2,.31,.8,1.63,2.08,1.75,1.17,.72,.47,.34,.23,17*0.,.33,1.,1.8,2.
     +67,5.33,6., 5.53,5.,17*0.,.17,.5,.9,1.83,2.67,3.0,2.77,2.5,3*0.,
     +0.,.0,0.,1.0,3.3,2.8,2.5,2.3,1.8,1.5,1.1,.8,.7,.55,.3,10*0.,
     +9*0.,.1,.4,1.,1.4,2.2,2.5,2.2,1.65,1.35,1.1,.8,.6,.4,3*0.,
     +9*0.,.15,.6,1.5,2.1,3.3,3.8,3.3,2.45,2.05,1.65,1.2,.9,.6,3*0.,
     +9*0.,.10,.2,.5,.7,1.3,1.55,1.9,1.8,1.55,1.35,1.15,.95,.7,3*0.,
     +9*0.,.00,.2,.5,.7,1.3,1.55,1.9,1.8,1.55,1.35,1.15,.95,.7,3*0.,
     +14*0.,.2,.5,.85,2.,2.15,2.05,1.75,1.,17*0.,.13,.33,.57,1.33,1.43,
     +1.36,1.17,.67,17*0.,.07,.17,.28,.67,.72,.69,.58,.33,3*0., 14*0.,.
     +4,.7,1.,1.6,1.8,2.3,1.9,1.7 /
     
c**** pi- n data
      data spikp4 / 0.,300.,40.,20.,13.,8.5,8.,9.5,12.,14.,15.5,20.,17.
     +,13.,10., 9.,8.5,8.,7.8,7.3,6.7, 9*0.,.23,.35,.7,.52,.4,.3,.2,.
     +15,.13,.11, .09,.07, 0.,0.033,.8,1.35,1.35,.5,15*0., 3*0.,.00,0.
     +80,2.2,3.6,4.6,4.7,3.5,2.4,1.8,1.4,.75,.47,.25,.13,.08, 6*0.,
     +0.,1.2,3.3,5.4,6.9,7.3,5.3,3.6,2.7,2.2,1.1,.73,.4,.22,.12, 9*0.,
     +.0,0.,2.0,4.4,6.8,9.9,7.9,6.0,3.8,2.5,2.,1.4,1.,.6,.35,10*0.,
     +.25,.55,.75,1.25,1.9,2.,1.8,1.5,1.25,1.,.8,6*0., 4*0.,.4,.85,1.1,
     +1.85,2.8,3.,2.7,2.2,1.85,1.5,1.2,6*0., 6*0.,.5,1.2,1.7,3.4,5.2,6.
     +4,6.1,5.6,5.2,6*0., 2*0.,.0,1.0,3.3,5.2,4.45,3.6,2.75,1.9,1.65,1.
     +3,.95,.6,.45,6*0., 3*0.,.0,.45,1.4,1.5,1.1,.85,.5,.3,.2,.15,8*0.,
     +5*0.,.0,.0,.6,.8,.95,.8,.7,.6,.5,.4,6*0., 5*0.,.0,.00,.85,1.2,1.
     +4,1.2,1.05,.9,.7,.55,6*0., 5*0.,.0,.00,1.,1.5,3.5,4.15,3.7,2.7,2.
     +3,1.75,6*0., 10*0.,.5,2.0,3.3,5.4,7. /
     
     
     
     
     
c**** k+ p data
      data spikp5/ 0.,20.,14.,12.,11.5,10.,8.,7.,6.,5.5,5.3,5.,4.5,4.4,
     +3.8,3.,2.8, 0.,.5,1.15,2.,1.3,.8,.45,10*0., 3*0.,0.9,2.5,3.,2.5,
     +2.3,2.,1.7,1.5,1.2,.9,.6,.45,.21,.2, 3*0.,0.9,2.5,3.,2.5,2.3,2.,
     +1.7,1.5,1.2,.9,.6,.45,.21,.2, 4*0.,1.0,2.1,2.6,2.3,2.1,1.8,1.7,1.
     +4,1.2,1.05,.9,.66, .5, 7*0.,.3,1.,1.,.9,.7,.4,.30,.2,.00,0.,
     +9*0.,.1,1.,2.2,3.5,4.20,4.55,4.85,4.9, 10*0.,.2,.7,1.6,2.5,2.2,1.
     +71,1.6, 6*0.,1.4,3.8,5.,4.7,4.4,4.,3.5,2.85,2.35,2.01,1.8,
     +12*0.,.1,.8,2.05,3.31,3.5, 12*0.,.034,.20,.75,1.04,1.24 /
     
     
     
     
     
c**** k+ n data
      data spikp6/ 0.,6.,11.,13.,6.,5.,3.,2.2,1.5,1.2,1.,.7,.6,.5,.45,.
     +35,.3, 0.,6.,11.,13.,6.,5.,3.,2.2,1.5,1.2,1.,.7,.6,.5,.45,.35,.3,
     +0.,.5,1.3,2.8,2.3,1.6,.9,10*0., 3*0.,0.9,2.5,3.,2.5,2.3,2.,1.7,1.
     +5,1.2,.9,.6,.45,.21,.2, 3*0.,0.9,2.5,3.,2.5,2.3,2.,1.7,1.5,1.2,.
     +9,.6,.45,.21,.2, 4*0.,1.0,2.1,2.6,2.3,2.0,1.8,1.7,1.4,1.2,1.15,.
     +9,.66, .5, 4*0.,1.0,2.1,2.6,2.3,2.1,1.8,1.7,1.4,1.2,1.15,.9,.66,
     +.5, 7*0.,.3,1.,1.,.9,.7,.4,.35,.2,.00,0., 7*0.,.3,1.,1.,.9,.7,.4,
     +.35,.2,.00,0., 9*0.,.1,1.,2.4,3.5,4.25,4.55,4.85,4.9, 9*0.,.1,1.,
     +2.4,3.5,4.25,4.55,4.85,4.9, 10*0.,.2,.7,1.6,2.5,2.2,1.71,1.6,
     +10*0.,.2,.7,1.6,2.5,2.2,1.71,1.6, 6*0.,1.4,3.8,5.,4.7,4.4,4.,3.5,
     +2.85,2.35,2.01,1.8, 6*0.,1.4,3.8,5.,4.7,4.4,4.,3.5,2.85,2.35,2.
     +01,1.8, 12*0.,.1,.8,2.05,3.31,3.5, 12*0.,.034,.20,.75,1.04,1.24,
     +.0,2.5,15.,21.5,15.3,3.,1.5,10*0./
     
     
     
     
     
c**** k- p data
      data skmpel/ 0.,35.,22.,25.,17.,9.,9.5,8.,7.,6.5,6.1,5.,4.8,4.6,
     +4.45,4.3,4.2, 0.,8.,3.5,8.,3.,1.9,1.7,1.,.9,.8,.75,.5,.42,.38,.
     +34,.25,.2, 0.,3.,3.2,3.5,1.5,1.4,1.1,.6,.5,.35,.28,.25,.18,.12,.
     +1,.08,.04, 0.,8.5,2.4,1.7,1.3,1.3,1.1,.5,.4,.4,.35,.3,.28,.2,.16,
     +.13,.11, 0.,7.,4.8,1.4,1.9,.9,.4,.2,.13,.1,.08,.06,.04,.02,.015,.
     +01,.01, 0.,5.5,1.,.8,.75,.32,.2,.1,.09,.08,.065,.05,.04,.022,.
     +017,2*.01/
      data spikp7/ 0.,.56,1.46,3.16,2.01,1.28,.74,10*0., 4*0.,1.13,2.
     +61,2.91,2.58,2.35,2.02,1.91,1.57,1.35,1.29,1.01,.74, .65, 4*0.,1.
     +13,2.61,2.91,2.58,2.35,2.02,1.91,1.57,1.35,1.29,1.01,.74, .65,
     +3*0.,1.00,3.03,3.36,2.8,2.58,2.24,1.91,1.68,1.35,1.01,.67,.5,.24,
     +.23, 3*0.,1.00,3.03,3.36,2.8,2.58,2.24,1.91,1.68,1.35,1.01,.67,.
     +5,.24, .23, 7*0.,.34,1.12,1.12,1.01,.78,.45,.39,.22,.07,0.,
     +7*0.,.34,1.12,1.12,1.01,.78,.45,.39,.22,.07,0., 6*0.,1.71,4.26,5.
     +6,5.57,4.93,4.48,3.92,3.19,2.63,2.25,2., 6*0.,1.71,4.26,5.6,5.57,
     +4.93,4.48,3.92,3.19,2.63,2.25,2., 10*0.,.22,.8,.75,1.,1.3,1.5,1.
     +3, 10*0.,.22,.8,.75,1.,1.3,1.5,1.3,13*0.,.1,.3,.7,1., 13*0.,.1,.
     +3,.7,1., 9*0.,.11,1.72,2.69,3.92,4.76,5.10,5.44,5.3, 9*0.,.11,1.
     +72,2.69,3.92,4.76,5.10,5.44,5.3, 4*0.,0.00,9.2,4.7,1.9,9*0.,
     +.0,2.5,15.,21.5,15.3,3.,1.5,10*0./
     
     
     
     
     
     
c**** k- n data
      data skmnel/ 0.,4.,9.5,20.,13.,9.5,6.,4.4,3.,2.4,2.,1.4,1.2,1.,.
     +9,.7,.6, 0.,4.5,6.,5.,2.5,2.,1.7,2.1,1.9,.9,.5,.3,.24,.2,.18,.1,.
     +09, 0.,1.8,2.,1.1,.9,.5,.5,.4,.4,.2,.1,.06,.05,.04,.03,.02,.02,
     +0.,1.5,2.,.9,1.1,.4,.6,.7,.65,.3,.17,.1,.08,.07,.06,.04,.03/
     
      data spikp8/ 0.,.56,1.29,2.26,1.01,.64,.37,10*0., 4*0.,1.13,2.61,
     +2.91,2.58,2.35,2.02,1.91,1.57,1.35,1.29,1.01,.74, .65, 3*0.,1.00,
     +3.03,3.36,2.8,2.58,2.24,1.91,1.68,1.35,1.01,.67,.5,.24, .23,
     +3*0.,1.00,3.03,3.36,2.8,2.58,2.24,1.91,1.68,1.35,1.01,.67,.5,.24,
     +.23, 7*0.,.34,1.12,1.12,1.01,.78,.45,.39,.22,.07,0., 6*0.,1.71,4.
     +26,5.6,5.57,4.93,4.48,3.92,3.19,2.63,2.25,2., 10*0.,.22,.8,.75,1.
     +,1.3,1.5,1.3, 13*0.,.1,.3,.7,1., 13*0.,.1,.3,.7,1., 9*0.,.11,1.
     +72,2.69,3.92,4.76,5.10,5.44,5.3, 4*0.,0.00,9.2,4.7,1.9,9*0.
     +/
     
     
     
     
     
     
     
c****  p p data
      data spikp9/ 0.,24.,25.,27.,23.,21.,20.,19.,17.,15.5,14.,13.5,13.
     +, 0.,3.6,1.7, 10*0., .0,0.0,8.7,17.7,18.8,15.9,11.7,8.,6.,5.3,4.
     +5,3.9,3.5, .0,.0,2.8,5.8,6.2,5.1,3.8,2.7,2.1,1.8,1.5,1.3,1.1,
     +4*0.,0.0,4.6,10.2,15.1,16.9,16.5,11.,5.5,3.5, 10*0.,4.3,7.6,9.,
     +10*0.,1.7,2.6,3., 6*0.,.3,.6,1.,1.6,1.3,.8,.6, 6*0.,.7,1.2,1.8,2.
     +5,1.8,1.3,1.2, 10*0.,.6,1.4,1.7, 10*0.,1.9,4.1,5.2/
     
     
     
     
     
     
c****  p n data
      data spikp0/ 0.,24.,25.,27.,23.,21.,20.,19.,17.,15.5,14.,13.5,13.
     +, 0.,1.8,.2, 10*0., .0,.0,3.2,6.05,9.9,5.1,3.8,2.7,1.9,1.5,1.4,1.
     +3,1.1, .0,.0,3.2,6.05,9.9,5.1,3.8,2.7,1.9,1.5,1.4,1.3,1.1,
     +4*0.,0.0,4.6,10.2,15.1,16.4,15.2,11.,5.4,3.5, 4*0.,0.0,4.6,10.2,
     +15.1,16.4,15.2,11.,5.4,3.5, 10*0.,.7,5.1,8., 10*0.,.7,5.1,8.,
     +10*.0,.3,2.8,4.7, 10*.0,.3,2.8,4.7, 7*0.,1.2,2.5,3.5,6.,5.3,2.9,
     +7*0.,1.7,3.6,5.4,9.,7.6,4.2, 3*0.,.0,0.0,7.7,6.1,2.9, 5*0./
     
     
     
     
     
     
     
c   nn - data
c
      data spkpv/ 0.,24.,25.,27.,23.,21.,20.,19.,17.,15.5,14.,13.5,
     +13., 0.,3.6,1.7,10*0.0, .0,0.0,8.7,17.7,18.8,15.9,11.7,8.,6.,5.3,
     +4.5,3.9,3.5, .0,.0,2.8,5.8,6.2,5.1,3.8,2.7,2.1,1.8,1.5,1.3,1.1,
     +4*0.,0.0,4.6,10.2,15.1,16.9,16.5,11.,5.5,3.5, 10*0.,4.3,7.6,9.,
     +10*0.,1.7,2.6,3., 6*0.,.3,.6,1.,1.6,1.3,.8,.6, 6*0.,.7,1.2,1.8,2.
     +5,1.8,1.3,1.2, 10*0.,.6,1.4,1.7, 10*0.,1.9,4.1,5.2/
     
     
     
     
     
     
c***************   ap - p - data
      data sappel/ 0.,176.,160.,105.,75.,68.,65.,50.,50.,43.,42.,40.5,
     +35.,30.,28., 25.,22.,21.,20.,18.,17., 11*0.,.05,.15,.18,.2,.2,.3,
     +.4,.6,.7,.85, 0.,1.,.9,.46,.3,.23,.18,.16,.14,.1,.08,.05,.02,.
     +015,4*.011,3*.005,0.,55.,50.,25.,15.,15.,14.,12.,10.,7.,6.,4.,3.
     +3,2.8,2.4,2.,1.8, 1.55,1.3,.95,.75, 0.,3.3,3.,1.5,1.,.7,.4,.35,.
     +4,.25,.18,.08,.04,.03,.023,.016,.014, .01,.008,.006,.005/
     
     
      data spikpe/ 0.,215.,193.,170.,148.,113.,97.,84., 78.,68.,64.,61.
     +,46.,36.,31.3,28.5,25.7,22.6,21.4,20.7,19.9, 8*0.,0.,2.,2.5,.2,
     +9*0., 8*0.,0.,0.,.3,1.4,2.2,1.2,1.1,1.,.8,.6,.5,.4,.3, 8*0.,0.,0.
     +,.3,1.4,2.2,1.2,1.1,1.,.8,.6,.5,.4,.3, 8*0.,0.,0.,.3,1.4,2.2,1.2,
     +1.1,1.,.8,.6,.5,.4,.3, 8*0.,0.,0.,.3,1.4,2.2,1.2,1.1,1.,.8,.6,.5,
     +.4,.3, 8*0.,0.,.6,2.5,5.,5.2,5.1,5.4,5.8,2.8,2.1,1.8,1.6,1.2,
     +8*0.,5*0.,1.3,1.5,2.,2.5,2.5,2.3,1.8,1.4, 8*0.,5*0.,1.3,1.5,2.,2.
     +5,2.5,2.3,1.8,1.4, 8*0.,5*0.,1.3,1.5,2.,2.5,2.5,2.3,1.8,1.4,
     +8*0.,5*0.,1.3,1.5,2.,2.5,2.5,2.3,1.8,1.4, 8*0.,6*0.,.2,.5,1.1,1.
     +6,1.4,1.1,.9, 8*0.,6*0.,.2,.5,1.1,1.6,1.4,1.1,.9, 8*0.,6*0.,.2,.
     +5,1.1,1.6,1.4,1.1,.9, 8*0.,6*0.,.2,.5,1.1,1.6,1.4,1.1,.9, 8*0.,9*
     +0.,.3,1.6,2.6,3.6,8*0.,9*0.,.3,1.6,2.6,3.6, 8*0.,9*0.,.3,1.6,2.6,
     +3.6,8*0.,9*0.,.3,1.6,2.6,3.6/
     
     
     
     
     
     
c***************   ap - n - data
      data sapnel/ 0.,176.,160.,105.,75.,68.,65.,50.,50.,43.,42.,40.5,
     +35.,30.,28., 25.,22.,21.,20.,18.,17., 11*0.,.05,.15,.18,.2,.2,.3,
     +.4,.6,.7,.85, 0.,1.,.9,.46,.3,.23,.18,.16,.14,.1,.08,.05,.02,.
     +015,4*.011,3*.005,0.,3.3,3.,1.5,1.,.7,.4,.35,.4,.25,.18,.08,.04,.
     +03,.023,.016,.014, .01,.008,.006,.005/
     
      data spikpz/ 0.,215.,193.,170.,148.,113.,97.,84., 78.,68.,64.,61.
     +,46.,36.,31.3,28.5,25.7,22.6,21.4,20.7,19.9, 8*0.,0.,2.4,.2,10*0.
     +, 8*0.,0.,0.,1.8,2.8,3.6,2.3,1.8,1.5,1.3,1.,.7,.5,.3, 8*0.,0.,0.,
     +1.8,2.8,3.6,2.3,1.8,1.5,1.3,1.,.7,.5,.3, 8*0.,0.,0.,1.8,2.8,3.6,
     +2.3,1.8,1.5,1.3,1.,.7,.5,.3, 8*0.,0.,0.,1.8,2.8,3.6,2.3,1.8,1.5,
     +1.3,1.,.7,.5,.3, 8*0.,5*0.,5.2,8.7,11.4,14.,11.9,7.6,6.,5.,
     +8*0.,5*0.,5.2,8.7,11.4,14.,11.9,7.6,6.,5., 8*0.,10*0.,1.,4.9,8.5,
     +8*0.,10*0.,1.,4.9,8.5, 8*0.,7*0.,1.9,2.3,4.,6.5,5.2,3.4, 8*0.,7*
     +0.,1.9,2.3,4.,6.5,5.2,3.4, 8*0.,7*0.,1.9,2.3,4.,6.5,5.2,3.4/
     
     
     
     
     
     
c
c
c***************   an - p - data
c
      data sanpel/ 0.,176.,160.,105.,75.,68.,65.,50.,50.,43.,42.,40.5,
     +35.,30.,28., 25.,22.,21.,20.,18.,17., 11*0.,.05,.15,.18,.2,.2,.3,
     +.4,.6,.7,.85, 0.,1.,.9,.46,.3,.23,.18,.16,.14,.1,.08,.05,.02,.
     +015,4*.011,3*.005,0.,3.3,3.,1.5,1.,.7,.4,.35,.4,.25,.18,.08,.04,.
     +03,.023,.016,.014, .01,.008,.006,.005/
     
      data spikpf / 0.,215.,193.,170.,148.,113.,97.,84., 78.,68.,64.,
     +61.,46.,36.,31.3,28.5,25.7,22.6,21.4,20.7,19.9, 8*0.,0.,2.4,.2,
     +10*0., 8*0.,0.,0.,1.8,2.8,3.6,2.3,1.8,1.5,1.3,1.,.7,.5,.3,
     +8*0.,0.,0.,1.8,2.8,3.6,2.3,1.8,1.5,1.3,1.,.7,.5,.3, 8*0.,0.,0.,1.
     +8,2.8,3.6,2.3,1.8,1.5,1.3,1.,.7,.5,.3, 8*0.,0.,0.,1.8,2.8,3.6,2.
     +3,1.8,1.5,1.3,1.,.7,.5,.3, 8*0.,5*0.,5.2,8.7,11.4,14.,11.9,7.6,6.
     +,5., 8*0.,5*0.,5.2,8.7,11.4,14.,11.9,7.6,6.,5., 8*0.,10*0.,1.,4.
     +9,8.5, 8*0.,10*0.,1.,4.9,8.5, 8*0.,7*0.,1.9,2.3,4.,6.5,5.2,3.4,
     +8*0.,7*0.,1.9,2.3,4.,6.5,5.2,3.4, 8*0.,7*0.,1.9,2.3,4.,6.5,5.2,3.
     +4/
     
     
     
     
     
c***  ko - n - data
      data spkp15/ 0.,20.,14.,12.,11.5,10.,8.,7.,6.,5.5,5.3,5.,4.5,4.4,
     +3.8,3.,2.8, 0.,.5,1.15,2.,1.3,.8,.45,10*0., 3*0.,0.9,2.5,3.,2.5,
     +2.3,2.,1.7,1.5,1.2,.9,.6,.45,.21,.2, 3*0.,0.9,2.5,3.,2.5,2.3,2.,
     +1.7,1.5,1.2,.9,.6,.45,.21,.2, 4*0.,1.0,2.1,2.6,2.3,2.1,1.8,1.7,1.
     +4,1.2,1.05,.9,.66, .5, 7*0.,.3,1.,1.,.9,.7,.4,.30,.2,.00,0.,
     +9*0.,.1,1.,2.2,3.5,4.20,4.55,4.85,4.9, 10*0.,.2,.7,1.6,2.5,2.2,1.
     +71,1.6, 6*0.,1.4,3.8,5.,4.7,4.4,4.,3.5,2.85,2.35,2.01,1.8,
     +12*0.,.1,.8,2.05,3.31,3.5, 12*0.,.034,.20,.75,1.04,1.24 /
     
     
     
     
     
c*** ako - p - data
      data spkp16/ 0.,4.,9.5,20.,13.,9.5,6.,4.4,3.,2.4,2.,1.4,1.2,1.,.
     +9,.7,.6, 0.,4.5,6.,5.,2.5,2.,1.7,2.1,1.9,.9,.5,.3,.24,.2,.18,.1,.
     +09, 0.,1.8,2.,1.1,.9,.5,.5,.4,.4,.2,.1,.06,.05,.04,.03,.02,.02,
     +0.,1.5,2.,.9,1.1,.4,.6,.7,.65,.3,.17,.1,.08,.07,.06,.04,.03,
     +0.,.56,1.29,2.26,1.01,.64,.37,10*0., 4*0.,1.13,2.61,2.91,2.58,2.
     +35,2.02,1.91,1.57,1.35,1.29,1.01,.74, .65, 3*0.,1.00,3.03,3.36,2.
     +8,2.58,2.24,1.91,1.68,1.35,1.01,.67,.5,.24, .23, 3*0.,1.00,3.03,
     +3.36,2.8,2.58,2.24,1.91,1.68,1.35,1.01,.67,.5,.24, .23, 7*0.,.34,
     +1.12,1.12,1.01,.78,.45,.39,.22,.07,0., 6*0.,1.71,4.26,5.6,5.57,4.
     +93,4.48,3.92,3.19,2.63,2.25,2., 10*0.,.22,.8,.75,1.,1.3,1.5,1.3,
     +13*0.,.1,.3,.7,1., 13*0.,.1,.3,.7,1., 9*0.,.11,1.72,2.69,3.92,4.
     +76,5.10,5.44,5.3, 4*0.,0.00,9.2,4.7,1.9,9*0. /
     
     
     
     
     
     
     
      data nurec/9,12,5*0,10,14,3*0,1,3,5,7,6*0,2,6,16,5*0, 10,13,5*0,
     +11,12,3*0,2,4,6,8,6*0,3,15,7,5*0/
c
c*****  data without names
c     datas     datas    datas      datas     datas
c******         *********
      data ikiic/0,15,41,67,82,93,111,134,149,160,173,184,208, 225,242,
     +253,268/
      data ieiic/0,21,46,71,92,109,126,143,160,173,186,199,220, 241,
     +262,279,296/
      data iriic/0,315,965,1615,1930,2117,2423,2814,3069,3212, 3381,
     +3524,4028,4385,4742,4929,5184/
c
c     particle masses in gev
      data amz/ 3*2.2 ,0.9576,3*1.887,2.4,2.03,2*1.43,2*1.7 ,3*0./
     
     
     
c     resonance width gamma in gev
      data gaz/ 3*.2 ,.1,.2,.2,.2,.2 ,.18,.2,.2,.15,.15, 3*0./
     
     
     
c     mean life time in seconds
      data tauz/ 16*0./
     
c     charge of particles and resonances
      data ichz/ 0,1,0 ,0,0,1,-1,0, 1,-1,0,0,1 ,3*0/
     
     
     
c     baryonic charge
      data ibarz/ 2,0,0, 5*0,1,-1,-1,1,1 ,3*0/
     
     
     
c     first number of decay channels used for resonances
c     and decaying particles
      data k1z/ 308,310,313 ,317,322,365,393,421,425,434,440,446,449
     +,3*460/
     
     
c     last number of decay channels used for resonances
c     and decaying particles
      data k2z/ 309,312,316 ,321,364,392,420,424,433,439,445,448,451
     +,3*460/
     
c     weight of decay channel
      data wtz/ .17,.83,.33,.33,.34,.17,.33,.33,.17, .01,.13,.36,.27,.
     +23, .0014,.0029,.0014,.0029,.0007,.0007,.0007,.0007,.0517,.0718,
     +.0144,.0431,.0359,.0718,.0014,.0273,.0014,.0431,.0129,.0129,
     +.0259,.0517,.0359,.0014,.0144,.0144,.0129,.0014,.0259,.0359,
     +.0072,.0474,.0948,.0259,.0072,.0144,.0287,.0431,.0144,.0287,
     +.0474,.0144,.0075, .0057,.0019,.0038,.0095,.0014,.0014,.0191,.
     +0572,.1430,.0029, .0029,.0477,.0477,.0477,.0477,.0477,.0019,.
     +0191,.0686,.0172, .0095,.1888,.0172,.0191,.0381,.0571,.0571,.
     +0190, .0057,.0019,.0038,.0095,.0014,.0014,.0191,.0572,.1430,.
     +0029, .0029,.0477,.0477,.0477,.0477,.0477,.0019,.0191,.0686,.
     +0172, .0095,.1888,.0172,.0191,.0381,.0571,.0571,.0190, 4*.25,.2,.
     +2,.12,.1,.07,.07,.14,.05,.05,.4,.2,.125,.075,.075,.125, .4,.075,.
     +125,.2,.125,.075,.3,.05,.65,.3,.05,.65, 9*1./
     
     
     
c     particle numbers in decay channel
      data nzk1/ 8,1,2,9,1,2,9,2,9, 7,13,31,15,24, 23,13,23,13,23,23,
     +14,13,23,31,98,33,33,32, 23,14,13,35,23,23,14,13,33,23,98,31,23,
     +14, 13,35,33,33,32,23,35,33,32,98,35,35,35,35,35, 13,13,13,13,23,
     +13,98,32,33,23,13,23,13,14, 13,32,13,98,23,13,32,32,13,33,32,98,
     +35,35, 14,14,14,14,23,14,98,34,34,23,14,23,14,14, 13,34,14,98,23,
     +14,34,34,14,33,32,98,35,35,104,61,105,62, 1,17,21,17,22,21,21,22,
     +21,2,67,68,69,2,9,9,68,69,70,2,9, 24,24,15,25,25,16, 9*0/
     
     
     
     
     
      data nzk2/ 8,8,1,8,9,8,8,1,1, 7,14,13,16,25, 23,14,23,14,31,33,
     +32,34,35,31,23,31,33, 34,31,32,34,31,33,32,33,33,35,31,33,31,33,
     +32,34,35,31,33,34,35,31,33,33,33,33,32,35, 35,35, 23,23,13,31,32,
     +33,13,31,32,31,31,32,33,32, 32,35,31,32,32,33,31,33,35,33,32,32,
     +32,35, 23,23,14,31,34,33,14,31,33,31,31,34,32,33, 34,35,31,34,34,
     +33,31,33,35,33,34,34,33,35,1,2,8,9, 25,13,35,32,32,33,31,13,23,
     +31,13,23,14,79,80,31,13,23,14,78,79, 8,1,8,1,8,1, 9*0/
     
     
     
     
     
     
     
      data nzk3/ 23,14,13,13,23,13,23,23,14, 0,7,14,0,0, 0, 0,23,23, 0,
     +0, 0, 0, 0, 0, 0, 0, 0, 0,33,31,31, 0,33,34,32,34, 0,35, 0,31,
     +35,35,35, 0,31,31,31,35,31,33,34,31,33,34,31,33,35, 0,23,14, 0,
     +0, 0, 0, 0, 0,32,33,33,33,32, 34, 0,35, 0,35,35,31,31,35,32,34,
     +31,33,32, 0,23,13, 0, 0, 0, 0, 0, 0,34,33,33,34,33, 34, 0,35, 0,
     +35,35,31,31,35,34,34,31,34,34,25*0,23,14,14,23,13,13, 9*0/
     
     
     
     
     
c
      runtes=100.
      eftes =100.
      do 10 i=1,17
         irii(i)=iriic(i)
         ikii(i)=ikiic(i)
         ieii(i)=ieiic(i)
   10 continue
      do 20 i=1,296
         siin(i)=0.
         plabf(i)=plabfc(i)
         umo(i)=umoc(i)
   20 continue
*
      iof=0
      num=315
      do 30 i=1,num
         wkh(iof+i)=spikp1(i)
   30 continue
      iof=iof+num
*
      num=278
      do 40 i=1,num
         wkh(iof+i)=spikpu(i)
   40 continue
      iof=iof+num
*
      num=372
      do 50 i=1,num
         wkh(iof+i)=spikpv(i)
   50 continue
      iof=iof+num
*
      num=278
      do 60 i=1,num
         wkh(iof+i)=spikpw(i)
   60 continue
      iof=iof+num
*
      num=372
      do 70 i=1,num
         wkh(iof+i)=spikpx(i)
   70 continue
      iof=iof+num
*
      num=315
      do 80 i=1,num
         wkh(iof+i)=spikp4(i)
   80 continue
      iof=iof+num
*
      num=187
      do 90 i=1,num
         wkh(iof+i)=spikp5(i)
   90 continue
      iof=iof+num
*
      num=306
      do 100 i=1,num
         wkh(iof+i)=spikp6(i)
  100 continue
      iof=iof+num
*
      num=102
      do 110 i=1,num
         wkh(iof+i)=skmpel(i)
  110 continue
      iof=iof+num
*
      num=289
      do 120 i=1,num
         wkh(iof+i)=spikp7(i)
  120 continue
      iof=iof+num
*
      num=68
      do 130 i=1,num
         wkh(iof+i)=skmnel(i)
  130 continue
      iof=iof+num
*
      num=187
      do 140 i=1,num
         wkh(iof+i)=spikp8(i)
  140 continue
      iof=iof+num
*
      num=143
      do 150 i=1,num
         wkh(iof+i)=spikp9(i)
  150 continue
      iof=iof+num
*
      num=169
      do 160 i=1,num
         wkh(iof+i)=spikp0(i)
  160 continue
      iof=iof+num
*
      num=143
      do 170 i=1,num
         wkh(iof+i)=spkpv(i)
  170 continue
      iof=iof+num
*
      num=105
      do 180 i=1,num
         wkh(iof+i)=sappel(i)
  180 continue
      iof=iof+num
*
      num=399
      do 190 i=1,num
         wkh(iof+i)=spikpe(i)
  190 continue
      iof=iof+num
*
      num=84
      do 200 i=1,num
         wkh(iof+i)=sapnel(i)
  200 continue
      iof=iof+num
*
      num=273
      do 210 i=1,num
         wkh(iof+i)=spikpz(i)
  210 continue
      iof=iof+num
*
      num=84
      do 220 i=1,num
         wkh(iof+i)=sanpel(i)
  220 continue
      iof=iof+num
*
      num=273
      do 230 i=1,num
         wkh(iof+i)=spikpf(i)
  230 continue
      iof=iof+num
*
      num=187
      do 240 i=1,num
         wkh(iof+i)=spkp15(i)
  240 continue
      iof=iof+num
*
      num=255
      do 250 i=1,num
         wkh(iof+i)=spkp16(i)
  250 continue
      iof=iof+num
*
      j=0
      do 260 i=1,268
         nrk(1,i)=nktemp(j+1)
         nrk(2,i)=nktemp(j+2)
         j=j+2
  260 continue
      do 270 i=1,30
         nure(i,1)=nurec(i)
         nure(i,2)=nurec(i+30)
  270 continue
*
      do 280 i=1,16
         l=i+94
         am (l)=amz(i)
         ga( l)=gaz(i)
         tau( l)=tauz(i)
         ich( l)=ichz(i)
         ibar( l)=ibarz (i)
         k1( l)=k1z(i)
         k2( l)=k2z(i)
  280 continue
      do 290 i=1,153
         l=i+307
         wt( l)=wt z(i)
         nzk( l,3)=nzk3(i)
         nzk( l,2)=nzk2(i)
         nzk( l,1)=nzk1(i)
  290 continue
      end