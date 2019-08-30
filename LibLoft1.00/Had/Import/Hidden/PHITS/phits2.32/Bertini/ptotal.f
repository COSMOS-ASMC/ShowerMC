************************************************************************
*                                                                      *
      subroutine ptotal(e,ind9,ind10,ind11,scatt,absor)
c ----------------------------------------------------------------------
c     given the kinetic energy in the c.m.routine will return the total
c     i-i or i-j pi-nucleon cross section for that energy and the
c     absorption cross section as well as the pi0-nucleon cross section
c     ind9=1 means an i-j collision, ind10=1 means an i-i collision, and
c     ind11=1 means an i-0 collisision.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c

      dimension scaii(67),scaij(67),scao(67)
      data scaii/0.0,2.,4.,6.,12.,19.,29.,42.,58.,77.,99.,129.,159.,179.
     1,191.,195.,185.,162.,143.,124.,109.,95.,83.,73.,64.,56.,50.,44.,
     240.,36.,33.,31.,29.,27.,25.,24.,22.,21.,20.,19.,18.,
     x    17.,17.,16.,16.,16.,16.,16.,17.,17.,18.,18.
     1,19.,20.,21.,22.,23.,24.,24.,25.,25.,26.,26.,26.,27.,27.,27./
      data scaij/0.0,2.,4.,6.,8.,11.,14.,18.,22.,29.,37.,47.,55.,62.,66.
     1,67.,65.,59.,52.,44.,39.,36.,33.,31.,29.,28.,27.,27.,27.,27.,27.,
     227.,27.,28.,29.,30.,31.,33.,34.,36.,38.,
     x    42.,45.,47.,48.,47.,45.,42.,40.,38.,38.,38.
     1,39.,40.,42.,45.,49.,52.,56.,58.,59.,58.,57.,54.,50.,47.,45./
      data scao/0.0,2.,4.,6.,10.,15.,22.,30.,40.,53.,68.,88.,107.,121.,
     1129.,131.,125.,111.,98.,84.,74.,66.,58.,52.,47.,42.,39.,36.,34.,
     232.,30.,29.,28.,28.,27.,27.,27.,27.,27.,28.,28.,
     x     30.,31.,32.,32.,32.,31.,29.,29.,28.,28.,28.,
     129.,30.,32.,34.,36.,38.,40.,42.,42.,42.,42.,40.,39.,37.,36./
c above tables are based on cern/hera 69-1 report
c above tables are good for c.m.k.e.e=660mev or lab.pi k.e.=990mev
c ----------------------------------------------------------------------
c
      if(e.ge.660.)go to 1
      absor=0.0
c the following is based on fact that above data are at equally
c spaced energy intervals of 10 mev
      rie=.1*e
      ie=rie
      ie=ie+1
      de=rie- float(ie-1)
      ae=1.-de
      if(ind9.eq.1)go to 10
      if(ind10.eq.1)go to 13
c zero ind(9),ind(10),and ind(11) before returning
      ind9=0
      ind10=0
      ind11=0
      scatt=ae*scao(ie)+de*scao(ie+1)
      return
10    continue
c zero ind(9),ind(10),and ind(11) before returning
      ind9=0
      ind10=0
      ind11=0
      scatt=ae*scaii(ie)+de*scaii(ie+1)
      return
13    continue
c zero ind(9),ind(10),and ind(11) before returning
      ind9=0
      ind10=0
      ind11=0
      scatt=ae*scaij(ie)+de*scaij(ie+1)
      return
c%%    1 print 100,e
c%%  100 format('1energy=',f10.4,' too large for ptotal')
c the following gives an automatic dump with correct control cards
c%%   x=0.0
c%%   y=1./x
c%%   z=y*y
    1 write(6,100) e
  100 format(/' *** error message from s.ptotal ***'
     &/' given energy was too large for ptotal.'
     &/' e =',1pe13.6)
      call parastop( 823 )
      end subroutine
************************************************************************
*                                                                      *
      subroutine gamma
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/switch/reat,vpion,isonsw,idistr,icb,nzman
      common/hiro/st(8000)
      common/palcom/npl(2),monin
      common/pinin/tcm,i9,i10,i11
      common/kstep/nst(30),ncol
      common/vegas/csa(14),icsb(16),ws(30),wsp(30),stp(200),den(9),
     1  denn(9),denp(9),ef(2,9),pf(2,9),poe1,cutofa(9),ind(40),
     2  amass(5),gws(15),cm(25),ela(20),cm1(21),cm2(14),wws(30),
     3  wwsp(30),aip(3),ain(3),aiscat,aidec,aicap,aitot,trick,fscat
      data qvpn,  qvnp  / 18.1203, 12.5867/
      data qalpn, qalnp /  5.8023,  1.8960/
      data q63pn, q63np /  4.1500,  0.0000/
      data q65pn, q65np /  2.1333,  1.3750/
      data qfepn, qfenp /  5.4447,  2.9662/
      data qpbpn, qpbnp /  3.6781,  4.2321/
c ----------------------------------------------------------------------
c
  108 format(24h interactions exceeds 12)
      notrys=0
      ind21=ind(21)
      if(ind(22)-1) 31,8,31
   31 go to (8,8,29,8,8,29),ind21
   29 if(ind(20)-1) 8,21,21
c
    8 ind(1)=1
      call look(gws(1),istep)
c
      if(ws(3)-10.0)40,50,50
c
   40 call cma
c
      if(ind(28).eq.1)go to 400
c
      if(ws(3).gt.2.)go to 500
c
      psel6=(cm(3)-ws(4)-wsp(4))*0.510976
c
    2 call anglen(ws(3),wsp(3),psel6,ans)
c
   79 call elas(ans)
c
   80 paul1=(cm(4)*cm(11)+cm(2)*cm(14))/cm(3)
c
  304 if(paul1-wsp(4)-wsp(24)) 3,4,4
c
    3 if(ind(29).eq.1) go to 3040
      if(ind(28).eq.1) go to 3040
      if(ind(27).eq.1) go to 3040
      if(ind(13).eq.1.and.fscat.ne.1) go to 3040
      iden1=ws(3)
      iden2=wsp(3)
      go to 3041
 3040 iden1=wws(3)
      iden2=wwsp(3)
 3041 ik=2
      if(iden1.le.2.and.iden2.le.2) ik=1
      npl(ik)=npl(ik)+1
c
 3043 if(ind(22)-1) 22,41,22
c
   41 wsp(3)=wwsp(3)
      wsp(4)=wwsp(4)
      wsp(5)=wwsp(5)
      ws(3)=wws(3)
      ws(4)=wws(4)
      ws(5)=wws(5)
      wsp(9)=wwsp(9)
      wsp(24)=wwsp(24)
      ws(24)=wws(24)
      ind(13)=0
      ind(14)=0
      ind(15)=0
      ind(26)=0
      ind(18)=0
      notrys=notrys+1
      if(notrys.lt.10)go to 16
c
      aidec=0.0
      aicap=1.
      random=unirn(dummy)
      wsp(15)=random
      wsp(15)=wsp(15)+wsp(15)-1.
      go to 16
c
   22 ind(20)=1
      if(ind(28).eq.1)call resto
      if(ind(27).eq.1)call resto
      if(ind(29).eq.1)call resto
      if(ind(15)-1)124,123,124
  123 call resto
  124 if(ind(14)-1)126,125,126
  125 call resto
  126 if(ind(13)-1)128,127,128
  127 if(fscat-1.0) 129,128,129
  129 call resto
  128 if(icsb(3)) 5,6,5
    6 icsb(3)=icsb(3)+1
    5 return
    4 paul1=(cm(2)*cm(10)+cm(4)*cm(7))/cm(3)
      if(ws(3)-2.0) 130,130,7
  130 if(paul1-ws(4)-ws(24)) 3,7,7
c
    7 ws(15)=ws(19)-ws(16)
      ncol = ncol +1
 2003 ind(22)=0
      ind(24)=0
      ind(25)=0
   81 call cmb
      call cmc(cm(7),cm(8),cm(9),cm(10))
      ws(5)=cm2(3)
      ws(6)=cm2(4)/cm2(2)
      ws(7)=cm2(5)/cm2(2)
      ws(8)=cm2(6)/cm2(2)
      call cmc(cm(11),cm(12),cm(13),cm(14))
      wsp(5)=cm2(3)
      wsp(6)=cm2(4)/cm2(2)
      wsp(7)=cm2(5)/cm2(2)
      wsp(8)=cm2(6)/cm2(2)
      ind(3)=0
      if(ind(13).eq.1)ws(30)=ws(30)+1.
      if(ind(28).eq.1)go to 28
      if(ind(27).eq.1)go to 28
      if(ind(29).eq.1)go to 28
      if(ind(13)-1) 86,92,86
   86 if(ind(15)-1)88,87,88
   87 call wwsst
      go to 38
c
   88 if(ind(14)-1) 38,28,38
   92 if(ind(26)-1) 82,28,82
   28 icsb(15)=icsb(15)+1
      icsb8=icsb(15)*icsb(5)
      st(icsb8+1)=-ws(1)
      st(icsb8+2)=5.0
      st(icsb8+3)=wwsp(3)
      st(icsb8+4)=wwsp(4)
      st(icsb8+5)=wwsp(5)
      st(icsb8+6)=0.0
      st(icsb8+7)=0.0
      st(icsb8+8)=0.0
      st(icsb8+9)=wwsp(9)
      st(icsb8+10)=ws(10)
      st(icsb8+18)=ws(10)
      st(icsb8+11)=ws(11)
      st(icsb8+19)=ws(11)
      st(icsb8+12)=ws(12)
      st(icsb8+20)=ws(12)
      st(icsb8+17)=ws(15)
      st(icsb8+13)=ws(16)
      st(icsb8+14)=ws(14)
      st(icsb8+15)=ws(25)
      st(icsb8+16)=wwsp(24)
      st(icsb8+21)=0.0
      call wwsst
      go to 38
c
   82 wsp(9)=wwsp(9)
   38 if(icsb(3)) 9,10,9
   10 icsb(3)=icsb(3)+2
      go to 11
    9 icsb(3)=icsb(3)+1
   11 icsb(1)=icsb(1)+1
      icsb(2)=icsb(2)+1
      wsp(1)=ws(1)
      wsp(2)=5.0
      i=wsp(3)
      if(ind(28).eq.1)go to 312
      if(ind(27).eq.1)go to 312
      if(ind(29).eq.1)go to 312
      if(ind(13)-1) 311,310,311
  310 if(fscat-1.0) 312,104,312
  312 j=wwsp(3)
      ws(j+21)=ws(j+21)+1.0
      monin=monin+1
      ws(15)=ws(15)+j*100.
      if(ind(27).eq.1)go to 14
      go to 103
  311 if(ind(14)-1) 102,103,102
  102 if(ind(15)-1)104,103,104
  104 ws(i+21)=ws(i+21)+1.0
      monin=monin+1
      ws(15)=ws(15)+i*100.
c
 103  if(cutofa(i+2)+wsp(24)+wsp(4)-wsp(5)) 14,15,15
c
   15 wsp(1)=-wsp(1)
      icsb(15)=icsb(15)+1
      icsb(1)=icsb(1)-1
      icsb(2)=icsb(2)-1
      icsb(3)=icsb(3)-1
      wsp(13)=ws(16)
      icsb8=icsb(15)*icsb(5)
   17 st(icsb8+1)=wsp(1)
      st(icsb8+2)=wsp(2)
      st(icsb8+3)=wsp(3)
      st(icsb8+4)=wsp(4)
      st(icsb8+5)=wsp(5)
      st(icsb8+6)=wsp(6)
      st(icsb8+7)=wsp(7)
      st(icsb8+8)=wsp(8)
      st(icsb8+9)=wsp(9)
      st(icsb8+18)=ws(10)
      st(icsb8+19)=ws(11)
      st(icsb8+20)=ws(12)
      st(icsb8+21)=wsp(5)
      st(icsb8+10)=ws(10)
      st(icsb8+11)=ws(11)
      st(icsb8+12)=ws(12)
      st(icsb8+13)=wsp(13)
      st(icsb8+14)=ws(14)
      st(icsb8+17)=ws(15)
      st(icsb8+15)=ws(25)
      st(icsb8+16)=wsp(24)
131   continue
c
      if(ws(3).gt.10.)go to 20
  105 if(ws(3)-3.0)23,20,20
   23 if(ws(3)-1.)18,19,18
   19 if(cutofa(3)+ws(24)+ws(4)-ws(5)) 20,21,21
   18 if(cutofa(4)+ws(24)+ws(4)-ws(5)) 20,21,21
c
   21 iws2=ws(2)
c
  116 icsb13=icsb(13)+iws2+16
c
      if(ind(20)-1) 46,39,46
   46 if(ind(14)-1) 47,48,47
c
   48 icsb13=icsb13-5
   47 st(icsb13+1)=ws(15)
      st(icsb13+2)=ws(10)
      st(icsb13+3)=ws(11)
      st(icsb13+4)=ws(12)
      st(icsb13+5)=ws(5)
      if(ind(14)-1) 62,59,62
   62 ws(2)=ws(2)+5.0
   59 icsb(3)=icsb(3)-1
   39 ws(1)=-ws(1)
      icsb(1)=icsb(1)-1
      icsb13=icsb(13)
      st(icsb13+1)=ws(1)
      st(icsb13+2)=ws(2)
      st(icsb13+3)=ws(3)
      st(icsb13+4)=ws(4)
      st(icsb13+5)=ws(5)
      st(icsb13+6)=ws(6)
      st(icsb13+7)=ws(7)
      st(icsb13+8)=ws(8)
      st(icsb13+9)=ws(9)
      st(icsb13+10)=ws(10)
      st(icsb13+11)=ws(11)
      st(icsb13+12)=ws(12)
      st(icsb13+13)=ws(13)
      st(icsb13+14)=ws(14)
      st(icsb13+15)=ws(25)
      st(icsb13+16)=ws(24)
      ind(20)=0
c
      if(ws(2)-60.0)107,106,106
  106 print 108
      ind(19)=1
      go to 24
  107 if(icsb(3))24,26,24
   26 icsb(11)=0
      icsb(16)=0
      icsb(10)=icsb(9)*icsb(11)
      ind(1)=0
      go to 24
c
   20 if(ind(15)-1)137,135,137
  137 if(ind(14)-1)136,135,136
  136 if(ind(26)-1) 49,135,49
49    continue
      if(ind(29).eq.1)go to 135
      if(ind(28).eq.1)go to 135
      iws2=ws(2)
  122 icsb13=icsb(13)+iws2+14
      st(icsb13+4)=ws(10)
      st(icsb13+5)=ws(11)
      st(icsb13+6)=ws(12)
      st(icsb13+3)=ws(15)
      st(icsb13+15)=ws(25)
      st(icsb13+7)=ws(5)
      ws(2)=ws(2)+5.0
      if(ws(3).lt.10.)go to 119
      icsb13=icsb(13)+77
      st(icsb13)=10.
      ws(26)=10.
  119 ind(20)=1
      if(ws(2)-60.0)24,106,106
   24 return
  135 icsb13=icsb(13)
      st(icsb13+2)=ws(2)
      st(icsb13+3)=ws(3)
      st(icsb13+4)=ws(4)
      st(icsb13+5)=ws(5)
      st(icsb13+6)=ws(6)
      st(icsb13+7)=ws(7)
      st(icsb13+8)=ws(8)
      st(icsb13+9)=ws(9)
      st(icsb13+10)=ws(10)
      st(icsb13+11)=ws(11)
      st(icsb13+12)=ws(12)
      st(icsb13+13)=ws(13)
      st(icsb13+14)=ws(14)
      st(icsb13+15)=ws(25)
      st(icsb13+16)=ws(24)
      st(icsb13+17)=ws(15)
      st(icsb13+18)=ws(10)
      st(icsb13+19)=ws(11)
      st(icsb13+20)=ws(12)
      st(icsb13+21)=ws(5)
      if(ws(3).lt.10.)go to 119
      st(icsb13+77)=10.
      ws(26)=10.
      if(ind(26).eq.1.)go to 119
      call updeca
      ws(29)=10.
      st(icsb13+77)=ws(26)
      st(icsb13+78)=ws(27)
      st(icsb13+79)=ws(28)
      st(icsb13+80)=ws(29)
      go to 119
c
   14 icsb(16)=icsb(16)+1
      icsb8=icsb(9)*(icsb(16)-1)
      stp(icsb8+1)=wsp(1)
      stp(icsb8+2)=wsp(2)
      stp(icsb8+3)=wsp(3)
      stp(icsb8+4)=wsp(4)
      stp(icsb8+5)=wsp(5)
      stp(icsb8+6)=wsp(6)
      stp(icsb8+7)=wsp(7)
      stp(icsb8+8)=wsp(8)
      stp(icsb8+9)=wsp(9)
      stp(icsb8+10)=ws(10)
      stp(icsb8+18)=ws(10)
      stp(icsb8+11)=ws(11)
      stp(icsb8+19)=ws(11)
      stp(icsb8+12)=ws(12)
      stp(icsb8+20)=ws(12)
      stp(icsb8+21)=wsp(5)
      stp(icsb8+13)=wsp(13)
      stp(icsb8+14)=ws(14)
      stp(icsb8+17)=ws(15)
      stp(icsb8+15)=ws(25)
      stp(icsb8+16)=wsp(24)
      if(ind(27).eq.1)go to 135
      go to 131
c
50    continue
      if(ind(22).eq.1)go to 35
c
      go to (13,37,37,13,37,37),ind21
c
37    if(ind(22)-1)13,35,13
   35 if(ind(12)-1)12,32,12
   32 if(ws(3)-13.0) 34,33,34
   33 if(wsp(3)-1.0) 12,43,12
   34 if(wsp(3)-2.0) 12,44,12
   43 wsp(3)=2.0
      wsp(4)=1838.67
      wsp(9)=pf(2,istep)
      wsp(24)=ef(2,istep)
      go to 45
   44 wsp(3)=1.0
      wsp(4)=1836.14
      wsp(9)=pf(1,istep)
      wsp(24)=ef(1,istep)
   45 random=unirn(dummy)
      gdc10=random*1.d0
      wsp(9)=wsp(9)*cubert(gdc10)
      wsp(5)=dsqrt(wsp(9)*wsp(9)+wsp(4)*wsp(4))
      wsp(19)=wsp(4)*wsp(4)
      ind(24)=1
      call ainter(gdc9,gdc10)
c
   12 aiscat=0.0
      aitot=aidec+aicap
      wwsp(3)=wsp(3)
      wwsp(4)=wsp(4)
      wwsp(5)=wsp(5)
      wws(3)=ws(3)
      wws(4)=ws(4)
      wws(5)=ws(5)
      wwsp(9)=wsp(9)
      wwsp(24)=wsp(24)
      wws(24)=ws(24)
c
   13 if(aitot.eq.0.0) go to 16
      aiscat=aiscat/aitot
      aidec=aidec/aitot
      aicap=aicap/aitot
   16 gws3=dsqrt(gws(5)*gws(5)-ws(4)*ws(4))
      wsp3=dsqrt(wsp(5)*wsp(5)-wsp(4)*wsp(4))
      e=((gws(5)*wsp(5)-wsp(15)*gws3*wsp3)/wsp(4)-ws(4))*0.510976
      u1=ws(4)*ws(4)+wsp(4)*wsp(4)
      u2=2.0*(gws(5)*wsp(5)-wsp(15)*gws3*wsp3)
      u=0.000510976*dsqrt(u1+u2)
      random=unirn(dummy)
      gf=random
      if(gf-aiscat)51,52,52
   52 if(gf-aiscat-aidec)54,53,53
   54 call saveis
      call decay(ans)
      if(ind(18)-1) 79,3,79
c
   51 ind(13)=1
      if(u-2.1)58,57,57
   58 u=2.1
   57 v=ws(4)*0.000510976
      if(v-1.11)60,61,61
   60 v=1.11
   61 call cma
      if(fscat-1.0) 56,55,56
c
   55 call aiangl(e,u,v,ans,ind(13),ind(14))
      go to 79
c
   56 call saveis
      wwsp(3)=wsp(3)
      wwsp(4)=wsp(4)
      wwsp(24)=wsp(24)
      wwsp(5)=wsp(5)
      wwsp(9)=wsp(9)
      wsp(9)=0.0
      if(fscat-0.33333)64,63,64
   63 ind(26)=1
      if(ws(3)-13.0)91,89,91
   89 ws(3)=14.0
      ws(24)=ef(1,istep)
  109 wsp(3)=1.0
      wsp(4)=amass(1)
      wsp(24)=ef(1,istep)
      delta=ef(1,istep)-ef(2,istep)+icb*csa(11)
      go to 93
   91 ws(3)=24.0
      ws(24)=ef(2,istep)
 101  wsp(3)=2.0
      wsp(4)=amass(2)
      wsp(24)=ef(2,istep)
      delta=ef(2,istep)-ef(1,istep)-icb*csa(11)
      go to 93
   64 random=unirn(dummy)
      gf=random*1.d0
      if(fscat-0.77778) 66,65,66
   65 if(gf*0.77778-0.44444) 97,94,94
   97 delta=0.0
      go to 93
c
   94 ind(26)=1
      if(ws(3)-14.0) 95,96,95
   96 ws(3)=13.0
      ws(24)=ef(1,istep)
      go to 101
   95 ws(3)=25.0
      ws(24)=ef(2,istep)
      go to 109
   66 if(gf*0.55555-0.44444) 98,97,97
   98 delta=0.0
      if(ws(3)-14.0) 100,99,100
   99 ws(3)=24.0
      wsp(3)=1.0
      wsp(4)=amass(1)
      wsp(24)=ef(1,istep)
      ws(24)=ef(2,istep)
      go to 93
  100 ws(3)=14.0
      wsp(3)=2.0
      wsp(4)=amass(2)
      wsp(24)=ef(2,istep)
      ws(24)=ef(1,istep)
c
93    continue
      wsp(19)=wsp(4)*wsp(4)
      call pih(delta)
      if(ind(18)-1) 55,3,55
c
   53 call saveis
      ind(14)=1
      if(u-2.1)67,68,68
   67 u=2.1
   68 wwsp(3)=wsp(3)
      wwsp(4)=wsp(4)
      wwsp(5)=wsp(5)
      wwsp(24)=wsp(24)
      wwsp(9)=wsp(9)
      call cma
      if(ind(12)-1)70,69,70
   69 if(ws(3)-13.0)71,71,72
   72 ws(3)=2.0
      wsp(3)=2.0
      wsp(24)=ef(2,istep)
      delta=ef(2,istep)-ef(1,istep)-icb*csa(11)
      go to 73
   71 ws(3)=1.0
      wsp(3)=1.0
      wsp(24)=ef(1,istep)
      delta=ef(1,istep)-ef(2,istep)+icb*csa(11)
      go to 73
   70 delta=0.0
      if(ws(3)-14.0)76,75,76
   76 ws(3)=2.0
      go to 74
   75 ws(3)=1.0
      go to 74
   73 i=wsp(3)
      wsp(4)=amass(i)
      wsp(19)=wsp(4)*wsp(4)
   74 i=ws(3)
      ws(4)=amass(i)
      if(ind21.gt.3)delta=delta-vpion
      call pih(delta)
      if(ind(18)-1) 320,3,320
  320 v=wws(4)*0.000510976
      if(v-1.11)77,78,78
   77 v=1.11
c
   78 call aiangl(e,u,v,ans,ind(13),ind(14))
      wsp(9)=0.0
      ws(9)=0.0
      go to 79
400   continue
      call saveis
      wwsp(3)=wsp(3)
      wwsp(4)=wsp(4)
      wwsp(5)=wsp(5)
      wwsp(24)=wsp(24)
      wwsp(9)=wsp(9)
      ecm=.510976*cm(3)
      call onepin(ecm,wwsp(3),wws(3),ws(3),wsp(3),ws(4),wsp(4))
      delta=0.0
      if(ws(3).eq.13.)delta=delta+ef(2,istep)-ef(1,istep)-icb*csa(11)
      if(ws(3).eq.25.)delta=delta-ef(2,istep)+ef(1,istep)+icb*csa(11)
      if(ind21.gt.3)delta=delta+vpion
      i=wsp(3)
      wsp(4)=amass(i)
      wsp(19)=wsp(4)*wsp(4)
      ws(4)=ws(4)/.510976
      call pih(delta)
      if(ind(18).eq.1)go to 22
      call angli1(ans)
      wsp(9)=0.0
      ws(9)=0.0
      if(wsp(3).eq.1.)wsp(24)=ef(1,istep)
      if(wsp(3).eq.2.)wsp(24)=ef(2,istep)
      if(ws(3).le.14.)ws(24)=ef(1,istep)
      if(ws(3).ge.24.)ws(24)=ef(2,istep)
      go to 79
c
500   continue
      call branch(tcm,i9,i10,i11,fin,finex)
      random=unirn(dummy)
      fint=random
      if(fint.le.finex)go to 510
      call pangle(tcm,i9,i10,i11,1.d0,1.d0,ans)
      go to 79
510   continue
      if(fint.le.fin)go to 550
      call saveis
      wwsp(3)=wsp(3)
      wwsp(4)=wsp(4)
      wwsp(5)=wsp(5)
      wwsp(24)=wsp(24)
      wwsp(9)=wsp(9)
      ind(29)=1
      delta=ef(2,istep)-ef(1,istep)-icb*csa(11)
      if(ws(3)-4.)515,520,530
515   continue
c
      ws(3)=4.
      wsp(3)=1.
      delta=-delta
      go to 540
520   continue
      if(wsp(3).eq.1.)go to 525
c
      delta=-delta
      ws(3)=5.
      wsp(3)=1.
      go to 540
525   continue
c
      ws(3)=3.
      wsp(3)=2.
      go to 540
530   continue
c
      ws(3)=4.
      wsp(3)=2.
540   continue
      i=ws(3)
      j=wsp(3)
      ws(4)=amass(i)
      wsp(4)=amass(j)
      wsp(19)=wsp(4)*wsp(4)
      call pih(delta)
      if(ind(18).eq.1)go to 22
      call pangle(tcm,i9,i10,i11,1.d0,2.d0,ans)
      if(wsp(3).eq.1.)wsp(24)=ef(1,istep)
      if(wsp(3).eq.2.)wsp(24)=ef(2,istep)
      wsp(9)=0.0
      ws(9)=0.0
      go to 79
550   continue
c
      ind(27)=1
      call saveis
      wwsp(3)=wsp(3)
      wwsp(4)=wsp(4)
      wwsp(5)=wsp(5)
      wwsp(24)=wsp(24)
      wwsp(9)=wsp(9)
      delta=ef(2,istep)-ef(1,istep)-icb*csa(11)
      ecm=.510976*cm(3)
      call pin1pi(ecm,wws(3),wwsp(3),ws(3),wsp(3),ws(4),wsp(4),dfact)
      delta=dfact*delta
      if(ind21.gt.3)delta=delta+vpion
      i=wsp(3)
      wsp(4)=amass(i)
      wsp(19)=wsp(4)*wsp(4)
      ws(4)=ws(4)/.510976
      call pih(delta)
      if(ind(18).eq.1)go to 22
      call angli1(ans)
      wsp(9)=0.0
      ws(9)=0.0
      if(ws(3).le.14.)ws(24)=ef(1,istep)
      if(ws(3).ge.24.)ws(24)=ef(2,istep)
      wsp(24)=0.0
      call elas(ans)
      go to 7
      end subroutine
      function cij(e)
c ----------------------------------------------------------------------
c     program calculates total cross sections for i-j(p-n)collisions
c     based on vegas polynomial fits below 400 mev and data of bugg et
c     al above.  cross sections are now up to 2.2 gev.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c ----------------------------------------------------------------------
c
      if(e.gt.40.)go to 10
      cij=6.9466+(9069.2-5057.4/e)/e
      go to 1000
10    continue
      if(e.gt.400.)go to 20
      cij=27.147+(1802.0+239380./e)/e
      go to 1000
20    continue
      if(e.gt.656.)go to 30
      cij=24.506+.021484*e
      go to 1000
30    continue
      if(e.gt.754.)go to 40
      cij=33.245+.0081633*e
      go to 1000
40    continue
      if(e.gt.1074.)go to 50
      cij=36.573+.00375*e
      go to 1000
50    continue
      if(e.gt.1343.)go to 60
      cij=35.410+.0048327*e
      go to 1000
60    continue
      if(e.gt.1685.)go to 70
      cij=38.758+.0023392*e
      go to 1000
70    continue
      cij=41.389+.00077821*e
1000  return
      end function

************************************************************************
*                                                                      *
      function cii(e)
c ----------------------------------------------------------------------
c     program calculates total cross-section for i-i (p-p) collisions
c     based on vegas polynomial fits below 450 mev and data of bugg et
c     al above.  cross sections are now up to 2.2 gev.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c ----------------------------------------------------------------------
c
      if(e.gt.40.)go to 10
      cii=5.3107+(3088.5-1174.8/e)/e
      go to 1000
10    continue
      if(e.gt.310.)go to 20
      cii=22.429-(11.148-93074./e)/e
      go to 1000
20    continue
      if(e.gt.450.)go to 30
      cii=3.5475+0.05331*e+887.37/e
      go to 1000
30    continue
      if(e.gt.656.)go to 40
      cii=-.42718+.066505*e
      go to 1000
40    continue
      if(e.gt.754.)go to 50
      cii=21.110+.033673*e
      go to 1000
50    continue
      if(e.gt.923.)go to 60
      cii=42.038+.0059172*e
      go to 1000
60    continue
      if(e.gt.1217.)go to 70
      cii=47.814-.00034014*e
      go to 1000
70    continue
      cii=50.994-.0029532*e
1000  return
      end function
