*********************************************************************\
***                                                                
      subroutine bert(finput,nopart,nomp,kind,eray,aray,bray,gray)
*
*                                                                      *
*       main routine of bertini model                                  *
*       last modified by K.Niita on 01/06/2001                         *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       finput  : input data                                           *
*        finput(1)  : target mass                                      *
*        finput(2)  : target charge                                    *
*        finput(3)  : projectile energy (MeV)                          *
*        finput(4)  : 0.0 fixed                                        *
*        finput(5)  : 1.0 fixed                                        *
*        finput(6)  : andit ,option of angular distribution of delta   *
*        finput(7)  : projectile id, ityp - 1                          *
*                                                                      *
*       nomp    : dimension size                                       *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       nopart  : number of out going particles                        *
*       kind(i) : particle type of i-th particle                       *
*       eray(i) : energy of i-th particle (MeV)                        *
*       aray(i),bray(i),gray(i) : unit momentum vector of i-th particle*
*                                                                      *
*      npcle,nhole : particle and hole number for pre-equ.             *
*      efermi      : fermi energy for pre-equ.                         *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'bert.inc'

      common/rtcom/tapcrs(6600)

*-----------------------------------------------------------------------

      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /preeq/  npcle, nhole, efermi, atar, ztar

*-----------------------------------------------------------------------

      dimension finput(7)
      dimension kind(nomp),eray(nomp),aray(nomp),bray(nomp),gray(nomp)

      dimension icc(12)
      data ifive,isix,zro,ne/5,6,0.,0/
      data eps9/1.0e-9/

*-----------------------------------------------------------------------
c//////////////
c      write(0,*) ' bert: finput'
c      write(0,*) finput
c///////////////

      amasno = finput(1)
      zee =    finput(2)
      einc =   finput(3)
      ctofe =  finput(4)
      casesn = finput(5)
      andit =  finput(6)
      prtin =  finput(7)
      ke = 0
      i18=0
      i19=0
c*pre
      npcle=1
      nhole=0
      call qfermi(efermi)
c*preq

      call zero

      do 11 i=1,12
   11 icc(i)=0
      space(13)=einc
      no=int(amasno)
      nmas=1+(no-1)*10

      do 33 i=1,10
      out(i)=crsc(nmas)
   33 nmas=nmas+1
   60 call rout1
      if(prtin-4.)105,105,135
  105 no= prtin+ 1.
  110 value1=einc+space(4)
      go to (1000,1000,121,3015,121),no
 3015 ne=3
      go to 135
  121 call rout2(tapcrs(1760))
      if(i1)200,3020,200
 3020 ne=4
      go to 135
  200 if(no-4)3005,3025,202
 3025 ne=5
      go to 135
  202 iv=0
  203 call xyi(1,14,30)
      ip=1
  201 call rout3
      if(begru)220,5090,220
  220 kk=i1
      xinc=xi(1)
      call rout4
      if(i1)1117,224,224
  224 i1=kk
  225 if(in)4480,226,4480
  226 if(ex-d(2))230,230,560
  230 curr(2)=out(13)
      wkrpn(3)=curr(2)
      wkrpn(6)=out(16)
      ifca=0
  231 call bg6ca(3,3)
  234 ifcc=0
23400 call abran(6)
      not=not
      knot=not
      go to (235,235,235,315,235,235),not
  235 call bg6c(isw(11))
      value1=rlke
      if(in)4830,23500,4830
23500 if(not-4)23501,315,23502
23501 any=space(not+13)
      go to 236
23502 any=s(not-4)
23503 if(not-5)236,540,550
  236 call rout5(tapcrs(5866),tapcrs(5740),tapcrs(5614))
  254 if(dabs(clsm-2.0).lt.eps9)go to 580
      if(clsm-2.0)640,580,255
  255 if(value1-value2)260,260,290
  260 if(isw(1))275,265,275
  265 ifc=ifcc+1
      if(in)4500,269,4500
  269 c(3)=0.0
  270 c(1)=curr(4)
      c(2)=curr(5)
      c(3)=c(3)+curr(6)+ex+d(1)
  274 it=it
      go to ( 670, 905, 910, 930, 955, 975, 980, 990, 995, 980,
     1       2000,4240,4245, 995,4250,4280,4285,4325,4330,4360,
     2       4365,4410,5050,5070,5075,5080, 980,5085),it
  275 if(isw(2))285,280,285
  280 ifc=2+ifcc
  281 if(in)4515,284,4515
  284 c(3)=d(2)+d(3)
      go to 270
  285 ifc=3+ifcc
      if(in)480,289,480
  289 c(3)=d(2)+d(3)+d(4)+d(5)
      go to 270
  290 call signex
  291 if(isw(1))295,225,295
  295 if(in)4535,296,4535
  296 if(ex-d(6))230,230,300
  300 if(isw(2))310,305,310
  305 ipec(7)=ipec(7)+1
      go to 201
  310 ipec(11)=ipec(11)+1
      go to 201
  345 i3=1
      go to 316
  315 i3=-1
  316 call rout6
      if(i3)3110,350,365
 3110 ne=22
      go to 135
  350 call rout6a
      if(dabs(clsm-2.0).lt.eps9) go to 660
      if(clsm-2.0)655,660,4870
  356 ifca=1
  357 if(isw(1))360,269,360
  360 if(isw(2))289,284,289
  365 call rout7
      if(i3)3110,425,425
  425 call rout7a
      i3=i3
      go to (480,3030,535,510,269,289,485,490,284),i3
 3030 ne=6
      go to 135
  480 value1=ex+d(4)+d(5)
  484 if(curr(10)-2.0)485,490,490
  485 c(1)=value1*curr(7)+curr(4)
      c(2)=value1*curr(8)+curr(5)
      c(3)=value1*curr(9)+curr(6)
      go to 274
  490 value1=value1+d(3)
      go to 500
  500 if(curr(10)-2.0)485,485,505
  505 value1=value1+d(2)
      go to 485
  510 if(inc)515,525,515
  515 c(3)=d(2)
      if(isw(3))520,270,520
  520   c(3)=c(3)+d(3)+d(4)
      go to 270
  525 if(isw(3))530,500,530
  530 value1=ex+d(4)
      go to 484
  535 value1=ex
      if(inc)284,484,284
  540 if(rlke-2500.0)5400,5400,3130

cKN
 3130 rlke = 2500.0d0
      write(*,*) '*** Warning ; BERT26, rlke is greater than 2500'
      goto 5400
c3130 ne=26
cKN   go to 135

 5400 if(rlke-180.0)541,541,545
  541 call signex
      if(dabs(clsm-2.0).lt.eps9) go to 601
      if(clsm-2.0)629,601,291
  545 value1=value1-180.0
      call crjab(1,tapcrs(4091))
      go to 254
  550 if(rlke-2500.0)554,554,3135

cKN
 3135 rlke = 2500.0d0
      write(*,*) '*** Warning ; BERT27, rlke is greater than 2500'
      goto 554
c3135 ne=27
cKN   go to 135

  554 if(rlke-180.0)541,541,555
  555 value1=value1-180.0
      call crjab(1,tapcrs(4325))
      go to 254
  560 if(d(3))570,565,570
  565 ipec(2)=ipec(2)+1
      go to 201
  570 isw(1)=1
      call spac32(31)
  574 if(in)4825,5740 ,4825
 5740 if(ex-d(3))575,575,615
 5750 if(in)576,575,576
  575 curr(2)=out(14)
      wkrpn(2)=curr(2)
      wkrpn(5)=out(17)
  576 call bg6ca(2,2)
      go to 234
  585 iv=-1
      go to 581
  580 iv=0
  581 call rout8
      i3=i3
      go to (601,4640,270,530),i3
  601 if(isw(3))605,574,605
  605 if(ex-d(5))5750 ,5750 ,610
  610 if(in)4840,6100 ,4840
 6100 call spac32(32)
  611 if(ex-d(6))230,230,310
  615 if(d(4))625,620,625
  620 call spac32(32)
  621 if(ex-d(6))230,230,305
  625 isw(2)=1
      isw(3)=1
      call spac32(30)
  629 if(in)4855,6290 ,4855
 6290 call rou10
      if(i3)605,636,636
  636 call bg6ca(1,1)
      go to 234
  640 if(value1-value2)645,645,650
  645 ifc=9+ifcc
      if(in)4620,284,4620
  650 call signex
      if(in)4855,6290 ,4855
  655 if(in)657,656,657
  656 ifca=6
      go to 284
  657 ifca=9*iabs(i6-2)+13*(i6-1)*(3-i6)
      go to 4620
  660 if(in)4885,661,4885
  661 ifca=7
      go to 515
  670 i3=1
      go to 673
 6700 i3=4
      go to 673
  671 i3=2
      go to 673
  672 i3=3
  673 call rou11(tapcrs(248))
      i3=i3
      go to (865,925,970,677),i3
  677 cst=crdt(2)- abs(snt*(crdt(2)-crdt(1)))
  680 snt=dsqrt(1.0-cst*cst)
  681 call rou12(ke)
      if(i3)3035,685,710
 3035 ne=7
      go to 135
  685 if(dabs((efrn-value1)/value1).lt.eps9) go to 686
      if(efrn-value1)720,686,690
  686 if(efrn)720,690,690
  690 fcn=fcn+1.0
  695 iv=-1
      go to 696
  705 iv=0
  696 i1=0
      call rou13
      if(i3)425,7090,350
 7090 ifc=ifc
      go to ( 225, 621, 611,1119,1236,1251, 574, 605, 629,1219,
     1       1249,1244,4480,4535,4535,4611,4544,4631,4825, 605,
     2       4480,4535,4535,4855,4480,4535,4535,4825,4955,4855),ifc
  710 if(dabs((efrp-value1)/value1).lt.eps9) go to 711
      if(efrp-value1)720,711,715
  711 if(efrp)720,715,715
  715 fcp=fcp+1.0
      go to 695
  850 i3=0
      go to 721
  861 i3=-1
      go to 721
  720 i3=1
  721 call rou14(iqstep,ke)
      i3=i3
      go to (710,685,3040,1340,4415),i3
 3040 ne=8
      go to 135
  865 i3=1
      go to 866
  876 i3=2
      go to 866
  882 i3=3
      go to 866
  883 i3=4
  866 call rou15(tapcrs(3401),tapcrs(3231))
      i3=i3
      go to (915,960,680,3045,920,965),i3
 3045 ne=9
      go to 135
  905 pt(2)=5.0
      ik=it
      pt(14)=2.0
      go to 671
  910 pt(2)=5.0
      pt(14)=1.0
      ik=it
      go to 671
  915 i3=1
      go to 926
  920 i3=2
      go to 926
  925 i3=3
  926 call rou16(tapcrs(3616),tapcrs(3446),tapcrs(626))
  927 if(i3)876,882,677
  930 pt(2)=3.0
      pt(14)=2.0
      ik=3
      go to 671
  955 pt(14)=2.0
 9550 ik=it
  956 pt(2)=4.0
      pm(3)=poms
      go to 672
  960 if(ik-23)961,5055,961
  961 i3=1
      go to 972
  966 i3=2
      go to 972
  971 i3=3
  972 call rou16(tapcrs(3831),tapcrs(3661),tapcrs(1004))
      go to 927
  965 if(ik-23)966,5060,966
  970 if(ik-23)971,5065,971
  975 pt(14)=1.0
      go to 9550
  980 pt(2)=1.0
  985 pt(14)=2.0
  986 call polt(iands,cst,snt)
      go to 681
  990 pt(2)=2.0
      go to 985
  995 pt(2)=1.0
      pt(14)=1.0
      go to 986
 2000 isw(9)=0
      isw(10)=0
 2040 i3=0
      go to 2002
 2001 i3=-1
 2002 call rou17(tapcrs(4793),tapcrs(4910),tapcrs(5011),tapcrs(4559),
     1tapcrs(5128))
      if(i3)3050,2084,2055
 3050 ne=10
      go to 135
 2055 k=3
 2056 if(dabs(pt(k-1)-1.0).lt.eps9) go to 2060
      if(pt(k-1)-1.0)2085,2060,2085
 2060 if(pt(k))2070,2070,2065
 2065 if(pt(k)-efrp)2070,2070,2075
 2070 fcp=fcp+1.0
 2071 pm(4)=dncms
      go to 705
 2075 m=pt(k-1)
      if(pt(k)-eco (m))2080,2080,2081
 2080 pt(k)=0.0
      pnbc(m)=pnbc(m)+1.0
 2081 if(dabs(col(15)-1.0).lt.eps9) go to 4035
      if(col(15)-1.0)2084,4035,2082
 2082 if(dabs(col(15)-3.0).lt.eps9) go to 4070
      if(col(15)-3.0)4071,4070,2083
 2083 if(dabs(col(15)-5.0).lt.eps9) go to 4145
      if(col(15)-5.0)4080,4145,3055
 3055 ne=11
      go to 135
 2084 call collm(0)
      if(pt(38))4015,4010,4015
 2085 if(dabs(pt(k-1)-2.0).lt.eps9) go to 2090
      if(pt(k-1)-2.0)4225,2090,4225
 2090 if(pt(k))4000,4000,2095
 2095 if(pt(k)-efrn)4000,4000,2075
 4000 fcn=fcn+1.0
      go to 2071
 4010 i3=1
      go to 4036
 4015 i3=2
      go to 4036
 4070 i3=4
      go to 4036
 4071 i3=5
      go to 4036
 4080 i3=6
      go to 4036
 4035 i3=3
 4036 call rou18
      i3=i3
      k=iv
      go to (2056,3060,2081,4235,3061),i3
 3060 ne=12
      go to 135
 3061 i18=i18+1
      go to 2071
 4145 call rou19
      if(i3)3065,4235,3066
 3065 ne=13
      go to 135
 3066 i19=i19+1
      go to 2071
 4225 if(col(15)-1.0)2084,3070,3070
 3070 ne=14
      go to 135
 4240 i3=2
      go to 4236
 4245 i3=3
      go to 4236
 4250 i3=4
      go to 4236
 4280 i3=5
      go to 4236
 4285 i3=6
      go to 4236
 4325 i3=7
      go to 4236
 4241 i3=8
      go to 4236
 4235 i3=1
 4236 call rou20(tapcrs(133),tapcrs(53),tapcrs(1862),tapcrs(1),
     1 tapcrs(1798))
      i3=i3
      go to (3075,861,2001,985,883,680,986,681),i3
 3075 ne=15
      go to 135
 4360 i3=2
      go to 4341
 4365 i3=3
      go to 4341
 4410 i3=4
      go to 4341
 4330 i3=1
 4341 call rou21(tapcrs(2678),tapcrs(2839),tapcrs(2940),tapcrs(3101),
     1 tapcrs(2502))
      go to 2002
 4830 iv=2
      go to 4481
 4500 iv=3
      go to 4481
 4515 iv=4
      go to 4481
 4535 iv=5
      go to 4481
 4825 iv=6
      go to 4481
 4640 iv=7
      go to 4481
 4840 iv=8
      go to 4481
 4855 iv=9
      go to 4481
 4620 iv=10
      go to 4481
 4885 iv=11
      go to 4481
 4611 iv=12
      go to 4481
 4544 iv=13
      go to 4481
 4631 iv=14
      go to 4481
 4955 iv=15
      go to 4481
 4415 iv=16
      go to 4481
 4696 iv=17
      go to 4481
 4670 iv=18
      go to 4481
 4479 iv=19
      go to 4481
 4650 iv=20
      go to 4481
 4610 iv=21
      go to 4481
 4870 iv=22
      go to 4481
 5005 iv=23
      go to 4481
 4480 iv=1
 4481 call rou22(npcle,nhole,
     &           tapcrs(1760),tapcrs(1779),tapcrs(5362),tapcrs(5614),
     &           tapcrs(5488))
      iv=iv
      if(i1.lt.0)go to 1117
 4482 go to( 5045, 3080,  605,  850, 1141,  485,  500,  484, 1290, 1161,
     1       1157, 1270,  576,23503,  636,  356,  480,  530,11410,  231,
     2       5035, 1170, 3105),iv
 3080 ne=16
      go to 135
 3105 ne=21
      go to 135
 5010 abx = 1.0
      value1=unirn(dummy)
      if(value1-ppnda)5015,365,365
 5015 it=27
      med=med
      absec=-hvp(med)
      go to 345
 5020 if(rlke-2500.0)5021,5021,3125

cKN
 3125 rlke = 2500.0d0
      write(*,*) '*** Warning ; BERT25, rlke is greater than 2500'
      goto 5021
c3125 ne=25
cKN   go to 135

 5021 if(rlke-180.0)5030,5030,5025
 5025 if(not-6)5026,5026,5040
 5026 value1=rlke-180.0
      call crjab(1,tapcrs(4208))
      go to 1170
 5030 if(dabs(clsm-2.0).lt.eps9) go to 1315
      if(clsm-2.0)1335,1315,1180
 5035 if(not-6)5005,5020,5020
 5040 value1=rlke-180.0
      call crjab(1,tapcrs(4442))
      go to 1170
 5045 isw(11)=0
      go to 4696
 5050 pt(14)=1.0
      go to 9550
 5060 i3=2
      go to 5056
 5065 i3=3
      go to 5056
 5055 i3=1
 5056 call rou16(tapcrs(4046),tapcrs(3876),tapcrs(1382))
 5057 go to 927
 5070 pt(2)=3.0
      pt(14)=2.0
      ik=it
      go to 671
 5075 pt(14)=2.0
      ik=23
      go to 956
 5080 pt(2)=5.0
      go to 6700
 5085 isw(9)=2
      go to 4241
 5090 itote  =ipec(2)+ipec(7)+ipec(11)
      if(itote-1)7004,7002,7001
 7001 call  exit
 7002 nopart = -1
 7003 continue
      return
 7004 nopart = esps(1)
      if(nopart.le.60)go to 7007
      write(6,7006) nopart
 7006 format(' nopart =',i3,' at stmt 7005 in bert, reduced to 60.')
      nopart = 60
 7007 continue
      do 7110 ndex = 1,nopart
      klmn = 8*(ndex-1) + 1
      kind(ndex) = esps(klmn+1)-1.
      eray(ndex) = esps(klmn+2)
      aray(ndex) = esps(klmn+3)
      bray(ndex) = esps(klmn+4)
 7110 gray(ndex) = esps(klmn+5)
      go to  7003
c%%  135 write(6,7125) ne
c%% 7125 format(' fatal error number ',i3,
c%%     1       ' in bert.  execution terminated.')
  135 write(6,7125) ne
 7125 format(/' *** error message from s.bert ***'
     &/' fatal error condition was found in bert.'
     &/' fatal error code =',i5,'  please trace by the code.')
      call parastop( 843 )
 1000 value2=einc+space(12)
      if(value1-160.0)1005,1005,1015
 1005 space(33)=2.9e-24
      fmax(2)=2.9e-24
      space(34)=.4e-24
      fmax(1)=.4e-24
      do 1010 i=9,12
 1010 s(i)=0.0
      go to 1090
 1015 call bovera(value2,dncms,ans)
      if(value1-560.0)1020,1020,1055
 1020 s(11)=0.0
      s(12)=0.0
      if(value1-400.0)1030,1025,1025
 1025 s(9)=20.9e-27*ans
      s(10)=10.5e-27*ans
      space(44)=58.0e-27
      space(45)=27.0e-27*ans
      go to 1054
 1030 if(value1-300.0)1040,1035,1035
 1035 s(9)=19.1e-27*ans
      s(10)=10.0e-27*ans
      space(44)=0.11e-24
      space(45)=38.0e-27
      go to 1054
 1040 if(value1-200.0)1050,1045,1045
 1045 s(9)=15.0e-27*ans
      s(10)=7.2e-27*ans
      space(44)=0.395e-24
      space(45)=0.107e-24
      go to 1054
 1050 s(9)=2.3e-27*ans
      s(10)=1.0e-27*ans
      space(44)=1.06e-24
      space(45)=0.183e-24
 1054 space(33)=space(44)
      space(34)=space(45)
      no=no
      go to 1090
 1055 if(value1-3600.0)1060,1060,3085

cKN
 3085 value1 = 3600.0d0
      write(0,*) '*** Warning ; BERT17, value1 is greater than 3600'
      write(0,*) 'Target m, z, Progj E (MeV)=', finput(1:3)
      write(0,*) ' proj. code ityp=finput(7)+1=',finput(7)+1
      goto 1060
c3085 ne=17
cKN   go to 135

 1060 s(9)=21.2e-27*ans
      s(10)=10.6e-27*ans
      if(value1-800.0)1065,1070,1070
 1065 s(11)=4.1e-27*ans
      s(12)=3.8e-27*ans
      space(46)=40.0e-27*ans
      space(47)=25.0e-27*ans
      go to 1089
 1070 if(value1-1680.0)1075,1085,1085
 1075 s(11)=13.0e-27*ans
      s(12)=11.7e-27*ans
      space(46)=33.0e-27*ans
      space(47)=25.0e-27*ans
      go to 1089
 1085 space(46)=26.0e-27*ans
      space(47)=19.5e-27*ans
      s(11)=15.0e-27*ans
      s(12)=12.3e-27*ans
 1089 space(33)=space(46)
      space(34)=space(47)
      no=no
 1090 go to (1091,2100),no
 1091 iv=1
 1092 call xyi(9,33,41)
      ip=2
 1095 if(no-2)1105,1100,1100
 1100 isw(4)=0
      go to 1110
 1105 isw(4)=1
 1110 call undis
      if(begru)1115,5090,1115
 1115 abx = 0.0
      xinc=xi(1)
      inc=1
      curr(1)=no
      curr(3)=dncms
      call  geo
      if(i1)1117,1118,1118
 1117 ne=20
      go to 135
 1118 call partin
      call spac32(43)
 1119 if(ex-d(2))1120,1120,1205
 1120 wkrpn(3)=out(13)
      wkrpn(6)=out(16)
      curr(2)=wkrpn(6)
      if(isw(4))1130,1135,1130
 1130 curr(2)=wkrpn(3)
 1135 call bg6ca(3,0)
      ifca=3
 1140 ifcc=3
 1141 ka=6
11410 call abran(ka)
      knot=not+6
      if(in.gt.0)knot=knot+6
      if(knot.eq.17)go to 5010
 1145 call bg6c(isw(4))
      if(rlke.gt.0.)go to 1146
 3120 ne=24
      go to 135
 1146 value1=rlke
      if(in.ne.0)go to 4670
      if(not.ge.5)go to 1290
 1150 if(not-2)1155,1160,1270
 1155 any=space(33)
 1157 call crjab(1,tapcrs(6168))
      go to 1170
 1160 any=space(34)
 1161 call crjab(1,tapcrs(5992))
 1170 if(dabs(clsm-2.0).lt.eps9) go to 1310
      if(clsm-2.0)1330,1310,1175
 1175 if(value1-value2)260,260,1180
 1180 call signex
      if(isw(1))1184,1183,1184
 1183 if(in)11830,1119,11830
11830 if(curr(1)-2.0)4480,4480,4479
 1184 if(in)4535,1185,4535
 1185 if(ex-d(6))1120,1120,1190
 1190 if(isw(2))1200,1195,1200
 1195 ipec(7)=ipec(7)+1
      go to 1095
 1200 ipec(11)=ipec(11)+1
      go to 1095
 1205 if(d(3))1215,1210,1215
 1210 ipec(2)=ipec(2)+1
      go to 1095
 1215 isw(1)=1
      call spac32(42)
 1219 if(ex-d(3))1220,1220,1230
 1220 wkrpn(2)=out(14)
      wkrpn(5)=out(17)
      curr(2)=wkrpn(5)
      if(isw(4))1221,1225,1221
 1221 curr(2)=wkrpn(2)
 1225 call bg6ca(2,0)
      go to 1140
 1230 if(d(4))1240,1235,1240
 1235 call spac32(43)
 1236 if(ex-d(6))1120,1120,1195
 1240 isw(2)=1
      isw(3)=1
      call spac32(41)
 1244 if(ex-d(4))1255,1255,1245
 1245 call spac32(42)
 1249 if(ex-d(5))1220,1220,1250
 1250 call spac32(43)
 1251 if(ex-d(6))1120,1120,1200
 1255 wkrpn(1)=out(15)
      wkrpn(4)=out(18)
      curr(2)=wkrpn(4)
      if(isw(4))1260,1265,1260
 1260 curr(2)=wkrpn(1)
 1265 call bg6ca(1,0)
      go to 1140
 1270 if(rlke-3500.0)1274,1274,3115
 3115 ne=23
      go to 135
 1274 if(rlke-360.0)1325,1325,1275
 1275 value1=rlke-360.0
      if(in)1278,1276,1278
 1276 any=s(knot)
 1278 if(not-4)1280,1285,3090
 3090 ne=18
      go to 135
 1280 call crjab(1,tapcrs(1926))
      go to 1170
 1285 call crjab(1,tapcrs(2214))
      go to 1170
 1290 if(rlke-3500.0)1294,1294,3115
 1294 if(rlke-920.0)1325,1325,1295
 1295 value1=rlke-920.0
      if(not-6)1300,1303,3095
 3095 ne=19
      go to 135
 1300 if(in)1302,1301,1302
 1301 any=s(11)
 1302 call crjab(1,tapcrs(2084))
      go to 1170
 1303 if(in)1305,1304,1305
 1304 any=s(12)
 1305 call crjab(1,tapcrs(2372))
      go to 1170
 1310 if(value1.le.value2)go to 585
 1315 call signex
      if(in)4650,1320,4650
 1320 if(isw(3))1249,1219,1249
 1325 if(dabs(clsm-2.0).lt.eps9) go to 1315
      if(clsm-2.0)1335,1315,1180
 1330 if(value1-value2)645,645,1335
 1335 call signex
      if(in)4610,1244,4610
 1340 if(esps(1))1345,1341,1345
 1341 nwds=1
      go to 1349
 1345 nwds=esps(1)*8.0+1.5
c     total no. of words(escaping particles)
 1349 continue
      nor=nor+1
 1360 in=0
      go to (201,1095),ip
 2100 iv=-1
      go to 1092
 3005 iv=2
      go to 203
c1370 return
      end subroutine


************************************************************************
*                                                                      *
      subroutine abran(k1)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
      value1=unirn(dummy)
      value1=value1*sigq
      value2=0.0
      not=1
      do 5 i=2,k1
      value2=ce(i)+value2
      if(value2-value1)5,15,15
    5 not=not+1
   15 continue
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine qfermi(efermi)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      dimension qrdius(4),qroh(4),qu(4)
      real*8 qrdius,qroh,qu,up,down,efermi
      data qrdius(1),qroh(1),qu(1)/0.,0.,0./
c ----------------------------------------------------------------------
c
      na=10*(amasno-1.)+1

      do 100 i=1,3
      qrdius(i+1)=crsc(na+i)
      qroh(i+1)=crsc(na+3+i)*1d30
      qu(i+1)=crsc(na+6+i)
  100 continue
      up=0.
      down=0.
      do 200 i=2,4
      up=up+(qrdius(i)**3-qrdius(i-1)**3)*qroh(i)*qu(i)
      down=down+(qrdius(i)**3-qrdius(i-1)**3)*qroh(i)
  200 continue
      efermi=(zee**(5./3.)+(amasno-zee)**(5./3.))/amasno*up/down
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine zero
c ----------------------------------------------------------------------
c     zero clear routine
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common /abcom/ amasno,andit,ab(9)
      common /ccom/  crsc(2500),ctofe,casesn,c(110)
      common /dcom/  dncms,d(9)
      common /ecom/  einc,e(491)
      common /fcom/  f(23)
      common /hcom/  h(6)
      common /icom/  i(52)
      common /klmn/  ln,nrt,kmn(8)
      common /opcom/ part(6),pnms,poms,prtin,op(1551)
      common /rcom/  rcpmv,r
      common /scom/  sf,sqnm,s(221)
      common /tcom/  trsym,t(9)
      common /uvwxz/ zee,uvwx(18)

c ----------------------------------------------------------------------
c
      do 1 n=1,9
        ab(n)=0.d0
    1 continue
c
      do 2 n=1,110
        c(n)=0.d0
    2 continue
c
      do 3 n=1,9
        d(n)=0.d0
    3 continue
c
      do 4 n=1,491
        e(n)=0.d0
    4 continue
c
      do 5 n=1,23
        f(n)=0.d0
    5 continue
c
      do 6 n=1,6
        h(n)=0.d0
    6 continue
c
      do 7 n=1,52
        i(n)=0
    7 continue
c
      do 8 n=1,8
        kmn(n)=0
    8 continue
c
      do 9 n=1,1551
        op(n)=0.d0
    9 continue
c
      do 10 n=1,221
        s(n)=0.d0
   10 continue
c
      do 11 n=1,9
        t(n)=0.d0
   11 continue
c
      do 12 n=1,18
        uvwx(n)=0.d0
   12 continue
c
      r=0.d0
c
csk --------------------------------- jun. 30 1994 -- blok 1 --  end  --
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine rout1
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      common/rtcom/tapcrs(6600)
c ----------------------------------------------------------------------
c
      out(11)=zee
      value2=zee**6.6666667d-1
      do 100 i=5,7
      space(i+2)=out(i)*out(11)
  100 space(i+5)=(out(i+3)*value2)+7.0d+0
      out(12)=amasno-out(11)
      value2=out(12)**(6.6666667d-1)
      do 110 i=5,7
      space(i-4)=out(i)*out(12)
  110 space(i-1)=(out(i+3)*value2)+7.0d+0
c
      do 120 i=1,3
      hvn(i)=0.5d+0*space(i+3)
      hvp(i)=0.5d+0*space(i+9)
      awd(i)=hvn(i)+hvp(i)
      fvnp(i)=0.5d+0*awd(i)
      vnvp(i)=space(i+3)-space(i+9)
      pmac(i)=vnvp(i)-hvn(i)
      ppan(i)=-vnvp(i)-hvp(i)
      thpn(i)=hvp(i)-vnvp(i)
      ffptfn(i)=-vnvp(i)+fvnp(i)
      tffn(i)=space(i+9)-fvnp(i)
  120 tffp(i)=vnvp(i)+tffn(i)
      pppda=(2.0d+0*zee)/(zee+amasno-1.0d+0)
!!!!   correction for 0/0:  KK  Mar.13.2014<<<<<<<<<<
      if(out(12) == 0.) then
!           happens when zee=1, out(11) = 1;  amasno=1
         ppmda = 1.0
         ppnda = 1.0
         ppnna = 1.0
      else
         ppmda=(2.0d+0*out(12))/(amasno+out(12)-1.0d+0)
         ppnda=(2.0d+0*zee*out(12))/(amasno*amasno-amasno)
         ppnna=
     *  (out(12)*out(12)-out(12))/(out(12)*out(12)+zee*zee-amasno)
      endif
!!!!!!!!!>>>>>>>>>>

      k=15
      do 130 i=4,6
      out(k)=space(i+6)+einc
      out(k+3)=space(i)+einc
  130 k=k-1
      out(30)=(zee/out(4))*1.4412d-13
      if (ctofe) 150,140,150
  140 ctofe=out(30)
  150 do 160 i=1,3
      cfepn(i+3)=space(i+3)+tapcrs(1)
  160 cfepn(i)=space(i+9)+ctofe
      in=0
      value1=6.28318531d+10*(1.19366207d-1)**3.3333333d-1
      do 170 i=1,3
      fmpn(i)=value1*(space(i+6))**3.3333333d-1
  170 fmpn(i+3)=value1*(space(i))**3.3333333d-1
c
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine rout2(t)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      dimension t(19)
c ----------------------------------------------------------------------
c
      i1=1
      value2=einc+space(12)
      call bovera(value2,pnms,ans)
      space(14)=.205e-24*ans
      space(15)=.0219e-24*ans
      space(16)=.0451e-24*ans
      if(value1-100.0)125,125,130
  125 fmax(1)=space(14)
      fmax(2)=space(15)
      fmax(3)=space(16)
      s(1)=0.0
      s(2)=0.0
12500 call crdet(1,t(1),einc)
       space(17)=crdt(1)
      if(i1)200,200,12501
12501 fmax(4)=space(17)
  200 return
  130 if(value2-2600.0)140,140,135
  135 i1=0
      go to 200
  140 space(17)=0.0
      if(value2-220.0)141,141,142
  141 s(1)=1.5e-27*ans
      s(2)=5.5e-27*ans
      go to 146
  142 if(value2-400.0)145,145,160
  145 space(14)=.205e-24
      space(16)=.0451e-24
      s(1)=5.5e-27*ans
      s(2)=20.0e-27*ans
  146 if(einc-360.0)12500,200,200
  160 if(value2-500.0)165,165,175
  165 space(14)=.113e-24
      space(15)=16.6e-27*ans
      space(16)=27.0e-27
      s(1)=7.6e-27*ans
      s(2)=20.3e-27*ans
      i1=-1
      go to 146
  175 s(2)=23.4e-27*ans
      space(15)=28.0e-27*ans
      if(value2-600.0)180,180,185
  180 space(14)=60.0e-27
      space(16)=17.5e-27
      s(1)=10.0e-27*ans
      go to 200
  185 if(value2-800.0)190,190,195
  190 space(14)=33.0e-27
      space(16)=12.0e-27*ans
      s(1)=22.0e-27*ans
      go to 200
  195 space(14)=16.5e-27*ans
      space(16)=9.5e-27*ans
      s(1)=26.0e-27*ans
      s(2)=25.0e-27*ans
      go to 200
      end subroutine


************************************************************************
*                                                                      *
      subroutine rout3
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      if(no-4)205,210,210
  205 isw(11)=1
      go to 215
  210 isw(11)=0
  215 call undis
      inc=1
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine rout4
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      call  geo
      if(i1)10,5,5
    5 curr(3)=pnms
      curr(1)=no
      call partin
      call spac32(32)
   10 return
      end subroutine


************************************************************************
*                                                                      *
      subroutine rout5(t,b,r)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      dimension t(126),b(126),r(126)
c ----------------------------------------------------------------------
c
      if(not-2)240,245,250
  240 call crjab(1,t(1))
      go to 254
  245 call crjab(1,b(1))
      go to 254
  250 call crjab(1,r(1))
  254 return
      end subroutine


************************************************************************
*                                                                      *
      subroutine rout6
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      if(i3)315,345,345
  315 abx=1.0
      med=clsm
      knot=not
      value1=unirn(dummy)
      if(isw(11))325,320,325
  320 if(value1-ppmda)330,365,365
  136 return
  365 i3=1
      go to 136
  325 if(value1-pppda)330,365,365
  330 if(isw(11))340,335,340
  335 it=13
      absec=pmac(med)
       go to 345
  340 it=14
      absec=-hvn(med)
  345 strkp=-1.0
      i1=0
      i2=med
      call bb
      i3=0
      go to 136
      end subroutine


************************************************************************
*                                                                      *
      subroutine rout6a
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
  350 i1=0
      call spisom
      strkp=-2.0
      i1=1
      call spisom
      strkp=-1.0
      com=(awd(med)-7.0)*2.0*rcpmv
      if(com-e(2))350,350,355
  355 pm(2)=2.0*dncms
      pm(3)=dncms
      e(2)=pm(2)+e(2)
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine rout7
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
       data eps9/1.0e-9/
c ----------------------------------------------------------------------
c
      i3=0
      if(abs(curr(1)-3.0).lt.eps9) go to 380
      if(curr(1)-3.0)135,380,370
  370 if(abs(curr(1)-5.0).lt.eps9) go to 375
      if(curr(1)-5.0)385,375,135
  135 i3=-1
  136 return
  375 it=7
      ifca=5
      absec=pmac(med)
      go to 400
  380 it=10
      ifca=3
       absec=ppan(med)
      go to 405
  385 value1=unirn(dummy)
      if(value1-ppnna)390,395,395
  390 it=8
      ifca=4
      absec=-hvn(med)
      go to 405
  395 it=9
      ifca=2
      absec=-hvp(med)
  400 strkp=-1.0
      e(1)=wkrpn(med)*rcpmv+pm(1)
      go to 410
  405 strkp=-2.0
      e(1)=wkrpn(med+3)*rcpmv+pm(1)
  410 if(inc)420,415,420
  415 call p1clc
      go to 136
  420 call p1cli
      go to 136
      end subroutine


************************************************************************
*                                                                      *
      subroutine rout7a
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
  425 i1=-1
      call spisom
      go to (430,435,430,430,435),ifca
  430 value1=space(med+3)-7.0
      go to 440
  435 value1=space(med+9)-7.0
  440 if(value1)135,445,445
  135 i3=2
  136 return
  445 if((value1*2.0*rcpmv)-e(2))425,425,450
  450 pm(3)=dncms
      pm(2)=2.0*dncms
      e(2)=pm(2)+e(2)
      value1=ex
       if(med-2)535,510,455
  535 i3=3
      go to 136
  510 i3=4
      go to 136
  455 if(inc)460,470,460
  460 if(isw(1))465,269,465
  269 i3=5
      go to 136
  465 if(isw(2))289,284,289
  289 i3=6
      go to 136
  284 i3=9
      go to 136
  470 if(isw(1))475,485,475
  485 i3=7
      go to 136
  475 if(isw(2))480,490,480
  480 i3=1
      go to 136
  490 i3=8
      go to 136
      end subroutine


************************************************************************
*                                                                      *
      subroutine rout8
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      i3=1
      if(iv)585,580,580
  580 if(value1-value2)585,585,600
  585 if(isw(3))595,590,595
  590 ifc=7+ifcc
      if(in)4640,591,4640
 4640 i3=2
 4641 return
  591 c(3)=d(2)
      go to 270
  270 i3=3
      go to 4641
  595 ifc=8+ifcc
      if(in)530,596,530
  530  i3=4
      go to 4641
  596 c(3)=d(2)+d(3)+d(4)
      go to 270
  600 call signex
      go to 4641
      end subroutine


************************************************************************
*                                                                      *
      subroutine rou10
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      i3=0
      if(ex-d(4))635,635,630
  630 call spac32(31)
      go to 605
  605 i3=-1
  606 return
  635 curr(2)=out(15)
      wkrpn(1)=out(15)
      wkrpn(4)=out(18)
      go to 606
      end subroutine


************************************************************************
*                                                                      *
      subroutine rou11(t)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      dimension t(378)
c ----------------------------------------------------------------------
c
      go to (670,671,672,6700),i3
  670 pt(2)=3.0
 6700 ik=it
      pt(14)=1.0
  671 pm(3)=pnms
  672 if(340.0-rlke)865,675,675
  865 i3=1
  866 return
  675 call mud(iands   ,snt,inpt)
      if(ik-3)676,925,970
  925 i3=2
      go to 866
  970 i3=3
       go to 866
  676 call crdet(21,t(1),rlke)
      i3=4
      go to 866
      end subroutine


************************************************************************
*                                                                      *
      subroutine rou12(ke)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      data eps9/1.0e-9/
c ----------------------------------------------------------------------
c
      i3=0
      call azit(iands   ,sopc,sops)
      call coll(0)
      if(col(15))135,682,135
  135 i3=-1
  136 return
  682 if(ke)135,683,681
  681 com = (e(4)-dncms)/rcpmv
      go to 684
  683 call  cole4
  684 i1= -1
      value1=com
       if(pt(14)-2.0)710,136,136
  710 i3=1
      go to 136
      end subroutine


************************************************************************
*                                                                      *
      subroutine rou13
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
       data eps9/1.0e-9/
c ----------------------------------------------------------------------
c
      i3=0
      if(iv)695,705,705
  695 if(abx)700,705,700
  700 if(ifca-2)701,703,702
  701 in=0
  350 i3=1
      go to 426
  705 call signex
      if(ifc-12)706,706,707
  702 if(ifca-6)703,701,710
  703 in=0
  425 i3=-1
  426 return
  710  if(ifca-8)701,704,704
  704 in=1
      go to  350
  706 in=0
      go to 426
  707 if(ifc-18)708,708,709
  708 in=-1
      go to 426
  709 in=1
      go to 426
      end subroutine


************************************************************************
*                                                                      *
      subroutine rou14(iqstep,ke)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      data eps9/1.0d-9/
c ----------------------------------------------------------------------
c
      if(i3) 861,850,720
  720 if(i1) 725,740,740
  725 i1=0
      value1=(e(3)-pm(3))/rcpmv
      if(abx) 730,735,730
  730 value1=value1+absec
  735 if(abs(pt(2)-2.0d+0).lt.eps9) go to 685
      if(pt(2)-2.0d+0) 710,685,740
  685 i3=2
  686 go to 1000
  710 i3=1
      go to 686
  740 call pinst
      if(i1)135,800,135
  135 i3=3
      go to 686
  800 i1=0
      m=pt(2)
      value2=value1
      if(m-3) 805,820,820
  805 continue

      if (iqstep.eq.3) then
          qeco1=unirn(dummy)
          qeco2=unirn(dummy)
          qeco =dmin1(qeco1, qeco2)
          qeco = eco(m) + 40.d+0 * qeco
          if (qeco-value2) 815, 806, 810
      else
           if (eco(m)-value2) 815, 806, 810
      endif

  806 if(eco(m)) 815,810,810
  810 pt(i1+3)=0.0d+0
      pnbc(m)=pnbc(m)+1.0d+0
      go to 834
  815 pt(i1+3)=value2
      if (i1) 860,835,860
c%%<02/07/97: avoid the machine dependence
c%%  820 if(abs((value2-clcfe)/clcfe).lt.eps9) go to 821
  820 if (clcfe.ne.0.0d+0) then
	 if (abs((value2-clcfe)/clcfe).lt.eps9) go to 821
      endif
c%%>
      if (value2-clcfe) 810,821,825
  821 if (value2) 810,825,825
  825 if (strkp+2.0) 830,830,815
  830 pt(3)=value1-space(med+3)+space(med+9)
  834 if (i1) 845,835,845
  835 m=pt(14)
      if (m-3) 840,135,135
  840 value2=com
      i1=12
      go to 805
  845 if (pt(3)) 860,850,860
  850 call punp
      if (i1) 135,1340,4415
 1340 i3=4
      go to 686
 4415 i3=5
      go to 686
  860 call collm(-1)
  861 if(ke.gt.0)go to 1000
      call stpr
      if (i1) 135,850,135
 1000 return
      end subroutine


************************************************************************
*                                                                      *
      subroutine rou15(t,b)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      dimension t(45),b(170)
c ----------------------------------------------------------------------
c
      go to (865,876,882,883),i3
  865 call snn
      if(i1)135,870,880
  135 i3=4
  136 go to 1000
  870 if(ik-3)875,915,960
  915 i3=1
      go to 136
  960 i3=2
      go to 136
  875 call dcintp(t(1))
  876 if(i2)680,680,135
  680 i3=3
       go to 136
  880 value1=rlke-340.0
      if(ik-3)881,920,965
  920 i3=5
      go to 136
  965 i3=6
      go to 136
  881 call crdet(5,b(1),value1)
  882 i1=6
      i2=0
  883 value1=unirn(dummy)
      if(value1-crdt(1))900,885,885
  885 value2=1.0
      value1=unirn(dummy)
      if(value1-crdt(4))890,895,895
  890 value1=unirn(dummy)
  891 cst=value1*value2
      go to 680
  895 com=0.0
      call big7t(iands   ,com,value1,i1)
      go to 891
  900 value2=-1.0
      value1=unirn(dummy)
      i1=i1+i2
      if(value1-crdt(2))890,895,895
 1000 return
      end subroutine


************************************************************************
*                                                                      *
      subroutine rou16(t,b,r)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      dimension t(45),b(170),r(378)
c ----------------------------------------------------------------------
c
      go to (915,920,925),i3
  915 call dcintp(t(1))
      i3=-1
  876 return
  920 call crdet(5,b(1),value1)
  882 i3=0
      go to 876
  925 call crdet(21,r(1),rlke)
      i3=1
      go to 876
      end subroutine


************************************************************************
*                                                                      *
      subroutine rou17(t,b,r,w,g)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      dimension t(117),b(101),r(117),w(234),g(234)
      data eps9/1.0e-9/
c ----------------------------------------------------------------------
c
      if(i3)2001,2001,2040
 2001 pt(38)=0.0
      value1=rlke-180.0
      call crdet(1,t(1),value1)
      com2=crdt(1)
      ftr=dncms*rlke*2.0*rcpmv+2.9877156e27
      univer=dsqrt(ftr)
 2005 value2=unirn(dummy)
      com=value2*com2
      call gene(b(1))
      com1=(com*com+ftr-.501264e26)/(2.0*univer)
      a=com1*com1-com*com
       if(a)2006,2009,2009
 2006 pacnt=pacnt+1.0
      go to 2005
 2009 unive=((univer-com1)*com1/univer)*dsqrt(a)
      call crdet(1,r(1),value1)
      com1  =unirn(dummy)
      if((unive/crdt(1))-com1)2005,2010,2010
 2010 call angid
      pm(3)=com
      pm(4)=poms
      pt(2)=3.0
      pt(4)=poms
      pt(14)=3.0
      pt(16)=poms
      pt(26)=1.0
      pt(28)=dncms
      if(isw(9))2020,2015,2020
 2015 if(isw(10))2030,2025,2030
 2020 if(isw(10))2035,135,2035
  135 i3=-1
  136 return
 2025 value1=.4
      value2=6.6666667e-1
      value3=0.0
      go to 2037
 2030 call crdet(2,w(1),value1)
      value3=3.3333333e-1
      go to 2036
 2035 call crdet(2,g(1),value1)
      value3=strkp
 2036 value1=crdt(1)
      value2=crdt(2)
 2037 call alpha
 2040 call ecpl
      if(i1)2045,2045,135
 2045 call coll(-1)
      if(col(15))135,2050,135
 2050 if(pt(38))2084,2055,2084
 2084 i3=0
      go to 136
 2055 pt(39)=0.0
      pt(3)=((e(4)-pm(4))/rcpmv)+pt(3)
      i3=1
      go to 136
      end subroutine


************************************************************************
*                                                                      *
      subroutine rou18
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
      ctofen=0  !cKN
c
      go to (4010,4015,4035,4095,4071,4080),i3
 4010 i=3
      col(15)=1.0
      k=27
      go to 4020
 4015 i=3
      col(15)=4.0
      k=15
 4020 pnidk(1)=pm(i)
      j=i
      do4025l=2,4
      pnidk(l)=pxyz(j)
 4025 j=j+4
       pnidk(5)=e(i)
      pnidk(6)=pt(k-11)
      call idk
      if(k-27)4031,4030,4031
 4030 pt(15)=pt(15)+((pnidk(12)-pnidk( 6))/rcpmv)
 4031 pt(k)=pt(k)+((pnidk(13)-dncms)/rcpmv)
      i3=1
 2057 iv=k
      return
 4035 k=3
      col(15)=2.0
      if(pt(2)-3.0)4071,4039,4039
 4039 if(pt(k)-2500.0)4040,4040,4055
 4055 i3=5
      go to 2057
 4040 if(pt(k))4050,4050,4045
 4045 ccofe =  eco(1)
      if (pt(k-1) - 4.0) 4047,4046,4046
 4046 ccofe = ccofe - ctofe + ctofen
 4047 if (pt(k) - ccofe ) 4050,4050,4070
 4050 m=pt(k-1)
      pnbc(m)=pnbc(m)+1.0
      pt(k)=0.0
      i3=3
      go to 2057
 4070 if(k-3)4071,4071,4095
 4071 col(15)=3.0
      k=15
      if(pt(14)-2.0)135,135,4039
  135 i3=2
      go to 2057
 4080 l=14
      do4085m=5,7
      pt(m)=pnidk(l)
      pt(m+12)=pnidk(l+3)
 4085 l=l+1
      pt(11)=pnidk(12)
      pt(12)=pnidk(6)
      i=4
      k=39
      col(15)=5.0
      go to 4020
 4095 i1=3
 4100 k=12*i1-33
      if(i1-4)4110,4115,4120
 4110 i2=-1
      go to 4130
 4115 i2=0
      go to 4130
 4120 if(i1-5)4115,4125,4235
 4235 i3=4
      go to 2057
 4125 i2=1
 4130 if(pt(k))4135,4140,4135
 4135 call pstor
 4140 i1=i1+1
      go to 4100
      end subroutine


************************************************************************
*                                                                      *
      subroutine rou19
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      data eps9/1.0e-9/
c ----------------------------------------------------------------------
      ctofen=0  !cKN
c
      pt(3)=pt(3)+((pt(11)-pt(12))/rcpmv)
      k=3
 4135 if(pt(k)-2500.0)4150,4150,4140
 4140 i3=1
      go to 136
 4150 if(pt(k))4160,4160,4155
 4155 ccofe =  eco(1)
      if (pt(k-1)-4.0) 4157,4156,4156
 4156 ccofe = ccofe - ctofe + ctofen
 4157 if(pt(k)-ccofe) 4160,4160,4200
 4160  pt(k)=0.0
      if(abs(pt(k-1)-3.0).lt.eps9) go to 4170
      if(pt(k-1)-3.0)135,4170,4165
  135 i3=-1
  136 return
 4165 if(pt(k-1)-5.0)4170,4170,135
 4170 m=pt(k-1)
      pnbc(m)=pnbc(m)+1.0
      go to 4185
 4175 i2=2
 4176 i1=(k/12)+3
      call pstor
 4185 if(k-15)4190,4195,4210
 4190 k=15
      if(pt(15))4195,4195,4175
 4195 k=27
      pt(27)=pt(27)+((pnidk(12)-pt(k+1))/rcpmv)
      go to 4135
 4200 if(k-15)4175,4205,4205
 4205 i2=0
      go to 4176
 4210 if(k-27)135,4215,4235
 4215 if(pt(39))4235,4235,4220
 4235 i3=0
      go to 136
 4220 i2=1
      k=39
      go to 4176
      end subroutine


************************************************************************
*                                                                      *
      subroutine rou20(t,b,r,w,g)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      dimension t(115),b(80),r(64),w(52),g(64)
c ----------------------------------------------------------------------
c
      go to (4235,4240,4245,4250,4280,4285,4325,4241),i3
 4235 pm(4)=dncms
      call pinst
      if(i1)135,861,135
  135 i3=1
  136 return
  861 i3=2
      go to 136
 4240 isw(9)=0
 4241 isw(10)=2
      i3=3
      go to 136
 4245 pt(2)=2.0
       i3=4
      go to 136
 4250 pt(2)=2.0
      pt(14)=1.0
 4251 pm(3)=dncms
      if(740.0-rlke)4270,4255,4255
 4255 if(300.0-rlke)4260,4265,4265
 4260 value1=rlke-300.0
      call crdet(5,t(1),value1)
      i1=3
      i2=3
  883 i3=5
      go to 136
 4265 call crdet(5,b(1),rlke)
      i1=3
      i2=1
      go to 883
 4270 if(3500.0-rlke)135,4271,4271
 4271 if(it-17)4275,135,135
 4275 call dcpr(r(1))
  679 if(i1)135,680,680
  680 i3=6
      go to 136
 4280 pt(2)=1.0
      pt(14)=2.0
      go to 4251
 4285 pt(2)=2.0
      pt(14)=2.0
 4286 pm(3)=dncms
      if(500.0-rlke)4290,986,986
  986 i3=7
      go to 136
 4290 if(1000.0-rlke)4315,4295,4295
 4295 value1=rlke-500.0
      call crdet(2,w(1),value1)
      value1=unirn(dummy)
      if(crdt(1)-value1)4310,4300,4300
 4300 cst=unirn(dummy)
 4301 snt=dsqrt(1.0-cst*cst)
      value1=unirn(dummy)
      if(0.5-value1)4305,681,681
  681 i3=8
      go to 136
 4305 cst=- abs(cst)
      go to 681
 4310 com=0.0
      call big7t(iands   ,com,cst,3)
      go to 4301
 4315 if(3500.0-rlke)135,4320,4320
 4320 call dcpr(g(1))
c     (p-p)diff.crs.sec.high en.
      go to 679
 4325 pt(2)=1.0
      pt(14)=1.0
      go to 4286
c     no pion production possible
      end subroutine


************************************************************************
*                                                                      *
      subroutine rou21(v,w,x,y,z)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      dimension v(161),w(101),x(161),y(130),z(176)
c ----------------------------------------------------------------------
c
      value2=rlke*4.81633308e24+9.0554256e27
      value3=dsqrt(value2)
      go to (4330,4360,4365,4410),i3
 4330 isw(12)=0
 4331 pt(38)=0.0
      i1=0
      ans=rlke
 4333 value1=ans-300.0
      call crdet(1,v(1),value1)
      ftr=crdt(1)
 4335 sn    =unirn(dummy)
      com=sn*ftr
      call gene(w(1))
       if(i1)4370,4336,4375
 4336 com1=(com*com-sqnm+value2)/(2.0*value3)
      a=com1*com1-com*com
      if(a)4337,4338,4338
 4337 pgcnt=pgcnt+1.0
      go to 4335
 4338 univer=dsqrt(a)*com1*(value3-com1)/value3
      call crdet(1,x(1),value1)
      com1=unirn(dummy)
      if(com1-(univer/crdt(1)))4340,4340,4335
 4340 pm(4)=dncms
      pm(3)=com
      call angid
      pt(4)=dncms
      pt(28)=dncms
      call alp19
 2040 return
 4360 isw(12)=2
      go to 4331
 4365 isw(13)=0
 4366 i1=-1
      ans=((value3-pnms)**2-9.0554256e27)/4.81633308e24
      go to 4333
 4370 com1=((value3+dncms-com)**2-9.0554256e27)/4.81633308e24
      com2=com
      ans=com1
      com4=ftr
      i1=1
      go to 4333
 4375 com1=(com2*com2-com*com+value2)/(2.0*value3)
      a=com1*com1-com2*com2
      if(a)4376,4377,4377
 4376 pecnt=pecnt+1.0
      go to 4380
 4377 univer=dsqrt(a)*com1*(value3-com1)/value3
      value1=rlke-920.0
      call crdet(1,y(1),value1)
      sn    =unirn(dummy)
      if(sn-(univer*ftr/(crdt(1)*com4)))4385,4385,4380
 4380 ftr=com4
      i1=-1
      go to 4335
 4385 value1=unirn(dummy)
      if(value1-.5)4390,4390,4395
 4390 pm(3)=com2
      pm(4)=com
      go to 4400
 4395 pm(3)=com
      pm(4)=com2
 4400 call angid
      pt(16)=dncms
      pt(40)=dncms
      if(isw(13))4401,4405,4401
 4401 call crdet(1,z(1),rlke)
      value1=crdt(1)
 4405 pt(2)=3.0
      pt(4)=poms
      pt(14)=1.0
      pt(26)=3.0
      pt(28)=poms
      pt(38)=1.0
      call alp28
      go to 2040
 4410 isw(13)=2
      go to 4366
      end subroutine


************************************************************************
*                                                                      *
      subroutine rou22(npcle,nhole,v,w,x,y,z)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      dimension v(19),w(19),x(126),y(126),z(126)
      data eps9/1.0e-9/
c ----------------------------------------------------------------------
c
c//////////////////  ! to avoid undef var.
      integer,save::lg
c//////////////
      inmiv = iv
      i5=curr(1)
      go to (4480,4830,4500,4515,4535,4825,4640,4840,4855,4620,
     1       4885,4611,4544,4638,4955,4415,4696,4670,4479,4650,
     2       4610,4870,5005),iv
 4415 do4420i=1,3
      xi(i)=curr(i+3)
 4420 dco(i)=curr(i+6)
      in=-1
      med=curr(10)
      i5=curr(1)
      call  geo
      if(i1)5046,4424,4424
 4424 if(abs(curr(1)-2.0).lt.eps9) go to 4690
       if(curr(1)-2.0)4430,4690,4425
 4425 if(abs(curr(1)-4.0).lt.eps9) go to 4696
      if(curr(1)-4.0)4695,4696,5045
 5045 iv=1
c*preq
 5046 continue 
!      kkrou = kkrou + 1  ! this line has no meaning; kkrou nowhere used

      if (inmiv.eq.16) then

         npcle = npcle + 1
         nhole = nhole + 1

      endif

      return
c*preq
c5046 return
 4430 isw(4)=1
 4431 abx=0.0
      i4=1
      i2=med
      if(i2-1)135,4436,4470
 4436 i4=6
      go to 4470
 4438 isw(6)=1
      isw(5)=1
      i3=0
      if(com-3600.0)4440,4440,135
  135 iv=2
      go to 5046
 4440 if(com-560.0)4445,4445,4450
 4445 if(com-160.0)4460,4460,4455
 4450 call dfmax
      i1=6
 4451 if(curr(1)-2.0)4452,4457,4457
 4457 i3=1
 4452 call store
      ex=0.0
      call signex
      go to (4479,4535,4610,4630,4544,4665,4810,4825,4855,605,4955),i4
  605 iv=3
      go to 5046
 4455 isw(6)=0
      call pfmax
      i1=4
      go to 4451
 4460 isw(5)=0
      isw(6)=0
      fmax(1)=.4e-24
      fmax(2)=2.9e-24
      fmax(3)=0.
      fmax(4)=0.
      fmax(5)=0.
      fmax(6)=0.
      i1=2
      go to 4451
 4470 isw(1)=0
      isw(2)=0
      isw(3)=0
 4471 m=med+idint(15.0-6.0*curr(1) )
      a=curr(2)-space(m)
 4472 do4475i=1,3
      wkrpn(i)=a+space(i+9)
 4475 wkrpn(i+3)=a+space(i+3)
 4476 m=4-3*isw(4)
      com=wkrpn(m)
      go to 4438
 4479 go to (4890,4890,4480),i2
 4480 if(ex-d(2))4495,4495,4485
 4485 if(d(3))4660,4490,4660
 4490 call ccpes
      if(i1)850,850,135
  850 iv=4
      go to 5046
 4495 if(in) 4951, 4951,4965
 4951 call bg6ca(3,0)
      ifcc=12
 4496 med=clsm
      iv=5
      go to 5046
 4500 value1=ex
      iv=6
      go to 5046
 4515 value1=ex+d(3)
  500 iv=7
      go to 5046
 4535 if(ex-d(6))4495,4495,4490
 4540 isw(1)=1
 4543 go to (4544,4544,4544,4825),i5
 4544 if(ex-d(3))4635,4635,4545
 4545 i2=3
      i4=2
      i3=1
      i1=2
      if(d(4))4550,4476,4550
 4550 isw(2)=1
      isw(3)=1
      i4=3
      i2=1
      go to 4476
 4610 go to (4611,4611,4611,4855),i5
 4611 if(ex-d(4))4615,4615,4625
 4615 call bg6ca(1,0)
      ifcc=7
      go to 4496
 4620 value1=ex
      iv=8
      go to 5046
 4625 i4=4
      i1=2
      i2=2
      i3=1
      go to 4476
 4630 if(i5-2)4638,4638,4634
 4634 if(i5-4)135,4955,135
 4638 if(ex-d(5))4635,4635,4655
 4635 call bg6ca(2,0)
      ifcc=10
      go to 4496
 4640 value1=ex
      go to 500
 4650 if(isw(3))4630,4543,4630
 4655 i4=2
      i2=3
 4656 i3=1
      i1=2
      go to 4476
 4660 isw(1)=1
      if(abs(curr(1)-3.0).lt.eps9) go to 4815
      if(curr(1)-3.0)4661,4815,4930
 4661 i4=5
      i2=2
      go to 4656
 4665 isw(1)=1
      isw(2)=1
      isw(3)=1
      go to 4479
 4670 any=fmax(not)
      if(not-5)4684,1290,4680
 1290 iv=9
      go to 5046
 4680 if(i5-4)1290,4685,1290
 4684 if(knot-15)4685,5000,5000
 4685 go to (4686,4686,4686,4985),i5
 4686 if(not-2)1161,1157,1270
 1161 iv=10
      go to 5046
 1157 iv=11
      go to 5046
 1270 iv=12
      go to 5046
 4690 isw(4)=0
      go to 4431
 4695 isw(11)=1
 4696 in=1
      isw(1)=0
      isw(2)=0
      isw(3)=0
      isw(5)=1
      isw(6)=0
      isw(7)=0
      isw(8)=1
      i6=i5-2
      i2=med
      com=curr(2)-space(med+9)
      go to (4706,4705,4707),med
 4705 isw(8)=0
 4706 isw(7)=1
 4707 do4708i=1,3
      wkrpn(i)=com+space(i+9)
 4708 wkrpn(i+3)=com+space(i+3)
      com=com+space(4)
 4710 if(com-2600.0)4715,4715,135
 4715 if(com-100.0)4716,4716,4717
 4716 lg=4
      go to 4910
 4717 lg=6
 4720 call spcn
      if(value1)4725,4725,4734
 4725 com=curr(2)-space(med+9)
      if(com-360.0)4730,4734,4734
 4730 go to (4731,4915,4731),i6
 4731 call crdet(1,v(1),com)
      fmax(4)=crdt(1)
 4734 go to (4735,4920,4735),i6
 4735 i1=6
 4736 i4=7
      i2=1
      i3=0
      go to (4755,4760,4755),i6
 4755 if(isw(11))4760,4790,4760
 4760 if(isw(7))4765,4785,4765
 4765 if(isw(8))4452,4780,4452
 4780 i2=2
      go to 4452
 4785 i2=3
      go to 4452
 4790 i3=1
      go to 4760
 4810 if(curr(1)-2.0)4811,4811,4479
 4811 a=curr(2)-space(med+9)
      go to 4472
 4815 i4=8
 4816 i2=2
 4818 i3=0
      i1=lg
      if(abs(curr(1)-4.0).lt.eps9) go to 4827
      if(curr(1)-4.0)4819,4827,4819
 4819 if(isw(11))4821,4820,4821
 4820 i3=1
 4821 if(curr(1)-3.0)4822,4823,4823
 4822 iv=23
      go to 5046
 4823 if(lg-4)4822,4452,4824
 4824 m=5-iabs(i5-4)
      univer=fmax(m)
      call spcn
      fmax(m)=univer
      go to 4452
 4827 i1=i1+1
      go to 4823
 4825 if(ex-d(3))4826,4826,4845
 4826 go to (576,4940,576),i6
  576 iv=13
      go to 5046
 4830 any=fmax(not)
      if(curr(1)-3.0)4831,4832,4832
 4831 ifc=12
23503 iv=14
      go to 5046
 4832 ifcc=(clsm-2.0)*((clsm*5.5)-8.5)+12.05
      go to 23503
 4840 if(i4-10)4841,4865,4841
 4841 if(abs(clsm-2.0).lt.eps9) go to 4865
      if(clsm-2.0)4842,4865,4842
 4842 i4=2
      go to 4816
 4845 i4=2
      i2=3
      if(d(4))4850,4846,4850
 4846 go to (4818,4960,4818),i6
 4850 isw(2)=1
      isw(3)=1
      i4=9
      i2=1
      go to 4846
 4855 if(ex-d(4))4856,4856,4860
 4856 go to (636,4970,636),i6
  636 iv=15
      go to 5046
 4860 i4=10
      go to (4816,4950,4816),i6
 4865 i4=2
      i2=3
      go to 4818
 4870 if(in)4875,356,4875
  356 iv=16
      go to 5046
 4875 ifca=8*iabs(i6-2)-11*(i6-1)*(i6-3)
      if(isw(1))4880,4500,4880
 4880 if(isw(2))480,4515,480
  480 iv=17
      go to 5046
 4885 ifca=10*iabs(i6-2)+12*(i6-1)*(3-i6)
      if(isw(3))530,4640,530
  530 iv=18
      go to 5046
 4890 if(curr(1)-3.0)4895,4900,4900
 4895 go to (4610,4540),med
 4900 isw(1)=1
      go to (4905,4825),med
 4905 isw(2)=1
      isw(3)=1
      go to 4855
 4910 isw(5)=0
      go to (4911,4980,4911),i6
 4911 call bovera(curr(2),pnms,ans)
      fmax(1)=.203e-24*ans
      fmax(2)=21.9e-27*ans
      fmax(3)=45.1e-27*ans
      com=curr(2)-space(med+9)
      call crdet(1,v(1),com)
      fmax(4)=crdt(1)
      i1=4
      go to 4736
 4915 call crdet(1,w(1),com)
      fmax(5)=crdt(1)
 4920 i1=7
      go to 4736
 4925 call bg6ca(3,4)
      ifcc=24
11410 ka=7
      med=clsm
      iv=19
      go to 5046
 4930 i4=8
 4931 i2=2
      go to 4818
 4940 call bg6ca(2,3)
 4945 ifcc=21
      go to 11410
 4950 i4=11
      go to 4931
 4955 if(ex-d(5))4940,4940,4975
 4960 i1=9-i5
      go to 4818
 4965 go to (231,4925,231),i6
  231 iv=20
      go to 5046
 4970 call bg6ca(1,2)
      go to 4945
 4975 i1=5
      go to 4865
 4980 call bovera(curr(2),poms,ans)
      fmax(1)=82.9e-27*ans
      fmax(2)=45.1e-27*ans
      fmax(3)=fmax(1)
      space(48)=fmax(1)
      fmax(4)=fmax(2)
      space(49)=fmax(2)
      com=curr(2)-space(med+9)
      call crdet(1,w(1),com)
      fmax(5)=crdt(1)
      space(50)=fmax(5)
      i1=5
      go to 4736
 4985 if(not-2)4990,4995,5035
 5035 iv=21
      go to 5046
 4990 call crjab(1,x(1))
 1170 iv=22
      go to 5046
 4995 call crjab(1,y(1))
      go to 1170
 5000 if(knot-16)5005,4995,4995
 5005 call crjab(1,z(1))
      go to 1170
      end subroutine


************************************************************************
*                                                                      *
      subroutine xyi(ii,jj,kk)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      w1=s(ii)*1.0e30
      w2=s(ii+1)*1.0e30
      if(iabs(iv)-1)5,20,5
    5 w3=space(jj)*1.0e30
      w4=space(jj+1)*1.0e30
      w5=space(jj+2)*1.0e30
      w6=space(jj+3)*1.0e30
      ll=ii+2
      if(iv)15,10,15
   20 w3=s(ii+2)*1.0e30
      w4=s(ii+3)*1.0e30
      w5=space(jj)*1.0e30
      w6=space(jj+1)*1.0e30
      ll=ii+4
       if(iv)35,25,25
   25 nm=7
      i7=1
   15 mm=7
      nn=1
      go to 40
   35 nm=1
      i7=7
   10 mm=1
      nn=7
   40 do55i=1,3
      s(ll)=w1*space(mm)
      s(ll+1)=w2*space(nn)
      if(iabs(iv)-1)50,45,50
   45 s(ll+2)=w3*space(nm)
      s(ll+3)=w4*space(i7)
      ll=ll+2
      nm=nm+1
      i7=i7+1
   50 mm=mm+1
      nn=nn+1
   55 ll=ll+2
      if(iv)60,65,70
   60 nm=7
      i7=1
      go to 80
   65 mm=1
      nn=7
      nm=7
      i7=7
      go to  72
   70 if(iv-1)75,75,71
   71 mm=7
      nn=1
      nm=1
      i7=7
   72 ll=jj+4
      go to 85
   75 nm=1
      i7=7
   80 ll=jj
   85 do100i=1,3
      space(ll+2)=w5*space(nm)
      space(ll+3)=w6*space(i7)
      if(iabs(iv)-1)90,95,90
   90 space(ll)=w3*space(mm)
      space(ll+1)=w4*space(nn)
      ll=ll+2
      mm=mm+1
      nn=nn+1
   95 nm=nm+1
      i7=i7+1
  100 ll=ll+2
      ll=kk+2
      if(iabs(iv)-1)105,120,105
  105 mm=jj+4
      nn=ii+2
      do110i=kk,ll
      space(i)=space(mm)+space(mm+1)+space(mm+2)+space(mm+3)
     1+s(nn)+s(nn+1)
      mm=mm+4
  110 nn=nn+2
  115 return
  120 mm=ii+4
      nn=jj+2
      do125i=kk,ll
      space(i)=space(nn)+space(nn+1)+s(mm)+s(mm+1)+s(mm+2)+s(mm+3)
      mm=mm+4
  125 nn=nn+2
      go to 115
      end subroutine


************************************************************************
*                                                                      *
      subroutine bg6b
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      do 2 i=1,21
    2 ce(i)=0.
      j=i2
      do5i=2,i4
      ce(i)=space(j)
    5 j=j+1
      j=i4+1
      do10i=j,6
      ce(i)=s(i3)
   10 i3=i3+1
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine bg6c(int1)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      if(knot-7)5,81,81
    5 abx=0.0
      if(knot-2)10,15,15
   10 if(int1)40,50,40
   15 if(knot-5)20,25,30
   20 if(int1)50,40,50
   25 it=11
      if(int1)45,55,45
   30 it=12
      if(int1)55,45,55
   40 it=2*knot-1
   45 strkp=-1.0
   46 i1=0
      go to 60
   50  it=2*knot
   55 strkp=-2.0
   56 i1=1
   60 i2=clsm
      call bb
   65 call isom
      if(knot.eq.4.or.knot.eq.17)rlke=0.
      if(rlke.le.2500.)return
      if(rlke.gt.3500..or.knot.lt.7.or.knot.gt.12)go to 65
   75 return
   81 if(knot.gt.12)go to 135
      if(in.ne.0)go to 155
   86 if(knot-8)90,115,136
   90 if(int1)95,105,95
   95 it=2*(knot+1)
  100 strkp=-2.0
      go to 46
  105 it=2*knot+1
  110 strkp=-1.0
      go to 56
  115 if(int1)125,130,125
  125 it=2*(knot+1)
  126 go to 45
  130 it=2*knot+1
      go to 55
  135 if(knot-18)136,195,195
  136 it=knot+10
  137 if(knot-10)140,145,150
  140 if(int1)45,55,45
  145 if(int1)100,110,100
  150 if(knot-12)140,145,185
  155 if(knot-8)160,170,136
  160 if(int1)165,136,165
  165 it=knot+11
      go to 45
  170 if(int1.ne.0)go to 180
      it=knot+knot-1
      go to 110
  180 it=knot+knot
      go to 100
  185 abx=0.0
      if(knot-15)45,55,190
  190 if(knot-18)55,45,55
  195 it=28
      go to 185
      end subroutine


************************************************************************
*                                                                      *
      subroutine bg6ca(k,l)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      clsm=k
      if(in)20,15,20
   15 curr(10)=clsm
      curr(11)=clsm
   20 efrp=space(k+9)-7.0
      efrn=space(k+3)-7.0
      pm(2)=dncms
      pm(1)=dncms
      if(k-l)25,5,10
   25 pm(1)=poms
      go to 10
    5 pm(1)=pnms
   10 return
      end subroutine


************************************************************************
*                                                                      *
      subroutine signex
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      expi=0.
    5 expa=unirn(dummy)
      expc=expa
   10 expb=unirn(dummy)
      if(expb.ge.expa)go to 20
      expa=unirn(dummy)
      if(expa.lt.expb)go to 10
      expi=expi+1.
         go to 5
   20 ex=ex+ (expi+expc)/sigq
      return
c     ex=distance in sampling routine
c     exponential random divided by sigma gi region i
      end subroutine


************************************************************************
*                                                                      *
      subroutine crjab(k1,pp)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c
      dimension pp(380)
c ----------------------------------------------------------------------
c
      call crdet(k1,pp(1),value1)
      value1=(pxyz(1)*pxyz(2)+pxyz(5)*pxyz(6)+pxyz(9)*pxyz(10))
     1 /e(1)
      value2=(value1/(p2*p2))*((e(2)/dncms)-1.0)-(1.0/dncms)
      value2=dncms*crdt(1)*dsqrt(p1oe1*p1oe1+2.0*value1
     1 *value2+p2*p2*value2*value2)/(e(2)*p1oe1*any)
      if(value2-1.0)10,10,5
    5 write(6,9)
    9 format(23h value2 .gt. 1 in crjab)
      value2=1.0
      go to 10
   10 value1=unirn(dummy)
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine spac32(i)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      ex=0.0
      i4=5
      sigq=9.99999e-1*space(i)
      if(i-31)3,5,10
    3 i2=18
      i3=3
      go to 15
    5 i2=22
      i3=5
      go to 15
   10 if(i-41)20,25,25
   20 i2=26
      i3=7
      go to 15
   25  i4=3
      if(i-42)30,35,40
   30 i2=35
      i3=13
      go to 15
   35 i2=37
      i3=17
      go to 15
   40 i2=39
      i3=21
   15 call bg6b
      call signex
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine idk
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
c       univ=pnidk(6)*pnidk(6)
      pnidk(7)=(pnidk(1)*pnidk(1)+univ-sqnm)/(2.0*pnidk(1))
      pnidk(8)=dsqrt(pnidk(7)*pnidk(7)-univ)
      call polt(iands   ,pnidk(20),pnidk(21))
      call azit(iands   ,pnidk(22),pnidk(23))
      pnidk(9)=pnidk(22)*pnidk(21)*pnidk(8)
      pnidk(10)=pnidk(21)*pnidk(23)*pnidk(8)
      pnidk(11)=pnidk(20)*pnidk(8)
      univ=pnidk(9)*pnidk(2)+pnidk(10)*pnidk(3)+pnidk(11)*pnidk(4)
      pnidk(12)=(pnidk(7)*pnidk(5)+univ)/pnidk(1)
      pnidk(13)=pnidk(5)-pnidk(12)
      univ=(((pnidk(5)/pnidk(1))-1.0)*univ)/(pnidk(2)*pnidk(2)+
     1pnidk(3)*pnidk(3)+pnidk(4)*pnidk(4))
      unive=pnidk(7)/pnidk(1)
      do 5 i=2,4
      pnidk(i+12)=pnidk(i)*(univ+unive)  +pnidk(i+7)
    5 pnidk(i+15)=pnidk(i)-pnidk(i+12)
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine polt(krn,cs,si)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c ----------------------------------------------------------------------
c
      cs=2.0*(unirn(dummy)-0.5)
      si= dsqrt(1.0-(cs*cs))
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine collm(m)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c
      data eps9/1.0e-9/
c ----------------------------------------------------------------------
c
      univ=e(2)+col(6)-col(11)
      unive=e(1)+col(11)
      univer=col(1)+col(6)
      k=16
      do5i=1,9,4
      col(k)=(pxyz(i)*univ-pxyz(i+1)*unive)/univer
      col(k+3)=(pxyz(i)+pxyz(i+1))/col(1)
    5 k=k+1
      col(22)=(pxyz(10)*pxyz(5)-pxyz(9)*pxyz(6))/col(1)
      col(23)=(pxyz(2)*pxyz(9)-pxyz(10)*pxyz(1))/col(1)
      col(24)=(pxyz(6)*pxyz(1)-pxyz(5)*pxyz(2))/col(1)
      a=snt/col(14)
      b=a*col(10)
      univ=col(10)*(cst-a*sopc*col(13))/col(12)
      unive=b*sops/col(12)
      univer=(sopc*b)+((e(3)+col(9))/(col(7)+1.0))
      k=19
      do10i=3,11,4
      pxyz(i)=col(k)*univer+col(k+3)*unive+col(k-3)*univ
   10 k=k+1
      if(m)11,15,11
   11 if(pt(15))15,25,15
   15 do20i=1,9,4
   20 pxyz(i+3)=pxyz(i)+pxyz(i+1)-pxyz(i+2)
      if(m)25,60,25
   25 if(pt(3))30,45,30
   30 pt(4)=pm(3)
      i1=3
   35 i2=-1
      call pstor
   40 if(i1-3)45,45,55
   45 if(pt(15))50,55,50
   50 pt(16)=dncms
      i1=4
      go to 35
   55 pt(27)=0.0
      pt(39)=0.0
   60 return
      end subroutine


************************************************************************
*                                                                      *
      subroutine  bovera(v,ve,ver)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      del=(ve/(v*rcpmv+ve))*(ve/(v*rcpmv+ve))
      ver=1.-.5*del
      if(del.gt..002)ver=dsqrt(1.-del)
      ver=dsqrt(((1.109*ver)+.691)*ver+.108)/ver
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine undis
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
c ----------------------------------------------------------------------
c
      if(begru-1.0)5,20,20
c     if negative, beginning

    5 begru=1.0

   10 ran1=2.*(unirn(dummy)-0.5)
      ran2=2.*(unirn(dummy)-0.5)
      if( ran1**2 + ran2**2 .ge. 1.00 ) goto 10

      xi(1) = ran1 * out(4)
      xi(2) = ran2 * out(4)
      xi(3) = -out(4)

      dco(1) = 0.0
      dco(2) = 0.0
      dco(3) = 1.0

      med = 4

      curr(4)  = xi(1)
      curr(5)  = xi(2)
      curr(6)  = xi(3)
      curr(7)  = 0.0
      curr(8)  = 0.0
      curr(9)  = 1.0
      curr(10) = med

c     x, y and z coordinates, also alpha, beta and gamma
c     direction cosines.  med=4, (no. of geom)

   15 return

   20 begru=begru+1.0
      if(casesn -begru)30,25,25
   25 frand=unirn(dummy)
      go to 10
   30 begru=0.0

      go to 15

c     final random in erand.  run completed
      end subroutine


************************************************************************
*                                                                      *
      subroutine geo
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      i1=0
      t1=out(2)*out(2)
      t2=out(3)*out(3)
      t3=out(4)*out(4)
      t4=2.0*t3
   95 t5=xi(1)*xi(1)+xi(2)*xi(2)+xi(3)*xi(3)
      go to(100,125,155,185),med
  100 t6=t5-t1
      if(t6)210,210,105
  105 temp=t1
  110 if((t6/temp)-5.0e-6)115,115,205
  115 do120i=1,3
      xi(i)=xi(i)*9.99995d-1
  120 curr(i+3)=xi(i)
       go to 95
  125 t6=t5-t1
      if(t6)130,130,150
  130 temp=t1
  135 if(5.0e-6+(t6/temp))205,140,140
  140 do145i=1,3
      xi(i)=xi(i)*10.00005d-1
  145 curr(i+3)=xi(i)
      go to 95
  150 t6=t5-t2
      temp=t2
      if(t6)210,210,110
  155 t6=t5-t2
      if(t6)160,160,165
  160 temp=t2
      go to 135
  165 t6=t5-t3
      if(t6)210,210,170
  170 temp=t3
      go to 110
c 175 if(xi(2))205,180,205
  180 if(curr(5))205,185,205
  185 t6=t5-t3
      if(t6)190,190,195
  190 temp=t3
      go to 135
  195 t6=t5-t4
      if(t6)210,210,200
  200 temp=t4
      go to 110
  205 i1=-1
      go to 10
  210 t4=xi(1)*dco(1)+xi(2)*dco(2)+xi(3)*dco(3)
      t6=t4*t4
      t6=t5-t6
      if(t3-t6)205,1000,1000
 1000 t3=dsqrt(t3-t6)
      temp=t2-t6
      t2=dsqrt( abs(temp))
      tempo=t1-t6
      t1=dsqrt( abs(tempo))
      do1005i=1,6
 1005 d(i)=0.0
      go to (5,15,39,65),med
    5 if(temp)205,6,6
    6 d(4)=t1-t4
      d(5)=t2-t1
      d(6)=t3-t2
   10 return
   15 if(temp)205,16,16
   16 d(6)=t3-t2
   19 if(t4)25,20,20
   20 d(3)=t2-t4
      go to 10
   25 if(tempo)20,30,30
   30 d(3)=-(t4+t1)
      d(4)=t1+t1
      d(5)=t2-t1
      go to 10
   39 if(t4)45,40,40
   40 d(2)=t3-t4
      go to 10
   45 if(temp)40,50,50
   50 d(2)=-(t4+t2)
      d(6)=t3-t2
      if(tempo)55,60,60
   55 d(3)=t2+t2
      go to 10
   60 d(3)=t2-t1
      d(5)=d(3)
      d(4)=t1+t1
      go to 10
   65 d(1)=-(t4+t3)
   69 if(temp)70,75,75
   70 d(2)=t3+t3
      go to 10
   75 d(2)=t3-t2
      d(6)=d(2)
      if(tempo)85,80,80
   80 d(3)=t2-t1
      d(5)=d(3)
      d(4)=t1+t1
      go to 10
   85 d(3)=t2+t2
      go to 10
      end subroutine


************************************************************************
*                                                                      *
      subroutine partin
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      data eps9/1.0e-9/
c ----------------------------------------------------------------------
c
      if(d(4))230,235,230
  230 ipec(10)=ipec(10)+1
      go to 250
  235 if(d(3))240,245,240
  240 ipec(6)=ipec(6)+1
      go to 250
  245 ipec(1)=ipec(1)+1
  250 do255i=1,3
  255 isw(i)=0
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine crdet(nodata,data,ener)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c
      dimension data(380)
c ----------------------------------------------------------------------
c
      ie= abs(ener/20.0)
      univ=(ener- dble(ie)*20.0)/20.0
   85 k=(nodata*ie)+1
   99 if(inpt)98,100,115
c%%98 write(6,97)
c%%97 format(' inpt .lt. 0 in crdet')
   98 write(6,97) inpt
   97 format(/' *** error message from s.crdet ***'
     &/' variable "inpt" was less than zero.'
     &/' inpt =',i5)
      call parastop( 838 )
  100 n=nodata
  150 l=k+nodata
      do110i=1,n
      crdt(i)=(data(l)-data(k))*univ+data(k)
      k=k+1
  110 l=l+1
      inpt=0
      return
  115 k=inpt-1+k
      n=2
      go to 150
      end subroutine


************************************************************************
*                                                                      *
      subroutine bb
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      i=i2
      if(i1.ne.0)i=i+3
      if(knot.le.6.or.knot.gt.12)clcfe=cfepn(i)
   25  e(1)=wkrpn(i)*rcpmv+pm(1)
      if(in.ne.0)go to 65
   50 call p1cli
      go to 60
   65 call p1clc
   60 return
      end subroutine


************************************************************************
*                                                                      *
      subroutine spisom
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
    5 call isom
      if(i1)10,10,15
   10 space(176)=pxyz(2)
      space(177)=pxyz(6)
      space(178)=pxyz(10)
      if(i1)25,20,20
   15 pxyz(2)=pxyz(2)+space(176)
      pxyz(6)=pxyz(6)+space(177)
      pxyz(10)=pxyz(10)+space(178)
      e(2)=(pxyz(10)*pxyz(10)+pxyz(6)*pxyz(6)+pxyz(2)*pxyz(2))/
     11.9032e14
   20 return
   25 i1=1
      go to 5
      end subroutine


************************************************************************
*                                                                      *
      subroutine p1clc
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      p1oe1=e(1)-.5*pm(1)*pm(1)/e(1)
      if(pm(1).gt..05*e(1))p1oe1=dsqrt(e(1)*e(1)-pm(1)*pm(1))
      pxyz(1)=p1oe1*curr(7)
      pxyz(5)=p1oe1*curr(8)
      pxyz(9)=p1oe1*curr(9)
      p1oe1=p1oe1/e(1)
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine p1cli
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      pxyz(1)=0.0
      pxyz(5)=0.0
      pxyz(9)=e(1)-.5*pm(1)*pm(1)/e(1)
      if(pm(1).gt..05*e(1))pxyz(9)=dsqrt(e(1)*e(1)-pm(1)*pm(1))
      p1oe1=pxyz(9)/e(1)
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine mud(krn,sine,inp)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c ----------------------------------------------------------------------
c
      sine=unirn(dummy)*20.
      inp = idint(     sine  + 1.)
      sine= dble(inp)-sine
c     sine=(.05n-r)/.05=n-r/.05   n=inpt   r/.05=(n-1)+x
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine azit(krn,sine,cosine)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c ----------------------------------------------------------------------
c
    5 r1=unirn(dummy)
       r1sq=r1*r1
      r2=2.*(unirn(dummy)-0.5)
      r2sq=r2*r2
      sum=.5*(r1sq+r2sq)
      if(sum.gt..5)go to 5
      cosine=(sum-r1sq)/sum
      sine=(r1*r2)/sum
   10 return
      end subroutine


************************************************************************
*                                                                      *
        subroutine  coll(m)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
 100  if(m)1,2,2
    1 a=pm(4)*pm(4)
      go to 3
    2 a=sqnm
    3  col(15)=0.0
      med=clsm
      eco (1)=cfepn(med)
      eco (2)=cfepn(med+3)
c     proton(neutron) energy cut-off
      col(1)=e(1)+e(2)
c     total energy particles 1 and 2
      do5i=1,3
    5 col(i+1)=pm(i)*pm(i)
c     mass particle i sqd.
      col(5)=col(3)+col(2)+2.0*(e(1)*e(2)-(pxyz(1)*pxyz(2)+pxyz(5)*
     1pxyz(6)+pxyz(9)*pxyz(10)))
      col(6)=dsqrt(col(5))
      col(7)=col(6)/col(1)
c     gam
      col(8)=2.0*col(6)
      col(9)=(col(4)+col(5)-a)/col(8)
      com2=col(9)*col(9)
   50 if(col(4)-2.9882156d27)51,51,52
c     gt,pm(3)=isobar--lte,test for roundoff range,(min)sqd+or-5d23
   51 if(col(4)-2.9872156d27)53,552,552
c     lt,pion or nucleon mass=pm(3)
  552 col(4)=2.9877156d27
      pm(3) = 5.466005d13
   52 if(com2-col(4))54,56,56
   53 if(col(4) - sqnm) 52,52,11
c     lte,have nucleon or pion--gt,go to error
   54 if(com2 - 9.9d-1 * col(4)) 56,55,55
   55 com2 = col(4)
      col(9) = pm(3)
   56 col(10)=dsqrt(com2-col(4))
c     p3 prime
      col(11)=(col(5)+col(2)-col(3))/col(8)
c     e1 prime
      col(12)=dsqrt(col(11)*col(11)-col(2))
c     p1 prime
      col(13)=(col(7)*e(1)-col(11))/col(12)
c     beta
      com=1.0-(col(13)*col(13)+col(7)*col(7))
      if(com-5.0d-6)10,25,25
   10 if(com+5.0d-6)11,20,20
   11 col(15)=1.0
c     error
   12 return
   20 col(14)=2.236067977d-3
      go to 30
   25 col(14)=dsqrt(com)
c     alpha
   30 e(3)=(col(9)+col(10)*(col(13)*cst+col(14)*sopc*snt))/col(7)
      e(4)=col(1)-e(3)
      go to 12
      end subroutine


************************************************************************
*                                                                      *
      subroutine cole4
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      k=clsm
      com=(e(4)-dncms)/rcpmv
      if(curr(1)-3.0)12,35,35
   35 if(abx)40,45,40
   40 com=com+absec
      go to 12
   45 if(it-5)50,65,50
   50 if(it-24)55,65,55
   55 if(it-6)60,70,60
   60 if(it-26)12,70,12
   65 univ=0.0
      go to 75
   70 univ=1.0
   75 unive=space(k+3)-space(k+9)
       if(univ)85,80,85
   80 com=com+unive
      go to 12
   85 com=com-unive
   12 return
      end subroutine


************************************************************************
*                                                                      *
      subroutine pinst
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      data eps9/1.0e-9/
c ----------------------------------------------------------------------
c
      i1=0
      med=clsm
      if(inc)5,70,5
    5 inc=0
      if(med-1)10,25,20
   10 i1=1
   15 return
   20 if(med-3)30,45,10
   25 ipec(12)=ipec(12)+1
      go to 15
   30 if(d(4))40,35,40
   35 ipec(8)=ipec(8)+1
      go to 15
   40  ipec(9)=ipec(9)+1
      go to 15
   45 if(d(3))50,60,50
   50 if(d(4))55,65,55
   55 ipec(5)=ipec(5)+1
      go to 15
   60 ipec(3)=ipec(3)+1
      go to 15
   65 ipec(4)=ipec(4)+1
      go to 15
   70 k=curr(11)
      k=3*(med-1)+k
      cc(k)=cc(k)+1.0
      go to 15
      end subroutine


************************************************************************
*                                                                      *
      subroutine punp
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      data eps9/1.0e-9/
c ----------------------------------------------------------------------
c
      if(pgvc(1))5,20,14
    5 i1=-1
   10 return
   14 pgvc(1)=pgvc(1)-11.0
      k=pgvc(1)+2.005
      do 15 i=1,11
      curr(i)=pgvc(k)
      pgvc(k)=0.0
   15 k=k+1
   16 i1=1
      go to 10
   20 if(plvc(1))5,55,25
   25 univ=0.0
      l=plvc(1)
      k=-10
      do 45 i=1,l
   30 k=k+12
      if(plvc(k))5,30,35
   35 if(plvc(k)-univ)45,40,40
   40 univ=plvc(k)
      m=k
   45 continue
      plvc(m)=0.0
      do 50 i=1,11
      m=m+1
      curr(i)=plvc(m)
   50 plvc(m)=0.0
      plvc(1)=plvc(1)-1.0
      go to 16
   55 i1=0
      go to 10
      end subroutine


************************************************************************
*                                                                      *
      subroutine stpr
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      data eps9/1.0e-9/
c ----------------------------------------------------------------------
c
      i1=0
      med=clsm
      do 5 i=3,39,12
      k=i
      if(pt(i))15,5,15
   15 if(abs(pt(k-1)-2.0).lt.eps9) go to 25
      if(pt(k-1)-2.0)20,25,20
   20 pt(k-2)=pt(k)-space(med+9)
      go to 40
   25 pt(k-2)=pt(k)-space(med+3)
   30 if(pt(k-2)-500.0)35,35,50
   35 call stpl(pt(k-2))
      if(i1)55,5,55
   40 if(abs(pt(k-1)-1.0).lt.eps9) go to 30
      if(pt(k-1)-1.0)45,30,45
   45 pt(k-2)=(dncms*pt(k-2))/pt(k+1)
      go to 30
   50 call stph(pt(k-1))
      if(i1)55,5,55
    5 continue
   55 return
      end subroutine


************************************************************************
*                                                                      *
      subroutine snn
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      if(2500.0-rlke)5,15,15
    5 i1=-1
   10 return
   15 if(1000.0-rlke)20,25,25
   20 i1=0
      go to 10
   25 i1=1
      go to 10
      end subroutine


************************************************************************
*                                                                      *
      subroutine dcintp(w)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c
      dimension w(50),z(18)
      data eps9/1.0e-9/
c ----------------------------------------------------------------------
c
      i2=0
      k=1
      do5i=1,5
      if(abs((w(k)-rlke)/rlke).lt.eps9) go to 4
      if(w(k)-rlke)5,4,10
    4 if(w(k))5,10,10
    5 k=k+9
    6 i2=1
      go to 46
   10 do15l=1,9
      z(l+9)=w(k)
      z(l)=w(k-9)
   15  k=k+1
      i=2
   20 unive=((z(i+ 9)-z(i)) * (rlke-z(1))/(z(10)-z(1)))+z(i)
      univ=unirn(dummy)
      if(univ-unive)25,90,90
   25 go to (6,30,35,40,35,35,35,35,35),i
   30 i=3
      go to 20
   35 univ=unirn(dummy)
      go to (45,50,55,60,65,70,75,80,85),i
   40 i=6
      go to 20
   45 cst=.25*univ-.25
   46 return
   50 cst=.9375+.0625*univ
      go to 46
   55 cst=.5*univ-1.0
      go to 46
   60 cst=.25+.25*univ
      go to 46
   65 cst=.25*univ-.5
      go to 46
   70 cst=.25*univ
      go to 46
   75 cst=.5+.25*univ
      go to 46
   80 cst=.75+.125*univ
      go to 46
   85 cst=.875+.0625*univ
      go to 46
   90 go to (6,95,95,105,110,115,120,120,130),i
   95 i=i+2
      go to 20
  105 i=7
      go to 20
  110 i=1
      go to 35
  115 i=4
      go to 35
  120 i=i+1
      go to 20
  130 i=2
      go to 35
      end subroutine


************************************************************************
*                                                                      *
       subroutine  big7t(krn,c,great,ix)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c ----------------------------------------------------------------------
c
      great=c
      i=ix+1
      do 5 k=1,i
      c=unirn(dummy)
      great=dmax1(c,great)
    5 continue
c     great is the largest of i random nos.
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine gene(z)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      dimension z(3800)
c ----------------------------------------------------------------------
c
      cd=com*100.0
      i=cd+1.0
      com=z(i)+(z(i+1)-z(i))*(cd- dble(i-1))
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine alp19
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
      univ=unirn(dummy)
      pt(2)=1.0
      pt(26)=1.0
      pt(14)=3.0
      pt(16)=poms
      if(isw(12))5,5,70
    5 if(univ-.25)10,10,55
   10 if(isw(4))20,15,20
   15 pt(2)=2.0
   20 univ  =unirn(dummy)
      if(univ-6.6666667e-1)25,25,40
   25 pt(14)=4.0
   26 if(isw(4))35,30,35
   30 pt(26)=2.0
   35  go to 1000
   40 pt(16)=pnms
      if(isw(4))30,50,30
   50 pt(14)=5.0
      go to 35
   55 if(isw(4))65,60,65
   60 pt(26)=2.0
      go to 40
   65 pt(2)=2.0
      pt(16)=pnms
      go to 35
   70 if(univ-.5)75,75,95
   75 if(isw(4))80,85,80
   80 pt(2)=2.0
   85 univ  =unirn(dummy)
      if(univ-3.3333333e-1)40,40,90
   90 pt(14)=4.0
      go to 26
   95 if(isw(4))105,100,105
  100 pt(2)=2.0
  105 univ  =unirn(dummy)
      if(univ-6.6666667e-1)110,110,115
  110 pt(14)=4.0
      if(isw(4))30,35,30
  115 pt(16)=pnms
      if(isw(4))50,30,50
 1000 return
      end subroutine


************************************************************************
*                                                                      *
      subroutine alp28
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      r=unirn(dummy)
      if(isw(13))95,5,95
    5 if(r-6.0e-1)10,10,50
   10 pt(4)=pnms
      r=unirn(dummy)
      if(isw(4))15,35,15
   15 if(r-3.3333333e-1)20,20,30
   20 pt(26)=5.0
   21 pt(28)=pnms
   25 return
   30 pt(26)=4.0
   31 pt(38)=2.0
      go to 25
   35 pt(2)=5.0
       pt(14)=2.0
      if(r-3.3333333e-1)40,40,45
   40 pt(28)=pnms
      go to 31
   45 pt(26)=4.0
      go to 25
   50 r=unirn(dummy)
      if(isw(4))55,75,55
   55 if(r-6.6666667e-1)60,60,70
   60 pt(2)=4.0
   65 r=unirn(dummy)
      if(r-6.6666667e-1)45,45,40
   70 pt(14)=2.0
   71 pt(4)=pnms
      go to 65
   75 if(r-6.6666667e-1)80,80,90
   80 pt(2)=4.0
   81 pt(14)=2.0
   85 r=unirn(dummy)
      if(r-6.6666667e-1)30,30,20
   90 pt(2)=5.0
      pt(4)=pnms
      go to 85
   95 if(r-value1)100,100,115
  100 pt(4)=pnms
      if(isw(4))110,105,110
  105 pt(2)=5.0
      pt(14)=2.0
      go to 21
  110 pt(38)=2.0
      go to 20
  115 r=unirn(dummy)
      if(isw(4))120,135,120
  120 if(r-3.3333333e-1)125,125,130
  125 pt(4)=pnms
      go to 81
  130 pt(2)=4.0
      go to 85
  135 if(r-3.3333333e-1)140,140,145
  140 pt(2)=5.0
      go to 71
  145 pt(14)=2.0
      go to 60
      end subroutine


************************************************************************
*                                                                      *
      subroutine dfmax
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      i=i2
      if(curr(1)-2.0)10,5,5
    5 i=i+3
   10 wk=wkrpn(i)
   15 call bovera(wk,dncms,univ)
      if(wk-560.0)20,20,30
   20 fmax(1)=25.3e-27*univ
      fmax(2)=40.0e-27*univ
      fmax(3)=21.0e-27*univ
      fmax(4)=10.5e-27*univ
      fmax(5)=0.0
      fmax(6)=0.0
   25 return
   30 if(wk-800.0)35,45,45
   35  fmax(2)=37.0e-27*univ
      fmax(5)=4.1e-27*univ
      fmax(6)=3.6e-27*univ
   40 fmax(1)=25.3e-27*univ
      fmax(3)=21.2e-27*univ
      fmax(4)=10.6e-27*univ
      go to 25
   45 if(wk-1680.0)50,55,55
   50 fmax(2)=33.0e-27*univ
      fmax(5)=12.3e-27*univ
      fmax(6)=11.3e-27*univ
      go to 40
   55 fmax(1)=19.4e-27*univ
      fmax(2)=26.0e-27*univ
      fmax(3)=21.2e-27*univ
      fmax(4)=10.6e-27*univ
      fmax(5)=15.0e-27*univ
      fmax(6)=12.3e-27*univ
      go to 25
      end subroutine


************************************************************************
*                                                                      *
      subroutine store
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      do 2 i=1,21
    2 ce(i)=0.
      n=i1
      if(i3)85,5,85
    5 k=i2+6
      l=i2
      go to (10,10,10,15,10),i5
   15 if(i1-5)20,100,20
   20 mm=-6
   19 l=k
      n=n-1
   10 do25i=1,n,2
      id=i
      ce(i+1)=fmax(i)*1.0e30*space(k)
       ce(i+2)=fmax(i+1)*1.0e30*space(l)
      go to (25,25,11,11,11),i5
   11 m=k
      k=l
      go to (25,25,30,40,35),i5
   35 if(id-2)25,25,36
   40 k=m+mm
      mm=-mm
      l=k
      go to 25
   36 k=i2
   30 l=m
   25 continue
      go to (26,26,26,95,26),i5
csib *****   changed by sib ('94.06.17)   *****
   26 sigq=0.0
csib6 sign=0.0
      do45i=2,8
   45 sigq=sigq+ce(i)
csib5 sign=sign+ce(i)
      go to (50,50,60,110,60),i5
   50 sigq=sigq*9.99999e-1
csib0 sign=sign*9.99999e-1
      return
   60 if(i1-4)65,65,50
   65 if(i3)70,70,75
   70 space(i2+87)=sigq
csib0 space(i2+87)=sign
      go to 80
   75 space(i2+106)=sigq
csib5 space(i2+106)=sign
   80 fmax(5)=0.0
      fmax(6)=0.0
      go to 50
   85 k=i2
      l=i2+6
      go to (10,10,10,90,10),i5
   90 mm=6
      go to 19
   95 ce(8)=fmax(7)*1.0e30*space(l)
      go to 26
  100 sub1=fmax(1)*1.0e30
      sub2=fmax(2)*1.0e30
      do105i=2,n,2
      ce(i)=sub1*space(k)
      ce(i+1)=sub2*space(k)
      k=l
  105 l=k+6
      ce(n+1)=space(k)*1.0e30*fmax(n)
      go to 26
  110 if(i1-5)50,115,50
  115 space(i2+68)=sigq
csib5 space(i2+68)=sign
csib ******************************************
      fmax(6)=0.0
      fmax(7)=0.0
      go to 50
      end subroutine


************************************************************************
*                                                                      *
      subroutine pfmax
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      i=i2
      if(curr(1)-2.0)10,5,5
    5 i=i+3
   10 wk=wkrpn(i)
      if(wk-160.0)20,20,30
   20 fmax(1)=0.183e-24
      fmax(2)=1.06e-24
      fmax(3)=0.0
      fmax(4)=0.0
   25 return
   30 call bovera(wk,dncms,univ)
      if(wk-400.0)35,60,60
   35 if(wk-300.0)40,55,55
   40 if(wk-200.0)45,50,50
   45  fmax(3)=2.3e-27*univ
      fmax(4)=1.0e-27*univ
      fmax(1)=0.107e-24
      fmax(2)=0.395e-24
      go to 25
   50 fmax(1)=0.09e-24
      fmax(2)=0.26e-24
      fmax(3)=14.8e-27*univ
      fmax(4)=7.3e-27*univ
      go to 25
   55 fmax(1)=30.0e-27*univ
      fmax(2)=0.073e-24
      fmax(3)=19.1e-27*univ
      fmax(4)=9.6e-27*univ
      go to 25
   60 fmax(1)=25.5e-27*univ
      fmax(2)=52.0e-27
      fmax(3)=21.0e-27*univ
      fmax(4)=10.5e-27*univ
      go to 25
      end subroutine


************************************************************************
*                                                                      *
      subroutine ccpes
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c
      data eps9/1.0d-9/
c ----------------------------------------------------------------------
c
      i1=0
      if(abs(curr(1)-2.0d+0).lt.eps9) go to 10
      if(curr(1)-2.0d+0)5,10,5
    5 k=curr(10)+9.05
      go to 15
   10 k=curr(10)+3.05
   15 if(esps(1)-60.0)30,20,20
   20 i1=1
   25 return
   30 l=esps(1)*8.0d+0+2.05d+0
      esps(l)=curr(1)
      esps(l+1)=curr(2)-space(k)
      m=13.05d+0-curr(11)
      cc(m)=cc(m)+1.0d+0
      m=4
      l=l+2
      n=l+2
      do 35 i=l,n
      esps(i)=curr(m+3)
      esps(i+3)=curr(m)
   35 m=m+1
      esps(1)=esps(1)+1.0d+0
      go to 25
      end subroutine


************************************************************************
*                                                                      *
      subroutine spcn
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      fmax(4)=0.0
      fmax(5)=0.0
      i6=curr(1)-1.95
      wk=wkrpn(i2)
      go to (5,10,5),i6
    5 univ=pnms
      go to 20
   10 univ=poms
   20 call bovera(wk,univ,unive)
      value1=0.0
      if(wk-220.0)25,25,55
   25 go to (30,45,30),i6
   30 fmax(5)=1.0e-27*unive
      fmax(6)=3.5e-27*unive
       fmax(1)=.205e-24*unive
      fmax(3)=.0451e-24*unive
   35 fmax(2)=.0219e-24*unive
   40 return
   45 fmax(6)=2.0e-27*unive
      fmax(1)=83.0e-27*unive
      fmax(2)=.0451e-24*unive
   50 fmax(3)=fmax(1)
      fmax(4)=fmax(2)
      fmax(7)=fmax(6)
      go to 40
   55 if(wk-400.0)60,60,75
   60 go to (65,70,65),i6
   65 fmax(5)=4.5e-27*unive
      fmax(6)=20.0e-27*unive
      fmax(1)=.205e-24
      fmax(3)=.0451e-24
      go to 35
   70 fmax(6)=11.0e-27*unive
      fmax(1)=83.0e-27
      fmax(2)=.0451e-24
      go to 50
   75 if(wk-500.0)80,80,100
   80 go to (85,90,85),i6
   85 fmax(1)=.113e-24
      fmax(2)=16.6e-27*unive
      fmax(3)=27.0e-27
      fmax(5)=6.6e-27*unive
      fmax(6)=20.0e-27*unive
      go to 40
   90 fmax(1)=48.0e-27
      fmax(2)=27.0e-27
      fmax(6)=12.0e-27*unive
      go to 50
  100 value1=1.0
      if(wk-600.0)105,105,120
  105 go to (110,115,110),i6
  110 fmax(1)=60.0e-27
      fmax(2)=26.0e-27*unive
      fmax(3)=17.0e-27*unive
      fmax(5)=9.0e-27*unive
      fmax(6)=23.0e-27*unive
      go to 40
  115 fmax(1)=27.0e-27*unive
      fmax(2)=17.0e-27*unive
      fmax(6)=15.5e-27*unive
      go to 50
  120 if(wk-800.0)125,125,140
  125 go to (130,135,130),i6
  130 fmax(1)=33.0e-27
      fmax(2)=28.0e-27*unive
      fmax(3)=12.0e-27*unive
      fmax(5)=21.0e-27*unive
      fmax(6)=23.4e-27*unive
      go to 40
  135 fmax(1)=15.0e-27*unive
      fmax(2)=12.0e-27*unive
      fmax(6)=19.0e-27*unive
      go to 50
  140 go to (145,150,145),i6
  145 fmax(1)=16.5e-27*unive
      fmax(2)=28.0e-27*unive
      fmax(3)=9.5e-27*unive
      fmax(5)=26.0e-27*unive
      fmax(6)=25.0e-27*unive
      go to 40
  150 fmax(1)=14.5e-27*unive
      fmax(2)=9.5e-27*unive
      fmax(6)=23.0e-27*unive
      go to 50
      end subroutine


************************************************************************
*                                                                      *
      subroutine isom
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
cKN   data ccc/.5/
c ----------------------------------------------------------------------
c
      call polt(iands,polc,pols)
      call azit(iands,sopc,sops)

      m=clsm+0.05
      if (strkp+2.0) 15,15,10
   10 fermn=fmpn(m)
      go to 60
   15 fermn=fmpn(m+3)

cKN
      ccc = 0.0
   60 call big7(iands,ccc,p2,2)
cKN
c  60 call big7(iands,ccc,p2,1)

      p2=fermn*p2
      a=p2*pols
      pxyz(2)=a*sopc
      pxyz(6)=a*sops
      pxyz(10)=p2*polc
      e(2)=dsqrt(p2*p2+sqnm)
      rlke=(((e(1)*e(2)-pxyz(1)*pxyz(2)-pxyz(5)*pxyz(6)-pxyz(9)*
     &      pxyz(10))/dncms)-pm(1))/rcpmv
c
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine pstor
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c ----------------------------------------------------------------------
c
      l=(i1*12)-28
      if(i2)3,14,23
    3 jj=0
      if(pm(3)-dncms)4,4,300
  300 i1=i1+1
      jj=1
    4 univ=dsqrt(pxyz(i1)*pxyz(i1)+pxyz(i1+4)*pxyz(i1+4)+pxyz(i1+8)
     1*pxyz(i1+8))
      k=i1+8
       do5i=i1,k,4
      pt(l)=pxyz(i)/univ
    5 l=l+1
      i1=i1-jj
    6 pt(l)=clsm
       pt(l+1)=curr(11)
      pt(l-6)=c(1)
      pt(l-5)=c(2)
      pt(l-4)=c(3)
      return
   14 k=14
      go to 25
   23 if(i2-2)24,35,35
   24 k=17
   25 univ=dsqrt(pnidk(k)*pnidk(k)+pnidk(k+1)*pnidk(k+1)+pnidk
     1(k+2)*pnidk(k+2))
      j=k+2
      do30i=k,j
      pt(l)=pnidk(i)/univ
   30 l=l+1
      go to 6
   35 univ=dsqrt(pt(l-3)*pt(l-3)+pt(l-2)*pt(l-2)+pt(l-1)*pt(l-1))
      k=l-1
      m=l-3
      do40i=m,k
      pt(l)=pt(i)/univ
   40 l=l+1
      go to 6
      end subroutine


************************************************************************
*                                                                      *
      subroutine dcpr(w)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c
      dimension w(64),z(16)
      data eps9/1.0e-9/
c ----------------------------------------------------------------------
c
      k=1
      do5i=1,8
      if(abs((w(k)-rlke)/rlke).lt.eps9) go to 500
      if(w(k)-rlke)5,500,10
  500 if(w(k))5,10,10
    5 k=k+8
    6 i1=-1
    7 return
   10 do15l=1,8
      z(l+8)=w(k)
      z(l)=w(k-8)
   15 k=k+1
       p=0.0
      univ=unirn(dummy)
      unive=(rlke-z(1))/(z(9)-z(1))
      do20i=2,7
      p=(z(i+8)-z(i))*unive+z(i)+p
      if(p-univ)20,20,25
   20 continue
      i=7
      if(p-9.99e-1)6,25,25
   25 if(.1001e+1-p)6,26,26
   26 i1=i
      univer=(z(16)-z(8))*unive+z(8)
      univ=unirn(dummy)
      go to (6,30,35,40,45,50,55),i1
   30 cst=1.0-univ*.5e-1
      go to 60
   35 cst=9.5e-1-univ*.5e-1
      go to 60
   40 cst=9.0e-1-univ*1.0e-1
      go to 60
   45 cst=8.0e-1-univ*2.0e-1
      go to 60
   50 cst=6.0e-1-univ*2.0e-1
      go to 60
   55 cst=4.0e-1-univ*4.0e-1
   60 univ=unirn(dummy)
      if(univ-univer)65,7,7
   65 cst=-cst
      go to 7
      end subroutine


************************************************************************
*                                                                      *
      subroutine stpl(w)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      dimension w(12)
      data eps9/1.0e-9/
c ----------------------------------------------------------------------
c
      do5k=2,950,12
      if(plvc(k))5,10,5
    5 continue
   20 i1=1
   25 return
   10 i=k
      j=i+11
      l=1
      do15k=i,j
      plvc(k)=w(l)
   15 l=l+1
      plvc(1)=plvc(1)+1.0
       go to 25
      end subroutine


************************************************************************
*                                                                      *
      subroutine stph(w)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      dimension w(11)
c ----------------------------------------------------------------------
c
      if(pgvc(1)-440.0)5,20,20
    5 i=pgvc(1)+2.005
      j=i+10
      l=1
      do15k=i,j
      pgvc(k)=w(l)
   15 l=l+1
      pgvc(1)=pgvc(1)+11.0
   25 return
   20 i1=1
      go to 25
      end subroutine


************************************************************************
*                                                                      *
       subroutine  big7(krn,c,great,ix)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c ----------------------------------------------------------------------
c
      great=c
      i=ix+1
      do 5 k=1,i
      c=unirn(dummy)
      great=dmax1(c,great)
    5 continue
c     great is the largest of i random nos.
      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine angid
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c
      data eps9/1.0e-9/
c ----------------------------------------------------------------------
c
      if(abs(andit-1.0).lt.eps9) go to 10
      if(andit-1.0)5,10,20
    5 r=unirn(dummy)
      if(r-.5)10,10,20
   10 call polt(iands   ,cst,snt)
   15 call azit(iands   ,sopc,sops)
      return
   20 snt=.003162
csk --------------------------------- jun. 30 1994 -- blok 2 -- start --

      a1=0.999995
      a2=unirn(dummy)-0.5
      aa=dabs(a1)
      if(a2.ge.0.) signk0=aa
      if(a2.lt.0.) signk0=-aa
      cst=signk0

csk --------------------------------- jun. 30 1994 -- blok 2 --  end  --
      go to 15
      end subroutine


************************************************************************
*                                                                      *
       subroutine alpha
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include 'bert.inc'
c
      data eps9/1.e-9/
c ----------------------------------------------------------------------
c
      univ=unirn(dummy)
      if(abs(value3).lt.eps9) go to 3
      if(value3)130,3,60
    3 if(univ-value1)5,5,50
    5 if(isw(11))15,10,15
   10 pt(2)=5.0
      pt(26)=2.0
   15 pt(4)=pnms
      pm(4)=pnms
      univ=unirn(dummy)
      if(univ-value2)20,20,30
   20 pt(14)=4.0
   25 return
   30 if(isw(11))45,35,45
   35 pt(26)=1.0
   36 pt(14)=5.0
   40 pt(16)=pnms
      go to 25
   45 pt(26)=2.0
      go to 40
   50 pt(2)=4.0
      if(isw(11))40,55,40
   55 pt(14)=5.0
      go to 45
   60 if(univ-value1)70,70,80
   70 pm(4)=pnms
      if(isw(11))65,75,65
   65 pt(2)=5.0
   66 pt(16)=pnms
   67 pt(4)=pnms
      go to 25
   75 pt(14)=5.0
      pt(26)=2.0
      go to 66
   80 if(univ-value2)85,85,105
   85 pt(2)=4.0
      univ=unirn(dummy)
      if(univ-value3)100,100,90
   90 if(isw(11))20,95,20
   95 pt(26)=2.0
      go to 20
  100 if(isw(11))45,36,45
  105 pm(4)=pnms
      pt(4)=pnms
      univ=unirn(dummy)
      if(univ-6.6666667e-1)110,110,120
  110 if(isw(11))95,115,95
  115 pt(2)=5.0
      go to 20
  120 if(isw(11))36,125,36
  125 pt(26)=2.0
      go to 65
  130 if(univ-value1)135,135,150
  135 pm(4)=pnms
      pt(4)=pnms
      univ=unirn(dummy)
      if(abs(value3+1.0).lt.eps9) go to 140
      if(value3+1.0)145,140,145
  140 if(univ-3.3333333e-1)36,36,95
  145 pt(2)=5.0
      if(univ-3.3333333e-1)45,45,20
  150 if(univ-value2)155,155,170
  155 pt(2)=4.0
      univ=unirn(dummy)
      if(abs(value3+1.0).lt.eps9) go to 160
      if(value3+1.0)165,160,165
  160 if(univ-6.6666667e-1)20,20,45
  165 if(univ-6.6666667e-1)95,95,36
  170 pm(4)=pnms
      pt(4)=pnms
      if(abs(value3+1.0).lt.eps9) go to 65
      if(value3+1.0)175,65,175
  175 pt(14)=5.0
      go to 45
      end subroutine


************************************************************************
*                                                                      *
      subroutine ecpl
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      include 'bert.inc'
c
      data eps9/1.0e-9/
c ----------------------------------------------------------------------
c
      i1=0
      med=clsm
      if(pt(38))240,5,240
    5 if(dabs(curr(1)-1.0).lt.eps9) go to 215
      if(curr(1)-1.0)135,215,10
   10 if(dabs(curr(1)-3.0).lt.eps9) go to 155
      if(curr(1)-3.0)170,155,15
   15 if(dabs(curr(1)-5.0).lt.eps9) go to 20
      if(curr(1)-5.0)110,20,135
   20 if(dabs(strkp+1.0)  .lt.eps9) go to 30
      if(strkp+1.0)25,30,135
   25 if(dabs(strkp+2.0)  .lt.eps9) go to 90
      if(strkp+2.0)135,90,135
   30  if(dabs(pt(2)-2.0)  .lt.eps9) go to 31
      if(pt(2)-2.0)32,31,32
   31 pt(3)=vnvp(med)
      go to 33
   32 pt(3)=0.0
   33 pt(15)=hvp(med)
      if(dabs(pt(26)-1.0) .lt.eps9) go to 75
      if(pt(26)-1.0)135,75,35
   35 if(dabs(pt(26)-2.0) .lt.eps9) go to 40
      if(pt(26)-2.0)135,40,135
   40 pt(27)=-ppan(med)
      if(curr(1)-2.0)235,41,41
   41 if(dabs(curr(1)-4.0) .lt.eps9) go to 135
      if(curr(1)-4.0)165,135,45
   45 if(dabs(pt(2)-3.0)  .lt.eps9) go to 55
      if(pt(2)-3.0)135,55,50
   50 if(dabs(pt(2)-5.0)   .lt.eps9) go to 70
      if(pt(2)-5.0)65,70,135
   55 if(dabs(pt(14)-5.0) .lt.eps9) go to 60
      if(pt(14)-5.0)135,60,135
   60 return
   65 if(dabs(pt(14)-4.0) .lt.eps9) go to 60
      if(pt(14)-4.0)135,60,135
   70 if(dabs(pt(14)-3.0) .lt.eps9) go to 60
      if(pt(14)-3.0)135,60,135
   75 pt(27)=hvp(med)
      if(curr(1)-2.0)225,76,76
   76 if(dabs(curr(1)-4.0).lt.eps9) go to 135
      if(curr(1)-4.0)145,135,80
   80 if(dabs(pt(2)-4.0)  .lt.eps9) go to 55
      if(pt(2)-4.0)135,55,85
   85 if(dabs(pt(2)-5.0)  .lt.eps9) go to 65
      if(pt(2)-5.0)135,65,135
   90 pt(3)=-vnvp(med)
   91 pt(15)=-pmac(med)
   92 if(dabs(pt(26)-1.0) .lt.eps9) go to 100
      if(pt(26)-1.0)135,100,95
   95 if(dabs(pt(26)-2.0) .lt.eps9) go to 105
      if(pt(26)-2.0)135,105,135
  100 pt(27)=-pmac(med)
      if(dabs(curr(1)-2.0).lt.eps9) go to 55
      if(curr(1)-2.0)55,55,101
  101 if(dabs(curr(1)-4.0).lt.eps9) go to 80
      if(curr(1)-4.0)45,80,102
  102 if(dabs(pt(2)-5.0)  .lt.eps9) go to 55
      if(pt(2)-5.0)135,55,135
  105 pt(27)=hvn(med)
      if(curr(1)-2.0)65,65,106
  106 if(dabs(curr(1)-4.0).lt.eps9) go to 45
      if(curr(1)-4.0)145,45,80
  110 if(dabs(strkp+1.0)  .lt.eps9) go to 120
      if(strkp+1.0)115,120,135
  115 if(dabs(strkp+2.0)  .lt.eps9) go to 90
      if(strkp+2.0)135,90,135
  120 pt(3)=0.0
  121 pt(15)=hvp(med)
      if(dabs(pt(26)-1.0) .lt.eps9) go to 130
      if(pt(26)-1.0)135,130,125
  125 if(dabs(pt(26)-2.0) .lt.eps9) go to 140
      if(pt(26)-2.0)135,140,135
  130 pt(27)=hvp(med)
      if(dabs(curr(1)-4.0).lt.eps9) go to 45
      if(curr(1)-4.0)65,45,135
  135 i1=1
      go to 60
  140 pt(27)=-ppan(med)
      if(dabs(curr(1)-4.0).lt.eps9) go to 145
      if(curr(1)-4.0)70,145,135
  145 if(dabs(pt(2)-3.0)  .lt.eps9) go to 65
      if(pt(2)-3.0)135,65,150
  150 if(dabs(pt(2)-4.0)  .lt.eps9) go to 70
      if(pt(2)-4.0)135,70,135
  155 if(dabs(strkp+1.0)  .lt.eps9) go to 30
      if(strkp+1.0)160,30,135
  160 if(dabs(strkp+2.0)  .lt.eps9) go to 90
      if(strkp+2.0)135,90,135
  165 if(dabs(pt(2)-3.0)  .lt.eps9) go to 70
      if(pt(2)-3.0)135,70,135
  170 if(dabs(strkp+1.0)  .lt.eps9) go to 180
      if(strkp+1.0)175,180,135
  175 if(dabs(strkp+2.0)  .lt.eps9) go to 190
      if(strkp+2.0)135,190,135
  180 pt(3)=0.0
      if(dabs(pt(2)-1.0)  .lt.eps9) go to 91
      if(pt(2)-1.0)135,91,185
  185 if(dabs(pt(2)-2.0)  .lt.eps9) go to 121
      if(pt(2)-2.0)135,121,135
  190 pt(15)=-pmac(med)
      if(dabs(pt(2)-1.0)  .lt.eps9) go to 200
      if(pt(2)-1.0)135,200,195
  195 if(dabs(pt(2)-2.0)  .lt.eps9) go to 210
      if(pt(2)-2.0)135,210,135
  200 if(dabs(pt(26)-2.0) .lt.eps9) go to 205
      if(pt(26)-2.0)135,205,135
  205 pt(3)=-vnvp(med)
      pt(27)=hvn(med)
      go to 55
  210 pt(3)=0.0
      go to 92
  215 if(dabs(strkp+1.0)  .lt.eps9) go to 30
      if(strkp+1.0)220,30,135
  220 if(dabs(strkp+2.0)  .lt.eps9) go to 180
      if(strkp+2.0)135,180,135
  225 if(dabs(pt(2)-1.0)  .lt.eps9) go to 65
      if(pt(2)-1.0)135,65,230
  230 if(dabs(pt(2)-2.0)  .lt.eps9) go to 70
      if(pt(2)-2.0)135,70,135
  235 if(dabs(pt(2)-1.0)  .lt.eps9) go to 70
      if(pt(2)-1.0)135,70,135
  240 if(dabs(curr(1)-1.0).lt.eps9) go to 250
      if(curr(1)-1.0)135,250,245
  245 if(dabs(curr(1)-2.0).lt.eps9) go to 260
      if(curr(1)-2.0)135,260,135
  250 if(dabs(strkp+1.0)  .lt.eps9) go to 360
      if(strkp+1.0)255,360,135
  255 if(dabs(strkp+2.0)  .lt.eps9) go to 270
      if(strkp+2.0)135,270,135
  260 if(dabs(strkp+1.0)  .lt.eps9) go to 270
      if(strkp+1.0)265,270,135
  265 if(dabs(strkp+2.0)  .lt.eps9) go to 405
      if(strkp+2.0)135,405,135
  270 if(pt(14))135,135,275
  275 if(pt(38))135,135,280
  280 if(dabs(pt(38)-2.0) .lt.eps9) go to 314
      if(pt(38)-2.0)285,314,135
  285 if(dabs(pt(14)-2.0) .lt.eps9) go to 340
      if(pt(14)-2.0)290,340,135
  290 pt(3)=tffn(med)
      pt(15)=tffn(med)
      pt(27)=tffn(med)
      pt(39)=tffn(med)
  295 if(pt(2)-3.0)135,135,300
  300 if(dabs(pt(2)-5.0)  .lt.eps9) go to 305
      if(pt(2)-5.0)310,305,135
  305 if(dabs(pt(26)-4.0) .lt.eps9) go to 60
      if(pt(26)-4.0)135,60,135
  310 if(dabs(pt(26)-5.0) .lt.eps9) go to 60
      if(pt(26)-5.0)135,60,135
  314 if(dabs(pt(14)-2.0) .lt.eps9) go to 345
      if(pt(14)-2.0)315,345,135
  315 pt(27)=ffptfn(med)
      pt(3)=fvnp(med)
  316 pt(39)=fvnp(med)
      pt(15)=fvnp(med)
  317 if(dabs(pt(2)-3.0)  .lt.eps9) go to 310
      if(pt(2)-3.0)135,310,320
  320 if(dabs(pt(2)-5.0)  .lt.eps9) go to 325
      if(pt(2)-5.0)305,325,135
  325 if(dabs(pt(26)-3.0) .lt.eps9) go to 60
      if(pt(26)-3.0)135,60,135
  340 pt(3)=ffptfn(med)
      pt(27)=fvnp(med)
      go to 316
  345 pt(3)=tffn(med)
      pt(27)=tffn(med)
      pt(15)=tffp(med)
      pt(39)=tffp(med)
  350 if(dabs(pt(2)-3.0)  .lt.eps9) go to 305
      if(pt(2)-3.0)135,305,355
  355 if(dabs(pt(2)-4.0)  .lt.eps9) go to 325
      if(pt(2)-4.0)135,325,135
  360 if(pt(14))135,135,365
  365 if(dabs(pt(38)-1.0) .lt.eps9) go to 375
      if(pt(38)-1.0)135,375,370
  370 if(dabs(pt(38)-2.0) .lt.eps9) go to 390
      if(pt(38)-2.0)135,390,135
  375 if(dabs(pt(14)-2.0) .lt.eps9) go to 385
      if(pt(14)-2.0)380,385,135
  380 pt(3)=hvp(med)
  381 pt(15)=hvp(med)
      pt(27)=hvp(med)
      pt(39)=hvp(med)
      if(dabs(pt(14)-2.0) .lt.eps9) go to 295
      if(pt(14)-2.0)317,295,135
  385 pt(27)=hvn(med)
  386 pt(3)=-pmac(med)
  387 pt(39)=hvn(med)
      pt(15)=hvn(med)
      if(dabs(pt(38)-2.0) .lt.eps9) go to 317
      if(pt(38)-2.0)350,317,135
  390 if(dabs(pt(14)-2.0) .lt.eps9) go to 400
      if(pt(14)-2.0)395,400,135
  395 pt(3)=hvn(med)
      pt(15)=hvn(med)
      pt(27)=-pmac(med)
      pt(39)=hvn(med)
      go to 350
  400 pt(15)=-ppan(med)
      pt(39)=-ppan(med)
      pt(3)=hvp(med)
      pt(27)=hvp(med)
      if(dabs(pt(2)-3.0)  .lt.eps9) go to 325
      if(pt(2)-3.0)135,325,135
  405 if(pt(14))135,135,410
  410 if(dabs(pt(38)-1.0) .lt.eps9) go to 420
      if(pt(38)-1.0)135,420,415
  415 if(dabs(pt(38)-2.0) .lt.eps9) go to 435
      if(pt(38)-2.0)135,435,135
  420 if(dabs(pt(14)-2.0) .lt.eps9) go to 430
      if(pt(14)-2.0)425,430,135
  425 pt(39)=-pmac(med)
      pt(15)=-pmac(med)
      pt(27)=-pmac(med)
      pt(3)=-pmac(med)
      if(dabs(pt(2)-5.0)  .lt.eps9) go to 310
      if(pt(2)-5.0)135,310,135
  430 pt(3)=thpn(med)
      go to 381
  435 if(dabs(pt(14)-2.0) .lt.eps9) go to 445
      if(pt(14)-2.0)440,445,135
  440 pt(3)=hvp(med)
      pt(15)=hvp(med)
      pt(39)=hvp(med)
      pt(27)=thpn(med)
      go to 295
  445 pt(27)=-pmac(med)
      go to 386
      end subroutine


************************************************************************
*                                                                      *
