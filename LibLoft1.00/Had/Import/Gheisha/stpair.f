*cmz :  3.14/16 04/07/89  09.06.11  by  nick van eijndhoven (cern)
*-- author :
      subroutine stpair
c
c *** strange particle pair production ***
c *** nve 14-mar-1988 cern geneva ***
c
c origin : h.fesefeldt 16-dec-1987
c
c the same formula for <k kb> vs available energy
c                  and <k y>  vs available energy
c for all reactions.
c choose charge combinations k+ k- , k+ k0b, k0 k0b or k0 k-
c                            k+ y0, k0 y+, k0 y-
c for antibaryon induced reactions half of the cross sections
c kb yb pairs are produced
c charge is not conserved , no experimental data available for
c exclusive reactions, therefore some average behaviour assumed.
c the ratio l/sigma is taken as 3:1 (from experimental low energy)
c
      common/consts/ pi,twpi,pibtw,mp,mpi,mmu,mel,mkch,mk0,smp,smpi,
     $               smu,ct,ctkch,ctk0,
     $               ml0,msp,ms0,msm,mx0,mxm,ctl0,ctsp,ctsm,ctx0,ctxm,
     $               rmass(35),rcharg(35)
c
                     real mp,mpi,mmu,mel,mkch,mk0,
     *                    ml0,msp,ms0,msm,mx0,mxm
c
      common/curpar/weight(10),ddeltn,ifile,irun,nevt,nevent,shflag,
     *              ithst,ittot,itlst,ifrnd,tofcut,cmom(5),ceng(5),
     *              rs,s,enp(10),np,nm,nn,nr,no,nz,ipa(200),
     *              atno2,zno2
c
      common/result/xend,yend,zend,rca,rce,amas,nch,tof,px,py,pz,
     *              userw,intct,p,en,ek,amasq,deltn,itk,ntk,ipart,ind,
     *              lcalo,icel,sinl,cosl,sinp,cosp,
     *              xold,yold,zold,pold,pxold,pyold,pzold,
     *              xscat,yscat,zscat,pscat,pxscat,pyscat,pzscat
                    real nch,intct
c
      common/mat   / lmat,
     *               den(21),radlth(21),atno(21),zno(21),absl(21),
     *               cden(21),mden(21),x0den(21),x1den(21),rion(21),
     *               matid(21),matid1(21,24),parmat(21,10),
     *               ifrat,ifrac(21),frac1(21,10),den1(21,10),
     *               atno1(21,10),zno1(21,10)
c
      common/event / nsize,ncur,next,ntot,eve(1200)
c
      common/prntfl/inbcd,newbcd,inbin,newbin,npevt,nevtp,lprt,nprt(10)
                    logical lprt,nprt
c
      common/errcom/ ier(100)
c
      common /vecuty/ pv(10,200)
c
c
c
      real kkb,ky
      dimension kkb(9),ky(12),ipakkb(2,9),ipaky(2,12),ipakyb(2,12)
      dimension avkkb(12),avky(12),avrs(12)
      dimension rndm(1)
      data kkb/0.2500,0.3750,0.5000,0.5625,0.6250,0.6875,0.7500,
     *         0.8750,1.000/
      data ky /0.200,0.300,0.400,0.550,0.625,0.700,0.800,0.850,
     *         0.900,0.950,0.975,1.000/
      data ipakkb/10,13,10,11,10,12,11,11,11,12,12,11,12,12,
     *            11,13,12,13/
      data ipaky /18,10,18,11,18,12,20,10,20,11,20,12,21,10,
     *            21,11,21,12,22,10,22,11,22,12/
      data ipakyb/19,13,19,12,19,11,23,13,23,12,23,11,24,13,
     *            24,12,24,11,25,13,25,12,25,11/
      data avrs/3.,4.,5.,6.,7.,8.,9.,10.,20.,30.,40.,50./
      data avkkb/0.0015,0.005,0.012,0.0285,0.0525,0.075,0.0975,
     *           0.123,0.28,0.398,0.495,0.573/
      data avky /0.005,0.03,0.064,0.095,0.115,0.13,0.145,0.155,
     *           0.20,0.205,0.210,0.212/
c
      if(ipa(3).le.0) go to 9999
      ier(50)=ier(50)+1
      ipa1=abs(ipa(1))
      ipa2=abs(ipa(2))
c --- protection against annihilation processes ---
      if ((ipa1 .eq. 0) .or. (ipa2 .eq. 0)) go to 9999
      eab=rs-abs(rmass(ipa1))-abs(rmass(ipa2))
      if(eab.lt.1.) go to 9999
c**
c** choose random replacement of produced kaons (16.12.87)
      do 111 i=1,60
      if(ipa(i).eq.0) goto 112
  111 continue
  112 i=i-3
      call grndm(rndm,1)
      i3=3+ifix(rndm(1)*i)
  114 call grndm(rndm,1)
      i4=3+ifix(rndm(1)*i)
      if(i.eq.1) i4=4
      if(i3.eq.i4) goto 114
c
c *** choose random replacement of produced kaons (16.12.87) ***
c --- get rs bin ---
      do 1 i=2,12
      if (rs .le. avrs(i)) go to 2
 1    continue
      i1=11
      i2=12
      go to 3
c
 2    continue
      i1=i-1
      i2=i
c
c *** use linear interpolation or extrapolation by y=rc*x+b ***
 3    continue
      dxnve=avrs(i2)-avrs(i1)
      dynve=log(avkkb(i2))-log(avkkb(i1))
      rcnve=dynve/dxnve
      bnve=log(avkkb(i1))-rcnve*avrs(i1)
      avk=rcnve*rs+bnve
      dxnve=avrs(i2)-avrs(i1)
      dynve=log(avky(i2))-log(avky(i1))
      rcnve=dynve/dxnve
      bnve=log(avky(i1))-rcnve*avrs(i1)
      avy=rcnve*rs+bnve
c
      avk=exp(avk)
      avy=exp(avy)
      if(avk+avy.le.0.) go to 9999
      if(ipa1.lt.14) avy=avy/2.
      if(ipa2.lt.14) avy=0.
      avy=avy+avk
      call grndm(rndm,1)
      ran=rndm(1)
      if(ran.lt.avk) goto 10
      if(ran.lt.avy) goto 20
      go to 9999
   10 if((eab-1.).lt.0.) go to 9999
      call grndm(rndm,1)
      ran=rndm(1)
      do 11 i=1,9
      if(ran.lt.kkb(i)) goto 12
   11 continue
      go to 9999
   12 ipa(i3)=ipakkb(1,i)
      ipa(i4)=ipakkb(2,i)
      goto 30
   20 if((eab-1.6).lt.0.) go to 9999
      call grndm(rndm,1)
      ran=rndm(1)
      do 21 i=1,12
      if(ran.lt.ky(i)) goto 22
   21 continue
      go to 9999
   22 if(ipa(1).lt.14) goto 23
      call grndm(rndm,1)
      if(rndm(1).lt.0.5) goto 23
      ipa1=abs(ipa(1))
      ipa(1)=ipaky(1,i)
      if(ipa1.eq.15) goto 25
      if(ipa1.eq.17) goto 25
      if(ipa1.eq.19) goto 25
      if(ipa1.gt.22) goto 25
      goto 24
   25 ipa(1)=ipakyb(1,i)
      ipa(i3)=ipakyb(2,i)
      goto 30
   23 ipa(2)=ipaky(1,i)
   24 ipa(i3)=ipaky(2,i)
c** check the available energy
   30 eab=rs
      ij=0
      do 31 i=1,60
      if(ipa(i).eq.0) goto 31
      ipa1=abs(ipa(i))
      eab=eab-abs(rmass(ipa1))
      ij=ij+1
      if(eab.lt.0.) goto 35
   31 continue
c      if (nprt(4)) write(newbcd,1003) (ipa(j),j=1,ij)
      go to 9999
   35 i=i-1
      l=i-1
      if(l.le.0) go to 9999
      do 36 j=i,60
   36 ipa(j)=0
c      if (nprt(4)) write(newbcd,1002) (ipa(j),j=1,l)
c
 1002 format(' *stpair* kkb/ky pair production not enough energy',
     $ ' reduce number of particles ',2x,20i3)
 1003 format(' *stpair* kkb/ky pair production energy sufficient',
     $ ' number of particles ',2x,20i3)
c
 9999 continue
c
      return
      end
