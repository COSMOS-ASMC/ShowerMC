*cmz :          10/01/91  19.07.05  by  federico carminati
*cmz :  3.14/16 11/10/90  18.30.06  by  rene brun
*-- author :
      subroutine twoclu(ippp,nfl,avern)
c
c *** generation of x- and pt- values for all produced particles ***
c *** nve 01-aug-1988 cern geneva ***
c
c origin : h.fesefeldt (11-oct-1987)
c
c a simple two cluster model is used
c this should be sufficient for low energy interactions
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
      common/genin /tecm,amass(18),npg,kgenev
      common/genout/pcm(5,18),wgt
c
c
      real nucsup
      dimension side(200),c1par(5),g1par(5),nucsup(5)
      dimension rndm(3)
      data c1par/0.6,0.6,0.35,0.15,0.10/
      data g1par/2.6,2.6,1.8,1.30,1.20/
      data nucsup/1.0,0.8,0.6,0.5,0.4/
c     data cb/3.0/
      data cb/0.01/
      bpp(x)=4.000+1.600*log(x)
c
      ek=enp(5)
      en=enp(6)
      p=enp(7)
      s=enp(8)
      rs=enp(9)
      cfa=0.025*((atno2-1.)/120.)*exp(-(atno2-1.)/120.)
      if(p.lt.0.001) goto 60
      nt=0
c**
c** check mass-indices for all particles
c**
      do 1 i=1,100
      if(ipa(i).eq.0) goto 1
      nt=nt+1
      ipa(nt)=ipa(i)
    1 continue
      call vzero(ipa(nt+1),200-nt)
c**
c** set the effectice 4-momentum-vector for interaction
c**
      pv(1,199)=p*px
      pv(2,199)=p*py
      pv(3,199)=p*pz
      pv(4,199)=en
      pv(5,199)=amas
      pv(6,199)=nch
      pv(7,199)=tof
      pv(8,199)=ipart
      pv(9,199)=0.
      pv(10,199)=userw
      ier(48)=ier(48)+1
c**
c** distribute particles in forward and backward hemisphere of cms
c** of the hadron nucleon interaction
c**
      side(1)= 1.
      side(2)=-1.
      targ=0.
      ifor=1
      iback=1
      do 3 i=1,nt
      if (i .le. 2) go to 78
      side(i)=1.
      call grndm(rndm,1)
      if (rndm(1) .lt. 0.5) side(i)=-1.
      if (side(i) .lt. 0.) go to 76
c
c --- particle in forward hemisphere ---
 77   continue
      ifor=ifor+1
      if (ifor .le. 18) go to 78
c
c --- change it to backward ---
      side(i)=-1.
      ifor=ifor-1
      iback=iback+1
      go to 78
c
c --- particle in backward hemisphere ---
 76   continue
      iback=iback+1
      if (iback .le. 18) go to 78
c
c --- change it to forward ---
      side(i)=1.
      iback=iback-1
      ifor=ifor+1
c**
c** suppression of charged pions for various reasons
c**
   78 if(ipart.eq.15.or.ipart.ge.17) goto 3
      if(abs(ipa(i)).ge.10) goto 3
      if(abs(ipa(i)).eq. 8) goto 3
      call grndm(rndm,1)
      if(rndm(1).gt.(10.-p)/6.) goto 3
      call grndm(rndm,1)
      if(rndm(1).gt.atno2/300.) goto 3
      ipa(i)=14
      call grndm(rndm,1)
      if(rndm(1).gt.zno2/atno2) ipa(i)=16
      targ=targ+1.
    3 continue
      tb=2.*iback
      call grndm(rndm,1)
      if(rs.lt.(2.0+rndm(1))) tb=(2.*iback+nt)/2.
c**
c** nucleons + some pions from intranuclear cascade
c**
      afc=0.312+0.200*log(log(s))
      xtarg=afc*(atno2**0.33-1.0)*tb
      if(xtarg.le.0.) xtarg=0.01
      call poisso(xtarg,ntarg)
      nt2=nt+ntarg
      if(nt2.le.170) goto 2
      nt2=170
      ntarg=nt2-nt
    2 continue
      if(nprt(4))
     *write(newbcd,3001) ntarg,nt
      nt1=nt+1
      if(ntarg.eq.0) goto 51
      ipx=ifix(p/3.)+1
      if(ipx.gt.5) ipx=5
      do 4 i=nt1,nt2
      call grndm(rndm,1)
      ran=rndm(1)
      if(ran.lt.nucsup(ipx)) goto 52
      call grndm(rndm,1)
      ipa(i)=-(7+ifix(rndm(1)*3.0))
      goto 4
   52 ipa(i)=-16
      pnrat=1.-zno2/atno2
      call grndm(rndm,1)
      if(rndm(1).gt.pnrat) ipa(i)=-14
      targ=targ+1.
    4 side(i)=-2.
      nt=nt2
c**
c** choose masses and charges for all particles
c**
   51 do 5 i=1,nt
      ipa1=abs(ipa(i))
      pv(5,i)=rmass(ipa1)
      pv(6,i)=rcharg(ipa1)
      pv(7,i)=1.
      if(pv(5,i).lt.0.) pv(7,i)=-1.
      pv(5,i)=abs(pv(5,i))
    5 continue
c**
c** mark leading strange particles
c**
      lead=0
      if(ipart.lt.10.or.ipart.eq.14.or.ipart.eq.16) goto 6
      ipa1=abs(ipa(1))
      if(ipa1.lt.10.or.ipa1.eq.14.or.ipa1.eq.16) goto 531
      lead=ipa1
      goto 6
  531 ipa1=abs(ipa(2))
      if(ipa1.lt.10.or.ipa1.eq.14.or.ipa1.eq.16) goto 6
      lead=ipa1
c**
c** check available kinetic energy , change hemisphere for particles
c** until it fits
c**
    6 if(nt.le.1) goto 60
      tavai=0.
      do 7 i=1,nt
      if(side(i).lt.-1.5) goto 7
      tavai=tavai+abs(pv(5,i))
    7 continue
      if(tavai.lt.rs) goto 12
      if(nprt(4))
     *write(newbcd,3002) (ipa(i),i=1,20),(side(i),i=1,20),tavai,rs
 3002 format(' *twoclu* check available energies'/
     *       1h ,20i5/1h ,20f5.0/1h ,'tavai,rs ',2f10.3)
      do 10 i=1,nt
      ii=nt-i+1
      if(side(ii).lt.-1.5) goto 10
      if(ii.eq.nt) goto 11
      nt1=ii+1
      nt2=nt
      do 8 j=nt1,nt2
      ipa(j-1)=ipa(j)
      side(j-1)=side(j)
      do 8 k=1,10
    8 pv(k,j-1)=pv(k,j)
      goto 11
   10 continue
   11 side(nt)=0.
      ipa(nt)=0
      nt=nt-1
      goto 6
   12 if(nt.le.1) goto 60
      b=bpp(p)
      if(b.lt.cb) b=cb
c**
c** choose masses for the 3 cluster: 1. forward cluster
c**   2. backward meson cluster  3. backward nucleon cluster
c**
      rmc0=0.
      rmd0=0.
      rme0=0.
      ntc=0
      ntd=0
      nte=0
      do 31 i=1,nt
      if(side(i).gt.0.) rmc0=rmc0+abs(pv(5,i))
      if(side(i).gt.0.) ntc =ntc +1
      if(side(i).lt.0..and.side(i).gt.-1.5) rmd0=rmd0+abs(pv(5,i))
      if(                  side(i).lt.-1.5) rme0=rme0+abs(pv(5,i))
      if(side(i).lt.0..and.side(i).gt.-1.5) ntd =ntd +1
      if(                  side(i).lt.-1.5) nte =nte +1
   31 continue
   32 call grndm(rndm,1)
      ran=rndm(1)
      rmc=rmc0
      if(ntc.le.1) goto 33
      ntc1=ntc
      if(ntc1.gt.5) ntc1=5
      rmc=-log(1.-ran)
      gpar=g1par(ntc1)
      cpar=c1par(ntc1)
      dumnve=gpar
      if (dumnve .eq. 0.0) dumnve=1.0e-10
      rmc=rmc0+rmc**cpar/dumnve
   33 rmd=rmd0
      if(ntd.le.1) goto 34
      ntd1=ntd
      if(ntd1.gt.5) ntd1=5
      call grndm(rndm,1)
      ran=rndm(1)
      rmd=-log(1.-ran)
      gpar=g1par(ntd1)
      cpar=c1par(ntd1)
      dumnve=gpar
      if (dumnve .eq. 0.0) dumnve=1.0e-10
      rmd=rmd0+rmd**cpar/dumnve
   34 if(rmc+rmd.lt.rs) goto 35
      if (rmc.le.rmc0.and.rmd.le.rmd0) then
         hnrmdc = 0.999*rs/(rmc+rmd)
         rmd = rmd*hnrmdc
         rmc = rmc*hnrmdc
      else
         rmc=0.1*rmc0+0.9*rmc
         rmd=0.1*rmd0+0.9*rmd
      end if
      goto 34
   35 if(nte.le.0) goto 38
      rme=rme0
      if(nte.eq.1) goto 38
      nte1=nte
      if(nte1.gt.5) nte1=5
      call grndm(rndm,1)
      ran=rndm(1)
      rme=-log(1.-ran)
      gpar=g1par(nte1)
      cpar=c1par(nte1)
      dumnve=gpar
      if (dumnve .eq. 0.0) dumnve=1.0e-10
      rme=rme0+rme**cpar/dumnve
c**
c** set beam , target of first interaction in cms
c**
   38 pv(1,181)=0.
      pv(2,181)=0.
      pv(3,181)=p
      pv(5,181)=abs(amas)
      pv(4,181)=sqrt(p*p+amas*amas)
      pv(1,182)=0.
      pv(2,182)=0.
      pv(3,182)=0.
      pv(4,182)=mp
      pv(5,182)=mp
     
c** transform into cms.
     
      call add(181,182,180)
      call lor(181,180,181)
      call lor(182,180,182)
      pf=(s+rmd*rmd-rmc*rmc)**2 - 4*s*rmd*rmd
      if(pf.lt.0.0001) pf=0.0001
      dumnve=2.0*rs
      if (dumnve .eq. 0.0) dumnve=1.0e-10
      pf=sqrt(pf)/dumnve
      if(nprt(4)) write(6,2002) pf,rmc,rmd,rs
c**
c** set final state masses and energies in cms
c**
      pv(5,183)=rmc
      pv(5,184)=rmd
      pv(4,183)=sqrt(pf*pf+pv(5,183)*pv(5,183))
      pv(4,184)=sqrt(pf*pf+pv(5,184)*pv(5,184))
c**
c** set |t| and |tmin|
c**
      t=-1.0e10
      call grndm(rndm,1)
      if (b .ne. 0.0) t=log(1.-rndm(1))/b
      call lengtx(181,pin)
      tacmin=(pv(4,181)-pv(4,183))**2-(pin-pf)**2
c**
c** caculate (sin(teta/2.)**2 and cos(teta), set azimuth angle phi
c**
      dumnve=4.0*pin*pf
      if (dumnve .eq. 0.0) dumnve=1.0e-10
      ctet=-(t-tacmin)/dumnve
      ctet=1.0-2.0*ctet
      if (ctet .gt. 1.0) ctet=1.0
      if (ctet .lt. -1.0) ctet=-1.0
      dumnve=1.0-ctet*ctet
      if (dumnve .lt. 0.0) dumnve=0.0
      stet=sqrt(dumnve)
      call grndm(rndm,1)
      phi=rndm(1)*twpi
c**
c** calculate final state momenta in cms
c**
      pv(1,183)=pf*stet*sin(phi)
      pv(2,183)=pf*stet*cos(phi)
      pv(3,183)=pf*ctet
      pv(1,184)=-pv(1,183)
      pv(2,184)=-pv(2,183)
      pv(3,184)=-pv(3,183)
c**
c** simulate backward nucleon cluster in lab. system and transform in
c** cms.
c**
      if(nte.eq.0) goto 28
      ga=1.2
      ekit1=0.04
      ekit2=0.6
      if(ek.gt.5.) goto 666
      ekit1=ekit1*ek**2/25.
      ekit2=ekit2*ek**2/25.
  666 a=(1.-ga)/(ekit2**(1.-ga)-ekit1**(1.-ga))
      do 29 i=1,nt
      if(side(i).gt.-1.5) goto 29
      call grndm(rndm,3)
      ran=rndm(1)
      ekit=(ran*(1.-ga)/a+ekit1**(1.-ga))**(1./(1.-ga))
      pv(4,i)=ekit+pv(5,i)
      dumnve=abs(pv(4,i)**2-pv(5,i)**2)
      pp=sqrt(dumnve)
      ran=rndm(2)
      cost=log(2.23*ran+0.383)/0.96
      if (cost .lt. -1.0) cost=-1.0
      if (cost .gt. 1.0) cost=1.0
      dumnve=1.0-cost*cost
      if (dumnve .lt. 0.0) dumnve=0.0
      sint=sqrt(dumnve)
      phi=twpi*rndm(3)
      pv(1,i)=pp*sint*sin(phi)
      pv(2,i)=pp*sint*cos(phi)
      pv(3,i)=pp*cost
      call lor(i,180,i)
   29 continue
c**
c** fragmentation of forward cluster and backward meson cluster
c**
   28 pv(1,1)=pv(1,183)
      pv(2,1)=pv(2,183)
      pv(3,1)=pv(3,183)
      pv(4,1)=pv(4,183)
      pv(1,2)=pv(1,184)
      pv(2,2)=pv(2,184)
      pv(3,2)=pv(3,184)
      pv(4,2)=pv(4,184)
      do 17 i=185,186
      do 16 j=1,3
   16 pv(j,i)=-pv(j,i-2)
      do 17 j=4,5
   17 pv(j,i)= pv(j,i-2)
      kgenev=1
      if(ntc.le.1) goto 26
      tecm= pv(5,183)
      npg=0
      do 18 i=1,nt
      if(side(i).lt.0.) goto 18
      npg=npg+1
      amass(npg)=abs(pv(5,i))
   18 continue
      if(nprt(4)) write(newbcd,2004) tecm,npg,(amass(i),i=1,npg)
      call phasp
      npg=0
      do 19 i=1,nt
      if(side(i).lt.0.) goto 19
      npg=npg+1
      pv(1,i)=pcm(1,npg)
      pv(2,i)=pcm(2,npg)
      pv(3,i)=pcm(3,npg)
      pv(4,i)=pcm(4,npg)
      if(nprt(4)) write(newbcd,2001) i,(pv(j,i),j=1,5)
      call lor(i,185,i)
      if(nprt(4)) write(newbcd,2001) i,(pv(j,i),j=1,10),ipa(i),side(i)
   19 continue
   26 if(ntd.le.1) goto 27
      tecm= pv(5,184)
      npg=0
      do 20 i=1,nt
      if(side(i).gt.0..or.side(i).lt.-1.5) goto 20
      npg=npg+1
      amass(npg)=abs(pv(5,i))
   20 continue
      if(nprt(4)) write(newbcd,2004) tecm,npg,(amass(i),i=1,npg)
      call phasp
      npg=0
      do 21 i=1,nt
      if(side(i).gt.0..or.side(i).lt.-1.5) goto 21
      npg=npg+1
      pv(1,i)=pcm(1,npg)
      pv(2,i)=pcm(2,npg)
      pv(3,i)=pcm(3,npg)
      pv(4,i)=pcm(4,npg)
      if(nprt(4)) write(newbcd,2001) i,(pv(j,i),j=1,5)
      call lor(i,186,i)
      if(nprt(4)) write(newbcd,2001) i,(pv(j,i),j=1,10),ipa(i),side(i)
   21 continue
c**
c** lorentz transformation in lab system
c**
   27 targ=0.
      do 36 i=1,nt
      if(pv(5,i).gt.0.5) targ=targ+1.
      call lor(i,182,i)
   36 continue
      if(targ.lt.0.5) targ=1.
c**
c** sometimes the leading strange particles are lost , set them back
c**
      if(lead.eq.0) goto 6085
      do 6081 i=1,nt
      if(abs(ipa(i)).eq.lead) goto 6085
 6081 continue
      i=1
      if(lead.ge.14.and.abs(ipa(2)).ge.14) i=2
      if(lead.lt.14.and.abs(ipa(2)).lt.14) i=2
      ipa(i)=lead
      ekin=pv(4,i)-abs(pv(5,i))
      pv(5,i)=rmass(lead)
      pv(7,i)=1.
      if(pv(5,i).lt.0.) pv(7,i)=-1.
      pv(5,i)=abs(pv(5,i))
      pv(6,i)=rcharg(lead)
      pv(4,i)=pv(5,i)+ekin
      call lengtx(i,pp)
      dumnve=abs(pv(4,i)**2-pv(5,i)**2)
      pp1=sqrt(dumnve)
c
      if (pp .ge. 1.0e-6) go to 8000
      call grndm(rndm,2)
      rthnve=pi*rndm(1)
      phinve=twpi*rndm(2)
      pv(1,i)=pp1*sin(rthnve)*cos(phinve)
      pv(2,i)=pp1*sin(rthnve)*sin(phinve)
      pv(3,i)=pp1*cos(rthnve)
      go to 8001
 8000 continue
      pv(1,i)=pv(1,i)*pp1/pp
      pv(2,i)=pv(2,i)*pp1/pp
      pv(3,i)=pv(3,i)*pp1/pp
 8001 continue
c
c** for various reasons, the energy balance is not sufficient,
c** check that,  energy balance, angle of final system e.t.c.
 6085 kgenev=1
      pv(1,184)=0.
      pv(2,184)=0.
      pv(3,184)=p
      pv(4,184)=sqrt(p*p+amas*amas)
      pv(5,184)=abs(amas)
      ekin0=pv(4,184)-pv(5,184)
      pv(1,185)=0.
      pv(2,185)=0.
      pv(3,185)=0.
      pv(4,185)=mp*targ
      pv(5,185)=pv(4,185)
      ekin=pv(4,184)+pv(4,185)
      i=184
      if(nprt(4)) write(newbcd,2001) i,(pv(j,i),j=1,5)
      i=185
      if(nprt(4)) write(newbcd,2001) i,(pv(j,i),j=1,5)
      call add(184,185,186)
      call lor(184,186,184)
      call lor(185,186,185)
      tecm=pv(4,184)+pv(4,185)
      npg=nt
      pv(1,188)=0.
      pv(2,188)=0.
      pv(3,188)=0.
      pv(4,188)=0.
      pv(5,188)=0.
      ekin1=0.
      do 598 i=1,npg
      if(nprt(4)) write(newbcd,2001) i,(pv(j,i),j=1,10),ipa(i),side(i)
      call add(188,i,188)
      ekin1=ekin1+pv(4,i)-pv(5,i)
      ekin=ekin-pv(5,i)
      if(i.gt.18) goto 598
      amass(i)=pv(5,i)
  598 continue
      if(npg.gt.18) goto 597
      call phasp
      ekin=0.
      do 599 i=1,npg
      pv(1,187)=pcm(1,i)
      pv(2,187)=pcm(2,i)
      pv(3,187)=pcm(3,i)
      pv(4,187)=pcm(4,i)
      pv(5,187)=amass(i)
      call lor(187,185,187)
  599 ekin=ekin+pv(4,187)-pv(5,187)
      call ang(188,184,cost,teta)
      if(nprt(4)) write(newbcd,2003) teta,ekin0,ekin1,ekin
c**
c** make shure, that  kinetic energies are correct
c** the 3. cluster is not produced within proper kinematics!!!
c** ekin= kinetic energy theoretically
c** ekin1= kinetic energy simulated
c**
  597 if(ekin1.eq.0.) goto 600
      pv(1,187)=0.
      pv(2,187)=0.
      pv(3,187)=0.
      pv(4,187)=0.
      pv(5,187)=0.
      wgt=ekin/ekin1
      ekin1=0.
      do 602 i=1,nt
      ekin=pv(4,i)-pv(5,i)
      ekin=ekin*wgt
      pv(4,i)=ekin+pv(5,i)
      dumnve=abs(pv(4,i)**2-pv(5,i)**2)
      pp=sqrt(dumnve)
      call lengtx(i,pp1)
c
      if (pp1 .ge. 1.0e-6) go to 8002
      call grndm(rndm,2)
      rthnve=pi*rndm(1)
      phinve=twpi*rndm(2)
      pv(1,i)=pp*sin(rthnve)*cos(phinve)
      pv(2,i)=pp*sin(rthnve)*sin(phinve)
      pv(3,i)=pp*cos(rthnve)
      go to 8003
 8002 continue
      pv(1,i)=pv(1,i)*pp/pp1
      pv(2,i)=pv(2,i)*pp/pp1
      pv(3,i)=pv(3,i)*pp/pp1
 8003 continue
c
      ekin1=ekin1+ekin
      call add(187,i,187)
  602 continue
      call ang(187,184,cost,teta)
      if(nprt(4)) write(newbcd,2003) teta,ekin0,ekin1
c**
c** rotate in direction of z-axis, see comments in 'genxpt'
c**
  600 pv(1,187)=0.
      pv(2,187)=0.
      pv(3,187)=0.
      pv(4,187)=0.
      pv(5,187)=0.
      do 596 i=1,nt
  596 call add(187,i,187)
*          call rannor(ran1,ran2)
      call grndm(rndm,2)
      ry=rndm(1)
      rz=rndm(2)
      rx=6.283185*rz
      a1=sqrt(-2.*log(ry))
      ran1=a1*sin(rx)
      ran2=a1*cos(rx)
      pv(1,187)=pv(1,187)+ran1*0.020*targ
      pv(2,187)=pv(2,187)+ran2*0.020*targ
      call defs(184,187,188)
      pv(1,187)=0.
      pv(2,187)=0.
      pv(3,187)=0.
      pv(4,187)=0.
      pv(5,187)=0.
      do 595 i=1,nt
      call trac(i,188,i)
  595 call add(187,i,187)
      call ang(187,184,cost,teta)
      if(nprt(4)) write(newbcd,2003) teta
c**
c** rotate in direction of primary particle
c**
      dekin=0.
      npions=0
      ek1=0.
      do 25 i=1,nt
      call defs1(i,199,i)
      if(nprt(4)) write(newbcd,2001) i,(pv(j,i),j=1,10),ipa(i),side(i)
      if(atno2.lt.1.5) goto 25
      call lengtx(i,pp)
      ekin=pv(4,i)-abs(pv(5,i))
      call normal(ran)
      ekin=ekin-cfa*(1.+0.5*ran)
      if (ekin .lt. 1.0e-6) ekin=1.0e-6
      call steeq(xxh,i)
      dekin=dekin+ekin*(1.-xxh)
      ekin=ekin*xxh
      if(abs(ipa(i)).ge.7.and.abs(ipa(i)).le.9) npions=npions+1
      if(abs(ipa(i)).ge.7.and.abs(ipa(i)).le.9) ek1=ek1+ekin
      pp1=sqrt(ekin*(ekin+2.*abs(pv(5,i))))
      pv(4,i)=ekin+abs(pv(5,i))
c
      if (pp .ge. 1.0e-6) go to 8004
      call grndm(rndm,2)
      rthnve=pi*rndm(1)
      phinve=twpi*rndm(2)
      pv(1,i)=pp1*sin(rthnve)*cos(phinve)
      pv(2,i)=pp1*sin(rthnve)*sin(phinve)
      pv(3,i)=pp1*cos(rthnve)
      go to 8005
 8004 continue
      pv(1,i)=pv(1,i)*pp1/pp
      pv(2,i)=pv(2,i)*pp1/pp
      pv(3,i)=pv(3,i)*pp1/pp
 8005 continue
c
   25 continue
      if(ek1.eq.0.) goto 23
      if(npions.le.0) goto 23
      dekin=1.+dekin/ek1
      do 22 i=1,nt
      if(abs(ipa(i)).lt.7.or.abs(ipa(i)).gt.9) goto 22
      call lengtx(i,pp)
      ekin=pv(4,i)-abs(pv(5,i))
      ekin=ekin*dekin
      if (ekin .lt. 1.0e-6) ekin=1.0e-6
      pp1=sqrt(ekin*(ekin+2.*abs(pv(5,i))))
      pv(4,i)=ekin+abs(pv(5,i))
c
      if (pp .ge. 1.0e-6) go to 8006
      call grndm(rndm,2)
      rthnve=pi*rndm(1)
      phinve=twpi*rndm(2)
      pv(1,i)=pp1*sin(rthnve)*cos(phinve)
      pv(2,i)=pp1*sin(rthnve)*sin(phinve)
      pv(3,i)=pp1*cos(rthnve)
      go to 8007
 8006 continue
      pv(1,i)=pv(1,i)*pp1/pp
      pv(2,i)=pv(2,i)*pp1/pp
      pv(3,i)=pv(3,i)*pp1/pp
 8007 continue
c
   22 continue
   23 igen=0
      if(atno2.lt.1.5) goto 40
c**
c** add black track particles
c**
      call selfab(sprob)
      tex=enp(1)
      spall=targ
      if(tex.lt.0.001) goto 445
      black=(1.5+1.25*targ)*enp(1)/(enp(1)+enp(3))
      call poisso(black,nbl)
      if(nprt(4))
     *write(newbcd,3003) nbl,tex
      if(ifix(targ)+nbl.gt.atno2) nbl=atno2-targ
      if(nt+nbl.gt.198) nbl=198-nt
      if(nbl.le.0) goto 445
      ekin=tex/nbl
      ekin2=0.
      call steep(xx)
      do 441 i=1,nbl
      call grndm(rndm,1)
      if(rndm(1).lt.sprob) goto 441
      if(nt.eq.198) goto 441
      if(ekin2.gt.tex) goto 443
      call grndm(rndm,1)
      ran1=rndm(1)
      call normal(ran2)
      ekin1=-ekin*log(ran1)-cfa*(1.+0.5*ran2)
      if(ekin1.lt.0.0) ekin1=-0.010*log(ran1)
      ekin1=ekin1*xx
      ekin2=ekin2+ekin1
      if(ekin2.gt.tex) ekin1=tex-(ekin2-ekin1)
      if (ekin1 .lt. 0.0) ekin1=1.0e-6
      ipa1=16
      pnrat=1.-zno2/atno2
      call grndm(rndm,3)
      if(rndm(1).gt.pnrat) ipa1=14
      nt=nt+1
      spall=spall+1.
      cost=-1.0+rndm(2)*2.0
      dumnve=1.0-cost*cost
      if (dumnve .lt. 0.0) dumnve=0.0
      sint=sqrt(dumnve)
      phi=twpi*rndm(3)
      ipa(nt)=-ipa1
      side(nt)=-4.
      pv(5,nt)=abs(rmass(ipa1))
      pv(6,nt)=rcharg(ipa1)
      pv(7,nt)=1.
      pv(4,nt)=ekin1+pv(5,nt)
      dumnve=abs(pv(4,nt)**2-pv(5,nt)**2)
      pp=sqrt(dumnve)
      pv(1,nt)=pp*sint*sin(phi)
      pv(2,nt)=pp*sint*cos(phi)
      pv(3,nt)=pp*cost
  441 continue
  443 if(atno2.lt.10.) goto 445
      if(ek.gt.2.0) goto 445
      ii=nt+1
      kk=0
      eka=ek
      if(eka.gt.1.) eka=eka*eka
      if(eka.lt.0.1) eka=0.1
      ika=3.6*exp((zno2**2/atno2-35.56)/6.45)/eka
      if(ika.le.0) go to 445
      do 444 i=1,nt
      ii=ii-1
      if(ipa(ii).ne.-14) goto 444
      ipa(ii)=-16
      ipa1  = 16
      pv(5,ii)=abs(rmass(ipa1))
      pv(6,ii)=rcharg(ipa1)
      kk=kk+1
      if(kk.gt.ika) goto 445
  444 continue
  445 tex=enp(3)
      if(tex.lt.0.001) goto 40
      black=(1.5+1.25*targ)*enp(3)/(enp(1)+enp(3))
      call poisso(black,nbl)
      if(nt+nbl.gt.198) nbl=198-nt
      if(nbl.le.0) goto 40
      ekin=tex/nbl
      ekin2=0.
      call steep(xx)
      if(nprt(4))
     *write(newbcd,3004) nbl,tex
      do 442 i=1,nbl
      call grndm(rndm,1)
      if(rndm(1).lt.sprob) goto 442
      if(nt.eq.198) goto 442
      if(ekin2.gt.tex) goto 40
      call grndm(rndm,1)
      ran1=rndm(1)
      call normal(ran2)
      ekin1=-ekin*log(ran1)-cfa*(1.+0.5*ran2)
      if(ekin1.lt.0.0) ekin1=-0.005*log(ran1)
      ekin1=ekin1*xx
      ekin2=ekin2+ekin1
      if(ekin2.gt.tex) ekin1=tex-(ekin2-ekin1)
      if (ekin1 .lt. 0.0) ekin1=1.0e-6
      call grndm(rndm,3)
      cost=-1.0+rndm(1)*2.0
      dumnve=1.0-cost*cost
      if (dumnve .lt. 0.0) dumnve=0.0
      sint=sqrt(dumnve)
      phi=twpi*rndm(2)
      ran=rndm(3)
      ipa(nt+1)=-30
      if(ran.gt.0.60) ipa(nt+1)=-31
      if(ran.gt.0.90) ipa(nt+1)=-32
      side(nt+1)=-4.
      pv(5,nt+1)=(abs(ipa(nt+1))-28)*mp
      spall=spall+pv(5,nt+1)*1.066
      if(spall.gt.atno2) goto 40
      nt=nt+1
      pv(6,nt)=1.
      if(ipa(nt).eq.-32) pv(6,nt)=2.
      pv(7,nt)=1.
      pv(4,nt)=pv(5,nt)+ekin1
      dumnve=abs(pv(4,nt)**2-pv(5,nt)**2)
      pp=sqrt(dumnve)
      pv(1,nt)=pp*sint*sin(phi)
      pv(2,nt)=pp*sint*cos(phi)
      pv(3,nt)=pp*cost
  442 continue
c**
c** store on event common
c**
   40 call grndm(rndm,1)
      if(rs.gt.(4.+rndm(1)*1.)) goto 42
      do 41 i=1,nt
      call lengtx(i,etb)
      if(etb.lt.p) goto 41
      etf=p
      pv(4,i)=sqrt(pv(5,i)**2+etf**2)
      dumnve=etb
      if (dumnve .eq. 0.0) dumnve=1.0e-10
      etf=etf/dumnve
      pv(1,i)=pv(1,i)*etf
      pv(2,i)=pv(2,i)*etf
      pv(3,i)=pv(3,i)*etf
   41 continue
   42 ekin=pv(4,200)-abs(pv(5,200))
      ekin1=pv(4,199)-abs(pv(5,199))
      ekin2=0.
      call tdelay(tof1)
      call grndm(rndm,1)
      ran=rndm(1)
      tof=tof-tof1*log(ran)
      do 44 i=1,nt
      ekin2=ekin2+pv(4,i)-abs(pv(5,i))
      if(pv(7,i).lt.0.) pv(5,i)=-pv(5,i)
      pv(7,i)=tof
      pv(8,i)=abs(ipa(i))
      pv(9,i)=0.
   44 pv(10,i)=0.
      if(nprt(4)) write(newbcd,2006) nt,ekin,enp(1),enp(3),ekin1,ekin2
      intct=intct+1.
      nmode=3
      if(spall.lt.0.5.and.atno2.gt.1.5) nmode=14
      call setcur(nt)
      ntk=ntk+1
      if(nt.eq.1) goto 300
      do 50 ii=2,nt
      i=ii-1
      if(ntot.lt.nsize/12) goto 43
      go to 9999
   43 call settrk(i)
   50 continue
 300  continue
      go to 9999
c**
c** it is not possible to produce a proper two cluster final state.
c** continue with quasi elastic scattering
c**
   60 if(nprt(4)) write(newbcd,2005)
      do 61 i=3,200
   61 ipa(i)=0
      ipa(1)=ipart
      ipa(2)=14
      if(nfl.eq.2) ipa(2)=16
      call twob(ippp,nfl,avern)
      go to 9999
c
 2000 format(' *twoclu* cms parameters of final state particles',
     $ ' after ',i3,' trials')
 2001 format(' *twoclu* track',2x,i3,2x,10f8.2,2x,i3,2x,f3.0)
 2002 format(' *twoclu* momentum ',f8.3,' masses ',2f8.4,' rs ',f8.4)
 2003 format(' *twoclu* teta,ekin0,ekin1,ekin ',4f10.4)
 2004 format(' *twoclu* tecm,npb,masses: ',f10.4,1x,i3,1x,8f10.4/
     $ 1h ,26x,15x,8f10.4)
 2005 format(' *twoclu* number of final state particles',
     $ ' less than 2 ==> continue with 2-body scattering')
 2006 format(' *twoclu*  comp.',1x,i5,1x,5f7.2)
 3001 format(' *twoclu* nuclear excitation ',i5,' particles produced',
     $ ' in addition to',i5,' normal particles')
 3003 format(' *twoclu* ',i3,' black track particles produced',
     $ ' with total kinetic energy of ',f8.3,' gev')
 3004 format(' *twoclu* ',i5,' heavy fragments with total energy of ',
     $ f8.4,' gev')
c
 9999 continue
      return
      end
