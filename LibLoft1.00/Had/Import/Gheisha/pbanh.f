*cmz :  3.14/16 13/03/89  14.48.40  by  nick van eijndhoven (cern)
*-- author :
      subroutine pbanh(nopt)
c *** anti proton annihilation at rest ***
c *** nve 04-mar-1988 cern geneva ***
c
c origin : h.fesefeldt (09-july-1987)
c
c nopt=0    no annihilation
c nopt=1    annih.in pi+ pi-
c nopt=2    annih.in pi0 pi0
c nopt=3    annih.in pi- pi0
c nopt=4    annih.in gamma gamma
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
      dimension brr(3)
      dimension rndm(3)
      data brr/0.125,0.25,0.5/
c
      pv(1,1)=0.
      pv(2,1)=0.
      pv(3,1)=0.
      pv(4,1)=abs(rmass(15))
      pv(5,1)=rmass(15)
      pv(6,1)=-1.
      pv(7,1)=tof
      pv(8,1)=ipart
      pv(9,1)=0.
      pv(10,1)=userw
      ier(86)=ier(86)+1
      isw=1
      call grndm(rndm,1)
      ran=rndm(1)
      if(ran.gt.brr(1)) isw=2
      if(ran.gt.brr(2)) isw=3
      if(ran.gt.brr(3)) isw=4
      nopt=isw
c**
c**  evaporation
c**
      call compo
      rmnve1=rmass(7)
      rmnve2=rmass(9)
      if (isw .eq. 2) rmnve1=rmass(8)
      if (isw .eq. 2) rmnve2=rmass(8)
      if (isw .eq. 3) rmnve1=rmass(9)
      if (isw .eq. 3) rmnve2=rmass(8)
      if (isw .eq. 4) rmnve1=rmass(1)
      if (isw .eq. 4) rmnve2=rmass(1)
      ek=rmass(14)+abs(rmass(15))-rmnve1-rmnve2
      tkin=exnu(ek)
      ek=ek-tkin
      if(ek.lt.0.0001) ek=0.0001
      ek=ek/2.
      en=ek+(rmnve1+rmnve2)/2.0
      s=amas*amas+rmass(14)**2+2.0*rmass(14)*en
      rs=sqrt(s)
      pcm=sqrt(abs(en*en-rmnve1*rmnve2))
      call grndm(rndm,2)
      phi=2.*pi*rndm(1)
      cost=-1.+2.*rndm(1)
      sint=sqrt(abs(1.-cost*cost))
      pv(1,2)=pcm*sint*cos(phi)
      pv(2,2)=pcm*sint*sin(phi)
      pv(3,2)=pcm*cost
      do 1 i=1,3
    1 pv(i,3)=-pv(i,2)
      pv(5,2)=rmnve1
      pv(5,3)=rmnve2
      if(isw.le.3) goto 2
      pv(5,2)=0.
      pv(5,3)=0.
    2 pv(4,2)=sqrt(pv(5,2)**2+pcm**2)
      pv(4,3)=sqrt(pv(5,3)**2+pcm**2)
      pv(7,2)=tof
      pv(7,3)=tof
      pv(9,2)=0.
      pv(9,3)=0.
      pv(10,2)=0.
      pv(10,3)=0.
      goto (21,22,23,24),isw
   21 pv(6,2)=1.
      pv(6,3)=-1.
      pv(8,2)=7.
      pv(8,3)=9.
      goto 25
   22 pv(6,2)=0.
      pv(6,3)=0.
      pv(8,2)=8.
      pv(8,3)=8.
      goto 25
   23 pv(6,2)=-1.
      pv(6,3)=0.
      pv(8,2)=9.
      pv(8,3)=8.
      goto 25
   24 pv(6,2)=0.
      pv(6,3)=0.
      pv(8,2)=1.
      pv(8,3)=1.
   25 nt=3
      if(atno2.lt.1.5) goto 40
      afc=0.312+0.200*log(log(s))+s**1.5/6000.
      targ=afc*(atno2**0.33 -1.0)
      cfa=0.025*((atno2-1.)/120.)*exp(-(atno2-1.)/120.)
      targ=1.
      tex=enp(1)
      if(tex.lt.0.001) goto 445
      black=(1.5+1.25*targ)*enp(1)/(enp(1)+enp(3))
      call poisso(black,nbl)
      if(ifix(targ)+nbl.gt.atno2) nbl=atno2-targ
      if(nt+nbl.gt.198) nbl=198-nt
      if(nbl.le.0) goto 445
      ekin=tex/nbl
      ekin2=0.
      call steep(xx)
      do 441 i=1,nbl
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
      if(ekin1.lt.0.) ekin1=0.001
      ipa1=16
      pnrat=1.-zno2/atno2
      call grndm(rndm,3)
      if(rndm(1).gt.pnrat) ipa1=14
      nt=nt+1
      cost=-1.+rndm(2)*2.
      sint=sqrt(abs(1.-cost*cost))
      phi=twpi*rndm(3)
      ipa(nt)=-ipa1
      pv(5,nt)=abs(rmass(ipa1))
      pv(6,nt)=rcharg(ipa1)
      pv(7,nt)=tof
      pv(8,nt)=ipa1
      pv(9,nt)=0.
      pv(10,nt)=0.
      pv(4,nt)=ekin1+pv(5,nt)
      pp=sqrt(abs(pv(4,nt)**2-pv(5,nt)**2))
      pv(1,nt)=pp*sint*sin(phi)
      pv(2,nt)=pp*sint*cos(phi)
      pv(3,nt)=pp*cost
  441 continue
  443 if(atno2.lt.230.) goto 445
      if(ek.gt.2.0) goto 445
      ii=nt+1
      kk=0
      eka=ek
      if(eka.gt.1.) eka=eka*eka
      if(eka.lt.0.1) eka=0.1
      ika=ifix(3.6/eka)
      do 444 i=1,nt
      ii=ii-1
      if(ipa(ii).ne.-14) goto 444
      ipa(ii)=-16
      ipa1  = 16
      pv(5,ii)=abs(rmass(ipa1))
      pv(6,ii)=rcharg(ipa1)
      pv(8,ii)=ipa1
      kk=kk+1
      if(kk.gt.ika) goto 445
  444 continue
c**
c** then also deuterons, tritons and alphas
c**
  445 tex=enp(3)
      if(tex.lt.0.001) goto 40
      black=(1.5+1.25*targ)*enp(3)/(enp(1)+enp(3))
      call poisso(black,nbl)
      if(nt+nbl.gt.198) nbl=198-nt
      if(nbl.le.0) goto 40
      ekin=tex/nbl
      ekin2=0.
      call steep(xx)
      do 442 i=1,nbl
      if(nt.eq.198) goto 442
      if(ekin2.gt.tex) goto 40
      call grndm(rndm,1)
      ran1=rndm(1)
      call normal(ran2)
      ekin1=-ekin*log(ran1)-cfa*(1.+0.5*ran2)
      if(ekin1.lt.0.0) ekin1=-0.010*log(ran1)
      ekin1=ekin1*xx
      ekin2=ekin2+ekin1
      if(ekin2.gt.tex) ekin1=tex-(ekin2-ekin1)
      if(ekin1.lt.0.) ekin1=0.001
      call grndm(rndm,3)
      cost=-1.+rndm(1)*2.
      sint=sqrt(abs(1.-cost*cost))
      phi=twpi*rndm(2)
      ran=rndm(3)
      ipa(nt+1)=-30
      if(ran.gt.0.60) ipa(nt+1)=-31
      if(ran.gt.0.90) ipa(nt+1)=-32
      inve=abs(ipa(nt+1))
      pv(5,nt+1)=rmass(inve)
      nt=nt+1
      pv(6,nt)=rcharg(inve)
      pv(7,nt)=tof
      pv(8,nt)=abs(ipa(nt))
      pv(9,nt)=0.
      pv(10,nt)=0.
      pv(4,nt)=pv(5,nt)+ekin1
      pp=sqrt(abs(pv(4,nt)**2-pv(5,nt)**2))
      pv(1,nt)=pp*sint*sin(phi)
      pv(2,nt)=pp*sint*cos(phi)
      pv(3,nt)=pp*cost
  442 continue
   40 intct=intct+1.
      call setcur(2)
      ntk=ntk+1
      if(nt.eq.2) go to 9999
      do 50 i=3,nt
      if(ntot.lt.nsize/12) goto 43
      go to 9999
   43 call settrk(i)
   50 continue
      call lengtx(3,pp)
      if(nprt(3))
     *write(newbcd,1001) xend,yend,zend,p,nch,pp,pv(6,3)
1001  format(1h0,'pb annihilation at rest  position',3(1x,f8.2),1x,
     * 'pi momenta,charges',2(1x,f8.4,1x,f4.1))
c
 9999 continue
      return
      end
