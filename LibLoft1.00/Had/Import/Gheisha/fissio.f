*cmz :  3.14/16 28/09/90  10.04.18  by  nick van eijndhoven (cern)
*-- author :
      function fissio(ek1)
c
c *** generation of photons and neutrons by fission ***
c *** nve 04-mar-1988 cern geneva ***
c
c origin : h.fesefeldt (21-mar-1987)
c
c the physics is based on u(238)
c for other materials extrapolations are used
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
c --- initialization flags for various gheisha routines ---
      common /kginit/ kginit(50)
c
c
      dimension spneut(10)
      dimension rndm(2)
      save spneut
      data spneut/10*0./
c
c --- initialization indicated by kginit(15) ---
      if (kginit(15) .ne. 0) go to 10
      kginit(15)=1
c
      xx=1.-0.5
      xxx=sqrt(2.29*xx)
      spneut(1)=exp(-xx/0.965)*(exp(xxx)-exp(-xxx))/2.
      do 1 i=2,10
      xx=i*1.-0.5
      xxx=sqrt(2.29*xx)
    1 spneut(i)=spneut(i-1)+exp(-xx/0.965)*(exp(xxx)-exp(-xxx))/2.
      do 2 i=1,10
    2 spneut(i)=spneut(i)/spneut(10)
c** in this routine we use mev as unit for energy and momentum
   10 nt=0
      ier(82)=ier(82)+1
      nd=ind+1
      pv(1,200)=px*p
      pv(2,200)=py*p
      pv(3,200)=pz*p
      pv(4,200)=en
      pv(5,200)=abs(amas)
      pv(6,200)=nch
      pv(7,200)=tof
      pv(8,200)=ipart
      pv(9,200)=0.
      pv(10,200)=userw
      pv(1,199)=0.
      pv(2,199)=0.
      pv(3,199)=0.
      pv(4,199)=atomas(atno(nd),zno(nd))
      pv(5,199)=pv(4,199)
      pv(6,199)=zno(nd)
      pv(7,199)=tof
      pv(8,199)=0.
      pv(9,199)=0.
      pv(10,199)=0.
      call add(200,199,198)
      pv(1,198)=-pv(1,198)
      pv(2,198)=-pv(2,198)
      pv(3,198)=-pv(3,198)
c** number of neutrons and photons
      fissio=0.
      e1=ek1*1000.
      if(e1.lt.1.0) e1=1.0
      avern=2.569+0.559*log(e1)
c**   take the following value if photofission is not included
      if(ifix(parmat(ind+1,8)).eq.0)
     *avern=2.569+0.900*log(e1)
      averg=9.500+0.600*log(e1)
      call normal(ran)
      nn=ifix(avern+ran*1.23+0.5)
      call normal(ran)
      ng=ifix(averg+ran*3.+0.5)
      if(nn.lt.1) nn=1
      if(ng.lt.1) ng=1
      exn=0.
      exg=0.
c** distribute kinetic energy
      do 15 i=1,nn
      call grndm(rndm,1)
      ran=rndm(1)
      do 11 j=1,10
      if(ran.lt.spneut(j)) goto 12
   11 continue
      j=10
   12 call grndm(rndm,1)
      ekin=(j-1)*1.+rndm(1)
      exn=exn+ekin
      pv(4,i)=ekin+rmass(16)*1000.
      pv(5,i)=rmass(16)*1000.
      pv(6,i)=0.
c** emission time for neutrons =0.
      pv(7,i)=tof
      pv(8,i)=16.
      pv(9,i)=0.
      pv(10,i)=0.
   15 continue
      nt=nn
      do 20 i=1,ng
      call grndm(rndm,1)
      ran=rndm(1)
      nt=nt+1
      pv(4,nt)=-0.87*log(ran)
      exg=exg+pv(4,nt)
      pv(5,nt)=0.
      pv(6,nt)=0.
c     ran=rndm(1)
c** emission time for photons= 2.5 e-8 sec
c     pv(7,nt)=tof-500.*log(ran)
c** changed 30.7.85
      pv(7,nt)=tof
      pv(8,nt)=1.
      pv(9,nt)=0.
      pv(10,nt)=0.
   20 continue
      if(nt.eq.0) go to 9999
      ex=exn+exg
      if(nprt(4))
     *write(newbcd,2000) atno(ind+1),nn,ng,ex
      fissio=ex/1000.
      do 49 i=1,nt
      pv(5,i)=pv(5,i)/1000.
      pv(4,i)=pv(4,i)/1000.
      call grndm(rndm,2)
      cost=-1.+2.*rndm(1)
      sint=sqrt(abs(1.-cost*cost))
      phi=rndm(2)*twpi
      pp=sqrt(abs(pv(4,i)**2-pv(5,i)**2))
      pv(1,i)=pp*sint*sin(phi)
      pv(2,i)=pp*sint*cos(phi)
      pv(3,i)=pp*cost
      call lor(i,198,i)
   49 continue
      intct=intct+1.
      nmode=6
      if(ipart.eq.1) nmode=7
      do 50 i=1,nt
      if(ntot.lt.nsize/12) goto 43
      ier(39)=ier(39)+1
      go to 9999
   43 call settrk(i)
   50 continue
c
 2000 format(1h ,'nuclear fission on material ',f6.1,
     *', neutrons, photons produced= ',2i3,' with ',f8.4,' mev total ene
     *rgy')
c
 9999 continue
      end
