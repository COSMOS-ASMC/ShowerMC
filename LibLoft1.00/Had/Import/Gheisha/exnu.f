*cmz :  3.14/16 10/01/91  19.15.42  by  federico carminati
*cmz :  3.14/16 13/03/89  14.48.39  by  nick van eijndhoven (cern)
*-- author :
      function exnu(ek1)
c
c *** nuclear evaporation as function of atomic number atno ***
c *** and kinetic energy ekin of primary particle ***
c *** nve 04-mar-1988 cern geneva ***
c
c origin : h.fesefeldt (10-dec-1986)
c
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
      dimension rndm(2)
c
      exnu=0.
      if(atno2.lt.1.5) go to 9999
      magic=0
      if(int(zno2+0.1).eq.82) magic=1
      ekin1=ek1
      if(ekin1.lt.0.1) ekin1=0.1
      if(ekin1.gt.4.) ekin1=4.
c**   0.35 value at 1 gev
c**   0.05 value at 0.1 gev
      cfa=(0.35-0.05)/2.3
      cfa= 0.35+cfa*log(ekin1)
      if(cfa.lt.0.15) cfa=0.15
      exnu=7.716*cfa*exp(-cfa)
      atno3=atno2
      if(atno3.gt.120.) atno3=120.
      cfa=((atno3-1.)/120.)*exp(-(atno3-1.)/120.)
      exnu=exnu*cfa
      fpdiv=1.-0.25*ekin1**2
      if(fpdiv.lt.0.50) fpdiv=0.50
      gfa=2.0*((atno2-1.)/70.)*exp(-(atno2-1.)/70.)
      enp(1)=exnu*fpdiv
      enp(3)=exnu-enp(1)
    4 call normal(ran1)
      call normal(ran2)
      if(magic.eq.1) then
         ran1=0.
         ran2=0.
      end if
      enp(1)=enp(1)*(1.+ran1*gfa)
      if(enp(1).lt.0.) enp(1)=0.
      enp(3)=enp(3)*(1.+ran2*gfa)
      if(enp(3).lt.0.) enp(3)=0.
    5 exnu=enp(1)+enp(3)
      if(exnu.lt.ek1) goto 10
      call grndm(rndm,2)
      enp(1)=enp(1)*(1.-0.5*rndm(1))
      enp(3)=enp(3)*(1.-0.5*rndm(2))
      goto 5
 10   continue
c
 9999 continue
      return
      end
