*cmz :  3.14/16 13/03/89  14.48.40  by  nick van eijndhoven (cern)
*-- author :
      subroutine pimabs(nopt)
c
c *** charged pion absorption by a nucleus ***
c *** nve 04-mar-1988 cern geneva ***
c
c origin : h.fesefeldt (09-july-1987)
c
c panofsky ratio (pi- p --> n pi0/pi- p --> n gamma) = 3/2
c for capture on proton (hydrogen),
c star production for heavier elements
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
      dimension rndm(4)
c
      call compo
      pv(1,1)=0.
      pv(2,1)=0.
      pv(3,1)=0.
      pv(4,1)=rmass(9)
      pv(5,1)=rmass(9)
      pv(6,1)=-1.
      pv(7,1)=tof
      pv(8,1)=ipart
      pv(9,1)=0.
      pv(10,1)=userw
      ier(87)=ier(87)+1
      if(atno2.gt.1.5) goto 30
      call grndm(rndm,2)
      ran=rndm(1)
      isw=1
      if(ran.lt.0.33) isw=2
      nopt=isw
      ran=rndm(2)
      tof1=-25.*log(ran)
      tof1=20.*tof1
      if(isw.eq.1) goto 20
      pcm=0.02
      call grndm(rndm,2)
      phi=2.*pi*rndm(1)
      cost=-1.+2.*rndm(2)
      sint=sqrt(abs(1.-cost*cost))
      pv(1,2)=pcm*sint*cos(phi)
      pv(2,2)=pcm*sint*sin(phi)
      pv(3,2)=pcm*cost
      pv(4,2)=pcm
      pv(5,2)=0.
      pv(7,2)=tof+tof1
      pv(9,2)=0.
      pv(10,2)=0.
      pv(6,2)=0.
      pv(8,2)=1.
      goto 21
   20 pv(1,2)=pv(1,1)
      pv(2,2)=pv(2,1)
      pv(3,2)=pv(3,1)
      pv(4,2)=pv(4,1)
      pv(5,2)=pv(5,1)
      pv(6,2)=0.
      pv(7,2)=tof+tof1
      pv(8,2)=8.
      pv(9,2)=0.
      pv(10,2)=0.
   21 intct=intct+1.
      call setcur(2)
      ntk=ntk+1
      go to 9999
c**
c** star production for pion absorption in heavy elements
c**
   30 enp(1)=0.0135
      enp(3)=0.0058
      nt=1
      tex=enp(1)
      black=0.5*log(atno2)
      call poisso(black,nbl)
      if(nbl.le.0) nbl=1
      if(nprt(3))
     *write(newbcd,3003) nbl,tex
      if(nt+nbl.gt.198) nbl=198-nt
      ekin=tex/nbl
      ekin2=0.
      do 31 i=1,nbl
      if(nt.eq.198) goto 31
      call grndm(rndm,4)
      ran2=rndm(1)
      ekin1=-ekin*log(ran2)
      ekin2=ekin2+ekin1
      ipa1=16
      pnrat=1.-zno2/atno2
      if(rndm(2).gt.pnrat) ipa1=14
      nt=nt+1
      cost=-1.+rndm(3)*2.
      sint=sqrt(abs(1.-cost*cost))
      phi=twpi*rndm(4)
      ipa(nt)=-ipa1
      pv(5,nt)=abs(rmass(ipa1))
      pv(6,nt)=rcharg(ipa1)
      pv(7,nt)=2.
      pv(4,nt)=ekin1+pv(5,nt)
      pp=sqrt(abs(pv(4,nt)**2-pv(5,nt)**2))
      pv(1,nt)=pp*sint*sin(phi)
      pv(2,nt)=pp*sint*cos(phi)
      pv(3,nt)=pp*cost
      if(ekin2.gt.tex) goto 33
   31 continue
   33 tex=enp(3)
      black=0.50*log(atno2)
      call poisso(black,nbl)
      if(nt+nbl.gt.198) nbl=198-nt
      if(nbl.le.0) nbl=1
      ekin=tex/nbl
      ekin2=0.
      if(nprt(3))
     *write(newbcd,3004) nbl,tex
      do 32 i=1,nbl
      if(nt.eq.198) goto 32
      call grndm(rndm,4)
      ran2=rndm(1)
      ekin1=-ekin*log(ran2)
      ekin2=ekin2+ekin1
      nt=nt+1
      cost=-1.+rndm(2)*2.
      sint=sqrt(abs(1.-cost*cost))
      phi=twpi*rndm(3)
      ran=rndm(4)
      ipa(nt)=-30
      if(ran.gt.0.60) ipa(nt)=-31
      if(ran.gt.0.90) ipa(nt)=-32
      inve=abs(ipa(nt))
      pv(5,nt)=rmass(inve)
      pv(6,nt)=rcharg(inve)
      pv(7,nt)=2.
      pv(4,nt)=pv(5,nt)+ekin1
      pp=sqrt(abs(pv(4,nt)**2-pv(5,nt)**2))
      pv(1,nt)=pp*sint*sin(phi)
      pv(2,nt)=pp*sint*cos(phi)
      pv(3,nt)=pp*cost
      if(ekin2.gt.tex) goto 40
   32 continue
c**
c** store on event common
c**
   40 call grndm(rndm,1)
      ran=rndm(1)
      tof1=-25.*log(ran)
      tof1=20.*tof1
      do 41 i=2,nt
      if(pv(7,i).lt.0.) pv(5,i)=-pv(5,i)
      pv(7,i)=tof+tof1
      pv(8,i)=abs(ipa(i))
      pv(9,i)=0.
   41 pv(10,i)=0.
      intct=intct+1.
      call setcur(2)
      ntk=ntk+1
      if(nt.eq.2) go to 9999
      do 50 i=3,nt
      if(ntot.lt.nsize/12) goto 43
      go to 9999
   43 call settrk(i)
   50 continue
c
 3003 format(1h ,i3,' black track particles produced with total kinetic
     * energy of ',f8.3,' gev')
 3004 format(1h ,i5,' heavy fragments with total energy of',f8.4,' gev')
c
 9999 continue
      return
      end