*cmz :  3.14/16 13/03/89  14.48.39  by  nick van eijndhoven (cern)
*-- author :
      subroutine captur(nopt)
c
c *** routine for capture of neutral baryons ***
c *** nve 04-mar-1988 cern geneva ***
c
c origin : h.fesefeldt (02-dec-1986)
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
      dimension rndm(3)
c
      nopt=1
      ier(81)=ier(81)+1
      pv(1,1)=px*p
      pv(2,1)=py*p
      pv(3,1)=pz*p
      pv(4,1)=en
      pv(5,1)=abs(amas)
      pv(6,1)=nch
      pv(7,1)=tof
      pv(8,1)=ipart
      pv(9,1)=0.
      pv(10,1)=userw
      nd=ind+1
      pv(1,2)=0.
      pv(2,2)=0.
      pv(3,2)=0.
      pv(4,2)=atomas(atno(nd),zno(nd))
      pv(5,2)=pv(4,2)
      pv(6,2)=zno(nd)
      pv(7,2)=tof
      pv(8,2)=0.
      pv(9,2)=0.
      pv(10,2)=0.
      call add(1,2,200)
      pv(1,200)=-pv(1,200)
      pv(2,200)=-pv(2,200)
      pv(3,200)=-pv(3,200)
      call normal(ran)
      p=0.0065+ran*0.0010
      call grndm(rndm,3)
      cost=-1.+rndm(1)*2.
      sint=sqrt(abs(1.-cost*cost))
      phi=twpi*rndm(2)
      pv(1,3)=p*sint*sin(phi)
      pv(2,3)=p*sint*cos(phi)
      pv(3,3)=p*cost
      pv(4,3)=p
      pv(5,3)=0.
      pv(6,3)=0.
      pv(8,3)=1.
      pv(9,3)=0.
      pv(10,3)=0.
      ran=rndm(3)
      tof=tof-480.*log(ran)
      pv(7,3)=tof
      call lor(3,200,3)
      nt=3
      xp=0.008-p
      if(xp.lt.0.) goto 9
      nt=4
      call grndm(rndm,2)
      cost=-1.+rndm(1)*2.
      sint=sqrt(abs(1.-cost*cost))
      phi=twpi*rndm(2)
      pv(1,4)=xp*sint*sin(phi)
      pv(2,4)=xp*sint*cos(phi)
      pv(3,4)=xp*cost
      pv(4,4)=xp
      pv(5,4)=0.
      pv(6,4)=0.
      pv(7,4)=tof
      pv(8,4)=1.
      pv(9,4)=0.
      pv(10,4)=0.
      call lor(4,200,4)
    9 intct=intct+1.
      call setcur(3)
      ntk=ntk+1
      if(nt.eq.4) call settrk(4)
c
      return
      end
