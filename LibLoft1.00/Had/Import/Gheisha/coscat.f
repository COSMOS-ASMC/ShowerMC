*cmz :  3.14/16 13/03/89  14.48.39  by  nick van eijndhoven (cern)
*-- author :
      subroutine coscat
c
c *** momentum generation for coherent elastic scattering ***
c *** nve 13-jul-1988 cern geneva ***
c
c origin : h.fesefeldt (03-dec-1986)
c
c approximation of bessel function for teta(lab)<=20 deg.
c is used . the nuclear radius is taken as r=1.25*e-13*(a)**1/3fm
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
      common/coscom/aa,bb,cc,dd,rr
c
c --- initialization flags for various gheisha routines ---
      common /kginit/ kginit(50)
c
c
      external fctcos
      dimension ff(20),atnox(3)
      dimension rndm(1)
c
      data atnox/9.,56.,207./
c
c --- initialization indicated by kginit(14) ---
      if (kginit(14) .ne. 0) go to 10
      kginit(14)=1
c
      if(.not.nprt(10)) goto 10
      write(newbcd,2001)
 2001 format(1h0,'ds/dt for coherent elastic scattering')
      do 3 l=1,3
      write(newbcd,2003) atnox(l),p
 2003 format(1h0,'calculated cross sections for a=',f5.1,' and p=',f8.2)
      do 2 i=1,20
      teta=(i-1)*pi/360.
      t=2.*p**2*(1.-cos(teta*1.d0))
      if(atnox(l).gt.62.) goto 4
      ff(i)=twpi*atnox(l)**1.63*exp(-14.5d0*atnox(l)**0.65*t)
     *     +twpi*1.4*atnox(l)**0.33*exp(-10.d0*t)
      goto 2
    4 ff(i)=twpi*atnox(l)**1.33*exp(-60.0d0*atnox(l)**0.33*t)
     *     +twpi*0.4*atnox(l)**0.40*exp(-10.d0*t)
    2 continue
      write(newbcd,2004) ff
 2004 format(1h ,10e12.3)
    3 continue
   10 if(p.lt.0.01) go to 9999
      if(atno2.lt.0.5) go to 9999
      ier(46)=ier(46)+1
      ran=ranres(dum)
      call vzero(ipa(1),100)
      ipa(1)=ipart
      if(atno2.gt.62.) goto 11
      aa=atno2**1.63
      bb=14.5*atno2**0.66
      cc=1.4*atno2**0.33
      dd=10.
      aa=aa/bb
      cc=cc/dd
      rr=(aa+cc)*ran
      goto 12
   11 aa=atno2**1.33
      bb=60.*atno2**0.33
      cc=0.4*atno2**0.40
      dd=10.
      aa=aa/bb
      cc=cc/dd
      rr=(aa+cc)*ran
   12 t1=-log(ran)/bb
      t2=-log(ran)/dd
      eps=0.001
      ind1=10
      call rtmi(t,val,fctcos,t1,t2,eps,ind1,ier1)
      if(ier1.eq.0) goto 14
      t=(3.*t1+t2)/4.
      ier(68)=ier(68)+1
   14 call grndm(rndm,1)
      phi=rndm(1)*twpi
      rr=t/(2.*p*p)
      if(rr.gt.1.) rr=0.
      cost=1.-rr
      dumnve=1.-cost*cost
      dumnve=max(dumnve,0.0)
      sint=sqrt(dumnve)
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
      pv(1,1)=p*sint*sin(phi)
      pv(2,1)=p*sint*cos(phi)
      pv(3,1)=p*cost
      pv(4,1)=en
      pv(5,1)=amas
      pv(6,1)=nch
      pv(7,1)=tof
      pv(8,1)=ipart
      pv(9,1)=0.
      pv(10,1)=0.
      call defs1(1,199,1)
      sinl1=sinl
      cosl1=cosl
      sinp1=sinp
      cosp1=cosp
      call setcur(1)
      if(nprt(4))
     *write(newbcd,1004) amas,p,sinl1,cosl1,sinp1,cosp1,sinl,cosl,
     *                   sinp,cosp,t1,t,t2,ier1
c
 1004 format(1h ,'coherent elastic scattering    mass ',f8.3,' momentum
     * ',f8.3/1h ,'direction ',4f10.4,' changed to ',4f10.4/
     *1h ,'t1,t,t2 ',3e10.3,' ier1 ',i2)
c
 9999 continue
      return
      end
