*cmz :  3.14/16 13/03/89  14.48.41  by  nick van eijndhoven (cern)
*-- author :
      subroutine setcur(ntr)
c
c *** storage of current track parameters ***
c *** nve 16-mar-1988 cern geneva ***
c
c origin : h.fesefeldt (26-jan-1984)
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
      dimension rndm(1)
c
      call lengtx(ntr,p)
      amas=pv(5,ntr)
      amasq=amas*amas
      nch=pv(6,ntr)
      tof=pv(7,ntr)
      ipart=ifix(pv(8,ntr)+0.1)
      if(pv(10,ntr).ne.0.) userw=pv(10,ntr)
      px=0.
      py=0.
      pz=0.
      if(p.lt.1.e-10) goto 4
      px=pv(1,ntr)/p
      py=pv(2,ntr)/p
      pz=pv(3,ntr)/p
    4 en=pv(4,ntr)
      ek=en-abs(amas)
      sinl=pz
      cosl=sqrt(abs(1.-sinl*sinl))
      if(abs(cosl).lt.1.e-10) goto 1
      sinp=py/cosl
      cosp=px/cosl
      goto 2
    1 call grndm(rndm,1)
      phi=rndm(1)*twpi
      sinp=sin(phi)
      cosp=cos(phi)
    2 if(nprt(3).or.nprt(4).or.nprt(5))
     *write(newbcd,1001) xend,yend,zend,rca,rce,amas,nch,tof,px,py,pz,
     *userw,intct,p,en,ek,amasq,deltn,itk,ntk,ipart,ind,lcalo,icel,
     *sinl,cosl,sinp,cosp
      return
 1001 format(1h ,'track parameter changed:',3f8.2,1x,2f7.0,1x,f8.3,1x,
     *f3.0,1x,f6.0,1x,3f6.3,1x,f10.0,1x,f5.0/10x,4f8.3,1x,f8.5,1x,6i5,
     *4f8.3)
      end
