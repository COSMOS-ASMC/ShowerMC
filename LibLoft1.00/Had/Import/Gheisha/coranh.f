*cmz :  3.14/16 13/03/89  14.48.43  by  nick van eijndhoven (cern)
*-- author :
      subroutine coranh(nihil,nfl)
c
c *** nuclear interactions for heavy fragments ***
c *** nve 06-may-1988 cern geneva ***
c
c origin : h.fesefeldt (09-july-1987)
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
      nihil=0
      if(amas.gt.0.)   go to 9999
      if(ipart.lt.14)  go to 9999
      if(ipa(1).ge.14) go to 9999
      if(ipa(2).ge.14) go to 9999
      nihil=1
c**
c**  do not be confused, this has nothing to do with relativistic
c**  kinematic
c
      tarmas=rmass(14)
      if (nfl .eq. 2) tarmas=rmass(16)
      ekcor=1.
      if(ek.gt.1.) ekcor=1./ek
      ek=2.*tarmas+ek*(1.+ekcor/atno2)
      en=ek+abs(amas)
      p =sqrt(abs(en*en-amas*amas))
      s =amas*amas+tarmas**2+2.0*tarmas*en
      rs=sqrt(s)
      enp(5)=ek
      enp(6)=en
      enp(7)=p
      enp(8)=s
      enp(9)=rs
c**
c**  evaporation
c**
      tkin=exnu(ek)
      enp(5)=ek-tkin
      if(enp(5).lt.0.0001) enp(5)=0.0001
      enp(6)=enp(5)+abs(amas)
      enp(7)=enp(6)*enp(6)-amasq
      enp(7)=sqrt(abs(enp(7)))
      enp(8)=amasq+rmass(14)**2+2.*rmass(14)*enp(6)
      enp(9)=sqrt(enp(8))
c**  check available energy for first interaction
      if(enp(5).gt.ceng(3)) go to 9999
      enp(5)=0.
      enp(6)=abs(amas)
      enp(7)=0.
      enp(8)=4.*rmass(14)**2
      enp(9)=2.*rmass(14)
c
 9999 continue
      return
      end
