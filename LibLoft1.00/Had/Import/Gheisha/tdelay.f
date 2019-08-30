*cmz :  3.14/16 28/09/90  10.22.58  by  nick van eijndhoven (cern)
*-- author :
      subroutine tdelay(x)
c
c *** time delay for nuclear reactions ***
c *** nve 16-mar-1988 cern geneva ***
c
c origin : h.fesefeldt (01-feb-1984)
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
c
      x=0.
      if(atno2.lt.1.5) return
      if(atno2.gt.230.) return
      if(ek.gt.0.2) return
      x=500.*exp(-ek/0.04)
      return
      end
