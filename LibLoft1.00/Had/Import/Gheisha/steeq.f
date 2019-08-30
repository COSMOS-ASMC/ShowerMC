*cmz :  3.14/16 28/09/90  10.13.20  by  nick van eijndhoven (cern)
*-- author :
      subroutine steeq(xxh,ipv)
c
c *** corrections for single particle spectra (shower particles) ***
c *** nve 16-mar-1988 cern geneva
c
c origin : h.fesefeldt (06-sep-1985)
c
      common/curpar/weight(10),ddeltn,ifile,irun,nevt,nevent,shflag,
     *              ithst,ittot,itlst,ifrnd,tofcut,cmom(5),ceng(5),
     *              rs,s,enp(10),np,nm,nn,nr,no,nz,ipa(200),
     *              atno2,zno2
c
      common /vecuty/ pv(10,200)
c
      common/result/xend,yend,zend,rca,rce,amas,nch,tof,px,py,pz,
     *              userw,intct,p,en,ek,amasq,deltn,itk,ntk,ipart,ind,
     *              lcalo,icel,sinl,cosl,sinp,cosp,
     *              xold,yold,zold,pold,pxold,pyold,pzold,
     *              xscat,yscat,zscat,pscat,pxscat,pyscat,pzscat
                    real nch,intct
c
c
      dimension alem(7),val0(7)
      dimension rndm(1)
c**   data   em/ 4.0 , 10.  , 15.  , 20.  ,  30. , 100. , 1000./
      data alem/ 1.40, 2.30 , 2.70 , 3.00 , 3.40 , 4.60 , 7.00 /
      data val0/ 0.00, 0.40 , 0.48 , 0.51 , 0.54 , 0.60 , 0.65 /
c
      xxh=1.
c
      if ((ipart .ne. 7) .and. (ipart .ne. 9)) go to 9999
      if (abs(ipa(ipv)) .ne. 8) go to 9999
      call grndm(rndm,1)
      if (rndm(1) .gt. log(atno2)) go to 9999
      ekw=pv(4,200)-abs(pv(5,200))
      alekw=log(ekw)
      if (alekw .le. alem(1)) go to 9999
c
c --- get energy bin ---
      do 1 i=2,7
      if (alekw .lt. alem(i)) go to 2
 1    continue
      xxh=val0(7)
      go to 3
c
c *** use linear interpolation or extrapolation by y=rc*x+b ***
 2    continue
      i1=i-1
      i2=i
      dxnve=alem(i2)-alem(i1)
      dynve=val0(i2)-val0(i1)
      rcnve=dynve/dxnve
      bnve=val0(i1)-rcnve*alem(i1)
      xxh=rcnve*alekw+bnve
c
 3    continue
      xxh=1.-xxh
c
 9999 continue
      return
      end
