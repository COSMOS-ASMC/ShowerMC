*cmz :          14/01/91  11.27.47  by  federico carminati
*cmz :  3.14/16 28/09/90  10.33.21  by  nick van eijndhoven (cern)
*-- author :
      function cinema(ek1)
c
c *** inelasticity in nuclear interactions as a function ***
c *** of atomic number atno2 and kinetic energy ek1 ***
c *** nve 12-jul-1988 cern geneva ***
c
c origin : h.fesefeldt (14-oct-1987)
c
c the functional dependence and the parameters have been obtained
c by study of various nuclear structure models.
c but: it is of course an interpolation as function of atomic
c      number, for certain nuclei a different description may be
c      more adequate. detailed tests have been performed for
c      fe, cu, pb ,u and some mixtures like nai, bgo, concrete.
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
c --- boundary limits for arguments of intrinsic functions ---
c --- xl denotes lower bound whereas xu denotes upper bound ---
      common /limits/ expxl,expxu
c
c
      cinema=0.0
      nd=ind+1
      ala=log(atno2)
      alek1=log(ek1)
      sig1=0.50
      sig2=0.50
      em=0.2390+0.0408*ala**2
      if (em. gt. 1.0) em=1.0
      cinem=0.0019*ala**3
      if(cinem.gt.0.15) cinem=0.15
      if (parmat(nd,10) .ge. 0.01) cinem=cinem*parmat(nd,10)
c
      if (alek1 .gt. em) go to 1
c
      corr=-(alek1-em)**2/(2.0*sig1**2)
      if (corr .lt. expxl) corr=expxl
      if (corr .gt. expxu) corr=expxu
      dum1=-ek1*cinem
      dum2=abs(dum1)
      dum3=exp(corr)
      cinema=0.0
      if (dum2 .ge. 1.0) cinema=dum1*dum3
      if ((dum2 .lt. 1.0) .and. (dum3 .gt. 1.0e-10)) cinema=dum1*dum3
      go to 2
c
 1    continue
      corr=-(alek1-em)**2/(2.0*sig2**2)
      if (corr .lt. expxl) corr=expxl
      if (corr .gt. expxu) corr=expxu
      dum1=-ek1*cinem
      dum2=abs(dum1)
      dum3=exp(corr)
      cinema=0.0
      if (dum2 .ge. 1.0) cinema=dum1*dum3
      if ((dum2 .lt. 1.0) .and. (dum3 .gt. 1.0e-10)) cinema=dum1*dum3
c
 2    continue
      if (cinema .lt. -ek1) cinema=-ek1
c
      return
      end
