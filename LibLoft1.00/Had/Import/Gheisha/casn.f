*cmz :  3.14/16 13/03/89  14.48.42  by  nick van eijndhoven (cern)
*-- author :
      subroutine casn(k,int,nfl)
c
c *** cascade of neutron ***
c *** nve 04-may-1988 cern geneva ***
c
c origin : h.fesefeldt (13-sep-1987)
c
c n  undergoes interaction with nucleon within nucleus.
c check if energetically possible to produce pions/kaons.
c if not assume nuclear excitation occurs and input particle
c is degraded in energy.    no other particles produced.
c if reaction is possible find correct number of pions/protons/
c neutrons produced using an interpolation to multiplicity data.
c replace some pions or protons/neutrons by kaons or strange baryons
c according to average multiplicity per inelastic reactions.
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
      common/prntfl/inbcd,newbcd,inbin,newbin,npevt,nevtp,lprt,nprt(10)
                    logical lprt,nprt
c
c --- boundary limits for arguments of intrinsic functions ---
c --- xl denotes lower bound whereas xu denotes upper bound ---
      common /limits/ expxl,expxu
c
c --- initialization flags for various gheisha routines ---
      common /kginit/ kginit(50)
c
c
      real n
      dimension pmul(2,1200),anorm(2,60),supp(10),cech(10),b(2)
      dimension rndm(1)
      save pmul,anorm
      data supp/0.,0.4,0.55,0.65,0.75,0.82,0.86,0.90,0.94,0.98/
      data cech/0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0./
      data b/0.35,0.0/,c/1.25/
c
c --- initialization indicated by kginit(17) ---
      if (kginit(17) .ne. 0) go to 10
      kginit(17)=1
c
c --- initialize pmul and anorm arrays ---
      do 9000 j=1,1200
      do 9001 i=1,2
      pmul(i,j)=0.0
      if (j .le. 60) anorm(i,j)=0.0
 9001 continue
 9000 continue
c
c** compute normalization constants
c** for n as target
c
      l=0
      do 1 np1=1,20
      np=np1-1
      npp1=np1+2
      do 1 nm1=np1,npp1
      nm=nm1-1
      do 1 nz1=1,20
      nz=nz1-1
      l=l+1
      if(l.gt.1200) goto 1
      nprot= -np+nm
      nneut=2-nprot
      nt=np+nm+nz
      if(nt.le.0.or.nt.gt.60) goto 1
      pmul(1,l)=pmltpc(np,nm,nz,nt,b(2),c)
      nprotf=nfac(nprot)
      nneutf=nfac(nneut)
      pmul(1,l)=pmul(1,l)/(nprotf*nneutf)
      anorm(1,nt)=anorm(1,nt)+pmul(1,l)
    1 continue
c** for p as target
      l=0
      do 2 np1=1,20
      np=np1-1
      nmm1=np1-1
      if(nmm1.le.1) nmm1=1
      npp1=np1+1
      do 2 nm1=nmm1,npp1
      nm=nm1-1
      do 2 nz1=1,20
      nz=nz1-1
      l=l+1
      if(l.gt.1200) goto 2
      nprot=1-np+nm
      nneut=2-nprot
      nt=np+nm+nz
      if(nt.le.0.or.nt.gt.60) goto 2
      pmul(2,l)=pmltpc(np,nm,nz,nt,b(1),c)
      nprotf=nfac(nprot)
      nneutf=nfac(nneut)
      pmul(2,l)=pmul(2,l)/(nprotf*nneutf)
      anorm(2,nt)=anorm(2,nt)+pmul(2,l)
    2 continue
      do 3 i=1,60
      if(anorm(1,i).gt.0.) anorm(1,i)=1./anorm(1,i)
      if(anorm(2,i).gt.0.) anorm(2,i)=1./anorm(2,i)
    3 continue
      if(.not.nprt(10)) goto 10
      write(newbcd,2001)
      do 4 nfl=1,2
      write(newbcd,2002) nfl
      write(newbcd,2003) (anorm(nfl,i),i=1,60)
      write(newbcd,2003) (pmul(nfl,i),i=1,1200)
    4 continue
c**  choose proton or neutron as target
   10 nfl=2
      call grndm(rndm,1)
      if(rndm(1).lt.zno2/atno2) nfl=1
      tarmas=rmass(14)
      if (nfl .eq. 2) tarmas=rmass(16)
      s=amasq+tarmas**2+2.0*tarmas*en
      rs=sqrt(s)
      enp(8)=amasq+tarmas**2+2.0*tarmas*enp(6)
      enp(9)=sqrt(enp(8))
      eab=rs-tarmas-rmass(16)
c**  elastic scattering
      np=0
      nm=0
      nz=0
      n=0.
      ncech=0
      if(int.eq.2) goto 20
c** introduce charge exchange reaction pn --> np
      if(nfl.eq.2) goto 100
      iplab=ifix(p*2.5)+1
      if(iplab.gt.10) iplab=10
      call grndm(rndm,1)
      if(rndm(1).gt.cech(iplab)/atno2**0.42) goto 100
      ncech=1
      goto 100
c**  check if energetically possible to produce one extra pion in react.
  20  if (eab .le. rmass(7)) goto 55
c**  suppression of high multiplicity events at low momentum
      ieab=ifix(eab*5.)+1
      if(ieab.gt.10) goto 22
      call grndm(rndm,1)
      if(rndm(1).lt.supp(ieab)) goto 22
      n=1.
      goto (24,23),nfl
 23   continue
      test=-(1+b(2))**2/(2.0*c**2)
      if (test .gt. expxu) test=expxu
      if (test .lt. expxl) test=expxl
      w0=exp(test)/2.0
      wm=exp(test)
      call grndm(rndm,1)
      ran=rndm(1)
      np=0
      nm=0
      nz=1
      if(ran.lt.w0/(w0+wm)) goto 100
      np=0
      nm=1
      nz=0
      goto 100
 24   continue
      test=-(1+b(1))**2/(2.0*c**2)
      if (test .gt. expxu) test=expxu
      if (test .lt. expxl) test=expxl
      w0=exp(test)
      wp=exp(test)/2.0
      test=-(-1+b(1))**2/(2.0*c**2)
      if (test .gt. expxu) test=expxu
      if (test .lt. expxl) test=expxl
      wm=exp(test)/2.0
      wt=w0+wp+wm
      wp=w0+wp
      call grndm(rndm,1)
      ran=rndm(1)
      np=0
      nm=0
      nz=1
      if(ran.lt.w0/wt) goto 100
      np=1
      nm=0
      nz=0
      if(ran.lt.wp/wt) goto 100
      np=0
      nm=1
      nz=0
      goto 100
c
   22 aleab=log(eab)
c** no. of total particles vs sqrt(s)-2*mp
      n=3.62567+0.665843*aleab+0.336514*aleab*aleab
     * +0.117712*aleab*aleab*aleab+0.0136912*aleab*aleab*aleab*aleab
      n=n-2.
c** normalization constant for  kno-distribution
      anpn=0.
      do 21 nt=1,60
      test=-(pi/4.0)*(nt/n)**2
      if (test .gt. expxu) test=expxu
      if (test .lt. expxl) test=expxl
      anpn=anpn+pi*nt*exp(test)/(2.0*n*n)
   21 continue
      anpn=1./anpn
c** p or n as target
      call grndm(rndm,1)
      ran=rndm(1)
      excs=0.
      goto (40,30),nfl
c** for n as target
   30 l=0
      do 31 np1=1,20
      np=np1-1
      npp1=np1+2
      do 31 nm1=np1,npp1
      nm=nm1-1
      do 31 nz1=1,20
      nz=nz1-1
      l=l+1
      if(l.gt.1200) goto 31
      nt=np+nm+nz
      if(nt.le.0.or.nt.gt.60) goto 31
      test=-(pi/4.0)*(nt/n)**2
      if (test .gt. expxu) test=expxu
      if (test .lt. expxl) test=expxl
      dum1=anpn*pi*nt*pmul(1,l)*anorm(1,nt)/(2.0*n*n)
      dum2=abs(dum1)
      dum3=exp(test)
      addnve=0.0
      if (dum2 .ge. 1.0) addnve=dum1*dum3
      if ((dum2 .lt. 1.0) .and. (dum3 .ge. 1.0e-10)) addnve=dum1*dum3
      excs=excs+addnve
      if(ran.lt.excs) goto 100
   31 continue
      goto 80
c** for p as target
   40 l=0
      do 41 np1=1,20
      np=np1-1
      nmm1=np1-1
      if(nmm1.le.1) nmm1=1
      npp1=np1+1
      do 41 nm1=nmm1,npp1
      nm=nm1-1
      do 41 nz1=1,20
      nz=nz1-1
      l=l+1
      if(l.gt.1200) goto 41
      nt=np+nm+nz
      if(nt.le.0.or.nt.gt.60) goto 41
      test=-(pi/4.0)*(nt/n)**2
      if (test .gt. expxu) test=expxu
      if (test .lt. expxl) test=expxl
      dum1=anpn*pi*nt*pmul(2,l)*anorm(2,nt)/(2.0*n*n)
      dum2=abs(dum1)
      dum3=exp(test)
      addnve=0.0
      if (dum2 .ge. 1.0) addnve=dum1*dum3
      if ((dum2 .lt. 1.0) .and. (dum3 .ge. 1.0e-10)) addnve=dum1*dum3
      excs=excs+addnve
      if(ran.lt.excs) goto 100
   41 continue
      goto 80
   50 if(nprt(4))
     *write(newbcd,1003) eab,n,nfl,np,nm,nz
      call stpair
      if(int.eq.1) call twob(16,nfl,n)
      if(int.eq.2) call genxpt(16,nfl,n)
      go to 9999
   55 if(nprt(4))
     *write(newbcd,1001)
      goto 53
c** exclusive reaction not found
   80 if(nprt(4))
     *write(newbcd,1004) rs,n
   53 int=1
      np=0
      nm=0
      nz=0
  100 do 101 i=1,60
  101 ipa(i)=0
      if(int.le.0) goto 131
      nprot=1-np+nm+(1-nfl)
      nneut=2-nprot
      goto (112,102),nfl
  102 goto (103,104),int
  103 ipa(1)=16
      ipa(2)=16
      nt=2
      goto 130
  104 if(nneut.eq.1) goto 105
      if(nneut.eq.2) goto 106
      ipa(1)=14
      ipa(2)=14
      goto 120
  105 ipa(1)=14
      ipa(2)=16
      call grndm(rndm,1)
      if(rndm(1).lt.0.5) goto 120
      ipa(1)=16
      ipa(2)=14
      goto 120
  106 ipa(1)=16
      ipa(2)=16
      goto 120
  112 goto (113,114),int
  113 ipa(1)=16
      ipa(2)=14
      nt=2
      if(ncech.eq.0) goto 130
      ipa(1)=14
      ipa(2)=16
      goto 130
  114 if(nneut.eq.1) goto 115
      if(nneut.eq.2) goto 116
      ipa(1)=14
      ipa(2)=14
      goto 120
  115 ipa(1)=14
      ipa(2)=16
      call grndm(rndm,1)
      if(rndm(1).lt.0.33) goto 120
      ipa(1)=16
      ipa(2)=14
      goto 120
  116 ipa(1)=16
      ipa(2)=16
  120 nt=2
      if(np.eq.0) goto 122
      do 121 i=1,np
      nt=nt+1
  121 ipa(nt)=7
  122 if(nm.eq.0) goto 124
      do 123 i=1,nm
      nt=nt+1
  123 ipa(nt)=9
  124 if(nz.eq.0) goto 130
      do 125 i=1,nz
      nt=nt+1
  125 ipa(nt)=8
  130 if(nprt(4))
     *write(newbcd,2004) nt,(ipa(i),i=1,20)
      goto 50
  131 if(nprt(4))
     *write(newbcd,2005)
c
1001  format('0*casn* cascade energetically not possible nuclear',
     * ' excitation',2x,f8.4,2x,'incident energy lost')
1003  format(' *casn* neutron-induced cascade,',
     $ ' avail. energy',2x,f8.4,
     $ 2x,'<ntot>',2x,f8.4,2x,'from',4(2x,i3),2x,'particles')
1004  format(' *casn* neutron-induced cascade,',
     $ ' exclusive reaction not found',
     $ ' try elastic scattering  avail. energy',2x,f8.4,2x,
     $ '<ntot>',2x,f8.4)
2001  format('0*casn* tables for mult. data neutron induced reaction',
     $ ' for definition of numbers see fortran coding')
2002  format(' *casn* target particle flag',2x,i5)
2003  format(1h ,10e12.4)
2004  format(' *casn* ',i3,2x,'particles , mass index array',2x,20i4)
2005  format(' *casn* no particles produced')
c
 9999 continue
      return
      end