*cmz :  3.14/16 13/03/89  14.48.44  by  nick van eijndhoven (cern)
*-- author :
      subroutine caspb(k,int,nfl)
c
c *** cascade of anti proton ***
c *** nve 04-may-1988 cern geneva ***
c
c origin : h.fesefeldt (13-sep-1987)
c
c pb  undergoes interaction with nucleon within nucleus.
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
c --- initialization flags for various gheisha routines ---
      common /kginit/ kginit(50)
c
c --- boundary limits for arguments of intrinsic functions ---
c --- xl denotes lower bound whereas xu denotes upper bound ---
      common /limits/ expxl,expxu
c
c
      real n
      dimension pmul1(2,1200),pmul2(2,400),anorm1(2,60),anorm2(2,60),
     $          supp(10),cech(20),anhl(29),b(2)
      dimension rndm(1)
      save pmul1,anorm1,pmul2,anorm2
      data supp/0.,0.4,0.55,0.65,0.75,0.82,0.86,0.90,0.94,0.98/
      data cech/0.14,0.17,0.18,0.18,0.18,0.17,0.17,0.16,0.155,0.145,
     *          0.11,0.082,0.065,0.050,0.041,0.035,0.028,0.024,0.010
     *         ,0.0/
      data anhl/1.00,1.00,1.00,1.00,1.0,1.00,1.0,1.00,1.00,0.90
     *         ,0.6,0.52,0.47,0.44,0.41,0.39,0.37,0.35,0.34,0.24
     *         ,0.19,0.15,0.12,0.10,0.09,0.07,0.06,0.05,0./
      data b/0.7,0.7/,c/1.25/
c
c --- initialization indicated by kginit(11) ---
      if (kginit(11) .ne. 0) go to 10
      kginit(11)=1
c
c --- initialize pmul and anorm arrays ---
      do 9000 j=1,1200
      do 9001 i=1,2
      pmul1(i,j)=0.0
      if (j .le. 400) pmul2(i,j)=0.0
      if (j .le. 60) anorm1(i,j)=0.0
      if (j .le. 60) anorm2(i,j)=0.0
 9001 continue
 9000 continue
c
c** compute normalization constants
c** for p as target
c
      l=0
      do 1 np1=1,20
      np=np1-1
      nmm1=np1-1
      if(nmm1.le.1) nmm1=1
      npp1=np1+1
      do 1 nm1=nmm1,npp1
      nm=nm1-1
      do 1 nz1=1,20
      nz=nz1-1
      l=l+1
      if(l.gt.1200) goto 1
      nt=np+nm+nz
      if(nt.le.0.or.nt.gt.60) goto 1
      pmul1(1,l)=pmltpc(np,nm,nz,nt,b(1),c)
      anorm1(1,nt)=anorm1(1,nt)+pmul1(1,l)
    1 continue
c** for n as target
      l=0
      do 2 np1=1,20
      np=np1-1
      npp1=np1+2
      do 2 nm1=np1,npp1
      nm=nm1-1
      do 2 nz1=1,20
      nz=nz1-1
      l=l+1
      if(l.gt.1200) goto 2
      nt=np+nm+nz
      if(nt.le.0.or.nt.gt.60) goto 2
      pmul1(2,l)=pmltpc(np,nm,nz,nt,b(2),c)
      anorm1(2,nt)=anorm1(2,nt)+pmul1(2,l)
    2 continue
      do 3 i=1,60
      if(anorm1(1,i).gt.0.) anorm1(1,i)=1./anorm1(1,i)
      if(anorm1(2,i).gt.0.) anorm1(2,i)=1./anorm1(2,i)
    3 continue
      if(.not.nprt(10)) goto 9
      write(newbcd,2001)
      do 4 nfl=1,2
      write(newbcd,2002) nfl
      write(newbcd,2003) (anorm1(nfl,i),i=1,60)
      write(newbcd,2003) (pmul1(nfl,i),i=1,1200)
    4 continue
c** do the same for annihilation channels
c** for p as target
c
    9 l=0
      do 5 np1=1,20
      np=np1-1
      nm=np
      do 5 nz1=1,20
      nz=nz1-1
      l=l+1
      if(l.gt.400) goto 5
      nt=np+nm+nz
      if(nt.le.1.or.nt.gt.60) goto 5
      pmul2(1,l)=pmltpc(np,nm,nz,nt,b(1),c)
      anorm2(1,nt)=anorm2(1,nt)+pmul2(1,l)
    5 continue
c** for n as target
      l=0
      do 6 np1=1,20
      np=np1-1
      nm=np+1
      do 6 nz1=1,20
      nz=nz1-1
      l=l+1
      if(l.gt.400) goto 6
      nt=np+nm+nz
      if(nt.le.1.or.nt.gt.60) goto 6
      pmul2(2,l)=pmltpc(np,nm,nz,nt,b(2),c)
      anorm2(2,nt)=anorm2(2,nt)+pmul2(2,l)
    6 continue
      do 7 i=1,60
      if(anorm2(1,i).gt.0.) anorm2(1,i)=1./anorm2(1,i)
      if(anorm2(2,i).gt.0.) anorm2(2,i)=1./anorm2(2,i)
    7 continue
      if(.not.nprt(10)) goto 10
      write(newbcd,3001)
      do 8 nfl=1,2
      write(newbcd,3002) nfl
      write(newbcd,3003) (anorm2(nfl,i),i=1,60)
      write(newbcd,3003) (pmul2(nfl,i),i=1,400)
    8 continue
c** choose proton or neutron as target
   10 nfl=2
      call grndm(rndm,1)
      if(rndm(1).lt.zno2/atno2) nfl=1
      tarmas=rmass(14)
      if (nfl .eq. 2) tarmas=rmass(16)
      s=amasq+tarmas**2+2.0*tarmas*en
      rs=sqrt(s)
      enp(8)=amasq+tarmas**2+2.0*tarmas*enp(6)
      enp(9)=sqrt(enp(8))
      eab=rs-tarmas-abs(rmass(15))
c**  elastic scattering
      ncech=0
      np=0
      nm=0
      nz=0
      n=0.
      if(int.eq.2) goto 20
c** introduce charge exchange reaction pb p --> nb n
      if(nfl.eq.2) goto 100
      iplab=ifix(p*10.)+1
      if(iplab.gt.10) iplab=ifix(p)+10
      if(iplab.gt.20) iplab=20
      call grndm(rndm,1)
      if(rndm(1).gt.cech(iplab)/atno2**0.75) goto 100
      ncech=1
      goto 100
c** annihilation channels
   20 iplab=ifix(p*10.)+1
      if(iplab.gt.10) iplab=ifix(p)+10
      if(iplab.gt.19) iplab=ifix(p/10.)+19
      if(iplab.gt.28) iplab=29
      call grndm(rndm,1)
      if(rndm(1).gt.anhl(iplab)) goto 19
      eab=rs
      if (eab .le. 2.0*rmass(7)) goto 55
      goto 222
c**  check if energetically possible to produce one extra pion in react.
   19 if (eab .le. rmass(7)) goto 55
c**  suppression of high multiplicity events at low momentum
      ieab=ifix(eab*5.)+1
      if(ieab.gt.10) goto 22
      call grndm(rndm,1)
      if(rndm(1).lt.supp(ieab)) goto 22
      n=1.
      goto (24,23),nfl
 23   continue
      test=-(1+b(1))**2/(2.0*c**2)
      if (test .lt. expxl) test=expxl
      if (test .gt. expxu) test=expxu
      w0=exp(test)
      test=-(-1+b(1))**2/(2.0*c**2)
      if (test .lt. expxl) test=expxl
      if (test .gt. expxu) test=expxu
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
      test=-(1+b(2))**2/(2.0*c**2)
      if (test .lt. expxl) test=expxl
      if (test .gt. expxu) test=expxu
      w0=exp(test)
      wp=exp(test)
      test=-(-1+b(2))**2/(2.0*c**2)
      if (test .lt. expxl) test=expxl
      if (test .gt. expxu) test=expxu
      wm=exp(test)
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
      if (test .lt. expxl) test=expxl
      if (test .gt. expxu) test=expxu
      dum1=pi*nt/(2.0*n*n)
      dum2=abs(dum1)
      dum3=exp(test)
      addnve=0.0
      if (dum2 .ge. 1.0) addnve=dum1*dum3
      if ((dum2 .lt. 1.0) .and. (dum3 .gt. 1.0e-10)) addnve=dum1*dum3
      anpn=anpn+addnve
   21 continue
      anpn=1./anpn
c** p or n as target
      call grndm(rndm,1)
      ran=rndm(1)
      excs=0.
      goto (30,40),nfl
c** for p as target
   30 l=0
      do 31 np1=1,20
      np=np1-1
      nmm1=np1-1
      if(nmm1.le.1) nmm1=1
      npp1=np1+1
      do 31 nm1=nmm1,npp1
      nm=nm1-1
      do 31 nz1=1,20
      nz=nz1-1
      l=l+1
      if(l.gt.1200) goto 31
      nt=np+nm+nz
      if(nt.le.0.or.nt.gt.60) goto 31
      test=-(pi/4.0)*(nt/n)**2
      if (test .lt. expxl) test=expxl
      if (test .gt. expxu) test=expxu
      dum1=anpn*pi*nt*pmul1(1,l)*anorm1(1,nt)/(2.0*n*n)
      dum2=abs(dum1)
      dum3=exp(test)
      addnve=0.0
      if (dum2 .ge. 1.0) addnve=dum1*dum3
      if ((dum2 .lt. 1.0) .and. (dum3 .gt. 1.0e-10)) addnve=dum1*dum3
      excs=excs+addnve
      if(ran.lt.excs) goto 100
   31 continue
      goto 80
c** for n as target
   40 l=0
      do 41 np1=1,20
      np=np1-1
      npp1=np1+2
      do 41 nm1=np1,npp1
      nm=nm1-1
      do 41 nz1=1,20
      nz=nz1-1
      l=l+1
      if(l.gt.1200) goto 41
      nt=np+nm+nz
      if(nt.le.0.or.nt.gt.60) goto 41
      test=-(pi/4.0)*(nt/n)**2
      if (test .lt. expxl) test=expxl
      if (test .gt. expxu) test=expxu
      dum1=anpn*pi*nt*pmul1(2,l)*anorm1(2,nt)/(2.0*n*n)
      dum2=abs(dum1)
      dum3=exp(test)
      addnve=0.0
      if (dum2 .ge. 1.0) addnve=dum1*dum3
      if ((dum2 .lt. 1.0) .and. (dum3 .gt. 1.0e-10)) addnve=dum1*dum3
      excs=excs+addnve
      if(ran.lt.excs) goto 100
   41 continue
      goto 80
c** annihilation channels
  222 ipa(1)=0
      ipa(2)=0
      aleab=log(eab)
c** no. of total particles vs sqrt(s)
      n=3.62567+0.665843*aleab+0.336514*aleab*aleab
     * +0.117712*aleab*aleab*aleab+0.0136912*aleab*aleab*aleab*aleab
      n=n-2.
c** normalization constant for  kno-distribution
      anpn=0.
      do 221 nt=2,60
      test=-(pi/4.0)*(nt/n)**2
      if (test .lt. expxl) test=expxl
      if (test .gt. expxu) test=expxu
      dum1=pi*nt/(2.0*n*n)
      dum2=abs(dum1)
      dum3=exp(test)
      addnve=0.0
      if (dum2 .ge. 1.0) addnve=dum1*dum3
      if ((dum2 .lt. 1.0) .and. (dum3 .gt. 1.0e-10)) addnve=dum1*dum3
      anpn=anpn+addnve
  221 continue
      anpn=1./anpn
c** p or n as target
      call grndm(rndm,1)
      ran=rndm(1)
      excs=0.
      goto (230,240),nfl
c** for p as target
  230 l=0
      do 231 np1=1,20
      np=np1-1
      nm=np
      do 231 nz1=1,20
      nz=nz1-1
      l=l+1
      if(l.gt.400) goto 231
      nt=np+nm+nz
      if(nt.le.1.or.nt.gt.60) goto 231
      test=-(pi/4.0)*(nt/n)**2
      if (test .lt. expxl) test=expxl
      if (test .gt. expxu) test=expxu
      dum1=anpn*pi*nt*pmul2(1,l)*anorm2(1,nt)/(2.0*n*n)
      dum2=abs(dum1)
      dum3=exp(test)
      addnve=0.0
      if (dum2 .ge. 1.0) addnve=dum1*dum3
      if ((dum2 .lt. 1.0) .and. (dum3 .gt. 1.0e-10)) addnve=dum1*dum3
      excs=excs+addnve
      if(ran.lt.excs) goto 120
  231 continue
      goto 80
c** for n as target
  240 l=0
      do 241 np1=1,20
      np=np1-1
      nm=np+1
      do 241 nz1=1,20
      nz=nz1-1
      l=l+1
      if(l.gt.400) goto 241
      nt=np+nm+nz
      if(nt.le.1.or.nt.gt.60) goto 241
      test=-(pi/4.0)*(nt/n)**2
      if (test .lt. expxl) test=expxl
      if (test .gt. expxu) test=expxu
      dum1=anpn*pi*nt*pmul2(2,l)*anorm2(2,nt)/(2.0*n*n)
      dum2=abs(dum1)
      dum3=exp(test)
      addnve=0.0
      if (dum2 .ge. 1.0) addnve=dum1*dum3
      if ((dum2 .lt. 1.0) .and. (dum3 .gt. 1.0e-10)) addnve=dum1*dum3
      excs=excs+addnve
      if(ran.lt.excs) goto 120
  241 continue
      goto 80
   50 if(nprt(4))
     *write(newbcd,1003) eab,n,nfl,np,nm,nz
      call stpair
      if(int.eq.1) call twob(15,nfl,n)
      if(int.eq.2) call genxpt(15,nfl,n)
      go to 9999
   55 if(nprt(4))
     *write(newbcd,1001)
      goto 53
c** exclusive reaction not found,assume elastic scattering
   80 if(nprt(4))
     *write(newbcd,1004)eab,n
   53 int=1
      np=0
      nm=0
      nz=0
  100 do 101 i=1,60
  101 ipa(i)=0
      if(int.le.0) goto 131
      goto (112,102),nfl
  102 goto (103,104),int
  103 ipa(1)=15
      ipa(2)=16
      nt=2
      goto 130
  104 if(np.eq.-1+nm) goto 105
      if(np.eq.   nm) goto 106
      ipa(1)=17
      ipa(2)=14
      goto 120
  105 ipa(1)=15
      ipa(2)=14
      call grndm(rndm,1)
      if(rndm(1).lt.0.5) goto 120
      ipa(1)=17
      ipa(2)=16
      goto 120
  106 ipa(1)=15
      ipa(2)=16
      goto 120
  112 goto (113,114),int
  113 ipa(1)=15
      ipa(2)=14
      nt=2
      if(ncech.eq.0) goto 130
      ipa(1)=17
      ipa(2)=16
      goto 130
  114 if(np.eq.  nm) goto 115
      if(np.eq.1+nm) goto 116
      ipa(1)=17
      ipa(2)=14
      goto 120
  115 ipa(1)=17
      ipa(2)=16
      call grndm(rndm,1)
      if(rndm(1).lt.0.33) goto 120
      ipa(1)=15
      ipa(2)=14
      goto 120
  116 ipa(1)=15
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
1001  format('0*caspb* cascade energetically not possible',
     $ ' continue with quasi-elastic scattering')
1003  format(' *caspb* antiproton-induced cascade,',
     $ ' avail. energy',2x,f8.4,
     $ 2x,'<ntot>',2x,f8.4,2x,'from',4(2x,i3),2x,'particles')
1004  format(' *caspb* antiproton-induced cascade,',
     $ ' exclusive reaction',
     $ ' not found  try elastic scattering  avail. energy',2x,f8.4,2x,
     $ ' <ntot>',2x,f8.4)
2001  format('0*caspb* tables for mult. data antiproton induced ',
     $ 'reaction  for definition of numbers see fortran coding')
2002  format(' *caspb* target particle flag',2x,i5)
2003  format(1h ,10e12.4)
2004  format(' *caspb* ',i3,2x,'particles , mass index array',2x,20i4)
2005  format(' *caspb* no particles produced')
3001  format('0*caspb* tables for mult. data antiproton induced ',
     $ ' annihilation reaction  for definition of numbers see fortran',
     $ ' coding')
3002  format(' *caspb* target particle flag',2x,i5)
3003  format(1h ,10e12.4)
c
 9999 continue
      return
      end
