*cmz :  3.14/16 13/03/89  14.48.45  by  nick van eijndhoven (cern)
*-- author :
      subroutine casasm(k,int,nfl)
c
c *** cascade of anti sigma- ***
c *** nve 04-may-1988 cern geneva ***
c
c origin : h.fesefeldt (13-sep-1987)
c
c s-b undergoes interaction with nucleon within nucleus.
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
     $          cech(10),anhl(25),iipa(10,2),b(2)
      dimension rndm(2)
      save pmul1,anorm1,pmul2,anorm2
      data cech/0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0./
      data anhl/1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,0.97,0.88
     $         ,0.85,0.81,0.75,0.64,0.64,0.55,0.55,0.45,0.47,0.40
     $         ,0.39,0.36,0.33,0.10,0.01/
      data iipa/14,14,14,14,14,19,24,16,14,14,
     $          25,25,25,25,25,14,14,25,19,24/
      data b/0.7,0.7/,c/1.25/
c
c --- initialization indicated by kginit(2) ---
      if (kginit(2) .ne. 0) go to 10
      kginit(2)=1
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
      nmm1=np1-2
      if(nmm1.le.1) nmm1=1
      npp1=np1
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
      nmm1=np1-1
      if(nmm1.lt.1) nmm1=1
      npp1=np1+1
      do 2 nm1=nmm1,npp1
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
      do 5 np1=3,20
      np=np1-1
      nm=np-2
      do 5 nz1=1,20
      nz=nz1-1
      if(nm+nz.eq.0) goto 5
      l=l+1
      if(l.gt.400) goto 5
      nt=np+nm+nz
      if(nt.le.2.or.nt.gt.60) goto 5
      pmul2(1,l)=pmltpc(np,nm,nz,nt,b(1),c)
      anorm2(1,nt)=anorm2(1,nt)+pmul2(1,l)
    5 continue
c** for n as target
      l=0
      do 6 np1=2,20
      np=np1-1
      nm=np-1
      do 6 nz1=1,20
      nz=nz1-1
      if(nm+nz.eq.0) goto 6
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
      eab=rs-tarmas-abs(rmass(25))
c**  elastic scattering
      ncech=0
      np=0
      nm=0
      nz=0
      ipa(1)=25
      ipa(2)=14
      if(nfl.eq.2) ipa(2)=16
      n=0.
      if(int.eq.2) goto 20
      iplab=ifix(p*2.5)+1
      if(iplab.gt.10) iplab=10
      call grndm(rndm,1)
      if(rndm(1).gt.cech(iplab)/atno2**0.42) goto 120
c** introduce charge and strangeness exchange reaction
c** s-b n --> lb p, s-b n --> s0b p,
c**                                  s-b p --> p s-b
c** s-b n --> n s-b, s-b n --> p lb, s-b n --> p s0b
      call grndm(rndm,1)
      ran=rndm(1)
      irn=ifix(ran/0.2)+1
      if(irn.gt.5) irn=5
      irn=irn+(nfl-1)*5
      ipa(1)=iipa(irn,1)
      ipa(2)=iipa(irn,2)
      ncech=1
      goto 120
c** annihilation channels
   20 iplab=ifix(p*10.)+1
      if(iplab.gt.10) iplab=ifix((p-1.)*5.)+11
      if(iplab.gt.15) iplab=ifix( p-2.    )+16
      if(iplab.gt.23) iplab=ifix((p-10.)/10.)+24
      if(iplab.gt.25) iplab=25
      call grndm(rndm,1)
      if(rndm(1).gt.anhl(iplab)) goto 19
      eab=rs
      if (eab .le. rmass(7)+rmass(10)) goto 55
      goto 222
c**  check if energetically possible to produce one extra pion in react.
   19 if (eab .le. rmass(7)) goto 55
      aleab=log(eab)
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
      if ((dum2 .lt. 1.0) .and. (dum3 .ge. 1.0e-10)) addnve=dum1*dum3
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
      nmm1=np1-2
      if(nmm1.le.1) nmm1=1
      npp1=np1
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
      if ((dum2 .lt. 1.0) .and. (dum3 .ge. 1.0e-10)) addnve=dum1*dum3
      excs=excs+addnve
      if(ran.lt.excs) goto 100
   31 continue
      goto 80
c** for n as target
   40 l=0
      do 41 np1=1,20
      np=np1-1
      nmm1=np1-1
      if(nmm1.lt.1) nmm1=1
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
      if (test .lt. expxl) test=expxl
      if (test .gt. expxu) test=expxu
      dum1=anpn*pi*nt*pmul1(2,l)*anorm1(2,nt)/(2.0*n*n)
      dum2=abs(dum1)
      dum3=exp(test)
      addnve=0.0
      if (dum2 .ge. 1.0) addnve=dum1*dum3
      if ((dum2 .lt. 1.0) .and. (dum3 .ge. 1.0e-10)) addnve=dum1*dum3
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
      if ((dum2 .lt. 1.0) .and. (dum3 .ge. 1.0e-10)) addnve=dum1*dum3
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
      do 231 np1=3,20
      np=np1-1
      nm=np-2
      do 231 nz1=1,20
      nz=nz1-1
      if(nm+nz.eq.0) goto 231
      l=l+1
      if(l.gt.400) goto 231
      nt=np+nm+nz
      if(nt.le.2.or.nt.gt.60) goto 231
      test=-(pi/4.0)*(nt/n)**2
      if (test .lt. expxl) test=expxl
      if (test .gt. expxu) test=expxu
      dum1=anpn*pi*nt*pmul2(1,l)*anorm2(1,nt)/(2.0*n*n)
      dum2=abs(dum1)
      dum3=exp(test)
      addnve=0.0
      if (dum2 .ge. 1.0) addnve=dum1*dum3
      if ((dum2 .lt. 1.0) .and. (dum3 .ge. 1.0e-10)) addnve=dum1*dum3
      excs=excs+addnve
      if(ran.lt.excs) goto 120
  231 continue
      goto 80
c** for n as target
  240 l=0
      do 241 np1=2,20
      np=np1-1
      nm=np-1
      do 241 nz1=1,20
      nz=nz1-1
      if(nm+nz.eq.0) goto 241
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
      if ((dum2 .lt. 1.0) .and. (dum3 .ge. 1.0e-10)) addnve=dum1*dum3
      excs=excs+addnve
      if(ran.lt.excs) goto 120
  241 continue
      goto 80
   50 if(nprt(4))
     *write(newbcd,1003) eab,n,nfl,np,nm,nz
      if(int.eq.1) call twob(25,nfl,n)
      if(int.eq.2) call genxpt(25,nfl,n)
      go to 9999
   55 if(nprt(4))
     *write(newbcd,1001)
      goto 53
c** exclusive reaction not found,assume elastic scattering
   80 if(nprt(4))
     *write(newbcd,1004) rs,n
   53 int=1
      np=0
      nm=0
      nz=0
      ipa(1)=25
      ipa(2)=14
      if(nfl.eq.2) ipa(2)=16
      goto 120
  100 do 101 i=1,60
  101 ipa(i)=0
      if(int.le.0) goto 131
      goto (102,112),nfl
  102 ncht=np-nm
      ncht=ncht+1
      if(ncht.le.0) ncht=1
      if(ncht.gt.3) ncht=3
      goto(103,104,105),ncht
  103 ipa(1)=25
      ipa(2)=14
      goto 120
  104 ipa(1)=25
      ipa(2)=16
      call grndm(rndm,2)
      if(rndm(1).lt.0.5) goto 120
      ipa(1)=24
      if(rndm(2).lt.0.5) ipa(1)=19
      ipa(2)=14
      goto 120
  105 ipa(1)=24
      call grndm(rndm,1)
      if(rndm(1).lt.0.5) ipa(1)=19
      ipa(2)=16
      goto 120
  112 ncht=np-nm
      ncht=ncht+2
      if(ncht.le.0) ncht=1
      if(ncht.gt.3) ncht=3
      goto(113,114,115),ncht
  113 ipa(1)=25
      ipa(2)=14
      goto 120
  114 ipa(1)=24
      call grndm(rndm,2)
      if(rndm(1).lt.0.5) ipa(1)=19
      ipa(2)=14
      if(rndm(2).lt.0.5) goto 120
      ipa(1)=25
      ipa(2)=16
      goto 120
  115 ipa(1)=24
      call grndm(rndm,1)
      if(rndm(1).lt.0.5) ipa(1)=19
      ipa(2)=16
  120 nt=2
      if(ipa(1).ne.0) goto 119
      if(nz.eq.0) goto 118
      if(nm.eq.0) goto 117
      call grndm(rndm,1)
      if(rndm(1).lt.0.5) goto 118
  117 ipa(3)=12
      nz=nz-1
      nt=3
      goto 119
  118 if(nm.eq.0) goto 119
      ipa(3)=13
      nm=nm-1
      nt=3
  119 if(np.eq.0) goto 122
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
1001  format('0*casasm* cascade energetically not possible',
     $ ' continue with quasi-elastic scattering')
1003  format(' *casasm* antisigma- -induced cascade',
     $ ' avail. energy',2x,f8.4,
     $ 2x,'<ntot>',2x,f8.4,2x,'from',4(2x,i3),2x,'particles')
1004  format(' *casasm* antisigma- -induced cascade,',
     $ ' exclusive reaction',
     $' not found  try elastic scattering  avail. energy',2x,f8.4,2x,
     $ '<ntot>',2x,f8.4)
2001  format('0*casasm* tables for mult. data antisigma- induced ',
     $'reaction  for definition of numbers see fortran coding')
2002  format(' *casasm* target particle flag',2x,i5)
2003  format(1h ,10e12.4)
2004  format(' *casasm* ',i3,2x,'particles , mass index array',2x,20i4)
2005  format(' *casasm* no particles produced')
3001  format('0*casasm* tables for mult. data antisigma- induced ',
     *'annihilation reaction  for definition of numbers see fortran',
     *' coding')
3002  format(' *casasm* target particle flag',2x,i5)
3003  format(1h ,10e12.4)
c
 9999 continue
      return
      end
