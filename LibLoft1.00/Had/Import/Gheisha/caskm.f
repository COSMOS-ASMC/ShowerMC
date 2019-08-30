*cmz :  3.14/16 13/03/89  14.48.42  by  nick van eijndhoven (cern)
*-- author :
      subroutine caskm(k,int,nfl)
c
c *** cascade of k- ***
c *** nve 04-may-1988 cern geneva ***
c
c origin : h.fesefeldt (13-sep-1987)
c
c k-  undergoes interaction with nucleon within nucleus.
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
      dimension pmul(2,1200),anorm(2,60),cech(10),cnk0(20),piy1(4),
     $          piy2(3),ipiy1(2,4),ipiy2(2,3),ipiy3(2,3),b(2)
      dimension rndm(1)
      save pmul,anorm
      data cech/1.,1.,1.,0.70,0.60,0.55,0.35,0.25,0.18,0.15/
      data cnk0/0.17,0.18,0.17,0.24,0.26,0.20,0.22,0.21,0.34,0.45
     $         ,0.58,0.55,0.36,0.29,0.29,0.32,0.32,0.33,0.33,0.33/
      data piy1/0.67,0.78,0.89,1.00/,piy2/0.68,0.84,1.00/
      data ipiy1/8,18,9,20,8,21,7,22/
      data ipiy2/9,18,9,21,8,22/,ipiy3/7,18,8,20,7,21/
      data b/0.7,0.7/,c/1.25/
c
c --- initialization indicated by kginit(4) ---
      if (kginit(4) .ne. 0) go to 10
      kginit(4)=1
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
      pmul(1,l)=pmltpc(np,nm,nz,nt,b(1),c)
      anorm(1,nt)=anorm(1,nt)+pmul(1,l)
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
      pmul(2,l)=pmltpc(np,nm,nz,nt,b(2),c)
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
      eab=rs-tarmas-rmass(13)
c
c**  elastic scattering
      np=0
      nm=0
      nz=0
      n=0.
      ipa(1)=13
      ipa(2)=14
      ncech=0
      if(nfl.eq.2) ipa(2)=16
      if(int.eq.2) goto 20
      goto 100
c**  check if energetically possible to produce one extra pion in react.
   20 iplab=ifix(p*5.)+1
      if(iplab.gt.10) goto 22
      call grndm(rndm,1)
      if(rndm(1).lt.cech(iplab)) goto 19
      if (eab .lt. rmass(7)) goto 55
      goto 22
c** charge exchange reaction (is included in inelastic cross section)
   19 iplab=ifix(p*10.)+1
      if(iplab.gt.20) iplab=20
      call grndm(rndm,1)
      if(rndm(1).gt.cnk0(iplab)) goto 24
      if(nfl.eq.1) goto 23
c** for k- n reaction no k n strangeness exchange possible
      int=1
      ipa(1)=13
      ipa(2)=16
      goto 100
   23 int=1
      ipa(1)=12
      ipa(2)=16
      goto 100
c** p l, p s reactions
   24 call grndm(rndm,1)
      ran=rndm(1)
      if(ran.lt.0.25) goto 25
      if(ran.lt.0.50) goto 26
      if(ran.lt.0.75) goto 27
c** k- p --> pi0 l or k- n --> pi- l
      ipa(1)=8
      if(nfl.eq.2) ipa(1)=9
      ipa(2)=18
      goto 100
c** k- p --> pi- s+
   25 ipa(1)=9
      ipa(2)=20
      if(nfl.eq.1) goto 100
      ipa(1)=13
      ipa(2)=16
      goto 100
c** k- p --> pi0 s0  or k- n --> pi- s0
   26 ipa(1)=8
      if(nfl.eq.2) ipa(1)=9
      ipa(2)=21
      goto 100
c** k- p --> pi+ s-  or k- n --> pi0 s-
   27 ipa(1)=7
      if(nfl.eq.2) ipa(1)=8
      ipa(2)=22
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
      dum1=anpn*pi*nt*pmul(1,l)*anorm(1,nt)/(2.0*n*n)
      dum2=abs(dum1)
      dum3=exp(test)
      addnve=0.0
      if (dum2 .ge. 1.0) addnve=dum1*dum3
      if ((dum2 .lt. 1.0) .and. (dum3 .ge. 1.0e-10)) addnve=dum1*dum3
      excs=excs+addnve
      if(ran.lt.excs) goto 50
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
      dum1=anpn*pi*nt*pmul(2,l)*anorm(2,nt)/(2.0*n*n)
      dum2=abs(dum1)
      dum3=exp(test)
      addnve=0.0
      if (dum2 .ge. 1.0) addnve=dum1*dum3
      if ((dum2 .lt. 1.0) .and. (dum3 .ge. 1.0e-10)) addnve=dum1*dum3
      excs=excs+addnve
      if(ran.lt.excs) goto 50
   41 continue
      goto 80
   50 goto (60,65),nfl
   60 if(np.eq.nm) goto 61
      if(np.eq.1+nm) goto 63
      ipa(1)=12
      ipa(2)=14
      goto 90
   61 call grndm(rndm,1)
      if(rndm(1).lt.0.75) goto 62
      ipa(1)=12
      ipa(2)=16
      goto 90
   62 ipa(1)=13
      ipa(2)=14
      goto 90
   63 ipa(1)=13
      ipa(2)=16
      goto 90
   65 if(np.eq.-1+nm) goto 66
      if(np.eq.nm) goto 68
      ipa(1)=12
      ipa(2)=16
      goto 90
   66 call grndm(rndm,1)
      if(rndm(1).lt.0.50) goto 67
      ipa(1)=12
      ipa(2)=16
      goto 90
   67 ipa(1)=13
      ipa(2)=14
      goto 90
   68 ipa(1)=13
      ipa(2)=16
c**  pi y production instead of k n
   90 call grndm(rndm,1)
      if(rndm(1).lt.0.5) goto 100
      if(ipa(1).eq.13.and.ipa(2).eq.16) goto 95
      if(ipa(1).eq.11.and.ipa(2).eq.14) goto 95
      if(ipa(1).eq.12.and.ipa(2).eq.14) goto 95
      call grndm(rndm,1)
      ran=rndm(1)
      do 91 i=1,4
      if(ran.lt.piy1(i)) goto 92
   91 continue
      goto 100
   92 ipa(1)=ipiy1(1,i)
      ipa(2)=ipiy1(2,i)
      goto 100
   95 call grndm(rndm,1)
      ran=rndm(1)
      do 96 i=1,3
      if(ran.lt.piy2(i)) goto 97
   96 continue
      goto 100
   97 if(ipa(2).eq.14) goto 98
      ipa(1)=ipiy2(1,i)
      ipa(2)=ipiy2(2,i)
      goto 100
   98 ipa(1)=ipiy3(1,i)
      ipa(2)=ipiy3(2,i)
      goto 100
   70 if(nprt(4))
     *write(newbcd,1003) eab,n,nfl,np,nm,nz
      call stpair
      if(int.eq.1) call twob(13,nfl,n)
      if(int.eq.2) call genxpt(13,nfl,n)
      go to 9999
c** nuclear excitation
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
      n=0.
      ipa(1)=13
      ipa(2)=14
      if(nfl.eq.2) ipa(2)=16
  100 do 101 i=3,60
  101 ipa(i)=0
      if(int.le.0) goto 131
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
      do 132 i=1,nt
      if(ipa(i).ne.12) goto 132
      call grndm(rndm,1)
      if(rndm(1).lt.0.5) goto 132
      ipa(i)=11
  132 continue
      goto 70
  131 if(nprt(4))
     *write(newbcd,2005)
c
1001  format('0*caskm* cascade energetically not possible',
     $ ' continue with quasi-elastic scattering')
1003  format(' *caskm* kaon- -induced cascade,',
     $ ' avail. energy',2x,f8.4,
     $ 2x,'<ntot>',2x,f8.4,2x,'from',4(2x,i3),2x,'particles')
1004  format(' *caskm* kaon- -induced cascade,',
     $ ' exclusive reaction not found',
     $ ' try elastic scattering  avail. energy',2x,f8.4,2x,
     $ '<ntot>',2x,f8.4)
2001  format('0*caskm* tables for mult. data kaon-  induced reaction',
     $ ' for definition of numbers see fortran coding')
2002  format(' *caskm* target particle flag',2x,i5)
2003  format(1h ,10e12.4)
2004  format(' *caskm* ',i3,2x,'particles , mass index array',2x,20i4)
2005  format(' *caskm* no particles produced')
c
 9999 continue
      return
      end
