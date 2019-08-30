*cmz :  3.14/16 13/03/89  14.48.47  by  nick van eijndhoven (cern)
*-- author :    nick van eijndhoven (cern)   02/02/89
      subroutine casx0(k,int,nfl)
c
c *** cascade of xi0 ***
c *** nve 20-jan-1989 cern geneva ***
c
c xi0  undergoes interaction with nucleon within nucleus.
c check if energetically possible to produce pions/kaons.
c if not, assume nuclear excitation occurs, degrade input particle
c in energy and no other particles are produced.
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
      dimension pmul(2,1200),anorm(2,60),cech(10),iipa(12,2),b(2)
      dimension rndm(1)
      save pmul,anorm
      data cech/0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0./
c --- array iipa denotes the strangeness and charge exchage reactions ---
c xi0 p --> s+ s0,  xi0 p --> s0 s+
c xi0 p --> s+ l0,  xi0 p --> l0 s+
c xi0 p --> p xi0
c xi0 n --> s0 s0
c xi0 n --> l0 l0
c xi0 n --> xi- p,  xi0 n --> p xi-
c xi0 n --> s+ s-,  xi0 n --> s- s+
c xi0 n --> n xi0
      data iipa/20,21,20,18,14, 21,18,27,14,20,22,16,
     *          21,20,18,20,26, 21,18,14,27,22,20,26/
      data b/0.7,0.7/,c/1.25/
c
c --- initialization indicated by kginit(20) ---
      if (kginit(20) .ne. 0) go to 10
      kginit(20)=1
c
c --- initialize pmul and anorm arrays ---
      do 9000 j=1,1200
      do 9001 i=1,2
      pmul(i,j)=0.0
      if (j .le. 60) anorm(i,j)=0.0
 9001 continue
 9000 continue
c
c *** compute normalization constants ***
c
c --- for p target ---
      l=0
      do 1 np1=1,20
      np=np1-1
      nmm1=np1-2
      if (nmm1 .le. 0) nmm1=1
      npp1=np1+1
      do 1 nm1=nmm1,npp1
      nm=nm1-1
      do 1 nz1=1,20
      nz=nz1-1
      l=l+1
      if (l .gt. 1200) go to 1
      nt=np+nm+nz
      if ((nt .le. 0) .or. (nt .gt. 60)) go to 1
      pmul(1,l)=pmltpc(np,nm,nz,nt,b(1),c)
      anorm(1,nt)=anorm(1,nt)+pmul(1,l)
 1    continue
c --- for n target ---
      l=0
      do 2 np1=1,20
      np=np1-1
      nmm1=np1-1
      if (nmm1 .le. 0) nmm1=1
      npp1=np1+2
      do 2 nm1=nmm1,npp1
      nm=nm1-1
      do 2 nz1=1,20
      nz=nz1-1
      l=l+1
      if (l .gt. 1200) go to 2
      nt=np+nm+nz
      if ((nt .le. 0) .or. (nt .gt. 60)) go to 2
      pmul(2,l)=pmltpc(np,nm,nz,nt,b(2),c)
      anorm(2,nt)=anorm(2,nt)+pmul(2,l)
 2    continue
c
      do 3 i=1,60
      if (anorm(1,i) .gt. 0.) anorm(1,i)=1./anorm(1,i)
      if (anorm(2,i) .gt. 0.) anorm(2,i)=1./anorm(2,i)
 3    continue
c
c      if (.not. nprt(10)) go to 10
c
      write(newbcd,2001)
 2001 format('0*casx0* tables for mult. data xi0 induced reaction',
     $ ' for definition of numbers see fortran coding')
      do 4 nfl=1,2
      write(newbcd,2002) nfl
 2002 format(' *casx0* target particle flag',2x,i5)
      write(newbcd,2003) (anorm(nfl,i),i=1,60)
      write(newbcd,2003) (pmul(nfl,i),i=1,1200)
 2003 format(1h ,10e12.4)
 4    continue
c
c --- select target nucleon ---
 10   continue
      nfl=2
      call grndm(rndm,1)
      if (rndm(1) .lt. (zno2/atno2)) nfl=1
      tarmas=rmass(14)
      if (nfl .eq. 2) tarmas=rmass(16)
      s=amasq+tarmas**2+2.0*tarmas*en
      rs=sqrt(s)
      enp(8)=amasq+tarmas**2+2.0*tarmas*enp(6)
      enp(9)=sqrt(enp(8))
      eab=rs-tarmas-rmass(26)
c
c --- reset strangeness fixing flag ---
      nvefix=0
c
c *** elastic scattering ***
      np=0
      nm=0
      nz=0
      n=0.
      ipa(1)=26
      ipa(2)=14
      if (nfl .eq. 2) ipa(2)=16
      ncech=0
c
      if (int .eq. 2) go to 20
c
c *** introduce charge and strangeness exchange reactions ***
      iplab=ifix(p*2.5)+1
      if (iplab .gt. 10) iplab=10
      call grndm(rndm,1)
      if (rndm(1) .gt. (cech(iplab)/atno2**0.42)) go to 120
      call grndm(rndm,1)
      ran=rndm(1)
      irn=ifix(ran*5.)+1
      if (nfl .eq. 2) irn=5+ifix(ran*7.)+1
      if (nfl .eq. 1) irn=max(irn,5)
      if (nfl .eq. 2) irn=max(irn,12)
      ipa(1)=iipa(irn,1)
      ipa(2)=iipa(irn,2)
      ncech=1
      go to 120
c
c --- check if energetically possible to produce one extra pion ---
 20   continue
      if (eab .le. rmass(7)) go to 55
c
c --- no. of total particles vs sqrt(s)-mp-msm ---
      aleab=log(eab)
      n=3.62567+0.665843*aleab+0.336514*aleab*aleab
     * +0.117712*aleab*aleab*aleab+0.0136912*aleab*aleab*aleab*aleab
      n=n-2.
c
c --- normalization constant for  kno-distribution ---
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
 21   continue
      anpn=1./anpn
c
c --- check for target nucleon type ---
      call grndm(rndm,1)
      ran=rndm(1)
      excs=0.
      go to (30,40),nfl
c
c --- proton target ---
 30   continue
      l=0
      do 31 np1=1,20
      np=np1-1
      nmm1=np1-2
      if (nmm1 .le. 0) nmm1=1
      npp1=np1+1
      do 31 nm1=nmm1,npp1
      nm=nm1-1
      do 31 nz1=1,20
      nz=nz1-1
      l=l+1
      if (l .gt. 1200) go to 31
      nt=np+nm+nz
      if ((nt .le. 0) .or. (nt .gt. 60)) go to 31
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
      if (ran .lt. excs) go to 100
   31 continue
      go to 80
c
c --- neutron target ---
 40   continue
      l=0
      do 41 np1=1,20
      np=np1-1
      nmm1=np1-1
      if (nmm1 .le. 0) nmm1=1
      npp1=np1+2
      do 41 nm1=nmm1,npp1
      nm=nm1-1
      do 41 nz1=1,20
      nz=nz1-1
      l=l+1
      if (l .gt. 1200) go to 41
      nt=np+nm+nz
      if ((nt .le. 0) .or. (nt .gt. 60)) go to 41
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
      if (ran .lt. excs) go to 100
   41 continue
      go to 80
c
 50   continue
c      if (nprt(4)) write(newbcd,1003) eab,n,nfl,np,nm,nz
 1003 format(' *casx0* xi0 -induced cascade,',
     $ ' avail. energy',2x,f8.4,
     $ 2x,'<ntot>',2x,f8.4,2x,'from',4(2x,i3),2x,'particles')
      if (int .eq. 1) call twob(27,nfl,n)
      if (int .eq. 2) call genxpt(27,nfl,n)
      go to 9999
c
c *** energetically not possible to produce one extra pion ***
 55   continue
c      if (nprt(4)) write(newbcd,1001)
 1001 format('0*casx0* cascade energetically not possible',
     $ ' continue with quasi-elastic scattering')
      go to 53
c
c *** exclusive reaction not found ***
 80   continue
c      if (nprt(4)) write(newbcd,1004) rs,n
 1004 format(' *casx0* xi0 -induced cascade,',
     $ ' exclusive reaction not found',
     $ ' try elastic scattering  avail. energy',2x,f8.4,2x,
     $ '<ntot>',2x,f8.4)
c
 53   continue
      int=1
      np=0
      nm=0
      nz=0
      ipa(1)=26
      ipa(2)=14
      if (nfl .eq. 2) ipa(2)=16
      go to 120
c
c *** inelastic interaction has occurred ***
c *** number of secondary mesons determined by kno distribution ***
 100  continue
      do 101 i=1,60
      ipa(i)=0
 101  continue
c
      if (int .le. 0) go to 131
c
c --- take target nucleon type into account ---
      go to (102,112),nfl
c
c --- proton target ---
 102  continue
c --- check for total charge of final state mesons to determine ---
c --- the kind of baryons to be produced taking into account    ---
c --- charge and strangeness conservation                       ---
      ncht=np-nm
      if (ncht .lt. 1) go to 103
      if (ncht .eq. 1) go to 104
      if (ncht .gt. 1) go to 105
c
 103  continue
c --- xi0 p ---
      ipa(1)=26
      ipa(2)=14
      if (ncht .eq. 0) go to 120
c --- charge mismatch ==> take a s+ and correct the strangeness ---
c --- by replacing a pi- by k- ---
c --- s+ p ---
      ipa(1)=20
      ipa(2)=14
      nvefix=1
      go to 120
c
 104  continue
c --- xi0 n ---
      ipa(1)=26
      ipa(2)=16
      call grndm(rndm,1)
      if (rndm(1) .lt. 0.5) go to 120
c --- xi- p ---
      ipa(1)=27
      ipa(2)=14
      go to 120
c
 105  continue
c --- xi- n ---
      ipa(1)=27
      ipa(2)=16
      go to 120
c
c --- neutron target ---
 112  continue
c --- check for total charge of final state mesons to determine ---
c --- the kind of baryons to be produced taking into account    ---
c --- charge and strangeness conservation                       ---
      ncht=np-nm
      if (ncht .lt. 0) go to 113
      if (ncht .eq. 0) go to 114
      if (ncht .gt. 0) go to 115
c
 113  continue
c --- xi0 p ---
      ipa(1)=26
      ipa(2)=14
      if (ncht .eq. -1) go to 120
c --- charge mismatch ==> take a s+ and correct the strangeness ---
c --- by replacing a pi- by k- ---
c --- s+ p ---
      ipa(1)=20
      ipa(2)=14
      nvefix=1
      go to 120
c
 114  continue
c --- xi0 n ---
      ipa(1)=26
      ipa(2)=16
      call grndm(rndm,1)
      if (rndm(1) .lt. 0.5) go to 120
c --- xi- p ---
      ipa(1)=27
      ipa(2)=14
      go to 120
c
 115  continue
c --- xi- n ---
      ipa(1)=27
      ipa(2)=16
c
c --- take pions for all secondary mesons ---
c --- replace pi by k in case of strangeness to be fixed ---
 120  continue
      nt=2
c
      if (np .eq. 0) go to 122
c
c --- pi+ ---
      do 121 i=1,np
      nt=nt+1
      ipa(nt)=7
 121  continue
c
 122  continue
      if (nm .eq. 0) go to 124
c
c --- pi- ---
      do 123 i=1,nm
      nt=nt+1
      ipa(nt)=9
      if (nvefix .ge. 1) ipa(nt)=13
c      if (nprt(4) .and. (nvefix .ge. 1)) print 3000
 3000 format(' *casx0* k- introduced')
      nvefix=nvefix-1
 123  continue
c
 124  continue
      if (nz .eq. 0) go to 130
c
c --- pi0 ---
      do 125 i=1,nz
      nt=nt+1
      ipa(nt)=8
 125  continue
c
c --- all secondary particles have been defined ---
c --- now go for momenta and x values ---
 130  continue
c      if (nprt(4)) write(newbcd,2004) nt,(ipa(i),i=1,60)
 2004 format(' *casx0* ',i3,' particles produced. mass index array : '/
     $ 3(1h ,20(i3,1x)/))
      go to 50
c
 131  continue
c      if (nprt(4)) write(newbcd,2005)
 2005 format(' *casx0* no particles produced')
c
 9999 continue
      return
      end
