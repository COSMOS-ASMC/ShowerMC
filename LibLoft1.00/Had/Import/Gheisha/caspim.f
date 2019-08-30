*cmz :  3.14/16 13/03/89  14.48.41  by  nick van eijndhoven (cern)
*-- author :
      subroutine caspim(k,int,nfl)
c
c *** cascade of pi- ***
c *** nve 04-may-1988 cern geneva ***
c
c origin : h.fesefeldt 13-sep-1987
c
c pi-  undergoes interaction with nucleon within nucleus.
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
      data cech/1.,0.95,0.79,0.32,0.19,0.16,0.14,0.12,0.10,0.08/
      data b/0.7,0.7/,c/1.25/
c
c --- initialization indicated by kginit(16) ---
      if (kginit(16) .ne. 0) go to 10
      kginit(16)=1
c
c --- initialize pmul and anorm arrays ---
      do 9000 j=1,1200
      do 9001 i=1,2
      pmul(i,j)=0.0
      if (j .le. 60) anorm(i,j)=0.0
 9001 continue
 9000 continue
c
c *** computation of normalization constants ***
c
c --- p target ---
      l=0
      do 1100 np1=1,20
      np=np1-1
      nmm1=np1-1
      if (nmm1 .le. 1) nmm1=1
      npp1=np1+1
c
      do 1101 nm1=nmm1,npp1
      nm=nm1-1
c
      do 1102 nz1=1,20
      nz=nz1-1
      l=l+1
      if (l .gt. 1200) goto 1199
      nt=np+nm+nz
      if (nt .le. 0) go to 1102
      if (nt .gt. 60) go to 1102
      pmul(1,l)=pmltpc(np,nm,nz,nt,b(1),c)
      anorm(1,nt)=anorm(1,nt)+pmul(1,l)
 1102 continue
c
 1101 continue
c
 1100 continue
c
 1199 continue
c
c --- n target ---
      l=0
      do 1200 np1=1,20
      np=np1-1
      npp1=np1+2
c
      do 1201 nm1=np1,npp1
      nm=nm1-1
c
      do 1202 nz1=1,20
      nz=nz1-1
      l=l+1
      if (l .gt. 1200) go to 1299
      nt=np+nm+nz
      if (nt .le. 0) go to 1202
      if (nt .gt. 60) go to 1202
      pmul(2,l)=pmltpc(np,nm,nz,nt,b(2),c)
      anorm(2,nt)=anorm(2,nt)+pmul(2,l)
 1202 continue
c
 1201 continue
c
 1200 continue
c
 1299 continue
c
      do 3 i=1,60
      if (anorm(1,i) .gt. 0.0) anorm(1,i)=1.0/anorm(1,i)
      if (anorm(2,i) .gt. 0.0) anorm(2,i)=1.0/anorm(2,i)
    3 continue
c
      if (.not. nprt(10)) go to 10
      write(newbcd,2001)
      do 4 nfl=1,2
      write(newbcd,2002) nfl
      write(newbcd,2003) (anorm(nfl,i),i=1,60)
      write(newbcd,2003) (pmul(nfl,i),i=1,1200)
    4 continue
c
c --- choose proton or neutron as target ---
 10   continue
      nfl=2
      call grndm(rndm,1)
      if (rndm(1) .lt. zno2/atno2) nfl=1
      tarmas=rmass(14)
      if (nfl .eq. 2) tarmas=rmass(16)
      s=amasq+tarmas**2+2.0*tarmas*en
      rs=sqrt(s)
      enp(8)=amasq+tarmas**2+2.0*tarmas*enp(6)
      enp(9)=sqrt(enp(8))
      eab=rs-tarmas-rmass(9)
c
c --- elastic scattering ---
      np=0
      nm=0
      nz=0
      n=0.0
      ncech=0
      ipa(1)=9
      ipa(2)=14
      if (nfl .eq. 2) ipa(2)=16
      if (int .eq. 2) goto 20
      goto 100
c
c --- check if energetically possible to produce one extra pion in react.
 20   continue
      if (eab .le. rmass(9)) go to 55
c
c --- suppression of high multiplicity events at low momentum ---
      ieab=ifix(eab*5.0)+1
      if (ieab .gt. 10) go to 22
      call grndm(rndm,1)
      if (rndm(1) .lt. supp(ieab)) go to 22
c
c --- charge exchange reaction (is included in inelastic cross section)
      iplab=ifix(p*5.0)+1
      if (iplab .gt. 10) iplab=10
      call grndm(rndm,1)
      if (rndm(1) .gt. cech(iplab)) go to 23
c
      if (nfl .eq. 1) goto 24
c
c --- n target ---
      int=1
      ipa(1)=9
      ipa(2)=16
      go to 100
c
c --- p target ---
 24   continue
      ipa(1)=8
      ipa(2)=16
      go to 100
c
 23   continue
      n=1.0
c
      if (nfl .eq. 1) go to 26
c
c --- n target ---
      dum=-(1+b(2))**2/(2.0*c**2)
      if (dum .lt. expxl) dum=expxl
      if (dum .gt. expxu) dum=expxu
      w0=exp(dum)
      dum=-(-1+b(2))**2/(2.0*c**2)
      if (dum .lt. expxl) dum=expxl
      if (dum .gt. expxu) dum=expxu
      wm=exp(dum)
      call grndm(rndm,1)
      ran=rndm(1)
      np=0
      nm=0
      nz=1
      if (ran .lt. w0/(w0+wm)) go to 50
      np=0
      nm=1
      nz=0
      go to 50
c
c --- p target ---
 26   continue
      dum=-(1+b(1))**2/(2.0*c**2)
      if (dum .lt. expxl) dum=expxl
      if (dum .gt. expxu) dum=expxu
      w0=exp(dum)
      wp=exp(dum)
      dum=-(-1+b(1))**2/(2.0*c**2)
      if (dum .lt. expxl) dum=expxl
      if (dum .gt. expxu) dum=expxu
      wm=exp(dum)
      wp=wp*10.
      wt=w0+wp+wm
      wp=w0+wp
      call grndm(rndm,1)
      ran=rndm(1)
      np=0
      nm=0
      nz=1
      if (ran .lt. w0/wt) go to 50
      np=1
      nm=0
      nz=0
      if (ran .lt. wp/wt) go to 50
      np=0
      nm=1
      nz=0
      goto 50
c
 22   continue
      aleab=log(eab)
c
c --- no. of total particles vs sqrt(s)-2*mp ---
      n=3.62567+0.665843*aleab+0.336514*aleab*aleab
     $ +0.117712*aleab*aleab*aleab+0.0136912*aleab*aleab*aleab*aleab
      n=n-2.0
c
c --- normalization constant for kno-distribution ---
      anpn=0.0
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
      anpn=1.0/anpn
c
      call grndm(rndm,1)
      ran=rndm(1)
      excs=0.0
      if (nfl .eq. 2) go to 40
c
c --- p target ---
      l=0
      do 310 np1=1,20
      np=np1-1
      nmm1=np1-1
      if (nmm1 .le. 1) nmm1=1
      npp1=np1+1
c
      do 311 nm1=nmm1,npp1
      nm=nm1-1
c
      do 312 nz1=1,20
      nz=nz1-1
      l=l+1
      if (l .gt. 1200) go to 80
      nt=np+nm+nz
      if (nt .le. 0) go to 312
      if (nt .gt. 60) go to 312
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
      if (ran .lt. excs) goto 50
 312  continue
c
 311  continue
c
 310  continue
      goto 80
c
c --- n target ---
 40   continue
      l=0
      do 410 np1=1,20
      np=np1-1
      npp1=np1+2
c
      do 411 nm1=np1,npp1
      nm=nm1-1
c
      do 412 nz1=1,20
      nz=nz1-1
      l=l+1
      if (l .gt. 1200) go to 80
      nt=np+nm+nz
      if (nt .le. 0) go to 412
      if (nt .gt. 60) go to 412
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
      if (ran .lt. excs) goto 50
 412  continue
c
 411  continue
c
 410  continue
      go to 80
c
 50   continue
      if (nfl .eq. 2) go to 65
c
c --- p target ---
      if (np .eq. nm) go to 61
      if (np .eq. 1+nm) go to 63
      ipa(1)=8
      ipa(2)=14
      go to 100
c
 61   continue
      call grndm(rndm,1)
      if (rndm(1) .lt. 0.75) go to 62
      ipa(1)=8
      ipa(2)=16
      go to 100
c
 62   continue
      ipa(1)=9
      ipa(2)=14
      go to 100
c
 63   continue
      ipa(1)=9
      ipa(2)=16
      go to 100
c
c --- n target ---
 65   continue
      if (np .eq. -1+nm) go to 66
      if (np .eq. nm) go to 68
      ipa(1)=8
      ipa(2)=16
      go to 100
c
 66   continue
      call grndm(rndm,1)
      if (rndm(1) .lt. 0.50) go to 67
      ipa(1)=8
      ipa(2)=16
      go to 100
c
 67   continue
      ipa(1)=9
      ipa(2)=14
      go to 100
c
 68   continue
      ipa(1)=9
      ipa(2)=16
      go to 100
c
 70   continue
      if (nprt(4)) write(newbcd,1003) eab,n,nfl,np,nm,nz
      call stpair
      if (int .eq. 1) call twob(9,nfl,n)
      if (int .eq. 2) call genxpt(9,nfl,n)
      go to 9999
c
c --- energetically not possible to produce cascade-particles ---
c --- continue with quasi-elastic scattering ---
 55   continue
      if (nprt(4)) write(newbcd,1001)
      go to 53
c
c --- exclusive reaction not found ---
 80   continue
      if (nprt(4)) write(newbcd,1004) rs,n
c
 53   continue
      int=1
      np=0
      nm=0
      nz=0
      n=0.0
      ipa(1)=9
      ipa(2)=14
      if (nfl .eq. 2) ipa(2)=16
c
 100  continue
      do 101 i=3,60
      ipa(i)=0
 101  continue
      if (int .le. 0) go to 131
c
 120  continue
      nt=2
      if (np .eq. 0) go to 122
      do 121 i=1,np
      nt=nt+1
      ipa(nt)=7
 121  continue
c
 122  continue
      if (nm .eq. 0) go to 124
      do 123 i=1,nm
      nt=nt+1
      ipa(nt)=9
 123  continue
c
 124  continue
      if (nz .eq. 0) go to 130
      do 125 i=1,nz
      nt=nt+1
      ipa(nt)=8
 125  continue
c
 130  continue
      if (nprt(4)) write(newbcd,2004) nt,(ipa(i),i=1,20)
      if (ipa(1) .eq. 7) np=np+1
      if (ipa(1) .eq. 8) nz=nz+1
      if (ipa(1) .eq. 9) nm=nm+1
      go to 70
c
 131  continue
      if (nprt(4)) write(newbcd,2005)
c
1001  format('0*caspim* cascade energetically not possible',
     $ ' continue with quasi-elastic scattering')
1003  format(' *caspim* pion- -induced cascade, avail. energy',2x,f8.4,
     $ 2x,'<ntot>',2x,f8.4,2x,'from',4(2x,i3),2x,'particles')
1004  format(' *caspim* pion- -induced cascade, exclusive reaction',
     $ ' not found try elastic scattering  avail. energy',2x,f8.4,2x,
     * '<ntot>',2x,f8.4)
2001  format('0*caspim* tables for multiplicity data pion- induced ',
     $ 'reaction for definition of numbers see fortran coding')
2002  format(' *caspim* target particle flag',2x,i5)
2003  format(1h ,10e12.4)
2004  format(' *caspim* ',i3,2x,'particles, mass index array',2x,20i4)
2005  format(' *caspim* no particles produced')
c
 9999 continue
      return
      end
