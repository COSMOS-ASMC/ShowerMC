c******************************************************************
c...A main to test some subroutines.
c******************************************************************
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      common/pjjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/pjdatr/mrpy(6),rrpy(100) 
      save /pjdatr/

c....Setup particle
      call jamsetpa
c...Random seed.
      mrpy(1)=190503

      itest=15
c     call xxx

c...Particle data.
c...(1) Output particle data defined by ludat
c...(2) Updat particle table.
      if(itest.eq.1)  call plist
      if(itest.eq.2)  call updat(2)
      if(itest.eq.3)  call code_tex

c...(11) determine resonance ID corresponding to mass. 
c...(12) Check ID of N* D*
c...(13) test particle charge.
c...(14) Test partial decay width.
c...(15) Test resonance decay.
      if(itest.eq.11) call idres_t
      if(itest.eq.12) call idnsds
      if(itest.eq.13) call kfch_t
      if(itest.eq.14) call test_pd
      if(itest.eq.15) call test_d

c...Cross sections.
c...(21) Test nn cross sections.
c...(22) Test cross sections.
c...(23) Test k-p cross sections.
c...(24) Test pi-N cross sections.
      if(itest.eq.21)  call cross_nn
      if(itest.eq.22)  call cross_t
      if(itest.eq.23)  call cross_kp
      if(itest.eq.24)  call cross_pin

c...PYHIA/JETSET
c...(31) Test fragmentation of jet system.
c...(32) Test hard interaction i.e. pythia routine.
      if(itest.eq.31)  call jetset_t
      if(itest.eq.32)  call jetset2_t
      if(itest.eq.33)  call test_py

c...Hard
c...(41) Structure function
c...(42) Calculate jet cross section
      if(itest.eq.41) call stfun
      if(itest.eq.42) call jet_t

c....Detbal. and integrals.
c...(51) absorption cross section by detailed balance; NR
c...(52) absorption cross section by detailed balance; RR
c...(53) absorption cross section by detailed balance; all
c...(54) Make data table for detailed balance.
      if(itest.eq.51) call detbalt1  ! NR
      if(itest.eq.52) call detbalt2  ! RR
      if(itest.eq.53) call detbalt3
      if(itest.eq.54) call integpar
      if(itest.eq.55) call bwres_t

c...Test resonance productions.
c...(61) t-channel resonance
      if(itest.eq.61) call test_res

      end

c**********************************************
      subroutine xxx

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      character cname*16

C...Define particle type codes.
        do kc=1,500
        kf=kchg(kc,4) 
        kfa=iabs(kf)
        kfl1=mod(kfa,10)
        kfl2=mod(kfa/10,10)
        kfl3=mod(kfa/100,10)
        kfl4=mod(kfa/1000,10)
        if(kfa.lt.10) then
          kfa=1
        else if(kfa.eq.21) then
          kfa=2
        else if(kfl2.eq.0) then  ! di-quark.
          kfa=100
        elseif(kfa.ge.91.and.kfa.le.93) then  ! clusters.
          kfa=3
        else if(kfa.le.100) then
          kfa=200
        else if(kfl4.eq.0) then   ! mesons
          kfa=4
        else  if(kfl4.ne.0) then  ! baryons
          kfa=5
        else
          kfa=999
        endif

c        kc1=0
c        kf1=0
c        if(kf.ne.0) then
c        kc1=jamcomp(kf)
c        kf1=kchg(kc1,4)
c        endif
c        print *,kc,kf,' ',cname,' ',kfa,'kc1 kf1',kc1,kf1
c        if(kc.ne.kc1) print *,'kc error',kc,kc1
c        if(kf.ne.kf1) print *,'kf error',kf,kf1

      end do

      end

c**********************************************

      subroutine plist

      implicit double precision(a-h, o-z)
      include 'jam2.inc'

c...(1) Print particle list.
      mstu(11)=11
      call pjlist(11)
      mstu(11)=12
      call pjlist(12)
      mstu(11)=13
      call pjlist(13)
      mstu(11)=16
      call pjlist(16)

      call jamupdat(1,30)

      do kc=1,500
       write(40,800)kc,kchg(kc,4),pmas(kc,1),pmas(kc,2),
     $    chaf(kc,(3-isign(1,kchg(kc,4)))/2)
      end do
 800  format(i3,1x,i9,1x,2(f11.3),1x,a16)
      end

c************************************************

      subroutine updat(mup)


c...Purpose: to update particle data table.
      implicit double precision(a-h, o-z)

      if(mup.eq.1) then
        open(30,file='part_new.dat',status='new')
        call jamupdat(1,30)
      else

c....Read new data
      open(40,file='part_new.dat',status='old')
      call jamupdat(2,40)

c....Write new data
      open(50,file='part_new.dat1',status='new')
      call jamupdat(1,50)

c...Make new block data.
      open(60,file='pydat.f',status='unknown')
      call jamupdat(4,60)
      endif

      end

c************************************************

      subroutine jet_t

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      common/hijdat/hidat0(10,10),hidat(10)
      common/hipyint/mint4,mint5,atco(200,20),atxs(0:200)
      save  /hiparnt/,/hijdat/,/hipyint/
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      save /pjpars/,/pjsubs/


c     isoft=0  ! No Soft part
      isoft=1  ! Soft cross section from Hijing

      ihpr2(10)=1  ! print warning messages
      ihnt2(1)=1   ! proj. mass
      ihnt2(2)=1   ! proj. z
      ihnt2(3)=1   ! targ. mass
      ihnt2(4)=1   ! targ. z
      ihnt2(5)=2212
      ihnt2(6)=2212
      hint1(8)=parc(28) ! nucl. mass
      hint1(9)=parc(28)
      ihpr2(6)=0  ! nuclear effect on the parton dist.

c...Second order running alpha_strong
      mstp(2)=2
c...Inclusion of K factor(d=0)
      mstp(33)=1
      parp(31)=hipr1(17)

c...Initial state radiation
      mstp(61)=1

c...Final state radiation
      mstp(71)=1


c...Effective minimum transverse momentum p_Tmin for multiple
c...interactions with mstp(82)=1  (D=1.40GeV/c)
      hipr1(8)=2.0d0
      parp(81)=hipr1(8)

        ckin(5)=hipr1(8)
        ckin(3)=hipr1(8)  ! minimum p_t transfer in (semi)hard scatt.
        ckin(4)=hipr1(9)
        if(hipr1(9).le.hipr1(8)) ckin(4)=-1.0d0
        ckin(9)=-10.0d0
        ckin(10)=10.0d0

      ih1=10
      ih2=11
      open(ih1,file='jet1.dat',status='unknown')
      open(ih2,file='jet2.dat',status='unknown')

c...Parton structure functions (d=9)
      mstc(82)=15  ! GRV 94L
c     mstc(82)=9   ! CTEQ2L
c     mstc(82)=3  ! Duke-Owens set 1 structure functions

      mstp(51)=mstc(82)
      ihpr2(7)=mstc(82)

      nmax=50
      smin=9.0d0
      parp(2)=smin
      parc(71)=smin
      smax=100000.d0
      b=1.d0/nmax*log10(smax/smin)
      do is=0,nmax

      hint1(1)=smin*10**(is*b)

c...energy dependent Pt cut off
c     hipr1(8)=2.5+0.12*log10(hint1(1)/50)**3
c     parp(81)=hipr1(8)
c...VNI
c     hipr1(8)=0.35d0*(hint1(1)**2)**0.14d0

      ptcut=hipr1(8)

      if(isoft.eq.1) then
c...Soft cross section from Hijing table fit
        i=0
20      i=i+1
        if(i.eq.10) go to 30
        if(hidat0(10,i).le.hint1(1)) go to 20
30      if(i.eq.1) i=2
        do 40 j=1,9
           hidat(j)=hidat0(j,i-1)+(hidat0(j,i)-hidat0(j,i-1))
     &     *(hint1(1)-hidat0(10,i-1))/(hidat0(10,i)-hidat0(10,i-1))
40      continue

c....Soft cross section
        hipr1(31)=hidat(5)
      else
        hipr1(31)=10.0d0
      endif

c     hipr1(30)=2.0*hidat(5)
      hipr1(30)=2.0d0*hipr1(31)

      call hijcrs

      sigjeta=hint1(10)
      sigin=hint1(12)
      sigtot=hint1(13)
      sigel=hint1(13)-hint1(12)
      sigelr=sigel/sigtot
      sig0=hipr1(31)
c     sigtrg=hint1(59)
c...diffractive to inel ratio   hipr1(33)
      sigjete=hint1(11)   ! average jet cross section with shadowing effect
      sigjet0=hint1(14)   ! jet cross section 
      sigjetp=hint1(15)   ! proj
      sigjett=hint1(16)   ! targ
      sigjetc=hint1(17)   ! cross term
c     sigjet =hint1(18)   ! effective jet x section

      write(ih1,800)hint1(1),ptcut,
     $   sigjete,sigjet0,sigjetp,sigjett,sigjetc

      write(ih2,800)hint1(1),sigtot,sigel,sigin,sigjeta,sig0,sigelr

      write(6,800)hint1(1),sigtot,sigel,sigin,sigjeta,sigelr

800   format(f12.3,8(1x,f9.3))

      end do

      end

c***********************************************************************

      subroutine stfun

c...Test structure functions.
      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      common/njet/n,ip_crs
      save  /hiparnt/,/njet/
      double precision f(2,7) ,qq,dx1,dx2

      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      save /pjpars/,/pjint1/
      dimension xpq(-25:25)

      q2=20.0d0
      qq=dble(q2)
      x=0.5d0
      dx1=dble(x)
      dx2=dble(x)
 
      kf=2212

c...Hijing
      ip_crs=0
      ihpr2(6)=0 ! nuclear shadowing effect off
      ihpr2(7)=1 ! set 1
      hipr1(15)=0.2d0
      hipr1(16)=2.0d0

c...Pythia
      mstp(51)=3

      dx=1.d0/100
      do i=1,100
        x=i*dx
        dx1=dble(x)
        dx2=dble(x)
        call parton(f,dx1,dx2,qq)
        call pjpdfu(kf,x,q2,xpq)
        write(6,800)
     $   x,f(1,7),xpq(0),f(1,3),xpq(1)-xpq(-1),f(1,1),xpq(2)-xpq(-2)
     $      ,f(1,2),xpq(-1)
      end do

800   format(9(f7.4,1x))
c     write(6,*)'x q2',x,q2
c     write(6,*)'f1= ',(f(1,i),i=1,7)
c     write(6,*)'f2= ',(f(2,i),i=1,7)
c     write(6,*)'xpq=',(xpq(i),i=-25,25)
c     write(6,*)' '
c     write(6,*)'f1= ',(f(1,i),i=1,7)
c     write(6,*)'xpq 0 1-6= ',(xpq(i),i=0,6)
c     write(6,*)'xpq sea= ',(xpq(i),i=-6,-1)

      end

c***********************************************************************

      subroutine cross_nn

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      dimension sig1(10),sig0(10),sigin(30),sigin1(9),sigin2(30)
      parameter ( emnuc = 0.9383d0, empion = 0.138d0 )
      parameter(nchnl=12)
      character simfile4(nchnl)*15,cfile4(nchnl)*15,simfile5(nchnl)*15
      dimension sigt(0:6,0:6,0:5)
c...Functions: momentum in two-particle cm.
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a) 
c...Function:lab. momentum.
      plabsr(a,b,c)=sqrt((a**2-b**2-c**2)**2/(4.d0*c**2)-b**2)

      data simfile4/'pp-ND.fit','pp-NNs.fit','pp-DD.fit','pp-NDs.fit',
     $              'pp-NsD.fit','pp-DDs.fit','pp-NsNs.fit',
     $              'pp-NsDs.fit','pp-DsDs.fit',
     $              'ppinel.fit','pp-NR.fit','pp-RR.fit'/

      data simfile5/'pn-ND.fit',  'pn-NNs.fit','pn-DD.fit','pn-NDs.fit',
     $              'pn-NsD.fit', 'pn-DDs.fit','pn-NsNs.fit',
     $              'pn-NsDs.fit','pn-DsDs.fit',
     $              'pninel.fit', 'pn-NR.fit', 'pn-RR.fit'/

      data cfile4/'N+D(1232)','NN*','DD','ND*',
     $              'N*D','DD*','N*N*',
     $              'N*D*','D*D*','Inel','NR','RR'/

      ipn=1 ! pp
c     ipn=2 ! pn

      ifile=10
      do i=1,nchnl
        ifile=ifile+1
        if(ipn.eq.1) open(ifile,file=simfile4(i),status='unknown')
        if(ipn.eq.2) open(ifile,file=simfile5(i),status='unknown')
        write(ifile,'(a)')'# '//cfile4(i)
        write(ifile,'(''# plab(gev/c) srt(GeV)  sigma(mb)  mult.'')')
      enddo

      kf1=2212
      kf2=2212
      if(ipn.eq.2) kf2=2112
      em1=emnuc
      em2=emnuc
      smin=em1+em2+0.0001d0
c     ds=0.1d0
      ds=0.05d0
      smax=20.0d0
      pare(3)=0.0d0
      icltyp=1
        nnn=nint((smax-smin)/ds)

      do i=0,nnn
        srt=smin+i*ds

        call jamxnnin(srt,sig1,ipn)
        if(ipn.eq.2) call jamxnnin(srt,sig0,0)

        plab=plabsr(srt,em1,em2)
        pr=pawt(srt,em1,em2)
        call jamcross(1,icltyp,srt,pr,kf1,kf2,em1,em2,
     $                 sig,sigel,sigin,mchanel,mabsrb,ijet,icon)

        sigint2=sig-sigel
        ppinel=30.9d0-28.9d0*plab**(-2.46d0)+0.192d0*log(plab)**2
     $    -0.835d0*log(plab)

        call jamcross(2,icltyp,srt,pr,kf1,kf2,em1,em2,
     $                 sig,sigel,sigin2,mchanel,mabsrb,ijet,icon)
         sigint3=0.0d0
         do ic=1,mchanel
           sigint3=sigint3+sigin2(ic)
         end do
         sigin1(1)=sigin2(1)+sigin2(2)
         sigin1(2)=sigin2(3)+sigin2(4)+sigin2(18)+sigin2(19)
         sigin1(3)=sigin2(5)+sigin2(6)
         sigin1(4)=sigin2(7)+sigin2(8)
         sigin1(5)=sigin2(9)+sigin2(10)
         sigin1(6)=sigin2(11)+sigin2(12)
         sigin1(7)=sigin2(13)
         sigin1(8)=sigin2(14)+sigin2(15)
         sigin1(9)=sigin2(16)+sigin2(17)

              sigdd1=sigin2(5)   ! d0d+
              sigdd2=sigin2(6)   ! d-d++

              write(60,811)srt,sigdd1,sigdd2

c...Diffractive
c       srt0=max(10,srt)
        call jamxtot(kf1,kf2,srt,pr,sigt)
        sigd=sigt(0,0,2)+sigt(0,0,3)

        if(ipn.eq.1) then
          sigin(1)=sig1(1)            ! ND
          sigin(2)=sig1(2)+sig1(10)  ! NN*+s-wave
          sigin(3)=sig1(3)            ! DD
          sigin(4)=sig1(4)            ! ND*
          sigin(5)=sig1(5)            ! N*D
          sigin(6)=sig1(6)            ! DD*
          sigin(7)=sig1(7)            ! N*N*
          sigin(8)=sig1(8)            ! N*D*
          sigin(9)=sig1(9)            ! D*D*
        else
          sigin(1)=0.5d0*sig1(1)                    ! ND
          sigin(2)=0.5d0*(sig1(2)+sig0(2))+sig0(10) ! NN*
          sigin(3)=0.5d0*(sig1(3)+sig0(3))         ! DD
          sigin(4)=0.5d0*sig1(4)                    ! ND*
          sigin(5)=0.5d0*sig1(5)                    ! N*D
          sigin(6)=0.5d0*(sig1(6)+sig0(6))         ! DD*
          sigin(7)=0.5d0*(sig1(7)+sig0(7))         ! n*p*
          sigin(8)=0.5d0*sig1(8)                    ! N*D*
          sigin(9)=0.5d0*(sig1(9)+sig0(9))         ! D*D*
        endif

        sigint=0.0d0
        do ic=1,9
          sigint=sigint+sigin(ic)
        end do

c        if(srt.le.3.8d0) then
c          sss=sigint
c          write(50,811)srt,sss,0.1
c        else if(plab.ge.11.0.and.plab.le.12.3) then
c          sss=0.8*sigint
c          write(50,811)srt,sss,0.01
c        else if(plab.ge.23.0.and.plab.le.25.0) then
c          sss=0.85*sigint
c          write(50,811)srt,sss,0.01
c        endif

c        write(50,*)srt,sigint3,sigint
c        write(50,*)(sigin1(ic),ic=1,9)
c        write(50,*)(sigin(ic),ic=1,9)
c        write(50,*)' '


        sigstr=sigint2-sigint
        write(6,811)srt,plab,sigint,sigint3,sigint2,sigstr,sigd
        sigstr=max(0.0d0,sigstr)
        if(srt.le.3) sigstr=0.0d0

        signr=sigin(1)+sigin(2)+sigin(4)
        sigrr=sigin(3)+sigin(5)+sigin(6)+sigin(7)+sigin(8)
     $    +sigin(9)

        if(srt.le.3.6d0) then
          write(33,811)plab,srt,max(0.0d0,ppinel-sigint)
        else if(srt.ge.5.5d0) then
          write(33,811)plab,srt,sigd
        endif

c       if(srt.le.4.0d0) then
c        write(34,811)srt,sigint2
c       else
c        write(34,811)srt,sigint
c       endif

        write(11,811)plab,srt,sigin(1)          ! ND
        write(12,811)plab,srt,sigin(2)          ! NN*+s-wave
        write(13,811)plab,srt,sigin(3)          ! DD
        write(14,811)plab,srt,sigin(4)          ! ND*
        write(15,811)plab,srt,sigin(5)          ! N*D
        write(16,811)plab,srt,sigin(6)          ! DD*
        write(17,811)plab,srt,sigin(7)          ! N*N*
        write(18,811)plab,srt,sigin(8)          ! N*D*
        write(19,811)plab,srt,sigin(9)          ! D*D*
        write(20,811)plab,srt,sigint,sigint2,sigstr ! inel
        write(21,811)plab,srt,signr             ! NR
        write(22,811)plab,srt,sigrr             ! RR
      end do

811   format(f10.3,1x,f10.3,1x,30(f10.4,1x))

      end

c***********************************************************************

      subroutine cross_t

      include 'jam1.inc'
      include 'jam2.inc'
      dimension sigin(50)
      character chaf1*16,chaf2*16
c     parameter(smax=10.,ds=0.05 )
      parameter(smax=30.d0,ds=0.02d0 )
      dimension sigt(0:6,0:6,0:5)
      logical liner
c...Functions: momentum in two-particle cm.
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a) 
c...Function:lab. momentum.
      plabsr(a,b,c)=sqrt((a**2-b**2-c**2)**2/(4.d0*c**2)-b**2)

      mstc(8)=4
c     mstc(65)=0 ! const. width
      msel=1
c     msel=2
      liner=.true.
      if(msel.eq.1) open(60,file='total.dat',status='unknown')
      if(msel.eq.2) open(60,file='inel.dat',status='unknown')

c...MB
      icltyp=2
c     kf1=211    !pi+
c     kf1=-211   !pi-
c     kf1=-321   !k-
c     kf2=2212   !p
c     kf2=2112   !n
      kf1=321    !K+
      kf2=32214  ! D(1600)+
c
c...MM
c     icltyp=3
c     kf1=211
c     kf2=-211

c...BB
      icltyp=1
c     kf1=2212
c     kf2=2112
      kf1=2112
      kf2=2114
c...AntiB-B
c     icltyp=4
c     kf1=2212
c     kf2=-2212
      
c...Compressed code.
      kc1=jamcomp(kf1)
      kc2=jamcomp(kf2)
c...Strangeness.
      istr1=kchg(kc1,7)*isign(1,kf1)
      istr2=kchg(kc2,7)*isign(1,kf1)
c...Charge.
      iz1=jamchge(kf1)
      iz2=jamchge(kf2)
      write(60,'(''# istr iz='',2i3,1x,2i3)')istr1,istr2,iz1,iz2
c...Mass.
      em1=pmas(kc1,1)
      em2=pmas(kc2,1)
c...Particle name.
      call pjname(kf1,chaf1)
      call pjname(kf2,chaf2)
      smin=em1+em2+0.0001d0
      write(60,800)kf1,chaf1,em1,kf2,chaf2,em2

c...Liner
      if(liner) then
        nnn=nint((smax-smin)/ds)
c...log
      else
        nnn=100
        b=1.d0/nnn*log10(smax/smin)
      endif

      pare(3)=0.0d0
      do i=0,nnn

        if(liner) then
          srt=smin+ds*i
        else
          srt=smin*10**(i*b)
        endif

        pr=pawt(srt,em1,em2)
        plab=plabsr(srt,em1,em2)
        call jamcross(1,icltyp,srt,pr,kf1,kf2,em1,em2,
     $                 sig,sigel,sigin,mchanel,mabsrb,ijet,icon)
        if(msel.eq.2)
     $  call jamcross(2,icltyp,srt,pr,kf1,kf2,em1,em2,
     $                 sig,sigel,sigin,mchanel,mabsrb,ijet,icon)

         call jamxtot(kf1,kf2,srt,pr,sigt)
c...sigt(0,0,0) total
c...sigt(0,0,1) elastic
c...sigt(0,0,2) diffractive a+b->x+b
c...sigt(0,0,3) diffractive a+b->a+x
c...sigt(0,0,4) double diffractive a+b->x+x
c...sigt(0,0,5) non diffractive
         sigd=sigt(0,0,2)+sigt(0,0,3)
         sigdd=sigt(0,0,4)

        if(i.eq.0) write(6,'(''# mc ma ijet='',3i3)')mchanel,mabsrb,ijet

        if(msel.eq.1) then
          sigint=sig-sigel
          write(60,810)srt,plab,sig,sigel,sigint,sigd,sigdd
        else
          sigtot=0.0d0
          do j=1,mchanel+mabsrb
           sigtot=sigtot+sigin(j)
          enddo

          write(60,810)srt,plab,sigtot,sigtot+sigel
     $           ,(sigin(j),j=1,mchanel+mabsrb)

c         call xsbw1(0,srt,pr,kf1,kf2,iz1,iz2,sigr,2)
c         write(60,815)mchanel,mabsrb,
c    $          srt,pr,sigtot,(sigin(j),j=1,2),sigr
c    $          srt,pr,sigtot,(sigin(j),j=1,mchanel+mabsrb)


        endif
      end do

800   format('# ',i5,1x,a8,1x,f8.5,' + ',i5,1x,a8,1x,f8.5)
810   format(g11.5,1x,g10.3,1x,g10.4,1x,20(g10.4,1x))
815   format(2i3,f7.3,f8.4,1x,f10.4,1x,20(f10.4,1x))

      end

c************************************************

      subroutine cross_kp

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      character chaf1*16,chaf2*16
c     parameter(smax=10.,ds=0.05 )
      parameter(smax=5.d0,ds=0.01d0 )
      common/bwdat1/sigbwp(30)
      parameter(mxchan=30)
      dimension sigin(mxchan),sigy(4)
c...Functions: momentum in two-particle cm.
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a) 
c...Function:lab. momentum.
      plabsr(a,b,c)=sqrt((a**2-b**2-c**2)**2/(4.d0*c**2)-b**2)

      mstc(8)=4
      open(60,file='totalkn.dat',status='unknown')

c...MB
      icltyp=2
      ikk=1     ! k-p
c     ikk=2     ! k-n
      if(ikk.eq.1) then
        kf1=-321
        kf2=2212
        open(61,file='k-p1.dat')
        open(51,file='k-p2.dat')
      else
       kf1=-321
       kf2=2112
        open(61,file='k-n1.dat')
        open(51,file='k-n2.dat')
      endif
      
c...Compressed code.
      kc1=jamcomp(kf1)
      kc2=jamcomp(kf2)
      istr1=kchg(kc1,7)*isign(1,kf1)
      istr2=kchg(kc2,7)*isign(1,kf1)
      istr=istr1+istr2
      iz1=jamchge(kf1)/3
      iz2=jamchge(kf2)/3
      write(60,'(''# istr iz='',2i3,1x,2i3)')istr1,istr2,iz1,iz2

      em1=pmas(kc1,1)
      em2=pmas(kc2,1)

      call pjname(kf1,chaf1)
      call pjname(kf2,chaf2)
      write(60,800)kf1,chaf1,em1,kf2,chaf2,em2
      write(60,801)

      smin=em1+em2+0.0001d0

      nnn=nint((smax-smin)/ds)
      do i=0,nnn

        srt=smin+ds*i
        pr=pawt(srt,em1,em2)
        plab=plabsr(srt,em1,em2)

c...Akaon-N total/background elastic.
        call jamxkp(kf1,kf2,srt,em1,em2,sig,sigelb,sigchb,sigy)

c...Breit-Wigner cross sections.
        call jamxbw1(istr,srt,pr,kf1,kf2,iz1,iz2,sigr,sige,2)

        call jamcross(1,icltyp,srt,pr,kf1,kf2,em1,em2,
     $                 siga,sigela,sigin,mchanel,mabsrb,ijet,icon)

c  2: k-p  => lambda + pi0
c  3: k- p => sigma- pi+                                       
c  4: k- p => sigma0 pi0                                        
c  5: k- p => sigma+ pi-                                        

c  2: k- n => lambda + pi-
c  3: k- n => sigma- pi0
c  4: k- n => sigma0 pi-

      sigtt=sigr+sigelb+sigchb+sigy(1)+sigy(2)+sigy(3)+sigy(4)
      sigeltt=sigelb+sige
      write(51,810)srt,plab
     $ ,sigbwp(1)+sigchb,sigchb,sigbwp(1)
     $ ,sigbwp(2)+sigy(1),sigy(1),sigbwp(2)
     $ ,sigbwp(3)+sigy(2),sigy(2),sigbwp(3)
     $ ,sigbwp(4)+sigy(3),sigy(3),sigbwp(4)
     $ ,sigbwp(5)+sigy(4),sigy(4),sigbwp(5)

      write(3,810)srt,plab,sigy(4)

      write(61,810)srt,plab,sig,sigtt,sigel,sigeltt,sigr,sige
     $     ,siga,sigela,sigelb
c       write(60,810)srt,plab
c    $   ,sigr+sigelb+sigchb+sigy
c    $   ,sige+sigelb
c    $   ,sigbwp(1)+sigchb
c    $   ,sigelb,sigchb
c    $   ,sigr,sige,sigbwp(1)
c    $   ,sigy

      end do

801   format('#  s  plab  total  elas  chex')
800   format('# ',i5,1x,a8,1x,f8.5,' + ',i5,1x,a8,1x,f8.5)
810   format(g10.3,1x,g10.3,1x,20(g10.4,1x))
815   format(2i3,f7.3,f8.4,1x,f10.4,1x,20(f10.4,1x))

      end

c************************************************

      subroutine cross_pin

c...Calculate pi-N cross sections.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      character chaf1*16,chaf2*16
c     parameter(smax=10.,ds=0.05 )
      parameter(smax=5.d0,ds=0.01d0 )
      common/bwdat1/sigbwp(30)
      parameter(mxchan=30)
      dimension sigin(mxchan)
c...Functions: momentum in two-particle cm.
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a) 
c...Function:lab. momentum.
      plabsr(a,b,c)=sqrt((a**2-b**2-c**2)**2/(4.d0*c**2)-b**2)

c...Job. mode.
      mstc(8)=4

c...Switch for Delta decay
      mstc(63)=3

c...MB
      icltyp=2
      ipi=1     ! pi+p
c     ipi=2     ! pi-p
      if(ipi.eq.1) then
        kf1=211
        kf2=2212
        open(61,file='pi+p1.dat')
        open(51,file='pi+p2.dat')
      else
       kf1=-211
       kf2=2212
        open(61,file='pi-p1.dat')
        open(51,file='pi-p2.dat')
      endif
      
c...Compressed code.
      kc1=jamcomp(kf1)
      kc2=jamcomp(kf2)
      istr1=kchg(kc1,7)*isign(1,kf1)
      istr2=kchg(kc2,7)*isign(1,kf1)
      istr=istr1+istr2
      iz1=jamchge(kf1)/3
      iz2=jamchge(kf2)/3
      write(60,'(''# istr iz='',2i3,1x,2i3)')istr1,istr2,iz1,iz2

      em1=pmas(kc1,1)
      em2=pmas(kc2,1)

      call pjname(kf1,chaf1)
      call pjname(kf2,chaf2)
      write(60,800)kf1,chaf1,em1,kf2,chaf2,em2
      write(60,801)

      smin=em1+em2+0.0001d0

      nnn=nint((smax-smin)/ds)
      do i=0,nnn

        srt=smin+ds*i
        pr=pawt(srt,em1,em2)
        plab=plabsr(srt,em1,em2)

c...Breit-Wigner cross sections.
        call jamxbw1(istr,srt,pr,kf1,kf2,iz1,iz2,sigr,sige,2)

        call jamcross(1,icltyp,srt,pr,kf1,kf2,em1,em2,
     $                 sig,sigel,sigin,mchanel,mabsrb,ijet,icon)

c...piN cross section
        call jamxpin(iz1,iz2,srt,pr,siga,sigela)

      sigelb=0.0d0
      sigchb=0.0d0
      sigtt=sigr+sigchb+sigel
      sigeltt=sige+sigel
      write(51,810)srt,plab
     $ ,sigbwp(1)+sigchb,sigchb,sigbwp(1)

      if(srt.le.1.3) then
       sigx=0.0d0
      else
       sigx=max(0.0d0,sigela-sige)
      endif
      write(3,810)srt,plab,sigx
      write(61,810)srt,plab,sig,sigtt,sigeltt,sigel,sigr,sige
     $     ,sigchb

c       write(60,810)srt,plab
c    $   ,sigr+sigelb+sigchb+sigy
c    $   ,sige+sigelb
c    $   ,sigbwp(1)+sigchb
c    $   ,sigelb,sigchb
c    $   ,sigr,sige,sigbwp(1)
c    $   ,sigy

      end do

801   format('#  s  plab  total  elas  chex')
800   format('# ',i5,1x,a8,1x,f8.5,' + ',i5,1x,a8,1x,f8.5)
810   format(g10.3,1x,g10.3,1x,20(g10.4,1x))
815   format(2i3,f7.3,f8.4,1x,f10.4,1x,20(f10.4,1x))

      end

c************************************************

      subroutine test_pd

c...Purpose: to test momentum dependent width.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

      parameter(maxbr=70)
      dimension ibranch(maxbr),pwid(maxbr)
c     parameter (emmin=1.07,emmax=2.5,dem=0.001)
      parameter (emmin=1.09d0,emmax=2.5d0,dem=0.01d0)

      nmax=nint((emmax-emmin)/dem)

c...No momentum dependent optin.
c     mstc(65)=0
c     mstc(63)=2  ! option for delta decay width formula

      de=0.01d0
      kf1=111
      kf2=111
      kfsp=1
      kfsn=0

c     kf=2226 ! D(1905)++
c     kf=2228 ! D(1950)++
c     kf=12212
c     kf=10221 ! f_0
c     kf=2224  ! delta
      kf=22112 ! N(1535)0
      kf=12128  ! N(1990)+


      kc=jamcomp(kf)
      write(6,800)chaf(kc,3-isign(1,kf)/2)
800   format('# ',a16)
      do i=1,nmax
        emcm=emmin+i*dem

c       mstc(63)=1
        call jamwidm(kc,kfsp,kfsn,kf1,kf2,emcm,ibranch,pwid,wid1,itag)

c       mstc(63)=2
c       call jamwidm(kc,kfsp,kfsn,kf1,kf2,emcm,ibranch,pwid,wid2,itag)

c       mstc(63)=3
c       call jamwidm(kc,kfsp,kfsn,kf1,kf2,emcm,ibranch,pwid,wid3,itag)

c       mstc(63)=4
c       call jamwidm(kc,kfsp,kfsn,kf1,kf2,emcm,ibranch,pwid,wid4,itag)

        print *,emcm,wid1

c       call jamwidm(kc,kfsp,kfsn,kf1,kf2,emcm,ibranch,pwid,wid,itag)
c       if(wid.lt.1e-8) then
c         tdec=1e+35
c       else
c         tdec=paru(3)/wid
c       endif
c       print *,emcm,wid,tdec

      end do
      

      end
c************************************************
      subroutine test_d

c...Purpose: to test the decay of unstable particles. 
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      parameter (emnuc=0.9383d0,empion=0.138d0)
      common/pjjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      save /pjjets/

c...Minimum kinetic energy in decays
c     parj(64)=parc(41)
c     parj(64)=0.001d0
      iang=0
      isang=0
      gpt=parc(43)

      mstj(21)=2
      do i=1,5
        k(1,i)=0
        p(1,i)=0
      end do
      ip=1
      n=1
      n=1

      kf=32114  ! D(1600)+
c     kf=11218  ! N(1990)0
      kf=12128  ! N(1990)+
c     kf=12212  ! N(1440)+
      kf=10221  ! f0
c     kf=-10211  ! a_0-

      k(1,1)=1
      k(1,2)=kf
      k(1,3)=0
      p(1,1)= 1.019d0
      p(1,2)= -2.600d0
      p(1,3)=-4.275d0

c     p(1,5)=1.0779
c     p(1,5)=emnuc+2*empion+parj(64)
c     p(1,5)=1.2083
c     p(1,5)=emnuc+3*empion+parj(64)

      p(1,5)=1.2d0
      p(1,4)= sqrt(p(1,5)**2+p(1,1)**2+p(1,2)**2+p(1,3)**2)

      do iev=1,100
      n=1
      k(1,1)=1
      k(1,2)=kf
      p(1,5)=pjmass(kf)
      call jamrdec(icon,iang,isang,gpt) 
      call pjlist(1)
      print *,'icon=',icon 
      end do

      end

c************************************************
      subroutine jetset_t

c...Purpose: to test the fragmentation of jet sysem.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      common/pjjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
c     common/pjdat1/mstu(200),paru(200),mstj(200),parj(200) 
c     common/pjdat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
c     common/pjdat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
c     save /pjjets/,/pjdat1/,/pjdat2/,/pjdat3/ 

      common/jampos1/jqconst(2),kfcq(4),icq(4),icms
      save /jampos1/

      icms=0
c...Choice of baryon production model.
      mstj(12)=2

c...Form of particle decay.
      mstj(21)=0

c....Hadron formation point.
      mstj(10)=3

c...Choice of mother pionters.
c     mstu(16)=2

c...Meson production par.
c     parj(14)=0.6222d0
c     parj(15)=0.12d0
c     parj(16)=0.25d0
c     parj(17)=0.36d0

c...Baryon resonance production par.
c     parj(27)=0.0d0
c     parj(14)=0.0d0
c     parj(41)=1.0d0
c     parj(42)=0.7d0

       do i=1,20
       write(31,8000)i,parf(i)
       enddo
 8000  format(i4,1x,f15.6)

      nevent=100
      ymin=-2.0D0
      ymax=2.0D0
      wy=0.1D0
      nymx=(ymax-ymin)/wy
      call vbook1(11,'dN/dy proton',nymx,ymin,ymax)
      call vbook1(12,'dN/dy pi+',nymx,ymin,ymax)

      npop=0
      jqconst(1)=1
      jqconst(2)=1
      em=2.2d0

      do iev=1,nevent

c     call pj4ent(0,2101,21,21,3,15.1,x1,x2,x4,x12,x14)
c     call pj2ent(0,1,2101,em)

      iii=0
      if(iii.eq.1) theN
      n=5
      do i=1,5
       do j=1,5
       k(i,j)=0
       p(i,j)=0
       v(i,j)=0
       end do
       k(i,1)=2
      end do
      k(5,1)=1

      k(1,2)=2
      k(2,2)=21
      k(3,2)=21
      k(4,2)=21
      k(5,2)=2203

      p(1,3)=0.82255
      p(1,4)=0.82257
      p(1,5)=0.00560
    
      p(2,1)=14.94472 
      p(2,2)=0.40513
      p(2,3)=29.13731
      p(2,4)=32.74892
      p(2,5)=0.00000

      p(3,1)=2.83848
      p(3,2)=0.72596
      p(3,3)=7.01564
      p(3,4)=7.60284
      p(3,5)=0.00000

      p(4,1)=-3.55951
      p(4,2)=-1.19170
      p(4,3)=-6.04311
      p(4,4)=7.11403
      p(4,5)=0.00000

      p(5,1)= -14.22370 
      p(5,2)=0.06061
      p(5,3)=-30.93240
      p(5,4)=34.05474
      p(5,5)=0.77133

      else

      n=2
      do i=1,5
      k(1,i)=0
      k(2,i)=0
      end do
      k(1,1)=2
      k(2,1)=1
      k(1,2)= -3
      k(2,2)= 3
      p(1,1)=0d0
      p(1,2)=0.00000
      p(1,3)=0.95269
      p(1,4)=0.97325
      p(1,5)=0.19900
      p(2,1)=0.00000 
      p(2,2)=0.00000
      p(2,3)=-0.95269
      p(2,4)=0.97325
      p(2,5)=0.19900
 

      endif


      write(6,'(///,''event='',i4)')iev
      call pjexec
      call pjlist(1)

        kfl4=mod(abs(k(n,2))/1000,10)
        if(kfl4.eq.0) npop=npop+1

       write(6,*)'kfcq',(kfcq(j),j=1,3)
       write(6,*)'icq ',(icq(j),j=1,3)
 
      do 390 i=1,n

c......This is already removed particle.
        if(k(i,1).ge.10.or.k(i,1).eq.0) go to 390
        kf=k(i,2)
        rap=0.5D0*log( max(p(i,4)+p(i,3),1.D-8)/max(p(i,4)-p(i,3), 
     & 1.D-8) )
        if(k(i,4).eq.2) mstd(198)=mstd(198)+1
        if(kf.eq.2212) call vfill1(11,rap,1.d0/wy)
        if(kf.eq.221) call vfill1(12,rap,1.d0/wy)

        gamma=p(i,4)/p(i,5)
        write(6,800)i,(v(i,jj),jj=1,4),v(i,4)/gamma
     $                       ,(k(i,jj),jj=1,5)

 800   format(i3,5(1x,f8.3),1x,i3,1x,i6,i3,1x,i6,i6)
        write(11,'(i3,10(f7.3,1x))')pjk(i,16),(v(i,jj),jj=3,4),
     $          (p(i,jj)/p(i,4),jj=1,3),(p(i,jj)/p(i,4)*v(i,4),jj=1,3)

         write(40,*)i,v(i,4),v(i,3),v(i,4)/gamma

 390  continue

c....Loop over event
      end do

c...Event weight
      fac=1.d0/dble(nevent)

c...Mass distributions.
      do i=1,2
       call vscale(10+i,fac)
       call vprint(10+i,0,0)
      end do

       do i=1,2000
       write(32,9000)i,parf(i)
       enddo
 9000  format(i4,1x,f15.6)
      print *,'mstd(198)',mstd(198)
      print *,'npop=',npop

 
c     mstu(16)=1
c     mrlu(1)=239873873
c     write(6,*)'mstu16=1'
c     do iev=1,10
c     call lu2ent(0,2101,3,3.1)
c     call lu2ent(0,1,-2,1.0)
c     call luexec
c     call lulist(1)
c     write(6,810)(i,i=0,23)
c     do i=1,n
c       write(6,810)i,(klu(i,j),j=1,22)
c     end do
c     end do
c     call lulist(21)

      end

c***********************************************************************

      subroutine idres_t

c...Purpose: to determine resonance ID corresponding to mass. 
      implicit double precision(a-h, o-z)
      include 'jam2.inc'


c....Find N* channel.
        write(6,*)'N*'
        iz0=0
        iof=iz0
        kcmin=mstc(22)+iof
        kcmax=mstc(23)+iof
        istep=2
        assign 10 to label
        goto 100
  10    continue

c...Delta*
        write(6,*)'D*'
        iz0=0
        iof=iz0+1
        kcmin=mstc(24)+iof
        kcmax=mstc(25)+iof
        istep=4

        assign 20 to label
        goto 100
  20    continue


c...Lambda*
        write(6,*)'L*'
        kcmin=mstc(26)
        kcmax=mstc(27)
        istep=1

        assign 30 to label
        goto 100
  30    continue


c...Sigma*
        write(6,*)'S*'
        iof=1+iz0
        kcmin=mstc(28)+iof
        kcmax=mstc(29)+iof
        istep=3
        assign 40 to label
        goto 100
  40    continue

c...Xi*
        write(6,*)'X*'
        iof=1+iz0
        kcmin=mstc(30)+iof
        kcmax=mstc(31)+iof
        istep=2
        assign 50 to label
        goto 100
  50    continue
        return

100     continue
        do i=kcmin,kcmax,istep
          write(6,*)kchg(i,4),pmas(i,3),
     $       pmas(i,1)-pmas(i,3),pmas(i,1)+pmas(i,3)
        end do


         goto label


      end

c*******************************************************
      subroutine idnsds

c...Check ID of N* D*
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      character chaf1*16

c      kc1=lucomp(kf1)
c      kc2=lucomp(kf2)

       kcmi=mstc(22)
       kcma=mstc(25)+3
       do kc1=kcmi,kcma

        jd1=idget1(kc1)
        kf=kchg(kc1,4)
	call pjname(kf,chaf1)
	write(6,*)chaf1,' ',kc1,kf,jd1

       end do

      end

c***********************************************************************

      subroutine kfch_t

c...Purpose: to test particle charge.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      character char*5


      do kc=1,500
       kf=kchg(kc,4)
       if(kf.eq.0) goto 100

       iz1=jamchge(kf)
       iz01=kchg(kc,1)

       kf2=0
       iz2=0
       iz02=0
       if(kchg(kc,3).eq.1) then
         kf2=-kf
         iz2=jamchge(kf2)
         iz02=kchg(kc,1)*isign(1,kf2)
       endif
        char='     '
        if((iz1.ne.iz01).or.(iz2.ne.iz02)) then
            char='error'
        endif
        write(6,'(7(i7,1x),5a)')kc,kf,iz1,iz01,kf2,iz2,iz02,char
       
100   end do


      end

c*****************************************************************
      subroutine test_py

      implicit double precision(a-h, o-z)
      common/pjjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/pjdat1/mstu(200),paru(200),mstj(200),parj(200) 
      common/pjdat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/pjdat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
      save /pjjets/,/pjdat1/,/pjdat2/,/pjdat3/ 

      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      common/pjint5/ngenpd,ngen(0:500,3),xsec(0:500,3)
      common/pjuppr/nup,kup(20,7),nfup,ifup(10,2),pup(20,5),q2up(0:10)

      parameter(mxpart=300)
      common/xhit/xo(2,4),srt,ecmc
      common/partdis1/rpart(mxpart,5),ppart(mxpart,5),kpart(mxpart,5)
      common/partdis2/tcoll0,npart0,npart1,npart2,npart3,npart4

      nev=20
      srt=2000.d0
      ecmc=2000.d0
      kf1=2112
      kf2=2112
      do i=1,4
        xo(1,i)=777.0d0
        xo(2,i)=999.0d0
      end do

c-----------------------------------------------------
      mstu(11)=8
c     mstj(12)=0
c mstj(21)=0: all particle decay are inhibited
c mstj(21)=2: 
      mstj(21)=0
      mstp(111)=1 ! master switch for fragmentation
c     mstp(81)=0  ! no multiple ineraction
      mstp(125)=2 ! list complete documentation of partnic process.

c...Second order running alpha_strong
        mstp(2)=2
c...Inclusion of K factor(d=0)
        mstp(33)=1
        parp(31)=2.5d0

c...Duke-Owens set 1 structure functions (d=9)
c       mstp(51)=3

c...Initial state radiation
        mstp(61)=1

c...Final state radiation
        mstp(71)=1

c...No multiple interaction
        mstp(81)=1  ! hijing(D)

c...Structure of mutliple interaction
        mstp(82)=1

c...Fragmentation off (have to be done by local call)
        mstp(111)=0

c...Effective minimum transverse momentum p_Tmin for multiple
c...interactions with mstp(82)=1  (D=1.40GeV/c)
        parp(81)=2.0d0

        ckin(5)=2.0d0
        ckin(3)=2.0d0  ! minimum p_t transfer in (semi)hard scatt.
        ckin(4)=-1.0d0
        ckin(9)=-10.0d0
        ckin(10)=10.0d0

c**** Switch on and off scattering channels

c...QCD subprocesses
        msel=0
        do 100 isub=1,200
        msub(isub)=0
100     continue
        msub(11)=1 ! f_i f_i(~) -> f_i f_i(~)
        msub(12)=1 ! f_i f_i~ -> f_k f_k~
        msub(13)=1 ! f_i f_i~ -> g g
        msub(28)=1 ! f_i g  -> f_i g
        msub(53)=1 ! g g -> f_k f_k~
        msub(68)=1 ! g g -> g g
        msub(81)=1 ! f_i f_i~ -> Q_i Q_i~
        msub(82)=1 ! g g -> Q_i Q_i~

c...Gluon decays are not allowed
        do 102 j=1,min(8,mdcy(21,3))
102        mdme(mdcy(21,2)+j-1,1)=0

        isel=4
c...Option to switch on B-quark production.
c       if(hint1(1).ge.20.0 .and. ihpr2(18).eq.1) isel=5
        mdme(mdcy(21,2)+isel-1,1)=1


c...Direct photon production allowed
      msub(14)=1   ! q + qbar -> g     + gamma
      msub(18)=1   ! q + qbar -> gamma + gamma
      msub(29)=1   ! q + g    -> q     + gamma
c-----------------------------------------------------

      call pjinit('cms',kf1,kf2,srt,icon)
c     call pjinit('five',ihnt2(5),ihnt2(6),hint1(1),icon)
c     call pjinit('CMS','p','p',200.)
c     call systat(2)

      do iev=1,nev
10      call pjevnt

c       write(6,*)'mint31=',mint(31)
c       if(mint(31).ne.1) then
c         go to 10
c       endif

        qmax=vint(55)
        isub=mint(1)
        if(iset(isub).eq.2) qmax=sqrt(parp(71))*qmax

        write(mstu(11),*)'iev',iev
        write(mstu(11),*)'mint31 qmax=',mint(31),qmax
        itag=0
        write(mstu(11),*)' '
        write(10,*)' '
        call pjlist(1)
        mstu11=mstu(11)
        mstu(11)=10
        call pjlist(2)
        mstu(11)=mstu11
        izero=0
        do i=21,n
         if(itag.eq.0.and.k(i,1).le.10) then
           write(mstu(11),*)'============================== '
           write(10,*)'============================== '
           itag=1
         endif
         write(mstu(11),'(i3,5(1x,g12.5))')i,(v(i,j),j=1,5)
         write(10,'(i3,5(1x,g12.5))')i,(v(i,j),j=1,5)
            vv=v(i,1)**2+v(i,2)**2+v(i,3)**2
            if(itag.eq.1.and.vv.le.0.0d0.and.abs(k(i,2)).lt.100) then
            if(k(i,1).ne.13.and.kchg(jamcomp(k(i,2)),2).ne.0) izero=1
            endif
        end do
        if(izero.eq.1) then
          write(mstu(11),*)'Zero'
          write(6,*)'Zero'
          write(9,*)'Zero'
        endif

c       write(mstu(11),*)' '
c       do i=21,n
c        write(mstu(11),'(i3,5(1x,g10.3))')i,(rpart(i,j),j=1,5)
c       end do

        mstu11=mstu(11)
        mstu(11)=9
        write(mstu(11),*)'mint31=',mint(31)
        call pjlist(2)
        mstu(11)=mstu11
      end do

      end

c*****************************************************************

      subroutine detbalt1

      implicit double precision(a-h, o-z)

c...Test for RN->NN cross sections.
      include 'jam2.inc'
c     dimension sigin(10)
      character chaf1*16,chaf2*16
c...Functions: momentum in two-particle cm.
c     pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a) 

c...Switch for delta decay width
      mstc(63)=2

      msel=2
      emnuc=0.938d0
      em1=1.18d0
      em2=0.938d0
      srt=2.8d0
      s0=em1+em2+0.001d0
      sn=100.0d0
      ds=1.0d0
      nn=(sn-s0)/ds

c...nd++ -> pp
      kf1=2224
      kf2=2112

c...N*
c      idd=1
c       kcmin=mstc(22)+idd-1
c       kcmax=mstc(23)+idd-1
c...D*
c       idd=3
c       kcmin=mstc(24)+idd-3
c       kcmax=mstc(25)+idd-3

c     do kc1=kcmin,kcmax
      kc1=jamcomp(kf1)
      kc2=jamcomp(kf2)
c     fac=4./(2.*4.)*0.5
c     kf1=kchg(kc1,4)

      call pjname(kf1,chaf1)
      call pjname(kf2,chaf2)
      write(6,800)kf1,chaf1,em1,kf2,chaf2,em2
800   format(//,'# ',i5,1x,a16,1x,f8.5,' + ',i5,1x,a16,1x,f8.5)


      do i=0,nn
       srt=s0+i*ds
c       call jamxnnin(srt,sigin,1)
c       sig=0.75*sigin(1)
c       pr=pawt(srt,em1,em2)
c       pr2new=0.25*srt*srt-emnuc**2

c....Giessen
c       mstc(62)=31
c       call jamdetb1(msel,kc1,srt,em1,em2,pr,db1,nhlf)
        mstc(62)=22
        call jamdetb1(msel,kc1,srt,em1,em2,pr,db2,nhlf)
c...Danielewicz
c       mstc(62)=33
c       call jamdetb1(msel,kc1,srt,em1,em2,pr,db3,nhlf)
c
c...Naive detbal
c       sig1=fac*sig*pr2new/(pr*pr)
c       sig2=fac*sig*db1*pr2new
c       sig3=fac*sig*db2*pr2new
c       sig4=fac*sig*db3*pr2new
c       print *,srt-em1,sig1,sig2,sig3,sig4
        write(6,810)srt,db2

      end do
c     end do
810   format(f9.3,1x,g13.6)

      end

c*****************************************************************

      subroutine detbalt2

      implicit double precision(a-h, o-z)

c...Test for RR->NN cross sections.
      include 'jam2.inc'

c     dimension sigin(10)
      parameter( pi=3.14159d0, emnuc=0.939d0, empion=0.138d0)
      character chaf1*16,chaf2*16
c...Functions: momentum in two-particle cm.
c     pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a) 

c...Switch for delta decay width
      mstc(63)=2
      msel=2   ! only integral value
      em1=1.07d0
      em2=1.07d0
      srt=2.8d0
      s0=em1+em2+0.001d0
      sn=10.0d0
      ds=0.2d0
      nn=(sn-s0)/ds
      srt=2.8d0

      kf1=2214
      kf2=2214
      kc1=jamcomp(kf1)
      kc2=jamcomp(kf2)
c....d+d+->pp
      fac=2.d0*2.d0/(4.d0*4.d0)*0.5d0


c...N*
       idd=1
        kcmin=mstc(22)+idd-1
        kcmax=mstc(23)+idd-1
        isp=2
c...D*
c       idd=3
c       kcmin=mstc(24)+idd-3
c       kcmax=mstc(25)+idd-3

      do kc1=kcmin,kcmax,isp
c     kc1=lucomp(kf1)
      kc2=jamcomp(kf2)
      kf1=kchg(kc1,4)
      call pjname(kf1,chaf1)
      call pjname(kf2,chaf2)
      write(6,800)kf1,chaf1,em1,kf2,chaf2,em2
800   format(//,'# ',i5,1x,a16,1x,f8.5,' + ',i5,1x,a16,1x,f8.5)

      do i=0,nn
       srt=s0+i*ds

c       call jamxnnin(srt,sigin,1)
c       sig=0.4*sigin(3)
c       pr=pawt(srt,em1,em2) 
c       pr2new=0.25*srt*srt-emnuc**2
c       pr2new=pawt(srt,emnuc,emnuc)
c
c...Naiv detailed balance
c      sig1=fac*sig*pr2new/(pr*pr)
c...Giessen
c      mstc(62)=21
c      call jamdetb2(msel,kc1,kc2,srt,em1,em2,pr,db1,nh1)
c...My
       mstc(62)=22
       call jamdetb2(msel,kc1,kc2,srt,em1,em2,pr,db2,nh2)
c...Danielewicz
c      mstc(62)=23
c      call jamdetb2(msel,kc1,kc2,srt,em1,em2,pr,db3,nh3)
c      call jamdetb3(msel,kc1,kc2,srt,em1,em2,pr,sum,nhlf)

c       sig2=fac*sig*db1*pr2new
c       sig3=fac*sig*db2*pr2new
c       sig4=fac*sig*db3*pr2new
c       sig5=fac*sig/sum*pr2new/(pr*pr)

c      print *,srt,db1,db2,db3
c      print *,srt,sig1,sig2,sig3,sig4,sig5
        write(6,810)srt,db2

      end do
      end do
810   format(f9.3,1x,g13.6)

      end

c*****************************************************************

      subroutine detbalt3

c...Test for RR->NN coef.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

      parameter( pi=3.14159d0, emnuc=0.939d0, empion=0.138d0)
      character chaf1*16,chaf2*16

c...Switch for decay width
      mstc(63)=1   ! 1:Frankfrut  2:Giessen
      mstc(62)=22
      em1=1.07d0
      em1=1.18d0
      em2=1.07d0

      kcmin1=1
      kcmax1=1
      kcmin2=1
      kcmax2=1
      isp1=1
      isp2=1
      idb=2
      ifn=10
      isel=1
      print *,'isel?'
      read(5,*)isel

c....(1)nd
      if(isel.eq.1) then
        kf1=2214
        kf2=2212
        kcmin1=jamcomp(kf1)
        kcmin2=jamcomp(kf2)
        kcmax1=kcmin1
        kcmax2=kcmin2
        em2=0.938d0
        idb=1
c....(2)nn*
      else if(isel.eq.2) then
        idd=1
        kcmin1=mstc(22)+idd-1
        kcmax1=mstc(23)+idd-1
        isp1=2
        kf2=2212
        kcmin2=jamcomp(kf2)
        kcmax2=kcmin2
        em2=0.938d0
        idb=1
c....(3)dd
      else if(isel.eq.3) then
        kf1=2214
        kf2=2214
        kcmin1=jamcomp(kf1)
        kcmin2=jamcomp(kf2)
        kcmax1=kcmin1
        kcmax2=kcmin2
c....(4)nd*
      else if(isel.eq.4) then
        idd=3
        kcmin1=mstc(24)+idd-3
        kcmax1=mstc(25)+idd-3
        isp1=4
        kf2=2212
        kcmin2=jamcomp(kf2)
        kcmax2=kcmin2
        em2=0.938d0
        idb=1
c....(5)n*d
      else if(isel.eq.5) then
        idd=1
        kcmin1=mstc(22)+idd-1
        kcmax1=mstc(23)+idd-1
        isp1=2
        kf2=2214
        kcmin2=jamcomp(kf2)
        kcmax2=kcmin2
c....(6)dd*
      else if(isel.eq.6) then
        idd=3
        kcmin1=mstc(24)+idd-3
        kcmax1=mstc(25)+idd-3
        isp1=4
        kf2=2214
        kcmin2=jamcomp(kf2)
        kcmax2=kcmin2
c....(7)n*n*
      else if(isel.eq.7) then
        idd=1
        kcmin1=mstc(22)+idd-1
        kcmax1=mstc(23)+idd-1
        isp1=2
        kcmin2=mstc(22)+idd-1
        kcmax2=mstc(23)+idd-1
        isp2=2
c....(8)n*d*
      else if(isel.eq.8) then
        idd=1
        kcmin1=mstc(22)+idd-1
        kcmax1=mstc(23)+idd-1
        isp1=2
        idd=3
        kcmin2=mstc(24)+idd-3
        kcmax2=mstc(25)+idd-3
        isp2=4
c....(9)d*d*
      else if(isel.eq.9) then
        idd=3
        kcmin1=mstc(24)+idd-3
        kcmax1=mstc(25)+idd-3
        isp1=4
        kcmin2=mstc(24)+idd-3
        kcmax2=mstc(25)+idd-3
        isp2=4
      endif

c...Min. energy
      s0=em1+em2+0.001d0
c...Max. energy
      sn=10.0d0
c..Energy bin
      ds=0.2d0
      nn=(sn-s0)/ds

      do 1000 kc1=kcmin1,kcmax1,isp1
      do 2000 kc2=kcmin2,kcmax2,isp2

      if((isel.eq.7.or.isel.eq.9).and.(kc2.lt.kc1)) goto 2000
      kf1=kchg(kc1,4)
      kf2=kchg(kc2,4)
      call pjname(kf1,chaf1)
      call pjname(kf2,chaf2)
      write(ifn,800)kf1,chaf1,em1,kf2,chaf2,em2
800   format(//,'# ',i5,1x,a16,1x,f8.5,' + ',i5,1x,a16,1x,f8.5)

      do i=0,nn
       srt=s0+i*ds
        
c...Numerical integration.
c      if(idb.eq.1) call jamdetb1(2,kc1,srt,em1,em2,1.0,db2,nhlf)
c      if(idb.eq.2) call jamdetb2(2,kc1,kc2,srt,em1,em2,1.0,db2,nhlf)

c...Fit by Chiba
       dbfit=jambwtbl(srt,kc1,kc2)

c      write(ifn,810)srt,db2,dbfit,nhlf
c      write(6,810)srt,db2,dbfit,nhlf
       write(ifn,810)srt,dbfit

      end do
2000  continue
1000  continue
810   format(f9.3,1x,2(1x,g13.6),1x,2i5)

      end

c*******************************************************
      subroutine integpar

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      character ch1*1,ckf1*8,ckf2*8
      dimension aa(7,250),ia(250)
      character achar(7)*23
      data achar/
     $ 'data (aa(1,i),i=1,231)/',
     $ 'data (aa(2,i),i=1,231)/',
     $ 'data (aa(3,i),i=1,231)/',
     $ 'data (aa(4,i),i=1,231)/',
     $ 'data (aa(5,i),i=1,231)/',
     $ 'data (aa(6,i),i=1,231)/',
     $ 'data (aa(7,i),i=1,231)/'/

      do i=1,250
	ia(i)=0
      do j=1,7
	aa(j,i)=0.0d0
      end do
      end do
      open(10,file='data',status='old')

      i=0
      icl=0
      do il=1,3000
	i=i+1
      read(10,800,err=999,end=900)
     $    i1,i2,kf1,ckf1,ch1,kf2,ckf2,a0,a1,a2,a3,a4,b0,b1

c     write(6,800)
c    $    i1,i2,kf1,ckf1,ch1,kf2,ckf2,a0,a1,a2,a3,a4,b0,b1

       kc1=jamcomp(kf1)
       kc2=jamcomp(kf2)
       if(kc1.le.0.or.kc2.le.0) then
	 write(6,*)'i kf1 kf2 kc1 kc2',i,kf1,kf2,kc1,kc2
	 stop
       endif
       jd1=idget1(kc1)
       jd2=idget1(kc2)
       ic=jamcpair(jd1,jd2)
c      print *,'jd1 jd2 ic',ic,jd1,jd2,ckf1,ckf2

       if(ic.ge.1.and.ic.le.250) then
	  icl=icl+1
	  ia(ic)= icl
	  aa(1,ic)=a0
	  aa(2,ic)=a1
	  aa(3,ic)=a2
	  aa(4,ic)=a3
	  aa(5,ic)=a4
	  aa(6,ic)=b0
	  aa(7,ic)=b1
       else
	 print *,'funny ic',ic
	 stop
       endif

      end do

800   format(2x,i2,1x,i2,1x,i6,1x,a8,a1,1x,i6,1x,a8,7(2x,f9.6))
999   continue
      write(6,*)'error',i
900   continue

      write(6,*)icl
      do i=1,250
	write(6,*)i,ia(i),aa(1,i)
      end do



      do ii=1,7
      write(11,810)ii,achar(ii)
810   format('c...a',i2,/,6x,a23)
      do i=1,32
      write(11,811)(aa(ii,j),j=7*(i-1)+1,i*7)
      end do
      i=33
      write(11,812)(aa(ii,j),j=7*(i-1)+1,i*7)
      end do
811   format(5x,'$',7(f8.4,','))
812   format(5x,'$',6(f8.4,','),f8.4,'/')


      write(11,820)
820   format(/'c...ia',/,6x,'data ia/')
      do i=1,32
      write(11,821)(ia(j),j=7*(i-1)+1,i*7)
      end do
      i=33
      write(11,822)(ia(j),j=7*(i-1)+1,i*7)
821   format(5x,'$',7(i4,','))
822   format(5x,'$',6(i4,','),i4,'/')


c...mstc(22) : (D=391) offset of KC code for N*
c...mstc(23) : (D=409) offset of KC code for N*
c...mstc(24) : (D=411) starting piont of KC code for D*
c...mstc(25) : (D=443) end point of KC code for D*

      end

c*******************************************************************

c     subroutine bwres_t(srt,idd1,idd2,kf1,kf2,kc1,kc2)
      subroutine bwres_t

c...Calculate individual resonance production probability.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      character chekc*80

      srt=3.0d0
      idd1=2
      idd2=id_nucl
      if(idd1.ge.1.and.idd1.le.2) then         ! N*
        kcmin1=mstc(22)+idd1-1
        kcmax1=mstc(23)+idd1-1
        isp1=2
      else if(idd1.ge.3.and.idd1.le.6) then    ! D*
        kcmin1=mstc(24)+idd1-3
        kcmax1=mstc(25)+idd1-3
        isp1=4
      else
        id1=kchg(kc1,5)
        if(id1.eq.id_nucl.or.id1.eq.id_delt) then
    	  isp1=1
          kcmin1=jamcomp(kf1)
          kcmax1=kcmin1
        else if(id1.eq.id_nucls) then
	  isp1=1
          kcmin1=jamcomp(kf1)
          kcmax1=kcmin1
        else
          write(chekc,*)idd1
          call jamerrm(1,0,'(bwres:)1 invalid idd ='//chekc)
        endif
      endif

      if(idd2.ge.1.and.idd2.le.2) then         ! N*
        kcmin2=mstc(22)+idd2-1
        kcmax2=mstc(23)+idd2-1
        isp2=2
      else if(idd2.ge.3.and.idd2.le.6) then    ! D*
        kcmin2=mstc(24)+idd2-3
        kcmax2=mstc(25)+idd2-3
        isp2=4
      else
        id2=kchg(kc2,5)
c        print *,'id2 kc2 kf2',id2,kc2,kf2
        if(id2.eq.id_nucl.or.id2.eq.id_delt) then
	  isp2=1
          kcmin2=jamcomp(kf2)
          kcmax2=kcmin2
        else if(id2.eq.id_nucls) then
	  isp2=1
          kcmin2=jamcomp(kf2)
          kcmax2=kcmin2
        else
          write(chekc,*)idd2
          call jamerrm(1,0,'(bwres:)2 invalid idd ='//chekc)
        endif
      endif

      bwtot=0.0d0
      do i=kcmin1,kcmax1,isp1
      do j=kcmin2,kcmax2,isp2
        spin1=mod(kchg(i,4),10)
        spin2=mod(kchg(j,4),10)
	bwtot=bwtot+jambwtbl(srt,i,j)*spin1*spin2
      end do
      end do

      do i=kcmin1,kcmax1,isp1
      do j=kcmin2,kcmax2,isp2
        spin1=mod(kchg(i,4),10)
        spin2=mod(kchg(j,4),10)
        write(6,*)jambwtbl(srt,i,j)*spin1*spin2/bwtot
      end do
      end do
      kc1=kcmax1
      kc2=kcmax2
 300  continue
      kf1=kchg(kc1,4)
      kf2=kchg(kc2,4)

      end

c*************************************************

      function idget1(kc1)

c...Get id of N* and D*
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

       id1=kchg(kc1,5)
       if(id1.eq.id_nucl) then
	 jd1=1
       else if(id1.eq.id_delt) then
	 jd1=2
       else if(id1.eq.id_nucls) then
	 jd1=3+(kc1-mstc(22))/2
       else if(id1.eq.id_delts) then
c        jd1=4+(mstc(23)-mstc(22))/2+(kc1-mstc(24))/4
	 jd1=4+(mstc(23)-mstc(22))/2+(kc1-mstc(24))/4
       else
	write(6,*)'error kc1=',kc1
	stop
       endif

       idget1=jd1

       end

c********************************************************************
      subroutine jetset2_t

C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/pjjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/pjdat1/mstu(200),paru(200),mstj(200),parj(200)
      common/pjdat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/pjdat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
C     save /pjjets/,/pjdat1/,/pjdat2/,/pjdat3/,/pjsubs/,/pjpars/
      common/jampos1/jqconst(2),kfcq(4),icq(4),icms
      save /jampos1/

      icms=0
      jqconst(1)=1
      jqconst(2)=1

c     call init
c...Choice of baryon production model.
      mstj(12)=2
c...Form of particle decay.
      mstj(21)=0
c....Hadron formation point.
      mstj(10)=3

c...Baryon resonance production parameter.
      parj(27)=0.0d0


      ymin=-2.0D0
      ymax=2.0D0
      wy=0.1D0
      nymx=(ymax-ymin)/wy
      call vbook1(11,'dN/dy proton',nymx,ymin,ymax)
      call vbook1(12,'dN/dy pi+',nymx,ymin,ymax)
      em=3.15d0

      nev=100
      do iev=1,nev
        call pj2ent(0,1,2101,em)
        call pjexec
c       call pjlist(1)

        do  i=1,n
c......This is already removed particle.
        if(k(i,1).ge.10.or.k(i,1).eq.0) go to 390
         kf=k(i,2)
         rap=0.5D0*log( max(p(i,4)+p(i,3),1.D-8)/max(p(i,4)-p(i,3), 
     & 1.D-8) )


          if(kf.eq.2212) call vfill1(11,rap,1.D0/wy)
          if(kf.eq.221) call vfill1(12,rap,1.D0/wy)
 390    end do

c....Loop over event
      end do

c...Event weight
      fac=1.d0/dble(nev)

c...Mass distributions.
      do i=1,2
       call vscale(10+i,fac)
       call vprint(10+i,0,0)
      end do

      end

c*******************************************************************
      subroutine init

      implicit double precision(a-h, o-z)
      common/pjdat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)

        mdcy(jamcomp(211),1)=0  ! pi+
        mdcy(jamcomp(-211),1)=0 ! pi-
        mdcy(jamcomp(111),1)=0  ! pi0
        mdcy(jamcomp(311),1)=0  ! k0
        mdcy(jamcomp(-311),1)=0 ! ak0
        mdcy(jamcomp(321),1)=0  ! k+
        mdcy(jamcomp(-321),1)=0 ! k-
        mdcy(jamcomp(411),1)=0    ! D+
        mdcy(jamcomp(-411),1)=0   ! D-
        mdcy(jamcomp(421),1)=0    ! D0
        mdcy(jamcomp(-421),1)=0   ! aD0
        mdcy(jamcomp(221),1)=0    ! eta
        mdcy(jamcomp(331),1)=0    ! eta'
        mdcy(jamcomp(441),1)=0    ! eta_c

        mdcy(jamcomp(310),1)=0
        mdcy(jamcomp(431),1)=0
        mdcy(jamcomp(-431),1)=0
        mdcy(jamcomp(511),1)=0
        mdcy(jamcomp(-511),1)=0
        mdcy(jamcomp(521),1)=0
        mdcy(jamcomp(-521),1)=0
        mdcy(jamcomp(531),1)=0
        mdcy(jamcomp(-531),1)=0
        mdcy(jamcomp(3122),1)=0
        mdcy(jamcomp(-3122),1)=0
        mdcy(jamcomp(3112),1)=0
        mdcy(jamcomp(-3112),1)=0
        mdcy(jamcomp(3212),1)=0
        mdcy(jamcomp(-3212),1)=0
        mdcy(jamcomp(3222),1)=0
        mdcy(jamcomp(-3222),1)=0
        mdcy(jamcomp(3312),1)=0
        mdcy(jamcomp(-3312),1)=0
        mdcy(jamcomp(3322),1)=0
        mdcy(jamcomp(-3322),1)=0
        mdcy(jamcomp(3334),1)=0
        mdcy(jamcomp(-3334),1)=0

      end


c*******************************************************************

      subroutine code_tex

      implicit double precision(a-h, o-z)
      include 'jam2.inc'


c...mstc(22)  offset of KC code for N*
c...mstc(23)  offset of KC code for N*
c...mstc(24)  starting piont of KC code for D*
c...mstc(25)  end point of KC code for D*
c...mstc(26)  offset of KC code for Lambda*
c...mstc(27)  offset of KC code for Lambda*
c...mstc(28)  offset of KC code for Sigma*
c...mstc(29)  offset of KC code for Sigma*
c...mstc(30)  offset of KC code for Xi*
c...mstc(31)  offset of KC code for Xi*
c...mstc(32)  offset of KC code for mesons
c...mstc(33)  offset of KC code for mesons

      character header(8)*80,endcom(3)*14
      character ckf*9
      character cnuc(80)*70,cdelt(80)*70,clambda(80)*70,csigma(80)*70,
     $  cxi(80)*70
      character cmeson1(80)*70,cmeson2(80)*70,cmeson3(80)*70
      character cmes*17,chafm*16
      character charb(80)*70,charm(80)*70

      data header/
     $ '\\begin{table}[ptb]',
     $ '\\captive{','\\protect\\label{t:codefive} }',
     $ '\\vspace{1ex}',
     $'\\begin{center}',
     $'\\begin{tabular}{|c|c|c||c|c|c|',
     $'@{\\protect\\rule{0mm}{\\tablinsep}}}',
     $'KF & Name & Printed & KF & Name & Printed'/
      data endcom/ '\\end{tabular}', '\\end{center}','\\end{table}'/

      do i=1,80
        cnuc(i)='          &              &'
      end do
      do i=1,80
        cdelt(i)='          &              &'
      end do

c...N and N*
      cnuc(1)='2112 & n & n0'
      cnuc(2)='2212 & p & p+'
      nn=2
      do i=mstc(22),mstc(23)+1
        nn=nn+1
        write(ckf(1:9),'(i9)')kchg(i,4)
        cnuc(nn)=ckf//' & '//'$'//chaf(i,1)(1:1)//chaf(i,1)(2:7)//
     $    '^'//chaf(i,1)(8:8)//'$ &  \\ttt{'//chaf(i,1)(1:8)//'}'
      end do

c...Delta
      cdelt(1)='1114 & $\\Delta^-$  & \\ttt{delta-}'
      cdelt(2)='2114 & $\\Delta^0$  & \\ttt{delta0}'
      cdelt(3)='2214 & $\\Delta^+$  & \\ttt{delta+}'
      cdelt(4)='2224 & $\\Delta^++$ & \\ttt{delta++}'
      nd=4
      do i=mstc(24),mstc(25)+3
        nd=nd+1
        write(ckf(1:9),'(i9)')kchg(i,4)
        cdelt(nd)=ckf//' & '//'$\\Delta'//chaf(i,1)(2:7)//
     $    '^{'//chaf(i,1)(8:9)//'}$ &  \\ttt{'//chaf(i,1)(1:9)//'}'
      end do

c...Write header.
      write(6,8000)(header(i),i=1,8)
      mm=max(nn,nd)
      do i=1,mm
       write(6,8100)cnuc(i),cdelt(i)
      end do

      write(6,8900)
      write(6,9000)(endcom(i),i=1,3)


c------------------------------------------------------------------------
c...Strange baryons.
c------------------------------------------------------------------------

      do i=1,80
        clambda(i)='          &              &'
      end do
      do i=1,80
        csigma(i)='          &              &'
      end do
      do i=1,80
        cxi(i)='          &              &'
      end do

c...Lambda.
      clambda(1)='3122 & $\\Lambda^0$ & \\ttt{Lambda0}'
      nl=1
      do i=mstc(26),mstc(27)
        nl=nl+1
        write(ckf(1:9),'(i9)')kchg(i,4)
        clambda(nl)=ckf//' & '//'$\\Lambda'//chaf(i,1)(2:7)//
     $    '^'//chaf(i,1)(8:8)//'$ &  \\ttt{'//chaf(i,1)(1:8)//'}'
      end do

c...Sigma.
      csigma(1)='3112 & $\\Sigma^-$     & \\ttt{Sigma-}'
      csigma(2)='3112 & $\\Sigma^0$     & \\ttt{Sigma0}'
      csigma(3)='3212 & $\\Sigma^+$     & \\ttt{Sigma+}'
      csigma(4)='3114 & $\\Sigma^{*-}$  & \\ttt{Sigma*-}'
      csigma(5)='3114 & $\\Sigma^{*0}$  & \\ttt{Sigma*0}'
      csigma(6)='3214 & $\\Sigma^{*+}$  & \\ttt{Sigma*+}'
      ns=6
      do i=mstc(28)+3,mstc(29)+2
        ns=ns+1
        write(ckf(1:9),'(i9)')kchg(i,4)
        csigma(ns)=ckf//' & '//'$\\Sigma'//chaf(i,1)(2:7)//
     $    '^'//chaf(i,1)(8:8)//'$ &  \\ttt{'//chaf(i,1)(1:9)//'}'
      end do

c...Xi.
      cxi(nl+2)='3312 & $\\Xi^-$     & \\ttt{Xi-}'
      cxi(nl+3)='3322 & $\\Xi^0$     & \\ttt{Xi0}'
      cxi(nl+4)='3312 & $\\Xi^{*-}$  & \\ttt{Xi*-}'
      cxi(nl+5)='3314 & $\\Xi^{*0}$  & \\ttt{Xi*0}'
      nx=nl+5
      do i=mstc(30)+2,mstc(31)+1
        nx=nx+1
        write(ckf(1:9),'(i9)')kchg(i,4)
        cxi(nx)=ckf//' & '//'$\\Xi'//chaf(i,1)(2:7)//
     $    '^'//chaf(i,1)(8:8)//'$ &  \\ttt{'//chaf(i,1)(1:9)//'}'
      end do
      nx=nx+2
      cxi(nx)='3334 & $\\Omega^-$  & \\ttt{Omega-}'

c...Write header.
      write(6,8010)(header(i),i=1,8)
      mm=min(nl,ns)
      do i=1,mm
       write(6,8200)clambda(i),csigma(i)
      end do
      mmx=max(nx,ns)
      do i=mm+1,mmx
       write(6,8200)cxi(i),csigma(i)
      end do

      write(6,8900)
      write(6,9000)(endcom(i),i=1,3)

8000  format(a19,/,a10,'Non-strange baryons',a30,/,a13,/,a15,/,a31,a37,
     $  '\\hline',/,a41,'\\\\ \\hline')
8010  format(a19,/,a10,'Strange baryons',a30,/,a13,/,a15,/,a31,a37,
     $  '\\hline',/,a41,'\\\\ \\hline')

8100  format(a45,' & ',a50,' \\\\')
8200  format(a50,' & ',a50,' \\\\')
8900  format('\\hline')
9000  format(a14,/,a13,/,a12)


c--------------------------------------------------------------------------
c...Meson codes.
c--------------------------------------------------------------------------
c...Loop over mesons.
      imes1=0
      imes2=0
      imes3=0
      imes4=0
      imes5=0
      do kc=mstc(32),mstc(33)
        id=kchg(kc,5)
        if(id.eq.id_pi.or.id.eq.id_light1) then
          if(imes1.eq.0) imes1=kc
          nmes1=kc
        else if(id.eq.id_light0) then
          if(imes2.eq.0) imes2=kc
          nmes2=kc
        else if(id.eq.id_str) then
          if(imes3.eq.0) imes3=kc
          nmes3=kc
        else if(id.eq.id_charm) then
          if(imes4.eq.0) imes4=kc
          nmes4=kc
        else
          if(imes5.eq.0) imes5=kc
          nmes5=kc
        endif
      end do

       print *,'imes1 nmes1',imes1,nmes1
       print *,'imes2 nmes2',imes2,nmes2
       print *,'imes3 nmes3',imes3,nmes3
 
      do i=1,80
      cmeson1(i)='          &               &'
      cmeson2(i)='          &               &'
      cmeson3(i)='          &               &'
      end do

c...Light meson t=1.
      nm1=0
      do i=imes1,nmes1
        nm1=nm1+1
        write(ckf(1:9),'(i9)')kchg(i,4)
        ml=0
        mu=0
        do il=1,16
         if(chaf(i,1)(il:il).eq.' ') then
            ml=il-1
            goto 100
         else if(chaf(i,1)(il:il).eq.'_') then
           mu=il
         endif
        end do
 100    continue
        if(chaf(i,1)(1:3).eq.'rho') then
          cmes='$\\rho'//chaf(i,1)(4:ml-1)//'^'//chaf(i,1)(ml:ml)//'$'
        else if(chaf(i,1)(1:2).eq.'pi') then
          cmes='$\\pi'//chaf(i,1)(3:ml-1)//'^'//chaf(i,1)(ml:ml)//'$'
        else if(chaf(i,1)(1:1).eq.'a') then
          cmes='a$'//chaf(i,1)(2:mu+1)//'^'//chaf(i,1)(ml:ml)//'$'
        else if(chaf(i,1)(1:1).eq.'b') then
          cmes='b$'//chaf(i,1)(2:mu+1)//'^'//chaf(i,1)(ml:ml)//'$'
        endif
        chafm=chaf(i,1)
        if(mu.ne.0) chafm=chaf(i,1)(1:mu-1)//'\\_'//chaf(i,1)(mu+1:9)

        cmeson1(nm1)=ckf//' & '//cmes//' &  \\ttt{'//chafm(1:9)//'}'
      end do

c...Light meson t=0.
      nm2=0
      do i=imes2,nmes2
        nm2=nm2+1
        write(ckf(1:9),'(i9)')kchg(i,4)
        ml=0
        mu=0
        do il=1,16
         if(chaf(i,1)(il:il).eq.' ') then
            ml=il-1
            goto 101
         else if(chaf(i,1)(il:il).eq.'_') then
           mu=il
         endif
        end do
 101    continue
        if(chaf(i,1)(1:5).eq.'omega') then
          cmes='$\\'//chaf(i,1)(1:ml)//'$'
        else if(chaf(i,1)(1:3).eq.'eta') then
          cmes='$\\'//chaf(i,1)(1:ml)//'$'
        else if(chaf(i,1)(1:3).eq.'phi') then
          cmes='$\\'//chaf(i,1)(1:ml)//'$'
        else if(chaf(i,1)(1:5).eq.'sigma') then
          cmes='$\\sigma$'
        else if(chaf(i,1)(1:1).eq.'h') then
          cmes='h$'//chaf(i,1)(2:ml)//'$'
        else if(chaf(i,1)(1:1).eq.'f') then
          cmes='f$'//chaf(i,1)(2:ml)//'$'
        endif
        chafm=chaf(i,1)
        if(mu.ne.0) chafm=chaf(i,1)(1:mu-1)//'\\_'//chaf(i,1)(mu+1:12)

        cmeson2(nm2)=ckf//' & '//cmes//' &  \\ttt{'//chafm(1:12)//'}'
      end do

c...Strange mesons.
      nm3=0
      do i=imes3,nmes3
        nm3=nm3+1
        write(ckf(1:9),'(i9)')kchg(i,4)
        ml=0
        mu=0
        do il=1,16
         if(chaf(i,1)(il:il).eq.' ') then
            ml=il-1
            goto 102
         else if(chaf(i,1)(il:il).eq.'_') then
           mu=il
         endif
        end do
 102    continue
        cmes='                 '
        if(kchg(i,4).eq.310) then
          cmes='K$_{S}^{0}$'
        else if(kchg(i,4).eq.130) then
          cmes='K$_{L}^{0}$'
        else if(kchg(i,4).eq.311) then
          cmes='K$^{0}$'
        else if(kchg(i,4).eq.321) then
          cmes='K$^{+}$'
        else if(kchg(i,4).eq.313) then
          cmes='K$^{*0}$'
        else if(kchg(i,4).eq.323) then
          cmes='K$^{*+}$'
        else if(chaf(i,1)(1:3).eq.'K*_') then
          cmes='K$^{*'//chaf(i,1)(ml:ml)//'}'//chaf(i,1)(3:4)//'$'
        else if(chaf(i,1)(1:3).eq.'K*(') then
          cmes='K'//chaf(i,1)(3:8)//'$^{'//chaf(i,1)(ml:ml)//'}$'
        else if(chaf(i,1)(1:2).eq.'K_') then
          cmes='K$'//chaf(i,1)(2:ml-1)//'^{'//chaf(i,1)(ml:ml)//'}$'
        endif
        chafm=chaf(i,1)
        if(mu.ne.0) chafm=chaf(i,1)(1:mu-1)//'\\_'//chaf(i,1)(mu+1:12)
        cmeson3(nm3)=ckf//' & '//cmes//' &  \\ttt{'//chafm(1:10)//'}'
      end do

c...Write strange mesons.
      write(6,800)
      mm=max(nm1,nm2,nm3)
      do i=1,mm
        write(6,810)cmeson1(i),cmeson2(i),cmeson3(i)
      end do
      write(6,8900)
      write(6,900)

 800  format('\\begin{table}[ptb]',/,
     $ '\\captive{light and strange mesons',
     $ '\\protect\\label{t:codefive} }',/,
     $ '\\vspace{1ex}',/,
     $ '\\begin{center}',/,
     $ '\\begin{tabular}{|c|c|c||c|c|c||c|c|c|',
     $ '@{\\protect\\rule{0mm}{\\tablinsep}}} \\hline',/,
     $ 'KF & Name & Printed & KF & Name & Printed',
     $ ' & KF & Name & Printed \\\\ \\hline')

 810  format(a50,' & ',a53,' & ',a50,' \\\\')
 900  format('\\end{tabular}'/'\\end{center}'/'\\end{table}')



c----------------------------------------------------------------------
c...Heavy mesons and hadrons.
c----------------------------------------------------------------------
      do i=1,80
        charb(i)='          &              &'
        charm(i)='          &              &'
      end do
      im=0
      ib=0
      do kc=101,mstu(6)
       id=kchg(kc,5)
       if(kchg(kc,4).eq.440) goto 111
       itag=0
c...Baryons.
       if(id.eq.id_charmb .or. id.eq.id_bottb) then
         ib=ib+1
         itag=1
c...Mesons.
       else if(id.eq.id_charm .or. id.eq.id_bott
     $   .or.id.eq.id_cc .or.id.eq.id_bb) then
         im=im+1
         itag=2
       endif
       if(itag.ne.0)then
        write(ckf(1:9),'(i9)')kchg(kc,4)
        ml=0
        mu=0
        do il=1,16
         if(chaf(kc,1)(il:il).eq.' ') then
            ml=il-1
            goto 104
         else if(chaf(kc,1)(il:il).eq.'_') then
           mu=il
         endif
        end do
 104    continue
        ll=ml
        chafm(1:ml)=chaf(kc,1)(1:ml)
        if(mu.ne.0) then
          chafm=chaf(kc,1)(1:mu-1)//'\\_'//chaf(kc,1)(mu+1:ml)
          ll=ml+1
        endif
        if(itag.eq.1) then
          charb(ib)=ckf//' & $'//chaf(kc,1)(1:ml)//'$ & '//
     $      '\\ttt{'//chafm(1:ll)//'}'
        else if(itag.eq.2) then
          charm(im)=ckf//' & $'//chaf(kc,1)(1:ml)//'$ & '//
     $      '\\ttt{'//chafm(1:ll)//'}'
        endif

       endif
 111  end do

      write(6,7000)
      mm=max(ib,im)
      do i=1,mm
        write(6,7200)charm(i),charb(i)
      end do
      write(6,8900)
      write(6,900)

7000  format('\\begin{table}[ptb]',/,
     $ '\\captive{Heavy mesons and baryons',
     $ '\\protect\\label{t:codefive} }',/,
     $ '\\vspace{1ex}',/,
     $ '\\begin{center}',/,
     $ '\\begin{tabular}{|c|c|c||c|c|c|',
     $ '@{\\protect\\rule{0mm}{\\tablinsep}}} \\hline',/,
     $ 'KF & Name & Printed',
     $ ' & KF & Name & Printed \\\\ \\hline')
7200  format(a40,' & ',a48,' \\\\')

      end

c*******************************************************************

      subroutine test_res


      include 'jam1.inc'
      include 'jam2.inc'

c...Commonblock for t-channel resonance productions.
      common/jamres1/kfo(2,20),noutpa

      parameter(mxchan=30)
      dimension sigin(mxchan)
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a)

      noutpa=0
      nv=2
      nbary=1
      nmeson=1
      call jamzero(1)
      call jamzero(2)
      kf1=311
      kf2=4214
      kc1=jamcomp(kf1)
      kc2=jamcomp(kf2)
      em1=pmas(kc1,1)
c     em2=pjmass(kf2)
      em2=2.9
      srt=4.5d0
      ibar1=isign(kchg(jamcomp(kf1),6),kf1)
      ibar2=isign(kchg(jamcomp(kf2),6),kf2)
      pr=pawt(srt,em1,em2)
      k(9,1)=ibar1
      k(9,2)=ibar2
      k(1,1)=1
      k(1,2)=1
      p(3,1)=-pr
      p(3,2)=-pr
      p(5,1)=em1
      p(5,2)=em2

      icltyp=jamcltyp(kf1,kf2,ibar1,ibar2)
      call jamcross(1,icltyp,srt,pr,kf1,kf2,em1,em2,
     $                 sig,sigel,sigin,mchanel,mabsrb,ijet,icon)

      pare(3)=(sig-sigel)*rn(0)
      call jamcross(3,icltyp,srt,pr,kf1,kf2,em1,em2,
     $                 sig,sigel,sigin,mchanel,mabsrb,ijet,icon)

      print *,'em1 em2',em1,em2
      call jamrmas2(kf1,kf2,kc1,kc2,srt,em1,em2,icon)
      print *,'icon',icon
      print *,'kf1 em1',kf1,em1
      print *,'kf2 em2',kf2,em2

      end

