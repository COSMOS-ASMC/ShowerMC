c....Analysis program produced for NA49

      include 'jam1.inc'
      include 'jam2.inc'
      common/myana/ispec

      open(33,file='phase.dat',form='unformatted',status='old')

c     read(33)nevent,eylab,ebeta,egamma,icm
c     mstc(2)=nevent
c     pard(17)=eylab
c     mstc(4)=icm

      read(33)nevent,ylab,beta,gamma,icm
      pard(5)=beta
      pard(6)=gamma
      mstc(2)=nevent
      pard(17)=ylab
      mstc(4)=icm

      print *,'nevent ylab beta gamma icm',nevent,ylab,beta,gamma,icm

      call ana1

      mevt=0
      do ie=1,1000000

        read(33,end=900,err=999)iev,nv,nbary,nmeson,b

c       read(33,end=900,err=999)iev,nv,pard(2)
        b=0.0d0

        mevt=mevt+1
	pard(2)=b
        nbary=0
        do i=1,nv
c         read(33)(k(j,i),j=1,11),(r(j,i),j=1,5),(p(j,i),j=1,5)
c         read(33)(k(j,i),j=1,9),(r(j,i),j=1,5),(p(j,i),j=1,5)

          read(33)(k(j,i),j=1,11),(r(j,i),j=1,5),(p(j,i),j=1,5)
     $              ,(v(j,i),j=1,5)

	   if(k(2,i).eq.3122.and.p(5,i).le.1.0) then
	     write(50,*)'lam???',ie,k(1,i),k(2,i),p(5,i)
	   endif
	   if(k(2,i).eq.2212.and.p(5,i).ge.1.0) then
	     write(50,*)'prot???',ie,k(1,i),k(2,i),p(5,i)
	   endif
          if(abs(k(9,i)).eq.3) nbary=nbary+1
        end do

c       call jamdeut
c       call jamlist(1)

        call ana2
        if(mod(ie,10).eq.0) write(6,*)iev,nv,b,ispec,mstd(92)

      end do
900   continue
      close(33)
      call ana3(mevt)
      stop
999   continue
      write(6,*)'error in reading file'
      end

C***********************************************************************

      subroutine ana1

C...Initial values and definitions for histograms.
      include 'jam1.inc'
      include 'jam2.inc'
      logical first
      dimension nyp(16),nyh(21),ptp(16),pth(21)
      save nyp,nyh,ptp,pth
      save wevt
      save ymin,ymax,wy,nymx
      save pmin,pmax,wp,npmx
      save emin,emax,we,nemx
      dimension nn(-500:500)
      save nn
      save first
      save ylab,beta,gamma
      data first/.true./
C.....Limits and bins for rapidity distributions.
      data ymin,ymax,wy/-3d0,7d0,0.3d0/
C.....Limits and bins for pT-distributions.
      data pmin,pmax,wp/0d0,5d0,0.1d0/
C.....Limits and bins for energy-momentum event fluctuations.
      data emin,emax,we/-3d0,7d0,0.15d0/

      if(first) then
        beta=pard(5)
        gamma=pard(6)
        do i=1,16
          nyp(i)=0
          ptp(i)=0d0
        end do 
        do i=1,20
          nyh(i)=0
          pth(i)=0d0
        end do 
        first=.false.
      endif 

      do i=-500,500
        nn(i)=0
      end do

C.....Event weight.
c     wevt=1d0/dble(mstc(2))
      wevt=1d0
      nymx=(ymax-ymin)/wy
      npmx=(pmax-pmin)/wp
      nemx=(emax-emin)/we
      ylab=pard(17)
      ylab=2.9

C...Initialize histogram booking.

c...rapidity
      call vbook1(1,'dN/dy negative(pi-k-pbar)',nymx,ymin,ymax)
      call vbook1(8, 'dE/dy     ',nemx,emin,emax)

c...CERN/SPS for Pb+Pb P.R.L 82(1999)2471 NA49 5%
      call vbook1(11,'1/ptdN/dpt 2.2<y<2.4 protons',npmx,pmin,pmax)
      call vbook1(12,'1/ptdN/dpt 2.4<y<2.6 protons',npmx,pmin,pmax)
      call vbook1(13,'1/ptdN/dpt 2.6<y<2.8 protons',npmx,pmin,pmax)
      call vbook1(14,'1/ptdN/dpt 2.8<y<3.0 protons',npmx,pmin,pmax)
      call vbook1(15,'1/ptdN/dpt 3.0<y<3.2 protons',npmx,pmin,pmax)
      call vbook1(16,'1/ptdN/dpt 3.2<y<3.4 protons',npmx,pmin,pmax)
      call vbook1(17,'1/ptdN/dpt 3.4<y<3.6 protons',npmx,pmin,pmax)
      call vbook1(18,'1/ptdN/dpt 3.6<y<3.8 protons',npmx,pmin,pmax)
      call vbook1(19,'1/ptdN/dpt 3.8<y<4.0 protons',npmx,pmin,pmax)
      call vbook1(20,'1/ptdN/dpt 4.0<y<4.2 protons',npmx,pmin,pmax)
      call vbook1(21,'1/ptdN/dpt 4.2<y<4.4 protons',npmx,pmin,pmax)
      call vbook1(22,'1/ptdN/dpt 4.4<y<4.6 protons',npmx,pmin,pmax)
      call vbook1(23,'1/ptdN/dpt 4.6<y<4.8 protons',npmx,pmin,pmax)
      call vbook1(24,'1/ptdN/dpt 4.8<y<5.0 protons',npmx,pmin,pmax)
      call vbook1(25,'1/ptdN/dpt 5.0<y<5.2 protons',npmx,pmin,pmax)
      call vbook1(26,'1/ptdN/dpt 5.2<y<5.4 protons',npmx,pmin,pmax)

      call vbook1(27,'1/ptdN/dpt 2.9<y<3.1 h-',npmx,pmin,pmax)
      call vbook1(28,'1/ptdN/dpt 3.1<y<3.3 h-',npmx,pmin,pmax)
      call vbook1(29,'1/ptdN/dpt 3.3<y<3.5 h-',npmx,pmin,pmax)
      call vbook1(30,'1/ptdN/dpt 3.5<y<3.7 h-',npmx,pmin,pmax)
      call vbook1(31,'1/ptdN/dpt 3.7<y<3.9 h-',npmx,pmin,pmax)
      call vbook1(32,'1/ptdN/dpt 3.9<y<4.1 h-',npmx,pmin,pmax)
      call vbook1(33,'1/ptdN/dpt 4.1<y<4.3 h-',npmx,pmin,pmax)
      call vbook1(34,'1/ptdN/dpt 4.3<y<4.5 h-',npmx,pmin,pmax)
      call vbook1(35,'1/ptdN/dpt 4.5<y<4.7 h-',npmx,pmin,pmax)
      call vbook1(36,'1/ptdN/dpt 4.7<y<4.9 h-',npmx,pmin,pmax)
      call vbook1(37,'1/ptdN/dpt 4.9<y<5.1 h-',npmx,pmin,pmax)
      call vbook1(38,'1/ptdN/dpt 5.1<y<5.3 h-',npmx,pmin,pmax)
      call vbook1(39,'1/ptdN/dpt 5.3<y<5.5 h-',npmx,pmin,pmax)
      call vbook1(40,'1/ptdN/dpt 5.5<y<5.7 h-',npmx,pmin,pmax)
      call vbook1(41,'1/ptdN/dpt 5.7<y<5.9 h-',npmx,pmin,pmax)
      call vbook1(42,'1/ptdN/dpt 5.9<y<6.1 h-',npmx,pmin,pmax)
      call vbook1(43,'1/ptdN/dpt 6.1<y<6.3 h-',npmx,pmin,pmax)
      call vbook1(44,'1/ptdN/dpt 6.3<y<6.5 h-',npmx,pmin,pmax)
      call vbook1(45,'1/ptdN/dpt 6.5<y<6.7 h-',npmx,pmin,pmax)
      call vbook1(46,'1/ptdN/dpt 6.7<y<6.9 h-',npmx,pmin,pmax)

C.....rapidity spectra of hadrons.
      call vbook1(51,'dN/dy net baryon',nymx,ymin,ymax)
      call vbook1(52,'dN/dy net proton',nymx,ymin,ymax)
      call vbook1(53,'dN/dy net lambda',nymx,ymin,ymax)
      call vbook1(54,'dN/dy net sigma-',nymx,ymin,ymax)
      call vbook1(55,'dN/dy net sigma0',nymx,ymin,ymax)
      call vbook1(56,'dN/dy net sigma+',nymx,ymin,ymax)
      call vbook1(57,'dN/dy net xi-',nymx,ymin,ymax)
      call vbook1(58,'dN/dy net xi0',nymx,ymin,ymax)
      call vbook1(59,'dN/dy net lam+sig0',nymx,ymin,ymax)

c...baryons
      call vbook1(61,'dN/dy total baryon',nymx,ymin,ymax)
      call vbook1(62,'dN/dy total proton',nymx,ymin,ymax)
      call vbook1(63,'dN/dy total lambda',nymx,ymin,ymax)
      call vbook1(64,'dN/dy total sigma-',nymx,ymin,ymax)
      call vbook1(65,'dN/dy total sigma0',nymx,ymin,ymax)
      call vbook1(66,'dN/dy total sigma+',nymx,ymin,ymax)
      call vbook1(67,'dN/dy total xi-',nymx,ymin,ymax)
      call vbook1(68,'dN/dy total xi0',nymx,ymin,ymax)
      call vbook1(69,'dN/dy total lam+sig0',nymx,ymin,ymax)

c....anti-baryons
      call vbook1(71,'dN/dy total anti-baryon',nymx,ymin,ymax)
      call vbook1(72,'dN/dy total anti-proton',nymx,ymin,ymax)
      call vbook1(73,'dN/dy total anti-lambda',nymx,ymin,ymax)
      call vbook1(74,'dN/dy total anti-sigma-',nymx,ymin,ymax)
      call vbook1(75,'dN/dy total anti-sigma0',nymx,ymin,ymax)
      call vbook1(76,'dN/dy total anti-sigma+',nymx,ymin,ymax)
      call vbook1(77,'dN/dy total anti-xi-',nymx,ymin,ymax)
      call vbook1(78,'dN/dy total anti-xi0',nymx,ymin,ymax)
      call vbook1(79,'dN/dy anti-lam+sig0',nymx,ymin,ymax)

c....Mesons.
      call vbook1(81,'dN/dy pi-',nymx,ymin,ymax)
      call vbook1(82,'dN/dy pi0',nymx,ymin,ymax)
      call vbook1(83,'dN/dy pi+',nymx,ymin,ymax)
      call vbook1(84,'dN/dy Kaon-',nymx,ymin,ymax)
      call vbook1(85,'dN/dy Kaon0/aKaon0',nymx,ymin,ymax)
      call vbook1(86,'dN/dy Kaon+',nymx,ymin,ymax)


      return

C***********************************************************************

      entry ana2

C...Do analysis of particle record.
      do 200 i=1,nv

C.......Skip all "dead" particles, photons, leptons, neutrinos, etc.
        if(k(1,i).le.0.or.k(1,i).gt.10) goto 200

        kf=k(2,i)
        if(kf.eq.0) goto 200
        kfa=abs(kf)

c       pmra=0d0
c       if(mstu(42).ge.2) pmra=max(0d0,p(5,i))
        pmra=max(0d0,p(5,i))
        pt=max(1d-10,sqrt(p(1,i)**2+p(2,i)**2))
        pt2=p(1,i)**2+p(2,i)**2
        pmt=max(1d-20,pmra**2+pt2)
        yp=sign(
     $    log(min( (sqrt(pmt+p(3,i)**2)+abs(p(3,i)))/sqrt(pmt),
     &  1d20 )),p(3,i) )
        yp=yp+ylab

C......VNIBOOK: fill histograms for dET/deta transerve energy.
        plab = gamma*( p(3,i) + beta * p(4,i) )
        elab = gamma*( p(4,i) + beta * p(3,i) )
	ppl=sqrt(plab**2+pt2)
	et_lab=elab*sin(pt/max(ppl,1d-8))
	eta_lab=0.5d0*log(max(ppl+plab,1.d-8)/max(ppl-plab,1.d-8))
	etl=eta_lab+ylab
        call vfill1(8,eta_lab,et_lab*wevt/we)

        kfl1=mod(kfa/1000,10)
        kfl2=mod(kfa/100,10)
        kfl3=mod(kfa/10,10)

      if(kfa.le.6.or.kfa.eq.21) then
c...Diquarks.
      else if(kfl3.eq.0) then
      else if(kfa.le.100) then
c....Hadrons.
      else

        nn(0)=nn(0)+1
        kch=jamchge(kf)
        if(kch.ne.0) then
          nn(1)=nn(1)+1
        endif
        kc=jamcomp(kf)
        if(kc.gt.100) nn(kc*isign(1,kf))=nn(kc*isign(1,kf))+1

c....h- distributions.
        if(kf.eq.-211.or.kf.eq.-321.or.kf.eq.-2212) then
          nn(2)=nn(2)+1
          pmra=0.1357d0
          pmt=max(1d-10,pmra**2+p(1,i)**2+p(2,i)**2)
          yp1=sign(log(min( (sqrt(pmt+p(3,i)**2)+abs(p(3,i)))/sqrt(pmt),
     &     1d20 )),p(3,i) )
          yp1=yp1+ylab
          call vfill1(1,yp1,wevt/wy)
          weit=wevt/(wp*pt*0.2d0)
          yh0=2.9d0
          ii=0
          do iy=27,46
            ii=ii+1
            if(yp1.ge.yh0.and.yp1.lt.yh0+0.2d0) then
              call vfill1(iy,pt,weit)
              nyh(ii)=nyh(ii)+1
              pth(ii)=pth(ii)+pt
            endif
            yh0=yh0+0.2d0
          end do
        endif

c...Mesons.
        if(kfl1.eq.0) then

          ii=0
          if(kf.eq.-211)     ii=81
          if(kf.eq.111)      ii=82
          if(kf.eq.211)      ii=83
          if(kf.eq.-321)     ii=84
          if(abs(kf).eq.311) ii=85
          if(kf.eq.321)      ii=86
          if(ii.ne.0) call vfill1(ii,yp,wevt/wy)

c...Baryons.
        else

        wei=isign(1,kf)*wevt/wy
        weip=isign(1,kf)*wevt/(wp*pt*0.2d0)
        call vfill1(51,yp,wei)
        if(kf.gt.0) call vfill1(61,yp,wevt/wy)
        if(kf.lt.0) call vfill1(71,yp,wevt/wy)

c.....protons.
        if(kfa.eq.2212) then
          call vfill1(52,yp,wei)
          if(kf.gt.0) call vfill1(62,yp,wevt/wy)
          if(kf.lt.0) call vfill1(72,yp,wevt/wy)

c         if(abs(k(7,i)).ne.1) then
          yh0=2.2d0
          ii=0
          do iy=11,26
             ii=ii+1
            if(yp.ge.yh0.and.yp.lt.yh0+0.2d0) then
              call vfill1(iy,pt,weip)
              nyp(ii)=nyp(ii)+1
              ptp(ii)=ptp(ii)+pt
            endif
            yh0=yh0+0.2d0
          end do
c         endif

c....lambda
        else if(kfa.eq.3122) then
          call vfill1(53,yp,wei)
          if(kf.gt.0) then
            call vfill1(63,yp,wevt/wy)
            call vfill1(69,yp,wevt/wy)
          else if(kf.lt.0) then
            call vfill1(73,yp,wevt/wy)
            call vfill1(79,yp,wevt/wy)
          endif
c.....sigma-
        else if(kfa.eq.3112) then
          call vfill1(54,yp,wei)
          if(kf.gt.0) call vfill1(64,yp,wevt/wy)
          if(kf.lt.0) call vfill1(74,yp,wevt/wy)
c.....sigma0
        else if(kfa.eq.3212) then
          call vfill1(55,yp,wei)
          if(kf.gt.0) then
            call vfill1(65,yp,wevt/wy)
            call vfill1(69,yp,wevt/wy)
          else if(kf.lt.0) then
            call vfill1(75,yp,wevt/wy)
            call vfill1(79,yp,wevt/wy)
          endif
c.....sigma+
        else if(kfa.eq.3222) then
          call vfill1(56,yp,wei)
          if(kf.gt.0) call vfill1(66,yp,wevt/wy)
          if(kf.lt.0) call vfill1(76,yp,wevt/wy)
c.....xi-
        else if(kfa.eq.3312) then
          call vfill1(57,yp,wei)
          if(kf.gt.0) call vfill1(67,yp,wevt/wy)
          if(kf.lt.0) call vfill1(77,yp,wevt/wy)
c.....xi0
        else if(kfa.eq.3322) then
          call vfill1(58,yp,wei)
          if(kf.gt.0) call vfill1(68,yp,wevt/wy)
          if(kf.lt.0) call vfill1(78,yp,wevt/wy)
        endif

      endif
      endif

C...End analysis of particle record.
  200 continue


      return

C***********************************************************************

      entry ana3(mevt)

C...VNIBOOK: output of histograms.

      fac=1d0/mevt
      mnorm=0
      mform=1

c....baryon rapidity distributions.
      call vscale(1,fac)
      call vprint(1,mnorm,mform)
      call vscale(8,fac)
      call vprint(8,mnorm,mform)

c...net p transverse.
      do i=1,16
      call vscale(10+i,fac)
      call vprint(10+i,mnorm,mform)
      end do

c...h- transverse.
      do i=1,20
      call vscale(26+i,fac)
      call vprint(26+i,mnorm,mform)
      end do

c...rapidity
c     do i=1,8
c       call vscale(30+i,fac)
c       call vprint(30+i,mnorm,mform)
c     end do

c....baryon rapidity distributions.
      do i=1,9
      call vscale(50+i,fac)
      call vprint(50+i,mnorm,mform)
      call vscale(60+i,fac)
      call vprint(60+i,mnorm,mform)
      call vscale(70+i,fac)
      call vprint(70+i,mnorm,mform)
      end do

      do i=1,6
        call vscale(80+i,fac)
        call vprint(80+i,mnorm,mform)
      end do

      iunit=10
      wevt=fac
      write(iunit,*)'total particle',nn(0)*wevt
      write(iunit,*)'charged particle',nn(1)*wevt
      write(iunit,*)'h-(pi-,k-p~) particle',nn(2)*wevt
      do i=101,500
       if(nn(i).ge.1.or.nn(-i).ge.1) then
         kf=kchg(i,4)
         if(kchg(i,3).ne.0) then
           write(iunit,'(i5,1x,i9,1x,2(a16,1x,g12.3))')
     $      i,kf,chaf(i,1),nn(i)*wevt,chaf(i,2),nn(-i)*wevt
         else
           write(iunit,'(i5,1x,i9,1x,a16,1x,g12.3)')
     $      i,kf,chaf(i,1),nn(i)*wevt
         endif
       endif
      end do

c...mean pt.
      write(iunit,*)'<pt> of net protn'
      yp=2.2d0
      do i=1,16
        ptpx=0d0
        if(nyp(i).ge.1) ptpx=ptp(i)/nyp(i)
        write(iunit,'(2(f10.4,1x))')yp+0.1d0,ptpx
        yp=yp+0.2d0
      end do
      write(iunit,*)'<pt> of h-'
      yp=2.9d0
      do i=1,20
        pthx=0d0
        if(nyh(i).ge.1) pthx=pth(i)/nyh(i)
        write(iunit,'(2(f10.4,1x))')yp+0.1d0,pthx
        yp=yp+0.2d0
      end do



      end

