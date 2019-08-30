c....Analysis program produced from mainR.f

      include 'jam1.inc'
      include 'jam2.inc'
      common/myana/ispec

      open(33,file='phase.dat',form='unformatted',status='old')


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
c       write(6,*) iev,nv,nbary,nmeson,b
	pard(2)=b
        nbary=0
        do i=1,nv
          read(33)(k(j,i),j=1,9),(r(j,i),j=1,5),(p(j,i),j=1,5)

c         write(6,*)(k(j,i),j=1,9),(r(j,i),j=1,5),(p(j,i),j=1,5)

c          if(ie.eq.1) then
c            write(51,*)(r(j,i),j=1,5) 
c            write(51,*)(v(j,i),j=1,5)
c          endif
c          if(v(4,i).lt.abs(v(3,i))) then
c            write(51,*)i,k(2,i),k(7,i),(v(j,i),j=1,5) 
c          endif

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

        mevt=mevt+1
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

c***********************************************************************

      subroutine ana1

      include 'jam1.inc'
      include 'jam2.inc'
      save wy,wp
      save ylab,yproj,ytarg

      ylab=pard(17)

c....Rapidity distribution.
      ylab=pard(17)
      dbetpr=pard(35)
      dbetta=pard(45)
      yproj=0.5d0*log((1.0d0+dbetpr)/(1.0d0-dbetpr))
      ytarg=0.5d0*log((1.0d0+dbetta)/(1.0d0-dbetta))

      ymin=ytarg*2.0D0
      ymax=yproj*2.0D0
      ymin=-10.D0
      ymax=10.D0
      nymx=30 ! RHIC
      if(mstc(4).eq.0) then
      else if(mstc(4).eq.3) then
      else
       ymax=ymax+ylab
       ymin=ymin+ylab
      endif
      wy=(ymax-ymin)/nymx
c     print *,'wy=',wy,ymin,ymax

      ymin=-10.0D0
      ymax=10.0D0
      wy=0.75D0
      nymx=(ymax-ymin)/wy

      call vbook1(1,'dN/dy - charged',nymx,ymin,ymax)
      call vbook1(2,'dN/dy - negative',nymx,ymin,ymax)
      call vbook1(3,'dN/dy - pi- k- pbar',nymx,ymin,ymax)
c...Transverse energy distributions.
      call vbook1(4,'dET/dy - all',nymx,ymin,ymax)
      call vbook1(5,'dN/d(eta) - charged',nymx,ymin,ymax)
      call vbook1(6,'P(n) - charged',50,0.D0,100.D0)
      call vbook1(7,'dET/d(eta) - all',nymx,ymin,ymax)
      call vbook1(8,'dET1/d(eta) - all',nymx,ymin,ymax)

      call vbook1(11,'dN/dy - net baryon',nymx,ymin,ymax)
      call vbook1(12,'dN/dy - p-p-bar',nymx,ymin,ymax)
      call vbook1(13,'dN/dy - net lambda',nymx,ymin,ymax)
      call vbook1(14,'dN/dy - net sigma',nymx,ymin,ymax)
      call vbook1(15,'dN/dy - net xi',nymx,ymin,ymax)
      call vbook1(16,'dN/dy - net omega',nymx,ymin,ymax)

      call vbook1(21,'dN/dy - baryon',nymx,ymin,ymax)
      call vbook1(22,'dN/dy - proton',nymx,ymin,ymax)
      call vbook1(23,'dN/dy - lambda',nymx,ymin,ymax)
      call vbook1(24,'dN/dy - sigma',nymx,ymin,ymax)
      call vbook1(25,'dN/dy - xi',nymx,ymin,ymax)
      call vbook1(26,'dN/dy - omega',nymx,ymin,ymax)

      call vbook1(31,'dN/dy - a-baryon',nymx,ymin,ymax)
      call vbook1(32,'dN/dy - a-proton',nymx,ymin,ymax)
      call vbook1(33,'dN/dy - a-lambda',nymx,ymin,ymax)
      call vbook1(34,'dN/dy - a-sigma',nymx,ymin,ymax)
      call vbook1(35,'dN/dy - a-xi',nymx,ymin,ymax)
      call vbook1(36,'dN/dy - a-omega',nymx,ymin,ymax)

      call vbook1(41,'dN/dy - pion- ',nymx,ymin,ymax)
      call vbook1(42,'dN/dy - pion0 ',nymx,ymin,ymax)
      call vbook1(43,'dN/dy - pion+ ',nymx,ymin,ymax)
      call vbook1(44,'dN/dy - k- ',nymx,ymin,ymax)
      call vbook1(45,'dN/dy - k0 ',nymx,ymin,ymax)
      call vbook1(46,'dN/dy - ak0 ',nymx,ymin,ymax)
      call vbook1(47,'dN/dy - k+ ',nymx,ymin,ymax)
      call vbook1(48,'dN/dy - D ',nymx,ymin,ymax)
      call vbook1(49,'dN/dy - aD ',nymx,ymin,ymax)

c...Transverse distributions.

c     if(pard(16).le.40.D0) then
c       pmin=0.0D0
c       pmax=3.0D0
c       npmx=30
c     else
        pmin=0.0D0
        pmax=10.0D0
        npmx=30
c     endif
      wp=(pmax-pmin)/npmx
      print *,'pmin pmax wp',pmin,pmax,wp
      call vbook1(51,'dN/dpt**2 - charged',npmx,pmin,pmax)
      call vbook1(52,'dN/dpt**2 - proton ',npmx,pmin,pmax)
      call vbook1(53,'dN/dpt**2 - lambda ',npmx,pmin,pmax)
      call vbook1(54,'dN/dpt**2 - sigma ',npmx,pmin,pmax)
      call vbook1(55,'dN/dpt**2 - xi ',npmx,pmin,pmax)
      call vbook1(56,'dN/dpt**2 - omega ',npmx,pmin,pmax)

      call vbook1(61,'dN/dpt**2 - pion',npmx,pmin,pmax)
      call vbook1(62,'dN/dpt**2 - kaon',npmx,pmin,pmax)
      call vbook1(63,'dN/dpt**2 - Dmeson',npmx,pmin,pmax)

      call vbook1(71,'dN/dpt**2 y<0.5- charged',npmx,pmin,pmax)
      call vbook1(72,'dN/dpt**2 y<0.5- proton ',npmx,pmin,pmax)
      call vbook1(73,'dN/dpt**2 y<0.5- lambda ',npmx,pmin,pmax)
      call vbook1(74,'dN/dpt**2 y<0.5- sigma ',npmx,pmin,pmax)
      call vbook1(75,'dN/dpt**2 y<0.5- xi ',npmx,pmin,pmax)
      call vbook1(76,'dN/dpt**2 y<0.5- omega ',npmx,pmin,pmax)

      call vbook1(81,'dN/dpt**2 y<0.5- pion',npmx,pmin,pmax)
      call vbook1(82,'dN/dpt**2 y<0.5- kaon',npmx,pmin,pmax)
      call vbook1(83,'dN/dpt**2 y<0.5- Dmeson',npmx,pmin,pmax)

c....Freaze-out point
      rmax=50.0D0
      wr=0.5D0
      nrmax=rmax/wr
      call vbook1(201,'dN/dr proton',nrmax,0.D0,rmax)
      call vbook1(202,'dN/dr lambda',nrmax,0.D0,rmax)
      call vbook1(203,'dN/dr sigma',nrmax,0.D0,rmax)
      call vbook1(204,'dN/dr pion+ ',nrmax,0.D0,rmax)
      call vbook1(205,'dN/dr pion- ',nrmax,0.D0,rmax)
      call vbook1(206,'dN/dr kaon+ ',nrmax,0.D0,rmax)
      call vbook1(207,'dN/dr kaon- ',nrmax,0.D0,rmax)
      call vbook1(208,'dN/dr D     ',nrmax,0.D0,rmax)

c....Freaze-out time
      tmax=90.0D0
      wt=2.0D0
      ntmax=tmax/wt
      call vbook1(211,'dN/dt proton',ntmax,0.D0,tmax)
      call vbook1(212,'dN/dt lambda',ntmax,0.D0,tmax)
      call vbook1(213,'dN/dt sigma ',ntmax,0.D0,tmax)
      call vbook1(214,'dN/dt pion+ ',ntmax,0.D0,tmax)
      call vbook1(215,'dN/dt pion- ',ntmax,0.D0,tmax)
      call vbook1(216,'dN/dt kaon+ ',ntmax,0.D0,tmax)
      call vbook1(217,'dN/dt kaon- ',ntmax,0.D0,tmax)
      call vbook1(218,'dN/dt D     ',ntmax,0.D0,tmax)



      call vbook1(101,'dN/dy - charged',nymx,ymin,ymax)
      call vbook1(102,'dE_t/dy',nymx,ymin,ymax)
      call vbook1(103,'dN/dy - net baryon',nymx,ymin,ymax)
      call vbook1(104,'dN/2piptdy - pion -0.5<y<0.5',npmx,pmin,pmax)

      return

c***********************************************************************

      entry ana2

      beta=pard(5)
      gamma=pard(6)

      ispec=0
      nch=0
      neg=0
c...Loop over all particles.
      do 3000 i=1,nv

       if(k(1,i).ge.10) goto 3000

c...Exclude spectetor.
c       if(abs(k(7,i)).eq.1) then
c         y=0.5*log(max(p(4,i)+p(3,i),1.e-8)/max(p(4,i)-p(3,i),1.e-8))
c         write(3,*)k(1,i),k(2,i),k(7,i),y
c         ispec=ispec+1
c         goto 3000
c       endif

        kf=k(2,i)
        kfa=abs(kf)
        kc=jamcomp(kf)
        if(kc.le.0.or.kc.gt.mstu(6)) then
           write(6,*)'Invalid code i kf kc',i,kf,kc,nv,nbary,nmeson
           goto 3000
        endif

        rap=0.5D0*log( max(p(4,i)+p(3,i),1.D-8)/max(p(4,i)-p(3,i), 
     & 1.D-8) )

        if(mstc(4).eq.0) then
        else if(mstc(4).eq.3) then
        else
         rap=rap+ylab
        endif

        ptsq=max(1.D-8,p(1,i)**2+p(2,i)**2)
        pt=sqrt(ptsq)
        pp=sqrt(ptsq+p(3,i)**2)
        pt=max(pt,1.D-8)
        eta=0.5D0*log( max(pp+p(3,i),1.D-8)/max(pp-p(3,i),1.D-8) )
        et=p(4,i)*sin(pt/max(pp,1.D-8))
        et1=p(4,i)*sin(pt/max(pp,1.D-8))+p(5,i)
c       et=sqrt(p(5,i)**2+ptsq)

c....Et didst.
        call vfill1(4,rap,et/wy)
        call vfill1(7,eta,et/wy)
        call vfill1(8,eta,et1/wy)
        call vfill1(102,rap,et/wy)

c...Charged particles.
        kch=jamchge(kf)
        if(kch.ne.0) then
          nch=nch+1
          call vfill1(101,rap,1.D0/wy)
          call vfill1(1,rap,1.D0/wy)
          call vfill1(5,eta,1.D0/wy)
          call vfill1(51,pt,1.D0/wp/pt/2) 
          if(abs(rap).le.0.5d0) call vfill1(71,pt,1.D0/wp/pt/2) 

c...Negative charged particles.
          if(kch.lt.0) then
            neg=neg+1
            call vfill1(2,rap,1.D0/wy)
c...h-(pi-,k-,p~)
            if(kf.eq.-211.or.kf.eq.-321.or.kf.eq.-2212)
     $       call vfill1(3,rap,1.D0/wy)
          endif
        endif


c...Transverse momentum at midrapidiy.
      if(abs(rap).le.0.5d0) then
         if(kf.eq.111.or.kfa.eq.211) then
           call vfill1(104,pt,1d0/(paru(2)*pt*wp))
         endif
      endif

c...Net baryons.
        ibar=kchg(kc,6)
        if(ibar.eq.3) then
          call vfill1(11,rap,isign(1,kf)/wy)
          call vfill1(103,rap,isign(1,kf)/wy)
        endif

c....Baryon dist.
      iii=0
      if(kfa.eq.2212) iii=2 ! proton
      if(kfa.eq.3122) iii=3 ! lambda
      if(kfa.eq.3212.or.kfa.eq.3222.or.kfa.eq.3112) iii=4 ! sigma
      if(kfa.eq.3312.or.kfa.eq.3322) iii=5 ! xi
      if(kfa.eq.3334) iii=6 ! Omega
      if(iii.ne.0) then
         call vfill1(10+iii,rap,isign(1,kf)/wy)
         if(kf.gt.0) then
           call vfill1(20+iii,rap,1.0d0/wy)
           call vfill1(50+iii,pt,1.D0/wp/pt/2) 
           if(abs(rap).le.0.5d0) call vfill1(70+iii,pt,1.D0/wp/pt/2) 
         else if(kf.lt.0) then
           call vfill1(30+iii,rap,1.0d0/wy)
         endif
      endif

c...meson rap. dist.
      iii=0
      if(kf.eq.-211) iii=1 ! pi
      if(kf.eq.111)  iii=2 ! pi0
      if(kf.eq.211)  iii=3 ! pi+
      if(kf.eq.-321) iii=4 ! k-
      if(kf.eq.-311) iii=5 ! ak0
      if(kf.eq. 321) iii=6 ! k+
      if(kf.eq. 311) iii=7 ! k0
      if(kf.eq.411.or.kf.eq.421) iii=8 ! D
      if(kf.eq.-411.or.kf.eq.-421) iii=8 ! aD
      if(iii.ne.0) call vfill1(40+iii,rap,1.0d0/wy)

c...meson pt dist.
      iii=0
      if(kfa.eq.211.or.kfa.eq.111) iii=1 ! pion
      if(kfa.eq.311.or.kfa.eq.321) iii=2 ! kaon
      if(kfa.eq.411.or.kfa.eq.421) iii=3 ! D
      if(iii.ne.0) then
        call vfill1(60+iii,pt,1.D0/wp/pt/2) 
        if(abs(rap).le.0.5d0) call vfill1(80+iii,pt,1.D0/wp/pt/2) 
      endif


c....Hyperons.
        if(kf.eq.3122.or.kf.eq.3212.or.kf.eq.3222.or.kf.eq.3112) then
            if(abs(k(7,i)).ne.1) then
              vry=r(1,i)**2+r(2,i)**2+r(3,i)**2
              if(vry.ge.0.0D0) then
                vr=sqrt(vry)
              else
                write(6,*)'hyperon v<0',kf,(r(j,i),j=1,5)
                goto 3000
              endif
              call vfill1(206,vr,1.D0) 
              call vfill1(216,r(4,i),1.D0) 
            endif
        endif


c...Freeze out.
      iii=0
      if(kf.eq.2212) iii=1
      if(kf.eq.3122) iii=2
      if(kfa.eq.3212.or.kfa.eq.3222.or.kfa.eq.3112) iii=3 ! sigma
      if(kf.eq. 211) iii=4
      if(kf.eq.-211) iii=5
      if(kf.eq. 321) iii=6
      if(kf.eq.-321) iii=7
      if(kfa.eq.411.or.kfa.eq.421) iii=8 ! D
      if(iii.ne.0) then
            if(abs(k(7,i)).ne.1) then
              vrl=r(1,i)**2+r(2,i)**2+r(3,i)**2
              if(vrl.ge.0.0D0) then
                vr=sqrt(vrl)
              else
                write(6,*)'v<0',kf,(r(j,i),j=1,5)
                goto 3000
              endif
              call vfill1(200+iii,vr,1.D0) 
              call vfill1(210+iii,r(4,i),1.D0) 
            endif
      endif


c------------------------------------------------------------
3000  continue

      return

c***********************************************************************

      entry ana3(mevt)

c...Output histograms.

c...Event weight
      print *,'event',mevt
      wevt=1.d0/dble(mevt)
c     fac=1.D0/dble(mstc(2))
      fac=wevt
      mnorm=0
      mform=1

c...ET distributions.
      do i=1,8
       call vscale(i,fac)
       call vprint(i,0,1)
      end do

c...Baryon Rapidity distributions.
      do i=1,6
       call vscale(10+i,fac)
       call vprint(10+i,0,0)
       call vscale(20+i,fac)
       call vprint(20+i,0,0)
       call vscale(30+i,fac)
       call vprint(30+i,0,0)
      end do

c...meson rap. dist.
      do i=1,9
       call vscale(40+i,fac)
       call vprint(40+i,0,0)
      end do

c...baryon Pt distributions.
      do i=1,6
       call vscale(50+i,fac)
       call vprint(50+i,0,1)
       call vscale(70+i,fac)
       call vprint(70+i,0,1)
      end do

c...meson Pt distributions.
      do i=1,3
       call vscale(60+i,fac)
       call vprint(60+i,0,1)
       call vscale(80+i,fac)
       call vprint(80+i,0,1)
      end do


c...Freaze-out points
c     do i=1,8
c      call vprint(200+i,1,0)
c      call vprint(210+i,1,0)
c     end do

      do i=1,4
        call vprint(100+i,0,fac)
      end do

      open(70,file='file70',status='unknown')
      write(70,*)'ylab yproj ytarg=',ylab,yproj,ytarg
      write(70,*)'average number of jet',pard(87)*fac
      close(70)


      end

