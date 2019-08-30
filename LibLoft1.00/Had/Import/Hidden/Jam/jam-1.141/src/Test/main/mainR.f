c...A main program for Au+Au RHIC

      include 'jam1.inc'
      include 'jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15
      real*8 jamemjet
      logical dump 
      data dump/.false./
c     data dump/.true./


      if(dump)
     $  open(33,file='phase.dat',form='unformatted',status='new')


c....Initialize JAM.
c     fname(1)='jam.cfg'  ! input file name.
      mstc(1) =12971  ! random seed.
      mevent=1        ! total simulation event
      bmin=0.0D0          ! minimum impact parameter
      bmax=3.0D0         ! maximum impact parameter
c     bmax=-100d0        ! maximum impact parameter
      dt=100.0D0          ! collision time(fm/c)
      nstep=1 
      cwin='200gev         '  ! incident energy
      frame='collider'        ! comp. frame
c     frame='nn'        ! comp. frame

      proj='197Au'         ! projectile
      targ='197Au'         ! target

c     mstc(8)=0          ! job mode.
c     mstc(156)=1        ! analysis of collision distribution
c     mstc(6)=-1         ! Glauber
c     mstc(51)=0         ! 0:only BB collision.

c     mstc(74)=1         ! soft rad.
c     mstc(76)=0         ! treatment of stirng deay
c     mstc(81)=0         ! 1:hard scattering on/off
c     mstc(84)=1         !  Option for nuclear shadowing effect.

c     mstc(155)=0        ! flow anal.
c     mstc(156)=1        ! analysis of collision distribution
c     mstc(156)=1        ! energy distribution of collisions
c     mstc(162)=1        ! Output collision histroy
c     mstc(165)=1        ! Output time evol. of particle yield
c     mstc(166)=1        ! Output time evol. of particle density
c     mstc(167)=1        ! Output time evol. of density (Gaussian)
c     parc(7)= 0.2D0     ! Output time interval (fm/c)

c     parc(54)=2.0d0     ! string tension

      call jaminit(mevent,bmin,bmax,dt,nstep,frame,proj,targ,cwin)
      nevent=mstc(2)
      if(dump)write(33)nevent,pard(17),pard(5),pard(6),mstc(4)


c...Initialize analysis.
      call anal1

c...Simulation start.
      do iev=1,nevent

c...Simulate one event.
        call jamevt(iev)

c...Dump phase space data.
        if(dump) then
          write(33)iev,nv,nbary,nmeson,pard(2)
          do i=1,nv
c           write(33)(k(j,i),j=1,11),(r(j,i),j=1,5),(p(j,i),j=1,5)
c    $              ,(v(j,i),j=1,5)
            write(33)(k(j,i),j=1,9),(r(j,i),j=1,5),(p(j,i),j=1,5)
          end do
        endif

        if(mod(iev,100).eq.0) write(6,*)'event=',iev

c...Data analysis.
        call anal2

c...List phase space data.
c       call jamlist(1)

      end do

      if(dump) close(33)

c...Final output.
      call jamfin

c...Print analysis results.
      call anal3

      write(6,*)'max res',mstd(199)

      end

c***********************************************************************

      subroutine anal1

      include 'jam1.inc'
      include 'jam2.inc'
      dimension npa(0:20)
      save wy,wp
      save ylab,yproj,ytarg
      save npa

c...npa(1) : average charged
c...npa(2) : average negative
c...npa(3) : pion-
c...npa(4) : pion0
c...npa(5) : pion+
c...npa(6) : antikaon
c...npa(7) : Kaon
c...npa(8) : lambda
c...npa(9) : a-lambda
c...npa(10): sigma
c...npa(11): a-sigma
c...npa(12): proton
c...npa(13): a-proton
c...npa(14):

      do i=0,20
       npa(i)=0
      end do


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

      if(pard(16).le.40.D0) then
        pmin=0.0D0
        pmax=3.0D0
        npmx=30
      else
        pmin=0.0D0
        pmax=10.0D0
        npmx=30
      endif
      wp=(pmax-pmin)/npmx
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

      entry anal2

      nch=0
      neg=0
c...Loop over all particles.
      do i=1,nv

        if(k(1,i).gt.10) goto 3000
        kf=k(2,i)
        kfa=abs(kf)
        kc=jamcomp(kf)
        if(kc.le.0.or.kc.gt.mstu(6)) then
           write(6,*)'Invalid code i kf kc',i,kf,kc,nv,nbary,nmeson
           goto 3000
        endif
        if(abs(k(7,i)).ne.1) then
          if(v(1,i).gt.1D+15.or.v(2,i).gt.1D+15.or.v(3,i).gt.1D+15)then
            write(6,*)'vv inf kf',kf,(v(j,i),j=1,5)
            goto 3000
          endif
          vv=v(1,i)**2+v(2,i)**2+v(3,i)**2
          if(vv.le.0.0D0) then
            write(6,*)'vv?? kf',kf,(v(j,i),j=1,5)
            goto 3000
          endif
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
c       et=sqrt(p(5,i)**2+ptsq)

        npa(0)=npa(0)+1
c....Et didst.
        call vfill1(4,rap,et/wy)
        call vfill1(102,rap,et/wy)

c...Charged particles.
        kch=jamchge(kf)
        if(kch.ne.0) then
          nch=nch+1
          npa(1)=npa(1)+1
          call vfill1(101,rap,1.D0/wy)
          call vfill1(1,rap,1.D0/wy)
          call vfill1(5,eta,1.D0/wy)
          call vfill1(51,pt,1.D0/wp/pt/2) 
          if(abs(rap).le.0.5d0) call vfill1(71,pt,1.D0/wp/pt/2) 

c...Negative charged particles.
          if(kch.lt.0) then
            neg=neg+1
            npa(2)=npa(2)+1
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
              vry=v(1,i)**2+v(2,i)**2+v(3,i)**2
              if(vry.ge.0.0D0) then
                vr=sqrt(vry)
              else
                write(6,*)'hyperon v<0',kf,(v(j,i),j=1,5)
                goto 3000
              endif
              call vfill1(206,vr,1.D0) 
              call vfill1(216,v(4,i),1.D0) 
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
              vrl=v(1,i)**2+v(2,i)**2+v(3,i)**2
              if(vrl.ge.0.0D0) then
                vr=sqrt(vrl)
              else
                write(6,*)'v<0',kf,(v(j,i),j=1,5)
                goto 3000
              endif
              call vfill1(200+iii,vr,1.D0) 
              call vfill1(210+iii,v(4,i),1.D0) 
            endif
      endif
 

3000  end do
      call vfill1(6,dble(nch),1.D0) 

      return

c***********************************************************************

      entry anal3

c...Output of histograms.

c...Event weight
      fac=1.D0/dble(mstc(2))

c...ET distributions.
      do i=1,6
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
        call vprint(100+i,0,1)
      end do

      open(70,file='file70',status='unknown')
      write(70,*)'ylab yproj ytarg=',ylab,yproj,ytarg
      write(70,*)'average number of jet',pard(87)*fac
      close(70)

      end

