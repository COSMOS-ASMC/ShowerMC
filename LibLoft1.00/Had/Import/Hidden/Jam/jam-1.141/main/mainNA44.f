c...A main program for Pb(158GeV/c)+Pb NA44

      include 'jam2.inc'
      include 'jam1.inc'
      character frame*8,proj*8,targ*8,cwin*15
c     data bmin/0.0/,bmax/-1.08/,dt/100.0/,mevent/1/


c....Initialize JAM.
c     fname(1)='jam.cfg'  ! input file name.
      mstc(1) =48827   ! random seed.
      mevent=100        ! total simulation event
      bmin=0.0          ! minimum impact parameter
c     bmax=-3.6         ! maximum impact parameter 6.4%
      bmax=-3.0         ! maximum impact parameter 6.4%
      dt=100.0          ! collision time(fm/c)
      nstep=1
      cwin='158gev         '  ! incident energy
      frame='nn      '        ! comp. frame
      proj='208Pb   '         ! projectile
      targ='208Pb   '         ! target

      mstc(8)=1   ! job mode.
c     mstc(51)=0  ! BB collisions only
      mstc(74)=0   ! dipole-approximated QCD radiation of the string
c     mstc(81)=0   ! 1:hard scattering off
      mstc(156)=0  ! analysis of collision distribution
      mstc(155)=1  ! flow anal.
      mstc(163)=1  ! time evolution of directed transverse flow
      parc(7)= 1.0    ! Output time interval (fm/c)

      call jaminit(mevent,bmin,bmax,dt,nstep,
     $                             frame,proj,targ,cwin)
      nevent=mstc(2)

c...Initialize analysis.
      call anal1

c...Simulation start.
      do iev=1,nevent

c...Simulate one event.
        call jamevt(iev)

        if(mod(iev,100).eq.0) write(6,*)'event=',iev

c...Data analysis.
        call anal2

c...List phase space data.
c       call jamlist(1)

      end do

c...Final output.
      call jamfin

c...Print analysis results.
      call anal3

      end

c***********************************************************************

      subroutine anal1

      include 'jam2.inc'
      include 'jam1.inc'
      save wy,wp
      save ylab,wevt

      wevt=1./real(mstc(2))

c....Rapidity distribution.
      ylab=pard(17)
      ymin=-2.0
      ymax=7.0
      wy=0.25
      nymx=(ymax-ymin)/wy

      pmin=0.0
      pmax=5.0
      wp=0.1
      npmx=(pmax-pmin)/wp

C...Trasnverse mass distributions.
      call vbook1(31,'dN/dy proton',npmx,pmin,pmax)
      call vbook1(32,'dN/dy K+',npmx,pmin,pmax)
      call vbook1(33,'dN/dy pi+',npmx,pmin,pmax)
      call vbook1(34,'dN/dy antiproton',npmx,pmin,pmax)
      call vbook1(35,'dN/dy K-',npmx,pmin,pmax)
      call vbook1(36,'dN/dy pi-',npmx,pmin,pmax)

c...all rapidity
      call vbook1(41,'dN/dy proton',npmx,pmin,pmax)
      call vbook1(42,'dN/dy K+',npmx,pmin,pmax)
      call vbook1(43,'dN/dy pi+',npmx,pmin,pmax)
      call vbook1(44,'dN/dy antiproton',npmx,pmin,pmax)
      call vbook1(45,'dN/dy K-',npmx,pmin,pmax)
      call vbook1(46,'dN/dy pi-',npmx,pmin,pmax)


c....Transverse flow.
      call vbook1(51,'px(y) nucleons',nymx,ymin,ymax)
      call vbook1(52,'pt(y) nucleons',nymx,ymin,ymax)
      call vbook1(53,'px(y) pions   ',nymx,ymin,ymax)
      call vbook1(54,'pt(y) pions   ',nymx,ymin,ymax)
      call vbook1(61,'n(y) nucleons   ',nymx,ymin,ymax)
      call vbook1(62,'n(y) pions      ',nymx,ymin,ymax)

      call vbook1(71,'<px(y)> nucleons',nymx,ymin,ymax)
      call vbook1(72,'<pt(y)> nucleons',nymx,ymin,ymax)
      call vbook1(73,'<px(y)> pions   ',nymx,ymin,ymax)
      call vbook1(74,'<pt(y)> pions   ',nymx,ymin,ymax)

c....Freaze-out point
      call vbook1(81,'dN/dr proton',40,0.,30.)
      call vbook1(82,'dN/dr pion+',40,0.,30.)
      call vbook1(83,'dN/dr pion-',40,0.,30.)
      call vbook1(84,'dN/dr kaon+',40,0.,30.)

c....Freaze-out time
      call vbook1(91,'dN/dt proton',40,0.,30.)
      call vbook1(92,'dN/dt pion+',40,0.,30.)
      call vbook1(93,'dN/dt pion-',40,0.,30.)
      call vbook1(94,'dN/dt kaon+',40,0.,30.)

      return

c***********************************************************************

      entry anal2

c...Loop over all particles.
      do i=1,nv

        kf=k(2,i)
        kc=lucomp(kf)
        if(kc.le.0.or.kc.gt.mstu(6)) then
           write(6,*)'Invalid code i kf kc',i,kf,kc,nv,nbary,nmeson
           goto 3000
        endif

        y=0.5*log( max(p(4,i)+p(3,i),1.e-8)/max(p(4,i)-p(3,i),1.e-8) )

        if(mstc(4).eq.0) then
        else if(mstc(4).eq.3) then
        else
         y=y+ylab
        endif

        ptsq=p(1,i)**2+p(2,i)**2
        pt=sqrt(ptsq)
        pp=sqrt(ptsq+p(3,i)**2)
        pt=max(pt,1.e-8)
        eta=0.5*log( max(pp+p(3,i),1.e-8)/max(pp-p(3,i),1.e-8) )
	et=p(4,i)*pt/max(pp,1.e-8)
        emt=sqrt(p(5,i)**2+ptsq)
	emt0=emt-p(5,i)
	px=p(1,i)

c......Transverse flow of nucleons and pions.
	if(kf.eq.2212.or.kf.eq.2112) then
	  call vfill1(51,y,px)
	  call vfill1(52,y,pt)
	  call vfill1(61,y,1.0)
        else if(kf.eq.111.or.abs(kf).eq.211) then
	  call vfill1(53,y,px)
	  call vfill1(54,y,pt)
	  call vfill1(62,y,1.0)
	endif

c......Protons
        if(kf.eq.2212) then
	  if(y.ge.2.4.and.y.le.2.8) call vfill1(31,emt0,wevt/wp/emt)
	  call vfill1(41,emt0,wevt/wp/emt)
c......Kaons
        else if(kf.eq.321) then
	  if(y.ge.2.5.and.y.le.3.4) call vfill1(32,emt0,wevt/wp/emt)
	  call vfill1(42,emt0,wevt/wp/emt)
c.....Positive pions
        else if(kf.eq.211) then
	  if(y.ge.3.0.and.y.le.4.1) call vfill1(33,emt0,wevt/wp/emt)
	  call vfill1(43,emt0,wevt/wp/emt)

        else if(kf.eq.-2212) then
	  if(y.ge.2.4.and.y.le.2.8) call vfill1(34,emt0,wevt/wp/emt)
	  call vfill1(44,emt0,wevt/wp/emt)
        else if(kf.eq.-321) then
	  if(y.ge.2.5.and.y.le.3.4) call vfill1(35,emt0,wevt/wp/emt)
	  call vfill1(45,emt0,wevt/wp/emt)
        else if(kf.eq.-211) then
	  if(y.ge.3.0.and.y.le.4.1) call vfill1(36,emt0,wevt/wp/emt)
	  call vfill1(46,emt0,wevt/wp/emt)
	endif

3000  end do

      return

c***********************************************************************

      entry anal3

c...Output of histograms.

      fac=1.0
      mnorm=0
      mform=1
      do i=1,6
	call vscale(30+i,fac)
	call vprint(30+i,mnorm,mform)
	call vscale(40+i,fac)
	call vprint(40+i,mnorm,mform)
      end do

c....Transverse flow.
      call vopera(51,'/',61,71,1.0,1.0)
      call vopera(52,'/',61,72,1.0,1.0)
      call vopera(53,'/',62,73,1.0,1.0)
      call vopera(54,'/',62,74,1.0,1.0)
      call vprint(61,mnorm,mform)
      call vprint(62,mnorm,mform)
      do i=1,4
	call vprint(70+i,mnorm,mform)
      end do

c...Freaze-out points
      do i=1,4
       call vprint(80+i,1,0)
       call vprint(90+i,1,0)
      end do

      end

