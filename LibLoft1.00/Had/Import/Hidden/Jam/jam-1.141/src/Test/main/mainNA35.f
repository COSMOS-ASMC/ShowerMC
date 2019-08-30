C***********************************************************************
 
      program mainNA35
 
c...NA35
c...CERN/SPS for S+S P.R.L 70(1994)1419 NA35 Veto trigger. 2%sigma_trig

C...Purpose: Generate fixed target S+S collisions at CERN SPS beam 
C...momentum P_beam = 200 A GeV (corresponding to root(s) = 17.8 A GeV).
C...Impact parametewr is chosen equal to zero, i.e. perfectly head-on.
C...Get histograms for dN/deta, dN/dpT^2, dET/deta, as well as for 
C...particle and energy fluctuations.

      include 'jam1.inc'
      include 'jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15

c....Initialize JAM.
      mstc(1)=1712915     ! random seed.
      mevent=300          ! total simulation event
      bmin=0.0d0          ! minimum impact parameter in fm
      bmax=-1.08d0        ! maximum impact parameter in fm
c....Eur.Phys.J 2.9%
c...WA80 collaboration. 25% 1450mb.
c     bmax=sqrt(0.25d0*1450*0.1/paru(1))
      dt=100.0d0          ! collision time(fm/c)
      nstep=1
      cwin='200gevc'      ! incident energy
      frame='nn'          ! comp. frame
      proj='32S'          ! projectile
      targ='32S'          ! target

c...Options.
      mstc(8)=0           ! job mode.
c     mstc(81)=1          ! 1:hard scattering on/off
c     mstc(156)=1         ! analysis of collision distribution
c     mstc(6)=-1          ! Glauber
c     mstc(51)= 1         ! only BB/MB collision.
c     mstc(51)= 0         ! 0:only BB collision.
c     mstc(74)= 1         ! soft rad.
c     parc(54)=2d0        ! string tension

      call jaminit(mevent,bmin,bmax,dt,nstep,frame,proj,targ,cwin)

c...Initialize analysis.
      call anal1

c...Simulation start.
      do iev=1,mevent

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
      
 
C***********************************************************************

      subroutine anal1

C...Initial values and definitions for histograms.
      include 'jam1.inc'
      include 'jam2.inc'
      save wevt,wy,wp,ylab

      wevt=1.0d0/mstc(2)
      ylab=pard(17)
      print *,'ylab',ylab

C.....Limits and bins for rapidity distributions.
      ymin=-2.d0
      ymax=+7.d0
      wy=0.25d0
      nymx=(ymax-ymin)/wy
C.....Limits and bins for pT-distributions.
      pmin=0.0d0
      pmax=5.0d0
      wp=0.1d0
      npmx=(pmax-pmin)/wp

C....Trasnverse energy distributions dET/dy.
      call vbook1(8, 'dET/dy',nymx,ymin,ymax)

c...Rapidity spectra.
      call vbook1(31,'dN/dy proton',nymx,ymin,ymax)
      call vbook1(32,'dN/dy net proton',nymx,ymin,ymax)
      call vbook1(33,'dN/dy h-(pi-k-pbar)',nymx,ymin,ymax)
      call vbook1(34,'dN/dy lambda',nymx,ymin,ymax)
      call vbook1(35,'dN/dy sigma-',nymx,ymin,ymax)
      call vbook1(36,'dN/dy sigma0',nymx,ymin,ymax)
      call vbook1(37,'dN/dy sigma+',nymx,ymin,ymax)
      call vbook1(38,'dN/dy lambda-like',nymx,ymin,ymax)
      call vbook1(39,'dN/dy K+',nymx,ymin,ymax)
      call vbook1(40,'dN/dy K-',nymx,ymin,ymax)

c...CERN/SPS for S+S P.R.L 70(1994)1419 NA35 Veto trigger. 2%sigma_trig
c...and EPJ.
      call vbook1(41,'1/ptdN/dpt 0.5<y<3.0 protons',npmx,pmin,pmax)
      call vbook1(42,'1/ptdN/dpt 0.8<y<2.0 h-',npmx,pmin,pmax)
      call vbook1(43,'1/ptdN/dpt 2.0<y<3.0 h-',npmx,pmin,pmax)
      call vbook1(44,'1/ptdN/dpt 3.0<y<4.0 h-',npmx,pmin,pmax)
      call vbook1(45,'1/ptdN/dpt 4.0<y<4.4 h-',npmx,pmin,pmax)

c...NA44 PRL 78(1997)2080  10%
c     call vbook1(101,'(s)1/mtdN/dmt 2.4<y<2.8 p',ntmx,tmin,tmax)
c     call vbook1(102,'(s)1/mtdN/dmt 2.5<y<3.4 k+',ntmx,tmin,tmax)
c     call vbook1(103,'(s)1/mtdN/dmt 3.0<y<4.1 pi+',ntmx,tmin,tmax)

c...WA80  EPJ  25%, nucl-ex/9805007
      call vbook1(121,'   1/ptdN/dpt 2.1<y<2.9 pi0',npmx,pmin,pmax)


c....Freaze-out point
      rmax=50.0D0
      wr=0.5D0
      nrmax=rmax/wr
      call vbook1(81,'dN/dr nucleon',nrmax,0.D0,rmax)
      call vbook1(82,'dN/dr Y      ',nrmax,0.D0,rmax)
      call vbook1(83,'dN/dr pion   ',nrmax,0.D0,rmax)
      call vbook1(84,'dN/dr kaon0  ',nrmax,0.D0,rmax)
      call vbook1(85,'dN/dr K+     ',nrmax,0.D0,rmax)
      call vbook1(86,'dN/dr K-     ',nrmax,0.D0,rmax)

c....Freaze-out time
      tmax=90.0D0
      wt=2.0D0
      ntmax=tmax/wt
      call vbook1(91,'dN/dt nucleon',ntmax,0.D0,tmax)
      call vbook1(92,'dN/dt Y      ',ntmax,0.D0,tmax)
      call vbook1(93,'dN/dt pion   ',ntmax,0.D0,tmax)
      call vbook1(94,'dN/dt kaon0  ',ntmax,0.D0,tmax)
      call vbook1(95,'dN/dt K+     ',ntmax,0.D0,tmax)
      call vbook1(96,'dN/dt K-     ',ntmax,0.D0,tmax)


      return

C***********************************************************************

      entry anal2

C...Do analysis of particle record.
      do 200 i=1,nv

C.......Skip all "dead" particles, photons, leptons, neutrinos, etc.
        if(k(1,i).le.0.or.k(1,i).gt.10) then
         print *,'k1 kf',k(1,i),k(2,i)
        endif
        if(k(1,i).le.0.or.k(1,i).gt.10) goto 200
        if(iabs(k(2,i)).le.100) goto 200

C.......Find charge.
        kch=jamchge(k(2,i))

        kf=k(2,i)
        kfa=abs(kf)
        y=0.5d0*log(max(p(4,i)+p(3,i),1d-8)/max(p(4,i)-p(3,i),1d-8))
        ycm=y
        if(mstc(4).eq.0) then
        else if(mstc(4).eq.3) then
        else
         y=y+ylab
        endif

        ptsq=p(1,i)**2+p(2,i)**2
        pt=sqrt(ptsq)
        pp=sqrt(ptsq+p(3,i)**2)
        pt=max(pt,1.d-8)
c       eta=0.5*log( max(pp+p(3,i),1.e-8)/max(pp-p(3,i),1.e-8) )
        et=p(4,i)*pt/max(pp,1.d-8)

c...Transverse energy
        call vfill1(8,y,et*wevt/wy)

c......h-(CERN/SPS definition)
	if(kf.eq.-211.or.kf.eq.-321.or.kf.eq.-2212) then
	  call vfill1(33,y,wevt/wy)
	  if(y.ge.0.8d0.and.y.le.2.0d0) call vfill1(42,pt,wevt/(wp*pt))
	  if(y.ge.2.0d0.and.y.le.3.0d0) call vfill1(43,pt,wevt/(wp*pt))
	  if(y.ge.3.0d0.and.y.le.4.0d0) call vfill1(44,pt,wevt/(wp*pt))
	  if(y.ge.4.0d0.and.y.le.4.4d0) call vfill1(45,pt,wevt/(wp*pt))
	endif

c...WA80  EPJ  25%, nucl-ex/9805007
          if(kf.eq.111) then
            if(y.ge.2.1.and.y.le.2.9) then
              call vfill1(121,pt,wevt/(wp*pt))
            endif
          endif

c......Protons
        if(kfa.eq.2212) then
	  call vfill1(32,y,isign(1,kf)*wevt/wy)
	  if(y.ge.0.5d0.and.y.le.3.0d0)
     $ call vfill1(41,pt,isign(1,kf)*wevt/(wp*pt))
	  if(kf.eq.2212) call vfill1(31,y,wevt/wy)
        else if(kfa.eq.3122) then
	  call vfill1(34,y,isign(1,kf)*wevt/wy)
	  call vfill1(38,y,isign(1,kf)*wevt/wy)
        else if(kfa.eq.3112) then
	  call vfill1(35,y,isign(1,kf)*wevt/wy)
        else if(kfa.eq.3212) then
	  call vfill1(36,y,isign(1,kf)*wevt/wy)
	  call vfill1(38,y,isign(1,kf)*wevt/wy)
        else if(kfa.eq.3222) then
	  call vfill1(37,y,isign(1,kf)*wevt/wy)
        else if(kfa.eq.3322) then ! xi0
	  call vfill1(38,y,isign(1,kf)*wevt/wy)
        else if(kf.eq.321) then
	  call vfill1(39,y,wevt/wy)
        else if(kf.eq.-321) then
	  call vfill1(40,y,wevt/wy)
	endif


c.....Freezeout.
        if(abs(ycm).le.1.0d0) then
        iii=0
        if(kf.eq.2212.or.kf.eq.2112) iii=1
        if(kf.eq.3112.or.kf.eq.3212.or.kf.eq.3122.or.kf.eq.3222) iii=2
        if(kf.eq.111.or.kfa.eq.211) iii=3
        if(kfa.eq.311) iii=4
        if(kf.eq.321) iii=5
        if(kf.eq.-321) iii=6
        if(iii.ne.0.and.abs(k(7,i)).ne.1) then
              vrp=v(1,i)**2+v(2,i)**2+v(3,i)**2
              if(vrp.ge.0.0D0)then
                vr=sqrt(vrp)
              else
                goto 200
              endif
              call vfill1(80+iii,vr,1.D0) 
              call vfill1(90+iii,v(4,i),1.D0) 
        endif
        endif


C...End analysis of particle record.
  200 continue

      return
C***********************************************************************

      entry anal3

C...VNIBOOK: output of histograms.

      fac=1.d0
      mnorm=0
      mform=1

c...Et.
      call vscale(8,fac)
      call vprint(8,mnorm,mform)

      do i=1,10
	call vscale(30+i,fac)
	call vprint(30+i,mnorm,mform)
      end do
      do i=1,5
	call vscale(40+i,fac)
	call vprint(40+i,mnorm,mform)
      end do

c...WA80  EPJ  25%, nucl-ex/9805007
      sig=paru(1)*parc(4)**2
      fac1=sig*10/paru(2)
      do i=1,1
        call vscale(120+i,fac1)
        call vprint(120+i,mnorm,mform)
      end do

c...Freaze-out points
      do i=1,6
       call vprint(80+i,1,0)
       call vprint(90+i,1,0)
      end do


      end
