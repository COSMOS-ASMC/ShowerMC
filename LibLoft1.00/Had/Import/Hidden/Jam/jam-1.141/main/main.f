c...A main program to use the initial condition of hadronic cascade
c...from statistical model.

      include 'jam1.inc'
      include 'jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15

c     do i=1,200
c      print *,'mstc',i,mstc(i),parc(i)
c     end do
c     stop

c=========Set input values and switches ========================
c....Initialize JAM
      mstc(1)=48127      ! random seed.
      mstc(8)=0          ! job mode.
      mstc(16)=0         ! display on/off.
      parc(6)=5.0        ! scale of display
      mstc(54)=0         !avoid first coll inside the same nucleus off

c....Switch on some analysis.
      mstc(156)=1        ! analysis of collision distribution
      mstc(162)=1        ! Output collision histroy
      mstc(165)=1        ! 
      parc(7)= 1.0D0     ! Output time interval (fm/c)
      mstc(81)=0         ! 1:hard scattering on/off

c....Initial setting for JAM.
      mevent=1           ! total simulation event
      frame='user'       ! comp. frame in this case, user defined 
      dt=100.d0          ! collision time(fm/c)
      nstep=1            ! time step (i.e. no time step)
      cwin='200gev'      ! initial c.m. energy per nucl, in this case,
                         ! most energitic two-body collisions expected
c...dummy in this case.
      bmin=0.0d0         ! minimum impact parameter (dummy)
      bmax=0.0d0         ! maximum impact parameter (dummy)
c     proj='209Pb'       ! projectile
c     targ='209Pb'       ! target
c     proj='179Au'       ! projectile
c     targ='179Au'       ! target
      proj='32S'         ! projectile
      targ='32S'         ! target

c================ end input section ==================================
      
c...Initialize jam.
      call jaminit(mevent,bmin,bmax,dt,nstep,frame,proj,targ,cwin)
      nevent=mstc(2)

c...Initialize analysis.
      call anal1

c...Simulation start.
      do iev=1,nevent

c...Sampling particle momentum and coordinate.
        call init_mom

        if(mod(iev,100).eq.0) write(6,*)'event=',iev

c...Simulate one event.
           call jamevt(iev)    !...generate hadronic cascade
           call anal2          !...Data analysis.

c....output phase space data.
cc      write(30,*)'iev=',iev,nv
cc      do i=1,nv
cc      write(30,800)k(1,i),k(2,i),(r(j,i),j=1,5),(v(j,i),j=1,5)
cc      if(k(1,10).ge.10) then
cc       write(20,*)'nev=',iev,' i=',i
cc       write(20,*)k(1,i),k(2,i),(r(j,i),j=1,4),(p(j,i),j=1,5)
cc      endif
cc      end do
cc800   format(i4,1x,i8,4(g12.4,1x),/,13x,5(g12.4))

      end do

c...Final output.
        call jamfin

c...Print analysis results.
        call anal3

c...print out initial phase space.
      fac=1d0/mstc(2)
c...Rapidty distributuons.
      do i=1,10
      call vscale(100+i,fac)
      call vprint(100+i,0,0)
      call vscale(120+i,fac)
      call vprint(120+i,0,1)
      end do

      end

c******************************************************************

      function isoget(ir)
      implicit double precision(a-h, o-z)
c...Convert hydro particle code to JAM code. Isospin is also generated.
      
      if(ir.eq.1) then        ! a0
        if(rn(0).lt.0.333333) then
          isoget=10111
        else if(rn(0).lt.0.666666) then
          isoget=10211
        else
          isoget=-10211
        endif
      else if(ir.eq.2.or.ir.eq.3) then   ! Delta
        if(rn(0).lt.0.25) then
          isoget=1114
        else if(rn(0).lt.0.5) then
          isoget=2114
        else if(rn(0).lt.0.75) then
          isoget=2214
        else
          isoget=2224
        endif
        if(ir.eq.3) isoget=-1*isoget
      else if(ir.eq.4) then   ! eta
        isoget=221
      else if(ir.eq.5) then   ! eta'
        isoget=331
      else if(ir.eq.6) then   ! f_0(980)
        isoget=10221
      else if(ir.eq.7) then   ! K0
        isoget=311
        if(rn(0).lt.0.5) isoget=-311
      else if(ir.eq.8) then   ! K*(892)
        isoget=313
        if(rn(0).lt.0.5) isoget=-313
      else if(ir.eq.9) then   ! Lambda
        isoget=3122
      else if(ir.eq.10) then   ! Lambda-bar
        isoget=-3122
      else if(ir.eq.11) then   ! omega
        isoget=223
      else if(ir.eq.12) then  ! phi(1020)
        isoget=333
      else if(ir.eq.13) then  ! rho(770)
        if(rn(0).lt.0.333333) then
          isoget=113
        else if(rn(0).lt.0.666666) then
          isoget=213
        else
          isoget=-213
        endif
      else if(ir.eq.14) then  ! f_0(600) i.e. sigma meson
        isoget=10220
      else if(ir.eq.15.or.ir.eq.16) then  ! Sigma
        if(rn(0).lt.0.333333) then
          isoget=3112
        else if(rn(0).lt.0.666666) then
          isoget=3212
        else
          isoget=3222
        endif
        if(ir.eq.16) isoget= -1*isoget
      else if(ir.eq.17) then  ! pion
        if(rn(0).lt.0.333333) then
          isoget=111
        else if(rn(0).lt.0.666666) then
          isoget=211
        else
          isoget=-211
        endif
      else if(ir.eq.18) then  ! K
        isoget=321
        if(rn(0).lt.0.5) isoget=-321
      else if(ir.eq.19.or.ir.eq.20) then  ! nucleon
        isoget=2212  ! proton
        if(rn(0).lt.0.5) isoget=2112 ! n
        if(ir.eq.20) isoget= -1*isoget
      else
        print *,'funny ir=',ir
        isoget=0
      endif

      end

c******************************************************************

      subroutine init_mom

c...Sample particle accroding to the statisitcal model initial
c...condition calculated in subr.init_part.
      include 'jam1.inc'
      include 'jam2.inc'
      real*8 jamdtim
      common/myana2/wy,wp,ylab
      character cdummy*80

      open(11,file='../data/particlesample_pos.dat',status='old')
      read(11,'(a)')cdummy
      ip=0
      nbary=0
      nmeson=0
1000  continue
      read(11,*,end=888,err=999) px,py,pz,e,em,ir,tau,rx,ry,eta
      px=px/1000
      py=py/1000
      pz=pz/1000
      kf=isoget(ir+1)   ! PDG particle code.
      if(kf.eq.0) goto 1000
      kc=jamcomp(kf)  ! internal particle code.
      ibary=isign(kchg(kc,6),kf)  ! baryon number
      if(ibary.eq.0) then
        nmeson=nmeson+1
      else
       nbary=nbary+1
      endif

      ip=ip+1

c...Zero the vector.
        call jamzero(ip)

c....Particle mass.
        pm=pjmass(kf)
        if(pmas(kc,2).le.1d-7.or.mdcy(kc,1).eq.0
     $              .or.mdcy(kc,2).eq.0.or.mdcy(kc,3).eq.0)then
          k(1,ip)=1
        else
          k(1,ip)=2
        endif
        k(2,ip)=kf
        k(3,ip)=0
        k(4,ip)=0
        k(5,ip)=-1
        k(6,ip)=0
        k(7,ip)=1
        k(8,ip)=1
        k(9,ip)=ibary
        k(10,ip)=0
        k(11,ip)=0
        p(1,ip)=px
        p(2,ip)=py
        p(3,ip)=pz
        p(4,ip)=sqrt(px**2+py**2+pz**2+pm**2)
        p(5,ip)=pm

        r(1,ip)=rx
        r(2,ip)=ry
        r(3,ip)=tau*sinh(eta)
        r(4,ip)=tau*cosh(eta)
        r(5,ip)=r(4,ip)           ! formation time

c...Vertex
        v(1,ip)=r(1,ip)
        v(2,ip)=r(2,ip)
        v(3,ip)=r(3,ip)
        v(4,ip)=r(4,ip)

c.....Set resonance decay time.
        v(5,ip)=1.d+35
        if(k(1,ip).eq.2) then
          v(5,ip)=r(4,ip)+jamdtim(1,kf,kc,k(1,ip),p(5,ip),p(4,ip))
        endif

      goto 1000
999   continue
      write(6,*) 'file read error'
888   continue
      write(6,*) 'read finished, file close'
      close(11)
      nv=ip   ! total number of particles
      print *,'total particle=',nv

c...C.M.correction.
      icm=0
      if(icm.eq.1) then
      cx=0.d0
      cy=0.d0
      cz=0.d0
      px=0.d0
      py=0.d0
      pz=0.d0
      s=0.d0
      do i=1,nv
        px=px+p(1,i)
        py=py+p(2,i)
        pz=pz+p(3,i)
        cx=cx+r(1,i)*p(5,i)
        cy=cy+r(2,i)*p(5,i)
        cz=cz+r(3,i)*p(5,i)
        s=s+p(5,i)
      end do

      cx=-cx/s
      cy=-cy/s
      cz=-cz/s
      px=-px/nv
      py=-py/nv
      pz=-pz/nv
      endif

      do i=1,nv
        if(icm.eq.1) then
          r(1,i)=r(1,i)+cx
          r(2,i)=r(2,i)+cy
          r(3,i)=r(3,i)+cz
          v(1,i)=r(1,i)
          v(2,i)=r(2,i)
          v(3,i)=r(3,i)
          p(1,i)=p(1,i)+px
          p(2,i)=p(2,i)+py
          p(3,i)=p(3,i)+pz
          p(4,i)=sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
        endif

        pt=sqrt(p(1,i)**2+p(2,i)**2)
        y=0.5d0*log(max(p(4,i)+p(3,i),1.d-8)/
     $              max(p(4,i)-p(3,i),1.d-8))

        call vfill1(101,y,1.0d0/wy)
        call vfill1(121,pt,1.d0/(pt*wp))

        kf=k(2,i)
        iii=0
        if(kf.eq.2212)  iii=2                       ! Protons.
        if(kf.eq.-2212) iii=3                       ! anti protons.
        if(abs(kf).eq.211.or.kf.eq.111)  iii=5      ! pi
        if(abs(kf).eq.321.or.abs(kf).eq.-311) iii=6 ! k
        if(kf.eq.3122)  iii=7                       ! lambda
        if(kf.eq.3112.or.kf.eq.3212.or.kf.eq.3222) iii=8  ! sigma
        if(kf.eq.3312.or.kf.eq.3312)  iii=9         ! xi
        if(kf.eq.3334)  iii=10                      ! omega

        if(iii.ne.0) then
          call vfill1(iii+100,y,1.0d0/wy)
          call vfill1(iii+120,pt,1.d0/(pt*wp))
        endif
        if(iii.eq.2.or.iii.eq.3) then
          call vfill1(104,y,isign(1,kf)/wy)
          call vfill1(124,pt,isign(1,kf)/(pt*wp))
        endif

 300  end do

      end

c***********************************************************************
c...This routine will be called in every collision ans decay.
      subroutine jamanaus1(indd,npart)

c...Analysize collision spectra.
      include 'jam1.inc'
      include 'jam2.inc'
      dimension indd(100)

c...Usuful vectors for analysis.
      ichanel=mste(1)
      icltyp=mste(2)

      kf1=kcp(1,1)   ! PDG code for colliding particle 1
      kf2=kcp(1,2)   ! PDG code for colliding particle 2
      kf3=kcp(2,1)   ! PDG code for outgoing  particle 1
      kf4=kcp(2,2)   ! PDG code for outgoing  particle 2

c...Decay of resonance or string fragmentation.
      if(icltyp.eq.-1) then
      endif

      end

c***********************************************************************

      subroutine anal1

      include 'jam1.inc'
      include 'jam2.inc'
      common/myana1/ispec
      common/myana2/wy,wp,ylab
      save xspec

      xspec=0

c....Rapidity distribution.
      ylab=pard(17)
      print *,'ylab',ylab
      ymin=-7.0d0
      ymax=7.0d0
      wy=0.25d0
      nymx=(ymax-ymin)/wy

      yminl=ymin+ylab
      ymaxl=ymax+ylab
      print *,'yminl ymaxl ylab',yminl,ymaxl,ylab

      pmin=0.0d0
      pmax=5.0d0
c     wp=0.1d0
      wp=0.2d0
      npmx=(pmax-pmin)/wp

      call vbook1(1,'dn/dy total',nymx,ymin,ymax)
      call vbook1(2,'dn/dy protons',nymx,ymin,ymax)
      call vbook1(3,'dn/dy antiprotons',nymx,ymin,ymax)
      call vbook1(4,'dn/dy net protons',nymx,ymin,ymax)
      call vbook1(5,'dn/dy pi',nymx,ymin,ymax)
      call vbook1(6,'dn/dy k',nymx,ymin,ymax)
      call vbook1(7,'dn/dy lambda',nymx,ymin,ymax)
      call vbook1(8,'dn/dy sigma',nymx,ymin,ymax)
      call vbook1(9,'dn/dy xi',nymx,ymin,ymax)
      call vbook1(10,'dn/dy Omega',nymx,ymin,ymax)

      call vbook1(19,'dETdy',nymx,ymin,ymax)
      call vbook1(20,'dETdy',nymx,ymin,ymax)

      call vbook1(21,'1/p/dpdy total',npmx,pmin,pmax)
      call vbook1(22,'1/p/dpdy protons',npmx,pmin,pmax)
      call vbook1(23,'1/p/dpdy antiprotons',npmx,pmin,pmax)
      call vbook1(24,'1/p/dpdy net potons',npmx,pmin,pmax)
      call vbook1(25,'1/p/dpdy pion',npmx,pmin,pmax)
      call vbook1(26,'1/p/dpdy k',npmx,pmin,pmax)
      call vbook1(27,'1/p/dpdy lambda',npmx,pmin,pmax)
      call vbook1(28,'1/p/dpdy sigma',npmx,pmin,pmax)
      call vbook1(29,'1/p/dpdy xi',npmx,pmin,pmax)
      call vbook1(30,'1/p/dpdy Omega',npmx,pmin,pmax)


c...for initial condition.
      call vbook1(101,'dn/dy total',nymx,ymin,ymax)
      call vbook1(102,'dn/dy protons',nymx,ymin,ymax)
      call vbook1(103,'dn/dy antiprotons',nymx,ymin,ymax)
      call vbook1(104,'dn/dy net protons',nymx,ymin,ymax)
      call vbook1(105,'dn/dy pi',nymx,ymin,ymax)
      call vbook1(106,'dn/dy k',nymx,ymin,ymax)
      call vbook1(107,'dn/dy lambda',nymx,ymin,ymax)
      call vbook1(108,'dn/dy sigma',nymx,ymin,ymax)
      call vbook1(109,'dn/dy xi',nymx,ymin,ymax)
      call vbook1(110,'dn/dy Omega',nymx,ymin,ymax)

      call vbook1(121,'1/p/dpdy total',npmx,pmin,pmax)
      call vbook1(122,'1/p/dpdy protons',npmx,pmin,pmax)
      call vbook1(123,'1/p/dpdy antiprotons',npmx,pmin,pmax)
      call vbook1(124,'1/p/dpdy net potons',npmx,pmin,pmax)
      call vbook1(125,'1/p/dpdy pion',npmx,pmin,pmax)
      call vbook1(126,'1/p/dpdy k',npmx,pmin,pmax)
      call vbook1(127,'1/p/dpdy lambda',npmx,pmin,pmax)
      call vbook1(128,'1/p/dpdy sigma',npmx,pmin,pmax)
      call vbook1(129,'1/p/dpdy xi',npmx,pmin,pmax)
      call vbook1(130,'1/p/dpdy Omega',npmx,pmin,pmax)

      return

c***********************************************************************

      entry anal2

      beta=pard(5)
      gamma=pard(6)
      ispec=0
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
        y=0.5d0*log( max(p(4,i)+p(3,i),1.d-8)/max(p(4,i)-p(3,i),1.d-8) )
        yl=y
        if(mstc(4).eq.0) then
        else if(mstc(4).eq.3.or.mstd(4).eq.100) then
        else
         yl=y+ylab
        endif

        ptsq=p(1,i)**2+p(2,i)**2
        pt=sqrt(ptsq)
        pp=sqrt(ptsq+p(3,i)**2)
	ee=sqrt(p(5,i)**2+pp**2)
	if(abs(ee-p(4,i)).gt.1e-7) then
	  print *,'ee p4',ee,p(4,i),kf
	endif
        pt=max(pt,1.d-8)
        eta=0.5d0*log( max(pp+p(3,i),1.d-8)/max(pp-p(3,i),1.d-8) )
	et=p(4,i)*pt/max(pp,1.d-8)
        emt=sqrt(p(5,i)**2+ptsq)
        emt0=emt-p(5,i)

        plab  = gamma*( p(3,i) + beta * p(4,i) )
        elab  = gamma*( p(4,i) + beta * p(3,i) )
	ppl=sqrt(plab**2+ptsq)
	etlab=elab*pt/max(ppl,1d-8)
	eta_lab=0.5d0*log(max(ppl+plab,1.d-8)/max(ppl-plab,1.d-8))
	etl=eta_lab+ylab
	call vfill1(19,etl,etlab/wy)
	call vfill1(20,y,et/wy)
c       print *,'gam bet',gamma,beta
c       print *,'plab elab',plab,elab,p(3,i),p(4,i)
c       print *,'etlab eta ylab',etlab,eta_lab,ylab
c       pause


        call vfill1(1,yl,1.0d0/wy)
        call vfill1(21,pt,1.d0/(pt*wp))

        iii=0
        if(kf.eq.2212)  iii=2    !...Protons.
        if(kf.eq.-2212) iii=3    !...anti protons.
        if(abs(kf).eq.211.or.kf.eq.111)  iii=5    ! pi
        if(abs(kf).eq.321.or.abs(kf).eq.-311)  iii=6  ! k
        if(kf.eq.3122)  iii=7  ! lambda
        if(kf.eq.3112.or.kf.eq.3212.or.kf.eq.3222)  iii=8  ! sigma
        if(kf.eq.3312.or.kf.eq.3312)  iii=9  ! xi
        if(kf.eq.3334)  iii=10  ! omega

        if(iii.ne.0) then
	  call vfill1(iii,yl,1.0d0/wy)
          call vfill1(iii+20,pt,1.d0/(pt*wp))
        endif
        if(iii.eq.2.or.iii.eq.3) then
	  call vfill1(4,yl,isign(1,kf)/wy)
          call vfill1(24,pt,isign(1,kf)/(pt*wp))
        endif

3000  continue
      xspec=xspec+ispec

      return

c***********************************************************************

      entry anal3

c...Output histograms.

      print *,'event',mstc(2)
      wevt=1.d0/dble(mstc(2))
      fac=wevt
      mnorm=0
      mform=1

      do i=1,2
      call vscale(18+i,fac)
      call vprint(18+i,0,1)
      end do

c...Rapidty distributuons.
      do i=1,10
      call vscale(i,fac)
      call vprint(i,0,0)
      call vscale(20+i,fac)
      call vprint(20+i,0,1)
      end do



      end

