c...A main program for calculation of et and mt at SPS.

      include 'jam1.inc'
      include 'jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15
      logical dump 
c     data bmin/0.0/,bmax/-1.08/,dt/100.0/,mevent/1/
c     data dump/.false./
      data dump/.true./

      if(dump)
     $  open(33,file='phase.dat',form='unformatted',status='new')

c....Initialize JAM.
      mstc(1) =38827   ! random seed.
      mevent=1000      ! total simulation event
      bmin=0.0D0          ! minimum impact parameter
      bmax=-3.2D0        ! NA49
c     bmax=-4.478d0      ! WA98 10%
      dt=100.0D0          ! collision time(fm/c)
      nstep=1
      cwin='158gev'  ! incident energy
      frame='nn'        ! comp. frame
      proj='208Pb'         ! projectile
      targ='208Pb'         ! target


c...WA80
c     proj='32S'         ! projectile
c     targ='32S'         ! target
c     cwin='200gev'      ! incident energy
c     bmax=-3.4          ! 25% central

      mstc(8)=0    ! job mode.

c     mstc(156)=1  ! analysis of collision distribution
c     mstc(155)=0  ! flow anal.
c     mstc(163)=0  ! time evolution of directed transverse flow
c     parc(7)= 1.0D0    ! Output time interval (fm/c)

      mstc(51)=0  ! BB collisions only

c     mstc(6)=-1  ! Glauber
c     mstc(81)= 0 ! <--- 1:Hard scattering on/off
c     mstc(156)=1        ! analysis of collision distribution
c     mstc(156)=1        ! energy distribution of collisions
c     mstc(162)=1        ! Output collision histroy
c     mstc(165)=1        ! Output time evol. of particle yield
c     mstc(166)=1        ! Output time evol. of particle density
c     mstc(167)=1        ! Output time evol. of density (Gaussian)
c     parc(7)= 0.5D0     ! Output time interval (fm/c)

c     parc(54)=1.6d0     ! string tension
c     parc(54)=2.0d0     ! string tension

      call jaminit(mevent,bmin,bmax,dt,nstep,
     $                             frame,proj,targ,cwin)
      nevent=mstc(2)
      if(dump)write(33)nevent,pard(17),pard(5),pard(6),mstc(4)

c...Initialize analysis.
      call ana1

c...Simulation start.
      do iev=1,nevent

c...Simulate one event.
        call jamevt(iev)

c...Dump phase space data.
        if(dump) then
          write(33)iev,nv,nbary,nmeson,pard(2)
          do i=1,nv
            write(33)(k(j,i),j=1,11),(r(j,i),j=1,5),(p(j,i),j=1,5)
     $              ,(v(j,i),j=1,5)
          end do
        endif

        if(mod(iev,100).eq.0) write(6,*)'event=',iev

c...Data analysis.
        call ana2

c...List phase space data.
c       call jamlist(1)

      end do

      if(dump) close(33)

c...Final output.
      call jamfin

c...Print analysis results.
      call ana3

      end

c***********************************************************************

      subroutine ana1

      include 'jam1.inc'
      include 'jam2.inc'
      common/myana/ispec
      save wy,wp,ylab
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

      call vbook1(1,'dn/dy protons',nymx,ymin,ymax)
      call vbook1(2,'dn/dy antiprotons',nymx,ymin,ymax)
      call vbook1(3,'dn/dy net protons',nymx,ymin,ymax)
      call vbook1(4,'dn/dy pi-',nymx,ymin,ymax)
      call vbook1(5,'dn/dy pi+',nymx,ymin,ymax)
      call vbook1(6,'dn/dy k-',nymx,ymin,ymax)
      call vbook1(7,'dn/dy k+',nymx,ymin,ymax)
      call vbook1(8,'dn/dy lambda',nymx,ymin,ymax)
      call vbook1(9,'dn/dy h-(k pi pbar)',nymx,ymin,ymax)

      call vbook1(19,'dETdy',nymx,ymin,ymax)
      call vbook1(20,'dETdy',nymx,ymin,ymax)

      call vbook1(21,'1/mtdn/dmtdy protons',npmx,pmin,pmax)
      call vbook1(22,'1/mtdn/dmtdy antiprotons',npmx,pmin,pmax)
      call vbook1(23,'1/mtdn/dmtdy net potons',npmx,pmin,pmax)
      call vbook1(24,'1/mtdn/dmtdy pion-',npmx,pmin,pmax)
      call vbook1(25,'1/mtdn/dmtdy pion+',npmx,pmin,pmax)
      call vbook1(26,'1/mtdn/dmtdy k-',npmx,pmin,pmax)
      call vbook1(27,'1/mtdn/dmtdy k+',npmx,pmin,pmax)
      call vbook1(28,'1/mtdn/dmtdy lambda',npmx,pmin,pmax)
      call vbook1(29,'1/mtdn/dmtdy h-',npmx,pmin,pmax)

      return

c***********************************************************************

      entry ana2

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

        if(mstc(4).eq.0) then
        else if(mstc(4).eq.3) then
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


        if(kf.eq.-211.or.kf.eq.-321.or.kf.eq.-2212) then
	  call vfill1(9,yl,1.0d0/wy)
	  if(abs(y).le.0.5d0) then
	    call vfill1(29,emt0,1.d0/(emt*wp*wy))
	  endif
	endif

c....Protons.
        if(kf.eq.2212) then
	  call vfill1(1,yl,1.0d0/wy)
	  call vfill1(3,yl,1.0d0/wy)
	  if(abs(y).le.0.5d0) then
	    call vfill1(21,emt0,1.d0/(emt*wp*wy))
	    call vfill1(23,emt0,1.d0/(emt*wp*wy))
	  endif
c.....anti protons.
        else if(kf.eq.-2212) then
	  call vfill1(2,yl,1.0d0/wy)
	  call vfill1(3,yl,-1.0d0/wy)
	  if(abs(y).le.0.5d0) then
	    call vfill1(22,emt0,1.d0/(emt*wp*wy))
	    call vfill1(23,emt0,-1.d0/(emt*wp*wy))
	  endif
c.....pion-
        else if(kf.eq.-211) then
          call vfill1(4,yl,1.0d0/wy)
	  if(abs(y).le.0.5d0) then
	    call vfill1(24,emt0,1.d0/(emt*wp*wy))
	  endif
c.....pion+
        else if(kf.eq.211) then
          call vfill1(5,yl,1.0d0/wy)
	  if(abs(y).le.0.5d0) then
	    call vfill1(25,emt0,1.d0/(emt*wp*wy))
	  endif
c.....kaon-
        else if(kf.eq.-321) then
          call vfill1(6,yl,1.0d0/wy)
	  if(abs(y).le.0.5d0) then
	    call vfill1(26,emt0,1.d0/(emt*wp*wy))
	  endif
c.....kaon+
        else if(kf.eq.321) then
          call vfill1(7,yl,1.0d0/wy)
	  if(abs(y).le.0.5d0) then
	    call vfill1(27,emt0,1.d0/(emt*wp*wy))
	  endif
c.....lambda
        else if(kf.eq.3122) then
          call vfill1(8,yl,1.0d0/wy)
	  if(abs(y).le.0.5d0) then
	    call vfill1(28,emt0,1.d0/(emt*wp*wy))
	  endif
        endif
3000  continue
      xspec=xspec+ispec

      return

c***********************************************************************

      entry ana3

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
      do i=1,9
      call vscale(i,fac)
      call vprint(i,0,0)
      call vscale(20+i,fac)
      call vprint(20+i,mnorm,mform)
      end do

      end

