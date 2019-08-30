c...Main program for Si(14.6GeV/c)+27Al,63Cu,197Au
c...P.R. C50(1994)1024 E-802 Collaboration
c...and for p(14.6GeV/c)+Be,27Al,63Cu,197Au
c...P.R. C45 (1992)2933 E-802 Collaboration

      include 'jam1.inc'
      include 'jam2.inc'

      character frame*8,proj*8,targ*8,cwin*15
      logical dump
c...ireac=1: Si+Al
c...ireac=2: Si+Cu
c...ireac=3: Si+Au
c...ireac=4: p+Be
c...ireac=5: p+Al
c...ireac=6: p+Cu
c...ireac=7: p+Au
      data dump/.false./
c     data dump/.true./
      data ireac/4/

      if(dump)
     $  open(33,file='phase.dat',form='unformatted',status='new')

c...Initialization.
c     mstc(1)=2882137  ! random seed.
      mstc(1)=2332337  ! random seed.
      bmin=0.0d0       ! minimum impact parameter
      dt=100.0d0       ! collision time(fm/c)
      nstep=1
      cwin='14.6gevc'  ! incident energy
      frame='nn'       ! comp. frame

      mstc(8)=3000     ! job mode.
      mstc(45)=0       ! deuteron
      mstc(156)=1      ! energy distribution of collisions
      parc(7)= 1.0d0   ! Output time interval (fm/c)

c     mstc(6)=-111     ! lead particle method.
c     mstc(6)=-1       ! Glauber
      parc(38)=0.10    ! low energy cut


      if(ireac.eq.1) then
c       mevent=10000      ! total simulation event
        mevent=1        ! total simulation event
        proj = '28Si'
        targ = '27Al'
        bmax=-1.797       ! 7% of sigma_INT
      else if(ireac.eq.2) then
        mevent=3000       ! total simulation event
        proj='28Si'
        targ='63Cu'
        bmax=-2.1988      ! 7% of 2.17barn E803
      else if(ireac.eq.3) then
        mevent=2000       ! total simulation event
        proj='28Si'
        targ='197Au'
        bmax=-2.89d0      ! 7% of 3.75barn E803
      else if(ireac.eq.4) then
c...p+9Be
        mevent=150000     ! total simulation event
        bmax=-4.0         ! maximum impact parameter
        proj='p'          ! projectile
        targ='9Be'        ! target

      else if(ireac.eq.5) then
c...p+27Al
        mevent=300000     ! total simulation event
        mevent=10000      ! total simulation event
c       bmax=-3.6d0       ! maximum impact parameter (407mb)Al
c       bmax=-3.767       ! maximum impact parameter (446mb)Al
        bmax=-4.5d0       ! maximum impact parameter (446mb)Al
        proj='p'          ! projectile
        targ='27Al'       ! target

      else if(ireac.eq.6) then
c...p+63Cu
        mevent=100000     ! total simulation event
        bmax=-4.8         ! maximum impact parameter (724mb)Cu
        bmax=-100.0
        proj='p'          ! projectile
        targ='63Cu'       ! target
      else if(ireac.eq.7) then
c...p+197Au
c       mevent=100000     ! total simulation event
        mevent=3000       ! total simulation event
ccc     bmax=-6.98        ! maximum impact parameter (1530mb)Au
        bmax=-7.50d0      ! maximum impact parameter (1530mb)Au
        bmax=-100.0
        proj='p'          ! projectile
        targ='197Au'      ! target
      else
        write(6,*)'Invalid ireac',ireac
        stop
      endif

c...Initialize JAM.
      call jaminit(mevent,bmin,bmax,dt,nstep,frame,proj,targ,cwin)

      if(dump)write(33)mevent,pard(17),pard(5),pard(6),mstc(4)

c...Initialize analysis.
      call ana1(ireac)

c...Loop over events.
      do iev=1,mevent

c....Simulate one event.
        call jamevt(iev)
c       if(mod(iev,500).eq.0)write(6,*)' event= ',iev
        if(mod(iev,50).eq.0)write(6,*)' event= ',iev

c...Analyzie data.
        call ana2

c...Dump phase space data if desired.
        if(dump) then
          write(33)iev,nv,nbary,nmeson,pard(2)
          do i=1,nv
            write(33)(k(j,i),j=1,11),(r(j,i),j=1,5),(p(j,i),j=1,5)
          end do
        endif

c...End loop over all events.
      end do

      if(dump) close(33)

c...Finish JAM.
      call jamfin

c..Output analysis results.
      call ana3(ireac)

      end

c***********************************************************************

      subroutine ana1(ireac)

      include 'jam1.inc'
      include 'jam2.inc'

      dimension dy(7,6),y0(7,6),ny(7,6)
      character cy1*3,cy2*3
      save wevt1,wevt2,ylab,wp,wevt
      save ihp,ihpip,ihpin,ihkp,ihkn,ihd
      save dyp,dypip,dypin,dykp,dykn,dyd
      save y0p,y0pip,y0pin,y0kp,y0kn,y0d
      save nyp,nypip,nypin,nykp,nykn,nyd
      save ihlast

c...Rapidity interval for central collisions.
c....                        p    pi+   pi-   k+    k-    d
      data (y0(1,i),i=1,6)/0.4d0,0.6d0,0.6d0,0.6d0,1.0d0,0.5d0/ !Si+Al
      data (y0(2,i),i=1,6)/0.4d0,0.6d0,0.6d0,0.6d0,1.0d0,0.4d0/ !Si+Cu
      data (y0(3,i),i=1,6)/0.5d0,0.7d0,0.7d0,0.7d0,0.9d0,0.5d0/ !Si+Au
      data (y0(4,i),i=1,6)/0.7d0,0.7d0,1.1d0,0.9d0,1.1d0,0.5d0/ !p+Be
      data (y0(5,i),i=1,6)/0.7d0,0.7d0,0.9d0,0.9d0,0.9d0,0.5d0/ !p+Al
      data (y0(6,i),i=1,6)/0.7d0,0.9d0,1.1d0,0.9d0,0.9d0,0.5d0/ !p+Cu
      data (y0(7,i),i=1,6)/0.7d0,0.7d0,0.7d0,0.7d0,0.9d0,0.5d0/ !p+Au

c...Rapidity bin for central collisions.
      data (dy(1,i),i=1,6)/0.2d0,0.2d0,0.2d0,0.4d0,0.4d0,0.2d0/
      data (dy(2,i),i=1,6)/0.2d0,0.2d0,0.2d0,0.2d0,0.4d0,0.2d0/
      data (dy(3,i),i=1,6)/0.2d0,0.2d0,0.2d0,0.2d0,0.4d0,0.2d0/
      data (dy(4,i),i=1,6)/0.2d0,0.2d0,0.2d0,0.2d0,0.2d0,0.2d0/
      data (dy(5,i),i=1,6)/0.2d0,0.2d0,0.2d0,0.2d0,0.2d0,0.2d0/
      data (dy(6,i),i=1,6)/0.2d0,0.2d0,0.2d0,0.2d0,0.2d0,0.2d0/
      data (dy(7,i),i=1,6)/0.2d0,0.2d0,0.2d0,0.2d0,0.2d0,0.2d0/

c...Number of bin for central collisions.
      data (ny(1,i),i=1,6)/9,12,12,5,3,5/  ! Si+Al
      data (ny(2,i),i=1,6)/9,12,12,7,3,6/  ! Si+Cu
      data (ny(3,i),i=1,6)/9,11,11,8,4,6/  ! Si+Au
      data (ny(4,i),i=1,6)/9,11, 9,7,5,5/  ! p+Be
      data (ny(5,i),i=1,6)/9,11,10,7,7,5/  ! p+Al
      data (ny(6,i),i=1,6)/9,10, 9,7,7,5/  ! p+Cu
      data (ny(7,i),i=1,6)/9,11,11,8,6,5/  ! p+Au


c Si+Al:    peripheral          central
c     p:  y=0.4-2.2 dy=0.2 #=10 y=0.4-2.0  dy=0.2 #= 9
c     pi: y=0.6-2.8 dy=0.2 #=12 y=0.6-2.8  dy=0.2 #=12
c     k+: y=0.8-2.0 dy=0.4 #= 4 y=0.6-2.2  dy=0.4 #= 5
c     k-: y=1.0-1.8 dy=0.4 #= 3 y=1.0-1.8  dy=0.4 #= 3
c     d:  y=0.4-1.0 dy=0.2 #= 4 y=0.5-1.2  dy=0.2 #= 5
c Si+Cu:    peripheral          central
c     p:  y=0.4-2.2 dy=0.2 #=10 y=0.4-2.0  dy=0.2 #= 9
c     pi: y=0.6-2.8 dy=0.2 #=12 y=0.6-2.8  dy=0.2 #=12
c     k+: y=0.8-2.0 dy=0.4 #= 4 y=0.6-1.8  dy=0.2 #= 7
c     k-: y=1.0-1.8 dy=0.4 #= 3 y=1.0-1.8  dy=0.4 #= 3
c     d:  y=0.4-1.0 dy=0.2 #= 4 y=0.4-1.4  dy=0.2 #= 6
c Si+Au:    peripheral          central
c     p:  y=0.5-2.1 dy=0.2 #= 9 y=0.5-2.1  dy=0.2 #= 9
c     pi: y=0.7-2.7 dy=0.2 #=11 y=0.7-2.7  dy=0.2 #=11
c     k+: y=0.9-2.1 dy=0.4 #= 4 y=0.7-2.1  dy=0.2 #= 8
c     k-: y=1.3-1.9 dy=0.6 #= 2 y=0.9-2.1  dy=0.4 #= 4
c     d:  y=0.5-1.1 dy=0.2 #= 4 y=0.5-1.5  dy=0.2 #= 6
     
      ylab=pard(17)
      wevt=1.0d0/dble(mstc(2))

      dyp=dy(ireac,1)
      dypip=dy(ireac,2)
      dypin=dy(ireac,3)
      dykp=dy(ireac,4)
      dykn=dy(ireac,5)
      dyd=dy(ireac,6)

      y0p=y0(ireac,1)
      y0pip=y0(ireac,2)
      y0pin=y0(ireac,3)
      y0kp=y0(ireac,4)
      y0kn=y0(ireac,5)
      y0d=y0(ireac,6)

      nyp=ny(ireac,1)
      nypip=ny(ireac,2)
      nypin=ny(ireac,3)
      nykp=ny(ireac,4)
      nykn=ny(ireac,5)
      nyd=ny(ireac,6)

c...dSig/dEt
      emax=4.0d0
      emin=0.0d0
      we=0.05d0
      nemx=nint((emax-emin)/we)
      call vbook1(1,'dSig/dE all    ',nemx,emin,emax)
      call vbook1(2,'dSig/dE charged',nemx,emin,emax)

c...dEt/d(eta)
      ymax=5.0d0
      ymin=-2.0d0
      wy=0.2d0
      nymx=nint((ymax-ymin)/wy)
      call vbook1(3,'dEt/deta all    ',nymx,ymin,ymax)
      call vbook1(4,'dEt/deta charged',nymx,ymin,ymax)

      call vbook1(11,'dN/dy proton',nymx,ymin,ymax)
      call vbook1(12,'dN/dy pi-',nymx,ymin,ymax)
      call vbook1(13,'dN/dy pi+',nymx,ymin,ymax)
      call vbook1(14,'dN/dy k-',nymx,ymin,ymax)
      call vbook1(15,'dN/dy k+',nymx,ymin,ymax)
      call vbook1(16,'dN/dy deuteon',nymx,ymin,ymax)

      pmin=0.0d0
      pmax=1.5d0
      wp=0.05d0
      npmx=nint((pmax-pmin)/wp)

c...Protons.
      ihist=20
      ihp=ihist
      y=y0p-dyp*3.d0/2.0d0
      do i=1,nyp
        ihist=ihist+1
        y=y+dyp
        write(cy1,'(f3.1)')y
        write(cy2,'(f3.1)')y+dyp
        call vbook1(ihist,'dN/dp proton'//cy1//'-'//cy2,npmx,pmin,pmax)
      end do

c...Positive pions.
      ihpip=ihist
      y=y0pip-dypip*3.d0/2.0d0
      do i=1,nypip
        ihist=ihist+1
        y=y+dypip
        write(cy1,'(f3.1)')y
        write(cy2,'(f3.1)')y+dypin
        call vbook1(ihist,'dN/dp pi+'//cy1//'-'//cy2,npmx,pmin,pmax)
      end do

c...Negative pions.
      ihpin=ihist
      y=y0pin-dypin*3.d0/2.0d0
      do i=1,nypin
        ihist=ihist+1
        y=y+dypin
        write(cy1,'(f3.1)')y
        write(cy2,'(f3.1)')y+dypin
        call vbook1(ihist,'dN/dp pi-'//cy1//'-'//cy2,npmx,pmin,pmax)
      end do

      pmin=0.0d0
      pmax=1.3d0
      npmx=nint((pmax-pmin)/wp)
      ihkp=ihist
c...Positive kaons.
      y=y0kp-dykp*3.d0/2.0d0
      do i=1,nykp
        ihist=ihist+1
        y=y+dykp
        write(cy1,'(f3.1)')y
        write(cy2,'(f3.1)')y+dykp
        call vbook1(ihist,'dN/dp k+'//cy1//'-'//cy2,npmx,pmin,pmax)
      end do

c...Negative kaons.
      ihkn=ihist
      y=y0kn-dykn*3.d0/2.0d0
      do i=1,nykn
        ihist=ihist+1
        y=y+dykn
        write(cy1,'(f3.1)')y
        write(cy2,'(f3.1)')y+dykn
        call vbook1(ihist,'dN/dp k-'//cy1//'-'//cy2,npmx,pmin,pmax)
      end do

c...Deuterons.
      ihd=ihist
      y=y0d-dyd*3.d0/2.0d0
      do i=1,nyd
        ihist=ihist+1
        y=y+dyd
        write(cy1,'(f3.1)')y
        write(cy2,'(f3.1)')y+dyd
        call vbook1(ihist,'dN/dp deut.'//cy1//'-'//cy2,npmx,pmin,pmax)
      end do

      ihlast=ihist

c...Eevent weight
      wevt1=1.0d0/dble(mstc(2))/we
      wevt2=1.0d0/dble(mstc(2))/wy

      return

c*************************************************************************

        entry ana2

c...Loop over all particles
        do i=1,nv
         
         if(k(1,i).ge.11) goto 3000
         kf=k(2,i)
         if(k(1,i).le.4) then
         kc=jamcomp(kf)
         if(kc.le.0.or.kc.gt.mstu(6))then
c          write(6,*)'Invalide code at i, kf, kc : ',i,kf,kc
           go to 3000
         end if
         endif

c...Rapidity cut.
           rap=0.5d0*log( max(p(4,i)+p(3,i),1.d-8)
     $             /max(p(4,i)-p(3,i),1.d-8) )
          if(mstc(4).eq.0) then
          else if(mstc(4).eq.3) then
          else
           rap=rap+ylab
          endif

        ptsq=p(1,i)**2+p(2,i)**2
        pt=sqrt(ptsq)
        pp=sqrt(ptsq+p(3,i)**2)
        pt=max(pt,1.d-8)
        eta=0.5d0*log( max(pp+p(3,i),1.d-8)/max(pp-p(3,i),1.d-8) )
	et=p(4,i)*pt/max(pp,1.d-8)
        emt=sqrt(p(5,i)**2+ptsq)
        emt0=emt-p(5,i)
        ppl=ptsq+p(3,i)**2
        if(mstc(4).eq.2.or.mstc(4).eq.3)then
          bet=pard(5)
          gam=pard(6)
          pl=gam*(p(3,i)+bet*p(4,i))
          el=gam*(p(4,i)+bet*p(3,i))
c         yl=0.5d0*log( max(el+pl,1.d-8)/max(el-pl,1.d-8) )
          ppl=ptsq+pl**2
        endif
        ppl=sqrt(ppl)
c...deuteron
         if(kf.eq.1001001000) goto 2000

        kch=jamchge(k(2,i))

        if(eta.ge.1.25d0.and.eta.le.2.5d0) then
          call vfill1(1,et,wevt1)
          if(kch.ne.0) call vfill1(2,et,wevt1)
        endif

c...Rapidity and transverse momentum distributions.
        if(kf.eq.2212) then
c         if(ppl.gt.0.3d0) then
          call vfill1(11,rap,wevt/wy)
          ih=ihp
          ytag=y0p-dyp*3.d0/2.0d0
          do ib=1,nyp
            ih=ih+1
            ytag=ytag+dyp
            if(rap.ge.ytag.and.rap.le.ytag+dyp)
     $              call vfill1(ih,emt0,wevt/(wp*emt*dyp))
          end do
c         endif
        else if(kf.eq.-211) then
c         if(ppl.gt.0.05d0) then
          call vfill1(12,rap,wevt/wy)
          ih=ihpin
          ytag=y0pin-dypin*3.d0/2.0d0
          do ib=1,nypin
            ih=ih+1
            ytag=ytag+dypin
            if(rap.ge.ytag.and.rap.le.ytag+dypin)
     $              call vfill1(ih,emt0,wevt/(wp*emt*dypin))
          end do
c         endif
        else if(kf.eq.211) then
c         if(ppl.gt.0.05d0) then
          call vfill1(13,rap,wevt/wy)
          ih=ihpip
          ytag=y0pip-dypip*3.d0/2.0d0
          do ib=1,nypip
            ih=ih+1
            ytag=ytag+dypip
            if(rap.ge.ytag.and.rap.le.ytag+dypip)
     $              call vfill1(ih,emt0,wevt/(wp*emt*dypip))
          end do
c         endif
        else if(kf.eq.-321) then
          call vfill1(14,rap,wevt/wy)
          ih=ihkn
          ytag=y0kn-dykn*3.d0/2.0d0
          do ib=1,nykn
            ih=ih+1
            ytag=ytag+dykn
            if(rap.ge.ytag.and.rap.le.ytag+dykn)
     $              call vfill1(ih,emt0,wevt/(wp*emt*dykn))
          end do
        else if(kf.eq.321) then
          call vfill1(15,rap,wevt/wy)
          ih=ihkp
          ytag=y0kp-dykp*3.d0/2.d0
          do ib=1,nykp
            ih=ih+1
            ytag=ytag+dykp
            if(rap.ge.ytag.and.rap.le.ytag+dykp)
     $              call vfill1(ih,emt0,wevt/(wp*emt*dykp))
          end do
        endif

        if(kf.eq.2212.or.kf.eq.2112) goto 3000
        if(eta.ge.1.25d0.and.eta.le.2.5d0) then
          call vfill1(1,et,wevt1)
          if(kch.ne.0) call vfill1(2,et,wevt1)
        endif
        call vfill1(3,eta,et*wevt2)
        if(kch.ne.0) call vfill1(4,eta,wevt2)
        goto 3000

c...deuteron
2000     continue
         call vfill1(16,rap,wevt/wy)
         ih=ihd
          ytag=y0d-dyd*3.d0/2.d0
          do ib=1,nyd
            ih=ih+1
            ytag=ytag+dyd
            if(rap.ge.ytag.and.rap.le.ytag+dykp)
     $        call vfill1(ih,emt0,wevt/(wp*emt*dyd))
          end do


3000     end do

         return

c*************************************************************************

      entry ana3(ireac)

c...Output hostogram

c...Rapidity
      if(ireac.le.3) then
        fac=1.0d0
      else
        fac=1.0d0
c       fac=paru(1)*10*(parc(4)**2-parc(3)**2)
      endif

      do j=1,6
      call vscale(10+j,fac)
      call vprint(10+j,0,1)
      end do

      do j=1,2
      call vscale(j,fac)
      call vprint(j,0,1)
      end do
      do j=3,4
      call vscale(j,1.0d0)
      call vprint(j,0,1)
      end do


c...Mt distributions.
      if(ireac.le.3) then
        fac=1/paru(2)
      else
        fac= paru(1)*10*(parc(4)**2-parc(3)**2)/paru(2)
      endif

      ihist=20
 40   ihist=ihist+1
      if(ihist.gt.ihlast) goto 50
      call vscale(ihist,fac)
      call vprint(ihist,0,1)
      goto 40
 50   continue

      end
