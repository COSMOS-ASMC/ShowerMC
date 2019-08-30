c...A main program for Au(11.6GeV/c)+Au P.R.C57 (1998) R466
c...calculate dn/dy, transverse momentum distritutions.

      include 'jam1.inc'
      include 'jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15

      mstc(1) =1766777   ! random seed.
      mevent=1           ! total simulation event
c     bmin=0.0D0         ! minimum impact parameter
c     bmax=-3.338D0      ! maximum impact parameter 350mb
      cwin='11.6gevc'    ! incident energy
      frame='nn'         ! comp. frame
      proj='197Au'       ! projectile
      targ='197Au'       ! target

c...cascade mode.
c     dt=100.0D0         ! collision time(fm/c)
c     nstep=1

c...BUU simulation
      dt=0.2d0
      nstep=150
      mstc(5)=30         ! number of test particle
      mstc(6)=12         ! buu mode
      mstc(101)=1        ! coulomb
      bmin=2.0d0
      bmax=2.0d0

      mstc(8)=1          ! job mode.
c     mstc(45)=1         ! deuteron coalecence

      mstc(155)=0        ! flow anal.
      mstc(156)=1        ! analysis of collision distribution
      mstc(156)=1        ! energy distribution of collisions
      mstc(162)=1        ! Output collision histroy
      mstc(165)=1        ! Output time evol. of particle yield
      mstc(166)=0        ! Output time evol. of particle density
      mstc(167)=0        ! Output time evol. of density (Gaussian)
      parc(7)= 1.0D0     ! Output time interval (fm/c)

c....Initialize JAM.
      call jaminit(mevent,bmin,bmax,dt,nstep,frame,proj,targ,cwin)

c...Initialize analysis.
      call anal1

c...Simulation loop start.
      do iev=1,mevent

c...Simulate one JAM event.
        call jamevt(iev)

        if(mod(iev,100).eq.0) write(6,*)'event=',iev

c...Data analysis.
        call anal2

      end do

c...Final output.
      call jamfin

c...Print analysis results.
      call anal3

      end

c***********************************************************************

      subroutine anal1

      include 'jam1.inc'
      include 'jam2.inc'
      dimension npa(0:20)
      save wy,wp
      save ylab,yproj,ytarg
      save npa
      save dely,dely2

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
      if(mstc(2).le.200) then
        ymin=ytarg*2.0D0
        ymax=yproj*2.0D0
        nymx=30
        wy=(ymax-ymin)/nymx
      else
        ymin=-4.0D0
        ymax=4.0D0
        wy=0.1D0
        nymx=(ymax-ymin)/wy
      endif

c     if(mstc(4).eq.0) then
c     else if(mstc(4).eq.3) then
c     else
c      ymax=ymax+ylab
c      ymin=ymin+ylab
c     endif


      call vbook1(11,'dN/dy - proton',nymx,ymin,ymax)
      call vbook1(12,'dN/dy - pion- ',nymx,ymin,ymax)
      call vbook1(13,'dN/dy - pion+ ',nymx,ymin,ymax)
      call vbook1(14,'dN/dy - lambda ',nymx,ymin,ymax)
      call vbook1(15,'dN/dy - kaon- ',nymx,ymin,ymax)
      call vbook1(16,'dN/dy - kaon+ ',nymx,ymin,ymax)
      call vbook1(17,'dN/dy - lambda/sigma0 ',nymx,ymin,ymax)
      call vbook1(18,'dN/dy - Y ',nymx,ymin,ymax)

c...Transverse distributions.

      if(pard(16).le.40.D0) then
        if(mstc(2).le.200) then
          pmin=0.0D0
          pmax=3.0D0
          npmx=30
          wp=(pmax-pmin)/npmx
        else
          pmin=0.0D0
          pmax=3.0D0
          wp=0.03D0
          npmx=(pmax-pmin)/wp
        endif
      else
        pmin=0.0D0
        pmax=20.0D0
        npmx=30
        wp=(pmax-pmin)/npmx
      endif

      dely=0.1D0
      call vbook1(21,'1/(2pim_t)dN/dm_tdy - p y=0.05D0',npmx,pmin,pmax)
      call vbook1(22,'1/(2pim_t)dN/dm_tdy - p y=0.15D0',npmx,pmin,pmax)
      call vbook1(23,'1/(2pim_t)dN/dm_tdy - p y=0.25D0',npmx,pmin,pmax)
      call vbook1(24,'1/(2pim_t)dN/dm_tdy - p y=0.35D0',npmx,pmin,pmax)
      call vbook1(25,'1/(2pim_t)dN/dm_tdy - p y=0.45D0',npmx,pmin,pmax)
      call vbook1(26,'1/(2pim_t)dN/dm_tdy - p y=0.55D0',npmx,pmin,pmax)
      call vbook1(27,'1/(2pim_t)dN/dm_tdy - p y=0.65D0',npmx,pmin,pmax)
      call vbook1(28,'1/(2pim_t)dN/dm_tdy - p y=0.75D0',npmx,pmin,pmax)
      call vbook1(29,'1/(2pim_t)dN/dm_tdy - p y=0.85D0',npmx,pmin,pmax)
      call vbook1(30,'1/(2pim_t)dN/dm_tdy - p y=0.95D0',npmx,pmin,pmax)
      call vbook1(31,'1/(2pim_t)dN/dm_tdy - p y=1.05D0',npmx,pmin,pmax)

      dely2=0.2D0
      call vbook1(32,'1/(2pim_t)dN/dm_tdy p 0<y<0.2D0',npmx,pmin,pmax)
      call vbook1(33,'1/(2pim_t)dN/dm_tdy pi- 0<y<0.2D0',npmx,pmin,pmax)
      call vbook1(34,'1/(2pim_t)dN/dm_tdy pi+ 0<y<0.2D0',npmx,pmin,pmax)


c....Freaze-out point
      call vbook1(81,'dN/dr proton',40,0.D0,80.D0)
      call vbook1(82,'dN/dr pion  ',40,0.D0,80.D0)
      call vbook1(83,'dN/dr akaon ',40,0.D0,80.D0)
      call vbook1(84,'dN/dr kaon  ',40,0.D0,80.D0)
      call vbook1(85,'dN/dr lambda',40,0.D0,80.D0)
      call vbook1(86,'dN/dr Y     ',40,0.D0,80.D0)

c....Freaze-out time
      call vbook1(91,'dN/dt proton',40,0.D0,80.D0)
      call vbook1(92,'dN/dt pion  ',40,0.D0,80.D0)
      call vbook1(93,'dN/dt akaon ',40,0.D0,80.D0)
      call vbook1(94,'dN/dt kaon  ',40,0.D0,80.D0)
      call vbook1(95,'dN/dt lambda',40,0.D0,80.D0)
      call vbook1(96,'dN/dt Y     ',40,0.D0,80.D0)

      return

c***********************************************************************

      entry anal2

      nch=0
      neg=0
c...Loop over all particles.
      do i=1,nv

        kf=k(2,i)
        kc=jamcomp(kf)
        if(kc.le.0.or.kc.gt.mstu(6)) then
           write(6,*)'Invalid code i kf kc',i,kf,kc,nv,nbary,nmeson
           goto 3000
        endif

        rap=0.5D0*log( max(p(4,i)+p(3,i),1.D-8)/max(p(4,i)-p(3,i), 
     & 1.D-8) )
        ycm=rap

        if(mstc(4).eq.0) then
        else if(mstc(4).eq.3) then
        else
         rap=rap+ylab
        endif

        ptsq=p(1,i)**2+p(2,i)**2
        pt=sqrt(ptsq)
        pp=sqrt(ptsq+p(3,i)**2)
        pt=max(pt,1.D-8)
        emt=sqrt(p(5,i)**2+ptsq)
        emt0=sqrt(p(5,i)**2+ptsq)-p(5,i)
        weiy=1.D0/(2*paru(1)*emt*dely*wp)
        weiy2=1.D0/(2*paru(1)*emt*dely2*wp)

        npa(0)=npa(0)+1

c...Charged particles.
        kch=jamchge(k(2,i))
        if(kch.ne.0) then
          nch=nch+1
          npa(1)=npa(1)+1

c...Negative charged particles.
          if(kch.lt.0) then
            neg=neg+1
            npa(2)=npa(2)+1
          endif
        endif

c------- Freeze out points -----------------------------------
c...Pions
      if(abs(kf).eq.211.or.kf.eq.111) then
          vr=sqrt(v(1,i)**2+v(2,i)**2+v(3,i)**2)
          call vfill1(82,vr,1.D0) 
          call vfill1(92,v(4,i),1.D0) 
c...Nucleons
      else if(kf.eq.2112.or.kf.eq.2212) then
         if(abs(k(7,i)).ne.1) then
              vr=sqrt(v(1,i)**2+v(2,i)**2+v(3,i)**2)
              call vfill1(81,vr,1.D0) 
              call vfill1(91,v(4,i),1.D0) 
         endif
c...Kaons
      else if(kf.eq.311.or.kf.eq.321) then
          vr=sqrt(v(1,i)**2+v(2,i)**2+v(3,i)**2)
          call vfill1(84,vr,1.D0) 
          call vfill1(94,v(4,i),1.D0) 
c...Anti-kaons
      else if(kf.eq.-311.or.kf.eq.-321) then
          vr=sqrt(v(1,i)**2+v(2,i)**2+v(3,i)**2)
          call vfill1(83,vr,1.D0) 
          call vfill1(93,v(4,i),1.D0) 
c....Hyperons.
      else if(kf.eq.3122.or.kf.eq.3212.or.kf.eq.3222.or.kf.eq.3112) then
          call vfill1(18,ycm,1.D0/wy) 
          if(kf.eq.3122.or.kf.eq.3212) call vfill1(17,ycm,1.D0/wy) 
         if(abs(k(7,i)).ne.1) then
              vr=sqrt(v(1,i)**2+v(2,i)**2+v(3,i)**2)
              call vfill1(86,vr,1.D0) 
              call vfill1(96,v(4,i),1.D0) 
              if(kf.eq.3122) then
                call vfill1(85,vr,1.D0) 
                call vfill1(95,v(4,i),1.D0) 
              endif
        endif
      endif
c------------------------------------------------------------

c.......Protons.
        if(abs(kf).eq.2212) then
          if(kf.eq.2212) then
            npa(12)=npa(12)+1
            call vfill1(11,ycm,1.D0/wy) 

            if(ycm.lt.0.0D0) then
            else if(ycm.le.0.1D0) then
              call vfill1(21,emt0,weiy) 
            else if(ycm.le.0.2D0) then
              call vfill1(22,emt0,weiy) 
            else if(ycm.le.0.3D0) then
              call vfill1(23,emt0,weiy) 
            else if(ycm.le.0.4D0) then
              call vfill1(24,emt0,weiy) 
            else if(ycm.le.0.5D0) then
              call vfill1(25,emt0,weiy) 
            else if(ycm.le.0.6D0) then
              call vfill1(26,emt0,weiy) 
            else if(ycm.le.0.7D0) then
              call vfill1(27,emt0,weiy) 
            else if(ycm.le.0.8D0) then
              call vfill1(28,emt0,weiy) 
            else if(ycm.le.0.9D0) then
              call vfill1(29,emt0,weiy) 
            else if(ycm.le.1.0D0) then
              call vfill1(30,emt0,weiy) 
            else if(ycm.le.1.1D0) then
              call vfill1(31,emt0,weiy) 
            endif

            if(ycm.ge.0.0D0.and.ycm.le.0.2D0) then
              call vfill1(32,emt0,weiy2) 
            endif
          else if(kf.eq.-2212) then
            npa(13)=npa(13)+1
          endif

c.......Pions.
        else if(kf.eq.-211) then
          npa(3)=npa(3)+1
          call vfill1(12,ycm,1.D0/wy) 
          if(ycm.ge.0.0D0.and.ycm.le.0.2D0) call vfill1(33,emt0,weiy2)
        else if(kf.eq.111) then
          npa(4)=npa(4)+1
        else if(kf.eq.211) then
          npa(5)=npa(5)+1
          call vfill1(13,ycm,1.D0/wy) 
          if(ycm.ge.0.0D0.and.ycm.le.0.2D0) call vfill1(34,emt0,weiy2)
c...Kaons
        else if(kf.eq.321) then
          call vfill1(16,ycm,1.D0/wy) 
        else if(kf.eq.-321) then
          call vfill1(15,ycm,1.D0/wy) 
c...Lambda
        else if(kf.eq.3122) then
          npa(8)=npa(8)+1
          call vfill1(14,ycm,1.D0/wy) 
c...Anti-Lambda
        else if(kf.eq.-3122) then
          npa(9)=npa(9)+1
        else if(kf.eq.3112) then ! Sigma-
          npa(10)=npa(10)+1
        else if(kf.eq.3212) then ! Sigma0
          npa(10)=npa(10)+1
        else if(kf.eq.3222) then ! Sigma+
          npa(10)=npa(10)+1
        else if(kf.eq.-3112) then ! a-Sigma-
          npa(11)=npa(11)+1
        else if(kf.eq.-3212) then ! a-Sigma0
          npa(11)=npa(11)+1
        else if(kf.eq.-3222) then ! a-Sigma+
          npa(11)=npa(11)+1
        endif

        if(kf.eq.321.or.kf.eq.311) then
          npa(7)=npa(7)+1
        else if(kf.eq.-321.or.kf.eq.-311) then
          npa(6)=npa(6)+1
        endif


3000  end do
      call vfill1(32,dble(nch),1.D0) 

      return

c***********************************************************************

      entry anal3

c...Output of histograms.

c...Event weight
      fac=1.D0/dble(mstc(2))

c...Rapidity distributions.
      do i=1,8
       call vscale(10+i,fac)
       call vprint(10+i,0,0)
      end do

c...Mt distributions.
      do i=1,14
       call vscale(20+i,fac)
       call vprint(20+i,0,1)
      end do

c...Freaze-out points
      do i=1,6
       call vprint(80+i,1,0)
       call vprint(90+i,1,0)
      end do

      open(70,file='file70',status='unknown')
      write(70,*)'ylab yproj ytarg=',ylab,yproj,ytarg
      write(70,*)'average mult',npa(0)*fac
      write(70,*)'charged',npa(1)*fac
      write(70,*)'negative',npa(2)*fac
      write(70,*)'pi- pi0 pi+',npa(3)*fac,npa(4)*fac,npa(5)*fac
      write(70,*)'pion  total',(npa(3)+npa(4)+npa(5))*fac
      write(70,*)'proton total',npa(12)*fac
      write(70,*)'a-proton total',npa(13)*fac
      write(70,*)'lambda a-lam total',npa(8)*fac,npa(9)*fac
      write(70,*)'sima a-sigma total',npa(10)*fac,npa(11)*fac

      write(70,*)'kaon   total',npa(7)*fac
      write(70,*)'akaon  total',npa(6)*fac
      write(70,*)'average number of jet',pard(87)*fac
      close(70)

      end

