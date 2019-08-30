c...A main program for S(200GeV/c)+S

      include 'jam1.inc'
      include 'jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15

c....Initialize JAM.
c     fname(1)='jam.cfg'  ! input file name.
      mstc(1) =48827      ! random seed.
      mevent=100          ! total simulation event
      bmin=0.0d0          ! minimum impact parameter
      bmax=-1000.0d0         ! maximum impact parameter
      dt=100.0d0          ! collision time(fm/c)
      nstep=1
      cwin='200gevc        '  ! incident energy
      frame='nn      '        ! comp. frame
      proj='p'         ! projectile
c     targ='p'         ! target
c     proj='32S'         ! target
c     proj='k0'         ! projectile
      targ='32S'         ! target
c     targ='p'         ! target
c     proj='197Au'         ! projectile
c     targ='197Au'         ! target
c     mstc(8)=2    ! job mode.
c     mstc(74)=0   ! dipole-approximated QCD radiation of the string
c     mstc(81)=1   ! 1:hard scattering on/off
c     mstc(156)=1  ! analysis of collision distribution

c     call xxx(cwin,proj,targ,frame)

      mstc(42)=0 ! allow weak decay after simulation.
      mdcy(jamcomp(111),1)=1    ! pi0
      mdcy(jamcomp(311),1)=1    ! k0
      mdcy(jamcomp(-311),1)=1   ! ak0 
      mdcy(jamcomp(310),1)=1    ! k0_S

      mdcy(jamcomp(3122),1)=1   ! Lambda0
      mdcy(jamcomp(-3122),1)=1
      mdcy(jamcomp(3112),1)=1   ! Sigma-
      mdcy(jamcomp(-3112),1)=1
      mdcy(jamcomp(3212),1)=1   ! Sigma0
      mdcy(jamcomp(-3212),1)=1
      mdcy(jamcomp(3222),1)=1   ! Sigma+
      mdcy(jamcomp(-3222),1)=1
      mdcy(jamcomp(3312),1)=1   ! Xi-
      mdcy(jamcomp(-3312),1)=1
      mdcy(jamcomp(3322),1)=1   ! Xi0
      mdcy(jamcomp(-3322),1)=1

      call jaminit(mevent,bmin,bmax,dt,nstep,frame,proj,targ,cwin)
      nevent=mstc(2)

c...Initialize analysis.
      call anal1

c...Simulation start.
      do iev=1,nevent

c...Simulate one event.
        call jamevt(iev)

        if(mod(iev,100).eq.0) write(6,*)'event=',iev

c...Dump phase space data.
c       do i=1,nv
c       write(30,800)i,k(2,i),k(3,i),k(4,i),k(5,i),k(7,i),(v(j,i),j=1,5)
c800    format(i4,1x,i7,1x,i6,1x,i11,1x,i11,i10,5(1x,g10.3))
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c        if(r(4,i).lt.abs(r(3,i))) then
c           write(77,*)i,k(1,i),k(2,i),k(4,i),r(4,i),r(3,i)
c        endif
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c       end do


c...Data analysis.
        call anal2

c...List phase space data.
c       call jamlist(1)

      end do

c...Final output.
      call jamfin

c...Print analysis results.
      call anal3

      write(6,*)'max res',mstd(199)
      write(6,*)'mstd(29)',mstd(29)

      end

c***********************************************************************

      subroutine anal1

      include 'jam1.inc'
      include 'jam2.inc'
      dimension npa(0:20)
      save wy,wp
      save ylab,yproj,ytarg
      save npa,ispec

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
      ispec=0

c....Rapidity distribution.
      ylab=pard(17)
      dbetpr=pard(35)
      dbetta=pard(45)
      yproj=0.5d0*log((1.0d0+dbetpr)/(1.0d0-dbetpr))
      ytarg=0.5d0*log((1.0d0+dbetta)/(1.0d0-dbetta))

      ymin=ytarg*2.0d0
      ymax=yproj*2.0d0

      if(pard(16).le.40.d0) then
        wy=0.25d0
c       wy=0.2d0
        nymx=(ymax-ymin)/wy
      else
        nymx=30
        ymin=-7.0d0
        ymax=7.0d0
        wy=(ymax-ymin)/nymx
      endif


      if(mstc(4).eq.0) then
      else if(mstc(4).eq.3) then
      else
       ymax=ymax+ylab
       ymin=ymin+ylab
      endif

      call vbook1(10,'dN/dy - pi- k- pbar',nymx,ymin,ymax)
      call vbook1(11,'dN/dy - proton',nymx,ymin,ymax)
      call vbook1(12,'dN/dy - pion- ',nymx,ymin,ymax)
      call vbook1(13,'dN/dy - pion+ ',nymx,ymin,ymax)
      call vbook1(14,'dN/dy - charged',nymx,ymin,ymax)
      call vbook1(15,'dN/dy - p-p-bar',nymx,ymin,ymax)
      call vbook1(16,'dN/dy - h-',nymx,ymin,ymax)
      call vbook1(17,'dN/dy - net baryon',nymx,ymin,ymax)


c...Transverse distributions.

      if(pard(16).le.40.d0) then
        pmin=0.0d0
        pmax=3.0d0
        npmx=30
      else
        pmin=0.0d0
        pmax=20.0d0
        npmx=30
      endif
      wp=(pmax-pmin)/npmx
      call vbook1(21,'dN/dpt**2 - proton ',npmx,pmin,pmax)
      call vbook1(22,'dN/dpt**2 - pion-  ',npmx,pmin,pmax)
      call vbook1(23,'dN/dpt**2 - pion+  ',npmx,pmin,pmax)
      call vbook1(24,'dN/dpt**2 - charged',npmx,pmin,pmax)

      call vbook1(31,'dN/d(eta) - charged',nymx,ymin,ymax)
      call vbook1(32,'P(n) - charged',50,0.d0,100.d0)
c...Transverse energy distributions.
      call vbook1(33,'dET/dy - all',nymx,ymin,ymax)
      call vbook1(34,'dET/dy - charged',nymx,ymin,ymax)


c....Freaze-out point
      call vbook1(81,'dN/dr proton',40,0.d0,30.d0)
      call vbook1(82,'dN/dr pion+',40,0.d0,30.d0)
      call vbook1(83,'dN/dr pion-',40,0.d0,30.d0)
      call vbook1(84,'dN/dr kaon+',40,0.d0,30.d0)

c....Freaze-out time
      call vbook1(91,'dN/dt proton',40,0.d0,30.d0)
      call vbook1(92,'dN/dt pion+',40,0.d0,30.d0)
      call vbook1(93,'dN/dt pion-',40,0.d0,30.d0)
      call vbook1(94,'dN/dt kaon+',40,0.d0,30.d0)

      return

c***********************************************************************

      entry anal2

      nch=0
      neg=0
c...Loop over all particles.
      do 3000 i=1,nv

        kf=k(2,i)
        if(kf.ge.1000000000) goto 3000
        kc=jamcomp(kf)
        if(kc.le.0.or.kc.gt.mstu(6)) then
           write(6,*)'Invalid code i kf kc',i,kf,kc,nv,nbary,nmeson
           goto 3000
        endif

        rap=0.5d0*log( max(p(4,i)+p(3,i),1.d-8)/max(p(4,i)-p(3,i), 
     & 1.d-8) )

c....Spect
         if(abs(k(7,i)).eq.1) ispec=ispec+1

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

c      if((k(2,i).eq.2212..or.k(2,i).eq.2112).and.p(5,i).ge.1.0) then
c        write(6,*)'prot???',i,k(1,i),k(2,i),p(5,i)
c        pause
c      endif

        npa(0)=npa(0)+1
        call vfill1(33,rap,et/wy)

        kch=jamchge(k(2,i))
c...Charged particles.
        if(kch.ne.0) then
          nch=nch+1
          npa(1)=npa(1)+1
          call vfill1(14,rap,1.d0/wy)
          call vfill1(31,eta,1.d0/wy)
          call vfill1(24,pt,1.d0/wp/pt/2) 
          call vfill1(34,rap,et/wy)

c...Negative charged particles.
          if(kch.lt.0) then
            neg=neg+1
            npa(2)=npa(2)+1
            call vfill1(16,rap,1.d0/wy)
          endif
        endif

c...Net baryons.
        ibar=kchg(kc,6)
        if(ibar.eq.3) then
          if(kf.gt.0) then
            call vfill1(17,rap,1.d0/wy) 
          else if(kf.lt.0) then
            call vfill1(17,rap,-1.d0/wy) 
          endif
        endif

c...h-(pi-,k-,p~)
        if(kf.eq.-211.or.kf.eq.-321.or.kf.eq.-2212) then
           call vfill1(10,rap,1.d0/wy)
        endif


c.......Protons.
        if(abs(kf).eq.2212) then
          if(kf.eq.2212) then
            npa(12)=npa(12)+1
            call vfill1(15,rap,1.d0/wy) 
            call vfill1(11,rap,1.d0/wy) 
            call vfill1(21,pt,1.d0/wp/pt/2) 

            if(abs(k(7,i)).ne.1) then
              vr=sqrt(v(1,i)**2+v(2,i)**2+v(3,i)**2)
              call vfill1(81,vr,1.d0) 
              call vfill1(91,v(4,i),1.d0) 
            endif

          else if(kf.eq.-2212) then
            npa(13)=npa(13)+1
            call vfill1(15,rap,-1.d0/wy) 
          endif

c.......Pions.
        else if(kf.eq.-211) then
          npa(3)=npa(3)+1
          call vfill1(12,rap,1.d0/wy) 
          call vfill1(22,pt,1.d0/wp/pt/2)

          vr=sqrt(v(1,i)**2+v(2,i)**2+v(3,i)**2)
          call vfill1(82,vr,1.d0) 
          call vfill1(92,v(4,i),1.d0) 

        else if(kf.eq.111) then
          npa(4)=npa(4)+1
        else if(kf.eq.211) then

          npa(5)=npa(5)+1
          call vfill1(13,rap,1.d0/wy) 
          call vfill1(23,pt,1.d0/wp/pt/2)

          vr=sqrt(v(1,i)**2+v(2,i)**2+v(3,i)**2)
          call vfill1(83,vr,1.d0) 
          call vfill1(93,v(4,i),1.d0) 

        else if(kf.eq.321) then

          vr=sqrt(v(1,i)**2+v(2,i)**2+v(3,i)**2)
          call vfill1(84,vr,1.d0) 
          call vfill1(94,v(4,i),1.d0) 

        endif

        if(kf.eq.3122) then
          npa(8)=npa(8)+1
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

        else if(kf.eq.321.or.kf.eq.311) then
          npa(7)=npa(7)+1
        else if(kf.eq.-321.or.kf.eq.-311) then
          npa(6)=npa(6)+1
        endif


3000  continue
      call vfill1(32,dble(nch),1.d0) 

      return

c***********************************************************************

      entry anal3

c...Output of histograms.

c...Event weight
      fac=1.d0/dble(mstc(2)*mstc(5))

c...ET distributions.
      do i=3,4
       call vscale(30+i,fac)
       call vprint(30+i,0,0)
      end do

c...Rapidity distributions.
      do i=0,7
       call vscale(10+i,fac)
       call vprint(10+i,0,0)
      end do

c...Pt distributions.
      do i=1,4
       call vscale(20+i,fac)
       call vprint(20+i,0,1)
      end do

c...dn/d(eta),P(n)
      do i=1,2
       call vscale(30+i,fac)
       call vprint(30+i,0,1)
      end do

c...Freaze-out points
      do i=1,4
       call vprint(80+i,1,0)
       call vprint(90+i,1,0)
      end do

      open(70,file='file70',status='unknown')
      write(70,*)'ylab yproj ytarg=',ylab,yproj,ytarg
      write(70,*)'# of initial particle',mstd(11),
     $  'spec',ispec*fac,' participant',mstd(11)-ispec*fac
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
