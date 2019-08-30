c...A main program for S(200GeV/c)+S

      include 'jam1.inc'
      include 'jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15

c....Initialize JAM.
      fname(1)='jam.cfg'  ! input file name.
      mstc(1) =48827      ! random seed.
      mevent=1            ! total simulation event
      bmin=0.0d0          ! minimum impact parameter
      bmax=-1.0d0         ! maximum impact parameter
      dt=100.0d0          ! collision time(fm/c)
      nstep=1
      cwin='200gevc        '  ! incident energy
      frame='nn      '        ! comp. frame
      proj='32S'         ! projectile
      targ='32S'         ! target
      mstc(8)=2   ! job mode.
      mstc(74)=0   ! dipole-approximated QCD radiation of the string
      mstc(81)=1   ! 1:hard scattering on/off
      mstc(156)=1  ! analysis of collision distribution

      call jaminit(mevent,bmin,bmax,dt,nstep,frame,proj,targ,cwin)
      nevent=mstc(2)

c...Initialize analysis.
      call anal1

       do i=1,2000
       write(31,8000)i,parf(i)
       enddo
 8000  format(i4,1x,f15.6)

c...Simulation start.
      do iev=1,nevent

c...Simulate one event.
        call jamevt(iev)

        if(mod(iev,100).eq.0) write(6,*)'event=',iev

c...Dump phase space data.
c       do i=1,nv
c       write(30,800)i,k(2,i),k(3,i),k(4,i),k(5,i),k(7,i),(v(j,i),j=1,5)
c       end do
c800    format(i4,1x,i7,1x,i6,1x,i11,1x,i11,i10,5(1x,g10.3))


c...Data analysis.
        call anal2

c...List phase space data.
c       call jamlist(1)

      end do

c...Final output.
      call jamfin

c...Print analysis results.
      call anal3

       do i=1,2000
       write(32,8000)i,parf(i)
       enddo

      write(6,*)'max res',mstd(199)
      write(6,*)'mstd(29)',mstd(29)
      print *,'mstd(198)',mstd(198)

      end

c***********************************************************************

      subroutine anal1

      include 'jam1.inc'
      include 'jam2.inc'
      dimension npa(0:20)
      save wy,wp,wm
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

c     ymin=ytarg*2.0d0
c     ymax=yproj*2.0d0
c     nymx=35
c     nymx=25 ! RHIC
c     nymx=30
c     wy=(ymax-ymin)/nymx
      ymin=-7.0d0
      ymax=7.0d0
      wy=0.15d0
c     wy=0.30d0
      nymx=(ymax-ymin)/wy

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

      call vbook1(18,'dN/dy - proj',nymx,ymin,ymax)
      call vbook1(19,'dN/dy - targ',nymx,ymin,ymax)

      call vbook1(51,'dN/dy - str',nymx,ymin,ymax)
      call vbook1(52,'dN/dy - res',nymx,ymin,ymax)


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

c...Mass distributins.
      rmin=0.8D0
      rmax=6.0D0
      wm=0.1D0
      nrmx=(rmax-rmin)/wm
      call vbook1(41,'dN/dmR proj.',nrmx,rmin,rmax)
      call vbook1(42,'dN/dmR targ.',nrmx,rmin,rmax)
      call vbook1(43,'dN/dmS proj.',nrmx,rmin,rmax)
      call vbook1(44,'dN/dmS targ.',nrmx,rmin,rmax)

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


c....Mass distibution
      kp1=mste(49)
      kt1=mste(50)
      em1=pare(49)
      em2=pare(50)
      if(kp1.eq.3) then
        call vfill1(43,em1,1.0d0/wm) 
      else
        call vfill1(41,em1,1.0d0/wm) 
      endif
      if(kt1.eq.3) then
        call vfill1(44,em2,1.0d0/wm) 
      else
        call vfill1(42,em2,1.0d0/wm) 
      endif

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

        rap=0.5d0*log( max(p(4,i)+p(3,i),1.d-8)/max(p(4,i)-p(3,i), 
     & 1.d-8) )

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

            if(k(7,i).le.-1) call vfill1(19,rap,1.D0/wy)
            if(k(7,i).ge.1) call vfill1(18,rap,1.D0/wy) 
            if(k(4,i).eq.92) then
              call vfill1(51,rap,1.d0/wy) 
            else
              call vfill1(52,rap,1.d0/wy) 
            endif


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


3000  end do
      call vfill1(32,dble(nch),1.d0) 

      return

c***********************************************************************

      entry anal3

c...Output of histograms.

c...Event weight
      fac=1.d0/dble(mstc(2))

c...Mass distributions.
      do i=1,4
       call vscale(40+i,fac)
       call vprint(40+i,0,0)
      end do

c...ET distributions.
      do i=3,4
       call vscale(30+i,fac)
       call vprint(30+i,0,0)
      end do

c...Rapidity distributions.
      do i=0,9
       call vscale(10+i,fac)
       call vprint(10+i,0,0)
      end do
      do i=1,2
       call vscale(50+i,fac)
       call vprint(50+i,0,0)
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
