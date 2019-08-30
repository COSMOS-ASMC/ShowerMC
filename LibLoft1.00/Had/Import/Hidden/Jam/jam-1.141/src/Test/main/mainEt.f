c...Main program for p(14.6GeV/c)+Be
c...P.R. C45 (1992)2933 E-802 Collaboration

        include 'jam1.inc'
        include 'jam2.inc'
        character frame*8,proj*8,targ*8,cwin*15

c...Initialization
      fname(1)='jam.cfg'  ! input file name.
      mstc(1) =78827   ! random seed.
      mevent=1000       ! total simulation event
      bmin=0.0D0          ! minimum impact parameter
      dt=50.0D0           ! collision time(fm/c)
      nstep=1
      cwin='14.6gevc       '  ! incident energy
      frame='lab     '        ! comp. frame
      proj='p       '         ! projectile

      mstc(8)=0   ! job mode.
c     mstc(42)=0    ! weak decay include

c     targ='9Be     '         ! target
c     bmax=-2.3125      ! maximum impact parameter (168mb)
c     bmax=-4.0         ! maximum impact parameter (168mb)

      targ='27Al    '         ! target
      bmax=-3.6D0         ! maximum impact parameter (407mb)Al


      call jaminit(mevent,bmin,bmax,dt,nstep,
     &               frame,proj,targ,cwin)
      nevent=mstc(2)

c...Initialize analysis.
      call anal1

      do iev=1,nevent
        call jamevt(iev)
        if(mod(iev,500).eq.0)write(6,*)' event= ',iev
        if(mod(iev,10).eq.0) call jamlist(1)
        call anal2
      end do

      call jamfin
      call anal3

      end

c***********************************************************************

      subroutine anal1

      include 'jam1.inc'
      include 'jam2.inc'
      parameter(vcut=0.6D0,ecut=0.05D0)
 
      save wevt1,wevt2,fac
     
c...dSig/dEt
      emax=3.0D0
      emin=0.0D0
      we=0.05D0
      nemx=nint((emax-emin)/we)
      call vbook1(1,'dSig/dE 1',nemx,emin,emax)
      call vbook1(2,'dSig/dE 2',nemx,emin,emax)
      call vbook1(3,'dEt/deta all mesons',nemx,emin,emax)
      call vbook1(4,'dEt/deta pi0',nemx,emin,emax)
      call vbook1(5,'dEt/deta pion+-',nemx,emin,emax)

c...dEt/d(eta)
      ymax=4.0D0
      ymin=-1.0D0
      wy=0.1D0
      nymx=nint((ymax-ymin)/wy)
      call vbook1(11,'dEt/deta 1',nymx,ymin,ymax)
      call vbook1(12,'dEt/deta 2',nymx,ymin,ymax)

c...Eevent weight
      fac= 0.001D0*paru(1)*10*(parc(4)**2-parc(3)**2) ! barn
      wevt1=1.0D0/dble(mstc(2))/we
      wevt2=1.0D0/dble(mstc(2))/wy

      return

c*************************************************************************

        entry anal2

c...Loop over all particles
        do i=1,nv
         
         kf=k(2,i)
         kc=jamcomp(kf)
         if(kc.le.0.or.kc.gt.mstu(6))then
           write(6,*)'Invalide code at i, kf, kc : ',i,kf,kc
           go to 3000
         end if

c...Rapidity cut.
c          rap=0.5*log( max(p(4,i)+p(3,i),1.e-8)
c    $             /max(p(4,i)-p(3,i),1.e-8) )
c         if(mstc(4).eq.0) then
c         else if(mstc(4).eq.3) then
c         else
c          rap=rap+ylab
c         endif

        ptsq=p(1,i)**2+p(2,i)**2
        pt=sqrt(ptsq)
        pp=sqrt(ptsq+p(3,i)**2)
        pt=max(pt,1.D-8)
        ee=p(4,i)
        eta=0.5D0*log( max(pp+p(3,i),1.D-8)/max(pp-p(3,i),1.D-8) )
	et=ee*pt/max(pp,1.D-8)
        vv=pp/ee
        kch=jamchge(k(2,i))

c...Experimental cut.
        if(ee.lt.ecut) goto 3000
        if(vv.lt.vcut) goto 3000

        if(eta.ge.1.25D0.and.eta.le.2.5D0) then
        if(k(9,i).eq.0)  call vfill1(3,et,wevt1)
        if(kf.eq.111)  call vfill1(4,et,wevt1)
        if(abs(kf).eq.211) call vfill1(5,et,wevt1)
        endif

c...Gamma and charged particles.
        if(kf.eq.111.or.kch.ne.0) then

          call vfill1(12,eta,et*wevt2)
c...Rapidity cut.
          if(eta.ge.1.25D0.and.eta.le.2.5D0) then
            call vfill1(2,et,wevt1)
          endif

        endif

3000     end do

         return

c*************************************************************************

      entry anal3

c...Output hostogram

      do j=1,5
      call vscale(j,fac)
      call vprint(j,0,1)
      end do
      do j=1,2
      call vprint(10+j,0,1)
      end do

      end
