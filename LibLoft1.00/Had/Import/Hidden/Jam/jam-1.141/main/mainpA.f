c...Main program for p(14.6GeV/c)+Be,27Al,63Cu,197Au
c...P.R. C45 (1992)2933 E-802 Collaboration

      include 'jam1.inc'
      include 'jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15

c...Initialization
c     mstc(1) =48827   ! random seed.
      mstc(1) =28821   ! random seed.
      bmin=0.0D0          ! minimum impact parameter
      dt=100.0D0          ! collision time(fm/c)
      nstep=1
      cwin='14.6gevc       '  ! incident energy
      frame='nn      '        ! comp. frame
      proj='p       '         ! projectile

c...9Be
c     mevent=300000     ! total simulation event
c     bmax=-3.0         ! maximum impact parameter
c     targ='9Be     '   ! target

c...27Al
c     mevent=300000     ! total simulation event
c     bmax=-3.6         ! maximum impact parameter (407mb)Al
c     bmax=-3.767       ! maximum impact parameter (446mb)Al
c     bmax=-100         ! min. bias
c     targ='27Al    '   ! target

c...63Cu
c     mevent=100000     ! total simulation event
c     bmax=-4.8         ! maximum impact parameter (724mb)Cu
c     targ='63Cu    '   ! target
c...197Au
      mevent=100000     ! total simulation event
ccc   bmax=-6.98        ! maximum impact parameter (1530mb)Au
      bmax=-7.50D0        ! maximum impact parameter (1530mb)Au
      targ='197Au   '   ! target

      mstc(8)=0   ! job mode.
      mstc(156)=1  ! energy distribution of collisions
      parc(7)= 1.0D0 ! Output time interval (fm/c)
c     mstc(55)= 1  ! 1:frozen resonance
c     mstc(51)= 1  ! 0:only BB collision.

cproj=28Si
ctarg=63Cu
cwin=14.6gevc
cbmin=0.0
cbmax=-2.1988  # 7% of 2.17barn E803
c
cproj=28Si
ctarg=197Au
cwin=14.6gevc
cbmin=0.0
cbmax=-2.89  # 7% of 3.75barn E803

c...Initialize JAM.
      call jaminit(mevent,bmin,bmax,dt,nstep,frame,proj,targ,cwin)
      nevent=mstc(2)

c...Initialize analysis.
      call anal1

      do iev=1,nevent
        call jamevt(iev)
        if(mod(iev,500).eq.0)write(6,*)' event= ',iev
        call anal2
      end do

      call jamfin
      call anal3

      end

c***********************************************************************

      subroutine anal1

      include 'jam1.inc'
      include 'jam2.inc'
 
      save wevt1,wevt2,ylab,wp,wevt
     
      ylab=pard(17)
      wevt=1.0D0/dble(mstc(2))

c...dSig/dEt
      emax=4.0D0
      emin=0.0D0
      we=0.05D0
      nemx=nint((emax-emin)/we)
      call vbook1(1,'dSig/dE all    ',nemx,emin,emax)
      call vbook1(2,'dSig/dE charged',nemx,emin,emax)

c...dEt/d(eta)
      ymax=5.0D0
      ymin=-2.0D0
      wy=0.2D0
      nymx=nint((ymax-ymin)/wy)
      call vbook1(3,'dEt/deta all    ',nymx,ymin,ymax)
      call vbook1(4,'dEt/deta charged',nymx,ymin,ymax)

      call vbook1(11,'dN/dy proton',nymx,ymin,ymax)
      call vbook1(12,'dN/dy pi-',nymx,ymin,ymax)
      call vbook1(13,'dN/dy pi+',nymx,ymin,ymax)
      call vbook1(14,'dN/dy k-',nymx,ymin,ymax)
      call vbook1(15,'dN/dy k+',nymx,ymin,ymax)

      pmin=0.0D0
      pmax=1.5D0
      wp=0.05D0
      npmx=nint((pmax-pmin)/wp)
      call vbook1(21,'dN/dp proton 0.6D0-0.8D0',npmx,pmin,pmax)
      call vbook1(22,'dN/dp proton 0.8D0-1.0D0',npmx,pmin,pmax)
      call vbook1(23,'dN/dp proton 1.0D0-1.2D0',npmx,pmin,pmax)
      call vbook1(24,'dN/dp proton 1.2D0-1.4D0',npmx,pmin,pmax)
      call vbook1(25,'dN/dp proton 1.4D0-1.6D0',npmx,pmin,pmax)
      call vbook1(26,'dN/dp proton 1.6D0-1.8D0',npmx,pmin,pmax)
      call vbook1(27,'dN/dp proton 1.8D0-2.0D0',npmx,pmin,pmax)
      call vbook1(28,'dN/dp proton 2.0D0-2.2D0',npmx,pmin,pmax)
      call vbook1(29,'dN/dp proton 2.2D0-2.4D0',npmx,pmin,pmax)

      call vbook1(31,'dN/dp pion-  0.6D0-0.8D0',npmx,pmin,pmax)
      call vbook1(32,'dN/dp pion-  0.8D0-1.0D0',npmx,pmin,pmax)
      call vbook1(33,'dN/dp pion-  1.0D0-1.2D0',npmx,pmin,pmax)
      call vbook1(34,'dN/dp pion-  1.2D0-1.4D0',npmx,pmin,pmax)
      call vbook1(35,'dN/dp pion-  1.4D0-1.6D0',npmx,pmin,pmax)
      call vbook1(36,'dN/dp pion-  1.6D0-1.8D0',npmx,pmin,pmax)
      call vbook1(37,'dN/dp pion-  1.8D0-2.0D0',npmx,pmin,pmax)
      call vbook1(38,'dN/dp pion-  2.0D0-2.2D0',npmx,pmin,pmax)
      call vbook1(39,'dN/dp pion-  2.2D0-2.4D0',npmx,pmin,pmax)
      call vbook1(40,'dN/dp pion-  2.4D0-2.6D0',npmx,pmin,pmax)
      call vbook1(41,'dN/dp pion-  2.6D0-2.8D0',npmx,pmin,pmax)

      call vbook1(42,'dN/dp pion+  0.6D0-0.8D0',npmx,pmin,pmax)
      call vbook1(43,'dN/dp pion+  0.8D0-1.0D0',npmx,pmin,pmax)
      call vbook1(44,'dN/dp pion+  1.0D0-1.2D0',npmx,pmin,pmax)
      call vbook1(45,'dN/dp pion+  1.2D0-1.4D0',npmx,pmin,pmax)
      call vbook1(46,'dN/dp pion+  1.4D0-1.6D0',npmx,pmin,pmax)
      call vbook1(47,'dN/dp pion+  1.6D0-1.8D0',npmx,pmin,pmax)
      call vbook1(48,'dN/dp pion+  1.8D0-2.0D0',npmx,pmin,pmax)
      call vbook1(49,'dN/dp pion+  2.0D0-2.2D0',npmx,pmin,pmax)
      call vbook1(50,'dN/dp pion+  2.2D0-2.4D0',npmx,pmin,pmax)
      call vbook1(51,'dN/dp pion+  2.4D0-2.6D0',npmx,pmin,pmax)
      call vbook1(52,'dN/dp pion+  2.6D0-2.8D0',npmx,pmin,pmax)

      pmin=0.0D0
      pmax=1.3D0
      npmx=nint((pmax-pmin)/wp)
      call vbook1(53,'dN/dp kaon-  0.8D0-1.0D0',npmx,pmin,pmax)
      call vbook1(54,'dN/dp kaon-  1.0D0-1.2D0',npmx,pmin,pmax)
      call vbook1(55,'dN/dp kaon-  1.2D0-1.4D0',npmx,pmin,pmax)
      call vbook1(56,'dN/dp kaon-  1.4D0-1.6D0',npmx,pmin,pmax)
      call vbook1(57,'dN/dp kaon-  1.6D0-1.8D0',npmx,pmin,pmax)
      call vbook1(58,'dN/dp kaon-  1.8D0-2.0D0',npmx,pmin,pmax)
      call vbook1(59,'dN/dp kaon-  2.0D0-2.2D0',npmx,pmin,pmax)

      call vbook1(60,'dN/dp kaon+  0.6D0-0.8D0',npmx,pmin,pmax)
      call vbook1(61,'dN/dp kaon+  0.8D0-1.0D0',npmx,pmin,pmax)
      call vbook1(62,'dN/dp kaon+  1.0D0-1.2D0',npmx,pmin,pmax)
      call vbook1(63,'dN/dp kaon+  1.2D0-1.4D0',npmx,pmin,pmax)
      call vbook1(64,'dN/dp kaon+  1.4D0-1.6D0',npmx,pmin,pmax)
      call vbook1(65,'dN/dp kaon+  1.6D0-1.8D0',npmx,pmin,pmax)
      call vbook1(66,'dN/dp kaon+  1.8D0-2.0D0',npmx,pmin,pmax)
      call vbook1(67,'dN/dp kaon+  2.0D0-2.2D0',npmx,pmin,pmax)

c...Eevent weight
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
           rap=0.5D0*log( max(p(4,i)+p(3,i),1.D-8)
     $             /max(p(4,i)-p(3,i),1.D-8) )
          if(mstc(4).eq.0) then
          else if(mstc(4).eq.3) then
          else
           rap=rap+ylab
          endif

        ptsq=p(1,i)**2+p(2,i)**2
        pt=sqrt(ptsq)
        pp=sqrt(ptsq+p(3,i)**2)
        pt=max(pt,1.D-8)
        eta=0.5D0*log( max(pp+p(3,i),1.D-8)/max(pp-p(3,i),1.D-8) )
	et=p(4,i)*pt/max(pp,1.D-8)
        emt=sqrt(p(5,i)**2+ptsq)
        emt0=emt-p(5,i)
        kch=jamchge(k(2,i))

        if(eta.ge.1.25D0.and.eta.le.2.5D0) then
          call vfill1(1,et,wevt1)
          if(kch.ne.0) call vfill1(2,et,wevt1)
        endif

        if(kf.eq.2212) then
          call vfill1(11,rap,wevt/wy)
        else if(kf.eq.-211) then
          call vfill1(12,rap,wevt/wy)
        else if(kf.eq.211) then
          call vfill1(13,rap,wevt/wy)
        else if(kf.eq.-321) then
          call vfill1(14,rap,wevt/wy)
        else if(kf.eq.321) then
          call vfill1(15,rap,wevt/wy)
        endif

        if(rap.lt.0.6D0) then
        else if(rap.le.0.8D0) then
          if(kf.eq.2212) call vfill1(21,emt0,wevt/(wp*emt))
          if(kf.eq.-211) call vfill1(31,emt0,wevt/(wp*emt))
          if(kf.eq. 211) call vfill1(42,emt0,wevt/(wp*emt))
          if(kf.eq. 321) call vfill1(60,emt0,wevt/(wp*emt))
        else if(rap.le.1.0D0) then
          if(kf.eq.2212) call vfill1(22,emt0,wevt/(wp*emt))
          if(kf.eq.-211) call vfill1(32,emt0,wevt/(wp*emt))
          if(kf.eq. 211) call vfill1(43,emt0,wevt/(wp*emt))
          if(kf.eq.-321) call vfill1(53,emt0,wevt/(wp*emt))
          if(kf.eq. 321) call vfill1(61,emt0,wevt/(wp*emt))
        else if(rap.le.1.2D0) then
          if(kf.eq.2212) call vfill1(23,emt0,wevt/(wp*emt))
          if(kf.eq.-211) call vfill1(33,emt0,wevt/(wp*emt))
          if(kf.eq. 211) call vfill1(44,emt0,wevt/(wp*emt))
          if(kf.eq.-321) call vfill1(54,emt0,wevt/(wp*emt))
          if(kf.eq. 321) call vfill1(62,emt0,wevt/(wp*emt))
        else if(rap.le.1.4D0) then
          if(kf.eq.2212) call vfill1(24,emt0,wevt/(wp*emt))
          if(kf.eq.-211) call vfill1(34,emt0,wevt/(wp*emt))
          if(kf.eq. 211) call vfill1(45,emt0,wevt/(wp*emt))
          if(kf.eq.-321) call vfill1(55,emt0,wevt/(wp*emt))
          if(kf.eq. 321) call vfill1(63,emt0,wevt/(wp*emt))
        else if(rap.le.1.6D0) then
          if(kf.eq.2212) call vfill1(25,emt0,wevt/(wp*emt))
          if(kf.eq.-211) call vfill1(35,emt0,wevt/(wp*emt))
          if(kf.eq. 211) call vfill1(46,emt0,wevt/(wp*emt))
          if(kf.eq.-321) call vfill1(56,emt0,wevt/(wp*emt))
          if(kf.eq. 321) call vfill1(64,emt0,wevt/(wp*emt))
        else if(rap.le.1.8D0) then
          if(kf.eq.2212) call vfill1(26,emt0,wevt/(wp*emt))
          if(kf.eq.-211) call vfill1(36,emt0,wevt/(wp*emt))
          if(kf.eq. 211) call vfill1(47,emt0,wevt/(wp*emt))
          if(kf.eq.-321) call vfill1(57,emt0,wevt/(wp*emt))
          if(kf.eq. 321) call vfill1(65,emt0,wevt/(wp*emt))
        else if(rap.le.2.0D0) then
          if(kf.eq.2212) call vfill1(27,emt0,wevt/(wp*emt))
          if(kf.eq.-211) call vfill1(37,emt0,wevt/(wp*emt))
          if(kf.eq. 211) call vfill1(48,emt0,wevt/(wp*emt))
          if(kf.eq.-321) call vfill1(58,emt0,wevt/(wp*emt))
          if(kf.eq. 321) call vfill1(66,emt0,wevt/(wp*emt))
        else if(rap.le.2.2D0) then
          if(kf.eq.2212) call vfill1(28,emt0,wevt/(wp*emt))
          if(kf.eq.-211) call vfill1(38,emt0,wevt/(wp*emt))
          if(kf.eq. 211) call vfill1(49,emt0,wevt/(wp*emt))
          if(kf.eq.-321) call vfill1(59,emt0,wevt/(wp*emt))
          if(kf.eq. 321) call vfill1(67,emt0,wevt/(wp*emt))
        else if(rap.le.2.4D0) then
          if(kf.eq.2212) call vfill1(29,emt0,wevt/(wp*emt))
          if(kf.eq.-211) call vfill1(39,emt0,wevt/(wp*emt))
          if(kf.eq. 211) call vfill1(50,emt0,wevt/(wp*emt))
        else if(rap.le.2.6D0) then
          if(kf.eq.-211) call vfill1(40,emt0,wevt/(wp*emt))
          if(kf.eq. 211) call vfill1(51,emt0,wevt/(wp*emt))
        else if(rap.le.2.8D0) then
          if(kf.eq.-211) call vfill1(41,emt0,wevt/(wp*emt))
          if(kf.eq. 211) call vfill1(52,emt0,wevt/(wp*emt))
        endif

        if(kf.eq.2212.or.kf.eq.2112) goto 3000
        if(eta.ge.1.25D0.and.eta.le.2.5D0) then
          call vfill1(1,et,wevt1)
          if(kch.ne.0) call vfill1(2,et,wevt1)
        endif
        call vfill1(3,eta,et*wevt2)
        if(kch.ne.0) call vfill1(4,eta,wevt2)

3000     end do

         return

c*************************************************************************

      entry anal3

c...Output hostogram

      fac=paru(1)*10*(parc(4)**2-parc(3)**2)
      do j=1,2
      call vscale(j,fac)
      call vprint(j,0,1)
      end do
      do j=3,4
      call vscale(j,1.0D0)
      call vprint(j,0,1)
      end do

c...Rapidity
      do j=1,5
      call vscale(10+j,fac)
      call vprint(10+j,0,1)
      end do

c...Mt distributions.
      fac= paru(1)*10*(parc(4)**2-parc(3)**2)/(2*paru(1)*0.2D0)
      do i=1,9
      call vscale(20+i,fac)
      call vprint(20+i,0,1)
      end do
      do i=31,67
      call vscale(i,fac)
      call vprint(i,0,1)
      end do

      end
