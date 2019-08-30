c...Main program for p+p~  UA1 srt=200GeV
c...P.R. D36(1987)16  NA24 Collaboration

      include 'jam1.inc'
      include 'jam2.inc'
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      character frame*8,proj*8,targ*8,cwin*15
      logical dump
      data dump/.false./
c     data dump/.true./


c     mstc(1)=143   ! random seed.
c     em=0.91
c     kf1=1
c     kf2=2101
c     call kfcnst(kf1,kf2,kf0,em)
c     kf0=10211
c     kf=kf0
c     em=0.349
c     call jamidres(kf,em,icon)
c     print *,'kf kf0',icon,kf0,kf,em

c     kfl10=4
c     kfl20=2101
c     print *,'em',emjetb(kfl10,kfl20)

c     kf=-10413
c     call jamdmass(kf,kfm,kfd,emin,emdn)
c     print *,'kf kfm kfd emin emdn',kf,kfm,kfd,emin,emdn
c     stop

      if(dump)
     $  open(33,file='phase.dat',form='unformatted',status='new')

c...Initialization
c     mstp(92)=3
c     parp(94)=1.5d0
c     parp(96)=1.5d0

c     mstp(92)=5
c     parp(98)=0.5d0
c     parp(96)=1.5d0

c     parp(94)=0.15d0*2
c     parp(96)=0.15d0*2

c     parp(95)=0.5d0
c     parp(97)=0.5d0

      mstc(1)=1418823   ! random seed.
      mevent=10000       ! total simulation event
      bmin=0.0D0          ! minimum impact parameter
      bmax=-1.0D0         ! maximum impact parameter
      dt=100.0D0          ! collision time(fm/c)
      nstep=1
      cwin='200gev         '  ! incident energy
      frame='collider'        ! comp. frame
      proj='p       '         ! projectile
      targ='p       '         ! target
      mstc(8)=0   ! job mode.
c     mstc(71)=0   ! 1:dffractive 
      mstc(72)=1   ! switch for hadron formation point from string decay.

c     mstc(81)=1   ! 1:hard scattering off
c     mstc(74)=1   ! soft radiaton
c     parc(56)=4.0

c     parp(91)=0.8 ! gaussian width of primordial kt dist.

      call jaminit(mevent,bmin,bmax,dt,nstep,
     &               frame,proj,targ,cwin)
      nevent=mstc(2)
      if(dump)write(33)nevent,pard(17),pard(5),pard(6),mstc(4)

c...Initialize analysis.
      call anal1

      do iev=1,nevent

        call jamevt(iev)
        if(mod(iev,500).eq.0)write(6,*)' event= ',iev
        call anal2
c...Dump phase space data.
        if(dump) then
          write(33)iev,nv,pard(2)
          do i=1,nv
            write(33)(k(j,i),j=1,9),(r(j,i),j=1,5),(p(j,i),j=1,5)
     $              ,(v(j,i),j=1,5)
          end do
        endif

      end do

      if(dump) close(33)
      call jamfin
      call anal3

      end

c***********************************************************************

      subroutine anal1

      include 'jam1.inc'
      include 'jam2.inc'
      save wp,wevt,fac,ycut1,ycut2,wy
     
      ycut1=0.25D0
      ycut2=-0.25D0
      dy=ycut1-ycut2

C...Initialize histogram booking.
C.....rapidity spectra.
      ymin=-7.D0
      ymax=+7.D0
      wy=0.25D0
      nymx=(ymax-ymin)/wy
      call vbook1(1,'dN/deta - all',nymx,ymin,ymax)
      call vbook1(2,'dN/deta - charged',nymx,ymin,ymax)


c...Momentum cut.
      pmax=10.0D0
      pmin=0.0D0
      wp=0.15D0
      npmx=nint((pmax-pmin)/wp)
c     npmx=50 
c     wp=(pmax-pmin)/float(npmx)
      call vbook1(11,'Eds/dp^3 pi0',npmx,pmin,pmax)
      call vbook1(12,'Eds/dp^3 pi+',npmx,pmin,pmax)
      call vbook1(13,'Eds/dp^3 pi-',npmx,pmin,pmax)
      call vbook1(14,'Eds/dp^3 all',npmx,pmin,pmax)
      call vbook1(15,'Eds/dp^3 charged',npmx,pmin,pmax)

c...Eevent weight
      fac=0.5*paru(1)*10*(parc(4)**2-parc(3)**2)/(2*paru(1)*dy)
      wevt=1.0D0/dble(mstc(2))

      call vbook1(101,'dN/dy - charged',nymx,ymin,ymax)
      call vbook1(102,'dE_t/dy',nymx,ymin,ymax)
      call vbook1(103,'dN/dy - net baryon',nymx,ymin,ymax)
      call vbook1(104,'dN/2piptdy - pion -0.5<y<0.5',npmx,pmin,pmax)

      call vbook1(111,'dN/dy - baryon',nymx,ymin,ymax)
      call vbook1(112,'dN/dy - anti-baryon',nymx,ymin,ymax)


      return

c***********************************************************************

        entry anal2

c...Loop over all particles
        do i=1,nv
         
         kf=k(2,i)
         kfa=abs(kf)
         kc=jamcomp(kf)
         if(kc.le.0.or.kc.gt.mstu(6))then
           write(6,*)'Invalide code at i, kf, kc : ',i,kf,kc
           go to 3000
         end if

c...Rapidity cut.
        yp=0.5*log( max(p(4,i)+p(3,i),1.d-8)/max(p(4,i)-p(3,i),1.d-8) )
        pt=sqrt(p(1,i)**2+p(2,i)**2)
        pp=sqrt(pt**2+p(3,i)**2)
        yeta=0.5D0*log( max(pp+p(3,i),1.D-8)/max(pp-p(3,i),1.D-8) )
        et=sqrt(max(0d0,p(5,i)**2+pt**2))

c       if(mstc(4).eq.0) then
c       else if(mstc(4).eq.3) then
c       else
c        yp=yp+ylab
c       endif

c...Et.
        call vfill1(102,yp,et*wevt/wy)

        kch=jamchge(k(2,i))
        call vfill1(1,yeta,wevt/wy)

c...Charged particles.
        if(kch.ne.0) then
          call vfill1(2,yeta,wevt/wy)
          call vfill1(101,yp,wevt/wy)
        endif

c...Net baryons.
        ibar=kchg(kc,6)
        if(ibar.eq.3) then
          call vfill1(103,yp,wevt*isign(1,kf)/wy)
          if(kf.gt.0) then
            call vfill1(111,yp,wevt/wy) 
          else if(kf.lt.0) then
            call vfill1(112,yp,wevt/wy) 
          endif
        endif

c...Transverse momentum at midrapidiy.
      if(abs(yp).le.0.5d0) then
         if(kf.eq.111.or.kfa.eq.211) then
           call vfill1(104,pt,wevt/(paru(2)*pt*wp))
         endif
      endif



c....Rapidity cut.
      if(yeta.lt.ycut2.or.yeta.gt.ycut1) goto 3000

       call vfill1(14,pt,wevt/(pt*wp))
       if(kch.ne.0) call vfill1(15,pt,wevt/(pt*wp))

c...pi0
         if(kf.eq.111) then
           call vfill1(11,pt,wevt/(pt*wp))
         else if(kf.eq.211) then
           call vfill1(12,pt,wevt/(pt*wp))
         else if(kf.eq.-211) then
           call vfill1(13,pt,wevt/(pt*wp))
         end if

3000     end do

         return

c***********************************************************************

        entry anal3

c...Output hostogram

      do j=1,2
        call vprint(j,0,1)
      end do

      do j=1,5
        call vscale(10+j,fac)
        call vprint(10+j,0,1)
      end do

      do i=1,4
        call vprint(100+i,0,1)
      end do
      call vprint(111,0,1)
      call vprint(112,0,1)

        end
