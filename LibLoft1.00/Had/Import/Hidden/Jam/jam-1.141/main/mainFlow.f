c...A main program for calculation of flow at SPS for WA98.

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
c     fname(1)='jam.cfg'  ! input file name.
      mstc(1) =48827   ! random seed.
      mevent=600        ! total simulation event
c     bmin=8.0D0          ! minimum impact parameter
c     bmax=-10.0D0        ! maximum impact parameter
      bmin=4.0D0          ! minimum impact parameter
      bmax=-8.0D0        ! maximum impact parameter
      dt=100.0D0          ! collision time(fm/c)
      nstep=1
      cwin='40.0gev         '  ! incident energy
c     cwin='158gev         '  ! incident energy
      frame='nn      '        ! comp. frame
      proj='197Au   '         ! projectile
      targ='197Au   '         ! projectile
c     proj='208Pb   '         ! projectile
c     targ='208Pb   '         ! target
      mstc(8)=0   ! job mode.
      mstc(156)=0  ! analysis of collision distribution
      mstc(155)=0  ! flow anal.
      mstc(163)=0  ! time evolution of directed transverse flow
      parc(7)= 1.0D0    ! Output time interval (fm/c)

c     mstc(51)=0  ! BB collisions only
c     mstc(6)=-1  ! Glauber

      call jaminit(mevent,bmin,bmax,dt,nstep,
     $                             frame,proj,targ,cwin)
      nevent=mstc(2)
      if(dump)write(33)nevent,pard(17),pard(5),pard(6),mstc(4)

c...Initialize analysis.
      call anal1

c...Simulation start.
      do iev=1,nevent

c...Simulate one event.
        call jamevt(iev)

c...Dump phase space data.
        if(dump) then
          write(33)iev,nv,nbary,nmeson,pard(2)
          do i=1,nv
            write(33)(k(j,i),j=1,11),(r(j,i),j=1,5),(p(j,i),j=1,5)
          end do
        endif

        if(mod(iev,100).eq.0) write(6,*)'event=',iev

c...Data analysis.
        call anal2

c...List phase space data.
c       call jamlist(1)

      end do

      if(dump) close(33)

c...Final output.
      call jamfin

c...Print analysis results.
      call anal3

      end

c***********************************************************************

      subroutine anal1

      include 'jam1.inc'
      include 'jam2.inc'
      save wy,wp
      save ylab,wevt

      wevt=1.D0/dble(mstc(2))

c....Rapidity distribution.
      ylab=pard(17)
c     ymin=-7.0D0
c     ymax=7.0D0
c     wy=0.25D0
c     nymx=(ymax-ymin)/wy

      ymin=-ylab*1.3
      ymax=ylab*1.3
      nymx=25
      wy=(ymax-ymin)/nymx

      pmin=0.0D0
      pmax=5.0D0
      wp=0.1D0
      npmx=(pmax-pmin)/wp

c....Directed and ellipse flow as a function of rapidity.
      call vbook1(11,'v1(y) nucleons',nymx,ymin,ymax)
      call vbook1(12,'v1(y) pions   ',nymx,ymin,ymax)
      call vbook1(13,'v2(y) nucleons',nymx,ymin,ymax)
      call vbook1(14,'v2(y) pions   ',nymx,ymin,ymax)

c....Directed and ellipse flow as a function of pt.
      call vbook1(21,'v2(pt) nucleons',npmx,pmin,pmax)
      call vbook1(22,'v2(pt) pions   ',npmx,pmin,pmax)
      call vbook1(23,'v2(pt) nucleons',npmx,pmin,pmax)
      call vbook1(24,'v2(pt) pions   ',npmx,pmin,pmax)

c....Directed and ellipse flow as a function of rapidity.
      call vbook1(31,'<v1(y)> nucleons',nymx,ymin,ymax)
      call vbook1(32,'<v1(y)> pions   ',nymx,ymin,ymax)
      call vbook1(33,'<v2(y)> nucleons',nymx,ymin,ymax)
      call vbook1(34,'<v2(y)> pions   ',nymx,ymin,ymax)

c....Directed and ellipse flow as a function of pt.
      call vbook1(41,'<v1(pt)> nucleons',npmx,pmin,pmax)
      call vbook1(42,'<v1(pt)> pions   ',npmx,pmin,pmax)
      call vbook1(43,'<v2(pt)> nucleons',npmx,pmin,pmax)
      call vbook1(44,'<v2(pt)> pions   ',npmx,pmin,pmax)

      call vbook1(61,'n(y) nucleons   ',nymx,ymin,ymax)
      call vbook1(62,'n(y) pions      ',nymx,ymin,ymax)
      call vbook1(63,'n(y) protons    ',nymx,ymin,ymax)
      call vbook1(64,'n(y) pion+      ',nymx,ymin,ymax)
      call vbook1(65,'n(pt) nucleons  ',npmx,pmin,pmax)
      call vbook1(66,'n(pt) pions     ',npmx,pmin,pmax)

c....Transverse flow.
      call vbook1(51,'px(y) nucleons',nymx,ymin,ymax)
      call vbook1(52,'pt(y) nucleons',nymx,ymin,ymax)

      call vbook1(53,'px(y) pions   ',nymx,ymin,ymax)
      call vbook1(54,'pt(y) pions   ',nymx,ymin,ymax)

      call vbook1(55,'px(y) protons ',nymx,ymin,ymax)
      call vbook1(56,'pt(y) protons ',nymx,ymin,ymax)

      call vbook1(57,'px(y) pion+   ',nymx,ymin,ymax)
      call vbook1(58,'pt(y) pion+   ',nymx,ymin,ymax)


      call vbook1(71,'<px(y)> nucleons',nymx,ymin,ymax)
      call vbook1(72,'<pt(y)> nucleons',nymx,ymin,ymax)

      call vbook1(73,'<px(y)> pions   ',nymx,ymin,ymax)
      call vbook1(74,'<pt(y)> pions   ',nymx,ymin,ymax)

      call vbook1(75,'<px(y)> protons ',nymx,ymin,ymax)
      call vbook1(76,'<pt(y)> protons ',nymx,ymin,ymax)

      call vbook1(77,'<px(y)> pion+   ',nymx,ymin,ymax)
      call vbook1(78,'<pt(y)> pion+   ',nymx,ymin,ymax)

      return

c***********************************************************************

      entry anal2

c...Loop over all particles.
      do 3000 i=1,nv

c...Exclude spectetor.
        if(abs(k(7,i)).eq.1) goto 3000

        kf=k(2,i)
        kc=jamcomp(kf)
        if(kc.le.0.or.kc.gt.mstu(6)) then
           write(6,*)'Invalid code i kf kc',i,kf,kc,nv,nbary,nmeson
           goto 3000
        endif

        y=0.5D0*log( max(p(4,i)+p(3,i),1.D-8)/max(p(4,i)-p(3,i),1.D-8) )

c       if(mstc(4).eq.0) then
c       else if(mstc(4).eq.3) then
c       else
c        y=y+ylab
c       endif

        ptsq=p(1,i)**2+p(2,i)**2
        pt=sqrt(ptsq)
        pp=sqrt(ptsq+p(3,i)**2)
        pt=max(pt,1.D-8)
        eta=0.5D0*log( max(pp+p(3,i),1.D-8)/max(pp-p(3,i),1.D-8) )
        et=p(4,i)*pt/max(pp,1.D-8)
        px=p(1,i)
        if(pt.gt.1D-8) then
          cos1=px/pt
          cos2 = (p(1,i)**2-p(2,i)**2)/ptsq
        else
          cos1=0.0D0
          cos2=0.0d0
        endif
c       phi=acos(cos1)
c       cos2=cos(2*phi)

c......Transverse flow of nucleons and pions.
        if(kf.eq.2212.or.kf.eq.2112) then
          call vfill1(61,y,1.0D0)
          call vfill1(65,pt,1.0D0)
          call vfill1(51,y,px)
          call vfill1(52,y,pt)
          call vfill1(11,y,cos1)
          call vfill1(13,y,cos2)
          call vfill1(21,pt,cos1)
          call vfill1(23,pt,cos2)
c....Protons.
          if(kf.eq.2212) then
            call vfill1(63,abs(y),0.5D0)
            call vfill1(55,abs(y),abs(px))
            call vfill1(56,abs(y),pt)
          endif
        else if(kf.eq.111.or.abs(kf).eq.211) then
          call vfill1(62,y,1.0D0)
          call vfill1(66,pt,1.0D0)
          call vfill1(53,y,px)
          call vfill1(54,y,pt)
          call vfill1(12,y,cos1)
          call vfill1(14,y,cos2)
          call vfill1(22,pt,cos1)
          call vfill1(24,pt,cos2)
c.....Positive pions.
          if(kf.eq.211) then
            call vfill1(64,abs(y),0.5D0)
            call vfill1(57,abs(y),abs(px))
            call vfill1(58,abs(y),pt)
          endif
        endif

3000  continue

      return

c***********************************************************************

      entry anal3

c...Output histograms.

      fac=1.0D0
      mnorm=0
      mform=1

c....Directed and ellipse flow.
      call vopera(11,'/',61,31,1.0D0,1.0D0)
      call vopera(12,'/',62,32,1.0D0,1.0D0)
      call vopera(13,'/',61,33,1.0D0,1.0D0)
      call vopera(14,'/',62,34,1.0D0,1.0D0)

      call vopera(21,'/',65,41,1.0D0,1.0D0)
      call vopera(22,'/',66,42,1.0D0,1.0D0)
      call vopera(23,'/',65,43,1.0D0,1.0D0)
      call vopera(24,'/',66,44,1.0D0,1.0D0)

c....Transverse flow.
      call vopera(51,'/',61,71,1.0D0,1.0D0)
      call vopera(52,'/',61,72,1.0D0,1.0D0)

      call vopera(53,'/',62,73,1.0D0,1.0D0)
      call vopera(54,'/',62,74,1.0D0,1.0D0)

      call vopera(55,'/',63,75,1.0D0,1.0D0)
      call vopera(56,'/',63,76,1.0D0,1.0D0)

      call vopera(57,'/',64,77,1.0D0,1.0D0)
      call vopera(58,'/',64,78,1.0D0,1.0D0)


      call vprint(61,mnorm,mform)
      call vprint(62,mnorm,mform)
      do i=1,4
        call vprint(30+i,mnorm,mform)
        call vprint(40+i,mnorm,mform)
      end do
      do i=1,8
        call vprint(70+i,mnorm,mform)
      end do

      end

