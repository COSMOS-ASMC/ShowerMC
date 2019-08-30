C***********************************************************************
 
      program main400gev
 
C...Purpose: Generate fixed target p+p collisions at momentum P_beam = 400 GeV.

      include 'jam1.inc'
      include 'jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15

c....Initialize JAM.
      mstc(1)=1712915     ! random seed.
      mevent=10000        ! total simulation event
      bmin=0.0d0          ! minimum impact parameter in fm
      bmax=-1.08d0        ! maximum impact parameter in fm

      dt=100.0d0          ! collision time(fm/c)
      nstep=1
      cwin='400gevc'      ! incident energy
      frame='nn'          ! comp. frame
      proj='p'            ! projectile
      targ='p'            ! target

c...Options.
      mstc(8)=0           ! job mode.
      mstc(17)=1          ! only inelastic collisions.
c     mstc(81)=1          ! 1:hard scattering on/off
c     mstc(6)=-1          ! Glauber
c     mstc(74)= 1         ! soft rad.

      call jaminit(mevent,bmin,bmax,dt,nstep,frame,proj,targ,cwin)

c...Initialize analysis.
      call anal1

c...Simulation start.
      do iev=1,mevent

c...Simulate one event.
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
      
 
C***********************************************************************

      subroutine anal1
C...Initial values and definitions for histograms.

      include 'jam1.inc'
      include 'jam2.inc'
      double precision wevt,fac,ycut1,ycut2
      real*8 ymin,ymax,wy
      real*8 pmin,pmax,wp
      real*8 dy
      integer nymx,npmx
      integer i,j,iii,kf,kfa,kch
      real*8 yp,pt,pp,yeta,et
      real*8 psum(5)
      save wp,wevt,fac,ycut1,ycut2,wy,mult,psum,ylab
     
      ylab=pard(17)
      ycut1=0.25D0
      ycut2=-0.25D0
      dy=ycut1-ycut2

c...Initialize particle multiplicities.
      do i=1,5
        psum(i)=0d0
      end do
 

C...Initialize histogram booking.
C.....rapidity spectra.
      ymin=-7.D0
      ymax=+7.D0
      wy=0.25D0
      nymx=(ymax-ymin)/wy
      call vbook1(1,'dN/deta - all',nymx,ymin,ymax)
      call vbook1(2,'dN/deta - charged',nymx,ymin,ymax)

      call vbook1(21,'dN/y - charged',nymx,ymin,ymax)
      call vbook1(22,'dN/y - negative',nymx,ymin,ymax)
      call vbook1(23,'dN/y - pi-',nymx,ymin,ymax)
      call vbook1(24,'dN/y - pi0',nymx,ymin,ymax)
      call vbook1(25,'dN/y - pi+',nymx,ymin,ymax)
      call vbook1(26,'dN/y - k-',nymx,ymin,ymax)
      call vbook1(27,'dN/y - k+',nymx,ymin,ymax)
      call vbook1(28,'dN/y - p',nymx,ymin,ymax)
      call vbook1(29,'dN/y - pbar',nymx,ymin,ymax)


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
c     fac=0.5*3.14159*10*(parc(4)**2-parc(3)**2)/(2*paru(1)*dy)
      fac=0.5*3.14159*10*3.0/(2*3.14159*dy)
      wevt=1.0d0/mstc(2)

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
         
         if(k(1,i).le.0.or.k(1,i).gt.10) goto 3000
         if(iabs(k(2,i)).le.100) goto 3000
         kf=k(2,i)
         kfa=abs(kf)
         kch=jamchge(kf)
         kc=jamcomp(kf)

         do j=1,4
           psum(j)=psum(j)+p(j,i)
         end do


c...Rapidity cut.
        yp=0.5*log( 
     $ max(p(4,i)+p(3,i),1.d-8)/max(p(4,i)-p(3,i),1.d-8) )
        pt=sqrt(p(1,i)**2+p(2,i)**2)
        pp=sqrt(pt**2+p(3,i)**2)
        yeta=0.5D0*log(
     $  max(pp+p(3,i),1.D-8)/max(pp-p(3,i),1.D-8) )
        et=sqrt(max(0d0,p(5,i)**2+pt**2))

cc      if(mstc(4).eq.0) then
cc      else if(mstc(4).eq.3) then
cc      else
cc       yp=yp+ylab
cc      endif

c...Et.
        call vfill1(102,yp,et*wevt/wy)
        call vfill1(1,yeta,wevt/wy)

c...Charged particles.
        if(kch.ne.0) then
          call vfill1(2,yeta,wevt/wy)
          call vfill1(21,yp,wevt/wy)
          call vfill1(101,yp,wevt/wy)
          if(kch.lt.0) call vfill1(22,yp,wevt/wy)
        endif


c...Net baryons.
c       ibar=kchg(kc,6)
c       if(ibar.eq.3) then
c         call vfill1(103,yp,wevt*isign(1,kf)/wy)
c         if(kf.gt.0) then
c           call vfill1(111,yp,wevt/wy) 
c         else if(kf.lt.0) then
c           call vfill1(112,yp,wevt/wy) 
c         endif
c       endif

      iii=0
      if(kf.eq.-211) iii=23
      if(kf.eq. 111) iii=24
      if(kf.eq. 211) iii=25
      if(kf.eq.-321) iii=26
      if(kf.eq. 321) iii=27
      if(kf.eq. 2212) iii=28
      if(kf.eq. -2212) iii=29
      if(iii.ne.0) call vfill1(iii,yp,wevt/wy)

c...Transverse momentum at midrapidiy.
      if(abs(yp).le.0.5d0) then
         if(kf.eq.111.or.kfa.eq.211) then
           call vfill1(104,pt,wevt/(2*3.14159d0*pt*wp))
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

c...Rapidity distributions.
      do j=1,9
        call vprint(20+j,0,1)
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

      write(6,*)(psum(j),j=1,4)


        end
