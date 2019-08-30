c...A example main program to illustrate how to use the option "user".
c...In that case, common blocks /jamevnt1/,/jamevnt2//jamevnt3/
c...should be specified by the user.
c...This program simulate K+ 40Ca. reaction. K+ is sampled inside 40Ca.

      include '../jam1.inc'
      include '../jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15

c....Initial parameters of JAM.
      mstc(1) =18827      ! random seed.
      mevent=3            ! total simulation event
      frame='user'        ! comp. frame
      dt=100.0d0          ! collision time(fm/c) =nstep*dt
      nstep=1             ! time step (shuld be =1 for cascade simulation)
      mstc(8)=0           ! job mode.
      mstc(16)=0          ! display on/off.
      parc(6)= 1.0        !  scale of display

      kfp=321             ! proj. PDG code.
c     kfp=-321            ! proj. PDG code.
      pproj=1.0           ! incident momentum.
      in=40               ! targ. neutron number.
      iz=20               ! targ. proton number.
      cwin='0.5gevc'      ! 2 body c.m energy expected in case of "user".
      proj='K+'           ! projectile
      targ='40Ca'         ! target

c...Dummy in the case of frame='user'.==============================
      bmin=0.0d0          ! minimum impact parameter
      bmax=-1.0d0         ! maximum impact parameter
c...End dummy ======================================================

c...Initialize JAM.
      call jaminit(mevent,bmin,bmax,dt,nstep,frame,proj,targ,cwin)
      nevent=mstc(2)

c...Simulation start.
      do iev=1,nevent

c...Sampling particle momentum and coordinate.
        call init_mine(kfp,pproj,in,iz)

c...Simulate one event.
        call jamevt(iev)

        if(mod(iev,100).eq.0) write(6,*)'event=',iev

c...List phase space data.
c       call jamlist(1)

      end do

c...Final output.
c      call jamfin


      end

c***********************************************************************

      subroutine init_mine(kf,pproj,in,iz)

c...A example to illustrate how to use the option "user".
c...Target is a nuleus (in,iz) and proj. is sampled randomly
c...inside a target.
c....kf: PDG code for proj.
c....pproj: incident momentum for proj.
c....in: neutron number of targ.
c....iz: proton number of targ.

      include '../jam1.inc'
      include '../jam2.inc'
c...My local common block for Wood-Saxon parameters.
      common /myjam1/rad0,rt00,radm,rwmax

      nv=in+iz+1   ! number of total initial particle.
      nmeson=1     ! number of mesons
      nbary=in+iz  ! number of baryons

c....Coordinate of targ.
      call woodsax(nbary,in)

c....Momentum of targ.
      call fermi_mom(nbary)

c...Set Proj.
      call jamzero(nv)
      p(1,nv)=0d0 
      p(2,nv)=0d0 
      p(3,nv)=pproj
      p(5,nv)=pjmass(kf)
      p(4,nv)=sqrt(p(5,nv)**2+p(3,nv)**2)
      k(1,nv)=1
      k(2,nv)=kf
      k(5,nv)=-1
      k(7,nv)=-1
      k(8,nv)=1
      kc=jamcomp(kf)
      k(9,nv)=kchg(kc,6)*isign(1,kf)  ! baryon number

c...Random distribution inside nucleus.
  10  rr = radm*rn(0)**1d0/3d0
      if(rn(0)*rwmax.gt.1d0/(1d0+exp((rr-rt00)/parc(11)))) goto 10
      cx=1.d0-2.d0*rn(0)
      sx=sqrt(1.d0-cx**2)
      phi=2*paru(1)*rn(0)
      r(1,nv)=rr*sx*cos(phi)
      r(2,nv)=rr*sx*sin(phi)
      r(3,nv)=rr*cx
      r(4,nv)=0.0
      r(5,nv)=0.0
      v(5,nv)=1d+35  ! life time.
      do j=1,4
       v(j,nv)=r(j,nv)
      end do

      end

c***********************************************************************
      subroutine woodsax(nnn,in)

      include '../jam1.inc'
      include '../jam2.inc'
c...My local common block for Wood-Saxon parameters.
      common /myjam1/rad0,rt00,radm,rwmax

c...Parameters for Wood Saxon sampling.
      rad0=1.19d0*nbary**0.333333d0-1.61d0*nbary**(-0.333333d0)
      rt00=rad0
      radm=rt00+2.5d0
      rwmax=1.0d0/(1.0d0+exp(-rt00/parc(11)))

c...Loop over target.
      do i=1,nnn

c...Zero the vector.
        call jamzero(i)

c....Sample coordinate accroding to the Wood-Saxon shape.
  10    rr = radm*rn(0)**1d0/3d0
        if(rn(0)*rwmax.gt.1d0/(1d0+exp((rr-rt00)/parc(11)))) goto 10
         cx=1.d0-2.d0*rn(0)
         sx=sqrt(1.d0-cx**2)
         phi=2*paru(1)*rn(0)
         r(1,i)=rr*sx*cos(phi)
         r(2,i)=rr*sx*sin(phi)
         r(3,i)=rr*cx
         r(4,i)=0.0
         r(5,i)=0.0
         v(5,i)=1d+35  ! life time.
    
c....Particle ID.
        if(i.le.in) then
          kf=2112   ! neutron
        else
          kf=2212   ! proton
        endif
        k(1,i)=1
        k(2,i)=kf
        k(5,i)=-1
        k(7,i)=1
        k(8,i)=1
        k(9,i)=3  ! baryon number
        p(5,i)=pjmass(kf)
      end do

c...C.M..correction.
      cx=0.d0
      cy=0.d0
      cz=0.d0
      s=0.d0
      do 100 i=1,nnn
        cx=cx+r(1,i)*p(5,i)
        cy=cy+r(2,i)*p(5,i)
        cz=cz+r(3,i)*p(5,i)
        s=s+p(5,i)
 100  continue
      cx=-cx/s
      cy=-cy/s
      cz=-cz/s
      do 101 i=1,nnn
        r(1,i)=r(1,i)+cx
        r(2,i)=r(2,i)+cy
        r(3,i)=r(3,i)+cz
        v(1,i)=r(1,i)
        v(2,i)=r(2,i)
        v(3,i)=r(3,i)
 101  continue

      end

c***********************************************************************
      subroutine fermi_mom(nnn)

c...Sample Fermi momentum inside nucleus.
      include '../jam1.inc'
      include '../jam2.inc'
      logical first
      data first /.true./

      if(first) then
        call hijwds(nnn,1,xmax)
        first=.false.
      endif

c....Sampling Fremi momenta of nucleus.
      do i=1,nnn

c...Calculate local Ferim momentum.
        dens=wdsax1(sqrt(r(1,i)**2+r(2,i)**2+r(3,i)**2))*nnn
        pfermi=paru(3)*(1.5d0*paru(1)**2*dens)**(1.0d0/3.0d0)
        pf=pfermi*rn(0)**(1.d0/3.d0)
        cth=1.d0-2.d0*rn(0)
        sth=sqrt(1.d0-cth**2)
        phi=paru(2)*rn(0)
        p(3,i)=pf*cth
        p(1,i)=pf*sth*cos(phi)
        p(2,i)=pf*sth*sin(phi)
      end do

c...C.M..correction.
      cx=0.d0
      cy=0.d0
      cz=0.d0
      do i=1,nnn
        cx=cx+p(1,i)
        cy=cy+p(2,i)
        cz=cz+p(3,i)
      end do

      cx=-cx/nnn
      cy=-cy/nnn
      cz=-cz/nnn
      do i=1,nnn
        p(1,i)=p(1,i)+cx
        p(2,i)=p(2,i)+cy
        p(3,i)=p(3,i)+cz
        p(4,i)=sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
      end do

      end

