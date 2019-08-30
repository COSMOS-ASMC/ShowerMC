c...Main program for p(14.6GeV/c)+Be,27Al,63Cu,197Au
c...P.R. C45 (1992)2933 E-802 Collaboration

      include '../jam1.inc'
      include '../jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15

c...Initialization
c     mstc(1) =48827   ! random seed.
      fname(1) = '0'
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

      do iev=1,nevent
        call jamevt(iev)
        if(mod(iev,500).eq.0)write(6,*)' event= ',iev
      end do

      call jamfin

      end
