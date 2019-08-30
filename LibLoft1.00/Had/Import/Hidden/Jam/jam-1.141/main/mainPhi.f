c...A main program for calculation of Pb+Pb at SPS with phi meson stable.

      include 'jam1.inc'
      include 'jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15
      logical dump 
c     data dump/.false./
      data dump/.true./

c     fname(1)='jam.cfg'  ! input file name.

      if(dump)
c    $  open(33,file='phase.dat',form='unformatted',status='new')
     $  open(33,file='phase.dat',status='new')

c...Make phi meson stable.
c     mdcy(jamcomp(333),1)=1   ! allowed to decay (default)
      mdcy(jamcomp(333),1)=0   ! stable

c...Change Phi decay width in GeV
       pmas(jamcomp(333),2)=0.00883d0
  

c....Initialize JAM.
      mstc(1) =328271  ! random seed.
      mevent=1000      ! total simulation event
      bmin=0.0D0       ! minimum impact parameter
      bmax=-3.2d0      ! max. impact parameter
      dt=100.0D0       ! collision time(fm/c)
      nstep=1
      cwin='158gev'    ! incident energy
      frame='nn'       ! comp. frame
      proj='208Pb'     ! projectile
      targ='208Pb'     ! target

c....RHIC
c     cwin='200gev'    ! incident energy
c     frame='collider' ! comp. frame
c     proj='197Au'     ! projectile
c     targ='197Au'     ! target

c     mstc(8)=1    ! job mode.
c     mstc(156)=1  ! analysis of collision distribution
c     mstc(155)=0  ! flow anal.
c     mstc(163)=0  ! time evolution of directed transverse flow
c     parc(7)= 1.0D0    ! Output time interval (fm/c)

c     mstc(51)=0  ! BB collisions only
c     mstc(6)=-1  ! Glauber
c     mstc(81)= 0 ! <--- 1:Hard scattering on/off
c     mstc(156)=1        ! analysis of collision distribution
c     mstc(156)=1        ! energy distribution of collisions
c     mstc(162)=1        ! Output collision histroy
c     mstc(165)=1        ! Output time evol. of particle yield
c     mstc(166)=1        ! Output time evol. of particle density
c     mstc(167)=1        ! Output time evol. of density (Gaussian)
c     parc(7)= 0.5D0     ! Output time interval (fm/c)
c     parc(54)=1.6d0     ! string tension
c     parc(54)=2.0d0     ! string tension

      call jaminit(mevent,bmin,bmax,dt,nstep,
     $                             frame,proj,targ,cwin)
      nevent=mstc(2)
c     if(dump)write(33)nevent,pard(17),pard(5),pard(6),mstc(4)
      if(dump)write(33,'(i9,1x,3(f12.5,1x),i4)')
     $ nevent,pard(17),pard(5),pard(6),mstc(4)

c...Simulation start.
      do iev=1,nevent

c...Simulate one event.
        call jamevt(iev)

c...Dump phase space data.
        if(dump) then
          write(33,'(4(i9,1x),f12.6)')iev,nv,nbary,nmeson,pard(2)
c         write(33)iev,nv,nbary,nmeson,pard(2)
          do i=1,nv
c           write(33)(k(j,i),j=1,11),(r(j,i),j=1,5),(p(j,i),j=1,5)
c    $              ,(v(j,i),j=1,5)
            write(33,'(7(i9,1x),15(e12.4,1x))')
     $  (k(j,i),j=1,7),(r(j,i),j=1,5),(p(j,i),j=1,5),(v(j,i),j=1,5)
          end do
        endif

        if(mod(iev,100).eq.0) write(6,*)'event=',iev

      end do

      if(dump) close(33)

c...Final output.
      call jamfin

      end
