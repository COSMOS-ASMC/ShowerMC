c...Main program for p(300GeV/c)+p 
c...P.R. D36(1987)16  NA24 Collaboration

        include 'jam1.inc'
        include 'jam2.inc'
        character frame*8,proj*8,targ*8,cwin*15

c...Initialization
      mstc(1) =48827   ! random seed.
      mevent=10000      ! total simulation event
      bmin=0.0D0          ! minimum impact parameter
      bmax=-1.0D0         ! maximum impact parameter
      dt=100.0D0          ! collision time(fm/c)
      nstep=1
      cwin='300gevc        '  ! incident energy
      frame='nn      '        ! comp. frame
      proj='p       '         ! projectile
      targ='p       '         ! target
      mstc(8)=0   ! job mode.
      mstc(81)=1   ! 1:hard scattering off
      mstc(17)=1       ! only inelastic collisions.

      call jaminit(mevent,bmin,bmax,dt,nstep,
     &               frame,proj,targ,cwin)
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
      save wp,wevt,fac,ymax,ymin
     
      ymax=0.52D0
      ymin=-0.65D0
      dy=ymax-ymin

c...Momentum cut.
      pmax=9.0D0
      pmin=0.0D0
      wp=0.5D0
      npmx=nint((pmax-pmin)/wp)
c     npmx=50 
c     wp=(pmax-pmin)/float(npmx)
      call vbook1(1,'Eds/dp^3 pi0',npmx,pmin,pmax)
      call vbook1(2,'Eds/dp^3 pi+',npmx,pmin,pmax)
      call vbook1(3,'Eds/dp^3 pi-',npmx,pmin,pmax)

c...Eevent weight
      fac= paru(1)*10*(parc(4)**2 - parc(3)**2)
      wevt=1.0D0/dble(mstc(2))/(2*paru(1)*wp*dy)

      return

c*************************************************************************

        entry anal2

c...Loop over all particles
        do i=1,nv
         
         kf=k(2,i)
         kc=jamcomp(kf)
         if(kc.le.0.or.kc.gt.mstc(6))then
           write(6,*)'Invalide code at i, kf, kc : ',i,kf,kc
           go to 3000
         end if

c...Rapidity cut.
           rap=0.5D0*log( max(p(4,i)+p(3,i),1.D-8)
     $             /max(p(4,i)-p(3,i),1.D-8) )
c         if(mstc(4).eq.0) then
c         else if(mstc(4).eq.3) then
c         else
c          rap=rap+ylab
c         endif
          if(rap.lt.ymin.or.rap.gt.ymax) goto 3000

          pt=sqrt(p(1,i)**2+p(2,i)**2)
c...pi0
         if(kf.eq.111) then
           call vfill1(1,pt,wevt/pt)
         else if(kf.eq.211) then
           call vfill1(2,pt,wevt/pt)
         else if(kf.eq.-211) then
           call vfill1(3,pt,wevt/pt)
         end if

3000     end do

         return

c*************************************************************************

        entry anal3

c...Output hostogram

        do j=1,3
        call vscale(j,fac)
        call vprint(j,0,1)
        end do

        end
