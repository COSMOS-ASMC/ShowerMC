c...A main program to use the initial condition of hadronic cascade
c...from statistical model.

      include 'jam1.inc'
      include 'jam2.inc'
      common/jamxml/mult(-500:500)
      character frame*8,proj*8,targ*8,cwin*15
      common/bolz1/iz,in,tch,potq,pots,gammas
      common/bolz2/nq,ns,nets,ispin,pmass,width
      common/bolz3/pmult(-500:500),pmult1(-500:500),xquark(-6:6)
      common/bolz4/vol,tbary,tdns,tstr,tchrg,tpart,rsd
     $ ,tdns2,tbary2,tstr2,tchrg2
      common/bolz5/xm,xb,mult1(-500:500)
      common/bopt1/icasc,ishape,istrange,imom,mtry
      common/bopt2/epss,frad,fz

      include 'jamntuple.inc'

c=========Set input values and switches ========================
c....Initialize JAM
      mstc(1)=48127      ! random seed.
      write(6,*) 'seed number'
      read(5,*) mstc(1)
      write(6,*) mstc(1)

      mstc(8)=0          ! job mode.
      mstc(16)=0         ! display on/off.
      parc(6)=5.0        ! scale of display
      mstc(54)=0         !avoid first coll inside the same nucleus off

c....Switch on some analysis.
      mstc(156)=1        ! analysis of collision distribution
      mstc(162)=1        ! Output collision histroy
      mstc(165)=1        ! 
      parc(7)= 1.0D0     ! Output time interval (fm/c)

c...Statistical model initial parameters.
      icasc=1            ! hadronic cascade on(1)/off(0)
      imom=1             ! initial momentum dist. 1:Bolzman 2:jet
      ishape=1           ! shape of fireball: 1:sphere, 2:cylinder
      frad=6.0d0         ! radius of the cylinder.
      fz=3.0d0           ! length of the cylinder.
      istrange=1         ! impose netstrangeness on(1)/not(0)
      mtry=100           ! max. number of iteration for strangeness pot.
      epss=1d-8          ! accuracy for net strangeness conservation.
      tch=0.17           ! chemical temperature
      potq=0.27/3.       ! baryon potential
      pots=0.0204496     ! strange potenital
      gammas=1.          ! strangness supression factor
      tpart=-1           ! particle multiplicity (before resonance decay)
                         ! required.
                         ! if tpart<0, multiplicity is determined from
                         ! baryon number conservation specified by
                         ! 'proj' and 'targ'.

      write(6,*) 'Tch [GeV]'
      read(5,*) tch
      write(6,*) tch
      write(6,*) 'mu_q [GeV]'
      read(5,*) potq
      write(6,*) potq
      write(6,*) 'gamma_s'
      read(5,*) gammas
      write(6,*) gammas
      write(6,*) 'multiplicity of total particle'
      read(5,*) tpart
      write(6,*) tpart

c....Initial setting for JAM.
      mevent=50          ! total simulation event
      write(6,*) 'number of event'
      read(5,*) mevent
      write(6,*) mevent

      frame='user'       ! comp. frame in this case, user defined 
      bmin=0.0d0         ! minimum impact parameter (dummy)
      bmax=0.0d0         ! maximum impact parameter (dummy)
      dt=100.d0          ! collision time(fm/c)
      nstep=1            ! time step (i.e. no time step)
      cwin='10gev'        ! initial c.m. energy per nucl, in this case,
                         ! most energitic two-body collisions expected

c      write(6,*)
c     + 'collision type, Pb+Pb=1, Au+Au=2, S+S=3, Si+Si=4, p+p=5'
c      read(5,*) iiiiiii
c      if     (iiiiiii.eq.1) then
c        proj='208Pb'                         ! projectile
c        targ='208Pb'                         ! target
c        write(6,*) 'Pb+Pb is selected'
c      elseif (iiiiiii.eq.2) then
c        proj='179Au'                         ! projectile
c        targ='179Au'                         ! target
c        write(6,*) 'Au+Au is selected'
c      elseif (iiiiiii.eq.3) then
c        proj='32S'                           ! projectile
c        targ='32S'                           ! target
c        write(6,*) 'S+S is selected'
c      elseif (iiiiiii.eq.4) then
c        proj='28Si'                          ! projectile
c        targ='28Si'                          ! target
c        write(6,*) 'Si+Si is selected'
c      elseif (iiiiiii.eq.5) then
c        proj='1p'                            ! projectile
c        targ='1p'                            ! target
c        write(6,*) 'p+p is selected'
c      else
c        write(6,*) 'I do not know the number. STOP.'
c        stop 
c      endif

      if(tpart.lt.0) then
        write(6,*) 'Projectile'
        read(5,*) proj
        write(6,*) proj

        write(6,*) 'Target'
        read(5,*) targ
        write(6,*) targ
      else
c....This is dummy.
        proj='32S'
        targ='32S'
      endif


c================ end input section ==================================
      
c...Initialize jam.
      call jaminit(mevent,bmin,bmax,dt,nstep,frame,proj,targ,cwin)
      nevent=mstc(2)

c....Total neutron number.
      in=mstd(12)-mstd(13)

c....Total charge.
      iz=mstd(13)

c...Initialize analysis.
      call anal1

c...Determine initial condition from a statistical model.
      call init_part

c...Simulation start.
      do iev=1,nevent
        nevt2 = iev   !.. for ntuple by MK
        nevt3 = iev   !.. for ntuple by MK

c...Sampling particle momentum and coordinate.
        call init_mom

        if(mod(iev,100).eq.0) write(6,*)'event=',iev

c...Simulate one event.
        if(icasc.eq.1)then
           call jamevt(iev)    !...generate hadronic cascade
           call anal2          !...Data analysis.
        endif

      end do

c...Final output.
      if(icasc.eq.1) then
        call jamfin

c...Print analysis results.
        call anal3
      endif

c...print out initial phase space.
      fac=1d0/mstc(2)
c...Rapidty distributuons.
      do i=1,10
      call vscale(100+i,fac)
      call vprint(100+i,0,0)
      call vscale(120+i,fac)
      call vprint(120+i,0,1)
      end do

C....Output hadron multiplicities.
      iunit=13
      open(unit=iunit,file='cascade.dat',status='unknown')
      write(iunit,6100)
 6100 format(//'#',78('=')/10x,
     $ 'Final hadron multiplicities after Cascade'
     $ /'#',78('=')/)
      write(iunit,'(''# total density='',g9.3,'' baryon number='',g9.3,
     $              ''Strangeness='',g9.3,
     $              ''total charge='',g9.3)')
     $ tdns2,tbary2*vol,tstr2*vol,tchrg2*vol
      write(iunit,'(''# volume'',g10.3)')vol

      write(iunit,'(''#'',3x,''before cascade'',2x,
     $  ''after cascade'',3x,''ratio(%)'',12x,''KF'',1x,''name'')')

      write(iunit,'(''#'',5x,3(g11.4,3x),7x,''Total particles'')')
     $ pmult1(0)*vol,mult(0)*fac,
     $ 100*(mult(0)*fac-pmult1(0)*vol)/pmult1(0)/vol
      write(iunit,'(''#'',5x,3(g11.4,3x),7x,''Charged particles'')')
     $ pmult1(41)*vol,mult(41)*fac,
     $ 100*(mult(41)*fac-pmult1(41)*vol)/pmult1(41)/vol
      write(iunit,'(''#'',5x,3(g11.4,3x),7x,''Negative particles'')')
     $ pmult1(42)*vol,mult(42)*fac
     $,100*(mult(42)*fac-pmult1(42)*vol)/pmult1(42)/vol
      write(iunit,'(''#'',5x,3(g11.4,3x),7x,''Positive particles'')')
     $ pmult1(43)*vol,mult(43)*fac,
     $ 100*(mult(43)*fac-pmult1(43)*vol)/pmult1(43)/vol

      inum=0
      do i=1,500
       if(i.eq.0.or.(i.ge.41.and.i.le.43)) goto 510
         ipart1=0
         ipart2=0
         if(pmult1(i).gt.0d0.or.mult(i).gt.0) ipart1=1
         if(pmult1(-i).gt.0d0.or.mult(-i).gt.0) ipart2=1 
         if(ipart1.ge.1.or.ipart2.ge.1) then
         kf=kchg(i,4)
         if(ipart1.eq.1) then
           inum=inum+1
           ratio=0d0
           part0=pmult1(i)*vol
           part1=mult(i)*fac
           if(part0.gt.1d-7) ratio=100*(part1-part0)/part0
           write(iunit,'(i5,1x,3(g11.4,3x),i9,1x,a10)')
     $    inum,part0,part1,ratio,kf,chaf(i,1)
         endif
         if(ipart2.eq.1) then
           inum=inum+1
           ratio=0d0
           part0=pmult1(-i)*vol
           part1=mult(-i)*fac
           if(part0.gt.1d-7) ratio=100*(part1-part0)/part0
           write(iunit,'(i5,1x,3(g11.4,3x),i9,1x,a10)')
     $    inum,part0,part1,ratio,-kf,chaf(i,2)
         endif
       endif
  510 end do

      close(iunit)

      iui=10
C....Monte Carlo sampled hadron multiplicities.
      write(iui,6300)
 6300 format(//1x,78('=')/10x,
     $  'Sampled hadron multiplicities'/1x,78('=')/)
      write(iui,*)'total particle',mult1(0)*fac
      write(iui,*)'charged particle',mult1(41)*fac
      write(iui,*)'negative particle',mult1(42)*fac
      write(iui,*)'positive particle',mult1(43)*fac
      write(iui,*)'baryons mesons',xb*fac,xm*fac
      do i=1,500
       if(i.eq.0.or.(i.ge.41.and.i.le.43)) goto 530
       if(mult1(i).ge.1.or.mult1(-i).ge.1) then
         kf=kchg(i,4)
         if(kchg(i,3).ne.0) then
           write(iui,'(i5,1x,i9,2x,2(a16,2x,g12.3,2x))')
     $      i,kf,chaf(i,1),mult1(i)*fac,chaf(i,2),mult1(-i)*fac
         else
           write(iui,'(i5,1x,i9,2x,a16,2x,g12.3)')
     $      i,kf,chaf(i,1),mult1(i)*fac
         endif
       endif
  530 end do

      end

c***********************************************************************

      subroutine init_mom

c...Sample particle accroding to the statisitcal model initial
c...condition calculated in subr.init_part.
      include 'jam1.inc'
      include 'jam2.inc'
      real*8 jamdtim
      common/bolz1/iz,in,tch,potq,pots,gammas
      common/bolz2/nq,ns,nets,ispin,pmass,width
      common/bolz3/pmult(-500:500),pmult1(-500:500),xquark(-6:6)
      common/bolz4/vol,tbary,tdns,tstr,tchrg,tpart,rsd
     $ ,tdns2,tbary2,tstr2,tchrg2
      common/myana2/wy,wp,ylab
      common/bolz5/xm,xb,mult1(-500:500)
      common/bopt1/icasc,ishape,istrange,imom,mtry
      common/bopt2/epss,frad,fz
      logical first
      data first/.true./
      save first

      include 'jamntuple.inc'

      rad=vol**(1d0/3d0)
      nv=int(tpart)
c     print *,'rad nv=',rad,nv
      mstd(12)=0
      mstd(13)=0
      mstd(14)=0
      nbary=0
      nmeson=0

c....For check.
      if(first) then
        xb=0d0
        xm=0d0
        do i=-500,500
          mult1(i)=0
        end do
        first=.false.
      endif

      do ip=1,nv

c....Fist determine particle species.
        xpart=tdns*rn(0)
        do i=-500,500
          if(i.eq.0.or.(i.ge.41.and.i.le.43)) goto 610
          if(pmult(i).le.0d0) goto 610
          xpart=xpart-pmult(i)
          if(xpart.le.0d0) then
            kc=abs(i)
            kf=kchg(kc,4)*isign(1,i)
            goto 520
          endif
  610   end do
  520   continue

        ibary=isign(kchg(kc,6),kf)
        mstd(12)=mstd(12)+ibary/3
        mstd(13)=mstd(13)+jamchge(kf)/3
        mstd(14)=mstd(14)+kchg(kc,7)*isign(1,kf)
        if(ibary.eq.0) then
          nmeson=nmeson+1
          xm=xm+1
        else
          nbary=nbary+1
          xb=xb+1
        endif

c...Zero the vector.
        call jamzero(ip)

c....Particle mass.
        pm=pjmass(kf)
        if(pmas(kc,2).le.1d-7.or.mdcy(kc,1).eq.0
     $              .or.mdcy(kc,2).eq.0.or.mdcy(kc,3).eq.0)then
          k(1,ip)=1
        else
          k(1,ip)=2
        endif
        k(2,ip)=kf
        k(3,ip)=0
        k(4,ip)=0
        k(5,ip)=-1
        k(6,ip)=0
        k(7,ip)=1
        k(8,ip)=1
        k(9,ip)=ibary
        k(10,ip)=0
        k(11,ip)=0

c.....Generate the energy from the local thermal distribution.
          if(imom.eq.1) then
1000     e=rn(0)
         e=e*rn(0)
         e=e*rn(0)
         if(e.le.0d0) goto 1000
         e=-tch*log(e)
         if(rn(0).gt.exp((e-sqrt(e**2+pm**2))/tch)) goto 1000

c         e=-tch*log((1.0d0-rn(0))*exp(-pm/tch))

c         pp=sqrt(abs(e**2-pm**2))

          cost=2d0*rn(0)-1d0
          sint=sqrt(1d0-cost**2)
          phi=paru(2)*rn(0)
          p(1,ip)=e*sint*cos(phi)
          p(2,ip)=e*sint*sin(phi)
          p(3,ip)=e*cost
        else if(imom.eq.2) then
          write(6,*)'sorry not yet implemented option: imom=',imom
          stop
        else
          write(6,*)'wrong option: imom=',imom
          stop
        endif

c       p(4,ip)=e
        p(4,ip)=sqrt(e**2+pm**2)
        p(5,ip)=pm
c       where variable e means momentum

c....Uniform fireball.
        if(ishape.eq.1) then
          rr=rad*rn(0)**(1d0/3d0)
          cx=1.d0-2.d0*rn(0)
          sx=sqrt(1.d0-cx**2)
          phi=paru(2)*rn(0)
          r(1,ip)=rr*sx*cos(phi)
          r(2,ip)=rr*sx*sin(phi)
          r(3,ip)=rr*cx
        else if(ishape.eq.2) then
          rr=frad*sqrt(rn(0))
          phi=paru(2)*rn(0)
          r(1,ip)=rr*cos(phi)
          r(2,ip)=rr*sin(phi)
          r(3,ip)=-fz/2+fz*rn(0)
        else
          write(6,*)'wrong option: ishape=',ishape
          stop
        endif
        r(4,ip)=0d0
        r(5,ip)=0d0

c...Vertex
        v(1,ip)=r(1,ip)
        v(2,ip)=r(2,ip)
        v(3,ip)=r(3,ip)
        v(4,ip)=r(4,ip)

c.....Set resonance decay time.
        v(5,ip)=1.d+35
        if(k(1,ip).eq.2)
     $  v(5,ip)=r(4,ip)+jamdtim(1,kf,kc,k(1,ip),p(5,ip),p(4,ip))


      end do



c...C.M.correction.
      cx=0.d0
      cy=0.d0
      cz=0.d0
      px=0.d0
      py=0.d0
      pz=0.d0
      s=0.d0
      do i=1,nv
        px=px+p(1,i)
        py=py+p(2,i)
        pz=pz+p(3,i)
        cx=cx+r(1,i)*p(5,i)
        cy=cy+r(2,i)*p(5,i)
        cz=cz+r(3,i)*p(5,i)
        s=s+p(5,i)
      end do

      cx=-cx/s
      cy=-cy/s
      cz=-cz/s
      px=-px/nv
      py=-py/nv
      pz=-pz/nv

      do i=1,nv
        r(1,i)=r(1,i)+cx
        r(2,i)=r(2,i)+cy
        r(3,i)=r(3,i)+cz
        v(1,i)=r(1,i)
        v(2,i)=r(2,i)
        v(3,i)=r(3,i)
        p(1,i)=p(1,i)+px
        p(2,i)=p(2,i)+py
        p(3,i)=p(3,i)+pz
        p(4,i)=sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)

        pt=sqrt(p(1,i)**2+p(2,i)**2)
        y=0.5d0*log(max(p(4,i)+p(3,i),1.d-8)/
     $              max(p(4,i)-p(3,i),1.d-8))

        call vfill1(101,y,1.0d0/wy)
        call vfill1(121,pt,1.d0/(pt*wp))

        kf=k(2,i)
        iii=0
        if(kf.eq.2212)  iii=2                              ! Protons.
        if(kf.eq.-2212) iii=3                              ! anti protons.
        if(abs(kf).eq.211.or.kf.eq.111)  iii=5             ! pi
        if(abs(kf).eq.321.or.abs(kf).eq.-311) iii=6        ! k
        if(kf.eq.3122)  iii=7                              ! lambda
        if(kf.eq.3112.or.kf.eq.3212.or.kf.eq.3222) iii=8   ! sigma
        if(kf.eq.3312.or.kf.eq.3312)  iii=9                ! xi
        if(kf.eq.3334)  iii=10                             ! omega

        if(iii.ne.0) then
	  call vfill1(iii+100,y,1.0d0/wy)
          call vfill1(iii+120,pt,1.d0/(pt*wp))
        endif
        if(iii.eq.2.or.iii.eq.3) then
	  call vfill1(104,y,isign(1,kf)/wy)
          call vfill1(124,pt,isign(1,kf)/(pt*wp))
        endif

c...Count hadron multiplicites.
        kc=jamcomp(kf)
        kch=jamchge(kf)
        mult1(0)=mult1(0)+1
        if(kch.ne.0) mult1(41)=mult1(41)+1
        if(kch.lt.0) mult1(42)=mult1(42)+1
        if(kch.gt.0) mult1(43)=mult1(43)+1
        mult1(kc*isign(1,kf))=mult1(kc*isign(1,kf))+1
 300  end do

c... for ntuple by MK
      nptcl2 = nv
      inptcl2 = (nptcl2-mod(nptcl2,5000))/5000 + 1
      do isub2=1,inptcl2
        if     (inptcl2.eq.1) then
          thtrk = nptcl2
        elseif (isub2.eq.inptcl2) then
          thtrk = mod(nptcl2,5000)
        else
          thtrk = 5000
        endif
        do iiii = 1,thtrk
          iiiii = iiii+(isub2-1)*5000
          thpdgid(iiii) =      k(2,iiiii)
          thpx(iiii)    = sngl(p(1,iiiii))
          thpy(iiii)    = sngl(p(2,iiiii))
          thpz(iiii)    = sngl(p(3,iiiii))
          thm(iiii)     = sngl(p(5,iiiii))

          thx(iiii)     = sngl(r(1,iiiii))
          thy(iiii)     = sngl(r(2,iiiii))
          thz(iiii)     = sngl(r(3,iiiii))
          tht(iiii)     = sngl(r(4,iiiii))
        enddo
        call HFNT(2)                       ! store data to ntuple, by MK
      enddo


      end

c***********************************************************************

      subroutine init_part

c...Calculate particle species from statistical model.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      common/bolz1/iz,in,tch,potq,pots,gammas
      common/bolz2/nq,ns,nets,ispin,pmass,width
      common/bolz3/pmult(-500:500),pmult1(-500:500),xquark(-6:6)
      common/bolz4/vol,tbary,tdns,tstr,tchrg,tpart,rsd
     $ ,tdns2,tbary2,tstr2,tchrg2
      common/bopt1/icasc,ishape,istrange,imom,mtry
      common/bopt2/epss,frad,fz
      dimension idlc(100),branch(100)
      dimension strtmp(2),potstmp(2)

      include 'jamntuple.inc'

      if(istrange.eq.0) then
        call init_pmult
        goto 5300
      endif

c....Initial strangeness chemical potentials.
      potstmp(1)=-100.5d0
      potstmp(2)= 100.5d0
      do ib=1,2
        pots=potstmp(ib)
        call init_pmult
        strtmp(ib)=tstr
      end do

c...Check initial condition.
      if(strtmp(1)*strtmp(2).gt.0d0) then
        write(6,*)'wrong initial condition for strangeness pot',
     $     'try again with different values of pots'
        write(6,*)'pots=',potstmp(1),potstmp(2)
        write(6,*)'str =',strtmp(1),strtmp(2)
        stop
      endif
      pots=(potstmp(1)+potstmp(2))/2

c...Loop over iteration for strangeness chemical potential to make
c...total strangeness number equals zero.
      itry=0
 5000 continue
      itry=itry+1
      if(itry.ge.mtry) then
        write(6,*)'strangness potential not converge after iteration of',ntry
        stop
      endif

      call init_pmult
      if(abs(tstr).lt.epss) goto 5300

       PRINT *,'total strangeness=',itry,pots,tstr

      if(tstr*strtmp(1).le.0d0) then
        potstmp(2)=pots
        strtmp(2)=tstr
        pots=(pots+potstmp(1))/2
        goto 5000
      else if(tstr*strtmp(2).le.0d0) then
        potstmp(1)=pots
        strtmp(1)=tstr
        pots=(pots+potstmp(2))/2
        goto 5000
      else
        write(6,*)'wrong initial strangeness potential'
        write(6,*)potstmp(1),potstmp(2)
        stop
      endif

 5300 continue

c....Calculate the volume of fireball from baryon number conservation.
      if(tpart.lt.0.0d0) then
        vol=(iz+in)/tbary
        tpart=tdns*vol
      else
        vol=tpart/tdns
      endif

c...Decay all resonances.
      do kc=-500,500
        pmult1(kc)=pmult(kc)
      end do

      ntry=0
 400  continue
      ntry=ntry+1
      if(ntry.ge.500) then
        write(6,*)'infinite loop in decay?'
        goto 410
      endif
      idec=0

c....Loop over all particles.
      do 2000 kc=1,500
        if(mdcy(kc,1).ne.1) goto 2000
        if(pmas(kc,2).lt.0.001d0) goto 2000
        if(pmult1(kc).le.0d0.and.pmult1(-kc).le.0d0) goto 2000
        kfa=kchg(kc,4)
        np=2
        if(kchg(kc,3).eq.0) np=1

c....Loop for particle/anti-particle.
        do 2100 ib=1,np

        if(ib.eq.1) then
          if(pmult1(kc).le.0d0) goto 2100
          dns0=pmult1(kc)
          kf=kfa
        else if(ib.eq.2) then
          if(pmult1(-kc).le.0d0) goto 2100
          dns0=pmult1(-kc)
          kf=-kfa
        endif
        kfs=isign(1,kf)

C...Check existence of decay channels. Particle/antiparticle rules.
        kca=kc
        if(mdcy(kc,2).gt.0) then
          mdmdcy=mdme(mdcy(kc,2),2)
          if(mdmdcy.gt.80.and.mdmdcy.le.90) kca=mdmdcy
        endif
        if(mdcy(kca,2).le.0.or.mdcy(kca,3).le.0) then
          goto 2000
        endif

        if(mod(kfa/1000,10).eq.0.and.kca.eq.85) kfs=-kfs
        if(kchg(kc,3).eq.0) then
          kfsp=1
          kfsn=0
          if(rn(0).gt.0.5d0) kfs=-kfs
        elseif(kfs.gt.0) then
          kfsp=1
          kfsn=0
        else
          kfsp=0
          kfsn=1
        endif

        brsu=0d0
        nope=0
        do 230 idl=mdcy(kc,2),mdcy(kc,2)+mdcy(kc,3)-1
          if(mdme(idl,1).ne.1.and.kfsp*mdme(idl,1).ne.2.and.
     &    kfsn*mdme(idl,1).ne.3) goto 230
          if(mdme(idl,2).gt.100) goto 230
          if(mdme(idl,1).ne.1.and.kfsp*mdme(idl,1).ne.2.and.
     &       kfsn*mdme(idl,1).ne.3) then
           if(idl.lt.mdcy(kca,2)+mdcy(kca,3)-1) goto 230
          elseif(mdme(idl,2).gt.100) then
           if(idl.lt.mdcy(kca,2)+mdcy(kca,3)-1) goto 230
          endif

C...Read out decay products. Convert to standard flavour code.
          jtmax=5
          if(mdme(idl+1,2).eq.101) jtmax=10
          do jt=1,jtmax
             if(jt.le.5) kp=kfdp(idl,jt)
             if(jt.ge.6) kp=kfdp(idl+1,jt-5)
             if(kp.eq.0) goto 2300
             kpa=iabs(kp)
             kcp=jamcomp(kpa)
             if(kchg(kcp,3).eq.0.and.kpa.ne.81.and.kpa.ne.82) then 
             elseif(kpa.ne.81.and.kpa.ne.82) then 
             else
              goto 2100
             endif
 2300     end do

          nope=nope+1
          if(nope.gt.100) then
            write(6,*)'decay branch over flow: dimension of branch()'
            stop
          endif
          brsu=brsu+brat(idl)
          branch(nope)=brat(idl)
          idlc(nope)=idl
  230   continue
        if(nope.eq.0) goto 2100

          ichg0=jamchge(kf)
          chg0=ichg0*dns0
          dnsxx=0d0
          chg1=0d0
c.....loop over all decay channels.
        do 2200 idl=1,nope
           ichg1=0
           idc=idlc(idl)
           dns1=branch(idl)/brsu*dns0
           dnsxx=dnsxx+dns1
           jtmax=5
           if(mdme(idc+1,2).eq.101) jtmax=10
           ikaon=0
           do 2400 jt=1,jtmax
             if(jt.le.5) kp=kfdp(idc,jt)
             if(jt.ge.6) kp=kfdp(idc+1,jt-5)
             if(kp.eq.0) goto 2400

c...Convert Ks, KL into k0, k0bar
             if(kp.eq.130.or.kp.eq.310) then
               if(ikaon.eq.0) then
                  kp=311
                  ikaon=1
                  if(rn(0).gt.0.5d0) then
                    kp=-311
                    ikaon=-1
                  endif
               else if(ikaon.eq.1) then
                 kp=-311
                 ikaon=0
               else if(ikaon.eq.-1) then
                 kp=311
                 ikaon=0
               endif
             endif

             kpa=iabs(kp)
             kcp=jamcomp(kpa)
             if(pmas(kcp,2).ge.0.001d0) idec=1
             if(kchg(kcp,3).eq.0.and.kpa.ne.81.and.kpa.ne.82) then 
              kfp=kp
             elseif(kpa.ne.81.and.kpa.ne.82) then 
               kfp=kfs*kp
             endif
             kcc=isign(kcp,kfp)
             pmult1(kcc)=pmult1(kcc)+dns1
             ichg1=ichg1+jamchge(kfp)
             chg1=chg1+jamchge(kfp)*dns1
 2400      continue

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            if(ichg1.ne.ichg0) then
              write(6,*)'ib kf0',ib,kf,ichg0,ichg1,chaf(kc,1)
              write(6,*)'jtmax',jtmax,brsu,branch(idl)

           do 2410 jt=1,jtmax
             if(jt.le.5) kp=kfdp(idc,jt)
             if(jt.ge.6) kp=kfdp(idc+1,jt-5)
             if(kp.eq.0) goto 2410
             kpa=iabs(kp)
             kcp=jamcomp(kpa)
             if(kchg(kcp,3).eq.0.and.kpa.ne.81.and.kpa.ne.82) then 
               kfp=kp
             elseif(kpa.ne.81.and.kpa.ne.82) then 
               kfp=kfs*kp
             endif
             write(6,*)jt,kfp,' ',chaf(abs(kcp),1)
 2410      continue

              stop
            endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 2200   continue

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
         if(abs(chg0-chg1).gt.1d-4.or.abs(dns0-dnsxx).gt.1d-4) then
             print *,'kf chrg0 chrg1',kf,chg0,chg1
             print *,'dns',dns0,dnsxx
            stop
         endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 

        if(ib.eq.1) pmult1(kc)=0d0
        if(ib.eq.2) pmult1(-kc)=0d0
 2100 continue
 2000 continue
      if(idec.eq.1) goto 400


  410 continue
c...Count hadron multiplicites for final hadrons.
      tdns2=0d0
      tbary2=0d0
      tstr2=0d0
      tchrg2=0d0
      pmult1(0)=0d0
      pmult1(41)=0d0
      pmult1(42)=0d0
      pmult1(43)=0d0
      do kc=-500,500
      if(kc.eq.0.or.kc.eq.41.or.kc.eq.42.or.kc.eq.43) goto 3000
      if(pmult1(kc).gt.0d0) then
        kca=abs(kc)
        kf=kchg(abs(kc),4)*isign(1,kc)
        dns=pmult1(kc)
        kch=jamchge(kf)
        pmult1(0)=pmult1(0)+dns
        if(kch.ne.0) pmult1(41)=pmult1(41)+dns
        if(kch.lt.0) pmult1(42)=pmult1(42)+dns
        if(kch.gt.0) pmult1(43)=pmult1(43)+dns
        tdns2=tdns2+dns
        tbary2=tbary2+kchg(kca,6)*dns/3*isign(1,kc)
        tstr2=tstr2+jamflav(kf,3)*dns
        tchrg2=tchrg2+jamchge(kf)*dns/3
      endif 
 3000 end do

      weit=vol
      iunit=10
      open(unit=iunit,file='statistic.dat',status='unknown')
      write(iunit,'(''Statisitcal model inputs:'',/
     $ ''Temperature='',g15.7,'' Baryon potential='',g15.7,
     $ ''Strangeness pot.='',g15.7)')tch, potq*3,pots

      write(iunit,'(''Initial proton '',i4,'' Initial neutron '',i4,
     $  '' Total number '',i4)')iz,in,iz+in

      write(iunit,'(''Outputs:''/)')
      write(iunit,'(''volume of fireball(fm^3)='',g15.7,
     $ '' Its radi='',g15.7''fm'')')vol,vol**(1./3.)
 
      write(iunit,'(/''total particle density'',g15.7,''1/fm^3'',/
     $''total baryon number density='',g15.7,/
     $''total charge density='',g15.7,/
     $''total strangness density='',g15.7/)')tdns,tbary,tchrg,tstr

      write(iunit,'(''total particle='',g15.7,''1/fm^3'',
     $''total baryon number='',g15.7,/
     $''total charge='',g15.7,/
     $''total strangness='',g15.7/)')tpart,tbary*vol,tchrg*vol,tstr*vol

      rsd=0.0d0
      totud=xquark(-1)+xquark(1)+xquark(-2)+xquark(2)
      if(totud.gt.1e-5) then
       rsd=200*(xquark(-3)+xquark(3))/totud
      endif
      write(iunit,820)((xquark(-j)+xquark(j))*weit,j=1,3),rsd
820   format('d+dbar',f10.5,/'u+ubar',f10.5,/'s+sbar',f10.5,/
     $  '2(s+sbar)/(u+ubar+d+dbar)',f10.5,' %')


C....Hadron multiplicities.
      write(iunit,6100)
 6100 format(//1x,78('=')/5x,
     $ 'Hadron multiplicities from statistical model',
     $ ' before resonane decays'/1x,78('=')/)


      write(iunit,'(''total particles'',g10.4,/
     $ ''total charged particles'',g10.4,/
     $ ''total negative particles'',g10.4,/
     $ ''total positive particles'',g10.4/)')
     $ pmult(0)*weit,pmult(41)*weit,pmult(42)*weit,pmult(43)*weit

c...Calculate baryon number.
      ybar=0
      ymes=0
      do i=1,500
        if(i.ge.41.and.i.le.43) then
        else
          kf=kchg(i,4)
          if(kchg(i,6).gt.0) then
            ybar=ybar+pmult(i)*kchg(i,6)-pmult(-i)*kchg(i,6)
          else
            ymes=ymes+pmult(i)+pmult(-i)
          endif
        endif
      end do

      write(iunit,'(''baryons:'',g11.4,/,''mesons:'',g11.4,/,)')
     $ ybar*weit/3,ymes*weit

      write(iunit,6110)
 6110 format(3x,'KC',8x,'KF',2x,'Name',8x,'Mult.',5x,'Prob.',
     $ 11x,'Name',4x,'Mult.',6x,'Prob.')


      do i=1,500
       if(i.eq.0.or.(i.ge.41.and.i.le.43)) goto 510
       if(pmult(i).gt.0d0.or.pmult(-i).gt.0d0) then
         kf=kchg(i,4)
         if(kchg(i,3).ne.0) then
           write(iunit,'(i5,1x,i9,2x,2(a10,2x,g9.3,1x,g9.3,2x))')
     $    i,kf,chaf(i,1),pmult(i)*weit,pmult(i)
     $    ,chaf(i,2),pmult(-i)*weit,pmult(i)
         else
           write(iunit,'(i5,1x,i9,2x,a10,2x,g9.3,1x,g9.3)')
     $      i,kf,chaf(i,1),pmult(i)*weit,pmult(i)
         endif
       endif
  510 end do

C....Final hadron multiplicities after decay.
      write(iunit,6200)
 6200 format(//1x,78('=')/10x,'Final hadron multiplicities after decay'
     $ /1x,78('=')/)

      write(iunit,'(''total particles'',g10.4,/
     $ ''total charged particles'',g10.4,/
     $ ''total negative particles'',g10.4,/
     $ ''total positive particles'',g10.4/)')
     $ pmult1(0)*weit,pmult1(41)*weit,pmult1(42)*weit,pmult1(43)*weit

      do i=1,500
       if(i.eq.0.or.(i.ge.41.and.i.le.43)) goto 710
       if(pmult1(i).gt.0d0.or.pmult1(-i).gt.0d0) then
         kf=kchg(i,4)
         if(kchg(i,3).ne.0) then
           write(iunit,'(i5,1x,i9,2x,2(a10,2x,g9.3,1x,g9.3,2x))')
     $    i,kf,chaf(i,1),pmult1(i)*weit,pmult1(i)
     $    ,chaf(i,2),pmult1(-i)*weit,pmult1(-i)
         else
           write(iunit,'(i5,1x,i9,2x,a10,2x,g9.3,1x,g9.3)')
     $      i,kf,chaf(i,1),pmult1(i)*weit,pmult1(i)
         endif
       endif
  710 end do

      close(10)

      tchem  = sngl(tch)
      muq    = sngl(potq) 
      mus    = sngl(pots)
      gams   = sngl(gammas)
      volume = sngl(vol)
      call HFNT(1)                       ! store  data  to ntuple, by MK

      end

c***********************************************************************

      subroutine init_pmult

c...Calculate particle species from statistical model.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      common/bolz1/iz,in,tch,potq,pots,gammas
      common/bolz2/nq,ns,nets,ispin,pmass,width
      common/bolz3/pmult(-500:500),pmult1(-500:500),xquark(-6:6)
      common/bolz4/vol,tbary,tdns,tstr,tchrg,tpart,rsd
     $ ,tdns2,tbary2,tstr2,tchrg2

c...Initialize particle multiplicities.
      tdns=0d0         ! total particle density
      tbary=0d0        ! total baryon number denstiy
      tstr=0d0         ! total strangeness number density
      tchrg=0d0        ! total charge number density
      do i=-500,500
        pmult(i)=0    ! hadron density.
        pmult1(i)=0   ! hadron density after resonance decay.
      end do

      do i=-6,6
        xquark(i)=0d0
      end do

      kcmin=132
      kcmax=423
c...Loop over particles.
      do kc=kcmin,kcmax
        kf=kchg(kc,4)
        if(kf.eq.130) goto 1000  ! K_L0
        if(kf.eq.310) goto 1000  ! K_S0
        kfl1=mod(kf/1000,10)
        kfl2=mod(kf/100,10)
        kfl3=mod(kf/10,10)

c....Exclude heavy flavour hadrons (current version of statistical model
c....in this program doese not handle those).
        if(kfl1.ge.4) goto 1000
        if(kfl2.ge.4) goto 1000
        if(kfl3.ge.4) goto 1000

        ispin=mod(abs(kf),10)
        if(ispin.eq.0.and.kf.ne.10220) goto 1000

c....Get hadron denisty.
        dns=density(kf,potq,pots,gammas,tch)
        tdns=tdns+dns
        tbary=tbary+kchg(kc,6)*dns/3
        tstr=tstr+jamflav(kf,3)*dns
        tchrg=tchrg+jamchge(kf)*dns/3

c...Count hadron multiplicites.
        pmult(kc)=pmult(kc)+dns
        kch=jamchge(kf)
        pmult(0)=pmult(0)+dns
        if(kch.ne.0) pmult(41)=pmult(41)+dns
        if(kch.lt.0) pmult(42)=pmult(42)+dns
        if(kch.gt.0) pmult(43)=pmult(43)+dns

c....Save quark contenet of the hadron.
        call attflv2(kf,ifla,iflb,iflc)
        if(ifla.ne.0) xquark(ifla)=xquark(ifla)+dns
        if(iflb.ne.0) xquark(iflb)=xquark(iflb)+dns
        if(iflc.ne.0) xquark(iflc)=xquark(iflc)+dns

c       write(6,'(i7,1x,a16,2x,g12.4,2(1x,g10.3),4(1x,i3))')
c    $ kf,chaf(kc,1),dns,pmass,width,ispin,nq,ns,nets

c......Anti-particle if it exists.
        if(kchg(kc,3).eq.1) then
          dns=density(-kf,potq,pots,gammas,tch)
          tdns=tdns+dns
          tbary=tbary-kchg(kc,6)*dns/3
          tstr=tstr+jamflav(-kf,3)*dns
          tchrg=tchrg+jamchge(-kf)*dns/3

c...Count hadron multiplicites.
          pmult(-kc)=pmult(-kc)+dns
          kch=jamchge(-kf)
          pmult(0)=pmult(0)+dns
          if(kch.ne.0) pmult(41)=pmult(41)+dns
          if(kch.lt.0) pmult(42)=pmult(42)+dns
          if(kch.gt.0) pmult(43)=pmult(43)+dns

c....Save quark contenet of the hadron.
          call attflv2(-kf,ifla,iflb,iflc)
          if(ifla.ne.0) xquark(ifla)=xquark(ifla)+dns
          if(iflb.ne.0) xquark(iflb)=xquark(iflb)+dns
          if(iflc.ne.0) xquark(iflc)=xquark(iflc)+dns

c         write(6,'(i7,1x,a16,2x,g12.4,2(1x,g9.4),4(1x,i3))')
c    $    -kf,chaf(kc,2),dns,pmass,width,ispin,nq,ns,nets
        endif

 1000 end do

      end

c************************************************************************

      real*8 function density(kf,potq,pots,gammas,tch)
c
c     Boltzmann approximation of hadron gas model
c
c     author: Masashi Kaneta
c
c      Input:
c
c    potq  ;  light   quark chemical potential [GeV]
c    pots  ;  strange quark chemical potential [GeV]
c    gammas;  strangeness supresion factor 
c    tch   ;  chemical freeze-put temperature  [GeV]
c    kf    ;  particle KF code  (Particle data group code)
c
c    Output:
c
c    density:  particle density [1/fm^3]
c 
c (width>1MeV)
c   =  (1/2 Pi^2) g tch^3 lambda_q lambda_s  gammas^nets
c  Integrate[ (x/tch)^2 K2(x/tch) x /((x^2-mass^2)^2+(x Gam)^2),
c                                  {x,mass-2*width,mass+2*width} ]
c /Integrate[                     x /((x^2-mass^2)^2+(x Gam)^2),
c                                  {x,mass-2*width,mass+2*width} ]
c
c (width<=1MeV)
c   = (1/2 Pi^2) g tch^3 lambda_q lambda_s  gammas^nets  (mass/tch)^2 K2(mass/tch)
c

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      parameter (fgev= 0.197327,pimass=0.14d0)
      parameter (n=500)
      real     BesK0, BesK1,w
      external BesK0, BesK1
      common/bolz2/nq,ns,nets,ispin,pmass,width

c...Preparation.
c...mass  ;  particle mass                    [GeV]
c...width ;  width of mass (B-W type)         [GeV]
c...ispin ;  spi freedum = 2J+1 
c...nq    ;  number of u,d quark, for u-bar,d-bar its shuld be negative value
c...ns    ;  number of s quark, for s-bar its shuld be negative value
c...nets  ;  net strangeness

      kc=jamcomp(kf)
      pmass=pmas(kc,1)                 ! particle mass
      width=pmas(kc,2)                 ! particle width
      kfa=abs(kf)
      ispin=max(1,mod(kfa,10))         ! spin  2J+1
      nq=jamflav(kf,1)+jamflav(kf,2)   ! net u+d
      ns=jamflav(kf,3)                 ! net s
      ibary=kchg(kc,6)                 ! baryon number times 3.

      nets=0
c....Mesons.
      if(ibary.eq.0) then
        if(kfa.eq.221.or.kfa.eq.331) then  ! eta, eta'
           nets=1
        else
          kf1=mod(kfa/100,10)
          kf2=mod(kfa/10,10)
          kfm=100*kf1+10*kf2
          kfms=kfm+ispin
          if(kfm.eq.220) then          ! omega mesons
          else if(kfm.eq.330) then     ! phi, f 
            nets=2
          else if(kfms.eq.223.or.kfms.eq.227) then
          else
            if(kf1.eq.3) nets=nets+1
            if(kf2.eq.3) nets=nets+1
          endif
        endif
c...Baryons.
      else if(abs(ibary).eq.3) then
        nets=abs(ns)
      endif

C-----|--1---------2---------3---------4---------5---------6---------7--
 
       s = 1d0
       xn = pmass-2*width
       xx = pmass+2*width
       dx = (xx-xn)/n
       if(xn.le.0.0d0) then
         if(kf.eq.10220) then   ! sigma meson
           xn = 2*pimass
         else
           write(6,*)'xn<0: mass width',pmass,width
           stop
         endif
       endif
       if (width.gt.0.001d0) then          ! calc. normarization factor
         xi = xn                      ! for B-W mass distribution
         s =  0.5d0*bm(xi,pmass,width)
         do i=2,n
           xi = xn + dx*(i-1)
           s = s + bm(xi,pmass,width)
         enddo
         xi = xx
         s = s + 0.5d0*bm(xi,pmass,width)
         s = dx * s 
 
         xi = xn                          ! calc. particle density
         dmk = 0.5d0*bwm(xi,pmass,width,tch)
         do i=2,n
           xi = xn + dx*(i-1)
           dmk = dmk + bwm(xi,pmass,width,tch)
         enddo
         xi = xx
         dmk = dmk + 0.5d0*bwm(xi,pmass,width,tch)
         area = dx * dmk /s 
       else
         w = real(pmass/tch)
         area = (w**2) * ( BesK0(w)+(2/w)*BesK1(w) )

       endif

 
       density=exp((nq*potq+ns*pots)/tch)*gammas**abs(nets)
     $   *area*(tch/fgev)**3*ispin*0.5d0/paru(1)**2
 
      end 

c    +-----------------------------------
      function bm(x,pmass,width)
      implicit double precision(a-h, o-z)

      bm =  x /( (x**2-pmass**2)**2+(x*width)**2 )

      end 

c    +-----------------------------------
      function bwm(x,pmass,width,tch)
      implicit double precision(a-h, o-z)
      real     besk0, besk1,w
      external besk0, besk1

      w = real(x/tch)
      bwm =  x *  (w**2)* ( BesK0(w)+(2/w)*BesK1(w) )
     +     / ( (x**2-pmass**2)**2 + (x*width)**2 )

      end 


c***********************************************************************

      subroutine anal1

      include 'jam1.inc'
      include 'jam2.inc'
      common/myana1/ispec
      common/myana2/wy,wp,ylab
      save xspec

      common/bolz1/iz,in,tch,potq,pots,gammas
      common/bolz4/vol,tbary,tdns,tstr,tchrg,tpart,rsd
     $ ,tdns2,tbary2,tstr2,tchrg2

      include 'jamntuple.inc'

c... nTuple initialization by MK
      write(6,*) 'CWN file name'
      read(5,*) oname

      call HLIMIT(NWPAWC)
      call HROPEN(2,'JAM',oname,'N',1024,istat)
      if( istat .ne. 0 ) then
        write(6,*) 'Error opening cwn file.'
        call HROUT(0,icycle,' ')
        call HREND('JAM')
        stop
      endif

      call HBNT(1, 'parameters' , ' ')
      call HBNT(2, 'Thermal'    , ' ')
      call HBNT(3, 'JAM cascade', ' ')

      chform1 = 'tch:R, muq:R, mus:R, gammas:R, vol:R '

      chform2 = 
     +  'nevt:I, nptcl:I, isub:I, trk[0,5000]:I, '//
     +  'pdgid(trk):I, '//
     +  'px(trk):R, py(trk):R, pz(trk):R, m(trk):R, '//
     +  'x(trk):R, y(trk):R, z(trk):R, t(trk):R '

      chform3 = 
     +  'nevt:I, nptcl:I, isub:I, trk[0,5000]:I, '//
     +  'pdgid(trk):I, '//
     +  'px(trk):R, py(trk):R, pz(trk):R, m(trk):R, '//
     +  'x(trk):R,  y(trk):R, z(trk):R, t(trk):R '

      call HBNAME(1, 'JAMNT1', tchem, chform1)
      call HBNAME(2, 'JAMNT2', nevt2, chform2)
      call HBNAME(3, 'JAMNT3', nevt3, chform3)

      xspec=0

c....Rapidity distribution.
      ylab=pard(17)
      print *,'ylab',ylab
      ymin=-7.0d0
      ymax=7.0d0
      wy=0.25d0
      nymx=(ymax-ymin)/wy

      yminl=ymin+ylab
      ymaxl=ymax+ylab
      print *,'yminl ymaxl ylab',yminl,ymaxl,ylab

      pmin=0.0d0
      pmax=5.0d0
c     wp=0.1d0
      wp=0.2d0
      npmx=(pmax-pmin)/wp

      call vbook1(1,'dn/dy total'      ,nymx,ymin,ymax)
      call vbook1(2,'dn/dy protons'    ,nymx,ymin,ymax)
      call vbook1(3,'dn/dy antiprotons',nymx,ymin,ymax)
      call vbook1(4,'dn/dy net protons',nymx,ymin,ymax)
      call vbook1(5,'dn/dy pi'         ,nymx,ymin,ymax)
      call vbook1(6,'dn/dy k'          ,nymx,ymin,ymax)
      call vbook1(7,'dn/dy lambda'     ,nymx,ymin,ymax)
      call vbook1(8,'dn/dy sigma'      ,nymx,ymin,ymax)
      call vbook1(9,'dn/dy xi'         ,nymx,ymin,ymax)
      call vbook1(10,'dn/dy Omega'     ,nymx,ymin,ymax)

      call vbook1(19,'dETdy' ,nymx,ymin,ymax)
      call vbook1(20,'dETdy' ,nymx,ymin,ymax)

      call vbook1(21,'1/p/dpdy total',npmx,pmin,pmax)
      call vbook1(22,'1/p/dpdy protons',npmx,pmin,pmax)
      call vbook1(23,'1/p/dpdy antiprotons',npmx,pmin,pmax)
      call vbook1(24,'1/p/dpdy net potons',npmx,pmin,pmax)
      call vbook1(25,'1/p/dpdy pion',npmx,pmin,pmax)
      call vbook1(26,'1/p/dpdy k',npmx,pmin,pmax)
      call vbook1(27,'1/p/dpdy lambda',npmx,pmin,pmax)
      call vbook1(28,'1/p/dpdy sigma',npmx,pmin,pmax)
      call vbook1(29,'1/p/dpdy xi',npmx,pmin,pmax)
      call vbook1(30,'1/p/dpdy Omega',npmx,pmin,pmax)


c...for initial condition.
      call vbook1(101,'dn/dy total',nymx,ymin,ymax)
      call vbook1(102,'dn/dy protons',nymx,ymin,ymax)
      call vbook1(103,'dn/dy antiprotons',nymx,ymin,ymax)
      call vbook1(104,'dn/dy net protons',nymx,ymin,ymax)
      call vbook1(105,'dn/dy pi',nymx,ymin,ymax)
      call vbook1(106,'dn/dy k',nymx,ymin,ymax)
      call vbook1(107,'dn/dy lambda',nymx,ymin,ymax)
      call vbook1(108,'dn/dy sigma',nymx,ymin,ymax)
      call vbook1(109,'dn/dy xi',nymx,ymin,ymax)
      call vbook1(110,'dn/dy Omega',nymx,ymin,ymax)

      call vbook1(121,'1/p/dpdy total',npmx,pmin,pmax)
      call vbook1(122,'1/p/dpdy protons',npmx,pmin,pmax)
      call vbook1(123,'1/p/dpdy antiprotons',npmx,pmin,pmax)
      call vbook1(124,'1/p/dpdy net potons',npmx,pmin,pmax)
      call vbook1(125,'1/p/dpdy pion',npmx,pmin,pmax)
      call vbook1(126,'1/p/dpdy k',npmx,pmin,pmax)
      call vbook1(127,'1/p/dpdy lambda',npmx,pmin,pmax)
      call vbook1(128,'1/p/dpdy sigma',npmx,pmin,pmax)
      call vbook1(129,'1/p/dpdy xi',npmx,pmin,pmax)
      call vbook1(130,'1/p/dpdy Omega',npmx,pmin,pmax)


      return

c***********************************************************************

      entry anal2

      beta=pard(5)
      gamma=pard(6)
      ispec=0

      catrk = 0                                                      ! by MK
c...Loop over all particles.
      do 3000 i=1,nv

       if(k(1,i).ge.10) goto 3000

        catrk = catrk + 1        ! count number of track after cascade, by MK

c...Exclude spectetor.
c       if(abs(k(7,i)).eq.1) then
c         y=0.5*log(max(p(4,i)+p(3,i),1.e-8)/max(p(4,i)-p(3,i),1.e-8))
c         write(3,*)k(1,i),k(2,i),k(7,i),y
c         ispec=ispec+1
c         goto 3000
c       endif

        kf=k(2,i)
        y=0.5d0*log( max(p(4,i)+p(3,i),1.d-8)/max(p(4,i)-p(3,i),1.d-8) )
        yl=y
        if(mstc(4).eq.0) then
        else if(mstc(4).eq.3.or.mstd(4).eq.100) then
        else
         yl=y+ylab
        endif

        ptsq=p(1,i)**2+p(2,i)**2
        pt=sqrt(ptsq)
        pp=sqrt(ptsq+p(3,i)**2)
	ee=sqrt(p(5,i)**2+pp**2)
	if(abs(ee-p(4,i)).gt.1e-7) then
	  print *,'ee p4',ee,p(4,i),kf
	endif
        pt=max(pt,1.d-8)
        eta=0.5d0*log( max(pp+p(3,i),1.d-8)/max(pp-p(3,i),1.d-8) )
	et=p(4,i)*pt/max(pp,1.d-8)
        emt=sqrt(p(5,i)**2+ptsq)
        emt0=emt-p(5,i)

        plab  = gamma*( p(3,i) + beta * p(4,i) )
        elab  = gamma*( p(4,i) + beta * p(3,i) )
	ppl=sqrt(plab**2+ptsq)
	etlab=elab*pt/max(ppl,1d-8)
	eta_lab=0.5d0*log(max(ppl+plab,1.d-8)/max(ppl-plab,1.d-8))
	etl=eta_lab+ylab
	call vfill1(19,etl,etlab/wy)
	call vfill1(20,y,et/wy)

c       print *,'gam bet',gamma,beta
c       print *,'plab elab',plab,elab,p(3,i),p(4,i)
c       print *,'etlab eta ylab',etlab,eta_lab,ylab
c       pause


        call vfill1( 1, yl, 1.0d0/wy    )
        call vfill1(21, pt, 1.d0/(pt*wp))

        iii=0
        if(kf.eq.2212)  iii=2    !...Protons.
        if(kf.eq.-2212) iii=3    !...anti protons.
        if(abs(kf).eq.211.or.kf.eq.111)  iii=5    ! pi
        if(abs(kf).eq.321.or.abs(kf).eq.-311)  iii=6  ! k
        if(kf.eq.3122)  iii=7  ! lambda
        if(kf.eq.3112.or.kf.eq.3212.or.kf.eq.3222)  iii=8  ! sigma
        if(kf.eq.3312.or.kf.eq.3312)  iii=9  ! xi
        if(kf.eq.3334)  iii=10  ! omega

        if(iii.ne.0) then
	  call vfill1(iii,yl,1.0d0/wy)
          call vfill1(iii+20,pt,1.d0/(pt*wp))
        endif
        if(iii.eq.2.or.iii.eq.3) then
	  call vfill1(4,yl,isign(1,kf)/wy)
          call vfill1(24,pt,isign(1,kf)/(pt*wp))
        endif

3000  continue

c... for ntuple by MK
      nptcl3 = nv
      inptcl3 = (nptcl3-mod(nptcl3,5000))/5000 + 1
      do isub3=1,inptcl3
        if     (inptcl3.eq.1) then
          catrk = nptcl3
        elseif (isub3.eq.inptcl3) then
          catrk = mod(nptcl3,5000)
        else
          catrk = 5000
        endif
        do iiii = 1,catrk
          iiiii = iiii+(isub3-1)*5000
          if (ifpd(p(1,iiiii)).eq.0) then
            write (*,*)  iiiii,":     ",
     +                   k(1,iiiii),k(2,iiiii),k(3,iiiii),k(4,iiiii),
     +                   p(1,iiiii),p(2,iiiii),p(3,iiiii),p(4,iiiii),
     +                   p(5,iiiii),
     +                   r(1,iiiii),r(2,iiiii),r(3,iiiii),r(4,iiiii),
     +                   r(5,iiiii)
          endif
          capdgid(iiii) =      k(2,iiiii)
          capx(iiii)    = sngl(p(1,iiiii))
          capy(iiii)    = sngl(p(2,iiiii))
          capz(iiii)    = sngl(p(3,iiiii))
          cam(iiii)     = sngl(p(5,iiiii))

          cax(iiii)     = sngl(r(1,iiiii))
          cay(iiii)     = sngl(r(2,iiiii))
          caz(iiii)     = sngl(r(3,iiiii))
          cat(iiii)     = sngl(r(4,iiiii))
        enddo
        call HFNT(3)                       ! store data to ntuple, by MK
      enddo



      xspec=xspec+ispec

      return

c***********************************************************************

      entry anal3

c...Output histograms.

      print *,'event',mstc(2)
      wevt=1.d0/dble(mstc(2))
      fac=wevt
      mnorm=0
      mform=1

      do i=1,2
      call vscale(18+i,fac)
      call vprint(18+i,0,1)
      end do

c...Rapidty distributuons.
      do i=1,10
      call vscale(i,fac)
      call vprint(i,0,0)
      call vscale(20+i,fac)
      call vprint(20+i,0,1)
      end do

c... for ntuple by MK
      call HROUT(0,icycle,' ')
      call HREND('JAM')


      end

