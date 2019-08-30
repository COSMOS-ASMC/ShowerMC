c...A main program for checking p+p reactions: energy dependence of
c...exclusive cross sections.

      include 'jam1.inc'
      include 'jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15
      character plabc*11
      common/pppp/plab
c...Function:lab. momentum.
      plabsr(a,b,c)=sqrt((a**2-b**2-c**2)**2/(4.d0*c**2)-b**2)

      bmin=0.0D0          ! minimum impact parameter
      bmax=-1.0D0         ! maximum impact parameter
      mstc(1) =254777     ! random seed.
      dt=10.0D0           ! collision time(fm/c)
      mevent=10000        ! total simulation event
      nstep=1
      frame='nn'          ! comp. frame
      proj='p'            ! projectile
      targ='p'            ! target

c...for analysis of resoance production.
      mstc(8)=0        ! job mode.
      mstc(55)= 1      ! 1:frozen resonance.
      mstc(76)= 1      ! string has a lifetime.
      mstc(41)=0       ! forbid unstable particle decay.
      mstc(13)=2       ! all warning are printed.
      mstc(17)=0       ! only inelastic collisions or not.
      parc(7)=100.0d0  ! output interval.

      em1=0.938D0
      em2=0.938D0

c     pmin=2.3D0
c     pmax=4.0D0
c     npmax=10

      pmin=1.9D0
      pmax=8.0D0
      npmax=20

c     pmin=1.9
c     pmax=2.6
c     npmax=15

      delp=(pmax-pmin)/npmax

c...Initialize anal for overall run.
      call anal0

c...Loop over incident momentum.
      do ip=1,npmax
c       plab=pmin*(pmax/pmin)**((ip-1)/dble(npmax-1))
c       plab=pmin+(ip-1)*delp
  
c       srt=pmin*(pmax/pmin)**((ip-1)/DBLE(npmax-1))
        srt=pmin+(ip-1)*delp
c       write(66,*)'srt=',srt

        plab=plabsr(srt,em1,em2)


c....Initialize JAM.
        write(plabc,'(f11.3)')plab
        cwin=plabc//'gevc'  ! incident energy
        call jaminit(mevent,bmin,bmax,dt,nstep,frame,proj,targ,cwin)
        nevent=mstc(2)

c...Initialize analysis.
      call anal1

c...Simulation start.
      do iev=1,nevent

c...Simulate one event.
        call jamevt(iev)
        if(mod(iev,1000).eq.0) write(6,*)'event=',iev

c...Data analysis.
        call anal2

      end do

c...Final output.
c     call jamfin

c...Print analysis results.
      call anal3

      end do  ! end loop over momentum

      call anal4

      end

c***********************************************************************

      subroutine anal0

c...Initialize analysis for overall run.

      include 'jam1.inc'
      include 'jam2.inc'
      parameter(ncdet=15,npdet=14) ! stable particles
      parameter(ncdet2=1)          ! unstable particles
      character simfile(ncdet)*15,simfile2(ncdet2)*15,simfile3(npdet)*15
      character cfile(ncdet)*15,cfile2(ncdet2)*15,cfile3(npdet)*15
      dimension ncount(ncdet),ncount2(ncdet2),ncount3(npdet)
      dimension ntyp(npdet), mchan(npdet,ncdet), ktyp(npdet)
      character chaf1*16,chaf2*16

c...For resonance production.
      parameter(nchnl=14)
      character simfile4(nchnl)*15,cfile4(nchnl)*15
      dimension numf4(15),numf41(15),ncount4(nchnl)

      common/pppp/plab

      save ncount,ncount2,ncount3,ncount4
      save sigel,siginel

      data simfile/
     &             'pi1a.dat',    'pi1b.dat'
     &            ,'pi2a.dat',   'pi2b.dat'
     &            ,'pi3a.dat',  'pi3b.dat'
     &            ,'pi4.dat', 'pi5.dat'
     &            ,'lambda1.dat', 'sigma1.dat'
     &            ,'sigma2.dat',  'lambda2.dat'
     &            ,'lambda3.dat',  'eta.dat'
     $            ,'lambda4.dat'/

      data cfile/
     &             'p n pi+',    'p p pion0'
     &            ,'p n pi^+ pi^0', 'p p pi+ pi-'
     &            ,'pn 2pi+ pi-',  'pp pi+ pi0 pi-'
     &            ,'pp 2pi+ 2pi-', 'pn 3pi+ 2pi-'
     &            ,'LpK+',   'S+nK+'
     &            ,'S0pK+',  'Ln+K+'
     &            ,'Lp0K+',  'pp eta'
     $            ,'Lambda p pi+ K0'/

      data simfile2/'omega.dat'/
      data cfile2/'pp-omega'/

      data simfile3/'pp-pX.dat','pp-nX.dat','pp-LambdaX.dat',
     $              'pp-S-X.dat','pp-S0X.dat','pp-S+X.dat',
     $              'pp-pi-X.dat','pp-pi0X.dat','pp-pi+X.dat',
     $              'pp-K-X.dat','pp-Kb0X.dat','pp-K0X.dat',
     $              'pp-K+X.dat','pp-etaX.dat'/

      data cfile3/'p + X','pp-n X','Lambda + X',
     $              'S^- +  X','S^0 + X.dat','S^+ + X',
     $              'pi^- X','pi^0 + X','pi^+ + X',
     $              'K^- X','Kb0 + X','K0 + X',
     $              'K^+ X','eta  + X'/

      data simfile4/'pp-ND.dat','pp-NNs.dat','pp-DD.dat','pp-NDs.dat',
     $              'pp-NsD.dat','pp-DDs.dat','pp-NsNs.dat',
     $              'pp-NsDs.dat','pp-DsDs.dat',
     $              'pp-NR.dat','pp-RR.dat',
     $              'pp-NStr.dat','pp-RStr.dat','pp-StrStr.dat'/
      data cfile4/'N+D(1232)','NN*','DD','ND*',
     $            'N*D','DD*','N*N*','N*D*','D*D*',
     $            'NR','RR',
     $            'N+String','R+String','String+String'/

      data numf4/0,1,3,2,5,7,4,6,8,9,0,0,0,0,0/
      data numf41/0,10,11,10,2*11,10,3*11,12,3*13,14/
c1 n 2 d 3 n* 4 d* 5 str

c2  1 2  nd    1  10
c3  2 2  dd    3  11
c4  3 1  nn*   2  10
c5  3 2  dn*   5  11
c6  3 3  n*n*  7  11
c7  1 4  nd*   4  10
c8  2 4  dd*   6  11
c9  4 3  n*d*  8  11
c10 4 4  d*d*  9  11
c11 1 5  ns       12
c12 2 5           13
c13 3 5           13
c14 4 5           13
c15 5 5  ss       14

c........... n  p  L  S- S0 S+ p- p0 p+ K- Kb K0 K+ eta
      data mchan/
     &       1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   ! np pi+
     &       0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,   ! pp pi0
     &       1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,   ! np pi0 pi+
     &       0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0,   ! pp pi- pi+
     &       1, 1, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0,   ! np pi- 2pi+
     &       0, 2, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0,   ! pp pi- pi0 pi+
     &       0, 2, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0,   ! pp 2pi- 2pi+
     &       1, 1, 0, 0, 0, 0, 2, 0, 3, 0, 0, 0, 0, 0,   ! np 2pi- 3pi+
     &       0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,   ! pL  K+
     &       1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,   ! nS+ K+
     &       0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0,   ! pS0 K+
     &       1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0,   ! nL  pi+ K+
     &       0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,   ! pL pi0 K+
     &       0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,   ! pp eta
     &       0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0 /  ! Lp pi+ K0
      data ktyp
     &      /2112,2212,3122,3112,3212,3222,
     $              -211,111,211,-321,-311,311,321,221/

      ifile=10
      do i=1,ncdet
        ifile=ifile+1
        open(ifile,file=simfile(i))
        write(ifile,'(a)')'# '//cfile(i)
        write(ifile,'(''# plab(gev/c) srt(GeV)  sigma(mb)  mult.'')')
      enddo

      do i=1,ncdet2
        ifile=ifile+1
        open(ifile,file=simfile2(i))
        write(ifile,'(a)')'# '//cfile2(i)
        write(ifile,'(''# plab(gev/c) srt(GeV)  sigma(mb)  mult.'')')
      enddo

      do i=1,npdet
        ifile=ifile+1
        open(ifile,file=simfile3(i))
        write(ifile,'(a)')'# '//cfile3(i)
        write(ifile,'(''# plab(gev/c) srt(GeV)  sigma(mb)  mult.'')')
      enddo

      do i=1,nchnl
        ifile=ifile+1
        open(ifile,file=simfile4(i))
        write(ifile,'(a)')'# '//cfile4(i)
        write(ifile,'(''# plab(gev/c) srt(GeV)  sigma(mb)  mult.'')')
      enddo

      open(1,file='sig.dat',status='unknown')

      return

c***********************************************************************

      entry anal1

      do i=1,ncdet
        ncount(i)=0
      enddo
      do i=1,ncdet2
        ncount2(i)=0
      enddo
      do i=1,npdet
        ncount3(i)=0
      enddo
      do i=1,nchnl
        ncount4(i)=0
      enddo

      sigel=0d0
      siginel=0d0

      return

c***********************************************************************

      entry anal2

c...Count data.
c----------------------------------------------------------------
c...Count resonance production.
      if(nv.ne.2) print *,'nv=',nv,(k(2,i),i=1,nv)
      em1=0.0d0
      em2=0.0d0
      kf1=0
      kf2=0
      kc1=0
      kc2=0

      if(nv.eq.2) then

        kf1=k(2,1)
        kf2=k(2,2)
        k1=k(1,1)
        k2=k(1,2)
        em1=p(5,1)
        em2=p(5,2)
        kc1=jamcomp(kf1)
        kc2=jamcomp(kf2)

        if(k1.eq.1) then
          id1=1
        else if(k1.eq.2) then
          id=kchg(kc1,5)
          if(id.eq.id_delt) id1=2
          if(id.eq.id_nucls) id1=3
          if(id.eq.id_delts) id1=4
        else
          id1=5
        endif

        if(k2.eq.1) then
          id2=1
        else if(k2.eq.2) then
          id=kchg(kc2,5)
          if(id.eq.id_delt)  id2=2
          if(id.eq.id_nucls) id2=3
          if(id.eq.id_delts) id2=4
        else
          id2=5
        endif

        idmin=min(id1,id2)
        idmax=max(id1,id2)
        icpair=(idmax*(idmax-1))/2+idmin
        i1=numf4(icpair)
        i2=numf41(icpair)
        if(i1.ne.0) ncount4(i1)=ncount4(i1)+1
        if(i2.ne.0) ncount4(i2)=ncount4(i2)+1
        call pjname(kf1,chaf1)
        call pjname(kf2,chaf2)
c       write(66,*)kf1,kf2,' ',chaf1,' ',chaf2

      endif
c----------------------------------------------------------------

c....Decay string
      call jamfdec

c...pp omega final
      if(nv.eq.3) then
         if(k(2,1).eq.2212.and.k(2,2).eq.2212.and.k(2,3).eq.223)
     $      ncount2(1)=ncount2(1)+1
      endif

      mstc(41)=1       ! unstable particle decay
      call jamfdec
      mstc(41)=0       ! forbid unstable particle decay

      if(mstc(17).eq.1.and.nv.eq.2) then
        print *,'elastic??? nv=',nv
        do i=1,nv
          write(6,*)k(1,i),k(2,i),p(5,i)
        end do
        stop
      endif

      do i=1,npdet
        ntyp(i)=0
      enddo

c...Loop over all particles.
      notdef=0
      do i=1,nv
        kf=k(2,i)
        itag=0
        do ipdet=1,npdet
          if (kf.eq.ktyp(ipdet)) then
            itag=1
            ntyp(ipdet)=ntyp(ipdet)+1
            ncount3(ipdet)=ncount3(ipdet)+1
          endif
        enddo
        if(itag.eq.0) notdef=1
3000  end do

      ichanel=-99
c...Loop over reaction channel.
      do icdet=1,ncdet
c....Loop over particle involed this reaction channel.
        do ipdet=1,npdet
          if (ntyp(ipdet).ne.mchan(ipdet,icdet)) goto 10
        enddo
        ichanel=icdet
        goto 4000
10      continue
      enddo

4000  continue
      if(notdef.eq.0.and.(ichanel.ge.1.and.ichanel.le.ncdet))
     $   ncount(ichanel)=ncount(ichanel)+1

c...Count elastic collision.
      if(nv.eq.2.and.k(2,1).eq.2212.and.k(2,2).eq.2212) then
        sigel=sigel+1d0
      else
        siginel=siginel+1d0
      endif

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c     if(ichanel.eq.2.or.ichanel.eq.3)then
c      if(nv.ne.3) then
c      write(90,*)' '
c      write(90,*)ichanel,mste(1)
c    $   ,kf1,kf2,em1,em2,chaf(kc1,1),' ',chaf(kc2,1)
c      do i=1,nv
c        kc=jamcomp(k(2,i))
c        write(90,*)k(1,i),k(2,i),p(5,i),' ',chaf(kc,1)
c      end do
c     endif
c     endif
      if(ichanel.eq.5) then
        write(66,*)kf1,kf2,' ',chaf1,' ',chaf2
      endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      return

c***********************************************************************

      entry anal3

c...Output results.

c...Event weight
      wei=1.D0/dble(mstc(2))
      fac=parc(4)**2*paru(1)*10*wei
      srt=pare(2)

c...Inclusive data
      ifn=10
      do i=1,ncdet
        ifn=ifn+1
        write(ifn,800)plab,srt,fac*ncount(i),ncount(i)*wei
      enddo

c...Final unstable particle
      do i=1,ncdet2
        ifn=ifn+1
        write(ifn,800)plab,srt,fac*ncount2(i),ncount2(i)*wei
      enddo

c...Exclusive data
      do i=1,npdet
        ifn=ifn+1
        write(ifn,800)plab,srt,fac*ncount3(i),ncount3(i)*wei
      enddo

c...Resonance productions.
      do i=1,nchnl
        ifn=ifn+1
        write(ifn,800)plab,srt,fac*ncount4(i),ncount4(i)*wei
      enddo

      write(1,810)plab,srt,pare(4),pare(5),sigel*fac
     $ ,pare(4)-pare(5),siginel*fac


800   format(2(f8.3,1x),f12.7,1x,f9.4)
810   format(2(f8.3,1x),10(f10.7,1x))

      return

c***********************************************************************

      entry anal4

      ifn=10
      do i=1,ncdet+ncdet2+npdet+nchnl
        ifn=ifn+1
        close(ifn)
      enddo

      end

