c...A main program for p + p

      include 'jam1.inc'
      include 'jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15

c....Initialize JAM.
      fname(1)='jampp.cfg'  ! input file name.
      mstc(17)=1          ! forgid elastic collision.
      call jaminit(mevent,bmin,bmax,dt,nstep,frame,proj,targ,cwin)
      nevent=mstc(2)

c...Initialize analysis.
      call anal1
      call anal11

c...Simulation start.
      do iev=1,nevent

c...Simulate one event.
        call jamevt(iev)

        if(mod(iev,1000).eq.0) write(6,*)'event=',iev

c...Data analysis.
        call anal2
        call anal21

      end do

c...Final output.
      call jamfin

c...Print analysis results.
      call anal3
      call anal31

      end

c***********************************************************************

      subroutine anal1

      include 'jam1.inc'
      include 'jam2.inc'
      dimension npa(0:20)
      save wy,wp
      save ylab,yproj,ytarg
      save npa

c...npa(1) : average charged
c...npa(2) : average negative
c...npa(3) : pion-
c...npa(4) : pion0
c...npa(5) : pion+
c...npa(6) : antikaon
c...npa(7) : Kaon
c...npa(8) : lambda
c...npa(9) : a-lambda
c...npa(10): sigma
c...npa(11): a-sigma
c...npa(12): proton
c...npa(13): a-proton
c...npa(14):

      do i=0,20
       npa(i)=0
      end do

c....Rapidity distribution.
      ylab=pard(17)
      dbetpr=pard(35)
      dbetta=pard(45)
      yproj=0.5d0*log((1.0d0+dbetpr)/(1.0d0-dbetpr))
      ytarg=0.5d0*log((1.0d0+dbetta)/(1.0d0-dbetta))
c     ymin=ytarg*1.5D0
c     ymax=yproj*1.5D0
      ymin=-2.2d0
      ymax=2.2d0
      wy=0.1d0
      nymx=(ymax-ymin)/wy
      if(mstc(4).eq.0) then
      else if(mstc(4).eq.3) then
      else
       ymax=ymax+ylab
       ymin=ymin+ylab
      endif

      call vbook1(10,'dN/dy - pi- k- pbar',nymx,ymin,ymax)
      call vbook1(11,'dN/dy - proton',nymx,ymin,ymax)
      call vbook1(12,'dN/dy - pion- ',nymx,ymin,ymax)
      call vbook1(13,'dN/dy - pion+ ',nymx,ymin,ymax)
      call vbook1(14,'dN/dy - charged',nymx,ymin,ymax)
      call vbook1(15,'dN/dy - p-p-bar',nymx,ymin,ymax)
      call vbook1(16,'dN/dy - h-',nymx,ymin,ymax)
      call vbook1(17,'dN/dy - net baryon',nymx,ymin,ymax)

c...Transverse distributions.
      if(pard(16).le.40.D0) then
        pmin=0.0D0
        pmax=3.0D0
        npmx=30
      else
        pmin=0.0D0
        pmax=20.0D0
        npmx=30
      endif
      wp=(pmax-pmin)/npmx
      call vbook1(21,'dN/dpt**2 - proton ',npmx,pmin,pmax)
      call vbook1(22,'dN/dpt**2 - pion-  ',npmx,pmin,pmax)
      call vbook1(23,'dN/dpt**2 - pion+  ',npmx,pmin,pmax)
      call vbook1(24,'dN/dpt**2 - charged',npmx,pmin,pmax)

      call vbook1(31,'dN/d(eta) - charged',nymx,ymin,ymax)
      call vbook1(32,'P(n) - charged',50,0.D0,100.D0)

      return

c***********************************************************************

      entry anal2

      nch=0
      neg=0
c...Loop over all particles.
      do i=1,nv

        kf=k(2,i)
        kc=jamcomp(kf)
        if(kc.le.0.or.kc.gt.mstu(6)) then
           write(6,*)'Invalid code i kf kc',i,kf,kc,nv,nbary,nmeson
           goto 3000
        endif

        rap=0.5D0*log( max(p(4,i)+p(3,i),1.D-8)/max(p(4,i)-p(3,i), 
     & 1.D-8) )

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

        npa(0)=npa(0)+1

c...Charged particles.
        kch=jamchge(k(2,i))
        if(kch.ne.0) then
          nch=nch+1
          npa(1)=npa(1)+1
          call vfill1(14,rap,1.D0/wy)
          call vfill1(31,eta,1.D0/wy)
          call vfill1(24,pt,1.D0/wp/pt/2) 

c...Negative charged particles.
          if(kch.lt.0) then
            neg=neg+1
            npa(2)=npa(2)+1
            call vfill1(16,rap,1.D0/wy)
          endif
        endif

c...Net baryons.
        ibar=kchg(kc,6)
        if(ibar.eq.3) then
          if(kf.gt.0) then
            call vfill1(17,rap,1.D0/wy) 
          else if(kf.lt.0) then
            call vfill1(17,rap,-1.D0/wy) 
          endif
        endif

c...h-(pi-,k-,p~)
        if(kf.eq.-211.or.kf.eq.-321.or.kf.eq.-2212) then
           call vfill1(10,rap,1.D0/wy)
        endif


c.......Protons.
        if(abs(kf).eq.2212) then
          if(kf.eq.2212) then
            npa(12)=npa(12)+1
            call vfill1(15,rap,1.D0/wy) 
            call vfill1(11,rap,1.D0/wy) 
            call vfill1(21,pt,1.D0/wp/pt/2) 


          else if(kf.eq.-2212) then
            npa(13)=npa(13)+1
            call vfill1(15,rap,-1.D0/wy) 
          endif

c.......Pions.
        else if(kf.eq.-211) then
          npa(3)=npa(3)+1
          call vfill1(12,rap,1.D0/wy) 
          call vfill1(22,pt,1.D0/wp/pt/2)

        else if(kf.eq.111) then
          npa(4)=npa(4)+1
        else if(kf.eq.211) then

          npa(5)=npa(5)+1
          call vfill1(13,rap,1.D0/wy) 
          call vfill1(23,pt,1.D0/wp/pt/2)

        else if(kf.eq.321) then
        endif

        if(kf.eq.3122) then
          npa(8)=npa(8)+1
        else if(kf.eq.-3122) then
          npa(9)=npa(9)+1

        else if(kf.eq.3112) then ! Sigma-
          npa(10)=npa(10)+1
        else if(kf.eq.3212) then ! Sigma0
          npa(10)=npa(10)+1
        else if(kf.eq.3222) then ! Sigma+
          npa(10)=npa(10)+1

        else if(kf.eq.-3112) then ! a-Sigma-
          npa(11)=npa(11)+1
        else if(kf.eq.-3212) then ! a-Sigma0
          npa(11)=npa(11)+1
        else if(kf.eq.-3222) then ! a-Sigma+
          npa(11)=npa(11)+1

        else if(kf.eq.321.or.kf.eq.311) then
          npa(7)=npa(7)+1
        else if(kf.eq.-321.or.kf.eq.-311) then
          npa(6)=npa(6)+1
        endif


3000  end do
      call vfill1(32,dble(nch),1.D0) 

      return

c***********************************************************************

      entry anal3

c...Output of histograms.

c...Event weight
      fac=1.D0/dble(mstc(2))

c...Rapidity distributions.
      do i=0,7
       call vscale(10+i,fac)
       call vprint(10+i,0,0)
      end do

c...Pt distributions.
      do i=1,4
       call vscale(20+i,fac)
       call vprint(20+i,0,1)
      end do

c...dn/d(eta),P(n)
      do i=1,2
       call vscale(30+i,fac)
       call vprint(30+i,0,1)
      end do

      open(70,file='file70',status='unknown')
      write(70,*)'ylab yproj ytarg=',ylab,yproj,ytarg
      write(70,*)'average mult',npa(0)*fac
      write(70,*)'charged',npa(1)*fac
      write(70,*)'negative',npa(2)*fac
      write(70,*)'pi- pi0 pi+',npa(3)*fac,npa(4)*fac,npa(5)*fac
      write(70,*)'pion  total',(npa(3)+npa(4)+npa(5))*fac
      write(70,*)'proton total',npa(12)*fac
      write(70,*)'a-proton total',npa(13)*fac
      write(70,*)'lambda a-lam total',npa(8)*fac,npa(9)*fac
      write(70,*)'sima a-sigma total',npa(10)*fac,npa(11)*fac

      write(70,*)'kaon   total',npa(7)*fac
      write(70,*)'akaon  total',npa(6)*fac
      write(70,*)'average number of jet',pard(87)*fac

      end

c***********************************************************************

      subroutine anal11

c...Initialize analysis for overall run.

      include 'jam1.inc'
      include 'jam2.inc'
      parameter(ncdet=16,npdet=14) ! stable particles
      character cfile(ncdet)*15,cfile3(npdet)*15
      dimension ncount(0:ncdet),ncount3(npdet)
      dimension ntyp(npdet), mchan(npdet,ncdet), ktyp(npdet)

      save ncount,ncount3

      data cfile/'elastic'
     &            ,'p n pi+',    'p p pion0'
     &            ,'p n pi^+ pi^0', 'p p pi+ pi-'
     &            ,'pn 2pi+ pi-',  'pp pi+ pi0 pi-'
     &            ,'pp 2pi+ 2pi-', 'pn 3pi+ 2pi-'
     &            ,'LpK+',   'S+nK+'
     &            ,'S0pK+',  'Ln+K+'
     &            ,'Lp0K+',  'pp eta'
     $            ,'Lambda p pi+ K0'/

      data cfile3/'p + X','pp-n X','Lambda + X',
     $              'S^- +  X','S^0 + X.dat','S^+ + X',
     $              'pi^- X','pi^0 + X','pi^+ + X',
     $              'K^- X','Kb0 + X','K0 + X',
     $              'K^+ X','eta  + X'/

c........... n  p  L  S- S0 S+ p- p0 p+ K- Kb K0 K+ eta
      data mchan
     1      /0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   ! pp (elastic)
     2       1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   ! np pi+ 1
     3       0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,   ! pp pi0 1
     4       1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,   ! np pi0 pi+ 2
     5       0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0,   ! pp pi- pi+ 2
     6       1, 1, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0,   ! np pi- 2pi+
     7       0, 2, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0,   ! pp pi- pi0 pi+
     8       0, 2, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0,   ! pp 2pi- 2pi+
     9       1, 1, 0, 0, 0, 0, 2, 0, 3, 0, 0, 0, 0, 0,   ! np 2pi- 3pi+
     $       0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,   ! pL  K+
     1       1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,   ! nS+ K+
     2       0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0,   ! pS0 K+
     3       1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0,   ! nL  pi+ K+
     4       0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,   ! pL pi0 K+
     5       0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,   ! pp eta
     6       0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0 /  ! Lp pi+ K0
      data ktyp
     &      /2112,2212,3122,3112,3212,3222,
     $              -211,111,211,-321,-311,311,321,221/


      do i=1,ncdet
        ncount(i)=0
      enddo
      do i=1,npdet
        ncount3(i)=0
      enddo

      return

c***********************************************************************

      entry anal21

c...Count data.

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
          if(ntyp(ipdet).ne.mchan(ipdet,icdet)) goto 10
        enddo
        ichanel=icdet
        goto 4000
10      continue
      enddo

4000  continue
      if(notdef.eq.0.and.(ichanel.ge.0.and.ichanel.le.ncdet))
     $   ncount(ichanel)=ncount(ichanel)+1

c     if(ichanel.eq.1.or.ichanel.eq.3) then
c       write(mstc(38),*)'itag=',itag,'notdef=',notdef,ichanel
c       write(mstc(38),*) (ntyp(i),i=1,npdet)
c       call jamlist(1)
c     endif

      return

c***********************************************************************

      entry anal31

c...Output results.

c...Event weight
      wei=1.d0/dble(mstc(2))
      fac=parc(4)**2*paru(1)*10*wei
      srt=pard(16)

      write(70,800)pard(15),srt
c...Exclusive data
      do i=1,ncdet
        write(70,810)cfile(i),fac*ncount(i),ncount(i)*wei
      enddo

c...Inclusive data
      do i=1,npdet
        write(70,810)cfile3(i),fac*ncount3(i),ncount3(i)*wei
      enddo

800   format('Beam energy C.M.eng=',2(f8.3,1x))
810   format(a15,2x,f12.7,1x,f9.4)

      end
