c   p(12GeV/c)+p  p(24GeV/c)+p
c
      include 'jam1.inc'
      include 'jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15

c....Initialize JAM.
      fname(1)='jam.cfg'  ! input file name.
      mevent=10
      bmin=0.0D0
      bmax=-1.0D0
      dt=100.0D0
      nstep=1
      win=200.D0
      frame='nn      '
      targ='32S     '
      proj='32S     '
      cwin='200gevc        '
      mstc(8)=1
      mstc(74)=0   ! dipole-approximated QCD radiation of the string
      mstc(81)=0  ! 1:hard scattering on
      call jaminit(mevent,bmin,bmax,dt,nstep,
     $                             frame,proj,targ,cwin)
      nevent=mstc(2)

c...Initialze event analysis
      call anal1

c...Loop over simulation
c==================================================
      do iev=1,nevent

c....Generate one event
        call jamevt(iev)

        if(mod(iev,100).eq.0) write(6,*)'event=',iev

c...Event analysis
        call anal2

c....End simulation
      end do
c==================================================

c...Finsih event
      call jamfin

c...Output event analysis
      call anal3

      end
c
c  Exp. data
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  exp. data from
c  V.Blobel, et.al.,  N.P.B69(1974) 454-492
c
c total inel cross section 12GeV/c=  29.75mb, 24GeV/c=  30.60mb
c
c inelastic charged particle multiplicity
c  12GeV/c: 3.43 +- 0.03
c  24GeV/c: 4.25 +- 0.03
c
c Integraed single inclusive cross sections per inelastic collision
c  note: Lambda and lambda-bar cross sections contain
c        lambda(lambda-bar) from sigma0(sigma-bar0) decays.
c
c pp 12GeV/c
c sig(mb) error   multiplicity   particle
c 21.100  +-0.4     .7092        pi-     
c 35.200  +-2.4    1.1832        pi0     
c 42.700  +-0.7    1.4353        pi+       pion total:3.3277
c  1.150  +-0.03    .0387        k0/a-k0 
c 37.500  +-0.6    1.2605        p       
c  1.120  +-0.03    .0376        lambda  
c   .003  +0.001    .0001        a-lambda
c         -0.002
c   .160  +-0.01    .0054        sigma-  
c   .490  +-0.02    .0165        sigma+  
c
c pp 24GeV/c
c sig(mb) error  multiplicity    particle
c 33.800  +-0.6    1.1046        pi-     
c 53.500  +-3.1    1.7484        pi0     
c 56.800  +-0.9    1.8562        pi+      pion total: 4.7092 
c  2.510  +-0.06    .0820        k0/a-k0 
c 37.900  +-0.6    1.2386        p       
c  1.760  +-0.06    .0575        lambda  
c   .021  +0.004    .0007        a-lambda
c         -0.010
c   .280  +-0.02    .0092        sigma-  
c   .850  +-0.03    .0278        sigma+  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c***********************************************************************

      subroutine anal1

c...Initialize all histgram.
      include 'jam1.inc'
      include 'jam2.inc'
      character chau*16
      save n_charge, n_pim, n_pi0, n_pip, n_kaon, n_kaon0,
     $  n_kaonp, n_kaonm, n_lamb, n_sigm, n_sig0, n_sigp,
     $ n_alamb, n_asigm, n_asig0, n_asigp, n_prot, n_aprot,
     $ n_100, n_mes, n_neg
      save wei,ylab,yproj
      save wp,wy

c....RHIC
c     data ymin,ymax,wy,pmin,pmax,wp/-10.0,10.0,0.5,0.0,9.0,0.5/
c....AGS/SPS
c     data ymin,ymax,wy,pmin,pmax,wp/-5.0,8.0,0.3,0.0,2.0,0.05/
c...Number of pions
      data n_charge,n_pim,n_pi0,n_pip,
     $ n_kaon,n_kaon0,n_kaonp,n_kaonm,
     $ n_lamb,n_sigm,n_sig0,n_sigp,
     $ n_alamb,n_asigm,n_asig0,n_asigp,n_prot,n_aprot,
     $ n_100, n_mes, n_neg/21*0/


c...Event weight
      wei=1.D0/dble(mstc(2))
      ylab=pard(17)

c...Kinematic bounds on rapidity.
      dbetpr=pard(35)
      dbetta=pard(45)
      yproj=0.5D0*log((1.0d0+dbetpr)/(1.0d0-dbetpr))
      ytarg=0.5D0*log((1.0d0+dbetta)/(1.0d0-dbetta))
      ymin=-yproj*1.6D0
      ymax=yproj*1.6D0
      nymx=30
      nymx1=nymx
      if(pard(16).le.40.D0) then
        pmin=0.0D0
        pmax=3.0D0
        npmx=30
      else
        pmin=0.0D0
        pmax=20.0D0
        npmx=30
      endif
      if(mstc(4).eq.0) then
      else if(mstc(4).eq.3) then
      else
       ymax=ymax+ylab
      endif
      wy=(ymax-ymin)/nymx
      wy1=(ymax-ymin)/nymx1
      wp=(pmax-pmin)/npmx

c     if(pard(16).le.40.) then
c       ymin=-5.0
c       ymax=8.0
c       wy1=0.3
c       wy=0.3
c       pmin=0.0
c       pmax=2.0
c       wp=0.05
c     else
c       ymin=-10.0
c       ymax=10.0
c       wy1=0.75
c       wy=0.5
c       pmin=0.0
c       pmax=9.0
c       wp=0.5
c     endif
c     nymx=(ymax-ymin)/wy
c     nymx1=(ymax-ymin)/wy1
c     npmx=(pmax-pmin)/wp

c...Inititalize histogram booking.

      if(nymx.le.0.or.nymx.gt.100) then
        write(6,*)'range error nymx',nymx
      endif

c...Rapidity distributions.
      call vbook1(10,'dN/dy - pi- k- pbar',nymx,ymin,ymax)
      call vbook1(11,'dN/dy - proton',nymx1,ymin,ymax)
      call vbook1(12,'dN/dy - pion- ',nymx,ymin,ymax)
      call vbook1(13,'dN/dy - pion+ ',nymx,ymin,ymax)
      call vbook1(14,'dN/dy - charged',nymx,ymin,ymax)
      call vbook1(15,'dN/dy - p-p-bar',nymx1,ymin,ymax)
      call vbook1(16,'dN/dy - h-',nymx,ymin,ymax)
      call vbook1(17,'dN/dy - net baryon',nymx1,ymin,ymax)
      call vbook1(18,'dN/dy - k+',nymx,ymin,ymax)
      call vbook1(19,'dN/dy - k-',nymx,ymin,ymax)
      call vbook1(20,'dN/dy - lambda',nymx,ymin,ymax)

      call vbook1(21,'dN/dpt**2 - proton ',npmx,pmin,pmax)
      call vbook1(22,'dN/dpt**2 - pion-  ',npmx,pmin,pmax)
      call vbook1(23,'dN/dpt**2 - pion+  ',npmx,pmin,pmax)
      call vbook1(24,'dN/dpt**2 - charged',npmx,pmin,pmax)

c...CERN/SPS for S+S NA35 P.R.L. 70(1994)1419.
      call vbook1(25,'1/ptdN/dpt - 0.5D0<y<3.0D0 proton',npmx,pmin,pmax)
      call vbook1(26,'1/ptdN/dpt - 0.8D0<y<2.0D0 h-',npmx,pmin,pmax)
      call vbook1(27,'1/ptdN/dpt - 2.0D0<y<3.0D0 h-',npmx,pmin,pmax)
      call vbook1(28,'1/ptdN/dpt - 3.0D0<y<4.0D0 h-',npmx,pmin,pmax)
      call vbook1(29,'1/ptdN/dpt - 4.2D0<y<4.4D0 h-',npmx,pmin,pmax)

c...dn/dm_t CERN/SPS transverse mass spectra of negative..
      call vbook1(101,'dN/m*dmt y=3.4D0 h-(pi k pbar)',npmx,pmin,pmax)
      call vbook1(102,'dN/m*dmt y=3.9D0 h-(pi k pbar)',npmx,pmin,pmax)
      call vbook1(103,'dN/m*dmt y=4.4D0 h-(pi k pbar)',npmx,pmin,pmax)
      call vbook1(104,'dN/m*dmt y=4.9D0 h-(pi k pbar)',npmx,pmin,pmax)
      call vbook1(105,'dN/m*dmt y=5.4D0 h-(pi k pbar)',npmx,pmin,pmax)

c...dn/dm_t CERN/SPS transverse mass spectra of p-pbar. central 5%
      call vbook1(111,'dN/m*dmt y=2.9D0 p-pbar',npmx,pmin,pmax)
      call vbook1(112,'dN/m*dmt y=3.4D0 p-pbar',npmx,pmin,pmax)
      call vbook1(113,'dN/m*dmt y=3.9D0 p-pbar',npmx,pmin,pmax)
      call vbook1(114,'dN/m*dmt y=4.4D0 p-pbar',npmx,pmin,pmax)
      call vbook1(115,'dN/m*dmt y=4.9D0 p-pbar',npmx,pmin,pmax)
      call vbook1(116,'dN/m*dmt y=5.4D0 p-pbar',npmx,pmin,pmax)
c...dn/mdm proton pion
      call vbook1(121,'dN/m*dmt 2.5D0<y<3.3D0 p',npmx,pmin,pmax)
      call vbook1(122,'dN/m*dmt 2.5D0<y<3.3D0 pbar',npmx,pmin,pmax)
      call vbook1(123,'dN/m*dmt 2.5D0<y<3.3D0 pi-',npmx,pmin,pmax)
      call vbook1(124,'dN/m*dmt 2.5D0<y<3.3D0 pi+',npmx,pmin,pmax)


c...Transverse flow
      call vbook1(31,'<p_x> nucleon',nymx1,ymin,ymax)
      call vbook1(32,'<p_x> pion',nymx,ymin,ymax)
      call vbook1(33,'<p_t> nucleon',nymx1,ymin,ymax)
      call vbook1(34,'<p_t> pion',nymx,ymin,ymax)

      call vbook1(41,'nucl dN/dy',nymx1,ymin,ymax)
      call vbook1(42,'pion dN/dy',nymx,ymin,ymax)

      call vbook1(51,'nucleon <p_x>',nymx1,ymin,ymax)
      call vbook1(52,'pion <p_x>',nymx,ymin,ymax)
      call vbook1(53,'nucleon <p_t>',nymx1,ymin,ymax)
      call vbook1(54,'pion <p_t>',nymx,ymin,ymax)

      call vbook1(55,'nucleon <p_x>-error',nymx1,ymin,ymax)
      call vbook1(56,'pion <p_x>-error',nymx,ymin,ymax)
      call vbook1(57,'nucleon <p_t>-error',nymx1,ymin,ymax)
      call vbook1(58,'pion <p_t>-error',nymx,ymin,ymax)

c...Nuclear density
c     call vbook1(71,'drho/dr -proj. density',100,0.0,12.)

      return

c***********************************************************************
      entry anal2
c...Anaysis of this event 

        nch=0
        npa=0
        do i=1,nv

         npa=npa+1


c...Select charged particle
         kch=jamchge(k(2,i))
         if(kch.ne.0) then
              nch=nch+1
              if(kch.lt.0) n_neg=n_neg+1
               n_charge=n_charge+1
         endif

c          rap=0.5*log( max(e(i)+p(3,i),1.e-8)/max(e(i)-p(3,i),1.e-8) )
           D1=sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
           rap=0.5D0*log( max(D1+p(3,i),1.D-8)/max(D1-p(3,i),1.D-8) )
           y=rap

          if(mstc(4).eq.0) then
          else if(mstc(4).eq.3) then
          else
           rap=rap+ylab
          endif

           ptsq=p(1,i)**2+p(2,i)**2
           pt=sqrt(ptsq)
           pt=max(pt,1.D-8)
           emt0=sqrt(p(5,i)**2+ptsq)
           kf=k(2,i)
           kc=jamcomp(kf)
           if(kc.le.0.or.kc.gt.mstu(6)) then
              write(6,*)'Invalid code i kf kc',i,kf,kc,nv,nbary,nmeson
              goto 3000
           endif

c          id=kchg(kc,5)
           ibar=kchg(kc,6)

c....Charged particles.
           if(kch.ne.0) then
             call vfill1(14,rap,wei/wy) 
             call vfill1(24,pt,wei/wp/pt/2) 

c.......Negative charged particles.
             if(kch.lt.0) then
               call vfill1(16,rap,wei/wy)
c...pi-,k-,pbar
           if(kf.eq.-211.or.kf.eq.-321.or.kf.eq.-2212) then
             call vfill1(10,rap,wei/wy) 
c....S+S CERN/SPS
             if(rap.lt.0.8D0) then
             else if(rap.le.2.0D0) then
               call vfill1(26,pt,wei/(wp*pt))
             else if(rap.ge.2.0D0.and.rap.le.3.0D0) then
               call vfill1(27,pt,wei/(wp*pt))
             else if(rap.ge.3.0D0.and.rap.le.4.0D0) then
               call vfill1(28,pt,wei/(wp*pt))
             else if(rap.ge.4.2D0.and.rap.le.4.4D0) then
               call vfill1(29,pt,wei/(wp*pt))
             endif

c.....Pb+Pb CERN/SPS h-  5%
             emt=emt0-0.138D0
	     if(rap.lt.3.15D0) then
             else if(rap.ge.3.15D0.and.rap.le.3.65D0) then
               call vfill1(101,emt,wei/(emt0*wp)) 
             else if(rap.le.4.15D0) then
               call vfill1(102,emt,wei/(emt0*wp)) 
             else if(rap.le.4.65D0) then
               call vfill1(103,emt,wei/(emt0*wp)) 
             else if(rap.le.5.15D0) then
               call vfill1(104,emt,wei/(emt0*wp)) 
             else if(rap.le.5.65D0) then
               call vfill1(105,emt,wei/(emt0*wp)) 
             endif

           else
c            write(99,*)kf
           endif

             endif

           endif



c...Net baryons.
           if(ibar.eq.3) then
             if(kf.gt.0) then
               call vfill1(17,rap,wei/wy1) 
             else if(kf.lt.0) then
               call vfill1(17,rap,-wei/wy1) 
             endif
           endif

c.......Protons.
           if(abs(kf).eq.2212) then
             emt=emt0-p(5,i)
             if(kf.eq.2212) then
               n_prot=n_prot+1
               call vfill1(15,rap,wei/wy1) 
               call vfill1(11,rap,wei/wy1) 
               call vfill1(21,pt,wei/wp/pt/2) 


c.............Pb+Pb
c              if(rap.ge.2.5.and.rap.le.3.3) then
c                call vfill1(121,emt,wei/(wp*emt0)) 
c              endif

	       if(rap.lt.2.65D0) then
               else if(rap.ge.2.65D0.and.rap.le.3.15D0) then !y=2.9D0
                 call vfill1(111,emt,wei/(wp*emt0)) 
               else if(rap.le.3.65D0) then            !y=3.4D0
                 call vfill1(112,emt,wei/(wp*emt0)) 
               else if(rap.le.4.15D0) then            !y=3.9D0
                 call vfill1(113,emt,wei/(wp*emt0)) 
               else if(rap.le.4.65D0) then            !y=4.4D0
                 call vfill1(114,emt,wei/(wp*emt0)) 
               else if(rap.le.5.15D0) then            !y=4.9D0
                 call vfill1(115,emt,wei/(wp*emt0)) 
               else if(rap.le.5.65D0) then            !y=5.4D0
                 call vfill1(116,emt,wei/(wp*emt0)) 
               endif

c.............S+S
               if(rap.ge.0.5D0.and.rap.le.3.0D0) then
                 call vfill1(25,pt,wei/(wp*pt)) 
               endif

             else if(kf.eq.-2212) then
               n_aprot=n_aprot+1
               call vfill1(15,rap,-wei/wy1) 
c.............S+S
               if(rap.ge.0.5D0.and.rap.le.3.0D0) then
                 call vfill1(25,pt,-wei/(wp*pt)) 
               endif

c              if(rap.ge.2.5.and.rap.le.3.3) then
c                call vfill1(122,emt,-wei/(wp*wp*emt0)) 
c              endif

               if(rap.lt.2.65D0) then
               else if(rap.ge.2.65D0.and.rap.le.3.15D0) then !y=2.9D0
                 call vfill1(111,emt,-wei/(wp*emt0)) 
               else if(rap.le.3.65D0) then            !y=3.4D0
                 call vfill1(112,emt,-wei/(wp*emt0)) 
               else if(rap.le.4.15D0) then            !y=3.9D0
                 call vfill1(113,emt,-wei/(wp*emt0)) 
               else if(rap.le.4.65D0) then            !y=4.4D0
                 call vfill1(114,emt,-wei/(wp*emt0)) 
               else if(rap.le.5.15D0) then            !y=4.9D0
                 call vfill1(115,emt,-wei/(wp*emt0)) 
               else if(rap.le.5.65D0) then            !y=5.4D0
                 call vfill1(116,emt,-wei/(wp*emt0)) 
               endif

             endif

c.......Pions.
           else if(kf.eq.-211) then
               call vfill1(12,rap,wei/wy) 
               call vfill1(22,pt,wei/wp/pt/2) 

c              if(rap.ge.2.5.and.rap.le.3.3) then
c                call vfill1(123,emt,wei/(wy*wp*emt0)) 
c              endif

               n_pim=n_pim+1
           else if(kf.eq.211) then
               call vfill1(13,rap,wei/wy) 
               call vfill1(23,pt,wei/wp/pt/2) 

c              if(rap.ge.2.5.and.rap.le.3.3) then
c                call vfill1(124,emt,wei/(wp*emt0)) 
c              endif

               n_pip=n_pip+1
           else if(kf.eq.111) then
               n_pi0=n_pi0+1
           else if(abs(kf).eq.311.or.abs(kf).eq.321) then
               if(kf.eq.321) then
                 n_kaonp=n_kaonp+1
                 call vfill1(18,rap,wei/wy) 
                endif
               if(kf.eq.-321) then
                 n_kaonm=n_kaonm+1
                 call vfill1(19,rap,wei/wy) 
               endif
               if(abs(kf).eq.311) n_kaon0=n_kaon0+1
               n_kaon=n_kaon+1
           else if(kf.eq.3122) then
               n_lamb=n_lamb+1
               call vfill1(20,rap,wei/wy) 

           else if(kf.eq.3112) then
               n_sigm=n_sigm+1
           else if(kf.eq.3212) then
               n_sig0=n_sig0+1
           else if(kf.eq.3222) then
               n_sigp=n_sigp+1

           else if(kf.eq.-3122) then
               n_alamb=n_alamb+1
           else if(kf.eq.-3112) then
               n_asigm=n_asigm+1
           else if(kf.eq.-3212) then
               n_asig0=n_asig0+1
           else if(kf.eq.-3222) then
               n_asigp=n_asigp+1

           else if(kc.le.100) then
               n_100=n_100+1
           else if(ibar.eq.0) then
c              call luname(kf,chau)
c              write(71,*)k(1,i),k(2,i),p(5,i),' ',chau
               n_mes=n_mes+1
           endif

c.....Transverse flow of nucleons and pions
           if(kf.eq.2112.or.kf.eq.2212) then
             call vfill1(31,y,p(1,i)) 
             call vfill1(33,y,pt) 
             call vfill1(55,y,p(1,i)**2) 
             call vfill1(57,y,pt**2) 
             call vfill1(41,y,1.0D0) 
           else if(kf.eq.111.or.abs(kf).eq.211) then
             call vfill1(32,y,p(1,i)) 
             call vfill1(34,y,pt) 
             call vfill1(56,y,p(1,i)**2) 
             call vfill1(58,y,pt**2) 
             call vfill1(42,y,1.0D0) 
          endif

3000    end do

      return

c...Ouptput.
c***********************************************************************
      entry anal3

      open(70,file='file70',status='unknown')
      write(70,*)'ylab yproj=',ylab,yproj,ytarg
      write(70,*)'charged',n_charge*wei
      write(70,*)'negative',n_neg*wei
      write(70,*)'pi- pi0 pi+',n_pim*wei,n_pi0*wei,n_pip*wei
      write(70,*)'pion  total',(n_pim+n_pi0+n_pip)*wei
      write(70,*)'proton total',n_prot*wei
      write(70,*)'a-proton total',n_aprot*wei
      write(70,*)'lambda a-lam total',n_lamb*wei,n_alamb*wei
      write(70,*)'sima- a-sigma- total',n_sigm*wei,n_asigm*wei
      write(70,*)'sima0 a-sigma0 total',n_sig0*wei,n_asig0*wei
      write(70,*)'sima+ a-sigma+ total',n_sigp*wei,n_asigp*wei
      sigtot=(n_sigm+n_sig0+n_sigp)*wei
      asigtot=(n_asigm+n_asig0+n_asigp)*wei
      write(70,*)'sigma a-sigma total',sigtot,asigtot

      write(70,*)'kaon+  total',n_kaonp*wei
      write(70,*)'kaon-  total',n_kaonm*wei
      write(70,*)'kaon  total',n_kaon*wei
      write(70,*)'kaon0 total',n_kaon0*wei
      write(70,*)'other mesons',n_mes*wei
      write(70,*)'lepton gamma etc',n_100*wei
      write(70,*)'average number of jet',mstd(55)*wei
      close(70)

c...Output of histograms.
      fac=1.D0
      mnorm=0
      mform=0
      do i=0,10
c      call vscale(10+i,fac)
       call vprint(10+i,mnorm,mform)
      end do

c...Transverse flow of nucl. and pions.
      call vopera(31,'/',41,51,1.0D0,1.0D0)
      call vopera(32,'/',42,52,1.0D0,1.0D0)
      call vopera(33,'/',41,53,1.0D0,1.0D0)
      call vopera(34,'/',42,54,1.0D0,1.0D0)

      call vopera(41,'M',31,55,1.0D0,1.0D0)
      call vopera(42,'M',32,56,1.0D0,1.0D0)
      call vopera(41,'M',33,57,1.0D0,1.0D0)
      call vopera(42,'M',34,58,1.0D0,1.0D0)
      do i=1,8
       call vprint(50+i,mnorm,mform)
      end do

      mform=1
      do i=1,9
       call vprint(20+i,mnorm,mform)
      end do

c...dn/dm_t CERN/SPS transverse mass spectra of p-pbar.
      do i=1,6
       call vprint(110+i,0,1)
      end do

c...dn/dm_t CERN/SPS transverse mass spectra of negative..
      do i=1,5
       call vprint(100+i,0,1)
      end do

c...dn/mdm proton pion
      do i=1,4
       call vprint(120+i,0,1)
      end do

      end
