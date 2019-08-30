c   p(12GeV/c)+p  p(24GeV/c)+p
c
      include 'jam1.inc'
      include 'jam2.inc'
      character chau*16
      character frame*8,proj*8,targ*8,cwin*15

c....Initialize JAM.
      fname(1)='jam.cfg'
      call jaminit(nev,b1,b2,dt,nstp,frame,proj,targ,cwin)
      nevent=mstc(2)

c...event weight
      wei=1.D0/dble(nevent)

      if(mstd(2).eq.1.and.mstd(5).eq.1) then
        ylab=0.0D0
      else
        ylab=pard(17)
      endif

c...Inititalize histogram booking.
      ymin=-4.5D0
      ymax=4.5D0
      wy1=0.1D0
      wy=0.1D0
      pmin=0.0D0
      pmax=3.0D0
      wp=0.05D0
      nymx=(ymax-ymin)/wy
      nymx1=(ymax-ymin)/wy1
      npmx=(pmax-pmin)/wp
c...String p_z
      pzmin=-2.D0
      pzmax=2.D0
      npzmx=(pzmax-pzmin)/wp

      call vbook1(31,'dN/dy - proton-proj',nymx,ymin,ymax)
      call vbook1(32,'dN/dy - proton-targ',nymx,ymin,ymax)

c...String Mass dist.
      rmin=0.8D0
      rmax=6.0D0
      wr=0.1D0
      nrmx=(rmax-rmin)/wr
      call vbook1(33,'dN/dm - proton-proj',nrmx,rmin,rmax)
      call vbook1(34,'dN/dm - proton-targ',nrmx,rmin,rmax)

c...String pt
      call vbook1(35,'dN/dp_t - proj',npmx,pmin,pmax)
      call vbook1(36,'dN/dp_t - targ',npmx,pmin,pmax)
      call vbook1(37,'dN/dp_z - proj',npzmx,pzmin,pzmax)
      call vbook1(38,'dN/dp_z - targ',npzmx,pzmin,pzmax)
      call vbook1(39,'dN/dm - mass',nrmx,rmin,rmax)
      call vbook1(40,'dN/dp_z - emnuc',npzmx,pzmin,pzmax)
      call vbook1(41,'dN/dp_t - prot',npmx,pmin,pmax)


      call vbook1(51,'dN/dy - p1 ',nymx,ymin,ymax)
      call vbook1(52,'dN/dy - p2 ',nymx,ymin,ymax)

c...Rapidity distributions.
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

c...Number of pions
      n_charge=0
      n_pim=0
      n_pi0=0
      n_pip=0
      n_kaon=0
      n_kaon0=0
      n_kaonp=0
      n_kaonm=0
      n_lamb=0
      n_sigm=0
      n_sig0=0
      n_sigp=0
      n_alamb=0
      n_asigm=0
      n_asig0=0
      n_asigp=0
      n_prot=0
      n_aprot=0
      n_100=0
      n_mes=0
      n_neg=0

c==================================================
      do iev=1,nevent

        call jamevt(iev)
        if(mod(iev,100).eq.0) write(6,*)'event=',iev

c       if(mod(iev,50).eq.0) then
c         mstc38=mstc(38)
c         mstc(38)=20
c         write(20,*)iev
c         call jamlist(1)
c         mstc(38)=mstc38
c       endif

          em1=pare(75)
          px1=pare(71)
          py1=pare(72)
          pz1=pare(73)
          em2=pare(99)
          px2=pare(96)
          py2=pare(97)
          pz2=pare(98)

c....Mass distibution
          call vfill1(33,em1,wei/wr) 
          call vfill1(34,em2,wei/wr) 
          if(abs(pz1).gt.2.8D0) then
            call vfill1(39,em1,wei/wr) 
          endif
          if(em1.le.1.0D0) then
              pt1=sqrt(px1**2+py1**2)
              call vfill1(40,pz1,wei/wp)
              call vfill1(41,pt1,wei/wp)
          endif
          if(em2.le.1.0D0) then
             pt2=sqrt(px2**2+py2**2)
             call vfill1(40,pz2,wei/wp)
             call vfill1(41,pt2,wei/wp)
          endif

c...Momentum distribution
c...p_t
          call vfill1(35,sqrt(pare(71)**2+pare(72)**2),wei/wp) 
          call vfill1(36,sqrt(pare(96)**2+pare(97)**2),wei/wp) 
c...p_z
          call vfill1(37,pare(73),wei/wp) 
          call vfill1(38,pare(98),wei/wp) 

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

         e1=sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
         rap=0.5D0*log( max(e1+p(3,i),1.D-8)/max(e1-p(3,i),1.D-8) )

          if(mstc(4).eq.0) then
          else if(mstc(4).eq.3) then
          else
           rap=rap+ylab
          endif

           pt=sqrt(p(1,i)**2+p(2,i)**2)
           pt=max(pt,1.D-8)

           if(kch.ne.0) then
             call vfill1(14,rap,wei/wy) 
             if(kch.le.0) call vfill1(16,rap,wei/wy)
             call vfill1(24,pt,wei/wp/pt/2) 
           endif

           kf=k(2,i)
           kc=jamcomp(kf)

           if(kc.le.0.or.kc.gt.mstu(6)) then
              write(6,*)'Invalid code i kf kc',i,kf,kc
              goto 3000
           endif

           id=kchg(kc,5)
           ibar=kchg(kc,6)

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
             if(kf.eq.2212) then
               n_prot=n_prot+1
               call vfill1(15,rap,wei/wy1) 
               call vfill1(11,rap,wei/wy1) 
               call vfill1(21,pt,wei/wp/pt/2) 
               if(i.eq.2) call vfill1(31,rap,wei/wy1) 
               if(i.eq.1) call vfill1(32,rap,wei/wy1) 
             else if(kf.eq.-2212) then
               n_aprot=n_aprot+1
               call vfill1(15,rap,-wei/wy1) 
             endif

c.......Pions.
           else if(kf.eq.-211) then
               call vfill1(12,rap,wei/wy) 
               call vfill1(22,pt,wei/wp/pt/2) 
               n_pim=n_pim+1
           else if(kf.eq.211) then
               call vfill1(13,rap,wei/wy) 
               call vfill1(23,pt,wei/wp/pt/2) 
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
               call jamname(kf,chau)
               write(71,*)k(1,i),k(2,i),p(5,i),' ',chau
               n_mes=n_mes+1
           endif

3000    end do

c....End simulation
      end do

      call jamfin
c     call jamlist(1)


      open(70,file='file70',status='unknown')
      write(70,*)'ylab=',ylab
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
      do i=1,10
       call vscale(10+i,fac)
       call vprint(10+i,mnorm,mform)
      end do
      do i=1,11
       call vscale(30+i,fac)
       call vprint(30+i,mnorm,mform)
      end do
      mform=1
      do i=1,4
       call vscale(20+i,fac)
       call vprint(20+i,mnorm,mform)
      end do

      mform=1
      do i=1,2
       call vscale(50+i,fac)
       call vprint(50+i,mnorm,mform)
      end do

      mform=1
      fac=wei*100.D0/12.D0
      do i=1,1
       call vscale(70+i,fac)
       call vprint(70+i,mnorm,mform)
      end do

c...Multiplicity distributions
c     do i=0,200,2
c      write(53,'(i4,5(e14.6,1x))')i,(mp(j,i)*wei,j=2,5)
c     end do
c     do i=0,200,2
c      write(54,'(i4,2(e14.6,1x))')i,mp(1,i)*wei
c     end do

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
