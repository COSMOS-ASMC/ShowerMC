c....Analysis program of WA80 Collaboration,  P.L.B361(1995)14.

      include 'jam1.inc'
      include 'jam2.inc'
      common/myana/ispec

      open(33,file='phase.dat',form='unformatted',status='old')

      read(33)nevent,ylab,beta,gamma,icm
      pard(5)=beta
      pard(6)=gamma
      mstc(2)=nevent
      pard(17)=ylab
      mstc(4)=icm

      print *,'nevent ylab beta gamma icm',nevent,ylab,beta,gamma,icm

      call ana1

      mevt=0
      do ie=1,1000000

        read(33,end=900,err=999)iev,nv,nbary,nmeson,b
c       read(33,end=900,err=999)iev,nv,pard(2)
c       b=0d0

	pard(2)=b
        nbary=0
        do i=1,nv
c         read(33)(k(j,i),j=1,11),(r(j,i),j=1,5),(p(j,i),j=1,5)
c         read(33)(k(j,i),j=1,9),(r(j,i),j=1,5),(p(j,i),j=1,5)
          read(33)(k(j,i),j=1,11),(r(j,i),j=1,5),(p(j,i),j=1,5)
     $              ,(v(j,i),j=1,5)


          if(abs(k(9,i)).eq.3) nbary=nbary+1
        end do

c       call jamdeut
c       call jamlist(1)

        mevt=mevt+1
        call ana2
        if(mod(ie,100).eq.0) write(6,*)iev,nv,b,ispec,mstd(92)

      end do
900   continue
      close(33)
      call ana3(mevt)
      stop
999   continue
      write(6,*)'error in reading file'
      end

c***********************************************************************

      subroutine ana1

      include 'jam1.inc'
      include 'jam2.inc'
      common/myana/ispec

      save ylab
      save nyp,nyh,ptp,pth
      save wevt
      save ymin,ymax,wy,nymx
      save pmin,pmax,wp,npmx
      save emin,emax,we,nemx
      dimension nn(-500:500)

C.....Limits and bins for rapidity distributions.
      data ymin,ymax,wy/-3d0,7d0,0.3d0/

C.....Limits and bins for pT-distributions.
c     data pmin,pmax,wp/0d0,5d0,0.1d0/
      data pmin,pmax,wp/0d0,5d0,0.2d0/

C.....Limits and bins for energy-momentum event fluctuations.
      data emin,emax,we/-3d0,7d0,0.15d0/


c....Rapidity distribution.
      ylab=pard(17)

C.....Event weight.
      wevt=1d0
      nymx=(ymax-ymin)/wy
      npmx=(pmax-pmin)/wp
      nemx=(emax-emin)/we

C...Initialize histogram booking.

c...rapidity
      call vbook1(1,'dN/dy negative(pi-k-pbar)',nymx,ymin,ymax)
      call vbook1(2,'(h)dN/dy negative(pi-k-pbar)',nymx,ymin,ymax)
      call vbook1(3,'(s1)dN/dy negative(pi-k-pbar)',nymx,ymin,ymax)
      call vbook1(4,'(s2)dN/dy negative(pi-k-pbar)',nymx,ymin,ymax)

      call vbook1(5,'total dN/dy net p',nymx,ymin,ymax)
      call vbook1(6,'(h)dN/dy net p',nymx,ymin,ymax)
      call vbook1(7,'(s)dN/dy net p',nymx,ymin,ymax)

      call vbook1(8, 'dE/dy     ',nemx,emin,emax)
      call vbook1(9, 'dE/dy hard',nemx,emin,emax)
      call vbook1(10,'dE/dy soft',nemx,emin,emax)


c....Mesons. 
      call vbook1(11,'dN/dy pi-',nymx,ymin,ymax)
      call vbook1(12,'dN/dy pi0',nymx,ymin,ymax)
      call vbook1(13,'dN/dy pi+',nymx,ymin,ymax)
      call vbook1(14,'dN/dy Kaon-',nymx,ymin,ymax)
      call vbook1(15,'dN/dy Kaon0/aKaon0',nymx,ymin,ymax)
      call vbook1(16,'dN/dy Kaon+',nymx,ymin,ymax)
      call vbook1(17,'dN/dy eta',nymx,ymin,ymax)

c....2.1<y<2.9
      call vbook1(21,'dN/dp 2.1<y<2.9 pi-',npmx,pmin,pmax)
      call vbook1(22,'dN/dp 2.1<y<2.9 pi0',npmx,pmin,pmax)
      call vbook1(23,'dN/dp 2.1<y<2.9 pi+',npmx,pmin,pmax)
      call vbook1(24,'dN/dp 2.1<y<2.9 Kaon-',npmx,pmin,pmax)
      call vbook1(25,'dN/dp 2.1<y<2.9 Kaon0/aKaon0',npmx,pmin,pmax)
      call vbook1(26,'dN/dp 2.1<y<2.9 Kaon+',npmx,pmin,pmax)
      call vbook1(27,'dN/dp 2.1<y<2.9 eta  ',npmx,pmin,pmax)

      call vbook1(31,'dN/dm 2.1<y<2.9 pi-',npmx,pmin,pmax)
      call vbook1(32,'dN/dm 2.1<y<2.9 pi0',npmx,pmin,pmax)
      call vbook1(33,'dN/dm 2.1<y<2.9 pi+',npmx,pmin,pmax)
      call vbook1(34,'dN/dm 2.1<y<2.9 Kaon-',npmx,pmin,pmax)
      call vbook1(35,'dN/dm 2.1<y<2.9 Kaon0/aKaon0',npmx,pmin,pmax)
      call vbook1(36,'dN/dm 2.1<y<2.9 Kaon+',npmx,pmin,pmax)
      call vbook1(37,'dN/dm 2.1<y<2.9 eta',npmx,pmin,pmax)



      call vbook1(41,'pt pi0 2.1<y<2.9',npmx,pmin,pmax)
      call vbook1(42,'pt eta 2.1<y<2.9',npmx,pmin,pmax)
      call vbook1(43,'pt eta/pi0 2.1<y<2.9',npmx,pmin,pmax)

      call vbook1(51,'mt pi0 2.1<y<2.9',npmx,pmin,pmax)
      call vbook1(52,'mt eta 2.1<y<2.9',npmx,pmin,pmax)
      call vbook1(53,'mt eta/pi0 2.1<y<2.9',npmx,pmin,pmax)

      call vbook1(61,'pt pi0 ',npmx,pmin,pmax)
      call vbook1(62,'pt eta ',npmx,pmin,pmax)
      call vbook1(63,'pt eta/pi0 ',npmx,pmin,pmax)

      call vbook1(71,'mt pi0 ',npmx,pmin,pmax)
      call vbook1(72,'mt eta ',npmx,pmin,pmax)
      call vbook1(73,'mt eta/pi0 ',npmx,pmin,pmax)




      return

c***********************************************************************

      entry ana2

      beta=pard(5)
      gamma=pard(6)
      ispec=0
c...Loop over all particles.
      do 3000 i=1,nv

       if(k(1,i).ge.10) goto 3000

c...Exclude spectetor.
c       if(abs(k(7,i)).eq.1) then
c         y=0.5*log(max(p(4,i)+p(3,i),1.e-8)/max(p(4,i)-p(3,i),1.e-8))
c         write(3,*)k(1,i),k(2,i),k(7,i),y
c         ispec=ispec+1
c         goto 3000
c       endif

        kf=k(2,i)
        kfa=abs(kf)
        yp=0.5d0*log(max(p(4,i)+p(3,i),1.d-8)/max(p(4,i)-p(3,i),1.d-8))

        ptsq=p(1,i)**2+p(2,i)**2
        pt=sqrt(ptsq)
        emt=sqrt(p(5,i)**2+ptsq)
        pp=sqrt(ptsq+p(3,i)**2)
	ee=sqrt(p(5,i)**2+pp**2)


        pt=max(pt,1.d-8)
        eta=0.5d0*log( max(pp+p(3,i),1.d-8)/max(pp-p(3,i),1.d-8) )
c       et=p(4,i)*pt/max(pp,1.d-8)
        et=sqrt(p(5,i)**2+ptsq)
        emt=sqrt(p(5,i)**2+ptsq)
        emt0=emt-p(5,i)

c       pmra=max(0d0,p(5,i))
c       pt=max(1d-10,sqrt(p(1,i)**2+p(2,i)**2))
c       pt2=p(1,i)**2+p(2,i)**2
c       pmt=max(1d-20,pmra**2+pt2)
c       yp=sign(
c    $    log(min( (sqrt(pmt+p(i,3)**2)+abs(p(i,3)))/sqrt(pmt),
c    &  1d20 )),p(i,3) )


        if(mstc(4).eq.0) then
        else if(mstc(4).eq.3) then
        else
         yp=yp+ylab
        endif

C......VNIBOOK: fill histograms for dET/deta transerve energy.
        pa=sqrt(pt**2+p(3,i)**2)
        et=p(4,i)*pt/max(pa,1.d-8)
        call vfill1(8,yp,et*wevt/we)

        plab = gamma*( p(3,i) + beta * p(4,i) )
        elab = gamma*( p(4,i) + beta * p(3,i) )
	ppl=sqrt(plab**2+ptsq)
	etlab=elab*pt/max(ppl,1d-8)
	eta_lab=0.5d0*log(max(ppl+plab,1.d-8)/max(ppl-plab,1.d-8))
	etl=eta_lab+ylab

c       print *,'gam bet',gamma,beta
c       print *,'plab elab',plab,elab,p(3,i),p(4,i)
c       print *,'etlab eta ylab',etlab,eta_lab,ylab
c       print *,'etlab et',etlab,et
c       pause

c       kch=jamchge(kf)

c....h- distributions.
        if(kf.eq.-211.or.kf.eq.-321.or.kf.eq.-2212) then
          pmra=0.1357d0
          pmt=max(1d-10,pmra**2+p(1,i)**2+p(2,i)**2)
          yp1=sign(log(min( (sqrt(pmt+p(3,i)**2)+abs(p(3,i)))/sqrt(pmt),
     &     1d20 )),p(3,i) )

          if(mstc(4).eq.0) then
          else if(mstc(4).eq.3) then
          else
           yp1=yp1+ylab
          endif
          call vfill1(1,yp1,wevt/wy)
        endif


          ii=0
          if(kf.eq.-211)     ii=11
          if(kf.eq.111)      ii=12
          if(kf.eq.211)      ii=13
          if(kf.eq.-321)     ii=14
          if(abs(kf).eq.311) ii=15
          if(kf.eq.321)      ii=16
          if(kf.eq.221)      ii=17

          if(ii.ne.0) then
	    call vfill1(ii,yp,wevt/wy)
	    if(yp.ge.2.1.and.yp.le.2.9) then
              call vfill1(ii+10,pt,wevt/(pt*wp*0.8))
              call vfill1(ii+20,emt,wevt/(paru(2)*pt*wp*0.8))
c             call vfill1(ii+20,emt,wevt/(paru(2)*emt*wp*0.8))
	      if(ii.eq.12) then
                call vfill1(41,pt,wevt)
                call vfill1(51,emt,wevt)
	      else if(ii.eq.17) then
                call vfill1(42,pt,wevt)
                call vfill1(52,emt,wevt)
	      endif
	    endif
	  endif

	  if(ii.eq.12) then
            call vfill1(61,pt,wevt)
            call vfill1(71,emt,wevt)
	  else if(ii.eq.17) then
            call vfill1(62,pt,wevt)
            call vfill1(72,emt,wevt)
	  endif



c.....protons.
        if(kfa.eq.2212) then
          wei=isign(1,kf)*wevt/wy
          weip=isign(1,kf)*wevt/(wp*pt*0.2d0)
          call vfill1(5,yp,wei)
	endif


3000  continue

      return

c***********************************************************************

      entry ana3(mevt)

c...Output histograms.

      print *,'event',mevt
      wevt=1.d0/dble(mevt)
      fac=wevt
      mnorm=0
      mform=1

c....baryon rapidity distributions.
      do i=1,10
      call vscale(i,fac)
      call vprint(i,mnorm,mform)
      end do


c...mt and pt.
      do i=1,7
        call vscale(10+i,fac)
        call vprint(10+i,mnorm,mform)
        call vscale(20+i,fac)
        call vprint(20+i,mnorm,mform)
        call vscale(30+i,fac)
        call vprint(30+i,mnorm,mform)
      end do

      call vopera(42,'/',41,43,1d0,1d0)
      call vopera(52,'/',51,53,1d0,1d0)
      call vopera(62,'/',61,63,1d0,1d0)
      call vopera(72,'/',71,73,1d0,1d0)
      call vprint(43,0,0)
      call vprint(53,0,0)
      call vprint(63,0,0)
      call vprint(73,0,0)

      end

