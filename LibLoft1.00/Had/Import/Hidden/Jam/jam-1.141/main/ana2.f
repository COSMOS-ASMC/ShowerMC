c....Analysis program produced from mainFlow.f

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

c          if(ie.eq.1) then
c            write(51,*)(r(j,i),j=1,5) 
c            write(51,*)(v(j,i),j=1,5)
c          endif
c          if(v(4,i).lt.abs(v(3,i))) then
c            write(51,*)i,k(2,i),k(7,i),(v(j,i),j=1,5) 
c          endif

	   if(k(2,i).eq.3122.and.p(5,i).le.1.0) then
	     write(50,*)'lam???',ie,k(1,i),k(2,i),p(5,i)
	   endif
	   if(k(2,i).eq.2212.and.p(5,i).ge.1.0) then
	     write(50,*)'prot???',ie,k(1,i),k(2,i),p(5,i)
	   endif
          if(abs(k(9,i)).eq.3) nbary=nbary+1
        end do

c       call jamdeut
c       call jamlist(1)

        mevt=mevt+1
        call ana2
        if(mod(ie,10).eq.0) write(6,*)iev,nv,b,ispec,mstd(92)

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
      save wy,wp,ylab
      save xspec

      xspec=0

c....Rapidity distribution.
      ylab=pard(17)
      ymin=-7.0d0
      ymax=7.0d0
c     wy=0.25d0
      wy=0.5d0
      nymx=(ymax-ymin)/wy

      yminl=ymin+ylab
      ymaxl=ymax+ylab
      print *,'yminl ymaxl ylab',yminl,ymaxl,ylab

      pmin=0.0d0
      pmax=10.0d0
c     wp=0.1d0
      wp=0.2d0
      npmx=(pmax-pmin)/wp

      call vbook1(1,'dn/dy protons',nymx,ymin,ymax)
      call vbook1(2,'dn/dy antiprotons',nymx,ymin,ymax)
      call vbook1(3,'dn/dy net protons',nymx,ymin,ymax)
      call vbook1(4,'dn/dy pi-',nymx,ymin,ymax)
      call vbook1(5,'dn/dy pi+',nymx,ymin,ymax)
      call vbook1(6,'dn/dy k-',nymx,ymin,ymax)
      call vbook1(7,'dn/dy k+',nymx,ymin,ymax)
      call vbook1(8,'dn/dy lambda',nymx,ymin,ymax)
      call vbook1(9,'dn/dy h-',nymx,ymin,ymax)
      call vbook1(10,'dn/dy antilambda',nymx,ymin,ymax)

      call vbook1(11,'dET/dy protons',nymx,ymin,ymax)
      call vbook1(12,'dET/dy antiprotons',nymx,ymin,ymax)
      call vbook1(13,'dET/dy net protons',nymx,ymin,ymax)
      call vbook1(14,'dET/dy pi-',nymx,ymin,ymax)
      call vbook1(15,'dET/dy pi+',nymx,ymin,ymax)
      call vbook1(16,'dET/dy k-',nymx,ymin,ymax)
      call vbook1(17,'dET/dy k+',nymx,ymin,ymax)
      call vbook1(18,'dET/dy lambda',nymx,ymin,ymax)
      call vbook1(19,'dETdy h-',nymx,ymin,ymax)
      call vbook1(20,'dETdy h-',nymx,ymin,ymax)

      call vbook1(21,'1/mtdn/dmtdy protons',npmx,pmin,pmax)
      call vbook1(22,'1/mtdn/dmtdy antiprotons',npmx,pmin,pmax)
      call vbook1(23,'1/mtdn/dmtdy net potons',npmx,pmin,pmax)
      call vbook1(24,'1/mtdn/dmtdy pion-',npmx,pmin,pmax)
      call vbook1(25,'1/mtdn/dmtdy pion+',npmx,pmin,pmax)
      call vbook1(26,'1/mtdn/dmtdy k-',npmx,pmin,pmax)
      call vbook1(27,'1/mtdn/dmtdy k+',npmx,pmin,pmax)
      call vbook1(28,'1/mtdn/dmtdy lambda',npmx,pmin,pmax)
      call vbook1(29,'1/mtdn/dmtdy h-',npmx,pmin,pmax)

      call vbook1(30,'1/mtdn/dmtdy phi',npmx,pmin,pmax)
      call vbook1(31,'1/mtdn/dmtdy eta',npmx,pmin,pmax)
      call vbook1(32,'1/mtdn/dmtdy xi',npmx,pmin,pmax)
      call vbook1(33,'1/mtdn/dmtdy sigma',npmx,pmin,pmax)
      call vbook1(34,'1/mtdn/dmtdy omega',npmx,pmin,pmax)

      call vbook1(41,'n(pt) net p',nymx,ymin,ymax)
      call vbook1(42,'<pt> net p',nymx,ymin,ymax)
      call vbook1(43,'n(pt) h-',nymx,ymin,ymax)
      call vbook1(44,'<pt> h-',nymx,ymin,ymax)

c....Freaze-out point
      call vbook1(81,'dN/dr proton',40,0.D0,80.D0)
      call vbook1(82,'dN/dr pion  ',40,0.D0,80.D0)
      call vbook1(83,'dN/dr akaon ',40,0.D0,80.D0)
      call vbook1(84,'dN/dr kaon  ',40,0.D0,80.D0)
      call vbook1(85,'dN/dr lambda',40,0.D0,80.D0)
      call vbook1(86,'dN/dr Y     ',40,0.D0,80.D0)

c....Freaze-out time
      call vbook1(91,'dN/dt proton',40,0.D0,80.D0)
      call vbook1(92,'dN/dt pion  ',40,0.D0,80.D0)
      call vbook1(93,'dN/dt akaon ',40,0.D0,80.D0)
      call vbook1(94,'dN/dt kaon  ',40,0.D0,80.D0)
      call vbook1(95,'dN/dt lambda',40,0.D0,80.D0)
      call vbook1(96,'dN/dt Y     ',40,0.D0,80.D0)

      call vbook1(101,'dN/dy - charged',nymx,ymin,ymax)
      call vbook1(102,'dE_t/dy',nymx,ymin,ymax)
      call vbook1(103,'dN/dy - net baryon',nymx,ymin,ymax)
      call vbook1(104,'dN/2piptdy - pion -0.5<y<0.5',npmx,pmin,pmax)

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
        y=0.5d0*log( max(p(4,i)+p(3,i),1.d-8)/max(p(4,i)-p(3,i),1.d-8) )
        ycm=y

        if(mstc(4).eq.0) then
        else if(mstc(4).eq.3) then
        else
         yl=y+ylab
        endif

        ptsq=p(1,i)**2+p(2,i)**2
        pt=sqrt(ptsq)
        pp=sqrt(ptsq+p(3,i)**2)
	ee=sqrt(p(5,i)**2+pp**2)
	if(abs(ee-p(4,i)).gt.1d-7) then
	  print *,'ee p4',ee,p(4,i),kf
	endif
        pt=max(pt,1.d-8)
        eta=0.5d0*log( max(pp+p(3,i),1.d-8)/max(pp-p(3,i),1.d-8) )
c       et=p(4,i)*pt/max(pp,1.d-8)
        et=sqrt(p(5,i)**2+ptsq)
        emt=sqrt(p(5,i)**2+ptsq)
        emt0=emt-p(5,i)

        plab = gamma*( p(3,i) + beta * p(4,i) )
        elab = gamma*( p(4,i) + beta * p(3,i) )
	ppl=sqrt(plab**2+ptsq)
	etlab=elab*pt/max(ppl,1d-8)
	eta_lab=0.5d0*log(max(ppl+plab,1.d-8)/max(ppl-plab,1.d-8))
	etl=eta_lab+ylab
	call vfill1(19,etl,etlab/wy)
	call vfill1(20,y,et/wy)
c       print *,'gam bet',gamma,beta
c       print *,'plab elab',plab,elab,p(3,i),p(4,i)
c       print *,'etlab eta ylab',etlab,eta_lab,ylab
c       print *,'etlab et',etlab,et
c       pause

c...Et.
        call vfill1(102,y,et/wy)

        kch=jamchge(kf)
c...Charged particles.
        if(kch.ne.0) then
          call vfill1(101,y,1d0/wy)
        endif

c...Net baryons.
        kc=jamcomp(kf)
        ibar=kchg(kc,6)
        if(ibar.eq.3) then
          call vfill1(103,y,isign(1,kf)/wy)
        endif

c...Transverse momentum at midrapidiy.
      if(abs(y).le.0.5d0) then
         if(kf.eq.111.or.kfa.eq.211) then
           call vfill1(104,pt,1d0/(paru(2)*pt*wp))
         endif
      endif


        if(kf.eq.-211.or.kf.eq.-321.or.kf.eq.-2212) then
	  call vfill1(9,y,1.0d0/wy)
	  call vfill1(43,y,pt/wy)
	  if(abs(y).le.0.5d0) then
	    call vfill1(29,emt,1.d0/(emt*wp*wy))
	  endif
	endif

c....Protons.
        if(kf.eq.2212) then
	  call vfill1(1,y,1.0d0/wy)
	  call vfill1(3,y,1.0d0/wy)
	  call vfill1(41,y,pt/wy)
	  call vfill1(11,y,et/wy)
	  call vfill1(13,y,et/wy)
	  if(abs(y).le.0.5d0) then
	    call vfill1(21,emt0,1.d0/(emt*wp*wy))
	    call vfill1(23,emt0,1.d0/(emt*wp*wy))
	  endif
c.....anti protons.
        else if(kf.eq.-2212) then
	  call vfill1(2,y,1.0d0/wy)
	  call vfill1(3,y,-1.0d0/wy)
	  call vfill1(41,y,-pt/wy)
	  call vfill1(12,y,et/wy)
	  call vfill1(13,y,-et/wy)
	  if(abs(y).le.0.5d0) then
	    call vfill1(22,emt0,1.d0/(emt*wp*wy))
	    call vfill1(23,emt0,-1.d0/(emt*wp*wy))
	  endif
c.....pion-
        else if(kf.eq.-211) then
          call vfill1(4,y,1.0d0/wy)
          call vfill1(14,y,et/wy)
	  if(abs(y).le.0.5d0) then
	    call vfill1(24,emt0,1.d0/(emt*wp*wy))
	  endif
c.....pion+
        else if(kf.eq.211) then
          call vfill1(5,y,1.0d0/wy)
          call vfill1(15,y,et/wy)
	  if(abs(y).le.0.5d0) then
	    call vfill1(25,emt0,1.d0/(emt*wp*wy))
	  endif
c.....kaon-
        else if(kf.eq.-321) then
          call vfill1(6,y,1.0d0/wy)
          call vfill1(16,y,et/wy)
	  if(abs(y).le.0.5d0) then
	    call vfill1(26,emt0,1.d0/(emt*wp*wy))
	  endif
c.....kaon+
        else if(kf.eq.321) then
          call vfill1(7,y,1.0d0/wy)
          call vfill1(17,y,et/wy)
	  if(abs(y).le.0.5d0) then
	    call vfill1(27,emt0,1.d0/(emt*wp*wy))
	  endif
c.....lambda
        else if(kf.eq.3122) then
          call vfill1(8,y,1.0d0/wy)
          call vfill1(18,y,et/wy)
	  if(abs(y).le.0.5d0) then
	    call vfill1(28,emt0,1.d0/(emt*wp*wy))
	  endif

        else if(kf.eq.-3122) then
          call vfill1(10,y,1.0d0/wy)

c...phi meson
        else if(kf.eq.333) then
	  if(abs(y).le.0.5d0) then
	    call vfill1(30,emt0,1.d0/(emt*wp*wy))
	  endif
c...eta meson
        else if(kf.eq.221) then
	  if(abs(y).le.0.5d0) then
	    call vfill1(31,emt0,1.d0/(emt*wp*wy))
	  endif
c...sigma
        else if(kf.eq.3112.or.kf.eq.3212.or.kf.eq.3222) then
	  if(abs(y).le.0.5d0) then
	    call vfill1(32,emt0,1.d0/(emt*wp*wy))
	  endif
c...xi
        else if(kf.eq.3312.or.kf.eq.3322) then
	  if(abs(y).le.0.5d0) then
	    call vfill1(33,emt0,1.d0/(emt*wp*wy))
	  endif
c...omega-
        else if(kf.eq.3334) then
	  if(abs(y).le.0.5d0) then
	    call vfill1(34,emt0,1.d0/(emt*wp*wy))
	  endif
        endif

c...Hyperons.
      if(kf.eq.3122.or.kf.eq.3212.or.kf.eq.3222.or.kf.eq.3112) then
          call vfill1(18,ycm,1.D0/wy) 
          if(kf.eq.3122.or.kf.eq.3212) call vfill1(17,ycm,1.D0/wy) 
      endif

c------- Freeze out points -----------------------------------
c...Pions
c     if(abs(ycm).le.1.0d0) then
c         vr=sqrt(v(1,i)**2+v(2,i)**2+v(3,i)**2)
c         vt=v(4,i)
c         vr=sqrt(v(1,i)**2+v(2,i)**2)
          vr=sqrt(r(1,i)**2+r(2,i)**2+r(3,i)**2)
          vt=r(4,i) 
          if(vt.lt.abs(r(3,i))) then
            print *,'tau',kf,v(4,i),v(3,i),r(4,i),r(3,i)
          endif
      if(abs(kf).eq.211.or.kf.eq.111) then
          call vfill1(82,vr,1.D0) 
          call vfill1(92,vt,1.D0) 
c...Nucleons
      else if(kf.eq.2112.or.kf.eq.2212) then
         if(abs(k(7,i)).ne.1) then
              call vfill1(81,vr,1.D0) 
              call vfill1(91,vt,1.D0) 
         endif
c...Kaons
      else if(kf.eq.311.or.kf.eq.321) then
          call vfill1(84,vr,1.D0) 
          call vfill1(94,vt,1.D0) 
c...Anti-kaons
      else if(kf.eq.-311.or.kf.eq.-321) then
          call vfill1(83,vr,1.D0) 
          call vfill1(93,vt,1.D0) 
c....Hyperons.
      else if(kf.eq.3122.or.kf.eq.3212.or.kf.eq.3222.or.kf.eq.3112) then
         if(abs(k(7,i)).ne.1) then
              call vfill1(86,vr,1.D0) 
              call vfill1(96,vt,1.D0) 
              if(kf.eq.3122) then
                call vfill1(85,vr,1.D0) 
                call vfill1(95,vt,1.D0) 
              endif
        endif
      endif
c     endif
c------------------------------------------------------------

3000  continue
      xspec=xspec+ispec

      return

c***********************************************************************

      entry ana3(mevt)

c...Output histograms.

      print *,'event',mevt
      wevt=1.d0/dble(mevt)
      fac=wevt
      mnorm=0
      mform=1

c....<pt(y)>
      call vopera(41,'/',3,42,1.0d0,1.0d0)
      call vopera(43,'/',9,44,1.0d0,1.0d0)
      call vprint(42,0,0)
      call vprint(44,0,0)

      call vscale(20,fac)
      call vprint(20,0,0)
c...Rapidty distributuons.
      do i=1,10
      call vscale(i,fac)
      call vprint(i,0,0)
      call vscale(10+i,fac)
      call vprint(10+i,0,0)
      call vscale(20+i,fac)
      call vprint(20+i,mnorm,mform)
      end do

      do i=0,4
      call vscale(30+i,fac)
      call vprint(30+i,mnorm,mform)
      enddo

      do i=1,4
        call vscale(100+i,fac)
        call vprint(100+i,0,1)
      end do

c...Freaze-out points
      do i=1,6
       call vprint(80+i,1,0)
       call vprint(90+i,1,0)
      end do

      end

