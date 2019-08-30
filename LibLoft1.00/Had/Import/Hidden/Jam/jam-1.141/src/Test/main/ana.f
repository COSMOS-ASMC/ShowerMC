c....Analysis program produced from mainFlow.f

      include 'jam1.inc'
      include 'jam2.inc'
      common/myana/ispec
c     real*4 eylab,ebeta,egamma,b,re(5),pe(5)
      dimension ptot(5),ptot1(5)

      open(33,file='phase.dat',form='unformatted',status='old')
      read(33)nevent,eylab,ebeta,egamma,icm
      mstc(2)=nevent
      pard(17)=eylab
      mstc(4)=icm
      print *,'nevent ylab beta gamma icm',nevent,ylab,beta,gamma,icm

      call ana1

      do ie=1,1000000

        read(33,end=900,err=999)iev,nv,b
        nbary=0
        do i=1,nv
          read(33)(k(j,i),j=1,11),(r(j,i),j=1,5),(p(j,i),j=1,5)

c         read(33)(k(j,i),j=1,11),(re(j),j=1,5),(pe(j),j=1,5)
c         do j=1,5
c           r(j,i)=re(j) 
c           p(j,i)=pe(j) 
c         end do

          if(abs(k(9,i)).eq.3) nbary=nbary+1
        end do

        do j=1,5
         ptot(j)=0.0d0
        end do
        do i=1,nv
         if(k(1,i).le.10) then
          do j=1,4
          ptot(j)=ptot(j)+p(j,i)
          end do
         endif
        end do

c       call jamdeut
c       call jamlist(1)
c       print *,'nclust',mstd(92)
c       call jamclust
c       write(6,*)'cluster',nbary,mstd(92)

        do j=1,5
         ptot1(j)=0.0d0
        end do
        do i=1,nv
         if(k(1,i).le.10) then
          do j=1,4
          ptot1(j)=ptot1(j)+p(j,i)
          end do
         endif
        end do

c       srt=sqrt(ptot(4)**2-ptot(1)**2-ptot(2)**2-ptot(3)**2)
c       srt1=sqrt(ptot1(4)**2-ptot1(1)**2-ptot1(2)**2-ptot1(3)**2)
c       write(6,*)ptot(1),ptot1(1),abs(ptot(1)-ptot1(1))
c       write(6,*)ptot(2),ptot1(2),abs(ptot(2)-ptot1(2))
c       write(6,*)ptot(3),ptot1(3),abs(ptot(3)-ptot1(3))
c       write(6,*)srt,srt1,abs(srt-srt1)

        call ana2
c       if(ie.eq.10) goto 900
        if(mod(ie,100).eq.0) write(6,*)iev,nv,b,ispec,mstd(92)

      end do
900   continue
      close(33)
      call ana3
      stop
999   continue
      write(6,*)'error in reading file'
      end

c***********************************************************************

      subroutine ana1

      include 'jam1.inc'
      include 'jam2.inc'
      common/myana/ispec
      save wy,wp,ylab,wevt
      save xspec

      xspec=0
      wevt=1.d0/dble(mstc(2))

c....Rapidity distribution.
      ylab=pard(17)
      ymin=-7.0d0
      ymax=7.0d0
      wy=0.25d0
      nymx=(ymax-ymin)/wy

      yminl=ymin+ylab
      ymaxl=ymax+ylab
      print *,'yminl ymaxl ylab',yminl,ymaxl,ylab

      pmin=0.0d0
      pmax=5.0d0
      wp=0.1d0
      npmx=(pmax-pmin)/wp

      call vbook1(1,'dn/dy nucleons',nymx,ymin,ymax)
      call vbook1(2,'dn/dy pions',nymx,ymin,ymax)
      call vbook1(3,'dn/dy deuteron',nymx,ymin,ymax)

c....Directed and ellipse flow as a function of rapidity.
      call vbook1(11,'v1(y) nucleons',nymx,ymin,ymax)
      call vbook1(12,'v1(y) pions   ',nymx,ymin,ymax)
      call vbook1(13,'v2(y) nucleons',nymx,ymin,ymax)
      call vbook1(14,'v2(y) pions   ',nymx,ymin,ymax)

c....Directed and ellipse flow as a function of pt.
      call vbook1(21,'v2(pt) nucleons',npmx,pmin,pmax)
      call vbook1(22,'v2(pt) pions   ',npmx,pmin,pmax)
      call vbook1(23,'v2(pt) nucleons',npmx,pmin,pmax)
      call vbook1(24,'v2(pt) pions   ',npmx,pmin,pmax)

c....Directed and ellipse flow as a function of rapidity.
      call vbook1(31,'<v1(y)> nucleons',nymx,ymin,ymax)
      call vbook1(32,'<v1(y)> pions   ',nymx,ymin,ymax)
      call vbook1(33,'<v2(y)> nucleons',nymx,ymin,ymax)
      call vbook1(34,'<v2(y)> pions   ',nymx,ymin,ymax)

c....Directed and ellipse flow as a function of pt.
      call vbook1(41,'<v1(pt)> nucleons',npmx,pmin,pmax)
      call vbook1(42,'<v1(pt)> pions   ',npmx,pmin,pmax)
      call vbook1(43,'<v2(pt)> nucleons',npmx,pmin,pmax)
      call vbook1(44,'<v2(pt)> pions   ',npmx,pmin,pmax)

      call vbook1(61,'n(y) nucleons   ',nymx,ymin,ymax)
      call vbook1(62,'n(y) pions      ',nymx,ymin,ymax)
      call vbook1(63,'n(y) protons    ',nymx,yminl,ymaxl)
      call vbook1(64,'n(y) pion+      ',nymx,yminl,ymaxl)

      call vbook1(65,'n(pt) nucleons  ',npmx,pmin,pmax)
      call vbook1(66,'n(pt) pions     ',npmx,pmin,pmax)

c....Transverse flow.
      call vbook1(51,'px(y) nucleons',nymx,ymin,ymax)
      call vbook1(52,'pt(y) nucleons',nymx,ymin,ymax)

      call vbook1(53,'px(y) pions   ',nymx,ymin,ymax)
      call vbook1(54,'pt(y) pions   ',nymx,ymin,ymax)

      call vbook1(55,'px(y) protons ',nymx,yminl,ymaxl)
      call vbook1(56,'pt(y) protons ',nymx,yminl,ymaxl)

      call vbook1(57,'px(y) pion+   ',nymx,yminl,ymaxl)
      call vbook1(58,'pt(y) pion+   ',nymx,yminl,ymaxl)


      call vbook1(71,'<px(y)> nucleons',nymx,ymin,ymax)
      call vbook1(72,'<pt(y)> nucleons',nymx,ymin,ymax)

      call vbook1(73,'<px(y)> pions   ',nymx,ymin,ymax)
      call vbook1(74,'<pt(y)> pions   ',nymx,ymin,ymax)

      call vbook1(75,'<px(y)> protons ',nymx,yminl,ymaxl)
      call vbook1(76,'<pt(y)> protons ',nymx,yminl,ymaxl)

      call vbook1(77,'<px(y)> pion+   ',nymx,yminl,ymaxl)
      call vbook1(78,'<pt(y)> pion+   ',nymx,yminl,ymaxl)



      call vbook1(81,'n(y) deuterons  ',nymx,ymin,ymax)
      call vbook1(82,'px(y) deuterons ',nymx,ymin,ymax)
      call vbook1(83,'<px(y)> deuterons ',nymx,ymin,ymax)

      call vbook1(84,'n(y) kaons      ',nymx,ymin,ymax)
      call vbook1(85,'px(y) kaons     ',nymx,ymin,ymax)
      call vbook1(86,'<px(y)> kaons     ',nymx,ymin,ymax)

      call vbook1(87,'n(y) antikaons  ',nymx,ymin,ymax)
      call vbook1(88,'px(y) antikaons ',nymx,ymin,ymax)
      call vbook1(89,'<px(y)> antikaons ',nymx,ymin,ymax)

      return

c***********************************************************************

      entry ana2

      ispec=0
c...Loop over all particles.
      do 3000 i=1,nv

       if(k(1,i).ge.10) goto 3000

c...Exclude spectetor.
        if(abs(k(7,i)).eq.1) then
c         y=0.5*log(max(p(4,i)+p(3,i),1.e-8)/max(p(4,i)-p(3,i),1.e-8))
c         write(3,*)k(1,i),k(2,i),k(7,i),y
          ispec=ispec+1
          goto 3000
        endif

c       if(k(1,i).eq.5) then
c         write(3,*)k(1,i),k(2,i)
c       endif

        kf=k(2,i)
        y=0.5d0*log( max(p(4,i)+p(3,i),1.d-8)/max(p(4,i)-p(3,i),1.d-8) )

        if(mstc(4).eq.0) then
        else if(mstc(4).eq.3) then
        else
         yl=y+ylab
        endif

        ptsq=p(1,i)**2+p(2,i)**2
        pt=sqrt(ptsq)
        pp=sqrt(ptsq+p(3,i)**2)
        pt=max(pt,1.d-8)
        eta=0.5d0*log( max(pp+p(3,i),1.d-8)/max(pp-p(3,i),1.d-8) )
	et=p(4,i)*pt/max(pp,1.d-8)
	px=p(1,i)
	if(pt.gt.1d-8) then
	  cos1=px/pt
	else
	  cos1=0.0d0
	endif
	phi=acos(cos1)
        cos2=cos(2*phi)

c...deuterons.
       if(kf.eq.1001001000) then
	  call vfill1(81,y,1.0d0)
	  call vfill1(82,y,px)
       else if(kf.eq.321.or.kf.eq.311) then
	  call vfill1(84,y,1.0d0)
	  call vfill1(85,y,px)
       else if(kf.eq.-321.or.kf.eq.-311) then
	  call vfill1(87,y,1.0d0)
	  call vfill1(88,y,px)

c......Transverse flow of nucleons and pions.
       else if(kf.eq.2212.or.kf.eq.2112) then
	  call vfill1(61,y,1.0d0)
	  call vfill1(65,pt,1.0d0)
	  call vfill1(51,y,px)
	  call vfill1(52,y,pt)
          call vfill1(11,y,cos1)
          call vfill1(13,y,cos2)
          call vfill1(21,pt,cos1)
          call vfill1(23,pt,cos2)

          call vfill1(1,y,wevt/wy)
c....Protons.
          if(kf.eq.2212) then
            yll=y+ylab
	    call vfill1(63,yll,1.0d0)
	    call vfill1(55,yll,px)
	    call vfill1(56,yll,pt)
          endif
        else if(kf.eq.111.or.abs(kf).eq.211) then
          call vfill1(2,y,wevt/wy)

	  call vfill1(62,y,1.0d0)
	  call vfill1(66,pt,1.0d0)
	  call vfill1(53,y,px)
	  call vfill1(54,y,pt)
          call vfill1(12,y,cos1)
          call vfill1(14,y,cos2)
          call vfill1(22,pt,cos1)
          call vfill1(24,pt,cos2)
c.....Positive pions.
          if(kf.eq.211) then
            yll=y+ylab
	    call vfill1(64,yll,1.0d0)
	    call vfill1(57,yll,px)
	    call vfill1(58,yll,pt)
          endif
	endif

3000  continue
      xspec=xspec+ispec

      return

c***********************************************************************

      entry ana3

      print *,'# of spectator',xspec*wevt
c...Output histograms.

      fac=1.0d0
      mnorm=0
      mform=1

c...Rapidty distributuons.
      call vprint(1,mnorm,mform)
      call vprint(2,mnorm,mform)

c....Directed and ellipse flow.
      call vopera(11,'/',61,31,1.0d0,1.0d0)
      call vopera(12,'/',62,32,1.0d0,1.0d0)
      call vopera(13,'/',61,33,1.0d0,1.0d0)
      call vopera(14,'/',62,34,1.0d0,1.0d0)

      call vopera(21,'/',65,41,1.0d0,1.0d0)
      call vopera(22,'/',66,42,1.0d0,1.0d0)
      call vopera(23,'/',65,43,1.0d0,1.0d0)
      call vopera(24,'/',66,44,1.0d0,1.0d0)

c....Transverse flow.
      call vopera(51,'/',61,71,1.0d0,1.0d0)
      call vopera(52,'/',61,72,1.0d0,1.0d0)

      call vopera(53,'/',62,73,1.0d0,1.0d0)
      call vopera(54,'/',62,74,1.0d0,1.0d0)

      call vopera(55,'/',63,75,1.0d0,1.0d0)
      call vopera(56,'/',63,76,1.0d0,1.0d0)

      call vopera(57,'/',64,77,1.0d0,1.0d0)
      call vopera(58,'/',64,78,1.0d0,1.0d0)


      call vprint(61,mnorm,mform)
      call vprint(62,mnorm,mform)
      do i=1,4
	call vprint(30+i,mnorm,mform)
	call vprint(40+i,mnorm,mform)
      end do
      do i=1,8
	call vprint(70+i,mnorm,mform)
      end do

      call vopera(82,'/',81,83,1.0d0,1.0d0)
      call vopera(85,'/',84,86,1.0d0,1.0d0)
      call vopera(88,'/',87,89,1.0d0,1.0d0)
      call vprint(83,mnorm,mform)
      call vprint(86,mnorm,mform)
      call vprint(89,mnorm,mform)

      end

