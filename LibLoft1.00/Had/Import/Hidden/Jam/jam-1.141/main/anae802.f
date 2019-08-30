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

      parc(151)=2.3d0 
      parc(152)=0.3d0 

c     parc(151)=4.0d0 
c     parc(152)=0.30d0 

      ireac=3
      call ana1(ireac)

      do ie=1,1000000

        read(33,end=900,err=999)iev,nv,nbary,nmeson,b
        if(mod(ie,100).eq.0) write(6,*)iev,nv,b,ispec
c       nbary=0
        do i=1,nv
          read(33)(k(j,i),j=1,11),(r(j,i),j=1,5),(p(j,i),j=1,5)
c         if(abs(k(9,i)).eq.3) nbary=nbary+1
c         write(8,800)k(1,i),k(2,i),k(7,i),(r(j,i),j=1,5)
        end do
c       nmeson=nv-nbary
  800   format(i3,1x,i9,1x,i4,3(1x,f8.3))

        call jamdeut
c       write(6,*)'cluster',nbary,mstd(92)
        if(mod(ie,50).eq.0)write(6,*)ie,mstd(92)
c       if(ie.eq.20) stop

        call ana2
c       if(ie.eq.10) goto 900

      end do
900   continue
      close(33)
      call ana3
      stop
999   continue
      write(6,*)'error in reading file'
      end

c***********************************************************************

      subroutine ana1(ireac)

      include 'jam1.inc'
      include 'jam2.inc'

      dimension dy(3,5),y0(3,5),ny(3,5)
      character cy1*3,cy2*3
      save wevt1,wevt2,ylab,wp,wevt
      save ihp,ihpip,ihpin,ihkp,ihkn,ihd
      save dyp,dypip,dypin,dykp,dykn,dyd
      save y0p,y0pip,y0pin,y0kp,y0kn,y0d
      save nyp,nypip,nypin,nykp,nykn,nyd
      save ihlast

c...Rapidity interval for central collisions.
      data (y0(1,i),i=1,5)/0.4d0,0.6d0,0.6d0,1.0d0,0.5d0/
      data (y0(2,i),i=1,5)/0.4d0,0.6d0,0.6d0,1.0d0,0.4d0/
      data (y0(3,i),i=1,5)/0.5d0,0.7d0,0.7d0,0.9d0,0.5d0/
c...Rapidity bin for central collisions.
      data (dy(1,i),i=1,5)/0.2d0,0.2d0,0.4d0,0.4d0,0.2d0/
      data (dy(2,i),i=1,5)/0.2d0,0.2d0,0.2d0,0.4d0,0.2d0/
      data (dy(3,i),i=1,5)/0.2d0,0.2d0,0.2d0,0.4d0,0.2d0/
c...Number of bin for central collisions.
      data (ny(1,i),i=1,5)/9,12,5,3,5/
      data (ny(2,i),i=1,5)/9,12,7,3,6/
      data (ny(3,i),i=1,5)/9,11,8,4,6/


c Si+Al:    peripheral          central
c     p:  y=0.4-2.2 dy=0.2 #=10 y=0.4-2.0  dy=0.2 #= 9
c     pi: y=0.6-2.8 dy=0.2 #=12 y=0.6-2.8  dy=0.2 #=12
c     k+: y=0.8-2.0 dy=0.4 #= 4 y=0.6-2.2  dy=0.4 #= 5
c     k-: y=1.0-1.8 dy=0.4 #= 3 y=1.0-1.8  dy=0.4 #= 3
c     d:  y=0.4-1.0 dy=0.2 #= 4 y=0.5-1.2  dy=0.2 #= 5
c Si+Cu:    peripheral          central
c     p:  y=0.4-2.2 dy=0.2 #=10 y=0.4-2.0  dy=0.2 #= 9
c     pi: y=0.6-2.8 dy=0.2 #=12 y=0.6-2.8  dy=0.2 #=12
c     k+: y=0.8-2.0 dy=0.4 #= 4 y=0.6-1.8  dy=0.2 #= 7
c     k-: y=1.0-1.8 dy=0.4 #= 3 y=1.0-1.8  dy=0.4 #= 3
c     d:  y=0.4-1.0 dy=0.2 #= 4 y=0.4-1.4  dy=0.2 #= 6
c Si+Au:    peripheral          central
c     p:  y=0.5-2.1 dy=0.2 #= 9 y=0.5-2.1  dy=0.2 #= 9
c     pi: y=0.7-2.7 dy=0.2 #=11 y=0.7-2.7  dy=0.2 #=11
c     k+: y=0.9-2.1 dy=0.4 #= 4 y=0.7-2.1  dy=0.2 #= 8
c     k-: y=1.3-1.9 dy=0.6 #= 2 y=0.9-2.1  dy=0.4 #= 4
c     d:  y=0.5-1.1 dy=0.2 #= 4 y=0.5-1.5  dy=0.2 #= 6
     
      ylab=pard(17)
      wevt=1.0d0/dble(mstc(2))

      dyp=dy(ireac,1)
      dypip=dy(ireac,2)
      dypin=dy(ireac,2)
      dykp=dy(ireac,3)
      dykn=dy(ireac,4)
      dyd=dy(ireac,5)

      y0p=y0(ireac,1)
      y0pip=y0(ireac,2)
      y0pin=y0(ireac,2)
      y0kp=y0(ireac,3)
      y0kn=y0(ireac,4)
      y0d=y0(ireac,5)

      nyp=ny(ireac,1)
      nypip=ny(ireac,2)
      nypin=ny(ireac,2)
      nykp=ny(ireac,3)
      nykn=ny(ireac,4)
      nyd=ny(ireac,5)

c...dSig/dEt
      emax=4.0d0
      emin=0.0d0
      we=0.05d0
      nemx=nint((emax-emin)/we)
      call vbook1(1,'dSig/dE all    ',nemx,emin,emax)
      call vbook1(2,'dSig/dE charged',nemx,emin,emax)

c...dEt/d(eta)
      ymax=5.0d0
      ymin=-2.0d0
      wy=0.2d0
      nymx=nint((ymax-ymin)/wy)
      call vbook1(3,'dEt/deta all    ',nymx,ymin,ymax)
      call vbook1(4,'dEt/deta charged',nymx,ymin,ymax)

      call vbook1(11,'dN/dy proton',nymx,ymin,ymax)
      call vbook1(12,'dN/dy pi-',nymx,ymin,ymax)
      call vbook1(13,'dN/dy pi+',nymx,ymin,ymax)
      call vbook1(14,'dN/dy k-',nymx,ymin,ymax)
      call vbook1(15,'dN/dy k+',nymx,ymin,ymax)
      call vbook1(16,'dN/dy deuteon',nymx,ymin,ymax)

      pmin=0.0d0
      pmax=1.5d0
      wp=0.05d0
      npmx=nint((pmax-pmin)/wp)

c...Protons.
      ihist=20
      ihp=ihist
      y=y0p-dyp*3.d0/2.0d0
      do i=1,nyp
        ihist=ihist+1
        y=y+dyp
        write(cy1,'(f3.1)')y
        write(cy2,'(f3.1)')y+dyp
        call vbook1(ihist,'dN/dp proton'//cy1//'-'//cy2,npmx,pmin,pmax)
      end do

c...Positive pions.
      ihpip=ihist
      y=y0pip-dypip*3.d0/2.0d0
      do i=1,nypip
        ihist=ihist+1
        y=y+dypip
        write(cy1,'(f3.1)')y
        write(cy2,'(f3.1)')y+dypin
        call vbook1(ihist,'dN/dp pi+'//cy1//'-'//cy2,npmx,pmin,pmax)
      end do

c...Negative pions.
      ihpin=ihist
      y=y0pin-dypin*3.d0/2.0d0
      do i=1,nypin
        ihist=ihist+1
        y=y+dypin
        write(cy1,'(f3.1)')y
        write(cy2,'(f3.1)')y+dypin
        call vbook1(ihist,'dN/dp pi-'//cy1//'-'//cy2,npmx,pmin,pmax)
      end do

      pmin=0.0d0
      pmax=1.3d0
      npmx=nint((pmax-pmin)/wp)
      ihkp=ihist
c...Positive kaons.
      y=y0kp-dykp*3.d0/2.0d0
      do i=1,nykp
        ihist=ihist+1
        y=y+dykp
        write(cy1,'(f3.1)')y
        write(cy2,'(f3.1)')y+dykp
        call vbook1(ihist,'dN/dp k+'//cy1//'-'//cy2,npmx,pmin,pmax)
      end do

c...Negative kaons.
      ihkn=ihist
      y=y0kn-dykn*3.d0/2.0d0
      do i=1,nykn
        ihist=ihist+1
        y=y+dykn
        write(cy1,'(f3.1)')y
        write(cy2,'(f3.1)')y+dykn
        call vbook1(ihist,'dN/dp k-'//cy1//'-'//cy2,npmx,pmin,pmax)
      end do

c...Deuterons.
      ihd=ihist
      y=y0d-dyd*3.d0/2.0d0
      do i=1,nyd
        ihist=ihist+1
        y=y+dyd
        write(cy1,'(f3.1)')y
        write(cy2,'(f3.1)')y+dyd
        call vbook1(ihist,'dN/dp deut.'//cy1//'-'//cy2,npmx,pmin,pmax)
      end do

      ihlast=ihist

c...Eevent weight
      wevt1=1.0d0/dble(mstc(2))/we
      wevt2=1.0d0/dble(mstc(2))/wy

      return

c*************************************************************************

        entry ana2

c...Loop over all particles
        do i=1,nv
         
         if(k(1,i).ge.11) goto 3000
         kf=k(2,i)
         if(k(1,i).le.4) then
         kc=jamcomp(kf)
         if(kc.le.0.or.kc.gt.mstu(6))then
c          write(6,*)'Invalide code at i, kf, kc : ',i,kf,kc
           go to 3000
         end if
         endif

c...Rapidity cut.
           rap=0.5d0*log( max(p(4,i)+p(3,i),1.d-8)
     $             /max(p(4,i)-p(3,i),1.d-8) )
          if(mstc(4).eq.0) then
          else if(mstc(4).eq.3) then
          else
           rap=rap+ylab
          endif

        ptsq=p(1,i)**2+p(2,i)**2
        pt=sqrt(ptsq)
        pp=sqrt(ptsq+p(3,i)**2)
        pt=max(pt,1.d-8)
        eta=0.5d0*log( max(pp+p(3,i),1.d-8)/max(pp-p(3,i),1.d-8) )
	et=p(4,i)*pt/max(pp,1.d-8)
        emt=sqrt(p(5,i)**2+ptsq)
        emt0=emt-p(5,i)
        ppl=ptsq+p(3,i)**2
        if(mstc(4).eq.2.or.mstc(4).eq.3)then
          bet=pard(5)
          gam=pard(6)
          pl=gam*(p(3,i)+bet*p(4,i))
          el=gam*(p(4,i)+bet*p(3,i))
c         yl=0.5d0*log( max(el+pl,1.d-8)/max(el-pl,1.d-8) )
          ppl=ptsq+pl**2
        endif
        ppl=sqrt(ppl)
c...deuteron
         if(kf.eq.1001001000) goto 2000

        kch=jamchge(k(2,i))

        if(eta.ge.1.25d0.and.eta.le.2.5d0) then
          call vfill1(1,et,wevt1)
          if(kch.ne.0) call vfill1(2,et,wevt1)
        endif

c...Rapidity and transverse momentum distributions.
        if(kf.eq.2212) then
          if(ppl.gt.0.3d0) then
          call vfill1(11,rap,wevt/wy)
          ih=ihp
          ytag=y0p-dyp*3.d0/2.0d0
          do ib=1,nyp
            ih=ih+1
            ytag=ytag+dyp
            if(rap.ge.ytag.and.rap.le.ytag+dyp)
     $              call vfill1(ih,emt0,wevt/(wp*emt*dyp))
          end do
          endif
        else if(kf.eq.-211) then
          if(ppl.gt.0.05d0) then
          call vfill1(12,rap,wevt/wy)
          ih=ihpin
          ytag=y0pin-dypin*3.d0/2.0d0
          do ib=1,nypin
            ih=ih+1
            ytag=ytag+dypin
            if(rap.ge.ytag.and.rap.le.ytag+dypin)
     $              call vfill1(ih,emt0,wevt/(wp*emt*dypin))
          end do
          endif
        else if(kf.eq.211) then
          if(ppl.gt.0.05d0) then
          call vfill1(13,rap,wevt/wy)
          ih=ihpip
          ytag=y0pip-dypip*3.d0/2.0d0
          do ib=1,nypip
            ih=ih+1
            ytag=ytag+dypip
            if(rap.ge.ytag.and.rap.le.ytag+dypip)
     $              call vfill1(ih,emt0,wevt/(wp*emt*dypip))
          end do
          endif
        else if(kf.eq.-321) then
          call vfill1(14,rap,wevt/wy)
          ih=ihkn
          ytag=y0kn-dykn*3.d0/2.0d0
          do ib=1,nykn
            ih=ih+1
            ytag=ytag+dykn
            if(rap.ge.ytag.and.rap.le.ytag+dykn)
     $              call vfill1(ih,emt0,wevt/(wp*emt*dykn))
          end do
        else if(kf.eq.321) then
          call vfill1(15,rap,wevt/wy)
          ih=ihkp
          ytag=y0kp-dykp*3.d0/2.d0
          do ib=1,nykp
            ih=ih+1
            ytag=ytag+dykp
            if(rap.ge.ytag.and.rap.le.ytag+dykp)
     $              call vfill1(ih,emt0,wevt/(wp*emt*dykp))
          end do
        endif

        if(kf.eq.2212.or.kf.eq.2112) goto 3000
        if(eta.ge.1.25d0.and.eta.le.2.5d0) then
          call vfill1(1,et,wevt1)
          if(kch.ne.0) call vfill1(2,et,wevt1)
        endif
        call vfill1(3,eta,et*wevt2)
        if(kch.ne.0) call vfill1(4,eta,wevt2)
        goto 3000

c...deuteron
2000     continue
         call vfill1(16,rap,wevt/wy)
         ih=ihd
          ytag=y0d-dyd*3.d0/2.d0
          do ib=1,nyd
            ih=ih+1
            ytag=ytag+dyd
            if(rap.ge.ytag.and.rap.le.ytag+dykp)
     $        call vfill1(ih,emt0,wevt/(wp*emt*dyd))
          end do


3000     end do

         return

c*************************************************************************

      entry ana3

c...Output hostogram

c...Rapidity
      fac=1.0d0
      do j=1,6
      call vscale(10+j,fac)
      call vprint(10+j,0,1)
      end do

      do j=1,2
      call vscale(j,fac)
      call vprint(j,0,1)
      end do
      do j=3,4
      call vscale(j,1.0d0)
      call vprint(j,0,1)
      end do


c...Mt distributions.
c     fac= paru(1)*10*(parc(4)**2-parc(3)**2)/(2*paru(1)*0.2d0)
      fac=1/paru(2)
      ihist=20
 40   ihist=ihist+1
      if(ihist.gt.ihlast) goto 50
      call vscale(ihist,fac)
      call vprint(ihist,0,1)
      goto 40
 50   continue

      end
