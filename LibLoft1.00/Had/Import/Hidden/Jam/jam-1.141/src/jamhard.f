c***********************************************************************
c  List of subprograms in order of appearance, with main purpose       *
c  (S = subroutine, F = function, B = block data)                      *
c                                                                      *
c s jamhard  to perform jets production                                *
c s jamhrdv  to extract space-time information                         *
c f jamprtim to calculate reaction time for partonic time evolutions   *
c                                                                      *
c***********************************************************************
c***********************************************************************

      subroutine jamhard(n_jet,jp,jt,jflg)

c...Purpose: to perform jets production.
c...JFLG=1 jets valence quarks,
c...JFLG=2 with gluon jet,
c...JFLG=3 with q-qbar prod for (JP,JT).
c...JFLG=0 jets can not be produced this time.
c...JFLG=-1, error occured abandon this event.

c...JAM common block.
      include 'jam1.inc'
      include 'jam2.inc'
      real*8 jamemjet
c...HIJING common block.
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      common/hipyint/mint4,mint5,atco(200,20),atxs(0:200)
      save  /hiparnt/,/hipyint/
c...PYHIA common block.
      common/jyjets/njet,npad,kjet(1000,5),pjet(1000,5),vjet(1000,5)
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      common/pjint5/ngenpd,ngen(0:500,3),xsec(0:500,3)
      save  /jyjets/,/pjsubs/,/pjpars/,/pjint1/,/pjint2/,/pjint5/

c...Commonblocks to use multiple scattering.
      common/jampyda1/qmult(10),imult(2,10)
      common/jamhrdev1/lead(2)
      save /jampyda1/,/jamhrdev1/

c...Used for space-time information.
      common/jamxhit/xo(2,5),srtc,ecmc

      save mstp122
      save kfcol1,kfcol2

c...Local values.
      real*8 jamdtim
      dimension pcm(5),psum(5)
      dimension ind(100),idel(100),kfv(2),k9v(2),kfq(2,2),jpq(2,2)
     $,ipa(2,2)
      dimension pjet0(2,5),pc(5),rc(4)
      dimension indd(100)

      data kfcol1, kfcol2/2212,2212/
      data mstp122/1/

      ih=mstc(38)

c...Itial settings.
      jflg=0
      isw1=1 ! =1: all produced particles are string
      isw2=1 ! =1: all produced particles are string
      nlq=0
      lead(1)=-1
      lead(2)=-1
      ksimul=k(8,jp)

c...Set incident momenta, energies and KF codes for proj. and targ.
      call jamcmom(jp,jt,pjet0,kfv,k9v,kfq,srt)

      kf01=kfv(1)
      kf02=kfv(2)
      do jm=1,2
        jpq(jm,1)=0
        jpq(jm,2)=0
        if(jm.eq.1) ip=jp
        if(jm.eq.2) ip=jt
        if(mod(abs(k(1,ip)),10).eq.4) then
          jpq(jm,1)=k(10,ip)
          jpq(jm,2)=k(11,ip)
        endif
      end do


      ihis0=min(1000,abs(k(6,jp))+abs(k(6,jt))+1)
      do j=1,5
        xo(1,j)=r(j,jp)
        xo(2,j)=r(j,jt)
      end do
   
c...c.m. energy and momentum.
      do j=1,4
        pcm(j)=pjet0(1,j)+pjet0(2,j)
      end do
      s=pcm(4)**2-(pcm(1)**2+pcm(2)**2+pcm(3)**2)
      if(s.le.0d0) then
         call jamerrm(1,0,'(jamhard:) s<0')
         return
      endif
      srtc=sqrt(s)
      ecmc=pcm(4)
      if(srtc.le.pjet0(1,5)+pjet0(2,5)) then
         call jamerrm(1,0,'(jamhard:) s<em1+em2')
         return
      endif
      pr=sqrt((s-(pjet0(1,5)+pjet0(2,5))**2)*
     $ (s-(pjet0(1,5)-pjet0(2,5))**2))/(2*srtc)

c...Must have enough energy to produce jets.
      d1=sqrt(pjet0(1,5)**2+pr**2)
      d2=sqrt(pjet0(2,5)**2+pr**2)
      ecut1=hipr1(1)+hipr1(8)+pjmass(kfq(1,1))+pjmass(kfq(1,2))
      ecut2=hipr1(1)+hipr1(8)+pjmass(kfq(2,1))+pjmass(kfq(2,2))
      if(d1.le.ecut1.or.d2.le.ecut2) then
        jflg=-1
        return
       endif

c...Add additional gluon to the total initial momentum.
      do jm=1,2
        jpq1=jpq(jm,1)
        jpq2=jpq(jm,2)
        if((jpq1.ne.0.and.jpq2.ne.0).and.(jpq2.gt.jpq1)) then
          do i=jpq1+1,jpq2-1
            do j=1,4
              pcm(j)=pcm(j)+p(j,i)
            end do
          end do
        endif
      end do
      pcm(5)=sqrt(max(0d0,pcm(4)**2-(pcm(1)**2+pcm(2)**2+pcm(3)**2)))

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      if(mstc(8).ge.3) then
        write(ih,*)'(jamhard:) nv nmeson',nv,nmeson
        write(ih,*)'srt pcm5',srt,pcm(5),kfv(1),kfv(2)
        write(ih,*)'k1 kf',jp,k(1,jp),k(2,jp),p(5,jp),kq(1,jp),kq(2,jp)
        write(ih,*)'k1 kf',jt,k(1,jt),k(2,jt),p(5,jt),kq(1,jt),kq(2,jt)
        write(ih,*)'jpq1 jpq2 a',jpq(1,1),jpq(1,2)
        write(ih,*)'jpq1 jpq2 b',jpq(2,1),jpq(2,2)
      endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      kst01=k(1,jp)
      kst02=k(1,jt)
      ihnt2(11)=jp
      ihnt2(12)=jt
      icltag=mstd(29)+1
      nmeson0=nmeson
      nv0=nv
      ncup=0
      ind(1)=0
      kfcol3=kf01
      kfcol4=kf02

c...Initialize Pythia.
      if((kfcol1.ne.kfcol3).or.(kfcol2.ne.kfcol4)) then
        do j=1,5
        pjet(1,j)=pjet0(1,j)
        pjet(2,j)=pjet0(2,j)
        end do

        mstp122=mstp122+1
        if(mstp122.ge.5) mstp(122)=0
        call pjinit('five',kfcol3,kfcol4,0.0d0,icon)
        if(icon.ne.0) then
           write(check(1),
     $  '(''icon srtc kfcol3,kfcol4'',i4,1x,g12.3,1x,i8,1x,i8)')
     $ icon,srtc,kcol3,kcol4
           call jamerrm(1,0,'iflg after pjinit')
          jflg=-2
          return
        endif
        kfcol1=kfcol3
        kfcol2=kfcol4
      endif

c----------------------
        miss=0
        miss2=0
        miss3=0
        misp=0
        mist=0
155   continue
      miss=miss+1
      if(miss.gt.200) then
        write(check(1),'(''n_jet='',i4)')n_jet
        call jamerrm(1,1,'(jamhrd) no maching jet number miss>200')
        n_jet=n_jet-1
c       miss=0
        miss2=miss2+1
        if(n_jet.eq.0.or.miss2.ge.50) then
         write(check(1),'(''n_jet='',i4)')n_jet
         call jamerrm(30,1,'(jamhrd:) no maching jet miss2')
        endif
      endif

152   continue

      do i=1,10
        imult(1,i)=0
        imult(2,i)=0
        qmult(i)=0.0d0
      end do

      do j=1,5
        pjet(1,j)=pjet0(1,j)
        pjet(2,j)=pjet0(2,j)
      end do

      do i=1,ncup
        k(1,ind(i))=30
      end do
      ncup=0
      nmeson=nmeson0
      nv=nv0

c... Perform hard scattering between jp and jt hadrons.
      call pjevnt

      if(msti(61).eq.1) then
        miss3=miss3+1
        if(miss3.lt.1000) goto 152
        write(check(1),'(''srt='',g13.3)')pcm(5)
        call jamerrm(1,1,'(jamhrd:) miss3>1000')
        jflg=-3
        return
      endif


c...Number of hard collision determined from HIJING in this time
      j_jet=mint(31)
      if(mstc(83).eq.1.and.j_jet.ne.n_jet) go to 155

      if(isw2.eq.1) then
        do i=9,njet
          kfa=abs(kjet(i,2))
          jc2=mod(kfa/10,10)
          if(kfa.eq.21) then
          else if(kfa.lt.10) then
          else if(kfa.gt.1000.and.jc2.eq.0) then
          else
           goto 155
          endif
        end do
      endif


c...q-qbar jets must has minimum mass HIPR1(1)
c       if(kjet(7,2).eq.-kjet(8,2)) then
c          qmass2=(pjet(7,4)+pjet(8,4))**2-(pjet(7,1)+pjet(8,1))**2
c    &           -(pjet(7,2)+pjet(8,2))**2-(pjet(7,3)+pjet(8,3))**2
c          qm=ulmass(kjet(7,2))
c          if(qmass2.lt.(2.0*qm+hipr1(1))**2) go to 155
c       endif

      if(isw1.eq.1) then
        pxp=p(1,jp)-pjet(3,1)
        pyp=p(2,jp)-pjet(3,2)
        pzp=p(3,jp)-pjet(3,3)
        pep=p(4,jp)-pjet(3,4)

        pxt=p(1,jt)-pjet(4,1)
        pyt=p(2,jt)-pjet(4,2)
        pzt=p(3,jt)-pjet(4,3)
        pet=p(4,jt)-pjet(4,4)

c...If the remain energy<ECUT,
c...the proj or targ can not produce jet anymore
        if(pep.le.ecut1) then
          misp=misp+1
          if(misp.lt.50) go to 155
          jflg=-4
          return
        endif
        if(pet.le.ecut2) then
          mist=mist+1
          if(mist.lt.50) go to 155
          jflg=-5
          return
        endif

c...The total W+, W- must be positive
        wp=pep+pzp+pet+pzt
        wm=pep-pzp+pet-pzt
        if(wp.lt.0.0d0 .or. wm.lt.0.0d0) then
          miss=miss+1
          if(miss.lt.50) go to 155
          jflg=-6
          return
        endif

c...The proj and targ remnants must have at least
c...a CM energy that can produce two strings
c...with minimum mass HIPR1(1)(see HIJSFT HIJFRG)
        sw=wp*wm
        ampx=sqrt((ecut1-hipr1(8))**2+pxp**2+pyp**2+0.01d0)
        amtx=sqrt((ecut2-hipr1(8))**2+pxt**2+pyt**2+0.01d0)
        sxx=(ampx+amtx)**2
        if(sw.lt.sxx.or.vint(43).lt.hipr1(1)) then
          miss=miss+1
          if(miss.lt.50) go to 155
          jflg=-7
          return
        endif  

      endif

c...1 is proj, 2 is targ.
c...3, 4 : before initial radiation
c...5, 6 : after initial radiation
c...7, 8 : jet
c
c...      x quark
c...   xx0x di-quark
c...     1x leptons
c...     xx gauge and Higgs bosons
c...    xxx meson
c...  x0xxx meson
c...   xxxx baryon
c...  xxxxx baryon

c...Hard scattering is successful.
      if(mstc(8).ge.3) then
        call pjlist(1)
        do i=9,njet
        write(mstc(38),'(i3,5(1x,g10.3))')i,(vjet(i,l),l=1,5)
        end do
      endif

c...Flag for participant.
      ipart1=0
      ipart2=0

c...Loop over produced partons.
      i=9
      ijet=0
3000  continue
      ijet=ijet+1
      if(ijet.ge.1000) then
        ih=mstc(38)
        write(ih,*)'Are there many particle? after hard scatt.'
        write(ih,*)'njet ijet i',njet,ijet,i
        write(ih,*)'jp jt',jp,jt
        write(ih,*)'k1 k2',k(2,jp),p(5,jp),k(2,jt),p(5,jt),pcm(5)
        call pjlist(1)

        write(ih,*)'Are there many particle? after hard scatt.'
        write(ih,*)'njet ijet i',njet,ijet,i
        write(ih,*)'jp jt',jp,jt
        write(ih,*)'k1 k2',k(2,jp),p(5,jp),k(2,jt),p(5,jt),pcm(5)
       call jamerrm(30,0,'(jamhrd:) many particle? after hard catt.')
      endif


c...End loop over partons.
      if(i.gt.njet) goto 2000

c...This is a hadron
c===========================================================
      if(kjet(i,1).eq.1) then
c===========================================================

        kf=kjet(i,2)
        kfa=abs(kf)
        jc2=mod(kfa/10,10)
        if((kfa.gt.10.and.kfa.le.100).and.kfa.ne.21) then
        else if(jc2.ne.0.and.kfa.gt.100) then
        else
          write(check(1),8000)i,(kjet(i,l),l=1,5)
 8000 format('i=',i9,'kjet',i4,1x,i9,1x,i4,1x,i4,1x,i4)
          call jamerrm(30,1,'(jamhrd:)something is wrong id')
        endif
            
c...Hadron remaines/leptons

c...Find index of JAM array.
      if(kjet(i,3).eq.1) then
        jtp=jp
        ks=kst01
        k6=k(6,jt)
        iorg=1
      else if(kjet(i,3).eq.2) then
        jtp=jt
        ks=kst02
        k6=k(6,jp)
        iorg=2
      else
       write(mstc(38),*)'(jamhrd:hadron)???',i,(kjet(i,l),l=1,5)
       call pjlist(1)
       goto 155
      endif

      nmeson=nmeson+1
      nv=nv+1
      if(nv.gt.mxv)
     $   call jamerrm(30,0,'particle too large (jamhrd)mxv->change')
      indx=nv
      ncup=ncup+1
      ind(ncup)=indx

c.....Copy pythia results to JAM array.
      kf=kjet(i,2)
      kc=jamcomp(kf)
      call jamkupda(5,indx,kf,kc,ks,k6,icltag,5)
      do j=1,4
        p(j,indx)=pjet(i,j)
        v(j,indx)=xo(iorg,j)
        r(j,indx)=xo(iorg,j)
      end do
      p(5,indx)=pjet(i,5)
      r(5,indx)=xo(iorg,4)
      v(5,indx)=r(5,indx)
     $             +jamdtim(1,kf,kc,k(1,indx),p(5,indx),p(4,indx))

      kq(1,indx)=0
      kq(2,indx)=0
      do kk=1,10
        vq(kk,indx)=0.0d0
      end do

      i=i+1
      goto 3000

c-- store parton system -----------------------------
c===========================================================
      else if(kjet(i,1).eq.2) then
c===========================================================

         lc=0
         idq=0
         jtp=0
         pc1=0.0d0
         pc2=0.0d0
         pc3=0.0d0
         pc4=0.0d0

c.....Get total momemtum of the color singlet system.
         do j=i,njet

           lc=lc+1
           if(abs(kjet(j,2)).ge.1000) then
              idq=idq+1
           endif
           pc1=pc1+pjet(j,1)
           pc2=pc2+pjet(j,2)
           pc3=pc3+pjet(j,3)
           pc4=pc4+pjet(j,4)
           if(kjet(j,1).eq.1) goto 22

         end do
   22    continue

         pc5sq=pc4**2-pc1**2-pc2**2-pc3**2
         if(pc5sq.le.0d0) then
           write(mstc(38),*)'(jamhrd:) pc5sq<0',pc5sq,i,i+lc-1
           write(mstc(38),*)'jp',jp,k(1,jp),k(2,jp),(p(j,jp),j=1,5)
           write(mstc(38),*)'jt',jp,k(1,jt),k(2,jt),(p(j,jt),j=1,5)
           call pjlist(1)
           goto 155
         endif
         pc5=sqrt(pc4**2-pc1**2-pc2**2-pc3**2)
         gg=pc4/pc5
         htime=0.0d0
         if(mstc(76).ge.1)
     $   htime=-parc(55)*log(max(1.d0-rn(0),1.d-35))*gg

c------------------------------------------------------------
c...Jet system
c------------------------------------------------------------

 2500   continue
        imax=i+lc-1
        k10=nv+1
        k11=nv+lc   ! max. of color flow

c.....Check this particle is participant or not.
c....XXX BB and BB~ collision only

      ix1=0
      ix2=0
      do jm=1,2
        ipa(jm,1)=0
        ipa(jm,2)=0
        if(jm.eq.1) jj=i
        if(jm.eq.2) jj=imax
        if(abs(kjet(jj,2)).ge.1000) then
          if(ipart1.eq.0.and.kjet(jj,3).eq.1) then
            if(jpq(1,1).gt.0.and.jpq(1,2).gt.0) then
              ipa(1,1)=jpq(1,1)+1
              ipa(1,2)=jpq(1,2)-1
              k11=k11+ipa(1,2)-ipa(1,1)+1
              ipart1=k11
              ix1=1
            else
              ipart1=k11
              ix1=2
            endif
            if(jm.eq.1) then
              ipart1=k10
              ix1=3
            endif
          else if(ipart2.eq.0.and.kjet(jj,3).eq.2)then
            if(jpq(2,1).gt.0.and.jpq(2,2).gt.0) then
              ipa(2,1)=jpq(2,1)+1
              ipa(2,2)=jpq(2,2)-1
              k11=k11+ipa(2,2)-ipa(2,1)+1
              ipart2=k11
              ix2=1
            else
              ipart2=k11
              ix2=2
            endif
            if(jm.eq.1) then
              ipart2=k10
              ix2=ix2+3
            endif
          endif
        endif
      end do

      if((ix1.eq.1.or.ix1.eq.4).and.(ix2.eq.1.or.ix2.eq.4)) then
        k11=nv+lc+jpq(1,1)-jpq(1,2)-1+jpq(2,1)-jpq(2,2)-1
        if(ix1.eq.1)ipart1=k11
        if(ix2.eq.1)ipart2=k11
      endif

      if(k11.gt.mxv) 
     $  call jamerrm(30,0,'particle too large (jamhrd)mxv')


c......Loop over color connected partons.
        do j=i,imax

          if(j.eq.imax) then
            do jm=1,2
              if(ipa(jm,1).ne.0) then
                do ip=ipa(jm,1),ipa(jm,2)
                 nv=nv+1
                 nmeson=nmeson+1
                 call jamexch(nv,ip)
                 k(10,nv)=k10
                 k(11,nv)=k11
                 ncup=ncup+1
                 ind(ncup)=nv
              end do
              endif
            end do
          endif

          nmeson=nmeson+1
          nv=nv+1
          indx=nv

          if(j.eq.1) then
            nlq=nlq+1
          endif

          kf=kjet(j,2)
          kc=jamcomp(kf)
          k(1,indx)=4
          k(2,indx)=kf
          k(3,indx)=1000*mste(1)+mste(2)
          k(4,indx)=1000*abs(kcp(2,1))+abs(kcp(2,2))
          k(5,indx)=mstd(29)+1
          k(6,indx)=ihis0
          k(7,indx)=2
          k(8,indx)=ksimul
          kfa=abs(kf)
          kflc=mod(kfa/10,10) 
          if(kfa.lt.10) then
            k(9,indx)=1*isign(1,kf)
          else if(kfa.gt.1000.and.kflc.eq.0) then
            k(9,indx)=2*isign(1,kf)
          else if(kfa.eq.21) then
            k(9,indx)=0
          endif

          k(10,indx)=k10                ! min. of color flow
          k(11,indx)=k11                ! max. of color flow

c...Set production vertex
          do im=1,5
            p(im,indx)=pjet(j,im)
            r(im,indx)=vjet(j,im)
            v(im,indx)=r(im,indx)
          end do
          r(5,indx)=r(4,indx)

c...Hadronization time
          v(5,indx)=htime+r(5,indx)

          kq(1,indx)=0
          kq(2,indx)=0
          do kk=1,10
            vq(kk,indx)=0.0d0
          end do

c...Save line number for updating collision array.
          ncup=ncup+1
          ind(ncup)=indx

        end do

        i=i+lc
        goto 3000

c===========================================================
      else
c===========================================================

        write(check(1),'(''kjet(i,1)'',i4)')kjet(i,1)
        call jamlist(1)
        call jamerrm(30,1,'(jamhrd)invalid jet system')

c===========================================================
      endif
c===========================================================

2000  continue



c======= reset the energy-momentum of incoming particles ========
900   continue

c....Check minimum string mass.
      ip=ind(1)-1
      imax=ind(ncup)
3400  continue
      ip=ip+1
      if(ip.ge.imax) goto 3440
      if(k(1,ip).ne.4) goto 3400

      k10=k(10,ip)
      k11=k(11,ip)
      kf1=k(2,k10)
      kf2=k(2,k11)
      i=ip
      ip=k11
      if(kf1.eq.21.or.kf2.eq.21) goto 3400

c.....Get total momemtum of the color singlet system.
      lc=0
      r4max=-100d0
      do j=1,4
        pc(j)=0.0d0
        rc(j)=0.0d0
      end do
      do ii=k10,k11
        lc=lc+1
        do j=1,4
        pc(j)=pc(j)+p(j,ii)
        rc(j)=rc(j)+r(j,ii)
        r4max=max(r4max,r(4,ii))
        end do
      end do
      pc(5)=sqrt(max(0d0,pc(4)**2-pc(1)**2-pc(2)**2-pc(3)**2))
      if(pc(5).le.0d0) goto 155

      if(mstc(8).ge.3) then
        write(ih,*)'check mass of parton system',k10,k11,pc(5)
     $  ,jamemjet(kf1,kf2),kf1,kf2
      endif

c-------------------------------------------------------------------
c...This parton system have enough mass to form string?

c...Check minimum string mass.
      if(pc(5).gt.jamemjet(kf1,kf2)) goto 3400

      call kfcnst(kf1,kf2,kf0,pc(5))
      kf=kf0
      call jamidres(kf,pc(5),icon)
      if(icon.ne.0) goto 155

c...Check minimum mass.
      call jamdmass(kf,kfm,kfd,emin,emdn)
      if(kf.eq.kfm.and.abs(pc(5)-emin).lt.1d-5) then
      else if(kf.eq.kfm.or.pc(5).lt.emdn) then
        write(check(1),8200)kf1,kf2,kf,pc(5),kf,kfm,emin,emdm
 8200   format('kf1 kf2 kfj',3i10,' pc5',g12.3,' kf kfm emin emdm',
     $  i7,1x,i7,1x,g10.3,1x,g10.3)
        call jamerrm(1,1,'(jamhard1:chh) mass small ')
        goto 155
      endif


      kc=jamcomp(kf)
      if(kc.eq.0) then
        write(check(1),8100)kf1,kf2,kf
 8100   format('kf1 kf2 kfj',3i10)
        call jamerrm(1,1,'(jamhrd:chh) no particle code')
        goto 155
      endif
 
cxxxxxxxxxxxxxxxxxxxxxxxxxxx
      if(mstc(8).ge.3) then
        write(ih,*)'convert parton systgem to hadron',kf1,kf2,kf,pc(5)
      endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxx

      nmeson=nmeson+1
      nv=nv+1
      if(nv.gt.mxv) call jamerrm(30,0,'(jamhrd:chh)particle too large')
      indx=nv
      do j=1,5
      p(j,indx)=pc(j)
      end do
      call jamkupda(6,indx,kf,kc,1,ihis0,icltag,5)
      kq(1,indx)=0
      kq(2,indx)=0
      do kk=1,10
        vq(kk,indx)=0.0d0
      end do
      do j=1,3
        r(j,indx)=rc(j)/lc
        v(j,indx)=r(j,indx)
      end do
      r(4,indx)=r4max
      r(5,indx)=r(4,indx)
      v(5,indx)=r(5,indx)
     $            +jamdtim(1,kf,kc,k(1,indx),p(5,indx),p(4,indx))

      do ix=k10,k11
        k(1,ix)=12
      end do
      ncup=ncup+1
      ind(ncup)=indx
      if(ipart1.eq.i.or.ipart1.eq.i+1) ipart1=indx
      if(ipart2.eq.i.or.ipart2.eq.i+1) ipart2=indx
      goto 3400
3440  continue


c=======================================================================
c...Collision are sucessful.

      mste(5)=j_jet
      jflg=1

c...Count jet number
      mstd(55)=mstd(55)+j_jet

c...Update collision counter
      mstd(29)=mstd(29)+1

c...Find leading particle.
      do jm=1,2
        if(jm.eq.1)ip=ipart1
        if(jm.eq.2)ip=ipart2
        if(ip.ne.0) then
          lead(jm)=ip
c         k10=k(10,ip)
c         k11=k(11,ip)
c         lead(jm)=k11
c         kfa=abs(k(2,k10))
c         if(kfa.ge.1000.and.mod(kfa/10,10).eq.0)lead(jm)=k10
        endif
      end do

c...Remove particle and collision list.
      call jamzero(jp)
      call jamzero(jt)
      k(1,jp)=51
      k(1,jt)=51

c...Remove collision
      if(mstc(6).ge.0) then
        call jamcupda(jp,-1,0)
        call jamcupda(jt,-1,0)
      end if

      do jm=1,2
        jpq1=jpq(jm,1)
        jpq2=jpq(jm,2)
        if(jpq1.ne.0.and.jpq2.ne.0)then
          do i=jpq1,jpq2
            k(1,i)=33
            if(mstc(6).ge.0) call jamcupda(i,-1,0)
          end do
        endif
      end do



c...Update parton number
      do i=1,5
       psum(i)=0.0d0
      end do
      do 5400 j=1,ncup
        i=ind(j)
        if(k(1,i).gt.10) goto 5400
        kfa=abs(k(2,i))
        kflc=mod(kfa/10,10) 
        psum(1)=psum(1)+p(1,i)
        psum(2)=psum(2)+p(2,i)
        psum(3)=psum(3)+p(3,i)
        psum(4)=psum(4)+p(4,i)

        if(kfa.lt.10) then
          mstd(83)=mstd(83)+1
        else if(kfa.gt.1000.and.kflc.eq.0) then
          mstd(83)=mstd(83)+2
        else if(kfa.eq.21) then
          mstd(84)=mstd(84)+1

c....Sort index of baryons
        else if(mstc(6).ge.0.and.abs(k(9,i)).eq.3) then
          if(k(1,jp).gt.10) then
            ind(j)=jp
            call jamexch(jp,i)
          else if(k(1,jt).gt.10) then
            ind(j)=jt
            call jamexch(jt,i)
          else
            mdel=0
            call jamindb(indx,ind,0,mdel,idel)
            call jamcupda(idel(1),-1,1)
            call jamexch(indx,i)
            ind(j)=indx
          endif
          call jamzero(i)
          k(1,i)=41
        endif
 5400 continue

c...Update collision
      if(mstc(6).ge.0.and.ncup.ge.1.and.mstc(52).ne.1
     $ .and.mstc(76).ge.1) then
        do i=1,ncup
          call jamcupda(ind(i),-4,1)
        end do
      endif

c...Print some information
      if(mstc(8).ge.2) then
        ih=mstc(38)

        if(mstc(8).ge.3) then
          write(ih,*)'**** JAMHRD time mste(5) mstd(55)',
     $            pard(1),mste(5),mstd(55),jp,jt
         write(ih,*)'quark gluon',mstd(83),mstd(84)
         write(ih,*)'HARD:'

         do j=1,2
         if(j.eq.1)i=jp
         if(j.eq.2)i=jt
         write(ih,'(2(1x,i5),1x,i6,1x,2i5,6(1x,g10.3))')
     $              i,(k(m,i),m=1,2),k(10,i),k(11,i),p(5,i)
     $              ,(p(m,i),m=1,5)
         enddo

         do j=1,ncup
           i=ind(j)
         write(ih,'(2(1x,i5),1x,i6,1x,2i5,6(1x,g10.3))')
     $              i,(k(m,i),m=1,2),k(10,i),k(11,i),p(5,i)
     $              ,(p(m,i),m=1,5)
         end do
         write(ih,*)'lead',lead(1),lead(2),ipart1,ipart2
        endif

         srt1=psum(4)**2-(psum(1)**2+psum(2)**2+psum(3)**2)
         if(srt1.le.0d0)
     $    call jamerrm(30,0,'srt<0 after hard scatt')
          srt1=sqrt(srt1)

         if(   (abs(pcm(1)-psum(1)).gt.0.1d0)
     $       .or.(abs(pcm(2)-psum(2)).gt.0.1d0)
     $         .or.(abs(pcm(3)-psum(3)).gt.0.1d0)
     $           .or.(abs(pcm(4)-psum(4)).gt.0.1d0)
     $             .or.(abs(pcm(5)-srt1).gt.0.1d0) ) then
            write(ih,*)'pcm  ',(pcm(l),l=1,5)
            write(ih,*)'psum ',(psum(l),l=1,4),srt1
            call jamerrm(30,0,'momentum not conserved after hard scatt')
         endif

         call jamcheck('<<After jamhrd>>')

      endif

c...Fragment strings.
      if(mstc(76).eq.0.and.nlq.ge.1) then
c       do ii=1,nlq
        do ii=nv0+1,nv
          if(k(1,ii).eq.4) then
            k10=k(10,ii)
            k11=k(11,ii)
            ilead1=0
            ilead2=0
            if(k10.eq.lead(1).or.k11.eq.lead(1)) ilead1=1
            if(k10.eq.lead(2).or.k11.eq.lead(2)) ilead2=1
            call jamsave(1,1,ip)
            call jamdec(ip,indd,nadd,icon)
            do i=1,nadd
             if(mod(abs(k(1,indd(i)))/10,10).eq.2) then
               if(ilead1.eq.1) lead(1)=indd(i)
               if(ilead2.eq.1) lead(2)=indd(i)
             endif
             if(mstc(6).ge.0) call jamcupda(indd(i),-1,1)
            end do
          endif
        end do 
      endif

      end

c***********************************************************************

      subroutine jamhrdv

c...Extract space-time information for hard parton-parton scattering.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      real*8 jamprtim
c....Pythia common block.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      save /jyjets/

      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      save /pjpars/,/pjint1/,/pjint2/

      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      common/hijcrdn/yp(3,300),yt(3,300)
      save  /hijcrdn/,/hiparnt/

      parameter(mxpart=300)
      common/partdis1/rpart(mxpart,5),ppart(mxpart,5),kpart(mxpart,5)
      common/partdis2/tcoll0,npart0,npart1,npart2,npart3,npart4
      common/partdis3/xorg(2,4)
      save /partdis1/,/partdis2/

c...Multiple scattering.
      common/jampyda1/qmult(10),imult(2,10)
      save /jampyda1/

      common/jamxhit/xo(2,5),srtc,ecmc

      dimension q2s(2),more(2)
      logical first
      data first/.true./
c     data tcoll0/0.0d0/


      if(first) then
c       call jamanap(1,0.0)
        tcoll0=0d0
        first=.false.
      endif

      x1=vint(141)
      x2=vint(142)
      do j=1,4
      xorg(1,j)=xo(1,j)
      xorg(2,j)=xo(2,j)
      end do

c...Check scattered partion is valence or not.
      do jt=1,2
       if(jt.eq.1) then
         kfb=mint(11)
         kfq=mint(13)
       else if(jt.eq.2) then
        kfb=mint(12)
        kfq=mint(14)
      endif
      j=abs(kfb)
      isg=isign(1,kfb)
      kfv1=mod(j/10,10)
      kfv2=mod(j/100,10)
      kfv3=mod(j/1000,10)
      if(kfv3.eq.0) then
        call attflv(kfb,kfv1,kfv2)
      else
        kfv1=isg*kfv1
        kfv2=isg*kfv2
        kfv3=isg*kfv3
      endif
      if(kfq.eq.kfv1.or.kfq.eq.kfv2.or.kfq.eq.kfv3) then
      else
c....Assign contracted distribution to the sea (anti)quarks
c...from uncertainty principle.
c...x,y direction is not yet implemented now!
        if(jt.eq.1) then
          ddz=2.d0*paru(3)/(x1*srtc)
          dz1=-ddz/2+rn(0)*ddz
          xorg(1,3)=xorg(1,3)+dz1
        else if(jt.eq.2) then
          ddz=2.d0*paru(3)/(x2*srtc)
          dz2=-ddz/2+rn(0)*ddz
          xorg(2,3)=xorg(2,3)+dz2
        endif
      endif
      end do


      q2mnc=parp(62)**2

c...Get gamma factor
      dbx=vint(8)
      dby=vint(9)
      dbz=vint(10)
      dga=ecmc/srtc
      db=sqrt(dbx**2+dby**2+dbz**2) 
      eps1=1d0-1d-12
c...Rescale boost vector if too close to unity.
      if(db.gt.eps1) then
        call jamerrm(3,0,'(jamhrd:) boost vector too large')
        dbx=dbx*(eps1/db)
        dby=dby*(eps1/db)
        dbz=dbz*(eps1/db)
        db=eps1
        dga=1d0/sqrt(1d0-db**2) 
      endif
      rtime=jamprtim(1,0,0.0d0,dga)

c....Calculate space-time point for space-like parton branching.
c...           ________
c... -------i3|        |i7--------
c... -------i4|________|i8--------
c
c     q2mx=vint(56)
c     isub=mint(1)
c     if(iset(isub).eq.2) q2mx=parp(67)*vint(56)

c....Reset vector and initialization.
      do i=mint(84)+1,n
      do j=1,5
        v(i,j)=0.0d0
      end do
      end do

      do jt=1,2
       more(jt)=1
       if(jt.eq.1) i=mint(84)+5
       if(jt.eq.2) i=mint(84)+6
       q2s(jt)=p(i,5)**2
      end do

c     v(mint(84)+1,4)=tcoll(1)
c     v(mint(84)+2,4)=tcoll(2)
c     do j=1,3
c       v(mint(84)+1,j)=xorg(1,j)+(xo(1,4)-tcoll0)*vz1
c       v(mint(84)+2,j)=xorg(2,j)+(xo(2,4)-tcoll0)*vz2
c     end do
      do j=1,4
        v(mint(84)+1,j)=xorg(1,j)
        v(mint(84)+2,j)=xorg(2,j)
      end do

c...Loop for initial state radiation.
      i=mint(84)+4
      ntry=0
100   continue

      ntry=ntry+1
      if(ntry.gt.n) then
        mstu(11)=8
        call pjlist(4)
        do i=mint(84)+1,n
         write(mstc(38),*)'Infinit loop??(spacelike)',
     $ i,(k(i,j),j=1,3),p(i,5)
        end do
        call jamerrm(30,0,'(jamhrdv:) Infinit looop at spacelike')
      endif

      jt=1
      i=i+1
      if(ntry.ne.1.and.q2s(2).gt.q2s(1) ) jt=2
      if(more(jt).eq.0) jt=3-jt

      ip1=mod(k(i,4),10000)
      ip2=mod(k(i,5),10000)

c.....Timelike partion branch
      if(k(i,1).eq.14.and.ip1.eq.ip2.and.ip1.eq.i+1) then
        nop1=1
        n1=i+1
        ip=min(mod(k(i,5)/10000,10000),mod(k(i,4)/10000,10000))
        im=k(i,3)
        gam=p(ip,4)/abs(p(ip,5))
c       v(i,4)=v(ip,4)
        v(i,4)=v(ip,4)-jamprtim(2,k(im,2),abs(p(im,5)),gam)
        v(i,5)=v(i,4)+p(i,4)/p(i,5)**2*paru(3)
        do j=1,3
         v(i,j)=v(ip,j)+p(ip,j)/p(ip,4)*(v(ip,4)-v(ip,5))
         v(n1,j)=v(i,j)
        end do
        v(n1,4)=v(i,4)
        v(n1,5)=v(i,5)
        n2=max(mod(k(i,5)/10000,10000),mod(k(i,4)/10000,10000))
        do 110 it=n1,n2-1
            if(k(it,1).le.13) goto 110
            it1=mod(k(it,4),10000)
            it2=mod(k(it,5),10000)
                if(it1.le.0.or.it2.le.0) then
                   mstu(11)=8
                   call pjlist(4)
                   call jamerrm(30,0,'(jamhrdv:) at timelike')
                endif
            v(it1,4)=v(it,5)   ! Production time.
            v(it2,4)=v(it,5)   ! Production time.
c........Production point.
            do j=1,3
              v(it1,j)=v(it,j)+p(it,j)/p(it,4)*(v(it,5)-v(it,4))
              v(it2,j)=v(it1,j)
            end do

c.........Dead time
            if(k(it1,1).eq.14) then
              gam=p(it1,4)/abs(p(it1,5))
              v(it1,5)=v(it1,4)+jamprtim(3,k(it1,2),abs(p(it1,5)),gam)
            else
              v(it1,5)=1.d+35
            endif
            if(k(it2,1).eq.14) then
              gam=p(it2,4)/abs(p(it2,5))
              v(it2,5)=v(it2,4)+jamprtim(3,k(it2,2),abs(p(it2,5)),gam)
            else
              v(it2,5)=1.d+35
            endif
110     continue
        v(i,4)=-1000.0d0
        v(i,5)=-1000.0d0
        i=n2
        if(i.le.n1+1) then
          write(mstc(38),*)'n1+1 n2 iorg',n1+1,i
          call jamerrm(30,0,'(jamhrdv:) 110')
        endif
      endif

c....Spacelike partion branching
      im=k(i,3)            ! Mother
      ip=min(ip1,ip2) ! daughter

      if(k(i,1).eq.14) then
        v(i,5)=v(ip,4)       ! dead time.
        if(im.eq.3.or.im.eq.4) then  ! initiator
          nop=2
          more(jt)=0
          q2s(jt)=0.0d0
c         v(i,4)=v(i,5)-p(i,4)/q2mnc*paru(3)
          v(i,4)=tcoll0     ! production time
          if(im.eq.3) jm=1
          if(im.eq.4) jm=2
          do j=1,3
            v(i,j)=xorg(jm,j)
          end do
        else
          nop=3
          q2s(jt)=p(i,5)**2
          gam=p(i,4)/abs(p(i,5))
          ptime=v(i,5)-jamprtim(2,k(i,2),abs(p(i,5)),gam)
          if(ptime.le.tcoll0) then
c            write(10,*)'Production time E Q',i,ptime,p(i,4),p(i,5)
             ptime=tcoll0
          endif
          v(i,4)=ptime ! Production time
          do j=1,3
            v(i,j)=v(ip,j)+p(ip,j)/p(ip,4)*(v(ip,4)-v(ip,5))
          end do
        endif
      else if(k(i,1).eq.13) then !produced particle i.e.will be on shell
        nop=4
        iq=min(mod(k(i,4)/10000,10000),mod(k(i,5)/10000,10000))
        if(iq.eq.0.or.iq.ge.i) then
              mstu(11)=8
              call pjlist(4)
              write(mstc(38),*)'space-like i iq=',
     $               i,iq,k(i,1),k(i,2),k(i,3),p(i,5),jt
              call jamerrm(30,0,'(jamhrdv:) space-like i iq')
        endif
        v(i,5)=1.d+35
        v(i,4)=v(iq,4)
        do j=1,3
          v(i,j)=v(iq,j)+p(iq,j)/p(iq,4)*(v(iq,4)-v(iq,5))
        end do
      else
          mstu(11)=8
          call pjlist(4)
          write(mstc(38),*)'invalid i k1=',i,k(i,1)
          call jamerrm(30,0,'(jamhrdv:) invalid i k1')
      endif

      if(more(1).eq.1.or.more(2).eq.1) goto 100

150   continue

c...Time-like branch
      icms=0
      irem=0
      do i=mint(84)+1,n
        if(k(i,2).eq.94) then
         icms=i
        endif
        if(k(i,3).eq.1.or.k(i,3).eq.2) then
         irem=i-1
         goto 200
        endif
      end do
200   continue

      npart2=icms+1
      npart3=irem
c...No time-like branching case.
      if(icms.eq.0) then
        k7=mint(84)+3
        k8=mint(84)+4
        v(k7,4)=rtime+xorg(1,4)
        v(k8,4)=rtime+xorg(2,4)
        v(k7,5)=1d+30
        v(k8,5)=1d+30
        do j=1,3
          v(k7,j)=xorg(1,j)+rtime*p(k7,j)/p(k7,4)
          v(k8,j)=xorg(2,j)+rtime*p(k8,j)/p(k8,4)
        end do
c       npart1=mint(84)+1
        npart1=k7
        goto 1000
      endif
      npart1=mint(84)+5
      v(icms,4)=-1000.d0
      v(icms,5)=-1000.d0

      i=icms+1
      v(i,4)=rtime+xorg(1,4)
      if(k(i,1).eq.14) then
        v(i,5)=v(i,4)+jamprtim(3,k(i,2),abs(p(i,5)),gam)
      else
        v(i,5)=1.d+20
      endif
      do j=1,3
        v(i,j)=xorg(1,j)
      end do

      i=icms+2
      v(i,4)=rtime+xorg(2,4)
      if(k(i,1).eq.14) then
        v(i,5)=v(i,4)+jamprtim(3,k(i,2),abs(p(i,5)),gam)
      else
        v(i,5)=1.d+20
      endif
      do j=1,3
c       v(i,5)=v(i,4)+p(i,4)/p(i,5)**2*paru(3)
        v(i,j)=xorg(2,j)
      end do

      n1=icms+1
      do 210 i=n1,irem
        if(k(i,1).le.13) goto 210
        j1=mod(k(i,4),10000)
        j2=mod(k(i,5),10000)

        if(p(i,5).lt.0.0001d0) then
           write(mstc(38),*)'time like Q im=',
     $            im,k(im,1),k(im,2),k(im,3),p(im,5),i
           mstu11=mstu(11)
           mstu(11)=mstc(38)
           call pjlist(2)
           mstu(11)=mstu11
           call jamerrm(30,0,'(jamhrdv:) p5<0')
        endif

c......Production time.
        v(j1,4)=v(i,5)
        v(j2,4)=v(i,5)

c....Production point.
        do j=1,3
          v(j1,j)=v(i,j)+p(i,j)/p(i,4)*(v(i,5)-v(i,4))
          v(j2,j)=v(j1,j)
        end do

c....Dead time.
        if(k(j1,1).eq.14) then
           gam=p(j1,4)/abs(p(j1,5))
           v(j1,5)=v(j1,4)+jamprtim(3,k(j1,2),abs(p(j1,5)),gam)
        else
          v(j1,5)=1.d+35
        endif
        if(k(j2,1).eq.14) then
          gam=p(j2,4)/abs(p(j2,5))
          v(j2,5)=v(j2,4)+jamprtim(3,k(j2,2),abs(p(j2,5)),gam)
        else
          v(j2,5)=1.d+35
        endif
210   continue

1000  continue


c...Calculate reaction time for multiple scattering.
      if(mint(31).ge.2) then
c       write(mstu(11),*)'Multiple int'
        do i=1,mint(31)-1
         i1=imult(1,i)
         i2=imult(2,i)
         qq=qmult(i)
         if(i1.eq.0.or.i2.eq.0) then
            write(mstc(38),*)'error (mult:)',i1,i2
            call jamerrm(30,0,'(jamhrdv:) i1 i2=0')
         endif
         if(v(i1,4).ne.0.0d0.or.v(i2,4).ne.0.0d0) then
            write(6,*)'error (mult:)v',i1,i2
         endif
         rtim=jamprtim(4,0,qq,gam)
         v(i1,4)=rtim+xorg(1,4)
         v(i2,4)=rtim+xorg(2,4)
         do j=1,3
          v(i1,j)=xorg(1,j)+rtim*p(i1,j)/p(i1,4)
          v(i2,j)=xorg(2,j)+rtim*p(i2,j)/p(i2,4)
         end do
c         write(22,*)i1,i2,qq,(v(i1,j),j=1,4),(v(i2,j),j=1,4)
c         write(mstu(11),*)i1,(v(i1,j),j=1,4)
c         write(mstu(11),*)i2,(v(i2,j),j=1,4)
        end do
      endif

      npart0=1
      npart4=n
      do i=1,n
        do j=1,5
          rpart(i,j)=v(i,j)
          ppart(i,j)=p(i,j)
          kpart(i,j)=k(i,j)
        end do
      end do

c...Copy results to final entry.
      do i=npart3+1,n
       im=k(i,3)
c...Remnant.
       im2=k(im,3)
       if(im2.eq.1.or.im2.eq.2) then
         do j=1,4
         v(i,j)=xorg(im2,j)
         end do
       else
        do j=1,5
         v(i,j)=v(im,j)
        end do
       endif
      enddo

cc    call jamanap(2,0.0)

      end

c***********************************************************************

      double precision function jamprtim(is,kf,qv,gam)

c....Calculate reaction time for partonic time evolutions.
c..is=1: parton-parton collision
c..is=2: space-like
c..is=3: time-like
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      common/pjint1/mint(400),vint(400)
      save /pjint1/
      dimension wdtp(0:200),wdte(0:200,0:5)
      parameter(parv44=0.1d0,parv63=0.1d0,parv64=1.0d0)
      parameter(mstv63=1)

c...Life time of space-like branch.
      if(is.eq.2) then
          if(qv.lt.1d-6) then
            jamprtim=1.0d0
            return
          endif
          jamprtim=min(1.0d0,parc(74)*gam/qv*paru(3))
c     else
c     tau=parc(74)/wid*paru(3)
c     jamprtim=-tau*log(rn(0))
c     xran=rn(0)
c     jamprtim=-tau*sqrt(xran/(1.-xran))
c     endif

        return
      endif

c...Life time for time-like branching.
      if(is.eq.3) then
        wpar=parc(73)
c       pms=min(qmax,ps(5))**2
        pms=qv**2
        if(mstj(44).eq.0) then
          alps=paru(111)
        else
          q2=max(paru(112),0.25d0*pms)
          alps=pjalem(q2)
        end if
        ig=0
        if(kf.eq.21) ig=1

c...Reaction time for parton-parton collision.
      else if(is.eq.1.or.is.eq.4) then

        wpar=parc(72)
c.....Multiple scattering.
        if(is.eq.4) then
          pms=qv*qv
          alps=vint(58)
          fps=1.d0
          fde=1.d0
          goto 100
        endif

        isub=mint(1)
        pms=vint(52)
        ig=0
        if(isub.eq.12.or.isub.eq.53) ig=1
        if(isub.eq.14.or.isub.eq.18.or.isub.eq.29) then
          alps=vint(57)
        else
          alps=vint(58)
        endif
      endif

      fps=1.d0
      if(ig.eq.1) then
        fde=2.d0
        wdtps=0.d0
        call pjwidt(21,pms,wdtp,wdte)
        do 275 i=1,mstj(45) 
        wdtps=wdtps+max(0.d0,wdtp(i))
  275   if(pms.gt.wdtps) fps=fps+1.d0
      else
        fde=2.d0/3.d0
      end if

100   continue
      width=fde*fps*alps*sqrt(pms)
      if(width.lt.1.d-10) then
       jamprtim=1d+30
       return
      endif

c...Sample life-time from probability distribution.
      taur=1.d0/width
      taul=gam*taur
      jamprtim=0.d0
      rd=rn(0)
      if(mstc(85).eq.1) jamprtim=-wpar*taul*alps*dlog(rd)
      if(mstc(85).eq.2) jamprtim=wpar*taul*alps*sqrt(rd/(1.d0-rd))

      end
