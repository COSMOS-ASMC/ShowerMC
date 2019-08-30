c the use
c file in order to avoid unnecessary long compilation times
c$$$$$$$$$$$$$$$$$$$$$$$
c       input parameter usage changed:
c             iproty=0 for nucleus is unchanded
c             iproty>0 is particle code listed in jetset6.2/6.3
c             iproty<0 is antiptcl                //
c$$$$$$$$$$$$$$$$$$$$$$$
c***********************************************************************

      subroutine ingeboC

c-----------------------------------------------------------------------
c-----administrates one event-------------------------------------------
c-----------------------------------------------------------------------
      common/indataC/elab,rots,nap,nzp,r0p,nat,nzt,r0t,iflspv,bmin,
     *                bmax,neve,iproty,ifermi,iflout
      common/evevecC/nevent,isppp,isppn,isptp,isptn,bimp,
     *                idi(2000),ipr(2000)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/ludat1C/mst(40),par(80)
      common/frstmaC/win(40),wta(200)
      common/frsttyC/idp(40),idt(250)
      common/frnynuC/bipa,nwp,nwt,ny,ni(500),nt(500)
      common/frstptC/pxsp(40),pysp(40),pxst(200),pyst(200)
      common/frstpmC/ppin(40),pmin(40),ppta(200),pmta(200)

c-----setting vectors to zero-------------------------------------------
c          kkkkkkkkkkkk
      integer discard
c       kkkkkkkkkkkkk
      data nevent,iflw/0,0/
c $$$$$$$$$ need not be cleared
c$$   do 1 i=1,2000
c$$     ipr(i)=0
c$$ 1 continue

      do 5 i=1,40
      pxsp(i)=0.
      pysp(i)=0.
    5 continue
      do 6 i=1,200
      pxst(i)=0.
      pyst(i)=0.
    6 continue
      n=0
      mst(12)=0
      mst(23)=0
      if (iflout.eq.1) mst(7)=0
c$    if(iflw.eq.0) then
c$    write(mst(20),770) 'the lund monte carlo - fritio version 1.6'
c$    write(mst(20),771) 'last date of change : 10 june  1986'
c$    iflw=1
c$770 format(' ',19x,a)
c$771 format(' ',22x,a)
c$    endif

c-----hildinC for identification of the participating nucleons----------a
      call hildinC

c-----angantC for creation of the nuclei and the subsequent nucleon----aa
c-----collisions--------------------------------------------------------
      call angantC

c-----ring for giving the excited nucleons masses-----------------------

c         kkkkkkkkkkk
      call ringC(discard)
      if(discard .ne. 0) then
         n = 0
         return
      endif
c         kkkkkkkkkkk

c-----updating of the event number, recording of impact parameter and---
c-----counting of spectators--------------------------------------------
      nevent=nevent+1
      bimp=bipa
      jcopp=0
      do 3 i=1,nwp
      if(idp(i).eq.41) jcopp=jcopp+1
    3 continue
      jcotp=0
      do 4 i=1,nwt
      if(idt(i).eq.41) jcotp=jcotp+1
    4 continue

      isppp=nzp-jcopp
      isppn=nap-nwp-isppp
      isptp=nzt-jcotp
      isptn=nat-nwt-isptp

      if(iproty.ne.0.or.nzp.le.-1) then
      isppp=0
      isppn=0
      endif

c-----------------------------------------------------------------------
c-----projectile loop. beleC:quarkflavours in string ends.---------------
c----------------------torsteC:pt for the gluon & way of fragmentation--a
      do 100 j=1,nwp
      call beleC(ifla,iflb,idp(j))
      call torsteC(win(j),idp(j),iflb,ifla,x1,x3,kt)
      nsav=n
      if(kt.eq.1) then
        call angurvC(1,j)
        goto 100
      elseif(kt.eq.2) then
        call lu2jetC(nsav+1,iflb,ifla,win(j))
        call luexecC
        mst(1)=nsav+1
      elseif(kt.eq.3) then
        call lu3jetC(nsav+1,iflb,ifla,win(j),x1,x3)
        call luexecC
        mst(1)=nsav+1
c-------pt--------------------------------------------------------------
c-------orientation of the threejet event-------------------------------
        x2=2.-x1-x3
c-------quarkjet axis minimizing pt ------------------------------------
        th=atan2(p(nsav+3,1),p(nsav+3,3))
        thp=3.1416-th
        psi=0.5*atan((x3*x3*sin(2*thp))/(x1*x1+x3*x3*cos(2*thp)))
        call luroboC(psi,0.,0.,0.,0.)
c-----------------------------------------------------------------------
        vink=2.*3.14159*rluC(0)
        call luroboC(0.,vink,0.,0.,0.)
      endif
c-----------------------------------------------------------------------
c-----boost to lab and pt for the whole string system ------------------
      beinx=2.*pxsp(j)/(ppin(j)+pmin(j))
      beiny=2.*pysp(j)/(ppin(j)+pmin(j))
      embein=2.*pmin(j)/(ppin(j)+pmin(j))
      call luroinC(0.,0.,beinx,beiny,embein)

      if(iflout.eq.3) then
        call lueditC(3)
      elseif(iflout.eq.2) then
        call lueditC(2)
      endif

c-----produced particles in fragmentation-------------------------------
      do 110 i=nsav+1,n
        ipr(i)=1
        idi(i)=0
  110 continue
      mst(1)=0
  100 continue

c-----------------------------------------------------------------------
c-----target loop  beleC: quarkflavours in string ends-------------------
c------------------torsteC: pt for the gluon & way of fragmentation-----a
      do 200 j=1,nwt
      call beleC(ifla,iflb,idt(j))
      call torsteC(wta(j),idt(j),ifla,iflb,x1,x3,kt)
      nsav=n
      if(kt.eq.1) then
        call angurvC(2,j)
        goto 200
      elseif(kt.eq.2) then
        call lu2jetC(nsav+1,ifla,iflb,wta(j))
        call luexecC
        mst(1)=nsav+1
      elseif(kt.eq.3) then
        call lu3jetC(nsav+1,ifla,iflb,wta(j),x1,x3)
        call luexecC
        mst(1)=nsav+1
        thet3=3.1416-pluC(nsav+3,13)
        call luroboC(thet3,0.,0.,0.,0.)
c-------pt global-------------------------------------------------------
c-------orientation of the threejet event ------------------------------
        x2=2.-x1-x3
c-------quarkjet axis minimizing pt ------------------------------------
        th=abs(atan2(p(nsav+1,1),p(nsav+1,3)))
        tht=th
        psi=0.5*atan((x1*x1*sin(2*tht))/(x3*x3+x1*x1*cos(2*tht)))
        call luroboC(-psi,0.,0.,0.,0.)
c-----------------------------------------------------------------------
        vink=2.*3.14159*rluC(0)
        call luroboC(0.,vink,0.,0.,0.)
      endif
c-----------------------------------------------------------------------
c-----boost to lab and pt for the whole string system ------------------
      betax=2.*pxst(j)/(ppta(j)+pmta(j))
      betay=2.*pyst(j)/(ppta(j)+pmta(j))
      epbeta=2.*ppta(j)/(ppta(j)+pmta(j))
      embeta=2.*pmta(j)/(ppta(j)+pmta(j))
      if (rots.ne.0.) then
        call lurotaC(0.,0.,betax,betay,epbeta)
      else
        call luroinC(0.,0.,betax,betay,embeta)
      endif

      if(iflout.eq.3) then
        call lueditC(3)
      elseif(iflout.eq.2) then
        call lueditC(2)
      endif

c-----produced particles in fragmentation-------------------------------
      do 210 i=nsav+1,n
        ipr(i)=0
        idi(i)=0
  210 continue
      mst(1)=0
  200 continue

c-----preparing the event as determined by the iflout setting-----------
      if(iflout.eq.1.or.iflout.eq.0) then
        nlines=0
        lastje=0
        do 111 i=1,n
          if(abs(k(i,2)).ge.500) then
            lastje=i
          else
            nlines=nlines+1
            do 112 j=1,5
              p(nlines,j)=p(i,j)
  112       continue
            iori=mod(k(i,1),10000)
            if (iori.le.lastje) then
              iori=0
            else
              iori=nlines+iori-i
            endif
            k(nlines,1)=(k(i,1)/10000)*10000+iori
            k(nlines,2)=k(i,2)
            idi(nlines)=idi(i)
            ipr(nlines)=ipr(i)
          endif
  111   continue
        n=nlines
      endif
      return
      end


      subroutine angurvC(l,j)
c-----------------------------------------------------------------------
c     taking care of diffractive nucleons(and diff pions,kaons)
c-----------------------------------------------------------------------
c        $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c        in original version /indataC/ is missing and iflout becomes
c        undefined. to avoid that, the next 2 lines inserted.
      common/indataC/elab,rots,nap,nzp,r0p,nat,nzt,r0t,iflspv,bmin,
     *                bmax,neve,iproty,ifermi,iflout
c        $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      common/frstpmC/ppin(40),pmin(40),ppta(200),pmta(200)
      common/frsttyC/idp(40),idt(250)
      common/evevecC/nevent,isppp,isppn,isptp,isptn,bimp,
     *                idi(2000),ipr(2000)
      common/lujetsC/n,k(2000,2),p(2000,5)
      common/frstptC/pxsp(40),pysp(40),pxst(200),pyst(200)
      nsav=n
c-----projectile--------------------------------------------------------
      if(l.eq.1) then
        e=(ppin(j)+pmin(j))/2.
        pz=(ppin(j)-pmin(j))/2.
        px=pxsp(j)
        py=pysp(j)
        kkod=idp(j)
c-----target------------------------------------------------------------
      elseif(l.eq.2) then
        e=(ppta(j)+pmta(j))/2.
        pz=(ppta(j)-pmta(j))/2.
        px=pxst(j)
        py=pyst(j)
        kkod=idt(j)
      endif
c-----------------------------------------------------------------------
      n=n+1
      if(l.eq.1) ipr(n)=1
      if(l.eq.2) ipr(n)=0
      idi(n)=1
      k(n,1)=0
      k(n,2)=kkod
      p(n,1)=px
      p(n,2)=py
      p(n,3)=pz
      p(n,4)=e
      p(n,5)=ulmassC(0,kkod)
      if (iflout.eq.3.and.kkod.eq.42) n=n-1
      return
      end




      subroutine beleC(ifla,iflb,j)
c-----------------------------------------------------------------------
c     giving spin and quarkflavour to the leading particles
c-----------------------------------------------------------------------
      spin=rluC(0)
c        pi
      if(j.eq.17) then
        if(spin.lt..5) then
          ifla=1
          iflb=-2
        else
          ifla=-2
          iflb=1
        endif
      elseif(j.eq.-17) then
        if(spin.lt..5) then
          ifla=-1
          iflb=2
        else
          ifla=2
          iflb=-1
        endif
c         k
      elseif(j.eq.18) then
        if(spin.lt..5) then
          ifla=1
          iflb=-3
        else
          ifla=-3
          iflb=1
        endif
      elseif(j.eq.-18) then
        if(spin.lt..5) then
          ifla=-1
          iflb=3
        else
          ifla=3
          iflb=-1
        endif
      elseif(j .eq. 19 .or. j .eq. 37 .or. j .eq. 38) then
c          k0
        if(spin .le. .5) then
          ifla=2
          iflb=-3
        else
          ifla=-3
          iflb=2
        endif
      elseif(j .eq.-19 .or. j .eq.-37 .or. j .eq.-38) then
c          k0_bar
        if(spin .le. .5) then
          ifla=-2
          iflb=3
        else
          ifla=3
          iflb=-2
        endif
      elseif(j.eq.41) then
c          p
        if(spin.lt.0.167) then
          ifla=1
          iflb=21
        elseif(spin.lt.0.5) then
          ifla=2
          iflb=11
        else
          ifla=1
          iflb=12
        endif
      elseif(j.eq.42) then
c          n
        if(spin.lt.0.167) then
          ifla=2
          iflb=21
        elseif(spin.lt.0.5) then
          ifla=1
          iflb=22
        else
          ifla=2
          iflb=12
        endif
      elseif(j.eq.-41) then
c           p_bar
        if(spin.lt.0.167) then
          ifla=-1
          iflb=-21
        elseif(spin.lt.0.5) then
          ifla=-2
          iflb=-11
        else
          ifla=-1
          iflb=-12
        endif
      elseif(j.eq.-42) then
c          n_bar
        if(spin.lt.0.167) then
          ifla=-2
          iflb=-21
        elseif(spin.lt.0.5) then
          ifla=-1
          iflb=-22
        else
          ifla=-2
          iflb=-12
        endif
      elseif(abs(j) .eq. 57  .or.( j .ge.43 .and. j .le. 45))then
c            lamda; sigma is treated as if lambda
          if(spin .lt. .167) then
              ifla=1
              iflb=23
          elseif(spin .le. .333) then
              ifla=1
              iflb=32
          elseif(spin .le. .5) then
              ifla=2
              iflb=13
          elseif(spin .le. .667) then
              ifla=2
              iflb=31
          elseif(spin .le. .833) then
              ifla=3
              iflb=21
          else
              ifla=3
              iflb=12
          endif
      elseif(abs(j) .eq. 23 .or.(  abs(j) .ge. 33
     *     .and. abs(j) .le. 35)) then
c          pi0 and rho 0/ omega ,fai
        u = rluC(0)
        if(u .lt. .5)then
c             uu_bar
           if(spin.lt. .5) then
              ifla=1
              iflb=-1
           else
              ifla=-1
              iflb=1
           endif
        else
c             d d_bar
           if(spin.lt. .5) then
              ifla=2
              iflb=-2
           else
              ifla=-2
              iflb=2
           endif
        endif
      endif

      return
      end




      subroutine torsteC(w,id,izpl,izmi,x1,x3,kt)
c-----------------------------------------------------------------------
c      twojet, threejet or diffractive nucleon. if threejet the routine
c      gives pt to the gluon
c-----------------------------------------------------------------------
      data alit /6./
c$$$$$$$$$$$$$$$$$$$$$$
      idabs=abs(id)
      if(idabs .eq. 17 .or. idabs .eq. 23) then
c           pi
           framin=0.8
      elseif(idabs .eq. 41 .or. idabs .eq. 42) then
c           n,p
           framin=1.2
      elseif(idabs .eq. 18) then
c           k+-
           framin=1.
      elseif(idabs .eq. 19 .or. idabs .eq. 37 .or. idabs .eq. 38) then
c           k 0
           framin=1.
      elseif(idabs .eq. 57) then
c             lambda
           framin=1.e37
      elseif(idabs .ge. 43 .and. idabs .le. 46) then
c            sigma
           framin=1.e37
      elseif(idabs .eq. 33 .or. idabs .eq. 34) then
c             rho,omega
           framin=1.e37
      elseif(idabs .eq. 35) then
c            fai
           framin=1.e37
      else
           write(*,*) ' ptcl id=',id, ' not supported in torsteC'
           stop
      endif
c     if(abs(id).eq.17) framin=0.8
c     if(abs(id).eq.18) framin=1.0
c     if(abs(id).eq.41.or.abs(id).eq.42) framin=1.2

      if(w.le.framin-0.001) then
        kt=1
      else
        slu=rluC(0)
c$$$$$$$$$$$$$$$$ originally pete=sqrt( .... ) is defined and used
        pete2= alit*((w**2/(16.*alit)+1.)**slu-1.)
        if(pete2.lt.0.25) then
          kt=2
        else
          pete=sqrt(pete2)
          kt=3
          iflag=0
          fakt=rluC(0)
          y=2.*(fakt-.5)*log(w/(2.*pete)+sqrt(w**2/(4.*pete2)-1.))
c$$$$$$$
          x1=1.-(pete/w)*exp(y)
          x3=1.-(pete/w)*exp(-y)
          x2=2.-x1-x3
c                               u,d
          if(abs(izpl).le.2)  zpmass=0.325
c                               s
          if(abs(izpl).eq.3)  zpmass=0.500
c                               c
          if(abs(izpl).eq.4)  zpmass=1.6
c                               uu1 etc
          if(abs(izpl).gt.10) zpmass=0.650

          if(abs(izmi).le.2)  zmmass=0.325
          if(abs(izmi).eq.3)  zmmass=0.500
          if(abs(izmi).gt.10) zmmass=0.650

c---------kinematical tests---------------------------------------------
          pa1=(0.5*x1*w)**2-zpmass**2
          pa3=(0.5*x3*w)**2-zmmass**2
          pa2=0.5*(2.-x1-x3)*w

          if(pa1 .le.0.0 .or. pa3 .le. 0.0  .or. pa2 .le. 0.) then
            iflag=1
            goto 13
          endif

          arg1=(pa3-pa1-pa2**2)/(2.*sqrt(pa1)*pa2)
          arg2=(pa2**2-pa1-pa3)/(2.*sqrt(pa1*pa3))
          if(abs(arg1).gt.1.0.or.abs(arg2).gt.1.0) iflag=1

   13     if(iflag.eq.1) kt=2

        endif
      endif
      return
      end




      subroutine ringC(discard)
c-----------------------------------------------------------------------
c     the routine gives masses to the excited nucleons
c-----------------------------------------------------------------------
      common/indataC/elab,rots,nap,nzp,r0p,nat,nzt,r0t,iflspv,bmin,
     *                bmax,neve,iproty,ifermi,iflout
      common/frstmaC/win(40),wta(200)
      common/frstptC/pxsp(40),pysp(40),pxst(200),pyst(200)
      common/frsttyC/idp(40),idt(250)
      common/frstpmC/ppin(40),pmin(40),ppta(200),pmta(200)
      common/frnynuC/bipa,nwp,nwt,ny,ni(500),nt(500)
      logical cms,lab
      dimension fmp(40),fmt(200)
      double precision dptsq,dmassq
c-----------------------------------------------------------------------
c    kkkkkkkkkkk
      integer endless, discard
      endless = 0
      discard = 0
c   kkkkkkkkkkkkkkk
   31 if(elab.eq.0.) then
        cms=.true.
        lab=.false.
      elseif(rots.eq.0.) then
        cms=.false.
        lab=.true.
      endif
c-----------------------------------------------------------------------
      do 6 i=1,nwt
        if (idt(i).eq.41) then
          fmt(i)=.93828
        else
          fmt(i)=.93957
        endif
    6 continue
c     kkkkkkkkkkkkkkkk
      endless = endless + 1
      if(endless .gt. 10) then
         discard = 1
         return
      endif
c     kkkkkkkkkkkkkkkkk
      wtamin=1.2
      wtdmin=0.94

      if(abs(idp(1)).eq.17) then
        fmp(1)=0.140
        winmin=1.0
        widmin=0.75
      elseif(abs(idp(1)).eq.18) then
        fmp(1)=0.495
        winmin=1.1
        widmin=0.85
      else
        do 5 i=1,nwp
          if (abs(idp(i)).eq.41) then
            fmp(i)=.93828
          else
            fmp(i)=.93957
          endif
    5   continue
        winmin=1.2
        widmin=0.94
      endif

      if(cms) then
        do 10 i=1,nwp
          plongp=rots*rots/4-(fmp(i)*fmp(i)+fmt(1)*fmt(1))/2+
     *    (fmp(i)*fmp(i)-fmt(1)*fmt(1))**2/(4*rots*rots)
          esen=sqrt(plongp+fmp(i)*fmp(i))
          ppin(i)=esen+sqrt(plongp)
          pmin(i)=fmp(i)*fmp(i)/ppin(i)
   10   continue
        do 20 i=1,nwt
          plongt=rots*rots/4-(fmp(1)*fmp(1)+fmt(i)*fmt(i))/2+
     *    (fmp(1)*fmp(1)-fmt(i)*fmt(i))**2/(4*rots*rots)
          esen=sqrt(plongt+fmt(i)*fmt(i))
          pmta(i)=esen+sqrt(plongt)
          ppta(i)=fmt(i)*fmt(i)/pmta(i)
   20   continue
      elseif(lab) then
        do 11 i=1,nwp
          ppin(i)=elab+sqrt(elab*elab-fmp(i)*fmp(i))
          pmin(i)=fmp(i)*fmp(i)/ppin(i)
   11   continue
        do 21 i=1,nwt
          ppta(i)=fmt(i)
          pmta(i)=fmt(i)
   21   continue
      endif
c---- helgeC gives fermi-motion to all nucleons in the nuclei -----------
      if (ifermi.eq.1) call helgeC
c-----------------------------------------------------------------------
      do 30 i=1,ny
      ppi=ppin(ni(i))
      pmi=pmin(ni(i))
      ppt=ppta(nt(i))
      pmt=pmta(nt(i))
      ppstox=pxsp(ni(i))
      ppstoy=pysp(ni(i))
      ptstox=pxst(nt(i))
      ptstoy=pyst(nt(i))
c-------transformation to nucleon-nucleon cms---------------------------
      bosfac=sqrt((ppi+ppt)/(pmi+pmt))
      ppicms=ppi/bosfac
      pmicms=pmi*bosfac
      pptcms=ppt/bosfac
      pmtcms=pmt*bosfac
c------- if the strings are backing in the cms => no effective ---------
c------- collision -----------------------------------------------------
      if(pmicms.gt.ppicms.or.pptcms.gt.pmtcms) goto 30
      nyrak=0
   40 nyrak=nyrak+1
c------- if we get here to often the collision better not take place----
      if (nyrak.gt.50) goto 31

      pmin(ni(i))=pmicms*(pmtcms/pmicms)**rluC(0)
      ppta(nt(i))=pptcms*(ppicms/pptcms)**rluC(0)
      ppin(ni(i))=ppicms+pptcms-ppta(nt(i))
      pmta(nt(i))=pmicms+pmtcms-pmin(ni(i))
      pxsp(ni(i))=ppstox
      pysp(ni(i))=ppstoy
      pxst(nt(i))=ptstox
      pyst(nt(i))=ptstoy
c------- halvdaC gives the strings pt ----------------------------------a
      call halvdaC(i,kfel)
      if (kfel.eq.1) goto 40

c-------calculation of excited masses-----------------------------------
      ptsq1=pxsp(ni(i))*pxsp(ni(i))+pysp(ni(i))*pysp(ni(i))
      ptsq2=pxst(nt(i))*pxst(nt(i))+pyst(nt(i))*pyst(nt(i))

      win(ni(i))=sqrt(ppin(ni(i))*pmin(ni(i))-ptsq1)
      wta(nt(i))=sqrt(ppta(nt(i))*pmta(nt(i))-ptsq2)


c-------no projectilemasses < mp and no targetmasses < mt --------------
      if(win(ni(i)).lt.widmin.or.wta(nt(i)).lt.wtdmin) goto 40
c-------no double diffractive collisions--------------------------------
      if(win(ni(i)).lt.winmin.and.wta(nt(i)).lt.wtamin) goto 40
c-------new disrtr.of p+ & p- if either win or wta < winmin or wtamin---
c-------keep p+ for the projectile & p- for the target -----------------
c-------some kinematical tests------------------------------------------
      if(win(ni(i)).lt.winmin) then
        pminny=(fmp(ni(i))*fmp(ni(i))+ptsq1)/ppin(ni(i))
        ppinny=ppin(ni(i))
        ploc=(ppinny-pminny)/2
        eloc=(ppinny+pminny)/2
        if (ptsq1.ge.ploc**2) goto 40
        ploc=abs(ploc)/ploc*sqrt(ploc**2-ptsq1)
        eloc=sqrt(eloc**2-ptsq1)
        ppinny=eloc+ploc
        pminny=eloc-ploc
        delpm=pmin(ni(i))-pminny
        delpp=ppin(ni(i))-ppinny
        pmin(ni(i))=pminny
        ppin(ni(i))=ppinny
        pmta(nt(i))=pmta(nt(i))+delpm
        ppta(nt(i))=ppta(nt(i))+delpp

        if (ptsq1.ge.ppin(ni(i))*pmin(ni(i))) goto 40
        win(ni(i))=sqrt(ppin(ni(i))*pmin(ni(i))-ptsq1)
        wta(nt(i))=sqrt(ppta(nt(i))*pmta(nt(i))-ptsq2)
        if(wta(nt(i)).lt.wtamin) goto 40
      endif
      if(wta(nt(i)).lt.wtamin) then
        pptany=(fmt(nt(i))*fmt(nt(i))+ptsq2)/pmta(nt(i))
        pmtany=pmta(nt(i))
        ploc=(pptany-pmtany)/2
        eloc=(pptany+pmtany)/2
        if (ptsq2.ge.ploc**2) goto 40
        ploc=abs(ploc)/ploc*sqrt(ploc**2-ptsq2)
        eloc=sqrt(eloc**2-ptsq2)
        pptany=eloc+ploc
        pmtany=eloc-ploc
        delpm=pmta(nt(i))-pmtany
        delpp=ppta(nt(i))-pptany
        pmta(nt(i))=pmtany
        ppta(nt(i))=pptany
        pmin(ni(i))=pmin(ni(i))+delpm
        ppin(ni(i))=ppin(ni(i))+delpp

        if (ptsq2.ge.ppta(nt(i))*pmta(nt(i))) goto 40
        win(ni(i))=sqrt(ppin(ni(i))*pmin(ni(i))-ptsq1)
        wta(nt(i))=sqrt(ppta(nt(i))*pmta(nt(i))-ptsq2)
        if(win(ni(i)).lt.winmin) goto 40
      endif
c-------no backing strings----------------------------------------------
      if(pmin(ni(i)).gt.ppin(ni(i)).or.ppta(nt(i)).gt.pmta(nt(i)))
     &goto 40
c-------transformation back to lab--------------------------------------
      ppin(ni(i))=ppin(ni(i))*bosfac
      pmin(ni(i))=pmin(ni(i))/bosfac
      ppta(nt(i))=ppta(nt(i))*bosfac
      pmta(nt(i))=pmta(nt(i))/bosfac

   30 continue

      do 331 i=1,nwp
        if (win(i).gt.winmin) then
          dptsq=dble(pxsp(i))*dble(pxsp(i))+
     *      dble(pysp(i))*dble(pysp(i))
          dmassq=dble(ppin(i))*dble(pmin(i))
          win(i)=sngl(dsqrt(dmassq-dptsq))
        endif
  331 continue
      do 332 i=1,nwt
        if (wta(i).gt.wtamin) then
          dptsq=dble(pxst(i))*dble(pxst(i))+
     *      dble(pyst(i))*dble(pyst(i))
          dmassq=dble(ppta(i))*dble(pmta(i))
          wta(i)=sngl(dsqrt(dmassq-dptsq))
        endif
  332 continue

      return
      end


      subroutine halvdaC(i,kfel)
c-----------------------------------------------------------------------
c     giving the excited nucleons pt
c-----------------------------------------------------------------------
      common/frstptC/pxsp(40),pysp(40),pxst(200),pyst(200)
      common/frnynuC/bipa,nwp,nwt,ny,nucp(500),nuct(500)
      common/frstpmC/ppin(40),pmin(40),ppta(200),pmta(200)

      kfel=0

      delpt=sqrt(-0.08*log(1.-rluC(0)))
      delfi=2*3.14159*rluC(0)
      delpx=delpt*cos(delfi)
      delpy=delpt*sin(delfi)

      pxsp(nucp(i))=pxsp(nucp(i))+delpx
      pysp(nucp(i))=pysp(nucp(i))+delpy
      pxst(nuct(i))=pxst(nuct(i))-delpx
      pyst(nuct(i))=pyst(nuct(i))-delpy

      ptsq1=pxsp(nucp(i))**2+pysp(nucp(i))**2
      ptsq2=pxst(nuct(i))**2+pyst(nuct(i))**2

      if (ptsq1.ge.ppin(nucp(i))*pmin(nucp(i)).or.
     *ptsq2.ge.ppta(nuct(i))*pmta(nuct(i))) kfel=1

      return
      end


      subroutine hildinC
c-----------------------------------------------------------------------
c the subroutine gives the identity to the involved particles
c-----------------------------------------------------------------------
      common/indataC/elab,rots,nap,nzp,r0p,nat,nzt,r0t,iflspv,bmin,
     *                bmax,neve,iproty,ifermi,iflout
      common/frsttyC/idp(40),idt(250)
      common/frnynuC/bipa,ninkol,ntakol,nkoll,ni(500),nt(500)

c$$$$$$$$$$$$$$$$
c---- projectiletype p,pbar,pi+,pi- ------------------------------------
c     if(abs(iproty).eq.2) then
c       idp(1)=iproty*18/2
c       goto 32
c     elseif(abs(iproty).eq.1) then
c       idp(1)=iproty*17
c       goto 32
c     elseif(iproty.eq.0) then
c       if (nzp.eq.-1.) then
c         idp(1)=-41
c         goto 32
c       endif
c     endif
c-----------------------------------------------------------------------
c$$$$$$$$$$$$$$$$
c$$$$$$$$$$$$$$$$
      if(iproty .ne. 0) then
          idp(1)=iproty
      else
c        p_bar,p, nucleus
c$$$$$$$$$$$$$$$$
c$$$$$$$$$$$$$$$$
          if (nzp.eq.-1.) then
              idp(1)=-41
          else
c$$$$$$$$$$$$$$$$
              iz=nzp
              ia=nap
              do 30 i=1,nap
                 slump=rluC(0)
                 qpos=float(iz)/ia
                 if (slump.lt.qpos) then
                    idp(i)=41
                    ia=ia-1
                    iz=iz-1
                 else
                    idp(i)=42
                    ia=ia-1
                 endif
   30         continue
c$$$$$$$$$$$$$$$$
          endif
c$$$$$$$$$$$$$$$$
c$$$$$$$$$$$$$$$$
      endif
c$$$$$$$$$$$$$$$$

   32 iz=nzt
      ia=nat
      do 40 i=1,nat
         slump=rluC(0)
         qpos=float(iz)/ia
         if (slump.lt.qpos) then
            idt(i)=41
            ia=ia-1
            iz=iz-1
         else
            idt(i)=42
            ia=ia-1
         endif
   40 continue

      return
      end

      subroutine angantC
c-----------------------------------------------------------------------
c the subroutine records the two particles involved in each binary
c collision
c-----------------------------------------------------------------------
c-------monte carlo generation of wounded nucleons----------------------
      common/indataC/elab,rots,nap,nzp,r0p,nat,nzt,r0t,iflspv,bmin,
     *                bmax,neve,iproty,ifermi,iflout
      common/frnynuC/bipa,nwp,nwt,ny,nucp(500),nuct(500)
      dimension sumco(3)
      dimension rgt(5),grt(5),rgtcub(5),coordt(250,3)
      dimension rgp(5),grp(5),rgpcub(5),coordp(40,3)
      dimension markt(250)
c-------the dimensions should be coordt(>nat,3),coordp(>nap,3),---------
c-------markt(>nat),nucp&nuct(<nat*nap)---------------------------------
      data el91,el92,el93,el94,el95/2.19722,4.59512,
     *6.90675,9.21024,11.51292/
      data c0,pi,rmin/.545,3.14159,1.12838/
c-------some initialization---------------------------------------------
c$$$$$$$$$$$$$$$
c     rint=1.00925
c-------pion,kaon-------------------------------------------------------
c     if (abs(iproty).eq.1) rint=.81759
c     if (abs(iproty).eq.2) rint=.75694
c$$$$$$$$$$$$$$$$$
c
c$$$$$$$$$$$$$$
      if(iproty .eq. 0) then
         rint=1.00925
      elseif(abs(iproty) .eq. 17 .or. abs(iproty) .eq. 23) then
c         pion
         rint=.81795
      elseif(abs(iproty) .eq. 18 .or. abs(iproty) .eq. 19) then
c         kaon
         rint=.75694
      else
         rint=.76
      endif
c$$$$$$$$$$$$$$$$$$$$$$
      sumtra=0
      pimp=0.
      rmsq=rmin*rmin
      risq=rint*rint
      ravt=r0t*nat**(1./3.)
      rgt(1)=ravt+c0*el91
      rgt(2)=ravt+c0*el92
      rgt(3)=ravt+c0*el93
      rgt(4)=ravt+c0*el94
      rgt(5)=ravt+c0*el95
      do 10 mcub=1,5
        rgtcub(mcub)=rgt(mcub)*rgt(mcub)*rgt(mcub)
   10 continue
      grt(1)=rgtcub(1)/3.
      grt(2)=rgtcub(2)/30.
      grt(3)=rgtcub(3)/300.
      grt(4)=rgtcub(4)/3000.
      grt(5)=rgtcub(5)/30000.
      grt(5)=grt(5)-grt(4)/10.
      grt(4)=grt(4)-grt(3)/10.
      grt(3)=grt(3)-grt(2)/10.
      grt(2)=grt(2)-grt(1)/10.
      grt(2)=grt(2)+grt(1)
      grt(3)=grt(3)+grt(2)
      grt(4)=grt(4)+grt(3)
      grt(5)=grt(5)+grt(4)
      grt(1)=grt(1)/grt(5)
      grt(2)=grt(2)/grt(5)
      grt(3)=grt(3)/grt(5)
      grt(4)=grt(4)/grt(5)
      grt(5)=1.
      ravp=r0p*nap**(1./3.)
      rgp(1)=ravp+c0*el91
      rgp(2)=ravp+c0*el92
      rgp(3)=ravp+c0*el93
      rgp(4)=ravp+c0*el94
      rgp(5)=ravp+c0*el95
      do 1010 mcub=1,5
        rgpcub(mcub)=rgp(mcub)*rgp(mcub)*rgp(mcub)
 1010 continue
      grp(1)=rgpcub(1)/3.
      grp(2)=rgpcub(2)/30.
      grp(3)=rgpcub(3)/300.
      grp(4)=rgpcub(4)/3000.
      grp(5)=rgpcub(5)/30000.
      grp(5)=grp(5)-grp(4)/10.
      grp(4)=grp(4)-grp(3)/10.
      grp(3)=grp(3)-grp(2)/10.
      grp(2)=grp(2)-grp(1)/10.
      grp(2)=grp(2)+grp(1)
      grp(3)=grp(3)+grp(2)
      grp(4)=grp(4)+grp(3)
      grp(5)=grp(5)+grp(4)
      grp(1)=grp(1)/grp(5)
      grp(2)=grp(2)/grp(5)
      grp(3)=grp(3)/grp(5)
      grp(4)=grp(4)/grp(5)
      grp(5)=1.
c-------start generate--------------------------------------------------
c-----------------------------------------------------------------------
c-------create the nat nucleons of the target---------------------------
 3333    if (nat.le.8) goto 950
c-------jump to special for alpha---------------------------------------
         do 101 j1=1,nat
  900       slu=rluC(0)
            if (slu.ge.grt(1)) goto 901
            slu=rluC(0)
            r=rgt(1)*slu**(1./3.)
            slu=rluC(0)
            rjam=(r-ravt)/c0
            rjam=1.+exp(rjam)
            rjam=1./rjam
            if (slu.ge.rjam) goto 900
            goto 905
  901       if (slu.ge.grt(2)) goto 902
            slu=rluC(0)
            r=rgtcub(2)-rgtcub(1)
            r=r*slu+rgtcub(1)
            r=r**(1./3.)
            slu=rluC(0)
            rjam=(r-ravt)/c0
            rjam=1.+exp(rjam)
            rjam=10./rjam
            if (slu.ge.rjam) goto 900
            goto 905
  902       if (slu.ge.grt(3)) goto 903
            slu=rluC(0)
            r=rgtcub(3)-rgtcub(2)
            r=r*slu+rgtcub(2)
            r=r**(1./3.)
            slu=rluC(0)
            rjam=(r-ravt)/c0
            rjam=1.+exp(rjam)
            rjam=100./rjam
            if (slu.ge.rjam) goto 900
            goto 905
  903       if (slu.ge.grt(4)) goto 904
            slu=rluC(0)
            r=rgtcub(4)-rgtcub(3)
            r=r*slu+rgtcub(3)
            r=r**(1./3.)
            slu=rluC(0)
            rjam=(r-ravt)/c0
            rjam=1.+exp(rjam)
            rjam=1000./rjam
            if (slu.ge.rjam) goto 900
            goto 905
  904       slu=rluC(0)
            r=rgtcub(5)-rgtcub(4)
            r=r*slu+rgtcub(4)
            r=r**(1./3.)
            slu=rluC(0)
            rjam=(r-ravt)/c0
            rjam=1.+exp(rjam)
            rjam=10000./rjam
            if (slu.ge.rjam) goto 900
  905       continue
            do 777 nerr=1,10
               slu=rluC(0)
               phi=slu*2.*pi
               slu=rluC(0)
               ctheta=1.-2*slu
               stheta=sqrt(1.-ctheta*ctheta)
               coordt(j1,1)=r*ctheta
               coordt(j1,2)=r*stheta*cos(phi)
               coordt(j1,3)=r*stheta*sin(phi)
c-------a nucleon is placed - check if it has a neighbour---------------
c-------too close-------------------------------------------------------
               if (j1.eq.1) goto 910
               do 911 j2=1,j1-1
                  avsx=coordt(j1,1)-coordt(j2,1)
                  avsy=coordt(j1,2)-coordt(j2,2)
                  avsz=coordt(j1,3)-coordt(j2,3)
                  avsq=avsx*avsx+avsy*avsy+avsz*avsz
                  if (avsq.le.rmsq) goto 777
  911          continue
               goto 910
  777       continue
            goto 900
  910       continue
c--------the nucleon (finally) passed the test--------------------------
c--------go on with next------------------------------------------------
  101    continue
c--------a target nucleus is created------------------------------------
c--------now the turn comes to the projectile---------------------------
c--------but first jump over the alpha stuff----------------------------
         goto 110
c--------special for alpha----------------------------------------------
  950    do 951 j1=1,nat
  952       do 953 j3=1,3
               sl1=rluC(0)
               sl1=-2.*log(sl1)
               if (nat.eq.1) then
                 sl1=0.
               else
                 sl1=sqrt(sl1*nat/(nat-1))*r0t
               endif
               sl2=rluC(0)
               sl2=cos(2.*pi*sl2)
               coordt(j1,j3)=sl1*sl2
  953       continue
            if (j1.eq.1) goto 951
            do 954 j2=1,j1-1
               avsx=coordt(j1,1)-coordt(j2,1)
               avsy=coordt(j1,2)-coordt(j2,2)
               avsz=coordt(j1,3)-coordt(j2,3)
               avsq=avsx*avsx+avsy*avsy+avsz*avsz
               if (avsq.le.rmsq) goto 952
  954       continue
  951    continue
c--------end of special for alpha---------------------------------------
  110    continue
c--------center the nucleons inside the nucleus-------------------------
         do 115 j2=1,3
            sumco(j2)=0.
  115    continue
         do 111 j1=1,nat
            do 112 j2=1,3
               sumco(j2)=sumco(j2)+coordt(j1,j2)
  112       continue
  111    continue
         do 116 j1=1,nat
            do 117 j2=1,3
               coordt(j1,j2)=coordt(j1,j2)-sumco(j2)/nat
  117       continue
  116    continue
c--------sort the nucleons on increasing z-coordinate-------------------
         do 778 j1=2,nat
            zco=coordt(j1,3)
            do 779 j2=1,j1-1
               if (zco.lt.coordt(j2,3)) then
                  xco=coordt(j1,1)
                  yco=coordt(j1,2)
                  do 780 j3=1,j1-j2
                     do 781 i=1,3
                        coordt(j1+1-j3,i)=coordt(j1-j3,i)
  781                continue
  780             continue
                  coordt(j2,1)=xco
                  coordt(j2,2)=yco
                  coordt(j2,3)=zco
                  goto 778
              endif
  779      continue
  778    continue
c----------end sorting--------------------------------------------------
         xmaxt=-999.
         xmint=999.
         ymaxt=-999.
         ymint=999.
         do 103 j1=1,nat
            if (coordt(j1,1).ge.xmaxt)
     *      xmaxt=coordt(j1,1)
            if (coordt(j1,1).le.xmint)
     *      xmint=coordt(j1,1)
            if (coordt(j1,2).ge.ymaxt)
     *      ymaxt=coordt(j1,2)
            if (coordt(j1,2).le.ymint)
     *      ymint=coordt(j1,2)
  103    continue
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
         if (nap.le.8) goto 1950
c-------jump to special for alpha---------------------------------------
         do 1101 l1=1,nap
 1900       slu=rluC(0)
            if (slu.ge.grp(1)) goto 1901
            slu=rluC(0)
            r=rgp(1)*slu**(1./3.)
            slu=rluC(0)
            rjam=(r-ravp)/c0
            rjam=1.+exp(rjam)
            rjam=1./rjam
            if (slu.ge.rjam) goto 1900
            goto 1905
 1901       if (slu.ge.grp(2)) goto 1902
            slu=rluC(0)
            r=rgpcub(2)-rgpcub(1)
            r=r*slu+rgpcub(1)
            r=r**(1./3.)
            slu=rluC(0)
            rjam=(r-ravp)/c0
            rjam=1.+exp(rjam)
            rjam=10./rjam
            if (slu.ge.rjam) goto 1900
            goto 1905
 1902       if (slu.ge.grp(3)) goto 1903
            slu=rluC(0)
            r=rgpcub(3)-rgpcub(2)
            r=r*slu+rgpcub(2)
            r=r**(1./3.)
            slu=rluC(0)
            rjam=(r-ravp)/c0
            rjam=1.+exp(rjam)
            rjam=100./rjam
            if (slu.ge.rjam) goto 1900
            goto 1905
 1903       if (slu.ge.grp(4)) goto 1904
            slu=rluC(0)
            r=rgpcub(4)-rgpcub(3)
            r=r*slu+rgpcub(3)
            r=r**(1./3.)
            slu=rluC(0)
            rjam=(r-ravp)/c0
            rjam=1.+exp(rjam)
            rjam=1000./rjam
            if (slu.ge.rjam) goto 1900
            goto 1905
 1904       slu=rluC(0)
            r=rgpcub(5)-rgpcub(4)
            r=r*slu+rgpcub(4)
            r=r**(1./3.)
            slu=rluC(0)
            rjam=(r-ravp)/c0
            rjam=1.+exp(rjam)
            rjam=10000./rjam
            if (slu.ge.rjam) goto 1900
 1905       continue
            do 1777 nerr=1,10
               slu=rluC(0)
               phi=slu*2.*pi
               slu=rluC(0)
               ctheta=1.-2*slu
               stheta=sqrt(1.-ctheta*ctheta)
               coordp(l1,1)=r*ctheta
               coordp(l1,2)=r*stheta*cos(phi)
               coordp(l1,3)=r*stheta*sin(phi)
c-------a nucleon is placed - check if it has a neighbour---------------
c-------too close-------------------------------------------------------
               if (l1.eq.1) goto 1910
               do 1911 l2=1,l1-1
                  avsx=coordp(l1,1)-coordp(l2,1)
                  avsy=coordp(l1,2)-coordp(l2,2)
                  avsz=coordp(l1,3)-coordp(l2,3)
                  avsq=avsx*avsx+avsy*avsy+avsz*avsz
                  if (avsq.le.rmsq) goto 1777
 1911          continue
               goto 1910
 1777       continue
            goto 1900
 1910       continue
c--------the nucleon (finally) passed the test--------------------------
c--------go on with next------------------------------------------------
 1101    continue
c--------a projectile nucleus is created--------------------------------
c--------now the colliding starts---------------------------------------
c--------but first jump over the alpha stuff----------------------------
         goto 1110
c--------special for alpha----------------------------------------------
 1950    do 1951 l1=1,nap
 1952       do 1953 l3=1,3
               sl1=rluC(0)
               sl1=-2.*log(sl1)
               if (nap.eq.1) then
                 sl1=0.
               else
                 sl1=sqrt(sl1*nap/(nap-1))*r0p
               endif
               sl2=rluC(0)
               sl2=cos(2.*pi*sl2)
               coordp(l1,l3)=sl1*sl2
 1953       continue
            if (l1.eq.1) goto 1951
            do 1954 l2=1,l1-1
               avsx=coordp(l1,1)-coordp(l2,1)
               avsy=coordp(l1,2)-coordp(l2,2)
               avsz=coordp(l1,3)-coordp(l2,3)
               avsq=avsx*avsx+avsy*avsy+avsz*avsz
               if (avsq.le.rmsq) goto 1952
 1954       continue
 1951    continue
c--------end of special for alpha---------------------------------------
 1110    continue
c--------center the nucleons inside the nucleus-------------------------
         do 1115 l2=1,3
            sumco(l2)=0.
 1115    continue
         do 1111 l1=1,nap
            do 1112 l2=1,3
               sumco(l2)=sumco(l2)+coordp(l1,l2)
 1112       continue
 1111    continue
         do 1116 l1=1,nap
            do 1117 l2=1,3
               coordp(l1,l2)=coordp(l1,l2)-sumco(l2)/nap
 1117       continue
 1116    continue
c--------sort the nucleons on increasing z-coordinate-------------------
         do 1778 l1=2,nap
            zco=coordp(l1,3)
            do 1779 l2=1,l1-1
               if (zco.lt.coordp(l2,3)) then
                  xco=coordp(l1,1)
                  yco=coordp(l1,2)
                  do 1780 l3=1,l1-l2
                     do 1781 i=1,3
                        coordp(l1+1-l3,i)=coordp(l1-l3,i)
 1781                continue
 1780             continue
                  coordp(l2,1)=xco
                  coordp(l2,2)=yco
                  coordp(l2,3)=zco
                  goto 1778
              endif
 1779      continue
 1778    continue
c----------end sorting--------------------------------------------------
         xmaxp=-999.
         xminp=999.
         ymaxp=-999.
         yminp=999.
         do 1103 l1=1,nap
            if (coordp(l1,1).ge.xmaxp)
     *      xmaxp=coordp(l1,1)
            if (coordp(l1,1).le.xminp)
     *      xminp=coordp(l1,1)
            if (coordp(l1,2).ge.ymaxp)
     *      ymaxp=coordp(l1,2)
            if (coordp(l1,2).le.yminp)
     *      yminp=coordp(l1,2)
 1103    continue
         xmax=xmaxt-xminp+rint
         xmin=xmint-xmaxp-rint
         ymax=ymaxt-yminp+rint
         ymin=ymint-ymaxp-rint
         ba=(xmax-xmin)*(ymax-ymin)
         ntra=0
         if (iflspv.eq.2.or.iflspv.eq.3) then
            slu=rluC(0)
            genus=sqrt(slu*(bmax*bmax-bmin*bmin)+bmin*bmin)
            slu=rluC(0)
            phgen=slu*2.*pi
            xgen=genus*cos(phgen)
            ygen=genus*sin(phgen)
         else
            slu=rluC(0)
            xgen=(xmax-xmin)*slu+xmin
            slu=rluC(0)
            ygen=(ymax-ymin)*slu+ymin
         endif
         bipa=sqrt(xgen*xgen+ygen*ygen)
         pimp=pimp+1.
         do 2010 j1=1,nat
            markt(j1)=0
 2010    continue
         ny=0
         nwp=0
         do 2020 l1=1,nap
            markp=0
            xlgen=coordp(l1,1)+xgen
            ylgen=coordp(l1,2)+ygen
            if (xlgen.gt.xmaxt+rint) goto 2020
            if (xlgen.lt.xmint-rint) goto 2020
            if (ylgen.gt.ymaxt+rint) goto 2020
            if (ylgen.lt.ymint-rint) goto 2020
            do 105 j1=1,nat
               ravsq=(coordt(j1,1)-xlgen)*
     *               (coordt(j1,1)-xlgen)+
     *               (coordt(j1,2)-ylgen)*
     *               (coordt(j1,2)-ylgen)
               if (risq.lt.ravsq) goto 105
               ny=ny+1
               markt(j1)=1
               markp=1
               nucp(ny)=l1
               nuct(ny)=j1
  105       continue
            nwp=nwp+markp
 2020    continue
         if ((iflspv.eq.1.or.iflspv.eq.3).and.nwp.lt.nap) goto 3333
         if (nwp.eq.0) goto 3333
         nwt=0
         do 2030 j1=1,nat
            nwt=nwt+markt(j1)
 2030    continue
         do 2700 i=1,nap
            ii=999
            do 2701 j=1,ny
               if (nucp(j).ge.i.and.nucp(j).lt.ii) ii=nucp(j)
 2701       continue
            if (ii.eq.i) goto 2700
            if (ii.eq.999) goto 2703
            do 2702 j=1,ny
               if (nucp(j).eq.ii) nucp(j)=i
 2702       continue
 2700    continue
 2703    do 2710 i=1,nat
            ii=999
            do 2711 j=1,ny
               if (nuct(j).ge.i.and.nuct(j).lt.ii) ii=nuct(j)
 2711       continue
            if (ii.eq.i) goto 2710
            if (ii.eq.999) goto 2713
            do 2712 j=1,ny
               if (nuct(j).eq.ii) nuct(j)=i
 2712       continue
 2710    continue
 2713    continue
      return
      end

      subroutine helgeC
c-----------------------------------------------------------------------
c---- helgeC gives fermi-motion to the nucleons in the nuclei -----------
c-----------------------------------------------------------------------
      common/indataC/elab,rots,nap,nzp,r0p,nat,nzt,r0t,iflspv,bmin,
     *                bmax,neve,iproty,ifermi,iflout
      common/frstptC/pxsp(40),pysp(40),pxst(200),pyst(200)
      common/frstpmC/ppin(40),pmin(40),ppta(200),pmta(200)
      common/frnynuC/bipa,nwp,nwt,ny,ni(500),nt(500)
      dimension fvect(3,250)
      data pi /3.14159/

      if (nap.le.1) goto 100
   19 do 10 i=1,3
        sum=0.
        do 11 j=1,nap
          sl1=rluC(0)
          sl2=rluC(0)
          sl1=sqrt(-2.*log(sl1)*nap/(nap-1))
          sl2=cos(2.*pi*sl2)
          fvect(i,j)=sl1*sl2*.1
          sum=sum+fvect(i,j)
   11   continue
        do 12 j=1,nwp
          fvect(i,j)=fvect(i,j)-sum/nap
   12   continue
   10 continue
      do 18 j=1,nwp
        sum=0.
        do 17 i=1,3
          sum=sum+fvect(i,j)**2
   17   continue
        if (sum.ge..3) goto 16
   18 continue
   16 if (sum.ge..3) goto 19
      do 13 j=1,nwp
        pxsp(j)=fvect(1,j)
        pysp(j)=fvect(2,j)
        bosfac=sqrt(ppin(j)/pmin(j))
        pplus=ppin(j)/bosfac
        pminus=pmin(j)*bosfac
        pplus=sqrt(pplus**2+fvect(3,j)**2)+fvect(3,j)
        pminus=sqrt(pminus**2+fvect(3,j)**2)-fvect(3,j)
        ppin(j)=pplus*bosfac
        pmin(j)=pminus/bosfac
   13 continue
  100 if (nat.le.1) goto 200
  119 do 110 i=1,3
        sum=0.
        do 111 j=1,nat
          sl1=rluC(0)
          sl2=rluC(0)
          sl1=sqrt(-2.*log(sl1)*nat/(nat-1))
          sl2=cos(2.*pi*sl2)
          fvect(i,j)=sl1*sl2*.1
          sum=sum+fvect(i,j)
  111   continue
        do 112 j=1,nwt
          fvect(i,j)=fvect(i,j)-sum/nat
  112   continue
  110 continue
      do 118 j=1,nwt
        sum=0.
        do 117 i=1,3
          sum=sum+fvect(i,j)**2
  117   continue
        if (sum.ge..3) goto 116
  118 continue
  116 if (sum.ge..3) goto 119
      do 113 j=1,nwt
        pxst(j)=fvect(1,j)
        pyst(j)=fvect(2,j)
        bosfac=sqrt(ppta(j)/pmta(j))
        pplus=ppta(j)/bosfac
        pminus=pmta(j)*bosfac
        pplus=sqrt(pplus**2+fvect(3,j)**2)+fvect(3,j)
        pminus=sqrt(pminus**2+fvect(3,j)**2)-fvect(3,j)
        ppta(j)=pplus*bosfac
        pmta(j)=pminus/bosfac
  113 continue
  200 return
      end


      subroutine ellidaC
c-----------------------------------------------------------------------
c-----ellidaC for eventlisting-------------------------------------------
c-----------------------------------------------------------------------
      common /indataC/ elab,rots,nap,nzp,r0p,nat,nzt,r0t,iflspv,bmin,
     *                bmax,neve,iproty,ifermi,iflout
      common /evevecC/ nevent,isppp,isppn,isptp,isptn,bimp,
     *                idi(2000),ipr(2000)
      common /lujetsC/ n,k(2000,2),p(2000,5)
      common /ludat1C/ mst(40),par(80)
      character chap*8

      if(nevent.eq.1) then
      write(mst(20),*)
      write(mst(20),*)
      write(mst(20),*)
      if(rots.eq.0.) write(mst(20),801)'energy',elab,' a gev'
      if(elab.eq.0.) write(mst(20),801)'roots ',rots,' a gev'
      write(mst(20),802)'projectile',' a=',nap,' z=',nzp
      write(mst(20),802)'target    ',' a=',nat,' z=',nzt
      write(mst(20),803)'this list contains',neve,' events'
      if(iflspv.eq.0) then
      write(mst(20),804)'minimum bias  '
      elseif(iflspv.eq.1) then
      write(mst(20),804)'spectator veto'
      elseif(iflspv.eq.2) then
      write(mst(20),804)'impact setting',bmin,' -',bmax,' fermi'
      else
      write(mst(20),804)'spect veto +impact setting',bmin,' -',
     *bmax,' fermi'
      endif
  801 format (' ',a,f6.0,a)
  802 format (' ',a,a,i4,a,i4)
  803 format (' ',a,i7,a)
  804 format (' ',a,f6.2,a,f6.2,a)
      endif

      write(mst(20),805) 'event',nevent
      write(mst(20),806) 'projectile spectators',
     *              '  protons',isppp,' neutrons',isppn
      write(mst(20),806) 'target     spectators',
     *              '  protons',isptp,' neutrons',isptn
      write(mst(20),807) 'impact parameter',bimp,' fermi'
      write(mst(20),808) 'the event contains',n,' lines'
  805 format (' ',a,i7)
  806 format (' ',a,a,i4,a,i4)
  807 format (' ',a,f6.2,a)
  808 format (' ',a,i5,a)

      if (n.eq.0) return
      if(iflout.eq.-1) then
        call lulistC(1)
        write(mst(20),*)
        write(mst(20),*)
      else
        write(mst(20),809) '    i  ori code  q d particle    px',
     *                '       py       pz     energy    mass'

        do 1 i=1,n
          if(iflout.eq.0) then
            iori=mod(k(i,1),10000)
          else
            iori=0
          endif
          ityp=k(i,2)
          icha=kluC(i,3)/3
          if(kluC(i,11).eq.0) then
            idec=1
          else
            idec=0
          endif
          call lunameC(ityp,chap)
          if (idec.eq.1) chap(8:8)='d'
          write(mst(20),810) i,iori,ityp,icha,idec,
     *                 chap,(p(i,j),j=1,5)
    1   continue
  809   format (' ',a,a)
  810   format (' ',3i5,i3,i2,1x,a8,5(1x,f8.3))
      endif
      return
      end


c-----routines same as luroboC except that the boost parameter beta------
c-----is changed to 1-beta and 1+beta respectively to make the boost----
c-----more acurate------------------------------------------------------
c***********************************************************************

      subroutine luroinC(the,phi,bex,bey,embez)
      common /lujetsC/ n,k(2000,2),p(2000,5)
      common /ludat1C/ mst(40),par(80)
      dimension rot(3,3),pv(3)
      double precision dp(4),dbex,dbey,dbez,dga,dbep,dgabep

      imax=n
      if(mst(2).gt.0) imax=mst(2)
      if(the**2+phi**2.lt.1e-20) goto 130
c...rotate (typically from z axis to direction theta,phi)
      rot(1,1)=cos(the)*cos(phi)
      rot(1,2)=-sin(phi)
      rot(1,3)=sin(the)*cos(phi)
      rot(2,1)=cos(the)*sin(phi)
      rot(2,2)=cos(phi)
      rot(2,3)=sin(the)*sin(phi)
      rot(3,1)=-sin(the)
      rot(3,2)=0.
      rot(3,3)=cos(the)
      do 120 i=max(1,mst(1)),imax
      if(mod(k(i,1)/10000,10).ge.6) goto 120
      do 100 j=1,3
  100 pv(j)=p(i,j)
      do 110 j=1,3
  110 p(i,j)=rot(j,1)*pv(1)+rot(j,2)*pv(2)+rot(j,3)*pv(3)
  120 continue

  130 if(bex**2+bey**2+(1.-embez)**2.lt.1e-20) return
c...lorentz boost (typically from rest to momentum/energy=beta)
      dbex=dble(bex)
      dbey=dble(bey)
      dbez=1d0-dble(embez)
      dga=1d0/dsqrt(1d0-dbex**2-dbey**2-dbez**2)
      do 150 i=max(1,mst(1)),imax
      if(mod(k(i,1)/10000,10).ge.6) goto 150
      do 140 j=1,4
  140 dp(j)=dble(p(i,j))
      dbep=dbex*dp(1)+dbey*dp(2)+dbez*dp(3)
      dgabep=dga*(dga*dbep/(1d0+dga)+dp(4))
      p(i,1)=sngl(dp(1)+dgabep*dbex)
      p(i,2)=sngl(dp(2)+dgabep*dbey)
      p(i,3)=sngl(dp(3)+dgabep*dbez)
      p(i,4)=sngl(dga*(dp(4)+dbep))
  150 continue

      return
      end

c*********************************************************************

      subroutine lurotaC(the,phi,bex,bey,epbez)
      common /lujetsC/ n,k(2000,2),p(2000,5)
      common /ludat1C/ mst(40),par(80)
      dimension rot(3,3),pv(3)
      double precision dp(4),dbex,dbey,dbez,dga,dbep,dgabep

      imax=n
      if(mst(2).gt.0) imax=mst(2)
      if(the**2+phi**2.lt.1e-20) goto 130
c...rotate (typically from z axis to direction theta,phi)
      rot(1,1)=cos(the)*cos(phi)
      rot(1,2)=-sin(phi)
      rot(1,3)=sin(the)*cos(phi)
      rot(2,1)=cos(the)*sin(phi)
      rot(2,2)=cos(phi)
      rot(2,3)=sin(the)*sin(phi)
      rot(3,1)=-sin(the)
      rot(3,2)=0.
      rot(3,3)=cos(the)
      do 120 i=max(1,mst(1)),imax
      if(mod(k(i,1)/10000,10).ge.6) goto 120
      do 100 j=1,3
  100 pv(j)=p(i,j)
      do 110 j=1,3
  110 p(i,j)=rot(j,1)*pv(1)+rot(j,2)*pv(2)+rot(j,3)*pv(3)
  120 continue

  130 if(bex**2+bey**2+(epbez-1.)**2.lt.1e-20) return
c...lorentz boost (typically from rest to momentum/energy=beta)
      dbex=dble(bex)
      dbey=dble(bey)
      dbez=dble(epbez)-1d0
      dga=1d0/dsqrt(1d0-dbex**2-dbey**2-dbez**2)
      do 150 i=max(1,mst(1)),imax
      if(mod(k(i,1)/10000,10).ge.6) goto 150
      do 140 j=1,4
  140 dp(j)=dble(p(i,j))
      dbep=dbex*dp(1)+dbey*dp(2)+dbez*dp(3)
      dgabep=dga*(dga*dbep/(1d0+dga)+dp(4))
      p(i,1)=sngl(dp(1)+dgabep*dbex)
      p(i,2)=sngl(dp(2)+dgabep*dbey)
      p(i,3)=sngl(dp(3)+dgabep*dbez)
      p(i,4)=sngl(dga*(dp(4)+dbep))
  150 continue

      return
      end
