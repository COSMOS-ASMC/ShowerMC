c***********************************************************************
c***********************************************************************
c                                                                      *
c        PART 3: Low energy h-h collisions                             *
c                                                                      *
c   List of subprograms in rough order of relevance with main purpose  *
c      (S = subroutine, F = function, B = block data, E = entry)       *
c                                                                      *
c s jamcross to administer the h-h collisions                          *
c s jamcbb   to handle non-strange BB collisions (S=0)                 *
c s jamcbb1  to handle N + N* or N* + N* collisions                    *
c s jamcbb2  to handle N(*) + D(*) collisions                          *
c s jamcbb3  to handle D(*) + D(*) collisions                          *
c s jamcbbs  to handle low energy Lambda-N/Sigma-N/Xi-N/LL collisions  *
c s jamcmbs0 to handle S=0 meson-baryon collisions                     *
c s jamcmbs1 to handle S=-1 meson-baryon collisions  Lambda-pi...      *
c s jamcmbs2 to handle S=-2 meson-baryon collisions  xi-N,...          *
c s jamckaon to treat Kaon-nonstrange baryon collisions                *
c s jamcpipi to treat low energy pi-pi collisions                      *
c s jamcabb  to handle antibaryon-baryon collisions                    *
c                                                                      *
c s jambmas  to find outgoing types in S=0 BB inel. collisions         *
c s jamrmas1 to generate masses accroding to the B-W for S=0 BB coll.  *
c s jamrmas2 to generate masses accroding to the B-W distribution      *
c s jamexpa  to find possible hadronic excitation states               *
c s jamdmas  to Generate an delta(1232)mass with BW                    *
c                                                                      *
c s jamxadq  to calculate cross section by additive quark model        *
c s jamxtot  to give total/elastic/diffractive cross sections          *
c s jamxnn   to calculate pp/pn total and elastic cross sections       *
c s jamxkp   to calculate K-bar N cross sections                       *
c s jamxpin  to give pi-N total and elastic cross sections             *
c s jamxbbar to give antibaryon-baryon cross sections                  *
c s jamxbw1  to calculate Berit-Wigner cross section for M-B           *
c s jamxbw2  to calculate Berit-Wigner cross section for M-M           *
c f jamxdelt to calculate pion plus nucleon to delta cross section     *
c                                                                      *
c s jambres1 to calculate individual resonance production probability  *
c s jambres2 to calculate non-strange baryonic resonance prob.         *
c s jamdetb1 to calculate integral for detailed balance for RN->NN     *
c f jambwf1  to provide B-W function for jamdetb1                      *
c s jamdetb2 to calculate integral for detailed balance for RR->NN     *
c f jambwf2  to provide B-W function for jamdetb2                      *
c f jambwf3  to provide B-W function for jambwf2                       *
c s jamdetb3 to calculate factor for RR->NN with constant width        *
c f jambwtbl to give integrated value of WB for non-strange resonances *
c f jamcpair to give order number for collision pair                   *
c                                                                      *
c f jamrgg92 to give total xsection fit by Regge theory based formula  *
c f jamrgg96 to give high energy cross sections by regge theory        *
c f jamchc96 to give fit formula for cross sections at high-energies   *
c f jamchc88 to give cross sections (particle data group  1988)        *
c s jamsighh to give h-h total and elastic x-section at low-energies   *
c s jamxnnin to give parametrization for nn-> x cross sections         *
c f jamsigkn to give kaon-nucleon one/two-pion production x-sections   *
c b jamsigda to contain data for low-energy hh cross sections          *
c s jamsigS1 to give S=-1 low energy Baryon-Baryon cross sections      *
c b jamss1da to contain data for S=-1 baryon - baryon cross sections   *
c s jamsigS2 to give S=-2 low energy Baryon-Baryon cross sections      *
c b jams2da  to contain S=-2 baryon - baryon cross sections data       *
c                                                                      *
c***********************************************************************

      subroutine jamcross(msel,icltyp,srt,pr,kf1,kf2,em1,em2,
     $                 sig,sigel,sigin,mchanel,mabsrb,ijet,icon)

c...Purpose: to administer h-h scattering.
c----------------------------------------------------------------------*
c...msel: =1: give total and elastic cross sections.
c...msel: =2: give inel. cross sections only.
c...msel: =3: Monte Calro evaluatin of inel. channel.
c...(Outputs)
c...mchanel: max. number of inelastic channel.
c...mabsrb : max. number of annihilation channel.
c...sig    : total cross section (mb)
c...sigel  : elastic cross section (mb)
c...sigin(): inelastic cross sections (mb)
c...Note
c...sigin() is orderd by 1...,mchanel, mchanel+1, ..., mchanel+mabsrb
c...icon   : Condition code.
c----------------------------------------------------------------------*
      
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

c...Commonblock for t-channel resonance productions.
      common/jamres1/kfo(2,20),noutpa

c...Local arrays.
      real*8 jamemjet
      parameter(mxchan=30)
      dimension sigin(mxchan)
      logical  anti

c...Check collision type.
      if(icltyp.le.0.or.icltyp.gt.6) then
        write(check(1),'(''icltyp='',i3)')icltyp
        write(check(2),'(''kf1 kf2='',i9,1x,i9)')kf1,kf2
        call jamerrm(30,2,'(jamcross:) invalid icltyp')
      endif

      if(msel.eq.3.and.(kf1.eq.92.or.kf2.eq.92)) then
        ijet=1
        return
      endif

c...Compressed code.
      kc1=jamcomp(kf1)
      kc2=jamcomp(kf2)

c...Check particle IDs.
      if((kc1.le.0.or.kc1.gt.mstu(6))
     $   .or. (kc2.le.0.or.kc2.gt.mstu(6))) then
        write(check(1),'(''kc1 kf1='',i3,1x,i9)')kc1,kf1
        write(check(2),'(''kc2 kf2='',i3,1x,i9)')kc2,kf2
        call jamerrm(30,2,'(jamcross:) invalid particle ID')
      endif

      ijet=0
      icon=0
      mchanel=0  ! # of t-channel
      mabsrb=0   ! # of s-channel
      mste(3)=0  ! final particle
      noutpa=0
      anti=.false.

c...Hadron-quark/parton-parton
      if(icltyp.ge.5) then
        if(kf1.eq.0.or.kf2.eq.0) then
          write(check(1),'(''icltyp='',i4)')icltyp
          write(check(2),'(''kf1 kf2'',i9,1x,i9)')kf1,kf2
          write(check(3),'(''kf1 kf2'',i3,1x,i3)')kc1,kc2
          call jamerrm(30,3,'(jamcross:)invalid icltyp kf1 kf2')
        endif
        call jamxadq(kf1,kf2,sig,sigel)
        mchanel=1
        return
      endif

      kfa=abs(kf1) 
      kf1a=mod(kfa/1000,10) 
      kf1b=mod(kfa/100,10) 
      kf1c=mod(kfa/10,10) 
      kfa=abs(kf2) 
      kf2a=mod(kfa/1000,10) 
      kf2b=mod(kfa/100,10) 
      kf2c=mod(kfa/10,10) 
      nheavy1=0
      nheavy2=0
      if(kf1a.ge.4) nheavy1=nheavy1+1
      if(kf1b.ge.4) nheavy1=nheavy1+1
      if(kf1c.ge.4) nheavy1=nheavy1+1
      if(kf2a.ge.4) nheavy2=nheavy2+1
      if(kf2b.ge.4) nheavy2=nheavy2+1
      if(kf2c.ge.4) nheavy2=nheavy2+1
      nheavy=nheavy1+nheavy2

      istr1=kchg(kc1,7)*isign(1,kf1)
      istr2=kchg(kc2,7)*isign(1,kf2)
      istr=istr1+istr2

c===================================
c...Baryon + baryon collision
c===================================
      if(icltyp.eq.1) then


c...If antibaryon-antibaryon collision,we should use the G-parity
c...transformation.
        if(kf1.lt.0.and.kf2.lt.0) then  ! antib-antib
          kf1=-kf1
          kf2=-kf2
          anti=.true.
        endif

        if(istr.eq.0.and.nheavy.eq.0) then

          call jamcbb(msel,srt,pr,kf1,kf2,kc1,kc2,em1,em2,
     $                  mchanel,sig,sigel,sigin,ic,mxchan,icon)
          if(icon.ne.0) goto 5000 ! skip collision.
          if(ic.eq.0) ijet=2      ! string formation.
          if(msel.eq.3.and.ijet.eq.2.and.srt.le.3.5d0) then
            io=mstc(38)
            write(io,*)'(jamcross:)low energy for string mchanl',mchanel
            write(io,*)'srt ijet icon ic',srt,ijet,icon,ic
            write(io,*)'pare(3) kf1 kf2 em1 em2',pare(3),kf1,kf2,em1,em2
            do i=1,mchanel
            write(io,*)sigin(i)
            end do
            icon=2
            goto 5000
          endif

c...Lambda N,etc.
        else if(nheavy.eq.0) then

c....Low energy Lambda-N/Sigma-N/Xi-N/LL...
          if((istr.eq.-1.and.srt.le.2.52d0)
     $                   .or.(istr.eq.-2.and.srt.le.2.73d0)) then
            call jamcbbs(msel,kf1,kf2,kc1,kc2,istr,srt,
     $                 sig,sigel,sigin,em1,em2,mchanel,mxchan,icon)
            if(icon.eq.0) goto 5000
          endif

          call jamxadq(kf1,kf2,sig,sigel)
          if(msel.eq.3) then
            if(srt.lt.parc(61)+0.3d0) then
              call jamrmas2(kf1,kf2,kc1,kc2,srt,em1,em2,icon)
            else
              ijet=1
            endif
          else
             mchanel=1
             sigin(1)=sig-sigel
          endif

        end if
        goto 5000

c===================================
c...Baryon + meson collisions part
c===================================
      else if(icltyp.eq.2) then

        ibar1=kchg(kc1,6)*isign(1,kf1)
        ibar2=kchg(kc2,6)*isign(1,kf2)

c...Check if this is a anti-B + M collision.
        if(ibar1.eq.-3.or.ibar2.eq.-3) then
         anti=.true.
         if(kchg(kc1,3).eq.1) kf1=-kf1
         if(kchg(kc2,3).eq.1) kf2=-kf2
        endif

        if(ibar1.ne.0) then
          kfb=kf1
          kfm=kf2
          embar=em1
          emmes=em2
          kcb=kc1
          kcm=kc2
          jswap=1
        else
          kfb=kf2
          kfm=kf1
          embar=em2
          emmes=em1
          kcb=kc2
          kcm=kc1
          jswap=0
        endif

        izb=kchg(kcb,1)*isign(1,kfb)/3
        izm=kchg(kcm,1)*isign(1,kfm)/3
        izt=izm+izb
        istr1=kchg(kcb,7)*isign(1,kfb)
        istr2=kchg(kcm,7)*isign(1,kfm)
        istr=istr1+istr2
        kfm0=kfm
        kfb0=kfb

c...Light collision systems.
        if(nheavy.eq.0) then
         ipath=0

c...piN/rhoN lambdaK+ etc.
          if(istr.eq.0) then
            call jamcmbs0(msel,srt,pr,istr,kfb,kfm,kcb,kcm,izb,izm
     $                ,emmes,embar,kf1,kf2,em1,em2
     $                ,sig,sigel,sigin,mabsrb,mchanel,mxchan
     $                ,ijet,jswap,icon)
           ipath=1

c....K-p,Lambda+pi etc
          else if(istr.eq.-1) then
            call jamcmbs1(msel,srt,pr,istr,kfb,kfm,kcb,kcm,izb,izm
     $               ,emmes,embar,kf1,kf2,em1,em2
     $               ,sig,sigel,sigin,mabsrb,mchanel,mxchan
     $               ,ijet,jswap,icon)
            ipath=2

c...Xi-pi etc.
          else if(istr.eq.-2) then
            call jamcmbs2(msel,srt,pr,istr,kfb,kfm,kcb,kcm,izb,izm
     $                   ,em1,em2,kf1,kf2,sig,sigel,sigin,mabsrb
     $                   ,mchanel,mxchan,ijet,icon)
            ipath=3

c....Kaon induced collisions.
          else if(istr.eq.1.and.istr1.eq.0) then
            call jamckaon(msel,kfm,kfb,srt,pr,ijet
     $                 ,emmes,embar,kf1,kf2,sig,sigel,sigin,mchanel
     $                 ,mxchan,ich,icon)
            if(icon.eq.999) goto 5000
            ipath=4

            if(msel.eq.3) then
              if(pare(3).lt.0.0d0) then
               if(ijet.ne.0) goto 5000
               if(jswap.eq.1) then
                 em1=embar
                 em2=emmes
                 itmp=kf1
                 kf1=kf2
                 kf2=itmp
               else
                 em1=emmes
                 em2=embar
               endif
              endif
            endif

          endif
   
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        if(msel.eq.3) then
        kcc1=jamcomp(kf1)
        kcc2=jamcomp(kf2)
        izf1=0 
        izf2=0 
        if(kcc1.ge.1) izf1=kchg(kcc1,1)*isign(1,kf1)/3
        if(kcc2.ge.1) izf2=kchg(kcc2,1)*isign(1,kf2)/3
        if(izf1+izf2.ne.izt) then
          print *,'charge not conserved ipath=',ipath,izt,izf1+izf2
          print *,'izm izb',izm,izb,' iz1 iz2',izf1,izf2
          print *,'kfm em1',kfm0,emmes,' kfb em2',kfb0,embar
          print *,'kf1 em1',kf1,em1,' kf2 em2',kf2,em2
          stop
        endif
        endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

          goto 5000
        endif

c....Sorry not implemented now!
          if(mstc(8).ge.3) then
           write(check(1),8000)srt,chaf(kc1,1),chaf(kc2,1)
 8000  format('srt=',g11.3,'GeV ',a10,' + ',a10)
           call jamerrm(1,1,'(jamcross:) Sorry not implemented')
          endif

          if(msel.eq.1) then
            call jamxadq(kf1,kf2,sig,sigel)
            sigel=sig
          else
            pare(3)=-10.d0
            kfb=kf2
            kfm=kf1
            call attflv(kfb,kfla1,kflb1)
            call attflv(kfm,kfla2,kflb2)
            sth=jamemjet(kfla1,kflb1)+jamemjet(kfla2,kflb2)
            if(srt.gt.sth) then
              ijet=1
            else
              call jamrmas2(kf1,kf2,kc1,kc2,srt,em1,em2,icon)
              if(icon.ne.0) ijet=-1
            endif
            sigin(1)=sig-sigel
            mchanel=1
            mabsrb=0
          endif

        goto 5000

c=============================
c...Meson + meson collisions
c=============================
      else if(icltyp.eq.3) then

        iz1=kchg(kc1,1)*isign(1,kf1)
        iz2=kchg(kc2,1)*isign(1,kf2)
        izt=iz1+iz2
        mchanel=1
        mabsrb=2
        sigin(1)=0.0d0  ! t-channel inel.
        sigin(2)=0.0d0  ! s-channel resonance formation
        sigin(3)=0.0d0  ! s-channel string formation
        sigres=0.0d0
        sigab=0.0d0
        sigabc=0.0d0

c.....Low energy pion - pion collisions.
c       id1=kchg(kc1,5)
c       id2=kchg(kc2,5)
c       ipair=jamcpair(id1,id2)
c       if(ipair.eq.jamcpair(id_pi,id_pi).and.srt.le.0.9) then
c          call jamcpipi(msel,srt,kf1,kf2,iz1,iz2,
c    $      em1,em2,sig,sigel,sigab,mchanel,mabsrb,ijet)
c          if(msel.eq.3.and.pare(3).le.0.0d0) return
c       endif

c....Calculate resonance cross sections.
        emrf=0d0
        if(srt.le.3.0d0)
     $    call jamxbw2(srt,pr,kf1,kf2,kc1,kc2,iz1,iz2,sigres,emrf,msel)
        sig=sigres+sigab
        if(msel.eq.3.and.pare(3).le.0.d0) then
          em1=srt
          em2=0.0d0
          return
        else if(msel.eq.2) then
          sigin(2)=sig
        endif

c...Additive quark cross section.
        if(srt.ge.emrf) then
          call jamxadq(kf1,kf2,siga,sigel)
c         sig=max(sig,siga)
          sig=max(sig,siga-sigel)  ! 2007/1/19
        endif
        sig=sig+sigel
        if(msel.eq.1) return

c...2010/8/24
        call attflv(kf1,ifla1,iflb1)
        call attflv(kf2,ifla2,iflb2)

c...Absorption impossible.
       if(izt.le.-6.or.izt.ge.6) goto 2000

        kfla1=max(ifla1,iflb1)
        kflb1=min(ifla1,iflb1)
        kfla2=max(ifla2,iflb2)
        kflb2=min(ifla2,iflb2)
        iann=0
        iann1=0
        iann2=0
        if(kfla1.eq.abs(kflb2)) iann1=1
        if(abs(kflb1).eq.kfla2) iann2=1
        if(iann1.eq.1.and.iann2.eq.1) then
          iann=1
          if(rn(0).gt.0.5d0) iann=2
        else if(iann1.eq.1) then
           iann=1
        else if(iann2.eq.1) then
           iann=2
        endif
        if(iann.eq.1) then
          call kfcnst(kfla2,kflb1,kfv1,0.0d0)
          sth=jamemjet(kfla2,kflb1)
          if(kfv1.eq.0)then
            call jamerrm(1,0,'(jamcross:)kfv1=0 at iann=1')
            goto 2000
          endif
        else if(iann.eq.2) then
          call kfcnst(kfla1,kflb2,kfv1,0.0d0)
          sth=jamemjet(kfla1,kflb2)
          if(kfv1.eq.0)then
            call jamerrm(1,0,'(jamcross:)kfv1=0 at iann=2')
            goto 2000
          endif
c         kfla=abs(kfla1)
c         kflb=abs(kflb2)
        else
          goto 2000
        endif
        if(jamcomp(kf1).eq.0) then
          write(mstc(38),*)'(jamcross:)m-m annihilation kf1=0??',kf1,kf2
          write(mstc(38),*)'kfla1 kflb1',kfla1,kflb1
          write(mstc(38),*)'kfla2 kflb2',kfla2,kflb2
          goto 2000
        endif

      sig1=max(0.0d0,sig-sigel-sigres)
c     mstr=0
c     mhvy=0
c     if(kfla.eq.3) mstr=mstr+1
c     if(kflb.eq.3) mstr=mstr+1
c     if(kfla.ge.4) mhvy=mhvy+1
c     if(kflb.ge.4) mhvy=mhvy+1
c.....Find s-channel contribution to the continuem.
c     sth=1.8d0+0.15d0*mstr+mhvy*1.0d0
c     sth=jamemjet(kfla1,kflb2)
      sths=sth+0.5d0
      if(srt.gt.sth) then
        if(srt.le.sths) then
          sigabc=sig1
        else
          sigabc=sig1*sths/srt
        endif
      endif
 
c....Monte Carlo for s-channel string formation.
        if(msel.eq.3) then
          xsig0=pare(3)
          pare(3)=pare(3)-sigabc
          if(pare(3).le.0.0d0) then
            kf1=kfv1
            kf2=0
            em1=srt
            em2=0.0d0
            ijet=2
            return
          endif
        else
          mabsrb=2
          sigin(3)=sigabc
        endif

2000   continue
      sigin(1)=max(0.0d0,sig-sigel-sigabc-sigres)

c....t-channel inelastic.
      if(msel.eq.3) then
         pare(3)=-10.d0
         sths=jamemjet(ifla1,iflb1)+jamemjet(ifla2,iflb2) !2010/8/24
         if(srt.ge.sths) then
           ijet=1
         else
           noutpa=0
           call jamrmas2(kf1,kf2,kc1,kc2,srt,em1,em2,icon)
           if(icon.ne.0) ijet=-1
         endif
         return
      endif
      return

c=============================
c...Anti-baryon + baryon  collisions.
c=============================
      else if(icltyp.eq.4) then

        call jamcabb(msel,srt,pr,kf1,kf2,kc1,kc2,
     $         em1,em2,sig,sigel,sigin,mchanel,mabsrb,mxchan,ijet)

      endif
      return

c...G-parity transformation.
 5000 continue
      if(anti) then
        kcx=jamcomp(kf1)
        kcy=jamcomp(kf2)
        if(kcx.ge.1) then
          if(kchg(kcx,3).eq.1) kf1=-kf1
        endif
        if(kcy.ge.1) then
          if(kchg(kcy,3).eq.1) kf2=-kf2
        endif
      endif

      end

c***********************************************************************

      subroutine jamcbb(msel,srt,pr,kf1,kf2,kc1,kc2,em1,em2,
     $                  mchanel,sig,sigel,sigin,ic,mxchan,icon)

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      dimension sigin(mxchan),sig1(10),sig0(10)
      parameter(rfac=1.0d0)

      ic=0
      iz1=kchg(kc1,1)*isign(1,kf1)
      iz2=kchg(kc2,1)*isign(1,kf2)
      snew=2*sqrt(parc(28)**2+pr**2)
      id1=kchg(kc1,5)
      id2=kchg(kc2,5)
      ipair=jamcpair(id1,id2)

c.....nucleon - nucleon
      if(ipair.eq.jamcpair(id_nucl,id_nucl)) then

          if( msel.eq.1) then
            call jamxnn(srt,iz1,iz2,sig,sigel)
          else

          if(kf1.eq.kf2) then
           ipn=1
          else
           ipn=2
          endif

c.......At least one pion production
c           if(kf1.eq.kf2) then
c             ipn=1
c             if(kf1.eq.2112) then
c               smin=parc(24)+parc(25)+parc(27)
c             else
c               smin=2*parc(25)+parc(26)
c             endif
c           else
c             ipn=2
c             smin=parc(24)+parc(25)+parc(26)
c           endif
            smin=2.0139999d0 ! thd0 in sub. jamxnnin
            if(srt.le.smin+parc(41)) then
              mchanel=0
              sigin(1)=0
              icon=-1
              return
            endif

c.....Resonance production cross sections.
            call jamxnnin(srt,sig1,ipn)
c 1) NN -> ND
c 2) NN -> NN*
c 3) NN -> DD
c 4) NN -> ND*
c 5) NN -> N*D
c 6) NN -> DD*
c 7) NN -> N*N*
c 8) NN -> N*D*
c 9) NN -> D*D*
            if(ipn.eq.1) then
              mchanel=17
              sigin(1)=0.25d0*sig1(1)   ! pd+
              sigin(2)=0.75d0*sig1(1)   ! nd++
              sigin(3)=sig1(2)          ! pp*
              sigin(4)=0.4d0*sig1(3)    ! d+d+
              sigin(5)=0.6d0*sig1(3)    ! d0d++
              sigin(6)=0.25d0*sig1(4)   ! pd*+
              sigin(7)=0.75d0*sig1(4)   ! nd*++
              sigin(8)=0.25d0*sig1(5)   ! p*d+
              sigin(9)=0.75d0*sig1(5)   ! n*d++
              sigin(10)=0.4d0*sig1(6)   ! d*d+
              sigin(11)=0.6d0*sig1(6)   ! d*d++
              sigin(12)=sig1(7)         ! p*p*
              sigin(13)=0.25d0*sig1(8)  ! p*d*+
              sigin(14)=0.75d0*sig1(8)  ! n*d*++
              sigin(15)=0.4d0*sig1(9)   ! d*+d*+
              sigin(16)=0.6d0*sig1(9)   ! d*0d*++
              sigin(17)=sig1(10)        ! s-wave
            else
              mchanel=19
              call jamxnnin(srt,sig0,0)
              sigin(1)=0.25d0*sig1(1)                   ! nd+
              sigin(2)=0.25d0*sig1(1)                   ! pd0
              sigin(3)=0.25d0*(sig1(2)+sig0(2))         ! np*
              sigin(4)=0.25d0*(sig1(2)+sig0(2))         ! n*p
              sigin(5)=0.05d0*sig1(3)+0.25d0*sig0(3)    ! d0d+
              sigin(6)=0.45d0*sig1(3)+0.25d0*sig0(3)    ! d-d++
              sigin(7)=0.25d0*sig1(4)                   ! nd*+
              sigin(8)=0.25d0*sig1(4)                   ! pd*0
              sigin(9)=0.25d0*sig1(5)                   ! n*d+
              sigin(10)=0.25d0*sig1(5)                  ! p*d0
              sigin(11)=0.05d0*sig1(6)+0.25d0*sig0(6) ! d0d*+
              sigin(12)=0.45d0*sig1(6)+0.25d0*sig0(6) ! d-d*++
              sigin(13)=0.5d0*(sig1(7)+sig0(7))        ! n*p*
              sigin(14)=0.25d0*sig1(8)                  ! n*d*+
              sigin(15)=0.25d0*sig1(8)                  ! p*d*0
              sigin(16)=0.05d0*sig1(9)+0.25d0*sig0(9) ! d*0d*+
              sigin(17)=0.45d0*sig1(9)+0.25d0*sig0(9) ! d*-d*++
              sigin(18)=0.5d0*sig0(10)                  ! s-wave
              sigin(19)=0.5d0*sig0(10)                  ! s-wave
            endif

c....pp->ppKaK (not implemeted now).
cc         if(srt.le.6.5d0) then
cc           sigin(mchanel+1)=sigkk(srt)
cc         endif

c.....Get outgoing masses and types
            if(msel.eq.3) then
               assign 10 to label
               goto 1000
 10            continue
               if(ic.eq.0) return

c......s-wave pion production
               mste(3)=0
               if(mstc(66).eq.1) then
                 if(ipn.eq.1.and.kf1.eq.2212.and.ic.eq.17) mste(3)=1
                 if(ipn.eq.1.and.kf1.eq.2112.and.ic.eq.17) mste(3)=2
                 if(ipn.eq.2.and.(ic.ge.18.and.ic.le.19)) mste(3)=3
               endif

               if(mste(3).eq.0) then
                 call jambmas(1,srt,pr,ic,ipn,kc1,kc2,kf1,kf2,iz1,iz2,
     $                               em1,em2,icon)
                 if(icon.ne.0) return
               endif
            endif

          endif

        return

c...Resonance cross sections.
      else

c.....Skip calculate resonance cross sections.
         if(srt.ge.100d0) then
           ic=0
           mchanel=0
           sigab=0.0d0
           goto 1200
         endif

c...First get inel. cross sections at equal c.m.s. momentum.
c.......N + N* or N* + N*  inelastic collisions
        if( ( ipair .eq. jamcpair(id_nucl,id_nucls) ) .or.
     &           ( ipair .eq. jamcpair(id_nucls,id_nucls) )  ) then

            call jamcbb1(srt,snew,em1,em2,pr,kf1,kf2,
     &          kc1,kc2,iz1,iz2,mchanel,sigin,ipn,mxchan)

c......Test: reduce resonance cross sectins.
            do i=1,mchanel-1
              sigin(i)=rfac*sigin(i)
            end do

            sigab=sigin(mchanel)
c.....Get outgoing masses and types
            if(msel.eq.3) then
               assign 20 to label
               goto 1000
 20            continue
               if(ic.ne.0) then
                 call jambmas(1,srt,pr,ic,ipn,kc1,kc2,kf1,kf2,iz1,iz2,
     $                                  em1,em2,icon)
                 if(icon.ne.0) return
               else
                 return
               endif
            endif

c.......N(*) + Delta(*) inelastic collisions
        else if( ( ipair .eq. jamcpair(id_nucl,id_delt)  ) .or.
     &         ( ipair .eq. jamcpair(id_nucl,id_delts) ) .or.
     &         ( ipair .eq. jamcpair(id_nucls,id_delt) ) .or.
     &         ( ipair .eq. jamcpair(id_nucls,id_delts))  ) then

          call jamcbb2(srt,snew,em1,em2,pr,kf1,kf2,kc1,kc2,
     &       id1,id2,iz1,iz2,sigin,mchanel,ipn,mxchan)
c......Test: reducce resonance cross sectins.
            do i=1,mchanel-1
              sigin(i)=rfac*sigin(i)
            end do

            sigab=sigin(mchanel)
c.....Get outgoing masses and types
            if(msel.eq.3) then
               assign 30 to label
               goto 1000
 30            continue
               if(ic.ne.0) then
                 call jambmas(2,srt,pr,ic,ipn,kc1,kc2,kf1,kf2,iz1,iz2,
     $                                      em1,em2,icon)
                 if(icon.ne.0) return
               else
                 return
               endif
            endif


c.....Delta(*) - Delta(*)  inelastic collisions
        else if( ( ipair .eq. jamcpair(id_delt,id_delt) ) .or.
     &         ( ipair .eq. jamcpair(id_delt,id_delts) ) .or.
     &         ( ipair .eq. jamcpair(id_delts,id_delts) )  ) then

          call jamccbb3(srt,snew,em1,em2,pr,kf1,kf2,kc1,kc2,id1,id2,
     &                   iz1,iz2,sigin,mchanel,ipn,mxchan)
c......Test: reducce resonance cross sectins.
            do i=1,mchanel-1
              sigin(i)=rfac*sigin(i)
            end do

            sigab=sigin(mchanel)
c.....Get outgoing masses and types
            if(msel.eq.3) then
               assign 40 to label
               goto 1000
 40            continue
               if(ic.ne.0) then
                 call jambmas(3,srt,pr,ic,ipn,kc1,kc2,kf1,kf2,iz1,iz2,
     $                                      em1,em2,icon)
                 if(icon.ne.0) return
               else
                 return
               endif
            endif

        else

        if(mstc(8).ge.1) then
          isn1=(3-isign(1,kf1))/2
          isn2=(3-isign(1,kf2))/2
          write(check(1),8000)srt,chaf(kc1,isn1),chaf(kc2,isn2)
 8000     format('srt=',g10.3,'GeV',a10,' + ',a10)
          call jamerrm(30,1,'(jamcbb) invalid code')
        endif

        endif
      endif


 1200 if(msel.eq.1) then

c...Elastic from NN cross sections.
c       if(iz1.eq.iz2) then
c         call jamxnn(snew,3,3,sig,sigel)
c       else
c         call jamxnn(snew,3,0,sig,sigel)
c       endif
        call jamxnn(snew,3,3,sigt1,sigel1)
        call jamxnn(snew,3,0,sigt2,sigel2)
        sig=0.5d0*(sigt1+sigt2)
        sigel=0.5d0*(sigel1+sigel2)

        if(snew.le.3.6d0) then
          sigint=0.0d0
          do i=1,mchanel
            sigint=sigint+sigin(i)
          end do
          sig=sigel+sigint

        else
c...Add absorption cross sections.
          sig=sig+sigab
        endif

      endif

      return
1000  continue

c========================================
c...Select possible inelastic channel.
c========================================

      ic=0
      itry=0
2100  continue
      itry=itry+1
      if(itry.ge.100) then
        write(check(1),8100)srt,pr,ic,mchanel
        write(check(2),8200)kf1,kf2,em1,em2
        write(check(3),8210)sigtt
        write(check(4),8300)(sigin(i),i=1,3)
 8100   format('srt',g10.3,' pr=',g10.3,' ic=',i3,' mchanel=',i3)
 8200   format('kf1',i10,' kf2=',i10,' em1=',g10.3,' em2=',g10.3)
 8210   format('sigtt=',g10.3)
 8300   format('sigin=',3(g10.3,1x))
        call jamerrm(1,3,'(jamcbb:) No inel. branch found')
        icon=1
        return
      endif

      ic=0
      sigtt=0.0d0
      do i=1,mchanel
        sigtt=sigtt+sigin(i)
        if(sigin(i).lt.0.0d0) then
        write(check(1),8400)kf1,em1,chaf(kc1,(3-isign(1,kf1))/2)
        write(check(2),8450)kf2,em2,chaf(kc2,(3-isign(1,kf2))/2)
        write(check(3),8460)i,sigin(i)
 8400   format('kf1 em1',i9,1x,g10.3,1x,a16)
 8450   format('kf2 em2',i9,1x,g10.3,1x,a16)
 8460   format('i sigin(i)',i9,1x,g10.3)
        sigin(i)=0.0d0
        call jamerrm(1,3,'inel. cross section has negative value')
        endif

        pare(3)=pare(3)-sigin(i)
        if(pare(3).le.0.0d0) then
          ic=i
          goto 300
        end if
      end do
 
c....Inelastic closed due to low energy.
      if(sigtt.lt.1d-6) then
        icon=1
        return
      endif

c....String impossible at low enrgy.
      if(snew.le.3.5d0) then
        pare(3)=sigtt*rn(0)
        goto 2100
      endif

300   continue
      goto label

      end

c**********************************************************************

      subroutine jamcbb1(srt,snew,em1,em2,pr,kf1,kf2,
     &       kc1,kc2,iz1,iz2,mchanel,sigin,icc,mxchan)

c...Purpose : to give cross sections for n + n* or n* + n* collisions.
c---------------------------------------------------------------------*
c...Outputs
c...sigin(): inel. cross sections.
c...mchanel: max. inel. channel.
c...at the moment dn -> x  cross sections are simply equated to
c...t=1 cross sections at equal cms momentum
c---------------------------------------------------------------------*
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      dimension sigin(mxchan),sig1(10),sig0(10)

c...See whether n + n* or n* + n* collision.
      if( kf1.eq.2112 .or. kf1.eq.2212 ) then
        isig=3
        call jamdetb1(1,kc2,srt,em1,em2,pr,detbal,nhlf)
      else if( kf2.eq.2112 .or. kf2.eq.2212 ) then
        isig=3
        call jamdetb1(1,kc1,srt,em1,em2,pr,detbal,nhlf)
      else
        call jamdetb2(1,kc1,kc2,srt,em1,em2,pr,detbal,nhlf)
        isig=7
      endif

c...Get cross sections from nn collision.
      if( iz1 .eq. iz2 ) then
        icc=1
        faci=0.5d0
        call jamxnnin(srt,sig1,1)
        sigab=sig1(isig)
        call jamxnnin(snew,sig1,1)
      else
        icc=2
        faci=1.d0
        call jamxnnin(srt,sig1,2)
        call jamxnnin(srt,sig0,0)
        sigab=0.25d0*(sig0(isig)+sig1(isig))
        if(isig.eq.3) sigab=0.5d0*sigab
        call jamxnnin(snew,sig1,2)
        call jamxnnin(snew,sig0,0)
      end if

            if(icc.eq.1) then
              mchanel=18
              sigin(1)=0.25d0*sig1(1)  ! pd+
              sigin(2)=0.75d0*sig1(1)  ! nd++
              sigin(3)=sig1(2)         ! pp*+
              sigin(4)=0.4d0*sig1(3)   ! d+d+
              sigin(5)=0.6d0*sig1(3)   ! d0d++
              sigin(6)=0.25d0*sig1(4)  ! pd*+
              sigin(7)=0.75d0*sig1(4)  ! nd*++
              sigin(8)=0.25d0*sig1(5)  ! p*d+
              sigin(9)=0.75d0*sig1(5)  ! n*d++
              sigin(10)=0.4d0*sig1(6)  ! d*d+
              sigin(11)=0.6d0*sig1(6)  ! d*d++
              sigin(12)=sig1(7)        ! p*p*
              sigin(13)=0.25d0*sig1(8) ! p*d*+
              sigin(14)=0.75d0*sig1(8) ! n*d*++
              sigin(15)=0.4d0*sig1(9)  ! d*+d*+
              sigin(16)=0.6d0*sig1(9)  ! d*0d*++
              sigin(17)=sig1(10)       ! s-wave
              sigin(18)=0.0d0          ! b*b*->nn
            else
              mchanel=20
              sigin(1)=0.25d0*sig1(1)                    ! nd+
              sigin(2)=0.25d0*sig1(1)                    ! pd0
              sigin(3)=0.25d0*(sig1(2)+sig0(2))          ! np*
              sigin(4)=0.25d0*(sig1(2)+sig0(2))          ! n*p
              sigin(5)=0.05d0*sig1(3)+0.25d0*sig0(3)   ! d0d+
              sigin(6)=0.45d0*sig1(3)+0.25d0*sig0(3)   ! d-d++
              sigin(7)=0.25d0*sig1(4)                    ! nd*+
              sigin(8)=0.25d0*sig1(4)                    ! pd*0
              sigin(9)=0.25d0*sig1(5)                    ! n*d+
              sigin(10)=0.25d0*sig1(5)                   ! p*d0
              sigin(11)=0.05d0*sig1(6)+0.25d0*sig0(6)  ! d0d*+
              sigin(12)=0.45d0*sig1(6)+0.25d0*sig0(6)  ! d-d*++
              sigin(13)=0.5d0*(sig1(7)+sig0(7))         ! n*p*
              sigin(14)=0.25d0*sig1(8)                   ! n*d*+
              sigin(15)=0.25d0*sig1(8)                   ! p*d*0
              sigin(16)=0.05d0*sig1(9)+0.25d0*sig0(9)  ! d*0d*+
              sigin(17)=0.45d0*sig1(9)+0.25d0*sig0(9)  ! d*-d*++
              sigin(18)=0.5d0*sig0(10)                   ! s-wave
              sigin(19)=0.5d0*sig0(10)                   ! s-wave
              sigin(20)=0.0d0                            ! b*b*->nn
            endif

c...Spin factor.
      spin=dble(mod(abs(kf1),10)*mod(abs(kf2),10))
      facsp=4.d0/spin

c...Calculate N(*)N(*) --> NN absorption cross section
c...from detailed balance.
      pr2new=0.25d0*srt*srt-parc(28)*parc(28)
      sigin(mchanel)=sigab*facsp*faci*pr2new*detbal

c...TEST
c     sigin(mchanel)=2.0*sigin(mchanel)
c     pare(4)=siginn(1)+siginn(2)+siginn(4)+siginn(10)
c     pare(5)=siginn(3)+siginn(5)+siginn(6)+siginn(7)
c    $       +sigin(8)+sigin(9)

      end

c***********************************************************************

      subroutine jamcbb2(srt,snew,em1,em2,pr,kf1,kf2,kc1,kc2,
     &       id1,id2,iz01,iz02,sigin,maxb,icc,mxchan)

c...Purpose : to give cross sections for n(*) + d(*) collisions.
c
c   srt (GeV) : c.m. energy of the colliding pair.
c   id1, id2:  particle ID
c   ires1, ires2: resonance number
c
c ... at the moment dn -> x  cross sections are simply equated to
c     t=1 cross sections at equal cms momentum
c=======================================================================

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      dimension sigin(mxchan),sig1(10)

      iz1=iz01/3
      iz2=iz02/3
      itz=iz1+iz2

c...1 is a nucleon(*) 2 is a delta(*)
      if(id1.eq.id_nucl .or. id1.eq.id_nucls) then
        izn=iz1
        izd=iz2
        iresn=1
        iresd=1
        if(id1.eq.id_nucl) iresn=0
        if(id2.eq.id_delt) iresd=0
        kcd=kc2
        emn=em1
        emd=em2
      else
        izn=iz2
        izd=iz1
        iresn=1
        iresd=1
        if(id1.eq.id_delt) iresd=0
        if(id2.eq.id_nucl) iresn=0
        kcd=kc1
        emn=em2
        emd=em1
      endif

      if(itz.ge.3.or.itz.le.-1) then
        idn=0
        faci=0.0d0
        icc=1
        detbal=0.d0
      else

        if(iresn.eq.0) then
          if(iresd.eq.0) then
            idn = 1   !   delta n ==> n n
            isig=1
            call jamdetb1(1,kcd,srt,emd,emn,pr,detbal,nhlf)

c.....delta* n ==> n n
          else
            idn=2
            isig=4
            call jamdetb1(1,kcd,srt,emd,emn,pr,detbal,nhlf)
          endif
        else

c.....d(d*) n* ==> n n
          call jamdetb2(1,kc1,kc2,srt,em1,em2,pr,detbal,nhlf)
          if(iresd.eq.0) then
            idn=3   !   delta n* ==> n n
            isig=5
          else
            idn=4   !   delta* n* ==> nn 
            isig=8
          end if
        endif

      end if 

      if(idn.eq.0) then
        sigab=0
c....pd- --> nn  16:nd++ --> pp
      else if( (izn.eq.1.and.izd.eq.-1) .or. 
     &         (izn.eq.0.and.izd.eq.2)  ) then

        faci=0.5d0
        icc=1
        call jamxnnin(srt,sig1,1)
        sigab=0.75d0*sig1(isig)

c....nd0- --> nn  12:pd+ --> pp
      else if( (izn.eq.0.and.izd.eq.0) .or. 
     &          (izn.eq.1.and.izd.eq.1)  ) then

        faci=0.5d0
        icc=1
        call jamxnnin(srt,sig1,1)
        sigab=0.25d0*sig1(isig)

c....pd0 --> pn  11:nd+ --> pn
      else if( (izn.eq.1.and.izd.eq.0) .or. 
     &          (izn.eq.0.and.izd.eq.1)  ) then

        faci=1.0d0
        icc=2
        call jamxnnin(srt,sig1,2)
        sigab=0.25d0*sig1(isig)
      end if

c...Get absorption cross sections
      if(idn.ne.0) then
c       if(vfd.le.0.0d0) vfd=1.0
c       sigab=facsp*faci*vfd*pr2new*sigab/(pr*pr)
        facsp=4.d0/dble(mod(abs(kf1),10)*mod(abs(kf2),10))
        pr2new=0.25d0*srt*srt-parc(28)*parc(28)
        sigab=facsp*faci*pr2new*sigab*detbal
      endif

c...Calculate nn cross section.
      call jamxnnin(snew,sig1,icc)
      sigin(11)=sigab
      do i=1,10
      sigin(i)=sig1(i)
      end do
c...Include s-wave pion into NN* branch
      sigin(2)=sigin(2)+sig1(10)
      sigin(10)=0.0d0      ! assume no s-wave pion
      if(idn.eq.0) then
        sigin(2)=0.0d0
        sigin(7)=0.0d0
      endif
      maxb=11

c...TEST
c     sigin(11)=2.0*sigin(11)
c     pare(4)=sig1(1)+sig1(2)+sig1(4)+sig1(10)
c     pare(5)=sig1(3)+sig1(5)+sig1(6)+sig1(7)
c    $       +sig1(8)+sig1(9)

      end

c**********************************************************************

      subroutine jamccbb3(srt,snew,em1,em2,pr,kf1,kf2,kc1,kc2,id1,id2,
     &                   iz01,iz02,sigin,maxb,icc,mxchan)

c...Purpose : to give cross sections for d(*) + d(*) collisions.
c...at the moment dn -> x  cross sections are simply equated to
c...t=1 cross sections at equal cms momentum

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      dimension sigin(mxchan),sig1(10),sig0(10)

      iz1=iz01/3
      iz2=iz02/3
      iztot=iz1+iz2
      izmin=min(iz1,iz2)
      izmax=max(iz1,iz2)

      if(id1.eq.id_delt.and.id2.eq.id_delt) then
        isig=3
      else if( (id1.eq.id_delt.and.id2.eq.id_delts)
     $     .or.(id1.eq.id_delts.and.id2.eq.id_delt) ) then
        isig=6
      else if(id1.eq.id_delts.and.id2.eq.id_delts) then
        isig=9
      else
        write(check(1),'(''id1 id2'',i4,1x,i4)')id1,id2
        call jamerrm(30,1,
     $     '(jamcbb3:) error in d+d collisions id1 id2:')
      endif

      nnyes=1
      ndyes=1
c...Identify type of collision.
      if(iztot.eq.4.or.iztot.eq.-2) then
        nnyes=0
        ndyes=0
      else if( iztot .ge. 3 .or. iztot.le.-1 ) then

        faci=0.0d0
        icc=1
        nnyes=0

c...d+d+->pp/d0d0->nn
      else if(izmin.eq.izmax) then  

        faci=0.5d0
        icc=1
        call jamxnnin(srt,sig1,1)
        sigab=0.4d0*sig1(isig)

c...d0d++ -> pp /d+d- -> nn
      else if( (izmin.eq.0.and.izmax.eq.2) .or.
     &             (izmin.eq.-1.and.izmax.eq.1) ) then  

        faci=0.5d0
        icc=1
        call jamxnnin(srt,sig1,1)
        sigab=0.6d0*sig1(isig)

c...d0d+ -> np
      else if(izmin.eq.0.and.izmax.eq.1) then  

        faci=1.0d0
        icc=2
        call jamxnnin(srt,sig1,2)
        call jamxnnin(srt,sig0,0)
        sigab=0.025d0*sig1(isig)+0.125d0*sig0(isig)

c...d-d++ -> np
      else if(izmin.eq.-1.and.izmax.eq.2) then  

        faci=1.0d0
        icc=2
        call jamxnnin(srt,sig1,2)
        call jamxnnin(srt,sig0,0)
        sigab=0.225d0*sig1(isig)+0.125d0*sig0(isig)

      else
        write(check(1),'(''iz1 iz2'',i3,1x,i3)')iz1,iz2
        call jamerrm(30,1,'(jamcbb3:) error in d+d collisions iz1 iz2')
      end if


c...dd -> nn  from the detailed valance 

      if(nnyes.eq.1) then
        call jamdetb2(1,kc1,kc2,srt,em1,em2,pr,detbal,nhlf)
        facsp=4.d0/dble(mod(abs(kf1),10)*mod(abs(kf2),10))
        pr2new=0.25d0*srt*srt-parc(28)*parc(28)
        sigab=facsp*faci*pr2new*sigab*detbal
      else
        sigab=0.0d0
      end if

      sigin(11)=sigab
      call jamxnnin(snew,sig1,1)
      do i=1,10
        sigin(i)=sig1(i)
      end do

c...Include s-wave pion into NN* branch
      sigin(2)=sigin(2)+sig1(10)
      sigin(10)=0.0d0  ! assume no s-wave pion 
      if(nnyes.eq.0) then
        sigin(2)=0.0d0
        sigin(7)=0.0d0
      endif
      if(ndyes.eq.0) then
        sigin(1)=0.0d0
        sigin(4)=0.0d0
        sigin(5)=0.0d0
        sigin(8)=0.0d0
      endif
      maxb=11

c...TEST
c     sigin(11)=2.0*sigin(11)
c     pare(4)=sig1(1)+sig1(2)+sig1(4)+sig1(10)
c     pare(5)=sig1(3)+sig1(5)+sig1(6)+sig1(7)
c    $       +sig1(8)+sig1(9)

      end

c***********************************************************************

      subroutine jamcbbs(msel,kf1,kf2,kc1,kc2,istr,srt,
     $ sig,sigel,sigin,em1,em2,mchanel,mxchan,icon)

c....Handle low energy Lambda-N/Sigma-N/Xi-N/LL... collisions.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      dimension sigin(mxchan),kfy(5,2),isig(5)

c...Initialize some values.
      icon=0
      mchanel=0
      id1=kchg(kc1,5)
      id2=kchg(kc2,5)
      ipair=jamcpair(id1,id2)
      iz01=kchg(kc1,1)/3
      iz02=kchg(kc2,1)/3
      izt=iz01+iz02
      iz1=iz01
      iz2=iz02
      isigt=0
      isige=0
      do i=1,5
        isig(i)=0
      end do


c...S=-1 BB collisions.
      if(istr.eq.-1) then

      isw=0
      if(id1.eq.id_nucl) then
        iz1=iz02
        iz2=iz01
        isw=1
      endif

c....Lambda + nucleon.
      if(ipair.eq.jamcpair(id_lamb,id_nucl)) then

        if(izt.eq.1) then           ! Lambda p
          if(msel.eq.1) then
            isigt=1
            isige=2
          else
            mchanel=2
            isig(1)=3
            isig(2)=4
            kfy(1,1)=3222  ! -> Sigma+ n
            kfy(1,2)=2112
            kfy(2,1)=3212  ! -> Sigma0 p
            kfy(2,2)=2212
          endif
        else if(izt.eq.0) then      ! Lambda n
          if(msel.eq.1) then
            isigt=5
            isige=6
          else
            mchanel=2
            isig(1)=7
            isig(2)=8
            kfy(1,1)=3212 ! -> Sigma0 n
            kfy(1,2)=2112
            kfy(2,1)=3112 ! -> Sigma- p
            kfy(2,2)=2212
          endif
        endif

c...Sigma + nuucleon.
      else if(ipair.eq.jamcpair(id_sigm,id_nucl)) then

        if(izt.eq.-1) then           ! Sigma- n
          if(msel.eq.1) then
            isigt=9
            isige=0
          else
            mchanel=0
          endif
        else if(izt.eq.0) then
          if(iz1.eq.0) then ! Sigma0 n
            if(msel.eq.1) then
              isigt=15
              isige=16
            else
              mchanel=2
              isig(1)=17
              isig(2)=18
              kfy(1,1)=3122 ! -> Lambda n
              kfy(1,2)=2112
              kfy(2,1)=3112 ! -> Sigma- p
              kfy(2,2)=2212
            endif
          else if(iz1.eq.-1) then  ! Sigma- p
            if(msel.eq.1) then
              isigt=11
              isige=12
            else
              mchanel=2
              isig(1)=13
              isig(2)=14
              kfy(1,1)=3122  ! -> Lambda n
              kfy(1,2)=2112
              kfy(2,1)=3212  ! -> Sigma0 n
              kfy(2,2)=2112
            endif
          endif
        else if(izt.eq.1) then
          if(iz1.eq.1) then          ! Sigma+ n
            if(msel.eq.1) then
              isigt=23
              isige=24
            else
              mchanel=2
              isig(1)=25
              isig(2)=26
              kfy(1,1)=3122  ! -> Lambda p
              kfy(1,2)=2212
              kfy(2,1)=3212  ! -> Sigma+ n
              kfy(2,2)=2212
            endif
          else if(iz1.eq.0) then       ! Sigma0 p
            if(msel.eq.1) then
              isigt=19
              isige=20
            else
              mchanel=2
              isig(1)=21
              isig(2)=22
              kfy(1,1)=3122  ! -> Lambda p
              kfy(1,2)=2212
              kfy(2,1)=3222  ! -> Sigma+ n
              kfy(2,2)=2112
            endif
          endif

        else if(izt.eq.2) then       ! Sigma+ p

          if(msel.eq.1) then
            isigt=10
            isige=0
          else
            mchanel=0
          endif
        endif

      endif


c...Get total and elastic cross sections.
      if(msel.eq.1) then
        if(isigt.eq.0.and.isige.eq.0) goto 1000
        call jamsigS1(sig,isigt,srt) 
        sigel=sig
        if(isige.ge.1) call jamsigS1(sigel,isige,srt) 

c...Get inelastic cross sections.
      else
        if(mchanel.eq.0) goto 1000
        sigin(1)=0d0
        sigin(2)=0d0
        do mch=1,mchanel
          if(isig(mch).ge.1) call jamsigS1(sigin(mch),isig(mch),srt) 
        end do

c...Monte carlo sampling.
        if(mchanel.ge.1.and.msel.eq.3) then
          ic=0
          do i=1,mchanel
            pare(3)=pare(3)-sigin(1)
            if(pare(3).le.0.0d0) then
              ic=i
              goto 100
            endif
          end do
          ic=mchanel
          if(mstc(8).ge.1.and.pare(3).gt.0.0d0)
     $  call jamerrm(1,0,'(jamcbbs:)S=-1 inel')
          pare(3)=-10d0
  100     continue
          if(isw.eq.0) then
            kf1=kfy(ic,1)
            kf2=kfy(ic,2)
          else
            kf1=kfy(ic,2)
            kf2=kfy(ic,1)
          endif
            em1=pmas(jamcomp(kf1),1)
            em2=pmas(jamcomp(kf2),1)
          return
        endif
      endif
      return


c...S=-2 BB collisions.
      else if(istr.eq.-2) then
 
        isw=0
        if(id1.eq.id_nucl.or.id1.eq.id_lamb) isw=1

c...Lambda Lambda ingoing.
        if(ipair.eq.jamcpair(id_lamb,id_lamb)) then

          if(msel.eq.1) then
            isigt=23
            isige=24
          else
            mchanel=4
            isig(1)=25
            isig(2)=26
            isig(3)=27
            isig(4)=28
            kfy(1,1)=3322  ! -> Xi0 n
            kfy(1,2)=2112
            kfy(2,1)=3312  ! -> Xi- p
            kfy(2,2)=2212
            kfy(3,1)=3212  ! -> Sigma0 Sigma0
            kfy(3,2)=3212
            kfy(4,1)=3222  ! -> Sigma+ Sigma-
            kfy(4,2)=3112
          endif

c.....Lambda Sigma ingoing.
        else if(ipair.eq.jamcpair(id_lamb,id_sigm)) then

          if(izt.eq.-1) then     ! Lambda Sigma-
            if(msel.eq.1) then
              isigt=29
              isige=30
            else
              mchanel=2
              isig(1)=31
              isig(2)=32
              kfy(1,1)=3312  ! Xi- n
              kfy(1,2)=2112
              kfy(2,1)=3212  ! Sigma0 Sigma-
              kfy(2,2)=3112
            endif
          else if(izt.eq.0) then ! Lambda Sigma0
            if(msel.eq.1) then
              isigt=33
              isige=34
            else
              mchanel=3
              isig(1)=35
              isig(2)=36
              isig(3)=37
              kfy(1,1)=3322  ! Xi0 n
              kfy(1,2)=2112
              kfy(2,1)=3312  ! Xi- p
              kfy(2,2)=2212
              kfy(3,1)=3222  ! S+S-
              kfy(3,2)=3112
            endif
          else if(izt.eq.1) then ! Lambda Sigma+
            if(msel.eq.1) then
              isigt=38
              isige=39
            else
              mchanel=2
              isig(1)=40
              isig(2)=41
              kfy(1,1)=3322  ! Xi0 p
              kfy(1,2)=2212
              kfy(2,1)=3222  ! S+S0
              kfy(2,2)=3212
            endif
          endif

c....Sigma Sigma ingoing.
        else if(ipair.eq.jamcpair(id_sigm,id_sigm)) then

          if(izt.eq.-2) then               ! Sigma- Sigma-
            if(msel.eq.1) then
              isigt=42
              isige=0
            else
              mchanel=0
            endif
          else if(izt.eq.-1) then          ! Sigma- Sigma0
            if(msel.eq.1) then
              isigt=43
              isige=44
            else
              mchanel=2
              isig(1)=45
              isig(2)=46
              kfy(1,1)=3312  ! -> Xi- n
              kfy(1,2)=2112
              kfy(2,1)=3122  ! -> LS-
              kfy(2,2)=3112
            endif
          else if(izt.eq.0) then  
            if(iz1.eq.0.and.iz2.eq.0) then ! Sigma0 Sigma0
              if(msel.eq.1) then
                isigt=47
                isige=48
              else
                mchanel=5
                isig(1)=49
                isig(2)=50
                isig(3)=51
                isig(4)=52
                isig(5)=53
                kfy(1,1)=3122  ! -> Lambda Lambda
                kfy(1,2)=3122
                kfy(2,1)=3322  ! -> Xi0 n
                kfy(2,2)=2112
                kfy(3,1)=3312  ! -> Xi- p
                kfy(3,2)=2212
                kfy(4,1)=3122  ! -> Lambda Sigma0
                kfy(4,2)=3212
                kfy(5,1)=3222  ! -> Sigma+ Sigma-
                kfy(5,2)=3112
              endif
            else                           ! Sigma+ Sigma-
              if(msel.eq.1) then
                isigt=54
                isige=55
              else
                mchanel=5
                isig(1)=56
                isig(2)=57
                isig(3)=58
                isig(4)=59
                isig(5)=60
                kfy(1,1)=3122  ! -> Lambda Lambda
                kfy(1,2)=3122
                kfy(2,1)=3322  ! -> Xi0 n
                kfy(2,2)=2112
                kfy(3,1)=3312  ! -> Xi- p
                kfy(3,2)=2212
                kfy(4,1)=3122  ! -> Lambda Sigma0
                kfy(4,2)=3212
                kfy(5,1)=3212  ! -> Sigma0 Sigma0
                kfy(5,2)=3212
              endif
            endif
          else if(izt.eq.1) then           ! Sigma+ Sigma0
            if(msel.eq.1) then
              isigt=61
              isige=62
            else
              mchanel=2
              isig(1)=63
              isig(2)=64
              kfy(1,1)=3322  ! -> Xi0 p
              kfy(1,2)=2212
              kfy(2,1)=3122  ! -> LS+
              kfy(2,2)=3222
            endif
          else if(izt.eq.2) then           ! Sigma+ Sigma+
            if(msel.eq.1) then
              isigt=65
              isige=0
            else
              mchanel=0
            endif
          endif

c....Xi N ingoing.
        else if(ipair.eq.jamcpair(id_xi,id_nucl)) then

          if(izt.eq.-1) then                ! Xi- n
            if(msel.eq.1) then
              isigt=1
              isige=2
            else
              mchanel=2
              isig(1)=3
              isig(2)=4
              kfy(1,1)=3122  ! -> Lambda Sigma-
              kfy(1,2)=3112
              kfy(2,1)=3212  ! -> Sigma0 Sigma-
              kfy(2,2)=3112
            endif
          else if(izt.eq.0) then
            if(iz1.eq.0.and.iz2.eq.0) then  !  Xi0 n
              if(msel.eq.1) then
                isigt=12
                isige=13
              else
                mchanel=5
                isig(1)=14
                isig(2)=15
                isig(3)=16
                isig(4)=17
                isig(5)=18
                kfy(1,1)=3122  ! -> Lambda Lambda
                kfy(1,2)=3122
                kfy(2,1)=3312  ! -> Xi- p
                kfy(2,2)=2212
                kfy(3,1)=3122  ! -> Lambda Sigma0
                kfy(3,2)=3212
                kfy(4,1)=3222  ! -> Sigma+ Sigma-
                kfy(4,2)=3112
                kfy(5,1)=3212  ! -> Sigma0 Sigma0
                kfy(5,2)=3212
              endif
            else                            !  Xi- p
              if(msel.eq.1) then
                isigt=5
                isige=6
              else
                mchanel=5
                isig(1)=7
                isig(2)=8
                isig(3)=9
                isig(4)=10
                isig(5)=11
                kfy(1,1)=3122  ! -> Lambda Lambda
                kfy(1,2)=3122
                kfy(2,1)=3322  ! -> Xi0 n
                kfy(2,2)=2112
                kfy(3,1)=3122  ! -> Lambda Sigma0
                kfy(3,2)=3212
                kfy(4,1)=3222  ! -> Sigma+ Sigma-
                kfy(4,2)=3112
                kfy(5,1)=3212  ! -> Sigma0 Sigma0
                kfy(5,2)=3212
              endif
            endif
          else if(izt.eq.1) then            ! Xi0 p
            if(msel.eq.1) then
              isigt=19
              isige=20
            else
              mchanel=2
              isig(1)=21
              isig(2)=22
              kfy(1,1)=3122  ! -> Lambda Sigma+
              kfy(1,2)=3222
              kfy(2,1)=3222  ! -> Sigma+ Sigma0
              kfy(2,2)=3212
            endif
          endif
        endif

c...Get total and elastic cross sections.
      if(msel.eq.1) then
        if(isigt.eq.0.and.isige.eq.0) goto 1000
        call jamsigS2(sig,isigt,srt) 
        sigel=sig
        if(isige.ge.1) call jamsigS2(sigel,isige,srt) 
        return

c...Get inelastic cross sections.
      else
        if(mchanel.eq.0) goto 1000
        do i=1,5
        sigin(i)=0d0
        end do
        do mch=1,mchanel
          call jamsigS2(sigin(mch),isig(mch),srt) 
        end do

c...Monte carlo sampling.
        if(mchanel.ge.1.and.msel.eq.3) then
          ic=0
          do i=1,mchanel
            pare(3)=pare(3)-sigin(i)
            if(pare(3).le.0d0) then
              ic=i
              goto 200
            endif
          end do
          pare(3)=-10d0
          ic=mchanel
200       continue
          if(isw.eq.0) then
            kf1=kfy(ic,1)
            kf2=kfy(ic,2)
          else
            kf1=kfy(ic,2)
            kf2=kfy(ic,1)
          endif
          em1=pmas(jamcomp(kf1),1)
          em2=pmas(jamcomp(kf2),1)
          return
        endif
        return

      endif
      endif

 1000 icon=1

      end

c***********************************************************************

      subroutine jamcmbs0(msel,srt,pr,istr,kfb,kfm,kcb,kcm,izb,izm,
     $ emmes,embar,kf1,kf2,em1,em2,
     $ sig,sigel,sigin,
     $ mabsrb,mchanel,mxchan,ijet,jswap,icon)

c...Treat strangeness S=0 meson-baryon cross section.

      implicit double precision(a-h, o-z)
      include 'jam2.inc'

c...Commonblock for t-channel resonance productions.
      common/jamres1/kfo(2,20),noutpa

      dimension sigin(mxchan)
      parameter(srt0=2.0d0)

      izt=izb+izm
      idb=kchg(kcb,5)
      idm=kchg(kcm,5)

c....Get total cross section.
c==================================================================
      if(msel.eq.1) then

c.....pion-N->Delta(1232)
        if( jamcpair(idb,idm) .eq. jamcpair(id_nucl,id_pi) ) then
          if(srt.le.1.85d0) then
            call jamxbw1(istr,srt,pr,kfm,kfb,izm,izb,sig,sige,msel)
          else
            call jamxpin(izm,izb,srt,pr,sig,sigel)
          endif
c....t-channel elastic.
          if(srt.le.4.0d0) then
            if(srt.le.1.4d0) then
              sigel=0.0d0
            else
              if(izt.eq.2.or.izt.eq.-1) then
                sigel=3.81d0*(srt-1.4d0)**2.156d0
     $                         /((srt-1.66d0)**2+0.0352d0)
              else
                sigel=45.49d0*(srt-1.4d0)**0.854d0
     $                           /((srt+1.22d0)**2-5.585d0)
              endif
            endif
          endif
          if(srt.le.1.85d0) sig=sig+sigel
          return
        endif

c....Higher resonance formation.
        call jamxbw1(istr,srt,pr,kfm,kfb,izm,izb,sig,sige,msel)
        if(srt.ge.2.0d0) then
          call jamxadq(kfm,kfb,siga,sigel)
          sig=max(sig,siga)
          sig=min(sig,1000.d0)
        endif

        if(izt.lt.-1.or.izt.gt.2) then
        if(srt.lt.parc(62)) then
          sig=sigel
        endif
        endif

        return 

c...End in case of total cross section.
      endif
c==================================================================

      sigin(1)=0.0d0  ! t-channel inel.
      sigin(2)=0.0d0  ! s-channel resonance formation
      sigin(3)=0.0d0  ! s-channel string formation
      sigabc=0.0d0
      sigres=0.0d0
      mchanel=1
      mabsrb=2

      if( jamcpair(idb,idm) .eq. jamcpair(id_nucl,id_pi) ) then
        call jamxpin(izm,izb,srt,pr,sig,sigel)
      else
        call jamxadq(kfm,kfb,sig,sigel)
      endif

c...Absorption impossible.
      if(izt.lt.-1.or.izt.gt.2) goto 2000

c.....Find resonance formation cross section.
      if(srt.lt.5.0d0) then
        call jamxbw1(istr,srt,pr,kfm,kfb,izm,izb,sigres,sige,msel)

c       if(srt.le.2.0d0.and.pare(3).gt.0.0d0) then
c        write(6,*)'funny! srt kfm kfb',srt,kfm,kfb,sigres,pare(3)
c        stop
c       endif

c....Monte Carlo for resonance formation.
        if(msel.eq.3.and.pare(3).le.0.0d0) then
          kf1=kfm
          kf2=0
          em1=srt
          em2=0
          return
        endif
        sigin(2)=sigres
      endif

      sig=max(sig,sigres)
      sig1=max(0.0d0,sig-sigel-sigres)

c.....Find s-channel contribution to the continuum.
      if(srt.le.srt0) then
c       sigabc=sig1
        sigabc=0.0d0
      else
        sigabc=sig1*srt0/srt
      endif
 
c....Monte Carlo for s-channel string formation.
        if(msel.eq.3) then
          pare(3)=pare(3)-sigabc
          if(pare(3).le.0.0d0) then
            if(izt.eq.-1) kf1=1114
            if(izt.eq.0) kf1=12112
            if(izt.eq.1) kf1=12212
            if(izt.eq.2) kf1=2224
            kf2=0
            em1=srt
            em2=0.0d0
            ijet=2
               if(srt.lt.2.0d0) then
                  write(6,*)'(jamcmbs0:) sigabc kf1 kf2 srt',kf1,kf2,srt
                  stop
                endif
            return
          endif

        else
          mabsrb=2
          sigin(3)=sigabc
        endif

2000  continue
      sigin(1)=max(0.0d0,sig-sigel-sigabc-sigres)
      if(msel.eq.3) then
         pare(3)=-10.d0
c        if(srt.ge.parc(62)) then
         if(srt.ge.3.0d0) then
           ijet=1
         else
           noutpa=0
           call jamrmas2(kfm,kfb,kcm,kcb,srt,emmes,embar,icon)
           if(icon.ne.0) ijet=-1
           if(jswap.eq.0) then
             kf1=kfm
             kf2=kfb
             em1=emmes
             em2=embar
           else
             kf1=kfb
             kf2=kfm
             em1=embar
             em2=emmes
           endif
c          if(izt.ge.-1.and.izt.le.2)
c    $         write(6,*)'(jamcmbs0): srt kf1 kf2',srt,kfm,kfb,izt
c          ijet=-1

         endif
         return
      endif

      end

c***********************************************************************

      subroutine jamcmbs1(msel,srt,pr,istr,kfb,kfm,kcb,kcm,izb,izm,
     $ emmes,embar,kf1,kf2,em1,em2,sig,sigel,sigin,
     $ mabsrb,mchanel,mxchan,ijet,jswap,icon)

c...Strangeness S=-1 meson-baryon(K-p,..) collisions.

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
c...Commonblock for t-channel resonance productions.
      common/jamres1/kfo(2,20),noutpa

c...Local arrays.
      dimension sigin(mxchan),sigy(4),sigyv(4)
      parameter(srt0=1.78d0)
      dimension kfy(4,2,4),kfkp(4,2)

c...Functions: momentum in two-particle cm.
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a) 

c...k-n
      data ((kfy(1,i,j),i=1,2),j=1,4)/
     $   -211,3122, 111,3112, -211,3212, 0,0/
c...k0p
      data ((kfy(2,i,j),i=1,2),j=1,4)/
     $    211,3122, 211,3212, 111,3222, 0,0/
c...k-p
      data ((kfy(3,i,j),i=1,2),j=1,4)/
     $     111,3122, 211,3112, 111,3212, -211,3222/
c...k0n
      data ((kfy(4,i,j),i=1,2),j=1,4)/
     $    111,3122, 211,3112, 111,3212, -211,3222/

      data (kfkp(1,i),i=1,2)/-321,2112/
      data (kfkp(2,i),i=1,2)/-311,2212/
      data (kfkp(3,i),i=1,2)/-321,2212/
      data (kfkp(4,i),i=1,2)/-311,2112/

      kf1=kfb
      kf2=kfm
      sigin(1)=0.0d0  ! t-channel inel.
      sigin(2)=0.0d0  ! s-channel resonance formation
      sigin(3)=0.0d0  ! s-channel string formation
      sigabc=0.0d0
      sigres=0.0d0
      mchanel=1
      mabsrb=2
      izt=izb+izm
      ik=0
c......K- n/antiK0p
      if(kfm.eq.-321.and.kfb.eq.2112) ik=1
      if(kfm.eq.-311.and.kfb.eq.2212) ik=2
      if(kfm.eq.-321.and.kfb.eq.2212) ik=3
      if(kfm.eq.-311.and.kfb.eq.2112) ik=4
      if(izt.eq.1) then
       ipk1=1
       ipk2=2
      else if(izt.eq.-1) then
       ipk1=2
       ipk2=2
      else if(izt.eq.0) then
       ipk1=3
       ipk2=4
      else
       ipk1=0
       ipk2=0
      endif

c...Initialize cross sections.
      noutpa=0
      siga=0.0d0
      sigela=0.0d0
      sigch=0.0d0
      do i=1,4
       sigy(i)=0.0d0
      end do

c....Get t-channel elastic, charge exchange and piY cross section.
      if(ik.ge.1) then
        call jamxkp(kfm,kfb,srt,emmes,embar,siga,sigela,sigch,sigy)

c...Save outgoing particle types.
        do i=1,4
          noutpa=noutpa+1
          kfo(1,noutpa)=kfy(ik,1,i)
          kfo(2,noutpa)=kfy(ik,2,i)
        enddo
      endif

c...Find resonance formation cross sections.
      sigr=0d0
      if(srt.le.5.0d0) then
        call jamxbw1(istr,srt,pr,kfm,kfb,izm,izb,sigr,sige,msel)
        if(msel.eq.3.and.pare(3).le.0.0d0) then
          kf1=kfm
          kf2=0
          em1=srt
          em2=0
          return
        endif
      endif

      sigin(2)=sigr
      sigin(1)=sigch+sigy(1)+sigy(2)+sigy(3)+sigy(4)
      sig=sigela+sigin(1)+sigin(2)
      sigel=sigela                ! background t-channel elastic.

c...t-channel piY->antKN cross section from detailed balance.
      sigyp=0.0d0
      if(ipk1.ne.0) then

c.....Find wether this is pi Y ingoing channel.
          do i=1,ipk1,ipk2
            do j=1,4
            kfi1=kfy(i,1,j)
            kfi2=kfy(i,2,j)
            if(kfi1.eq.kfm.and.kfi2.eq.kfb) then
              ikf=i
              icg=j
              goto 200
            endif
            end do
          end do
          goto 220
 200      continue
c.....Outgoing particle codes and masses.
          kfi1=kfkp(ikf,1)
          kfi2=kfkp(ikf,2)
          emi1=pmas(jamcomp(kfi1),1)
          emi2=pmas(jamcomp(kfi2),1)
          if(srt.lt.emi1+emi2+0.0001d0) goto 220
          prf=pawt(srt,emi1,emi2)
          call jamxkp(kfi1,kfi2,srt,emi1,emi2,tmp1,tmp2,tmp3,sigyv)
          spini=max(1,mod(abs(kfm),10))*max(1,mod(abs(kfb),10))
          sigyp=sigyv(icg)*(prf/pr)**2*2.0d0/spini
        endif
 220  continue
      sigin(1)=sigin(1)+sigyp

c...Additive quark cross section except KN ingoing above resonance
c...region.
c     if(ik.eq.0.and.srt.ge.srt0) then
      if(ik.eq.0) then
        call jamxadq(kfm,kfb,siga,sigela)
      endif

c...Calculate total and elastic cross sections.
      if(msel.eq.1) then
        if(ik.ge.1) then
         if(srt.le.1.8d0) then
           sig=sigr+sigela+sigch+sigy(1)+sigy(2)+sigy(3)+sigy(4)
         else
           sig=siga
         endif
        else
          sig=max(sig+sigyp,siga)
          sigel=sigela
        endif
        sig=min(sig,300.d0)
        return
      endif
c=================================================================

c...Monte Carlo for t-channel back ground reactions. Now KN only.
      if(msel.eq.3.and.ik.ne.0) then

c...First charge exchange reaction.

        pare(3)=pare(3)-sigch
        if(pare(3).le.0.0d0) then
          if(izm.eq.-1.and.izb.eq.1) then
            kfm=-311
            kfb=2112
          else if(izm.eq.0.and.izb.eq.0) then
            kfm=-321
            kfb=2212
          endif
          emmes=pmas(jamcomp(kfm),1)
          embar=pmas(jamcomp(kfb),1)
          if(srt.le.em1+em1+0.001d0) then
            pare(3)=pare(3)+sigch
            goto 10
          endif
          goto 4000
        endif

c...Save outgoing particle types for charge exchange reaction.
        noutpa=noutpa+1
        kfo(1,noutpa)=kfm
        kfo(2,noutpa)=kfb

c....Monte-Calro for KN->piY.
 10     do i=1,4
          pare(3)=pare(3)-sigy(i)
          if(pare(3).le.0.0d0) then
            kfm=kfy(ik,1,i)
            kfb=kfy(ik,2,i)
            emmes=pmas(jamcomp(kfm),1)
            embar=pmas(jamcomp(kfb),1)
            goto 4000
          endif
        end do
      endif

c....Monte-Calro for piY->KN.
      if(msel.eq.3) then
        pare(3)=pare(3)-sigyp
        if(pare(3).le.0.0d0) then
          kfm=kfi1
          kfb=kfi2
          emmes=emi1
          embar=emi2
          noutpa=noutpa+1
          kfo(1,noutpa)=kfi1
          kfo(2,noutpa)=kfi2
          goto 4000
        endif
      endif

c...Gap of cross section from experimental total xsection.
      sig1=max(0.0d0,siga-sig)

c...t-channel resonance (string impossible due to low energy).
      if(srt.ge.srt0.and.srt.le.2.5d0) then
        sigtc=sig1
      else
        sigtc=sig1*2.5d0/srt
      endif
      sigin(1)=sigin(1)+sigtc
      sig1=max(0.0d0,sig1-sigtc)
      if(msel.eq.3) then
        pare(3)=pare(3)-sigtc
        if(pare(3).le.0.0d0)  goto 3000
      endif

c...Absorption impossible.
      if(izt.lt.-1.or.izt.gt.1) goto 2000

c.....Find s-channel contribution to the continuum.
      if(srt.le.3.0d0) then
        sigabc=sig1
      else
        sigabc=sig1*3.0d0/srt
      endif
      sig1=max(0.0d0,sig1-sigabc)

c....Monte Carlo for s-channel string formation.
      if(msel.eq.3) then
          pare(3)=pare(3)-sigabc
          if(pare(3).le.0.0d0) then
            if(izt.eq.-1) kf1=3112
            if(izt.eq.0) kf1=3212
            if(izt.eq.1) kf1=3222
            kf2=0
            em1=srt
            em2=0.0d0
            ijet=2
            return
          endif
      else
          sigin(3)=sigabc
      endif

c...t-channel string formation.
 2000 continue
      sigin(1)=sigin(1)+sig1
      if(msel.eq.3) then
         pare(3)=-10.d0
         ijet=1
         goto 4000
      endif
      return

c...t-channel resonance formation.
 3000 continue
      call jamrmas2(kfm,kfb,kcm,kcb,srt,emmes,embar,icon)
      if(icon.ne.0) ijet=-1

 4000 if(jswap.eq.0) then
        kf1=kfm
        kf2=kfb
        em1=emmes
        em2=embar
      else
        kf1=kfb
        kf2=kfm
        em1=embar
        em2=emmes
      endif

      end

c***********************************************************************

      subroutine jamcmbs2(msel,srt,pr,istr,kfb,kfm,kcb,kcm,izb,izm,
     $ em1,em2,kf1,kf2,sig,sigel,sigin,
     $ mabsrb,mchanel,mxchan,ijet,icon)

c...Strangeness S=-2 meson-baryon cross section. xi-pi,...

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      dimension sigin(mxchan)

      izt=izb+izm
c....Get total cross section.
      if(msel.eq.1) then
            call jamxbw1(istr,srt,pr,kfm,kfb,izm,izb,sig,sige,msel)
            call jamxadq(kf1,kf2,siga,sigel)
            if(srt.lt.parc(62).and.(izt.lt.-1.or.izt.gt.1)) then
              sig=sigel
            else
              sigel=0.0d0
            endif
            return 
      endif

c.....Find resonance formation.
          if(srt.lt.3.0d0) then
           call jamxbw1(istr,srt,pr,kfm,kfb,izm,izb,sigres,sige,msel)
              if(msel.eq.3) then
              if(pare(3).le.0.0d0) then
                 kf1=kfm
                 kf2=0
                 em1=srt
                 em2=0
                 return
              endif
              endif

          endif

          if(msel.eq.3) then
                pare(3)=-10.d0
                if(srt.ge.parc(62)) then
                 ijet=1
                else
                 ijet=-1
                endif
              return
          endif

              mabsrb=1
              sigin(2)=sigres
              icon=2

      end

c***********************************************************************

      subroutine jamckaon(msel,kfm,kfb,srt,pr,ijet
     $        ,em1,em2,kf1,kf2,sig,sigel,sigin,mchanel,mxchan,ich,icon)

c...Purpose: to treat Kaon-nonstrange baryon collisions.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

c...Commonblock for t-channel resonance productions.
      common/jamres1/kfo(2,20),noutpa

c...msel=1: calculate total and elastic cross sections.
c...msel=2: calculate inelastic cross sections.
c...msel=3: Monte Calro evaluation of inelastic channel.
      real*8 jamsigkn
      real*8 jamchc96,jamrgg96
      parameter(emkc=.49360d0, emk0=.49770d0,emk=.495650d0,widk=0.05d0)
      parameter(empion=0.139d0,emdelt=1.232d0,widdlt=0.12d0)
      parameter(emp=.93830d0 ,emn=.93960d0,emnuc=0.93895d0)
      parameter(ekinmi=0.001d0,emdmin=1.08d0)
      parameter (bhad=2.3d0,eps=0.0808d0,facel=0.0511d0)
      dimension sigin(mxchan)
      dimension kfdelt(2,2,4)
c...Functions: momentum in two-particle cm.
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a) 
c...Lab. momentum.
      plabsr(a,b,c)=sqrt((a**2-b**2-c**2)**2/(4.d0*c**2)-b**2)

c....Data for one delta production.
      data kfdelt/
     $  321,2214, ! k+ p->k+ d+ 
     $  311,2224, ! k+ p->k0 d++
     $  311,2114, ! k0 n->k0 d0
     $  321,1114, ! k0 n->k+ d-
     $  321,2114, ! k+ n->k+ d0
     $  311,2214, ! k+ n->k0 d+
     $  311,2214, ! k0 p->k0 d+
     $  321,2114/ ! k0 p->k+ d0

c...Check particle code.
      if(kfm.le.0.or.kfb.le.0) then
         write(check(1),8000)kfm,kfb,msel,srt
         call jamerrm(1,1,'(jamckaon:)invalid particle')
         icon=999
      endif

      mchanel=0
      ijet=0
      xsig0=pare(3)

      kfma=abs(kfm)
      kflr1=mod(kfma/10000,10) 
      kflb1=mod(kfma/100,10) 
      kflc1=mod(kfma/10,10) 
      kfls1=mod(kfma,10) 
      kfm1=10*kflb1+kflc1

      kcb=jamcomp(kfb)
      kcm=jamcomp(kfm)
      izm=kchg(kcm,1)*isign(1,kfm)/3
      izb=kchg(kcb,1)*isign(1,kfb)/3
      id2=kchg(kcb,5)
      izt=izm+izb

      ibra=0
      if(izm.eq.1.and.izb.eq.1)  ibra=1  ! k+ p
      if(izm.eq.0.and.izb.eq.0)  ibra=2  ! k0 n
      if(izm.eq.1.and.izb.eq.0)  ibra=3  ! k+ n
      if(izm.eq.0.and.izb.eq.1)  ibra=4  ! k0 p

      if(izm.eq.1.and.izb.eq.2)  ibra=5  ! k+ d++
      if(izm.eq.0.and.izb.eq.-1) ibra=6  ! k0 d-
      if(izm.eq.1.and.izb.eq.-1) ibra=7  ! k+ d-
      if(izm.eq.0.and.izb.eq.2)  ibra=8  ! k0 d++

      if((kfm.eq.311.or.kfm.eq.321)
     $           .and. (kfb.eq.2112.or.kfb.eq.2212)) then
         if(kfm.eq.311) ema=emk0
         if(kfm.eq.321) ema=emkc
         if(kfb.eq.2112) emb=emn
         if(kfb.eq.2212) emb=emp
         snew=srt
         if(srt.gt.ema+emb) then
           plab=sqrt( ((srt*srt-ema*ema-emb*emb)/(2.d0*emb) )**2
     $              -ema*ema)
         else
          plab=0.0d0
         endif
      else
        if(kfm1.eq.32) then
          snew=sqrt(emnuc**2+pr**2)+sqrt(emkc**2+pr**2)
          if(snew.gt.emkc+emnuc) then
          plab=plabsr(snew,emkc,emnuc)
c         plab=sqrt(((snew*snew-emkc*emkc-emnuc*emnuc)
c    $           /(2.d0*emnuc))**2-emkc*emkc)
          else
            write(check(1),8100)srt,kfm,em1,kfb,em2
            call jamerrm(1,1,'(jamckaon:1) snew<emkc+emp')
            plab=0.0d0
            sig=0d0
            sigel=0d0
            return
          endif

        else if(kfm1.eq.31) then
          snew=sqrt(emnuc**2+pr**2)+sqrt(emkc**2+pr**2)
          if(snew.gt.emk0+emnuc) then
            plab=plabsr(snew,emk0,emnuc)
          else
            write(check(1),8100)srt,kfm,em1,kfb,em2
            call jamerrm(1,1,'(jamckaon:2) snew<emkc+emp')
            plab=0.0d0
            sig=0d0
            sigel=0d0
            return
          endif
        else
          write(check(1),8100)srt,kfm,em1,kfb,em2
          call jamerrm(1,1,'(jamckaon:3) no path')
          plab=0.0d0
          sig=0d0
          sigel=0d0
          return
        endif
      endif

c------------------------------------------------
c.....Get total and elastic cross sections.
c------------------------------------------------

      if(msel.eq.1) then

c......K+ p/K0n i.e. isospin =1 channel.
          if(izm.eq.izb) then

c...... K+ p Total/elastic
            if (plab.le.3.0d0) then
              call jamsighh(sig,15,snew)
              call jamsighh(sigel,16,snew)
            else if(snew.lt.30.d0) then
              sig=jamchc96(8,plab)
              sigel=jamchc96(9,plab)
            else
              sig=jamrgg96(snew,10)
              bel=2*(2.3d0+0.8d0)+4.d0*(srt*srt)**0.079d0-4.2d0
              sigel=facel*sig**2/bel
            endif

c.......K+n/K0p
          else
 
c... K+ n Total/elastic
            if(plab.le.4.0d0)then
              if(plab.le.1.05d0) then
                sig=18.6d0*plab**(0.86d0)
              else
                call jamsighh(sig,17,snew)
              endif
              if(plab.le.1.05d0) then
                sigel=7.6761d0*plab**(0.566017d0)
              else if(plab.le.4.5d0)then
                sigel=4.03297d0+9.71556d0
     $           /(1.10417d0+exp((plab-0.674619d0)/0.835265d0))
              else
                call jamsighh(sigel,18,snew)
              endif
            else if(snew.lt.20.d0) then
              sig=jamchc96(10,plab)
              sigel=jamchc96(9,plab)
            else
              sig=jamrgg96(snew,11)
              bel=2*(2.3d0+0.8d0)+4.d0*(srt*srt)**0.079d0-4.2d0
              sigel=facel*sig**2/bel
            endif

          endif

         return
      endif

c------------------------------------------------
c.....Get charge exchange cross sections.
c------------------------------------------------
      ich=0

      if(izt.ne.1) goto 10
      if(srt.lt.emn+emkc+ekinmi) goto 10

c...Charge exchange.
      mchanel=1
      snew1=sqrt(emn**2+pr**2)+sqrt(emkc**2+pr**2)
      if(snew.gt.emkc+emn) then
      plab1=sqrt( ((snew*snew-emkc*emkc-emn*emn)/(2.d0*emp) )**2
     $              -emkc*emkc)
      else
        goto 10
      endif

      if(plab1.le.2.0d0) then
        call jamsighh(sigin(1),20,snew1)
      else
       pl=min(6.3d0,plab1)
       sigin(1)=jamchc96(26,pl)
      endif

c...Monte Calro for charge exchange.
      if(msel.eq.3) then
        ich=1
        pare(3)=pare(3)-sigin(1)
        if(pare(3).le.0.0d0) then
          if(kfm.eq.321) then
            kf1=311
            em1=emk0
          else if(kfm.eq.311) then
            kf1=321 
            em1=emkc
          else if(kfm1.eq.31) then
            kf1=kflr1*10000+320+kfls1
          else if(kfm1.eq.32) then
            kf1=kflr1*10000+310+kfls1
          endif 

          if(kfb.eq.2112) then
            kf2=2212
            em2=emp
          else if(kfb.eq.2212) then
            kf2=2112
            em2=emn
          else
            kc2=jamcomp(kfb)
            if(izb.eq.1) then
              kf2=kchg(kc2-1,4)
            else if(izb.eq.0) then
              kf2=kchg(kc2+1,4)
            endif
          endif
          return
        endif
      endif

c------------------------------------------------
c.....Get inelastic cross sections.
c------------------------------------------------
c...(1) K^+ + p --> K^+  + D+   1/4*sig(t=1)
c...(2) K^+ + p --> K^0  + D++  3/4*sig(t=1)
c...(3) K^+ + p --> K^+* + p (k+* p --> k+ p k* absorption)

c...1. K^+ + n --> K^0 + p                          charge exchange
c...2. K^+ + n --> K^+ + D0  0.25*sig(t=1) isig=1   delta production
c...3. K^+ + n --> K^0 + D+  0.25*sig(t=1) isig=1   delta production
c...4. K^+ + n --> K^+* + n  isig=3                 k* production
c...5. K^+ + n --> K^0* + p  isig=4                 k* production

c...1. K^0 + p --> K^+ + n    charge exchange
c...2. K^0 + p --> K^0 + D+   0.25*sig(t=1)
c...3. K^0 + p --> K^+ + D0   0.25*sig(t=1)
c...4. K^0 + p --> K^0* + p   3
c...5. K^0 + p --> K^+* + n   4

c...1. K^0 + n --> K^0 + D0   0.25
c...2. K^0 + n --> K^+ + D-   0.75
c...3. K^0 + n --> K^0* + n   t=1  equated to k+p->k+* p
c     (K^0* + n --> K^0 + n   k* absorption )

10    continue

c...Check if there is enough energy to produce delta
      if(srt.le.(em1+emdmin+empion+ekinmi)) then
        goto  150
      endif

      if(id2.eq.id_delt.or.id2.eq.id_delts) goto 100

c...(2) D(1232) produciton.
      mchanel=mchanel+1
      sigin(mchanel)=jamsigkn(1,snew,plab)

c....Monte Calro evaluation of delta channel.
      if(msel.eq.3) then

          ich=2
          if(ibra.eq.1.or.ibra.eq.2) then
           sigin1=sigin(mchanel)*0.25d0
           sigin2=sigin(mchanel)*0.75d0
          else
           sigin1=sigin(mchanel)*0.5d0
           sigin2=sigin(mchanel)*0.5d0
          endif

          pare(3)=pare(3)-sigin1
          iselect=0
          if(pare(3).le.0.0d0) then
            iselect=1
          else
            pare(3)=pare(3)-sigin2
            if(pare(3).le.0.0d0) iselect=2
          endif

          if(iselect.ne.0) then
            kf1=kfdelt(1,iselect,ibra)
            kf2=kfdelt(2,iselect,ibra)
            kc1=jamcomp(kf1)
            kc2=jamcomp(kf2)
            em1=pmas(kc1,1)

            emmin=max(pmas(kc2,1)-pmas(kc2,3),emdmin+parc(41))
            emmax=srt-em1-ekinmi
            emr=pmas(kc2,1)
            wid=pmas(kc2,2)
            call jambwmas(emmin,emmax,emr,wid,em2,icon)
            return
          endif
     
      endif

100   continue
c...(3) D(1232) absorption using detaild balance.

      if(id2.eq.id_delt.and.(izt.ge.0.and.izt.le.2)) then

        widcof=max(1.0d0,paru(1)/(atan(2.d0*(srt-emk-emdelt)/widdlt)-
     a         atan(2.d0*(emnuc+empion-emdelt)/widdlt)))
        if(srt.gt.em1+emnuc) then
            pfinal=pawt(srt,em1,emnuc)
            mchanel=mchanel+1
            sigin(mchanel)=0.5d0*(pfinal/pr)**2
     $                     *jamsigkn(1,srt,plab)*widcof
        else
          write(check(1),'(''kfm kfb srt'',i9,1x,i9,1x,g10.3)')
     $       kfm,kfb,srt
          call jamerrm(3,1,
     $       '(jamckaon:) srt<emk+emnuc D(1232)absroption')
          pfinal=0.0d0
        endif

c....Monte Calro evaluation of delta absorption channel.
      if(msel.eq.3) then
          ich=3
          pare(3)=pare(3)-sigin(mchanel)
          iselect=0
          if(pare(3).le.0.0d0) then
            if(kfb.eq.1114) then ! d-
              kf2=2112
              kf1=10000*kflr1+310+kfls1
            else if(kfb.eq.2114) then ! d0
             if(izm.eq.1) then        ! d0 k+
               if(rn(0).le.0.5d0) then
                 kf2=2112
                 kf1=10000*kflr1+320+kfls1
               else
                 kf2=2212
                 kf1=10000*kflr1+310+kfls1
               endif
             else if(izm.eq.0) then  ! d0 k0
                kf2=2112
                kf1=kfm
             endif
            else if(kfb.eq.2214) then ! d+
             if(izm.eq.0) then       ! d+ k0
               if(rn(0).le.0.5d0) then
                 kf2=2112
                 kf1=10000*kflr1+320+kfls1
               else
                 kf2=2212
                 kf1=10000*kflr1+310+kfls1
               endif
             else if(izm.eq.1) then  ! d+ k+
                kf2=2212
                kf1=kfm
             endif
            else if(kfb.eq.2224) then ! d++
               kf2=2212
               kf1=10000*kflr1+320+kfls1
            endif

            if(kfm.eq.311.or.kfm.eq.321) em1=pmas(jamcomp(kf1),1)
            em2=pmas(jamcomp(kf2),1)
            return
          endif
     
      endif
   
      endif

150   continue
c...Check if there is enough energy to produce delta
      if(srt.le.(em2+empion+emk+empion+ekinmi)) then
        goto  5000
      endif

c....(4) K(892) produciton.
      mchanel=mchanel+1
      if(ibra.le.2.or.ibra.eq.5.or.ibra.eq.6) then
        sigin1=jamsigkn(2,snew,plab)
        sigin(mchanel)=sigin1
        jch=1
      else
        jch=2
        sigin1=jamsigkn(3,snew,plab)
        sigin(mchanel)=sigin1
        mchanel=mchanel+1
        sigin2=jamsigkn(4,snew,plab)
        sigin(mchanel)=sigin2
      endif

c....Monte Calro evaluation of K(892) production channel.
      if(msel.eq.3) then
        ich=4
        pare(3)=pare(3)-sigin1
        iselect=0
        if(pare(3).le.0.0d0) then
          iselect=1
        else if(jch.eq.2) then
          pare(3)=pare(3)-sigin2
          if(pare(3).lt.0.0d0) iselect=2
        endif

        if(iselect.ne.0) then
            if(iselect.eq.1) then
             if(izm.eq.0) kf1=313
             if(izm.eq.1) kf1=323
             kf2=kfb
            else
             if(izm.eq.0) then
                kf1=323
                kf2=kchg(kcb-1,4)
             else if(izm.eq.1) then
                kf1=313
                kf2=kchg(kcb+1,4)
             endif
            endif
            kc1=jamcomp(kf1)
            emmin=pmas(kc1,1)-pmas(kc1,3)
            emmax=srt-em2-ekinmi
            emr=pmas(kc1,1)
            wid=pmas(kc1,2)
            call jambwmas(emmin,emmax,emr,wid,em1,icon)
            return
        endif
     
      endif

      if(kfm.eq.311.or.kfm.eq.321) goto 200

c...(5) K*(892) absorption using detaild balance.

        widcof=max(1.0d0,paru(1)/(atan(2.d0*(srt-emnuc-0.892d0)/widk)-
     a         atan(2.d0*(emk+empion-0.892d0)/widk)))
        if(srt.gt.emk+emnuc) then
          pfinal=pawt(srt,emk,emnuc)
        else
        write(check(1),'(''kfm kfb srt'',i9,1x,i9,1x,g10.3)')kfm,kfb,srt
        call jamerrm(30,1,'(jamckaon:) srt<emk+emnuc K*(892)abs.')
        endif

      mchanel=mchanel+1
      if(ibra.le.2.or.ibra.eq.5.or.ibra.eq.6) then
        sigin1=jamsigkn(2,srt,plab)*(pfinal/pr)**2*widcof
        sigin(mchanel)=sigin1
        jch=1
      else
        jch=2
        sigin1=jamsigkn(3,srt,plab)*(pfinal/pr)**2*widcof
        sigin(mchanel)=sigin1
        sigin2=jamsigkn(4,srt,plab)*(pfinal/pr)**2*widcof
        mchanel=mchanel+1
        sigin(mchanel)=sigin2
      endif

c....Monte Calro evaluation of K(892) absorption channel.
      if(msel.eq.3) then
        ich=5
        pare(3)=pare(3)-sigin1
        iselect=0
        if(pare(3).le.0.0d0) then
          iselect=1
        else if(jch.eq.2) then
          pare(3)=pare(3)-sigin2
          if(pare(3).le.0.0d0) iselect=2
        endif

        if(iselect.ne.0) then
            if(iselect.eq.1) then
             if(izm.eq.0) kf1=311
             if(izm.eq.1) kf1=321
             kf2=kfb
            else
             if(izm.eq.0) then
                kf1=321
                kf2=kchg(kcb-1,4)
             else if(izm.eq.1) then
                kf1=311
                kf2=kchg(kcb+1,4)
             endif
            endif
            if(kf1.eq.311) em1=emk0
            if(kf1.eq.321) em1=emkc
            return
        endif
     
      endif




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
200   continue

      if(id2.eq.id_delt.or.id2.eq.id_delts) goto 300
      if(srt.le.(emk+2*empion+emnuc+ekinmi)) goto 300

c...(6) K(892)+D(1232) produciton.

      mchanel=mchanel+1
      sigin(mchanel)=jamsigkn(5,snew,plab)

      if(msel.eq.3) then
          ich=6
          if(ibra.eq.1.or.ibra.eq.2) then
           sigin1=sigin(mchanel)*0.25d0
           sigin2=sigin(mchanel)*0.75d0
          else
           sigin1=sigin(mchanel)*0.5d0
           sigin2=sigin(mchanel)*0.5d0
          endif

          pare(3)=pare(3)-sigin1
          iselect=0
          if(pare(3).le.0.0d0) then
            iselect=1
          else
            pare(3)=pare(3)-sigin2
            if(pare(3).le.0.0d0) iselect=2
          endif

          if(iselect.ne.0) then
            kf1=kfdelt(1,iselect,ibra)+2  ! K*=312/323
            kf2=kfdelt(2,iselect,ibra)
            kc1=jamcomp(kf1)
            kc2=jamcomp(kf2)

            emmin1=pmas(kc1,1)-pmas(kc1,3)
            emmin2=max(pmas(kc2,1)-pmas(kc2,3),emdmin+parc(41))
            if(srt.le.emmin1+emmin2+ekinmi) then
              pare(3)=xsig0
              kf1=kfm
              kf2=kfb
              goto 300
            endif

            emmax1=srt-emmin2-ekinmi
            emr1=pmas(kc1,1)
            wid1=pmas(kc1,2)
            itry=0
210         continue
            itry=itry+1
            call jambwmas(emmin1,emmax1,emr1,wid1,em1,icon)

            emmax2=srt-em1-ekinmi
            if(emmax2.le.emmin2) then
               if(itry.le.50) goto 210
               em1=emmin1
               em2=emmin2
            else
              emr2=pmas(kc2,1)
              wid2=pmas(kc2,2)
              call jambwmas(emmin2,emmax2,emr2,wid2,em2,icon)
            endif
            return
          endif
     
      endif


300   continue

      if(kfm.eq.311.or.kfm.eq.321) goto 5000 
      if(id2.ne.id_delt)  goto 5000
      if(izt.lt.0.or.izt.gt.2)  goto 5000

c...(7) K*(892)D(1232) absorption using detaild balance.
      widcof1=max(1.0d0,paru(1)/(atan(2.d0*(srt-emnuc-0.892d0)/widk)-
     a         atan(2.d0*(emk+empion-0.892d0)/widk)))

      widcof2=max(1.0d0,paru(1)/(atan(2.d0*(srt-emk-1.232d0)/widdlt)-
     a         atan(2.d0*(emnuc+empion-1.232d0)/widdlt)))

      if(srt.gt.emk+emnuc) then
          pfinal=pawt(srt,emk,emnuc)
      else
        write(check(1),'(''kfm kfb srt'',i9,1x,i9,1x,g10.3)')kfm,kfb,srt
        call jamerrm(30,1,'(jamckaon:) srt<emk+emnuc at K*D abs.')
      endif

      mchanel=mchanel+1
      sigin(mchanel)=0.5d0*(pfinal/pr)**2
     $ *jamsigkn(5,srt,plab)*widcof1*widcof2
      if(msel.eq.3) then
          ich=7
          pare(3)=pare(3)-sigin(mchanel)
          iselect=0
          if(pare(3).le.0.0d0) then
            if(kfb.eq.1114) then ! d-
              kf1=311
              kf2=2112
            else if(kfb.eq.2114) then ! d0
             if(izm.eq.1) then        ! d0 k+
               if(rn(0).le.0.5d0) then
                 kf2=2112
                 kf1=321
               else
                 kf2=2212
                 kf1=311
               endif
             else if(izm.eq.0) then  ! d0 k0
                kf2=2112
                kf1=311
             endif
            else if(kfb.eq.2214) then ! d+
             if(izm.eq.0) then       ! d+ k0
               if(rn(0).le.0.5d0) then
                 kf2=2112
                 kf1=321
               else
                 kf2=2212
                 kf1=311
               endif
             else if(izm.eq.1) then  ! d+ k+
                kf2=2212
                kf1=321
             endif
            else if(kfb.eq.2224) then ! d++
               kf2=2212
               kf1=321
            endif

            if(kf1.eq.311) em1=emk0
            if(kf1.eq.321) em1=emkc
            if(kf2.eq.2112) em2=emn
            if(kf2.eq.2212) em2=emp
            return
          endif
     
      endif

c------------------------------------------------------------
5000  continue
      if(msel.eq.3.and.pare(3).gt.0.0d0) then
         pare(3)=-120.d0
         if(srt.ge.4.0d0) then
            ijet=1
         else
           noutpa=0
           call jamrmas2(kfm,kfb,kcm,kcb,srt,em1,em2,icon)
           if(icon.ne.0) then
             ijet=-1
             return
           endif
           kf1=kfm
           kf2=kfb
         endif
      endif

 8000 format('kfm kfb msel srt=',i9,1x,i9,1x,i2,g10.3)
 8100 format('srt',g10.3,' kfm em1',i9,1x,g11.4,' kfb em2',i9,1x,g11.4)

      end

c***********************************************************************

      subroutine jamcpipi(msel,srt,kf1,kf2,iz1,iz2,
     $     em1,em2,sig,sigel,sigab,mchanel,mabsrb,ijet)

c...Purpose: to treat low energy(rho region) pi-pi cross sections.

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      parameter(emrhot=0.770d0,empi=0.138d0,emrh=5.8d0*empi)
      parameter(srtpi=1.3d0)

      iztot=iz1+iz2
c...Calculate pi-pi ispspin cross sections.
c...Ref. G.Bertsch, P.R.D37 (1988) 1202
      sig0=0.0d0
      sig1=0.0d0
      sig2=0.0d0
c...Check energy.
      if(srt.le.em1+em2+0.001d0) then
        write(check(1),8000)srt,em1,em2
 8000   format('srt',g12.3,'em1 em2',2(g12.3,1x))
        call jamerrm(3,1,'(jamcpipi:) Energy too small')
        goto 100
      endif
      prel=sqrt((srt**2-(em1+em2)**2)*(srt**2-(em1-em2)**2))/(2.d0*srt)
      gamrh=0.095d0*prel*(prel/empi/(1.0d0+(prel/emrhot)**2))**2
      if(abs(emrh-srt).gt.0.00001d0) then
        del00=atan(0.5d0*2.06d0*prel/(emrh-srt))
      else
        del00=1.57d0
      endif
      if(abs(emrhot-srt).gt.0.00001d0) then
        del11=atan(0.5d0*gamrh/(emrhot-srt))
      else
        del11=1.57d0
      endif
      del20=-0.12d0*prel/empi
      fact=80*paru(1)*paru(3)**2/(prel*prel)
      sig0=fact*sin(del00)**2
      sig1=3*fact*sin(del11)**2
      sig2=fact*sin(del20)**2

c...Calculate total and elastic cross sections.
 100  if(msel.eq.1) then
        if(srt.le.srtpi) then
          if(abs(iztot).eq.6) then       ! pi- pi-/pi+ pi+ 
            sig=sig2
            sigel=sig
            sigab=0.0d0
          else if(abs(iztot).eq.3) then  ! pi- pi0/pi0 pi+
            sigab=0.5d0*sig1
            sigel=0.5d0*sig2
            sig=sigab+sigel
          else if(abs(iztot).eq.0) then
            if(abs(iz1).eq.3) then       ! pi- pi+ 
              sigab=0.5d0*sig1+sig0/3.d0
              sigel=sig2/6.d0
              sig=sigel+sigab
            else                         ! pi0 pi0
              sigab=sig0/3.d0
              sigel=2.d0/3.d0*sig2
              sig=sigab+sigel
            end if
          endif
        else
          call jamxadq(kf1,kf2,sig,sigel)
        endif
        return
      endif

c....(2) Calculation of inel. channel.
      if(abs(iztot).eq.6) then  ! pi-pi-/pi+pi+
        mchanel=0
        mabsrb=0
        sigab=0.0d0
        if(msel.eq.3) then
          pare(3)=-10.d0
          if(srt.le.1.8d0) then
            ijet=-1
          else
            ijet=1
          endif
          return
        endif
      else if(iz1.eq.0.and.iz2.eq.0) then  ! pi0pi0 -> sigma
        mchanel=0
        mabsrb=1
        sigab=1.d0/3.d0*sig0
        if(msel.eq.3) then
          pare(3)=pare(3)-sigab
          if(pare(3).le.0.0d0) then
            kf1=10220
            kf2=0
            em1=srt
            em2=0.0d0
            return
          else
            pare(3)=-10.0d0
            if(srt.le.1.8d0) then
              ijet=-1
             else
              ijet=1
             endif
          endif
          return
        endif
      else if(abs(iz1+iz2).eq.3) then  ! pi-pi0/pi+pi0 -> rho-/rho+
        mchanel=0
        mabsrb=1
        sigab=0.5d0*sig1
        if(msel.eq.3) then
          pare(3)=pare(3)-sigab
          if(pare(3).le.0.0d0) then
            jj=isign(1,kf1)*isign(1,kf2)
            kf1=213*jj
            kf2=0
            em1=srt
            em2=0.0d0
            return
          else
            pare(3)=-10.0d0
            if(srt.le.1.8d0) then
              ijet=-1
             else
              ijet=1
             endif
          endif
          return
        endif
      else if(iz1+iz2.eq.0) then  ! pi-pi+ -> rho/sigma
        mchanel=0
        mabsrb=2
        sigab1=0.5d0*sig1
        sigab2=1.d0/3.d0*sig0
        sigab=sigab1+sigab2
        if(msel.eq.3) then
          pare(3)=pare(3)-sigab1
          if(pare(3).le.0.0d0) then
            kf1=113
            kf2=0
            em1=srt
            em2=0.0d0
            return
          endif
          pare(3)=pare(3)-sigab2
          if(pare(3).le.0.0d0) then
            kf1=10220
            kf2=0
            em1=srt
            em2=0.0d0
            return
          else
            pare(3)=-10.0d0
            if(srt.le.1.8d0) then
              ijet=-1
             else
              ijet=1
             endif
          endif
          return
        endif
      else
        write(check(1),'(i9,1x,i9,1x,i3,1x,i3)')kf1,kf2,iz1,iz2
        call jamerrm(30,1,'(jamcpipi:) ??? kf1 kf2 iz1 iz2=')
      endif

      end

c***********************************************************************

      subroutine jamcabb(msel,srt,pr,kf1,kf2,kc1,kc2,
     $            em1,em2,sig,sigel,sigin,mchanel,mabsrb,mxchan,ijet)

c...Purpose: to treat baryon-antibaryon cross sections.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      dimension sigin(mxchan),kfl1(3),kfl2(3)
   
      iz1=kchg(kc1,1)*isign(1,kf1)
      iz2=kchg(kc2,1)*isign(1,kf2)
c     id1=kchg(kc1,5)
c     id2=kchg(kc2,5)
      sig=0.0d0
      sigel=0.0d0
      snew=2*sqrt(parc(28)**2+pr**2)
      if(snew.gt.2*parc(28)) then
        plab=sqrt((snew**2-parc(28)**2-parc(28)**2)**2/
     $                             (4.d0*parc(28)**2)-parc(28)**2)
      else
        write(check(1),'(''plab pr='',2(g14.3,1x))')plab,pr
        call jamerrm(1,1,'(jamcabb:)plab<0')
        return
      endif

c...Get total and elastic cross sections.
      call jamxbbar(snew,plab,iz1,iz2,sig,sigel)

c...Rescale the cross section.
      call jamxadq(2212,2212,sig1,sigel1)
      call jamxadq(kf1,kf2,sig2,sigel2)
      sig=sig2/sig1*sig
      sigel=sigel2/sigel1*sigel
      sigann=min(sig-sigel,sig2/sig1*(24/plab**1.1d0+38/plab**0.5d0))


c...Calculate annihilation cross section.
      ianti=1
      if(kf2.lt.0) ianti=2

c...Quark contents.
      kf1a=abs(kf1)
      kfl1(3)=mod(kf1a/1000,10) 
      kfl1(2)=mod(kf1a/100,10) 
      kfl1(1)=mod(kf1a/10,10) 

      kf2a=abs(kf2)
      kfl2(3)=mod(kf2a/1000,10) 
      kfl2(2)=mod(kf2a/100,10) 
      kfl2(1)=mod(kf2a/10,10) 

      kfl1a=max(kfl1(1),kfl1(2),kfl1(3))
      kfl1c=min(kfl1(1),kfl1(2),kfl1(3))
      kfl1b=kfl1(1)+kfl1(2)+kfl1(3)-kfl1a-kfl1c

      kfl2a=max(kfl2(1),kfl2(2),kfl2(3))
      kfl2c=min(kfl2(1),kfl2(2),kfl2(3))
      kfl2b=kfl2(1)+kfl2(2)+kfl2(3)-kfl2a-kfl2c

      if(kfl1a.eq.kfl2a.and.kfl1b.eq.kfl2b) then
        kfla=kfl1c
        kflb=kfl2c
      else if(kfl1a.eq.kfl2a.and.kfl1c.eq.kfl2c) then
        kfla=kfl1b
        kflb=kfl2b
      else if(kfl1b.eq.kfl2b.and.kfl1c.eq.kfl2c) then
        kfla=kfl1a
        kflb=kfl2a
      else
        sig=sig-sigann
        sigann=0.0d0
      endif

      if(msel.eq.1) return

      mchanel=1
      mabsrb=1
      sigin(1)=max(0.0d0,sig-sigel-sigann)
      sigin(2)=sigann

      if(msel.eq.3) then

       pare(3)=pare(3)-sigann
       if(pare(3).le.0.0d0) then

         if(ianti.eq.1) then
          itmp=kfla
          kfla=kflb
          kflb=-itmp
         else
          kflb=-kflb
         endif
         call kfcnst(kfla,kflb,kf,0.0d0)
         if(kf.eq.0)then
           call jamerrm(30,0,'(jamcabb:)kf=0')
         endif
         kf1=kf
         kf2=0
         em1=srt
         em2=0
         ijet=2
         return
       endif

       pare(3)=-10.d0
       ijet=1
       return    
      endif


c      sigint=sig-sigel
c      pare(3)=pare(3)-sigint
c      if(pare(3).le.0.0) then
c        ijet=1
c        return
c      endif

c...Parametrization by
c...J.Gugnon and J. Vandermeulen, Ann. Phys. (Paris) 14, (1989)49.

c...Annihilation crosss section.
c     if(iann(1).ne.0.or.iann(2).ne.0.or.iann(3).ne.0) then
c       ian1=iann(1)*iann(2)
c       ian2=iann(2)*iann(3)
c       ian3=iann(3)*iann(1)
c       if(ian1.ne.0.and.iann(1).ne.iann(2)) then
c         kfan1=kfl2(iann(1))
c         kfan2=kfl2(iann(3))
c
c       if(iann(1).ne.iann(2).or. 
c       mabsrb=1
c       sigin(mabsrb)=24/plab**1.1+38/plab**0.5
c     endif
c
c....Charge exchange.
c     if(id1.eq.id_nucl.and.id2.eq.id_nucl) then
c       if(kfa1.eq.kfa2) then
c        mchanel=mchanel+1
c        if(plab.le.0.5) then
c          sigin(mchanel)=10.9*(plab-0.1)*p**(-1.6)
c        else
c          sigin(mchanel)=7.1*plab**(-0.9)
c        endif
c     endif

c...Non-strange production cross section.
c     sigin()=30*(plab-0.793)**1.5/(2+(plab-0.793)**1.5) 
c...Strange production ppbar-> yybar
c     sigin()=3*(plab-1.435)/(10.+(plab-1.435))


      end

c***********************************************************************

      subroutine jambmas(icl,srt,pr,ic,icc,kc1,kc2,kf1,kf2,iz01,iz02,
     $                              em1,em2,icon)

c...Purpose: to find outgoing types in BB inel. collisions.

c...icl     : =1:nn/n*n/n*n*,  =2:n(*)/d(*), =3:d(*)d(*)
c...srt     : invariant mass in GeV                     (input)
c...pr      : relative momentum(GeV/c) in two-body c.m. (input)
c...ic      : collision channel                         (input)
c...kc1,kc2 : ingoing particle KC codes                 (input)
c...kf1,kf2 : outgoing particle ID                      (input/output)
c...em1,em2 : outgoing particle masses                  (input/output)

      implicit double precision(a-h, o-z)
      include 'jam2.inc'

      dimension jz1(2,20),ipout1(2,20) ! pn
      dimension jz2(2,18),ipout2(2,18) ! pp/nn
      dimension iddn(2,11),jz3(3,11),jz4(11),jz5(2,11)
      dimension emind(4,2)
      character chap(4)*16
      parameter (emnuc=0.9383d0,empion=0.138d0)

c...1:nucl.    2: delta    3: ncul*   4: delta*
c...Outgoing particle IDs for pn collisions.
      data ipout1/2,1,  2,1,  3,1,  3,1,  2,2,
     a            2,2,  4,1,  4,1,  3,2,  3,2,
     a            2,4,  2,4,  3,3,  3,4,  3,4,
     a            4,4,  4,4,  3,1,  3,1,  1,1/

c...Outgoing charges for pn collisions.
      data jz1/   1,0,  0,1, 1,0,  0,1,  0,1,
     a           -1,2,  1,0, 0,1,  0,1,  1,0,
     a            0,1, -1,2, 0,1,  0,1,  1,0,
     a            0,1, -1,2, 1,0,  0,1,  0,1/

c...Outgoing particle IDs for pp collisions.
      data ipout2/2,1,  2,1,  3,1,  2,2,  2,2,
     a            4,1,  4,1,  3,2,  3,2,  2,4,
     a            2,4,  3,3,  3,4,  3,4,  4,4,
     a            4,4,  3,1,  1,1/

c...Outgoing charges for pp collisions.
      data jz2/   1,1,  2,0,  1,1, 1,1,  0,2,
     a            1,1,  2,0,  1,1, 0,2,  1,1,
     a            0,2,  1,1,  1,1, 0,2,  1,1,
     a            0,2,  1,1,  1,1/


c...Outgoing min. masses.
      data emind/
c    $ 1.0792000d0,1.0773d0,1.0792d0,1.0779000d0,
     $ 1.08000d0,1.08d0,1.08d0,1.0800d0,
     $ 1.08d0,1.08d0,1.08d0,1.08d0/
c    $ 1.22,1.22,1.22,1.22/

c    $ 1.07920003,1.07459998,1.0733,1.07790005,
c================================================
c 1)  np => d+(1232) + n   2) np => d0(1232) + p 
c 3)  np => np*            4) np => pn*
c 5)  np => d0d+           6) np => d-d++            
c 7)  np => nd*+           8) np => pd*0
c 9)) np=>n* d+           10) np=>p* d0
c 11) np=>d0 d*+          12) np=>d- d*++
c 13) np=>n* p*
c 14) np=>n* d*+          15) np=>p* d*0
c 16) np=>d*0 d*+         17) np=>d*- d*++
c 18,19) s-wave
c 20) -> nn 
c=======================================================
c 1) p+p => d+ + p /nn=>nd0  2) p+p => d++ + n /nn=>pd-
c 3) pp=>pp* /nn=>nn*
c 4) pp=>d+d+/nn=>d0d0       5) pp=>d0d++/nn=>d+d-
c 6) pp=>p d*+/nn=>nd*0 c    7) pp=>n d*++/nn=>pd*-
c 8) pp=>p* d+/nn=>n*d0 c    9) pp=>n* d++/nn=>p*d-
c 10)pp=>d+ d*+/nn=>d0d*0   11) pp=>d0 d*++/nn=>d+d*-
c 12)pp=>p* p*/nn=>n*n*
c 13)pp=>p* d*+/nn=>n*d0*   14) pp=>n* d*++/nn=>p*d*-
c 15)pp=>d*+d*+/nn=>d*0d*0  16) pp=>d*0d*++/nn=>d*+d*-
c 17) s-wave
c 18) ->nn


c...(1)BB->ND (2)BB->NN* (3)BB->DD (4)BB->ND* (5)BB->N*D
c...(6)BB->DD* (7)BB->N*N* (8)BB->N*D* (9)BB->D*D*
c...(10)s-wave (11)BB->N string (12)BB->NN
c
c...Outgoing states for n(*)n(*)/dd
      data iddn/2,1, 3,1, 2,2, 4,1, 2,3, 2,4, 3,3, 4,3, 4,4,
     $ 0,0,  1,1/

c...Outgoing charges for DN in total charge -1 or 3
      data jz3/-1,0,0,   99,99,0,  -1,0,1,    -1,0,0,  -1,0,0,
     a         -1,0,1,   99,99,0,  -1,0,0,  -1,0,1,
     a         99,99,0,  99,99,0/

c...Outgoing charges for dn total charge = 2 or 0
      data jz4/1,0,1,1,1,1,0,1,1,0,0/

c...Outgoing charges for dn total charge = 1.
      data jz5/0,1, 0,1, -1,2, 0,1, 0,1, -1,2, 0,1, 0,1, -1,2,
     $   0,1,  0,1/

c=======================================================================

      icon=0
      iswave=0
      iz1=iz01/3
      iz2=iz02/3
      id1=kchg(kc1,5)
      id2=kchg(kc2,5)
      izi1 = iz1
      izi2 = iz2
      iztot = iz1 + iz2
      ipair = jamcpair(id1,id2)
      id01=id1
      id02=id2
      kf01=kf1
      kf02=kf2
      em1o=em1
      em2o=em2

c======================================================================*
c...n+n /n+n*/ n*+n*  collsions
c======================================================================*

      if(icl.eq.1) then

c....pn
        if( iz1 .ne. iz2 ) then

           if(ic.ge.1.and.ic.le.20) then
              id1=ipout1(1,ic)
              id2=ipout1(2,ic)
              iz1=jz1(1,ic)
              iz2=jz1(2,ic)
              if(ic.ge.18.and.ic.le.19) iswave=1
            else
              write(check(1),'(''kf1 kf2'',i9,1x,i9)')kf1,kf2
              write(check(2),'(''ic='',i5)')ic
              call jamerrm(1,2,'(jambmas:)ic error nuc-nuc pn section')
              icon=-1
              return
            endif
c...pp
        else

          if(ic.ge.1.and.ic.le.18) then
             id1=ipout2(1,ic)
             id2=ipout2(2,ic)
             iz1=jz2(1,ic)
             iz2=jz2(2,ic)
             if(ic.eq.17) iswave=1
           else
             write(check(1),'(''kf1 kf2'',i9,1x,i9)')kf1,kf2
             write(check(2),'(''ic='',i5)')ic
             call jamerrm(1,2,'(jambmas:)ic error nuc-nuc pp section')
             icon=-1
             return
          endif

        endif

c...Convert charge in case of nn.
        if(iztot.eq.0) then
          iz1=1-iz1  
          iz2=1-iz2
        endif  

c======================================================================*
c      nucleon - delta ;  inelastic collisions
c
c...1) dn -> nn 2) dn -> nn*  3) dn -> dd 4) dn -> nd* 5) dn -> n*d
c   6) dn -> dd* 7) dn -> n*n* 8) dn -> n*d* 9) dn -> d*d*
c======================================================================*

      else if(icl.eq.2) then

        if(ic.eq.10) then
          write(check(1),'(''ic='',i5)')ic
          write(check(2),'(''kf1 kf2'',i9,1x,i9)')kf1,kf2
          write(check(3),'(''id1 id2'',i4,1x,i4)')id1,id2
          call jamerrm(1,3,'(jambmas:)invalid ic at dn1')
          icon=-1
          return
        else if(ic.ge.1.and.ic.le.11) then
          id1=iddn(1,ic)
          id2=iddn(2,ic)
        else
          write(check(1),'(''ic='',i5)')ic
          write(check(2),'(''kf1 kf2'',i9,1x,i9)')kf1,kf2
          write(check(3),'(''id1 id2'',i4,1x,i4)')id1,id2
          call jamerrm(1,3,'(jambmas:)invalid ic at dn1')
          icon=-1
          return
        endif

        randx=rn(0)

c......n + d-  or  p + d++ 
        if( iztot.eq.-1 .or. iztot.eq.3 ) then
             iz1=jz3(1,ic)
             iz2=jz3(2,ic)

            if(iz1.eq.99.and.iz2.eq.99) then
                call jamerrm(1,0,'invalid channel dn->nn???')
                icon=-1
                return
            endif

            if(jz3(3,ic).eq.1) then
              if( randx .ge. 0.5d0 ) then
                iztm = iz1
                iz1 = iz2
                iz2 = iztm
              endif 
            endif

c........  p + d+  or  n + d0 
        else if( iztot.eq.0 .or. iztot.eq.2 ) then
            iz1=1
            iz2=1
            if(jz4(ic).eq.1) then
              if(randx.le.0.5d0) then
                iz1 = 2
                iz2 = 0
              endif
            endif

c.....p  +  d0    or    n  +  d+
        else if(iztot.eq.1) then

          iz1=1
          iz2=0
          if(randx.le.0.5d0) then
             iz1=jz5(1,ic)
             iz2=jz5(2,ic)
          endif

        endif

        if(iztot.eq.1 .or.iztot.eq.0 .or. iztot.eq.3) then
          iz1 = 1 - iz1
          iz2 = 1 - iz2
        endif

c======================================================================*
c      delta - delta ;  inelastic collisions
c======================================================================*

      else if(icl.eq.3) then

        if(ic.eq.10) then
          write(check(1),'(''ic='',i5)')ic
          write(check(2),'(''kf1 kf2'',i9,1x,i9)')kf1,kf2
          write(check(3),'(''id1 id2'',i4,1x,i4)')id1,id2
          call jamerrm(1,3,'(bbtyp:) fatal error at DD channel(1)')
          icon=-1
          return
        else if(ic.ge.1.and.ic.le.11) then
            id1=iddn(1,ic)
            id2=iddn(2,ic)
        else
          write(check(1),'(''ic='',i5)')ic
          write(check(2),'(''kf1 kf2'',i9,1x,i9)')kf1,kf2
          write(check(3),'(''id1 id2'',i4,1x,i4)')id1,id2
          call jamerrm(1,3,'(bbtyp:) fatal error at DD channel(2)')
          icon=-1
          return
        endif

        if((iztot.eq.0).or.(iztot.eq.2)) then

          if(ic.eq.2 .or. ic.eq.7 .or.ic.eq.11) then
            iz1 = 0
            iz2 = 0
          else if(ic.eq.1.or.ic.eq.4 .or. ic.eq.5 .or. ic.eq.8 ) then 
            if(rn(0).le.0.5d0) then
              iz1 = 0
              iz2 = 0
            else
              iz1 = -1
              iz2 = 1
            endif
          else if( ic.eq.3 .or. ic.eq.6 .or. ic.eq.9 ) then 
            if(rn(0).le.0.5d0) then
              iz1 = 0
              iz2 = 0
            else
              iz1 = 1
              iz2 = -1
            endif
          endif
          if( iztot .eq. 2 ) then
            iz1 = 1 - iz1
            iz2 = 1 - iz2
          endif

        else if( iztot .eq. 1 ) then

          if(ic.eq.1.or.ic.eq.2.or.ic.eq.4.or.ic.eq.5
     $                .or.ic.eq.7.or.ic.eq.8.or.ic.eq.11 ) then
            if(rn(0).le.0.5d0) then
              iz1 = 0
              iz2 = 1
            else
              iz1 = 1
              iz2 = 0
            endif
          else if( ic.eq.3 .or. ic.eq.6 .or. ic.eq.9 ) then  ! DD
            if(rn(0).le.0.5d0) then
              iz1 = -1
              iz2 = 2
            else
              iz1 = 0
              iz2 = 1
            endif
          endif

        else if(iztot.eq.3) then  ! d++ d+
          iz1=2
          iz2=1
        else if(iztot.eq.-1) then  ! d- d0
          iz1=-1
          iz2=0
c       else if( iztot .ge. 3  .or. iztot .le. -1 ) then
        else

          if( ic.eq.3 .or. ic.eq.6 .or. ic.eq.9 ) then
          else
            write(check(1),'(''ic='',i5,''iztot='',i4)')ic,iztot
            write(check(2),'(''kf1 kf2'',i9,1x,i9)')kf1,kf2
            write(check(3),'(''iz1 iz2'',i4,1x,i4)')iz1,iz2
            call jamerrm(30,3,'(bbtyp:)error at dd charge=3/-1')
            icon=-1
            return
          endif

        endif

      else
        write(check(1),'(''ipair='',i4)')ipair
        write(check(2),'(''id1 id2='',i4,1x,i4)')id1,id2
        call jamerrm(30,2,'(bbtyp:)error invalid ipair')
      endif

c....Nucleon
      if(id1.eq.1) then
         iexc1=0
         if(iz1.eq.0) then
           kf1=2112
           em1min=0.9396d0
         else if(iz1.eq.1) then
           kf1=2212
           em1min=0.9383d0
         endif
         kc1=jamcomp(kf1)
c...Delta(1232)
      else if(id1.eq.2) then
         iexc1=100
         em1min=emind(iz1+2,1)+parc(41)
         if(iz1.eq.-1) kf1=1114
         if(iz1.eq.0) kf1=2114
         if(iz1.eq.1) kf1=2214
         if(iz1.eq.2) kf1=2224
         kc1=jamcomp(kf1)
c...N*
      else if(id1.eq.3) then
         if(iswave.eq.0) then
           iexc1=1
           if(iz1.eq.1) iexc1=2
c.....n* thresholds is asuumed to be 2*empion+emnuc
           em1min=emnuc+2*empion+parc(41)
         else
           iexc1=1000
           em1min=emnuc+empion+parc(41)
           kf1=12212
           if(iz1.eq.0) kf1=12112
           kc1=jamcomp(kf1)
         endif

c...Delta*
      else if(id1.eq.4) then
         if(iz1.eq.-1) iexc1=3
         if(iz1.eq. 0) iexc1=4
         if(iz1.eq. 1) iexc1=5
         if(iz1.eq. 2) iexc1=6
c        em1min=emind(iz1+2,2)+parc(41)
c.....d* thresholds is asuumed to be 3*empion+emnuc
         em1min=emnuc+3*empion+parc(41)
      endif


c...Same statement for 2
c....Nucleon
      if(id2.eq.1) then
         iexc2=0
         if(iz2.eq.0) then
           kf2=2112
           em2min=0.9396d0
         else if(iz2.eq.1) then
           kf2=2212
           em2min=0.9383d0
         endif
         kc2=jamcomp(kf2)

c...Delta(1232)
      else if(id2.eq.2) then
         iexc2=100
         em2min=emind(iz2+2,1)+parc(41)
         if(iz2.eq.-1) kf2=1114
         if(iz2.eq.0) kf2=2114
         if(iz2.eq.1) kf2=2214
         if(iz2.eq.2) kf2=2224
         kc2=jamcomp(kf2)
c...N*
      else if(id2.eq.3) then
         if(iswave.eq.0) then
           iexc2=1
           if(iz2.eq.1) iexc2=2
           em2min=emnuc+2*empion+parc(41)
         else
           em2min=emnuc+empion+parc(41)
           iexc2=1000
           kf2=12212
           if(iz2.eq.0) kf2=12112
           kc2=jamcomp(kf2)
         endif
c...Delta*
      else if(id2.eq.4) then
         if(iz2.eq.-1) iexc2=3
         if(iz2.eq. 0) iexc2=4
         if(iz2.eq. 1) iexc2=5
         if(iz2.eq. 2) iexc2=6
c        em2min=emind(iz2+2,2)+parc(41)
         em2min=emnuc+3*empion+parc(41)
      endif

c...Check minimum mass.
      if(srt.le.em1min+em2min+0.0001d0) then
        write(check(1),8100)srt,em1min,em2min
 8100   format('srt<em1+em2 srt em1min em2min',3(g12.3,1x))
        write(check(2),8110)ic,kf1,kf2
 8110   format('ic kf1 kf2=',i3,1x,i9,1x,i9)
        call jamerrm(1,2,'(jambmas:)minimum mass violate')
        kf1=kf01
        kf2=kf02
        icon=-1
        return
      endif

      iex=iexc1+iexc2
      em1=em1min
      em2=em2min
c...N+N final
      if(iex.eq.0) then

c....D(1232)+N final or  s-wave N(1440)+N or D(1232)+D(1232) final
      else if(iex.eq.100.or.iex.eq.200.or.iex.eq.1000) then
        call jamrmas1(iexc1,iexc2,kf1,kf2,kc1,kc2
     $     ,srt,em1,em2,em1min,em2min,icon)
        if(icon.ne.0) then
          kf1=kf01
          kf2=kf02
          em1=em1o
          em2=em2o
          return
        endif

      else

c...p(m)=breit-wigner
        call jamrmas1(iexc1,iexc2,kf1,kf2,kc1,kc2
     $                       ,srt,em1,em2,em1min,em2min,icon)
        if(icon.ne.0) then
          kf1=kf01
          kf2=kf02
          em1=em1o
          em2=em2o
          return
        endif

      endif

cxxxxxxxx Check minimum mass. xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      if(srt.le.em1+em2+0.0001d0) then
        write(6,*)'srt<em1+em2 srt em1 em2',srt,em1,em2
        write(6,*)'ic kf1 kf2=',ic,kf1,kf2
        kf1=kf01
        kf2=kf02
        icon=-1
        return
      endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
c...Make a random choice for the nucleon that will become a resonance.

3000  if(rn(0).gt.0.5d0) then
        itmp=kf1
        kf1=kf2
        kf2=itmp
        tmp=em1
        em1=em2
        em2=tmp
        itmp=iexc1
        iexc2=iexc1
        iexc1=itmp
      endif

c...Check charge conservation.
      izt2=iz1+iz2
      if(iztot.ne.izt2) then
        kc1f=jamcomp(kf1)
        kc2f=jamcomp(kf2)
        call pjname(kf01,chap(1))
        call pjname(kf02,chap(2))
        call pjname(kf1,chap(3))
        call pjname(kf1,chap(4))
        write(check(1),'(''ic='',i3,''iztot'',i4)')ic,iztot
        write(check(2),8000)chap(1),chap(2),chap(3),chap(4)
        write(check(3),'(''initial:izi1='',i3,''izi2='',i3)')izi1,izi2
        write(check(4),'(''final:  iz1 ='',i3,''iz2 ='',i3)')iz1,iz2
 8000   format(a16,1x,a16,' ==> ',a16,1x,a16)
        call jamerrm(30,4,'(jambmas:) Charge not conserved')
      endif 

      end

c***********************************************************************

      subroutine jamrmas1(iex1,iex2,kf1,kf2,kc1,kc2
     $ ,srt,em1,em2,emin1,emin2,icon)

c...Generate masses accroding to the Breit-Wigner distribution for
c...non-strange BB collisions.

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      parameter(maxbr=70)
      dimension ibranch(maxbr),pwid(maxbr)
      logical rel
      data isw/1/
      data rel /.true./

      icon=0
      emax2=0.0d0
      s=srt*srt
      pf0=sqrt((s-(emin1+emin2)**2)*(s-(emin1-emin2)**2))/(2.d0*srt)
      ntry=0
1000  continue

c...Particle1: find pole mass and max. mass.
      if(iex1.eq.0) then
        em01=em1
        emax1=em1
      else
        emax1=min(parc(52),srt-emin2-0.001d0)
        if(emax1.lt.emin1) then
          write(check(1),'(''(jamrmas1:)srt emax1<emin1'',
     $    3(1x,g11.3))')srt,emax1,emin1
          write(check(2),'(''emin2 emax2'',2(1x,g11.3),
     $    '' kf1 kf2'',2i10)')emin2,emax2,kf1,kf2
          call jamerrm(1,2,'(jamrmas1:)(1)invalid mass')
          icon=-1
          return
        endif
        
        if(mstc(64).le.1.and.(iex1.ge.1.and.iex1.le.6)) then
          call jambres2(srt,iex1,emax1,kf1,kc1)
        endif
        em01=pmas(kc1,1)     ! pole mass
      endif

c...Particle2: find pole mass and max. mass.
      if(iex2.eq.0) then
        em02=em2
        emax2=em2
      else
        emax2=min(parc(52),srt-emin1-0.001d0)
        if(emax2.lt.emin2) then
          write(check(1),'(''(jamrmas1:)srt emax2<emin2'',
     $    3(1x,g11.3))')srt,emax2,emin2
          write(check(2),'(''emin1 emax1'',2(1x,g11.3),
     $    '' kf1 kf2'',2i10)')emin1,emax1,kf1,kf2
          call jamerrm(1,2,'(jamrmas1:)(2)invalid mass')
          icon=-1
          return
        endif
        if(mstc(64).eq.1.and.(iex2.ge.1.and.iex2.le.6)) then
          call jambres2(srt,iex2,emax2,kf2,kc2)
        endif
        em02=pmas(kc2,1)     ! pole mass
      endif

c...Resonance prob. from integrated Breit-Wigner.
      if(mstc(64).ge.2) then
        call jambres1(srt,iex1,iex2,kf1,kf2,kc1,kc2)
        if(iex1.ge.1.and.iex1.le.6) em01=pmas(kc1,1)
        if(iex2.ge.1.and.iex2.le.6) em02=pmas(kc2,1)
      endif

      if(isw.eq.0) then
        gam1=pmas(kc1,2)
        gam2=pmas(kc2,2)
      endif
  
 100  continue
      ntry=ntry+1
      if(ntry.gt.100) goto 400

      if(iex1.ne.0) then
        itry1=0
 110    em1=emin1+(emax1-emin1)*rn(0) 
        itry1=itry1+1
        if(itry1.ge.300) then
          call jamerrm(1,0,'(jamrmas1:1) infinit loop? emin1 emax1')
          goto 400
        endif
        if(isw.eq.1)call jamwidm(kc1,1,0,0,0,em1,ibranch,pwid,gam1,itag)
        if(emax1.ge.em01) then
          if(rel) then
            bwmax1=2/gam1
          else
            bwmax1=4/gam1
          endif
        else
          if(rel) then
            bwmax1=2*gam1*em01**2/((emax1**2-em01**2)**2+(gam1*em01)**2)
          else
            bwmax1=gam1/((emax1-em01)**2+gam1**2/4.d0)
          endif
        endif
        if(rel) then
          bw=2*em1*em01*gam1/((em1**2-em01**2)**2+(em01*gam1)**2)
        else
          bw=gam1/((em1-em01)**2+gam1**2/4.d0)
        endif
        if(rn(0)*bwmax1.gt.bw) goto 110
      endif

      if(iex2.ne.0) then
        itry2=0
 200    em2=emin2+(emax2-emin2)*rn(0) 
        itry2=itry2+1
        if(itry2.ge.300) then
           call jamerrm(1,0,'(jamrmas1:2) infinit loop? emin2 emax2')
           goto 400
        endif
        if(isw.eq.1)call jamwidm(kc2,1,0,0,0,em2,ibranch,pwid,gam2,itag)
        if(emax2.ge.em02) then
          if(rel) then
            bwmax2=2/gam2
          else
            bwmax2=4/gam2
          endif
        else
          if(rel) then
            bwmax2=2*em02**2*gam2/((emax2**2-em02**2)**2+(em02*gam2)**2)
          else
            bwmax2=gam2/((emax2-em02)**2+gam2**2/4.d0)
          endif
        endif
        if(rel) then
          bw=2*em2*em02*gam2/((em2**2-em02**2)**2+(em02*gam2)**2)
        else
          bw=gam2/((em2-em02)**2+gam2**2/4.d0)
        endif
        if(rn(0)*bwmax2.gt.bw) goto 200
      endif

300   if(em1+em2+0.001d0.gt.srt) then
         if(mstc(64).ge.2) goto 100
         if(mstc(64).le.1) goto 1000
      endif
      pf=sqrt((s-(em1+em2)**2)*(s-(em1-em2)**2))/(2.d0*srt)
      if(rn(0).gt.pf/pf0) then
        if(mstc(64).ge.2) goto 100
        if(mstc(64).le.1) goto 1000
      endif
      return

400   continue
      if(iex1.ne.0.and.iex2.ne.0) then
        em1=emin1+(emax1-emin1)*rn(0) 
        em2=emin2+(min(parc(52),srt-em1-0.001d0)-emin2)*rn(0) 
      else if(iex1.ne.0) then
        em1=emin1+(emax1-emin1)*rn(0) 
      else if(iex2.ne.0) then
        em2=emin2+(emax2-emin2)*rn(0) 
      endif

      end

c***********************************************************************

      subroutine jamrmas2(kf1,kf2,kc1,kc2,srt,em1,em2,icon)

c...Generate masses accroding to the Breit-Wigner distribution.

c     implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'

c...Commonblock for t-channel resonance productions.
      common/jamres1/kfo(2,20),noutpa

c...Local arrays.
      parameter(maxbr=70)
      dimension ibranch(maxbr),icheck(50,50),ik(2)
      dimension ncount(50,2),icount(2),kfsign(2),kcv(2),kfv(2),iex(2)
      dimension pwid(maxbr),spin(2),xspin(2),emr(2),em0(2)
     $ ,emax(2),gamw(2),emin(2)
      character cerror*80
      data isw/1/

      izt0=kchg(kc1,1)*isign(1,kf1)+kchg(kc2,1)*isign(1,kf2)

      icon=0
      ncheck=0
      do i=1,50
      do j=1,50
        icheck(i,j)=0
      end do
      end do

c...Flag for particle-antiparticle rule.
      kfsign(1)=isign(1,kf1)
      kfsign(2)=isign(1,kf2)
      if(kchg(kc1,3).eq.0) kfsign(1)=1
      if(kchg(kc2,3).eq.0) kfsign(2)=1

c...Find minimum and maximum masses.
      call jamdmass(kf1,kfm1,kfd1,emin1,emdn1)
      call jamdmass(kf2,kfm2,kfd2,emin2,emdn2)
      emax1=srt-emdn2
      emax2=srt-emdn1
      iext=1
      if(kf1.eq.kfm1.and.kf2.eq.kfm2) iext=0

c...Find final states.
      call jamexpa(kf1,emax1,icount(1),ncount(1,1),ispin1,ig1)
      call jamexpa(kf2,emax2,icount(2),ncount(1,2),ispin2,ig2)

c...Max. spins in final states.
      spin(1)=dble(ispin1)
      spin(2)=dble(ispin2)

c...Energically not allowed.
      if(icount(1).eq.1.and.icount(2).eq.1) then
        kf1=kfm1
        kf2=kfm2
        em1=pmas(ncount(1,1),1)
        em2=pmas(ncount(1,2),1)
        return
      endif

c...Final max. c.m. momentum.
      s=srt*srt
      if(iext.eq.1) then
        if(srt.lt.emin1+emin2+0.001d0) then
          icon=1
          return
        endif
        pf0=sqrt((s-(emin1+emin2)**2)*(s-(emin1-emin2)**2))
     $    /(2.d0*srt)
      else
        ema=emin1+emdn2
        emb=emin2+emdn1
        if(ema.ge.emb) then
          if(srt.lt.ema+0.001d0) then
            icon=1
            return
          endif
          pf0=sqrt((s-ema**2)*(s-(emin1-emdn2)**2))/(2d0*srt)
        else
          if(srt.lt.emb+0.001d0) then
            icon=1
            return
          endif
          pf0=sqrt((s-emb**2)*(s-(emin2-emdn1)**2))/(2d0*srt)
        endif
      endif

c...Exclude the final states which are already checked before.
c...(e.g. explicit parametization of cross sections)
      mcheck=icount(1)*icount(2)
      if(noutpa.ge.1) then
      do ik1=1,icount(1)
      do ik2=1,icount(2)
        kcv1=ncount(ik1,1)
        kcv2=ncount(ik2,2)
        kfv1=kchg(kcv1,4)*kfsign(1)
        kfv2=kchg(kcv2,4)*kfsign(2)
        do i=1,noutpa
          if(kfo(1,i).eq.kfv1.and.kfo(2,i).eq.kfv2) then
            icheck(ik1,ik2)=-1
            mcheck=mcheck-1
          endif
        end do
      end do
      end do
      endif

c...Forbid elastic scattering.
      if(kf1.eq.kfm1.and.kf2.eq.kfm2) then
        icheck(ig1,ig2)=-1
        mcheck=mcheck-1
      endif

c...Check number of final states.
      if(mcheck.lt.1) then
       write(mstc(38),*)'(jamrmas2:a)There are no inelastic??'
     $,kf1,kf2,srt
       icon=1
       return
      endif

      mstd(199)=max(mstd(199),icount(1))
      mstd(199)=max(mstd(199),icount(2))

c...Start Monte Calro sampling.
      ntry=0
1000  continue
      ntry=ntry+1
      if(ntry.gt.100) then
        if(iex(1).ge.2.and.iex(2).ge.2) then
         emr(1)=emin(1)+(emax(1)-emin(1))*rn(0) 
         emr(2)=emin(2)+(srt-emr(1)-0.001d0-emin(2))*rn(0) 
        else if(iex(1).ge.2) then
         emr(1)=emin(1)+(emax(1)-emin(1))*rn(0) 
         emr(2)=em0(2)
        else if(iex(2).ge.2) then
         emr(2)=emin(2)+(emax(2)-emin(2))*rn(0) 
         emr(1)=em0(1)
        endif
        if(emr(1)+emr(2).gt.srt) then
          icon=1
        endif
        goto 600
      endif

      itry1=0
 1100 continue
      itry1=itry1+1
      if(itry1.gt.200) then
        ierror=2
        write(check(1),'(''kf1 kf2 srt'',2i10,1x,g10.3)')kf1,kf2,srt
        write(check(2),'(''xspin1 xspin2'',2(1x,g10.3))')
     $  xspin(1),xspin(2)
        cerror='(jamrmas2:)infinit loop at itry1'
        goto 9000
      endif

      itry2=0
  10  continue
      itry2=itry2+1
      if(itry2.ge.2500) then
        ierror=0
        cerror='(jamrmas2:B)There are no branches??'
        goto 9000
      endif
      ik(1)=1+rn(0)*icount(1)
      ik(2)=1+rn(0)*icount(2)
      if(icheck(ik(1),ik(2)).eq.-1) goto 10

c...Particle: find pole mass and max. mass.
      do jt=1,2
        kcv(jt)=ncount(ik(jt),jt)
        kfv(jt)=kchg(kcv(jt),4)*kfsign(jt)
        xspin(jt)=max(1,mod(abs(kfv(jt)),10))/spin(jt)
        if(kfv(jt).eq.221) xspin(jt)=xspin(jt)*0.3d0  ! eta
        if(kfv(jt).eq.331) xspin(jt)=xspin(jt)*0.4d0  ! eta'
        if(rn(0).gt.xspin(jt).and.icount(jt).ne.1) goto 10
        iex(jt)=1
        em0(jt)=pmas(kcv(jt),1) ! pole mass
        emr(jt)=em0(jt)
        emax(jt)=em0(jt)
        if(pmas(kcv(jt),2).gt.0.001d0) then
          iex(jt)=2
          emin(jt)=pmas(kcv(jt),1)-pmas(kcv(jt),3)
        endif
      end do

      do jt=1,2
        if(iex(jt).eq.2) then
          if(iex(3-jt).eq.1) then
            emin2=srt-em0(3-jt)-0.0001d0
          else if(iex(3-jt).eq.2) then
            emin2=srt-emin(3-jt)-0.0001d0
          endif
          emax(jt)=max(emin(jt),emin2)
        endif
      end do

c...Check energy.
      if(emax(1)+emax(2).gt.srt) then
        icheck(ik(1),ik(2))=-1
        ncheck=ncheck+1
        if(ncheck.eq.mcheck) then
          ierror=0
          cerror='(jamrmas2:)There are no inelastic??'
          goto 9000
        endif
        goto 1100
      endif


c...Get width of the resonances.
      if(isw.eq.0) then
        gamw(1)=pmas(kcv(1),2)
        gamw(2)=pmas(kcv(2),2)
      endif

c...Genarate resonance mass accroding to Breit-Wignwer distribution.
      do jt=1,2
      if(iex(jt).ge.2) then
        jtry1=0
 110    emr(jt)=emin(jt)+(emax(jt)-emin(jt))*rn(0) 
        jtry1=jtry1+1
        if(jtry1.ge.300) goto 1000

c.....Momentum dependent width.
        if(isw.eq.1)
     $    call jamwidm(kcv(1),1,0,0,0,emr(1),ibranch,pwid,gamw(jt),itag)
        if(gamw(jt).le.0.0d0) goto 1000

        if(emax(jt).ge.em0(jt)) then
          bwmax1=4/gamw(jt)
        else
          bwmax1=gamw(jt)/((emax(jt)-em0(jt))**2+gamw(jt)**2/4.d0)
        endif
        bw=gamw(jt)/((emr(jt)-em0(jt))**2+gamw(jt)**2/4.d0)
        if(rn(0)*bwmax1.gt.bw) goto 110
      endif
      end do

c...Check masses.
300   if(emr(1)+emr(2)+0.001d0.gt.srt) goto 1000

c...Check final phase space.
      pf=sqrt((s-(emr(1)+emr(2))**2)*(s-(emr(1)-emr(2))**2))/(2d0*srt)
      if(rn(0).gt.pf/pf0) goto 1000

c...Set final flavors and masses genarated.
600   continue
      kf1=kfv(1)
      kf2=kfv(2)
      kc1=kcv(1)
      kc2=kcv(2)
      em1=emr(1)
      em2=emr(2)

cxxxxxxxxxxxxxxx
c...Check total charge conservation.
      iz1=kchg(kc1,1)*isign(1,kf1)
      iz2=kchg(kc2,1)*isign(1,kf2)
      izt=iz1+iz2
      if(izt.ne.izt0) then
      write(6,*)'res2',srt,kcp(1,2),kcp(2,2),kf1,kf2,em1,em2
      write(6,*)'res2:',srt,kf1,kf2,'->',kfv1,kfv2,' '
     $ ,chaf(kc1,(3-isign(1,kf1))/2),' ',chaf(kc2,(3-isign(1,kf2))/2),
     $ ' -> ',chaf(kcv1,(3-isign(1,kf1))/2),
     $ ' ',chaf(kcv2,(3-isign(1,kf2))/2)
     $,em1,em2
      endif
cxxxxxxxxxxxxxxx

      return

 9000 continue
      icon=1

      if(mstc(8).eq.0.or.mstc(13).eq.0) return
      if(mstc(13).eq.1.and.(mstd(25).gt.mstc(14))) return
      mstd(25)=mstd(25)+1
      call jamerrm(1,ierror,cerror)
      ih=mstc(38)
      write(ih,*)'ntry kf1 kf2 srt',ntry,kf1,kf2,srt
      write(ih,*)'emin1 emax1',emin1,emax1
      write(ih,*)'emin2 emax2',emin2,emax2
      write(ih,*)'kfv1,kfv2',kfv(1),kfv(2)
      write(ih,*)'spin1 spin2',spin(1),spin(2)
      write(ih,*)'xspin1 xspin2',xspin(1),xspin(2)
      write(ih,*)'emin1 emax1 em0',emin(1),emax(1),em0(1)
      write(ih,*)'emin2 emax2 em0',emin(2),emax(2),em0(2)
      write(ih,*)'icount1 icount2',icount(1),icount(2)
      do jt=1,2
      write(ih,*)'jt=',jt
      do i=1,icount(jt)
        ii=ncount(i,jt)
        write(ih,*)i,kchg(ii,4),pmas(ii,1),pmas(ii,1)-pmas(ii,3)
     $   ,' ',chaf(ii,1)
      end do
      enddo
      do i=1,icount(1)
      do j=1,icount(2)
        i1=ncount(i,1)
        i2=ncount(j,2)
        write(ih,*)'icheck i j',i,j,icheck(i,j),kchg(i1,4),kchg(i2,4)
      end do
      end do


      end

c***********************************************************************

      subroutine jamexpa(kf0,em,icount,ncount,ispin,igr)

c...Purpose: to find possible hadronic excitation states.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'


c...Local arrays.
      dimension ncount(50)

      igr=0
      kfsign=isign(1,kf0)
      kc0=jamcomp(kf0)
      ispin=1
      if(kc0.ge.1) then
        id0=kchg(kc0,5)
        iz0=kchg(kc0,1)/3
        ibary=kchg(kc0,6)
      else
       write(check(1),'(''kf0='',i9)')kf0
       call jamerrm(30,1,'(jamexpa:) Unrecognized KF code')
      endif
      if(ibary.ne.3.and.ibary.ne.0) then
        write(check(1),'(''kf0='',i9)')kf0
        call jamerrm(30,1,'(jamexpa:) kf0=')
      endif

      if(ibary.eq.3) goto 200


c....Mesons.
        kfa=abs(kf0)
        kf1=mod(kfa/100,10)
        kf2=mod(kfa/10,10)
        kfm=100*kf1+10*kf2
        if(kchg(kc0,3).eq.0) then
          kfsign=1
        endif
        kfg=(kfm+1)*kfsign

        icount=0
        itry=0
1000    continue
        itry=itry+1

        do ir=0,5
        do is=1,7,2
          kfmes=(ir*10000+kfm+is)*kfsign
          kcm=jamcomp(kfmes)
c         do ip=1,noutpa
c           if(kfo(iopa,ip).eq.kfmes) goto 110
c         end do
          if(kcm.ne.0) then
          if((em.ge.pmas(kcm,1)-pmas(kcm,3)).or.kfmes.eq.kfg) then
            icount=icount+1
            ncount(icount)=kcm
            ispin=max(ispin,max(1,mod(abs(kfmes),10)))
            if(kfmes.eq.kfg) igr=icount
          endif
          endif
 110    end do
        end do

        jtry=0
        if(kfm.eq.220) then
           kfm=330
           jtry=1
        else if(kfm.eq.110) then
           kfm=220
           jtry=2
        else if(kfm.eq.330) then
           kfm=110
           jtry=3
        endif
        if(itry.le.2.and.jtry.ge.1) goto 1000

      return

c...Baryons.
 200  continue
      icount=0
      igr=1
      if(id0.eq.id_nucl.or.id0.eq.id_nucls
     $      .or.id0.eq.id_delt.or.id0.eq.id_delts) then
c....Find N* D* channel.

        if(iz0.eq.-1) kfd=1114
        if(iz0.eq.0)  kfd=2112
        if(iz0.eq.1)  kfd=2212
        if(iz0.eq.2)  kfd=2224

c       kfb=kfd*kfsign
c       icheck=0
c       do ip=1,noutpa
c         if(kfo(iopa,ip).eq.kfb) goto 210
c       end do

        icount=icount+1
        ncount(icount)=jamcomp(kfd)

 210   continue
        ispin=2
        if(iz0.eq.-1.or.iz0.eq.2) ispin=4

        if(iz0.eq.0.or.iz0.eq.1) then
          kcmin=mstc(22)+iz0
          kcmax=mstc(23)+iz0
          do kc=kcmin,kcmax,2
c           kfb=kchg(kc,4)*kfsign
c           do ip=1,noutpa
c             if(kfo(iopa,ip).eq.kfb) goto 220
c           end do
            if(em.ge.pmas(kc,1)-pmas(kc,3)) then
            icount=icount+1
            ncount(icount)=kc
c           ncount(icount)=kchg(kc,4)*kfsign
            ispin=max(ispin,max(1,mod(kchg(kc,4),10)))
            endif
  220     end do
        endif

c...Delta.
         if(iz0.eq.1) then
           icount=icount+1 
           ncount(icount)=jamcomp(2214)
         else if(iz0.eq.0) then
           icount=icount+1 
           ncount(icount)=jamcomp(2114)
         endif
         ispin=max(ispin,4)

c...Delta*
        kcmin=mstc(24)+iz0+1
        kcmax=mstc(25)+iz0+1
        istep=4
        do kc=kcmin,kcmax,istep
          if(em.ge.pmas(kc,1)-pmas(kc,3)) then
            icount=icount+1
            ncount(icount)=kc
            ispin=max(ispin,max(1,mod(kchg(kc,4),10)))
          endif
 230    end do

c...Lambda/Sigma
      else if(id0.eq.id_lamb.or.id0.eq.id_lambs
     $    .or.id0.eq.id_sigm.or.id0.eq.id_sigms) then

        icount=1
        if(iz0.eq.-1) ncount(icount)=jamcomp(3112)
        if(iz0.eq. 0) ncount(icount)=jamcomp(3122)
        if(iz0.eq. 1) ncount(icount)=jamcomp(3222)
        ispin=max(ispin,2)

c...Lambda*
        if(iz0.eq.0) then
          kcmin=mstc(26)
          kcmax=mstc(27)
          istep=1
          do kc=kcmin,kcmax,istep
            if(em.ge.pmas(kc,1)-pmas(kc,3)) then
            icount=icount+1
            ncount(icount)=kc
            ispin=max(ispin,max(1,mod(kchg(kc,4),10)))
            endif
 240      end do
        endif

        if(id0.eq.id_lamb.or.id0.eq.id_lambs) then
          icount=icount+1
          ncount(icount)=jamcomp(3212)
          ispin=max(ispin,2)
        endif

c...Sigma*
        kcmin=mstc(28)+iz0+1
        kcmax=mstc(29)+iz0+1
        istep=3
        icount=icount+1
c...S(1353)
        if(iz0.eq.-1) ncount(icount)=jamcomp(3114)
        if(iz0.eq. 0) ncount(icount)=jamcomp(3214)
        if(iz0.eq. 1) ncount(icount)=jamcomp(3224)
        ispin=max(ispin,4)
        do kc=kcmin,kcmax,istep
          if(em.ge.pmas(kc,1)-pmas(kc,3)) then
          icount=icount+1
          ncount(icount)=kc
          ispin=max(ispin,max(1,mod(kchg(kc,4),10)))
          endif
  250   end do

c...Xi*
      else if(id0.eq.id_xi.or.id0.eq.id_xis) then
        kcmin=mstc(30)+iz0+1
        kcmax=mstc(31)+iz0+1
        istep=2
        ispin=2
        icount=1
        if(iz0.eq.-1) then
          ncount(1)=jamcomp(3312)
          kc=jamcomp(3314)
          if(em.le.pmas(kc,1)-pmas(kc,3)) return
          ncount(2)=kc
        else if(iz0.eq. 0) then
          ncount(1)=jamcomp(3322)
          kc=jamcomp(3324)
          if(em.le.pmas(kc,1)-pmas(kc,3)) return
          ncount(2)=kc
        endif
        icount=2
        ispin=4
        do kc=kcmin,kcmax,istep
          if(em.gt.pmas(kc,1)-pmas(kc,3)) then
          icount=icount+1
          ncount(icount)=kc
          ispin=max(ispin,max(1,mod(kchg(kc,4),10)))
          endif
 260    end do

      else
        icount=1
        ncount(icount)=kc0
        ispin=max(1,mod(abs(kf0),10))
      endif


 2000 continue


      end


c*******************************************************************

      subroutine jamdmas(srt,em1,em2,imod,idlt2,emmax0,ic)

c...Generation of an allowed delta(1232)-mass
c
c    srt: c.m. enrgy in gev (input)
c    em1: delta mass (output)
c    em2: outgoing particle mass with delta (input/output)
c         in the case of nn->dd, em2 means the min.mass of
c         outgoing other delta mass.
c         i.e. em2=emnuc+empion should be used.
c   emmax0: limit of mass
c   idlt2 =0: 1delta production
c   idlt2 =1: 2delta production
c
c     ic: =0: fail to generate delta mass  (output)
c 
c imod=   0: constant delta width           (input)
c     others= momentum dependent delta width 
c         1: Randrap  ref: Randrup, NP A314 (1979) 429.
c                 Rittenberg, Rev.Mod.Phys. 43 (1971) s1.
c         2: Giessen
c         3: Kitazoe
c         4: Barz/Iwe
c note: all parametization is the same at low energy,
c       but at higher energy region, difference can be seen:
c       Giessen < Randrap < Kitazoe = Barz/Iwe
c       Kitazoe or Barz/Iwe parametrization can produce higher mass
c       delta than others.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision(a-h, o-z)
      parameter(bet2=0.090d0, qqr2=0.051936d0, gamr=0.11d0)
c...Functions: momentum in two-particle cm.
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a) 
      data  emnuc,emdelt,empion,widdlt/
     $      0.938d0,1.232d0,.138d0,   0.120d0 /
      data  const,const1, const2, ekinmi /
     $    -2.26d0,-1.154215d0,2.725011d0, 0.001d0/
      data  pscal1,pscal2,p0ref/0.238d0, 0.318d0, 0.227d0/

      ic=0
      emmax=min(emmax0,srt-em2-ekinmi)
      emmin=emnuc+empion+ekinmi
      if(emmax.le.emmin) then
        em1=emnuc+empion+ekinmi
        return
      end if

c....Delta mass with constant width
      if(imod.eq.0) then
c       pi=4.0d0*atan(1.0d0)
c       const=2.d0*(emnuc+empion-emdelt)/widdlt
c       const1=atan(const)
c       const2=pi/2.d0-const1
        xmax=(atan(2.d0*(emmax-emdelt)/widdlt)-const1)/const2
        xmax=min(1.0d0,xmax)
        x=xmax*rn(0)
        t=tan(x*const2)
        em1=emdelt+.5d0*widdlt*(const+t)/(1.d0-const*t)
        if(em1.gt.emmax) em1=emmax
        if(em1.lt.emmin) em1=emmin

c...Generation of 2 delta
        if(idlt2.ne.0) then
          emmax=min(emmax0,srt-em1-ekinmi)
          if(emmax.gt.3.0d0) emmax=3.0d0
          if(emmax.le.emmin) then
            em2=emnuc+empion+ekinmi
            return
          end if
          xmax=(atan(2.d0*(emmax-emdelt)/widdlt)-const1)/const2
          xmax=min(1.0d0,xmax)
          x=xmax*rn(0)
          t=tan(x*const2)
          em2=emdelt+.5d0*widdlt*(const+t)/(1.d0-const*t)
          if(em2.gt.emmax) em2=emmax
          if(em2.lt.emmin) em2=emmin
        endif
        if(em1.le.emnuc+empion) then
          call jamerrm(1,0,'(jamdmas:)delta mass small')
        endif
        return
      endif

c....Generation of delta mass with momentum dependent width
      itry=0
  100 continue
      itry=itry+1
      if(itry.ge.500.and.idlt2.eq.0) then
        write(6,*)'delta1 mass itry>500 imod',imod,em1,em2
        if(em1.le.emnuc+empion) then
          write(6,*)'(jamdlmas:)delta2 mass small em1 em2=',em1,em2
        endif
        return
      endif

      em1=emmin+(emmax-emmin)*rn(0)
      prel=pawt(em1,emnuc,empion)
      pp2=prel*prel

      if(imod.eq.1) then
c...Randrup
        gamd=widdlt*(prel**3/(1.d0+(prel/pscal1)**2+(prel/pscal2)**4))
     a     /(p0ref**3/(1.d0+(p0ref/pscal1)**2+(p0ref/pscal2)**4))

      else if(imod.eq.2) then
c...Giessen
            form= (1.d0+qqr2/bet2)/(1.d0+pp2/bet2)
            gamd=sqrt(pp2/qqr2)**3*emdelt/em1*gamr*form**2

      else if(imod.eq.3) then
c...Kitazoe
           gamd=0.47d0/(1.0d0+0.6d0*pp2/empion**2)*pp2/empion**2*prel
      else if(imod.eq.4) then
c...Barz/Iwe
         gamd=29.d0*prel**3/(1.d0+40.d0*pp2)
      else
           gamd=0.0d0
      endif

      prob=0.25d0*gamd*gamd/((em1-emdelt)**2+0.25d0*gamd*gamd)
      probmx=1.0d0
      if(emmax.lt.emdelt)
     a  probmx=0.25d0*gamd*gamd/((em1-emmax)**2+0.25d0*gamd*gamd)

      if(prob/probmx.gt.1.0d0) then
         write(6,*) 'jamdlmas-warning: 1.0d0< prob=',prob,probmx
      end if

      if(rn(0)*probmx.gt.prob) goto 100
c     write(6,*)'imod idlt2',imod,idlt2
c     write(6,*)'srt,em1,em2 emmax em1',srt,em1,em2,emmax,em1
      if(em1.le.emnuc+empion) then
        write(6,*)'(jamdlmas:)3delta mass small em1 em2=',em1,em2
      endif
      if(idlt2.eq.0) return


c...Generation of 2-delta mass
      emmax=min(emmax0,srt-em1-ekinmi)
        if(emmax.gt.3.0d0) emmax=3.0d0
      if(emmax.le.emmin) then
         em2=emnuc+empion+ekinmi
         return
      end if
      em2=emmin+(emmax-emmin)*rn(0)
      prel=pawt(em2,emnuc,empion)
      pp2=prel*prel
      if(imod.eq.1) then
c...Randrup
        gamd=widdlt*(prel**3/(1.d0+(prel/pscal1)**2+(prel/pscal2)**4))
     a     /(p0ref**3/(1.d0+(p0ref/pscal1)**2+(p0ref/pscal2)**4))

      else if(imod.eq.2) then
c...Giessen
            form= (1.d0+qqr2/bet2)/(1.d0+pp2/bet2)
            gamd=sqrt(pp2/qqr2)**3*emdelt/em2*gamr*form**2

      else if(imod.eq.3) then
c...Kitazoe
           gamd=0.47d0/(1.0d0+0.6d0*pp2/empion**2)*pp2/empion**2*prel
      else if(imod.eq.4) then
c...Barz/Iwe
         gamd=29.d0*prel**3/(1.d0+40.d0*pp2)
      else
           gamd=0.0d0
      endif

       prob=0.25d0*gamd*gamd/((em2-emdelt)**2+0.25d0*gamd*gamd)
       probmx=1.0d0
       if(emmax.lt.emdelt)
     a probmx=0.25d0*gamd*gamd/((em2-emmax)**2+0.25d0*gamd*gamd)

      if(prob/probmx.gt.1.0d0) then
         write(6,*) 'jamdlmass-warning: 1.0d0< prob=',prob,probmx
      end if

      if(rn(0)*probmx.gt.prob) then
       if(itry.ge.2000) then
        write(6,*)'2delta mass itry>200 imod',imod,em2
        return
       endif
       goto 100
      endif

      end


c***********************************************************************

      subroutine jamxadq(kf1,kf2,sig,sigel)

c...Calculate cross section from additive quark model.
c...Ref;K.Goulianos,Phys.Rep.101,169(1983),H.Sorge,Z.P.C59(1993)85.

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      parameter(signn=40.0d0,sigelnn=10.0d0)

      kfa=abs(kf1) 
      kfla=mod(kfa/1000,10) 
      kflb=mod(kfa/100,10) 
      kflc=mod(kfa/10,10) 
      if(kfa.le.10) kflc=kfa
      str1=0
      nq1=0
      if(kfla.eq.3) str1=str1+1
      if(kflb.eq.3) str1=str1+1
      if(kflc.eq.3) str1=str1+1
      if(kfla.ne.0) nq1=nq1+1
      if(kflb.ne.0) nq1=nq1+1
      if(kflc.ne.0) nq1=nq1+1
      nq1=max(1,nq1)

      kfa=abs(kf2) 
      kfla=mod(kfa/1000,10) 
      kflb=mod(kfa/100,10) 
      kflc=mod(kfa/10,10) 
      if(kfa.le.10) kflc=kfa
      str2=0
      nq2=0
      if(kfla.eq.3) str2=str2+1
      if(kflb.eq.3) str2=str2+1
      if(kflc.eq.3) str2=str2+1
      if(kfla.ne.0) nq2=nq2+1
      if(kflb.ne.0) nq2=nq2+1
      if(kflc.ne.0) nq2=nq2+1
      if(kf2.eq.21) nq2=1
      nq2=max(1,nq2)

      sig=signn*(nq1/3.d0)*(nq2/3.d0)*
     &        (1.d0-0.4d0*str1/nq1)*(1.d0-0.4d0*str2/nq2)
      sigel=sig**1.5*(sigelnn/signn**1.5)

      end

c*********************************************************************
 
      subroutine jamxtot(kf01,kf02,srt,pr,sigt)
 
C...Parametrizes total, elastic and diffractive cross-sections
C...for different energies and beams. Donnachie-Landshoff for
C...total and Schuler-Sjostrand for elastic and diffractive.
c...sigt(0,0,0) total
c...sigt(0,0,1) elastic
c...sigt(0,0,2) a+b-> x+b
c...sigt(0,0,3) a+b-> a+x
c...sigt(0,0,4) double diffractive
c...sigt(0,0,5) non-diffractive
C...Process code IPROC:
C...=  1 : p + p;
C...=  2 : pbar + p;
C...=  3 : pi+ + p;
C...=  4 : pi- + p;
C...=  5 : pi0 + p;
C...=  6 : phi + p;
C...=  7 : J/psi + p;
C...= 11 : rho + rho;
C...= 12 : rho + phi;
C...= 13 : rho + J/psi;
C...= 14 : phi + phi;
C...= 15 : phi + J/psi;
C...= 16 : J/psi + J/psi;
C...= 21 : gamma + p (DL);
C...= 22 : gamma + p (VDM).
C...= 23 : gamma + pi (DL);
C...= 24 : gamma + pi (VDM);
C...= 25 : gamma + gamma (DL);
C...= 26 : gamma + gamma (VDM).

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      dimension sigt(0:6,0:6,0:5)
      save /pjpars/,/pjint1/
      dimension nproc(30),xpar(30),ypar(30),ihada(20),ihadb(20),
     & pmhad(4),bhad(4),betp(4),ifitsd(20),ifitdd(20),ceffs(10,10),
     & ceffd(10,10),sigtmp(6,0:5)

c...Function:lab. momentum.
      plabsr(a,b,c)=sqrt((a**2-b**2-c**2)**2/(4.d0*c**2)-b**2)
 
C...Common constants.
      data eps/0.0808d0/, eta/-0.4525d0/, alp/0.25d0/, cres/2.d0/, 
     &  pmrc/1.062d0/,
     &smp/0.880d0/, facel/0.0511d0/, facsd/0.0336d0/, facdd/0.0084d0/
 
C...Number of multiple processes to be evaluated (= 0 : undefined).
      data nproc/7*1,3*0,6*1,4*0,4*3,2*6,4*0/

C...X and Y parameters of sigmatot = X * s**epsilon + Y * s**(-eta).
      data xpar/2*21.70d0,3*13.63d0,10.01d0,0.970d0,3*0.d0,
     &8.56d0,6.29d0,0.609d0,4.62d0,0.447d0,0.0434d0,4*0.d0,
     &0.0677d0,0.0534d0,0.0425d0,0.0335d0,2.11d-4,1.31d-4,4*0.d0/
      data ypar/56.08d0,98.39d0,27.56d0,36.02d0,31.79d0,-1.51d0, 
     & -0.146d0,3*0.d0,
     &13.08d0,-0.62d0,-0.060d0,0.030d0,-0.0028d0,0.00028d0,4*0.d0,
     &0.129d0,0.115d0,0.081d0,0.072d0,2.15d-4,1.70d-4,4*0.d0/
 
C...Beam and target hadron class:
C...= 1 : p/n ; = 2 : pi/rho/omega; = 3 : phi; = 4 : J/psi.
      data ihada/2*1,3*2,3,4,3*0,3*2,2*3,4,4*0/
      data ihadb/7*1,3*0,2,3,4,3,2*4,4*0/
C...Characteristic class masses, slope parameters, beta = sqrt(X).
      data pmhad/0.938d0,0.770d0,1.020d0,3.097d0/
      data bhad/2.3d0,1.4d0,1.4d0,0.23d0/
      data betp/4.658d0,2.926d0,2.149d0,0.208d0/
 
C...Fitting constants used in parametrizations of diffractive results.
      data ifitsd/2*1,3*2,3,4,3*0,5,6,7,8,9,10,4*0/
      data ifitdd/2*1,3*2,3,4,3*0,5,6,7,8,9,10,4*0/
      data ((ceffs(j1,j2),j2=1,10),j1=1,10)/
     & 0.213d0, 0.0d0, -0.47d0, 150.d0, 0.213d0, 0.0d0, -0.47d0, 150.d0,
     &  0.d0, 0.d0,
     & 0.213d0, 0.0d0, -0.47d0, 150.d0, 0.267d0, 0.0d0, -0.47d0, 100.d0,
     &  0.d0, 0.d0,
     & 0.213d0, 0.0d0, -0.47d0, 150.d0, 0.232d0, 0.0d0, -0.47d0, 110.d0,
     &  0.d0, 0.d0,
     & 0.213d0, 7.0d0, -0.55d0, 800.d0, 0.115d0, 0.0d0, -0.47d0, 110.d0,
     &  0.d0, 0.d0,
     & 0.267d0, 0.0d0, -0.46d0,  75.d0, 0.267d0, 0.0d0, -0.46d0,  75.d0,
     &  0.d0, 0.d0,
     & 0.232d0, 0.0d0, -0.46d0,  85.d0, 0.267d0, 0.0d0, -0.48d0, 100.d0,
     &  0.d0, 0.d0,
     & 0.115d0, 0.0d0, -0.50d0,  90.d0, 0.267d0, 6.0d0, -0.56d0, 420.d0,
     &  0.d0, 0.d0,
     & 0.232d0, 0.0d0, -0.48d0, 110.d0, 0.232d0, 0.0d0, -0.48d0, 110.d0,
     &  0.d0, 0.d0,
     & 0.115d0, 0.0d0, -0.52d0, 120.d0, 0.232d0, 6.0d0, -0.56d0, 470.d0,
     &  0.d0, 0.d0,
     & 0.115d0, 5.5d0, -0.58d0, 570.d0, 0.115d0, 5.5d0, -0.58d0, 570.d0,
     &  0.d0, 0.d0/
      data ((ceffd(j1,j2),j2=1,10),j1=1,10)/
     & 3.11d0, -7.34d0,  9.71d0, 0.068d0, -0.42d0, 1.31d0, -1.37d0, 
     &   35.0d0,  118.d0, 0.d0,
     & 3.11d0, -7.10d0,  10.6d0, 0.073d0, -0.41d0, 1.17d0, -1.41d0, 
     &   31.6d0,   95.d0, 0.d0,
     & 3.12d0, -7.43d0,  9.21d0, 0.067d0, -0.44d0, 1.41d0, -1.35d0, 
     &   36.5d0,  132.d0, 0.d0,
     & 3.13d0, -8.18d0, -4.20d0, 0.056d0, -0.71d0, 3.12d0, -1.12d0, 
     &   55.2d0, 1298.d0, 0.d0,
     & 3.11d0, -6.90d0,  11.4d0, 0.078d0, -0.40d0, 1.05d0, -1.40d0, 
     &   28.4d0,   78.d0, 0.d0,
     & 3.11d0, -7.13d0,  10.0d0, 0.071d0, -0.41d0, 1.23d0, -1.34d0, 
     &   33.1d0,  105.d0, 0.d0,
     & 3.12d0, -7.90d0, -1.49d0, 0.054d0, -0.64d0, 2.72d0, -1.13d0, 
     &   53.1d0,  995.d0, 0.d0,
     & 3.11d0, -7.39d0,  8.22d0, 0.065d0, -0.44d0, 1.45d0, -1.36d0, 
     &   38.1d0,  148.d0, 0.d0,
     & 3.18d0, -8.95d0, -3.37d0, 0.057d0, -0.76d0, 3.32d0, -1.12d0, 
     &   55.6d0, 1472.d0, 0.d0,
     & 4.18d0, -29.2d0,  56.2d0, 0.074d0, -1.36d0, 6.67d0, -1.14d0, 
     &  116.2d0, 6532.d0, 0.d0/
 

C...Parameters. Combinations of the energy.
      aem=paru(101)  ! alpha_em: electromagnetic fine structure const.
      pmth=parp(102) ! D=0.28GeV: the mass of diffractive states.

c     s=vint(2)
c     srt=vint(1)

      s=srt*srt
      seps=s**eps
      seta=s**eta
      slog=log(s)
 
C...Ratio of gamma/pi (for rescaling in structure functions).
c     vint(281)=(xpar(22)*seps+ypar(22)*seta)/
c    &                    (xpar(5)*seps+ypar(5)*seta)
c
c     if(mint(50).ne.1) return

C...Find process number (for lookup tables).
      kf1=kf01
      kf2=kf02
      kc1=jamcomp(kf1)
      kc2=jamcomp(kf2)
      ibar1=kchg(kc1,6)*isign(1,kf1)
      ibar2=kchg(kc2,6)*isign(1,kf2)
      kf1a=abs(kf1)
      kf2a=abs(kf2)
      iord=1
      icross=0
      factor=1.0d0
      kfla1=mod(kf1a/1000,10) 
      kflb1=mod(kf1a/100,10) 
      kflc1=mod(kf1a/10,10) 
      kfla2=mod(kf2a/1000,10) 
      kflb2=mod(kf2a/100,10) 
      kflc2=mod(kf2a/10,10) 

      if(ibar1*ibar2.eq.9) then
        iproc=1
        if(srt.le.10.d0) then
          call jamxnn(srt,3,3,sig,sigel)
          icross=1
        endif
      else if( ibar1 .eq. 0 .and. ibar2 .eq. 0 ) then
        kfm1=10*kflb1+kflc1
        kfm2=10*kflb2+kflc2
        iproc=11
        if(kfm1.eq.33.and.kfm2.le.30) then
          iproc=12
          iord=2
        else if(kfm2.eq.33.and.kfm1.le.30) then
          iproc=12
          iord=2
        else if(kfm1.eq.33.and.kfm2.eq.33) then
          iproc=14
        else if(kfm1.eq.44.and.kfm2.eq.44) then
          iproc=16
        endif
        if(srt.le.10.0d0) then
          call jamxadq(kf1,kf2,sig,sigel)
          icross=1
        endif

      else if( ibar1 * ibar2 .eq. 0 ) then
        if(ibar1.lt.0.or.ibar2.lt.0) then
          if(kchg(kc1,3).ne.0) kf1=-kf1
          if(kchg(kc2,3).ne.0) kf2=-kf2
        endif
        if(ibar1.eq.0) then
         iz1=jamchge(kf1)/3
         iz2=jamchge(kf2)/3
         kfm1=10*kflb1+kflc1
         kfm=kfm1
        else
         iz1=jamchge(kf2)/3
         iz2=jamchge(kf1)/3
         iord=2
         kfm2=10*kflb2+kflc2
         kfm=kfm2
        endif
        iproc=3
c.....pi- n/pi+ p
        if((iz1.eq.-1.and.iz2.eq.0).or.(iz1.eq.1.and.iz2.eq.1))iproc=3
c.....pi+ n/pi- p
        if((iz1.eq.-1.and.iz2.eq.1).or.(iz1.eq.1.and.iz2.eq.0))iproc=4
c.....pi0 n/pi0 p
        if((iz1.eq.0.and.iz2.eq.0).or.(iz1.eq.0.and.iz2.eq.1))iproc=5
        if(kfm.eq.33) iproc=6
        if(kfm.eq.44) iproc=7
        if(srt.le.10.0d0) then
          call jamxadq(kf1,kf2,sig,sigel)
          icross=1
        endif
c       if(kfm1.eq.32.or.kfm2.eq.31) factor=0.82

      else if(ibar1*ibar2.eq.-9) then
        iproc=2
        if(kf1.gt.0) iord=2
        snew=2*sqrt(parc(28)**2+pr**2)
        plab=plabsr(snew,parc(28),parc(28))
        if(snew.lt.10.0d0) then
          iz1=jamchge(kf1)
          iz2=jamchge(kf2)
          call jamxbbar(snew,plab,iz1,iz2,sig,sigel)
          icross=1
        endif
      else
       call jamerrm(30,0,'(jamxtot:) invalid kf1 kf2')
      endif

      if(iproc.le.0.or.iproc.gt.26) then
       call jamerrm(30,0,'(jamxtot:) invalid iproc')
      endif
 
C...Order flavours of incoming particles: KF1 < KF2.
c     if(iabs(mint(11)).le.iabs(mint(12))) then
c       kf1=iabs(mint(11))
c       kf2=iabs(mint(12))
c       iord=1
c     else
c       kf1=iabs(mint(12))
c       kf2=iabs(mint(11))
c       iord=2
c     endif

c     isgn12=isign(1,mint(11)*mint(12))
C...Find process number (for lookup tables).
c     if(kf1.gt.1000) then
c       iproc=1
c       if(isgn12.lt.0) iproc=2
c     elseif(kf1.gt.100.and.kf2.gt.1000) then
c       iproc=3
c       if(isgn12.lt.0) iproc=4
c       if(kf1.eq.111) iproc=5
c     elseif(kf1.gt.100) then
c       iproc=11
c     elseif(kf2.gt.1000) then
c       iproc=21
c       if(mint(123).eq.2) iproc=22
c     elseif(kf2.gt.100) then
c       iproc=23
c       if(mint(123).eq.2) iproc=24
c     else
c       iproc=25
c       if(mint(123).eq.2) iproc=26
c     endif
 

C... Number of multiple processes to be stored; beam/target side.
      npr=nproc(iproc)
      mint(101)=1
      mint(102)=1
      if(npr.eq.3) then
        mint(100+iord)=4
      elseif(npr.eq.6) then
        mint(101)=4
        mint(102)=4
      endif
      n1=0
      if(mint(101).eq.4) n1=4
      n2=0
      if(mint(102).eq.4) n2=4
 
C...Do not do any more for user-set or undefined cross-sections.
c     if(mstp(31).le.0) return
c
c     if(npr.eq.0) call luerrm(26,
c    &'(PYXTOT:) cross section for this process not yet implemented')
 
C...Parameters. Combinations of the energy.
c     aem=paru(101)
c     pmth=parp(102)
c     s=vint(2)
c     srt=vint(1)
c     seps=s**eps
c     seta=s**eta
c     slog=log(s)
 
C...Loop over multiple processes (for VDM).
      do 110 i=1,npr
      if(npr.eq.1) then
        ipr=iproc
      elseif(npr.eq.3) then
        ipr=i+4
        if(kf2.lt.1000) ipr=i+10
      elseif(npr.eq.6) then
        ipr=i+10
      endif
 
C...Evaluate hadron species, mass, slope contribution and fit number.
      iha=ihada(ipr)
      ihb=ihadb(ipr)
      pma=pmhad(iha)
      pmb=pmhad(ihb)
      bha=bhad(iha)
      bhb=bhad(ihb)
      isd=ifitsd(ipr)
      idd=ifitdd(ipr)
 
C...Skip if energy too low relative to masses.
      do 100 j=0,5
      sigtmp(i,j)=0.d0
  100 continue
      if(srt.lt.pma+pmb+0.2d0) then
        goto 110
      endif
 
C...Total cross-section. Elastic slope parameter and cross-section.
      bel=2.d0*bha+2.d0*bhb+4.d0*seps-4.2d0
      if(icross.eq.0) then
        sigtmp(i,0)=xpar(ipr)*seps+ypar(ipr)*seta*factor
        sigtmp(i,1)=facel*sigtmp(i,0)**2/bel
      else
        sigtmp(i,0)=sig
        sigtmp(i,1)=sigel
      endif
 
C...Diffractive scattering A + B -> X + B.
      bsd=2.d0*bhb
      sqml=(pma+pmth)**2
      sqmu=s*ceffs(isd,1)+ceffs(isd,2)
      sum1=log((bsd+2.d0*alp*log(s/sqml))/
     &(bsd+2.d0*alp*log(s/sqmu)))/(2.d0*alp)
      bxb=ceffs(isd,3)+ceffs(isd,4)/s
      sum2=cres*log(1.d0+((pma+pmrc)/(pma+pmth))**2)/
     &(bsd+2.d0*alp*log(s/((pma+pmth)*(pma+pmrc)))+bxb)
      sigtmp(i,2)=facsd*xpar(ipr)*betp(ihb)*max(0.d0,sum1+sum2)
 
C...Diffractive scattering A + B -> A + X.
      bsd=2.d0*bha
      sqml=(pmb+pmth)**2
      sqmu=s*ceffs(isd,5)+ceffs(isd,6)
      sum1=log((bsd+2.d0*alp*log(s/sqml))/
     &(bsd+2.d0*alp*log(s/sqmu)))/(2.d0*alp)
      bax=ceffs(isd,7)+ceffs(isd,8)/s
      sum2=cres*log(1.d0+((pmb+pmrc)/(pmb+pmth))**2)/
     &(bsd+2.d0*alp*log(s/((pmb+pmth)*(pmb+pmrc)))+bax)
      sigtmp(i,3)=facsd*xpar(ipr)*betp(iha)*max(0.d0,sum1+sum2)
 
C...Order single diffractive correctly.
      if(iord.eq.2) then
        sigsav=sigtmp(i,2)
        sigtmp(i,2)=sigtmp(i,3)
        sigtmp(i,3)=sigsav
      endif
 
C...Double diffractive scattering A + B -> X1 + X2.
      yeff=log(s*smp/((pma+pmth)*(pmb+pmth))**2)
      deff=ceffd(idd,1)+ceffd(idd,2)/slog+ceffd(idd,3)/slog**2
      sum1=deff+yeff*(log(max(1d-10,yeff/deff))-1.d0)/(2.d0*alp)
      if(yeff.le.0) sum1=0.d0
      sqmu=s*(ceffd(idd,4)+ceffd(idd,5)/slog+ceffd(idd,6)/slog**2)
      slup=log(max(1.1d0,s/(alp*(pma+pmth)**2*(pmb+pmth)*(pmb+pmrc))))
      sldn=log(max(1.1d0,s/(alp*sqmu*(pmb+pmth)*(pmb+pmrc))))
      sum2=cres*log(1.d0+((pmb+pmrc)/(pmb+pmth))**2)*log(slup/sldn)/
     &(2.d0*alp)
      slup=log(max(1.1d0,s/(alp*(pmb+pmth)**2*(pma+pmth)*(pma+pmrc))))
      sldn=log(max(1.1d0,s/(alp*sqmu*(pma+pmth)*(pma+pmrc))))
      sum3=cres*log(1.d0+((pma+pmrc)/(pma+pmth))**2)*log(slup/sldn)/
     &(2.d0*alp)
      bxx=ceffd(idd,7)+ceffd(idd,8)/srt+ceffd(idd,9)/s
      slrr=log(s/(alp*(pma+pmth)*(pma+pmrc)*(pmb+pmth)*(pmb*pmrc)))
      sum4=cres**2*log(1.d0+((pma+pmrc)/(pma+pmth))**2)*
     &log(1.d0+((pmb+pmrc)/(pmb+pmth))**2)/max(0.1d0,2.d0*alp*slrr+bxx)
      sigtmp(i,4)=facdd*xpar(ipr)*max(0.d0,sum1+sum2+sum3+sum4)
 
C...Non-diffractive by unitarity.
      sigtmp(i,5)=sigtmp(i,0)-sigtmp(i,1)-sigtmp(i,2)-sigtmp(i,3)-
     &sigtmp(i,4)
  110 continue
 
C...Put temporary results in output array: only one process.
      if(mint(101).eq.1.and.mint(102).eq.1) then
        do 120 j=0,5
        sigt(0,0,j)=sigtmp(1,j)
  120   continue
 
C...Beam multiple processes.
      elseif(mint(101).eq.4.and.mint(102).eq.1) then
        do 140 i=1,4
        conv=aem/parp(160+i)
        i1=max(1,i-1)
        do 130 j=0,5
        sigt(i,0,j)=conv*sigtmp(i1,j)
  130   continue
  140   continue
        do 150 j=0,5
        sigt(0,0,j)=sigt(1,0,j)+sigt(2,0,j)+sigt(3,0,j)+sigt(4,0,j)
  150   continue
 
C...Target multiple processes.
      elseif(mint(101).eq.1.and.mint(102).eq.4) then
        do 170 i=1,4
        conv=aem/parp(160+i)
        iv=max(1,i-1)
        do 160 j=0,5
        sigt(0,i,j)=conv*sigtmp(iv,j)
  160   continue
  170   continue
        do 180 j=0,5
        sigt(0,0,j)=sigt(0,1,j)+sigt(0,2,j)+sigt(0,3,j)+sigt(0,4,j)
  180   continue
 
C...Both beam and target multiple processes.
      else
        do 210 i1=1,4
        do 200 i2=1,4
        conv=aem**2/(parp(160+i1)*parp(160+i2))
        if(i1.le.2) then
          iv=max(1,i2-1)
        elseif(i2.le.2) then
          iv=max(1,i1-1)
        elseif(i1.eq.i2) then
          iv=2*i1-2
        else
          iv=5
        endif
        do 190 j=0,5
        jv=j
        if(i2.gt.i1.and.(j.eq.2.or.j.eq.3)) jv=5-j
        sigt(i1,i2,j)=conv*sigtmp(iv,jv)
  190   continue
  200   continue
  210   continue
        do 230 j=0,5
        do 220 i=1,4
        sigt(i,0,j)=sigt(i,1,j)+sigt(i,2,j)+sigt(i,3,j)+sigt(i,4,j)
        sigt(0,i,j)=sigt(1,i,j)+sigt(2,i,j)+sigt(3,i,j)+sigt(4,i,j)
  220   continue
        sigt(0,0,j)=sigt(1,0,j)+sigt(2,0,j)+sigt(3,0,j)+sigt(4,0,j)
  230   continue
      endif
 
C...Scale up uniformly for Donnachie-Landshoff parametrization.
      if(iproc.eq.21.or.iproc.eq.23.or.iproc.eq.25) then
        rfac=(xpar(iproc)*seps+ypar(iproc)*seta)/sigt(0,0,0)
        do 260 i1=0,n1
        do 250 i2=0,n2
        do 240 j=0,5
        sigt(i1,i2,j)=rfac*sigt(i1,i2,j)
  240   continue
  250   continue
  260   continue
      endif
 
      return
      end

c***********************************************************************

      subroutine jamxnn(srt,iz1,iz2,sig,sigel)

c...pp/pn total and elastic cross sections.

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      real*8 jamchc96,jamrgg96
      parameter (bhad=2.3d0,eps=0.0808d0,facel=0.0511d0)
      parameter(smin=2.0139999d0) ! thd0 in sub. jamxnnin


c....nn/pp collision.
      if(iz1.eq.iz2) then

        if(srt.lt.4.5d0) then
          call jamsighh(sig,1,srt)
          call jamsighh(sigel,2,srt)
c       else if(srt.le.30.0d0) then
        else
          if(iz1.eq.3) then
            plab=sqrt((srt**2-parc(24)**2-parc(24)**2)**2/
     $           (4.d0*parc(24)**2)-parc(24)**2)
          else
            plab=sqrt((srt**2-parc(25)**2-parc(25)**2)**2/
     $           (4.d0*parc(25)**2)-parc(25)**2)
          endif
          sig=jamchc96(16,plab)
          sigel=jamchc96(17,plab)
c       else
c         sig=jamrgg96(srt,1)
c         bel=4.d0*bhad+4.d0*(srt*srt)**eps-4.2d0
c         sigel=facel*sig**2/bel
        endif

c...np collision.
      else

        if(srt.lt.4.5d0) then
          call jamsighh(sig,3,srt)
          call jamsighh(sigel,4,srt)
c       else if(srt.le.30.d0) then
        else
          plab=sqrt((srt**2-parc(24)**2-parc(25)**2)**2/
     $          (4.d0*parc(24)**2)-parc(25)**2)
          sig=jamchc96(18,plab)
          sigel=jamchc96(17,plab)
c       else
c         sig=jamrgg96(srt,3)
c         bel=4.d0*bhad+4.d0*(srt*srt)**eps-4.2d0
c         sigel=facel*sig**2/bel
        endif
      end if

      if(srt.le.smin+parc(41)) then
        sigel=sig
      endif

      end


c***********************************************************************

      subroutine jamxkp(kfm,kfb,srt,emmes,embar,sig,sigel,sigch,sigy)

c....Purpose: to give anti-K-nucleon total/background elastic
c...background charge exchange and background KN->Ypi cross sections.
      implicit double precision(a-h, o-z)
      real*8 jamchc96,jamrgg96
      dimension sigy(4)
      parameter (bhad1=2.3d0,bhad2=1.4d0,eps=0.0808d0,facel=0.0511d0)

      sig=0.0d0
      sigel=0.0d0
      sigch=0.0d0
      do i=1,4
       sigy(i)=0.0d0
      end do
      if(srt.gt.emmes+embar) then
           plab=sqrt(((srt*srt-emmes*emmes-embar*embar)/(2.d0*embar) 
     &  )**2
     $              -emmes*emmes)
      else
        return
      endif

c......K- n/antiK0 p i.e. isospin =1 channel.
      if((kfm.eq.-321.and.kfb.eq.2112)
     $   .or.(kfm.eq.-311.and.kfb.eq.2212))then

c...... K- n Total/elastic
            if (plab.le.2.5d0) then
              call jamsighh(sig,13,srt)
              call jamsighh(sigel,14,srt)
              sig=0.9d0*sig
            else if(srt.lt.30.d0) then
              sig=jamchc96(14,plab)
              sigel=jamchc96(13,plab)
            else
              sig=jamrgg96(srt,13)
              bel=2*(2.3d0+0.8d0)+4.d0*(srt*srt)**0.079d0-4.2d0
              sigel=facel*sig**2/bel
            endif

          if(srt.le.10.0d0) sigel=min(100.0d0,8.5d0*plab**(-1.0d0))
          sigy(1)=min(60.0d0,2*0.617537d0*plab**(-1.98735d0)) ! k-n->Lpi
          sigy(2)=min(60.0d0,0.09d0*plab**(-2.42d0))          ! k-n->Spi
          sigy(3)=sigy(2)                               ! k-n->Spi

c.......K-p/antiK0n
      else if((kfm.eq.-321.and.kfb.eq.2212)
     $        .or.(kfm.eq.-311.and.kfb.eq.2112))then
 
            if (plab.le.3.0d0) then
              call jamsighh(sig,11,srt)
              call jamsighh(sigel,12,srt)
            else if(srt.lt.30.d0) then
              sig=jamchc96(12,plab)
              sigel=jamchc96(13,plab)
            else
              sig=jamrgg96(srt,12)
              bel=2*(2.3d0+0.8d0)+4.d0*(srt*srt)**0.079d0-4.2d0
              sigel=facel*sig**2/bel
            endif

c...Background elastic.
        if(srt.le.1.55d0) then
          sigel=min(100.0d0,9.0d0*plab**(-1.3d0))
        else if(srt.le.1.69d0) then
          sigel=6.0d0*plab**(-1.7d0)
        else if(srt.le.1.9d0) then
          sigel=13.00d0*plab**(-1.21d0)
        else if(srt.le.10.0d0) then
          sigel=9.88d0*plab**(-0.637d0)
        endif

c...Background charge exchange.
        if(srt.le.1.65d0) then
          sigch=min(100.0d0,1.0d0*plab**(-1.547d0))
c         sigch=0.6*plab**(-1.9)
        else if(srt.le.1.87d0) then
         sigch=0.6d0
        else
         sigch=0.005d0
        endif

c...t-channel hyperon productions(K+N->Lam/Sig+pi).
       sigy(1)=min(60.0d0,0.617537d0*plab**(-1.98735d0))
       sigy(2)=min(60.0d0,0.18475d0*plab**(-3.06927d0))
       sigy(3)=min(60.0d0,0.172513d0*plab**(-2.83577d0))
c      if(srt.le.1.63d0) then
c        sigy(4)=min(60.0d0,1.44898d0*plab**(-1.62197d0))
         sigy(4)=min(60.0d0,0.23d0*plab**(-3.32d0))
c      else if(srt.le.1.9d0) then
c        sigy(4)=0.6d0*plab**(-1.6d0)
c      else
c        sigy(4)=1.10938d0*plab**(-1.68098d0)
c      endif

      endif

      end

c***********************************************************************

      subroutine jamxpin(izpi,izn,srt,pr,sig,sigel)

c...Purpose: to give pi-N total and elastic cross sections.
      implicit double precision(a-h, o-z)
      real*8 jamxdelt,jamchc96,jamrgg96
      parameter(emnuc=0.938d0,empi=0.139d0)
      parameter (bhad1=2.3d0,bhad2=1.4d0,eps=0.0808d0,facel=0.0511d0)
      character chekc*80

c...Too slow. no reaction
      if(srt.lt.emnuc+empi+0.0001d0) then
        sig=0.0d0
        sigel=0.0d0
        return
      end if

c....Cross section for delta formation from fitted Breit-Wigner.
      if(srt.le.1.3d0) then
       sig=jamxdelt(izpi,izn,srt,pr)
       sigel=0.0d0
       return
      endif
 
c...pi- n/pi+ p
      if((izpi.eq.-1.and.izn.eq.0).or.(izpi.eq.1.and.izn.eq.1))then
        icha=1
c...pi+ n/pi- p
      elseif((izpi.eq.-1.and.izn.eq.1).or.(izpi.eq.1.and.izn.eq.0))then
        icha=2
c...pi0 n/pi0 p
      elseif((izpi.eq.0.and.izn.eq.0).or.(izpi.eq.0.and.izn.eq.1))then
        icha=3
      else
        write(chekc,'(i4,1x,i4)')izpi,izn
        call jamerrm(30,0,'(jamxpin:)error izpi izn='//chekc)
        return
      end if

c...Cross section from table fit.
      if(srt.le.3.d0) then
        if(icha.eq.1) then
          call jamsighh(sig,7,srt)
          call jamsighh(sigel,8,srt)
        else if(icha.eq.2) then
          call jamsighh(sig,9,srt)
          call jamsighh(sigel,10,srt)
        else
          call jamsighh(sig1,7,srt)
          call jamsighh(sig2,9,srt)
          call jamsighh(sigel1,8,srt)
          call jamsighh(sigel2,10,srt)
          sig=0.5d0*(sig1+sig2)
          sigel=0.5d0*(sigel1+sigel2)
        endif
c...Cross section from HERA fit.
      else if(srt.le.30.d0) then
        plab=((srt*srt-emnuc**2-empi**2)/(2.d0*emnuc))**2-empi**2
        if(plab.le.0.0d0) then
          sig=0.0d0
          sigel=0.0d0
          return
        else
          plab=sqrt(plab)
        endif
        if(icha.eq.1) then
          sig=jamchc96(3,plab)
          sigel=jamchc96(4,plab)
        else if(icha.eq.2) then
          sig=jamchc96(5,plab)
          sigel=jamchc96(6,plab)
        else
          sig=0.5d0*(jamchc96(5,plab)+jamchc96(3,plab))
          sigel=0.5d0*(jamchc96(6,plab)+jamchc96(4,plab))
        endif

c...Regge fit.
      else
        if(icha.eq.1) then
          sig=jamrgg96(srt,7)
        else if(icha.eq.2) then
          sig=jamrgg96(srt,8)
        else
          sig=0.5d0*(jamrgg96(srt,7)+jamrgg96(srt,8))
        endif
        bel=2.d0*bhad1+2.d0*bhad2+4.d0*(srt*srt)**eps-4.2d0
        sigel=facel*sig**2/bel
      end if

      end

c***********************************************************************

      subroutine jamxbbar(srt,plab,iz1,iz2,sig,sigel)

c...Calculate baryon-antibaryon total and elastic cross sections.
      implicit double precision(a-h, o-z)
      real*8 jamchc96,jamrgg96
      parameter (bhad=2.3d0,eps=0.0808d0,facel=0.0511d0)

      if(srt.lt.5.0d0) then
        call jamsighh(sig,5,srt)
        call jamsighh(sigel,6,srt)
      else if(srt.lt.20.d0) then
        if(iz1.eq.iz2) then
          sig=jamchc96(21,plab)
          sigel=jamchc96(22,plab)
        else
          sig=jamchc96(23,plab)
          sigel=jamchc96(22,plab)
        endif
      else
        if(iz1.eq.iz2) then
          sig=jamrgg96(srt,2)
        else
          sig=jamrgg96(srt,4)
        endif
        bel=4.d0*bhad+4.d0*(srt*srt)**eps-4.2d0
        sigel=facel*sig**2/bel
      endif

      end

c***********************************************************************

      subroutine jamxbw1(istr,srt,pr,kf1,kf2,iz1,iz2,sig,sigel,msel)

c...Purpose: to calculate Berit-Wigner cross section for M-B.
c -------------------------                                            *
c  meson-baryon absorptions for no total strangeness                   *
c --------------------------                                           *
c                                                                      *
c    pion   + n        ==> d(1232)/n*/d*                               *
c    pion   + d(1232)  ==> d(1232)/n*/d*                               *
c    pion   + n*       ==> d*                                          *
c    rho    + n        ==> n*/d*                                       *
c    eta    + n        ==> n*                                          *
c    eta'   + n        ==> n*                                          *
c    omega  + n        ==> n*                                          *
c    kaon   + lambda   ==> n*                                          *
c    kaon   + sigma    ==> n*/d*                                       *
c                                                                      *
c -------------------------                                            *
c  meson-baryon absorptions for S=-1 strangeness                       *
c --------------------------                                           *
c                                                                      *
c    Kbar   + n        ==> L*/S*                                       *
c    pion   + Sigma    ==> Lambda*/Sigma*                              *
c    .....etc.                                                         *
c                                                                      *
c iz1,iz2: charge
c                                                                    *
c--------------------------------------------------------------------*

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      real*8 jamxdelt
      parameter(maxbr=70)
      dimension ibranch(maxbr),pwid(maxbr)
      common/bwdat1/sigbwp(30)
      data idelta/0/
      data itest/0/

      izt=iz1+iz2
      isel=msel-2

c     wei1 = cofcgi(itb,itm,3,itzb,itzm,itz)**2
c     wei2 = cofcgi(itb,itm,1,itzb,itzm,itz)**2

      itry=0
10    continue
      itry=itry+1
      if(itry.ge.10) then
        write(check(1),'(''xsig sig'',g10.3,1x,g10.3)')pare(3),sig
        call jamerrm(30,1,'(jamxbw1:)No branch found m-B itry')
      endif
c============================================================
c...n + pi => d(1232)
c============================================================

c...Get N + pi => D(1232) cross section.
      sig=0.0d0
      sigel=0.0d0
      do i=1,30
      sigbwp(i)=0.0d0
      end do
      ipin=0
      if((kf1.eq.111.or.abs(kf1).eq.211)
     $.and.(kf2.eq.2112.or.kf2.eq.2212) )  ipin=1

      if(idelta.eq.1.and.ipin.eq.1) then
        if(srt.le.3.0d0) then
          sig=jamxdelt(iz1,iz2,srt,pr)
          if(izt.eq.-1) then
            kfd=1114
            sigel=sig
          else if(izt.eq.0) then
            kfd=2114
            if(iz2.eq.-1)sigel=1.d0/3.d0*sig
            if(iz2.eq.0)sigel=2.d0/3.d0*sig
          else if(izt.eq.1) then
            kfd=2214
            if(iz2.eq.1)sigel=1.d0/3.d0*sig
            if(iz2.eq.0)sigel=2.d0/3.d0*sig
          else if(izt.eq.2) then
            kfd=2224
            sigel=sig
          else
            write(check(1),'(''izt='',i3)')izt
            write(check(2),'(''kf1 kf2'',i9,1x,i9)')kf1,kf2
            call jamerrm(30,2,'(jamxbw1:)Charge invaild izt')
          endif
        endif
c....Monte Calro.
        if(isel.eq.1) then
           pare(3)=pare(3)-sig
          if(pare(3).le.0.d0) then
            kf1=kfd
            kf2=0
            return
          endif
        endif
      end if

      ispin1=max(1,mod(abs(kf1),10))
      ispin2=max(1,mod(abs(kf2),10))
      factbw=paru(1)*10.0d0*(paru(3)/pr)**2/(ispin1*ispin2)

      if(istr.eq.0) then

*============================================================
*  (1)  Delta* formation cross sections,baryon + meson => d*
*============================================================
c...offset: 1:n* 2:p*  3:d*- 4:d*0 5:d*+ 6:d*++
        if(idelta.eq.0.and.ipin.eq.1)then
          if(srt.le.1.077d0) goto 101
          if(izt.eq.-1) kfd=1114
          if(izt.eq.0)  kfd=2114
          if(izt.eq.1)  kfd=2214
          if(izt.eq.2)  kfd=2224
          kcmin=jamcomp(kfd)
          kcmax=kcmin
          isp=1
          assign 101 to label
          goto 5000
        endif
101     kcmin=mstc(24)+izt+1
        kcmax=mstc(25)+izt+1
        isp=4
        assign 100 to label
        goto 5000
100     continue
c============================================================
c (2) Get N* formation cross sections,   n + pi => n*
c============================================================
        if(izt.ge.0.and.izt.le.1) then
          kcmin=mstc(22)+izt
          kcmax=mstc(23)+izt
          isp=2
        else
          goto 7000
        endif
        assign 200 to label
        goto 5000
200     continue
        goto 7000

      else if(istr.eq.-1) then

c============================================================
c  (3)  Sigma* formation cross sections,baryon + meson => d*
c============================================================
        if(izt.lt.-1.or.izt.gt.1) return

c...Sigma(1385)
        if(izt.eq.-1) kcmin=jamcomp(3114)
        if(izt.eq. 0) kcmin=jamcomp(3214)
        if(izt.eq. 1) kcmin=jamcomp(3224)
        kcmax=kcmin
        isp=1
        assign 300 to label
        goto 5000
300     continue

c...Sigma*
        kcmin=mstc(28)+izt+1
        kcmax=mstc(29)+izt+1
        isp=3
        assign 400 to label
        goto 5000
400     continue
        if(izt.ne.0) goto 7000
c============================================================
c  (4) Get L* formation cross sections
c============================================================
        kcmin=mstc(26)
        kcmax=mstc(27)
        isp=1
        assign 500 to label
        goto 5000
500     continue
        goto 7000

      else if(istr.eq.-2) then
c============================================================
c  (5) Get Xi* formation cross sections
c============================================================
        if(izt.lt.-1.or.izt.gt.0) return
        kcmin=mstc(30)+izt+1
        kcmax=mstc(31)+izt+1
        isp=2
        assign 600 to label
        goto 5000
600     continue
        goto 7000

      else
        write(check(1),'(''istr='',i4)')istr
        call jamerrm(30,1,'(jamxbw1:)Invalid strangness istr')
      endif

5000  continue

      do ir=kcmin,kcmax,isp
        
       call jamwidm(ir,1,0,kf1,kf2,srt,ibranch,pwid,totwid,itag)

        if(itag.ne.0) then
          emr=pmas(ir,1)
          kfr=kchg(ir,4)
          ispin=max(1,mod(abs(kfr),10))
          sigbw=factbw*ispin*pwid(itag)*totwid
     &         /( (srt-emr)**2 + 0.25d0*totwid**2 )
          sigbwel=factbw*ispin*pwid(itag)**2
     &         /( (srt-emr)**2 + 0.25d0*totwid**2 )

          sigbw=min(sigbw,200.d0)
          if(isel.eq.1) then
            pare(3)=pare(3)-sigbw
            if(pare(3).le.0.0d0) then
              kf1=kfr
              kf2=0
              return
            endif 
          endif

          sig=sig+sigbw
          sigel=sigel+sigbwel

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c......Calculate Exclusive 2-body final Cross section (optional).
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(itest.eq.1) then
          ibc=0
          ibr=0
          do ib=mdcy(ir,2),mdcy(ir,2)+mdcy(ir,3)-1
            ibr=ibr+1
            kfd1=kfdp(ib,1)
            kfd2=kfdp(ib,2)
            kfd3=kfdp(ib,3)
            ibp=0

c....For pi- p reaction.
c.........Charge exchange.
            if((kfd1.eq.111.and.kfd2.eq.2112.and.kfd3.eq.0).or.
     $         (kfd2.eq.111.and.kfd1.eq.2112.and.kfd3.eq.0)) then
              ibp=1
            endif
c....For k-p reaction.
c............Charge exchange.
            if((kfd1.eq.-311.and.kfd2.eq.2112.and.kfd3.eq.0).or.
     $         (kfd2.eq.-311.and.kfd1.eq.2112.and.kfd3.eq.0)) then
              ibp=1
c...........k-p ->pi0 Lambda
            else if((kfd1.eq.111.and.kfd2.eq.3122.and.kfd3.eq.0).or.
     $         (kfd2.eq.111.and.kfd1.eq.3122.and.kfd3.eq.0)) then
              ibp=2
c...........k-p ->pi+ Sigma-
            else if((kfd1.eq.211.and.kfd2.eq.3112.and.kfd3.eq.0).or.
     $              (kfd2.eq.211.and.kfd1.eq.3112.and.kfd3.eq.0)) then
              ibp=3
c...........k-p ->pi0 Sigma0
            else if((kfd1.eq.111.and.kfd2.eq.3212.and.kfd3.eq.0).or.
     $              (kfd2.eq.111.and.kfd1.eq.3212.and.kfd3.eq.0)) then
              ibp=4
c...........k-p ->pi- Sigma+
            else if((kfd1.eq.-211.and.kfd2.eq.3222.and.kfd3.eq.0).or.
     $              (kfd2.eq.-211.and.kfd1.eq.3222.and.kfd3.eq.0)) then
              ibp=5
            endif
 
c....For k-n induced reaction.
c........k-n -> pi- lambda
            if((kfd1.eq.-211.and.kfd2.eq.3122.and.kfd3.eq.0).or.
     $         (kfd2.eq.-211.and.kfd1.eq.3122.and.kfd3.eq.0)) then
              ibp=2
c........k-n -> pi0 sigma-
            else if((kfd1.eq.111.and.kfd2.eq.3112.and.kfd3.eq.0).or.
     $              (kfd2.eq.111.and.kfd1.eq.3112.and.kfd3.eq.0)) then
              ibp=3
c........k-n -> pi- sigma0
            else if((kfd1.eq.-211.and.kfd2.eq.3212.and.kfd3.eq.0).or.
     $              (kfd2.eq.-211.and.kfd1.eq.3212.and.kfd3.eq.0)) then
              ibp=4
            endif
    
            if(ibp.ne.0)
     $      sigbwp(ibp)=sigbwp(ibp)+factbw*ispin*pwid(itag)*pwid(ibr)
     &                  /( (srt-emr)**2 + 0.25d0*totwid**2 )
 
          end do
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        endif
      end do

      goto label

7000  continue

c...No resonance formation case like S(1600)+rho.
c     if(sig.le.1d-4) return
c     if(msel.eq.3) then
c       if((pare(3).gt.0.0d0.and.srt.le.2.0d0)
c    $                       .or.(pare(3).le.0.0d0.and.kf2.ne.0)) then
c       write(check(1),8000)srt,kf1,kf2,pare(3)
c8000   format('srt=',g12.3,'kf1 kf2',i9,1x,i9,' xsig',g12.3)
c       call jamerrm(3,1,'(jamxbw1:) warning resonance formation false')
c       pare(3)=rn(0)*sig
c       goto 10
c       endif
c     endif

      end

c***********************************************************************

      subroutine jamxbw2(srt,pr,kf1,kf2,kc1,kc2,iz1,iz2
     $                                  ,sig,emrf,isel)

c...Purpose: to calculate Berit-Wigner cross section for M-M.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      parameter(maxbr=70)
      dimension ibranch(maxbr),pwid(maxbr)
      parameter(xnorm=1.0d0)

      isg=1
      sig=0.0d0
      emrf=0.0d0
      izt=iz1+iz2
      kf01=kf1
      kf02=kf2
      if(izt.le.-6.or.izt.ge.6) then
        sig=0.0d0
        return
      endif 
      if(kf01.lt.0.and.kchg(kc2,3).eq.0) then
        kf01=-kf01
        isg=-1
        izt=3
      else if(kf02.lt.0.and.kchg(kc1,3).eq.0) then
        kf02=-kf02
        isg=-1
        izt=3
      endif
      ispin1=max(1,mod(abs(kf1),10))
      ispin2=max(1,mod(abs(kf2),10))
      factbw=paru(1)*10.0d0*(paru(3)/pr)**2/(ispin1*ispin2)
      factbw=xnorm*factbw
*============================================================
      kcmin=mstc(32)
      kcmax=mstc(33)
      isp=1
      assign 100 to label
      goto 5000
100   continue
      return

5000  continue

c     irr=0
      do 10 ir=kcmin,kcmax,isp

c...No particle
       if(kchg(ir,4).eq.0) goto 10

       if(ir.eq.221.or.ir.eq.222) goto 10 ! K0_L/K0_S

c      if(ir.eq.223) goto 10 ! sigma-m
c.....rho+/rho0/sigma-m
c      if(ir.eq.121.or.ir.eq.131.or.ir.eq.223) goto 10

c...Charge desn't match
       if(kchg(ir,1).ne.izt) goto 10
       if(srt.lt.pmas(ir,1)-pmas(ir,3)) goto 10

       call jamwidm(ir,1,0,kf01,kf02,srt,ibranch,pwid,totwid,itag)

        if(itag.ne.0) then
          emr=pmas(ir,1)
          emrf=max(emrf,emr)
          kfr=kchg(ir,4)
          ispin=max(1,mod(abs(kfr),10))
          sigbw=factbw*ispin*pwid(itag)*totwid
     &         /( (srt-emr)**2 + 0.25d0*totwid**2 )
          sigbw=min(120.0d0,sigbw)
          sig=sig+sigbw
          if(isel.eq.3) then
            pare(3)=pare(3)-sigbw
            if(pare(3).le.0.0d0) then
              kf1=kchg(ir,4)*isg
              kf2=0
              goto 3000
            endif 
          endif

        endif

  10  continue

      goto label

3000  continue

      end

c***********************************************************************

      real*8 function jamxdelt(izpi,izn,srt,pr)

c...Purpose: to calculate pion plus nucleon to delta cross section.
c... ref. Gy.Wolf et al., Nucl.Phys. A517 (1990) 615
c...      A.Engel et al., Nucl.Phys. A572 (1994) 657

      implicit double precision(a-h, o-z)
      parameter(sig0=135.d0, sig1=200.d0, sig2=70.d0)
      parameter(gmr=0.11d0, be=0.3d0, be2=be*be)
      parameter(empion=0.138d0,emnuc=0.9383d0,emdelt=1.232d0)
      parameter(qr=0.227894932d0)
cccc  parameter(qr=sqrt((emdelt**2-(emnuc+empion)**2)
cccc &    *(emdelt**2-(emnuc-empion)**2))/(2*emdelt) )

      izzz=(izn+izpi+2)/2
      weight=sig0
      if(izzz.ne.1) weight=sig1
      if((izzz.eq.1).and.(izpi.ne.0)) weight=sig2
      v1=be2/(be2+pr**2)
      v2=be2/(be2+qr**2)
      gamma=((pr/qr)**3 *(emdelt/srt)*(v1/v2)**2*gmr)**2/4
      jamxdelt=weight*(qr/pr)**2*gamma/((srt-emdelt)**2+gamma)

      end


c*******************************************************************

      subroutine jambres1(srt,idd1,idd2,kf1,kf2,kc1,kc2)

c...Calculate individual resonance production probability
c...in non-strange BB collisions.
      implicit double precision(a-h, o-z)
      real*8 jambwtbl
      include 'jam2.inc'

      if(idd1.ge.1.and.idd1.le.2) then         ! N*
        kcmin1=mstc(22)+idd1-1
        kcmax1=mstc(23)+idd1-1
        isp1=2
      else if(idd1.ge.3.and.idd1.le.6) then    ! D*
        kcmin1=mstc(24)+idd1-3
        kcmax1=mstc(25)+idd1-3
        isp1=4
      else
        id1=kchg(kc1,5)
        if(id1.eq.id_nucl.or.id1.eq.id_delt) then
          isp1=1
          kcmin1=jamcomp(kf1)
          kcmax1=kcmin1
        else if(id1.eq.id_nucls) then
          isp1=1
          kcmin1=jamcomp(kf1)
          kcmax1=kcmin1
        else
          write(check(1),'(''idd1='',i3)')idd1
          call jamerrm(1,1,'(jambres1:)1 invalid idd =')
        endif
      endif

      if(idd2.ge.1.and.idd2.le.2) then         ! N*
        kcmin2=mstc(22)+idd2-1
        kcmax2=mstc(23)+idd2-1
        isp2=2
      else if(idd2.ge.3.and.idd2.le.6) then    ! D*
        kcmin2=mstc(24)+idd2-3
        kcmax2=mstc(25)+idd2-3
        isp2=4
      else
        id2=kchg(kc2,5)
        if(id2.eq.id_nucl.or.id2.eq.id_delt) then
          isp2=1
          kcmin2=jamcomp(kf2)
          kcmax2=kcmin2
        else if(id2.eq.id_nucls) then
          isp2=1
          kcmin2=jamcomp(kf2)
          kcmax2=kcmin2
        else
          write(check(1),'(''idd2='',i3)')idd2
          call jamerrm(1,1,'(jambres1:)2 invalid idd =')
        endif
      endif

      bwtot=0.0d0
      do i=kcmin1,kcmax1,isp1
      do j=kcmin2,kcmax2,isp2
        spin1=mod(kchg(i,4),10)
        spin2=mod(kchg(j,4),10)
        bwtot=bwtot+jambwtbl(srt,i,j)*spin1*spin2
      end do
      end do

      xbw=bwtot*rn(0)
      do i=kcmin1,kcmax1,isp1
      do j=kcmin2,kcmax2,isp2
        spin1=mod(kchg(i,4),10)
        spin2=mod(kchg(j,4),10)
        xbw=xbw-jambwtbl(srt,i,j)*spin1*spin2
        if(xbw.le.0.0d0) then
         kc1=i
         kc2=j
         goto 300
        end if             
      end do
      end do
      kc1=kcmax1
      kc2=kcmax2
 300  continue
      kf1=kchg(kc1,4)
      kf2=kchg(kc2,4)

      end


c***********************************************************************

      real*8 function jambwtbl(srt,kc1,kc2)

c...Give integrated value of Wreit-Bigner for non-strange resonances.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      dimension aa(7,231)

c...a 1
      data (aa(1,i),i=1,231)/
     $  0.0000d0, 15.1526d0,  7.0270d0, 14.1900d0,  6.1684d0,  5.2539d0,
     &   7.3682d0,
     $  4.2570d0,  3.5854d0,  3.1261d0,  6.2091d0,  3.6691d0,  3.0932d0,
     &   2.7471d0,
     $  2.4169d0,  4.2858d0,  2.7302d0,  2.3019d0,  2.1727d0,  1.9234d0,
     &   1.5593d0,
     $  3.9654d0,  2.5230d0,  2.1215d0,  2.0241d0,  1.7956d0,  1.4575d0,
     &   1.3562d0,
     $  3.8810d0,  2.4932d0,  2.0980d0,  2.0170d0,  1.7901d0,  1.4571d0,
     &   1.3571d0,
     $  1.3584d0,  3.7163d0,  2.4313d0,  2.0480d0,  1.9903d0,  1.7679d0,
     &   1.4441d0,
     $  1.3460d0,  1.3477d0,  1.3384d0,  3.5986d0,  2.3638d0,  1.9917d0,
     &   1.9441d0,
     $  1.7280d0,  1.4128d0,  1.3157d0,  1.3178d0,  1.3090d0,  1.2800d0,
     &   3.3815d0,
     $  2.2007d0,  1.8530d0,  1.8009d0,  1.6008d0,  1.3060d0,  1.2142d0,
     &   1.2161d0,
     $  1.2073d0,  1.1801d0,  1.0881d0,  1.2751d0,  0.8783d0,  0.7433d0,
     &   0.7659d0,
     $  0.6863d0,  0.5666d0,  0.5214d0,  0.5241d0,  0.5221d0,  0.5093d0,
     &   0.4670d0,
     $  0.1950d0, 11.2787d0,  6.0546d0,  4.6944d0,  3.8161d0,  3.3199d0,
     &   2.5416d0,
     $  2.3602d0,  2.3405d0,  2.2955d0,  2.2390d0,  2.0816d0,  0.8685d0,
     &   4.4536d0,
     $ 10.0835d0,  5.5980d0,  4.5949d0,  3.9670d0,  3.4773d0,  2.7440d0,
     &   2.5668d0,
     $  2.5552d0,  2.5205d0,  2.4650d0,  2.2872d0,  0.9881d0,  7.6614d0,
     &   9.6594d0,
     $  6.8717d0,  3.9583d0,  3.2906d0,  2.9127d0,  2.5621d0,  2.0454d0,
     &   1.9151d0,
     $  1.9091d0,  1.8864d0,  1.8463d0,  1.7122d0,  0.7461d0,  3.4337d0,
     &   3.6323d0,
     $  2.6945d0,  3.5180d0,  2.3437d0,  1.9748d0,  1.9535d0,  1.7376d0,
     &   1.4365d0,
     $  1.3458d0,  1.3489d0,  1.3423d0,  1.3152d0,  1.2135d0,  0.5373d0,
     &   2.2211d0,
     $  2.4543d0,  1.8513d0,  1.3491d0,  3.2559d0,  2.0863d0,  1.7499d0,
     &   1.7018d0,
     $  1.5134d0,  1.2396d0,  1.1544d0,  1.1567d0,  1.1491d0,  1.1245d0,
     &   1.0370d0,
     $  0.4503d0,  1.9602d0,  2.1521d0,  1.6160d0,  1.1559d0,  0.9879d0,
     &   3.2960d0,
     $  2.1938d0,  1.8490d0,  1.8298d0,  1.6284d0,  1.3447d0,  1.2582d0,
     &   1.2612d0,
     $  1.2544d0,  1.2288d0,  1.1339d0,  0.5001d0,  2.0834d0,  2.3034d0,
     &   1.7371d0,
     $  1.2632d0,  1.0810d0,  1.1821d0,  3.2169d0,  2.1677d0,  1.8283d0,
     &   1.8256d0,
     $  1.6265d0,  1.3473d0,  1.2611d0,  1.2648d0,  1.2593d0,  1.2338d0,
     &   1.1375d0,
     $  0.5030d0,  2.0685d0,  2.2982d0,  1.7356d0,  1.2693d0,  1.0850d0,
     &   1.1880d0,
     $  1.1947d0,  3.0227d0,  1.9658d0,  1.6514d0,  1.6229d0,  1.4446d0,
     &   1.1872d0,
     $  1.1060d0,  1.1087d0,  1.1021d0,  1.0785d0,  0.9946d0,  0.4324d0,
     &   1.8600d0,
     $  2.0528d0,  1.5448d0,  1.1115d0,  0.9486d0,  1.0394d0,  1.0438d0,
     &   0.9113d0,
     $  2.7213d0,  1.7839d0,  1.4979d0,  1.4854d0,  1.3237d0,  1.0897d0,
     &   1.0140d0,
     $  1.0166d0,  1.0114d0,  0.9894d0,  0.9113d0,  0.3949d0,  1.6955d0,
     &   1.8819d0,
     $  1.4169d0,  1.0226d0,  0.8700d0,  0.9553d0,  0.9600d0,  0.8360d0,
     &   0.7669d0/
c...a 2
      data (aa(2,i),i=1,231)/
     $  0.0000d0,  2.0140d0,  2.1520d0,  2.1520d0,  2.2900d0,  2.4280d0,
     &   2.1520d0,
     $  2.2900d0,  2.4280d0,  2.4280d0,  2.1520d0,  2.2900d0,  2.4280d0,
     &   2.4280d0,
     $  2.4280d0,  2.1520d0,  2.2900d0,  2.4280d0,  2.4280d0,  2.4280d0,
     &   2.4280d0,
     $  2.1520d0,  2.2900d0,  2.4280d0,  2.4280d0,  2.4280d0,  2.4280d0,
     &   2.4280d0,
     $  2.1520d0,  2.2900d0,  2.4280d0,  2.4280d0,  2.4280d0,  2.4280d0,
     &   2.4280d0,
     $  2.4280d0,  2.1520d0,  2.2900d0,  2.4280d0,  2.4280d0,  2.4280d0,
     &   2.4280d0,
     $  2.4280d0,  2.4280d0,  2.4280d0,  2.1520d0,  2.2900d0,  2.4280d0,
     &   2.4280d0,
     $  2.4280d0,  2.4280d0,  2.4280d0,  2.4280d0,  2.4280d0,  2.4280d0,
     &   2.1520d0,
     $  2.2900d0,  2.4280d0,  2.4280d0,  2.4280d0,  2.4280d0,  2.4280d0,
     &   2.4280d0,
     $  2.4280d0,  2.4280d0,  2.4280d0,  2.1520d0,  2.2900d0,  2.4280d0,
     &   2.4280d0,
     $  2.4280d0,  2.4280d0,  2.4280d0,  2.4280d0,  2.4280d0,  2.4280d0,
     &   2.4280d0,
     $  2.4280d0,  2.2900d0,  2.4280d0,  2.5660d0,  2.5660d0,  2.5660d0,
     &   2.5660d0,
     $  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,
     &   2.7040d0,
     $  2.2900d0,  2.4280d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,
     &   2.5660d0,
     $  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.7040d0,
     &   2.7040d0,
     $  2.2900d0,  2.4280d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,
     &   2.5660d0,
     $  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.7040d0,
     &   2.7040d0,
     $  2.7040d0,  2.2900d0,  2.4280d0,  2.5660d0,  2.5660d0,  2.5660d0,
     &   2.5660d0,
     $  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,
     &   2.7040d0,
     $  2.7040d0,  2.7040d0,  2.7040d0,  2.2900d0,  2.4280d0,  2.5660d0,
     &   2.5660d0,
     $  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,
     &   2.5660d0,
     $  2.5660d0,  2.7040d0,  2.7040d0,  2.7040d0,  2.7040d0,  2.7040d0,
     &   2.2900d0,
     $  2.4280d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,
     &   2.5660d0,
     $  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.7040d0,  2.7040d0,
     &   2.7040d0,
     $  2.7040d0,  2.7040d0,  2.7040d0,  2.2900d0,  2.4280d0,  2.5660d0,
     &   2.5660d0,
     $  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,
     &   2.5660d0,
     $  2.5660d0,  2.7040d0,  2.7040d0,  2.7040d0,  2.7040d0,  2.7040d0,
     &   2.7040d0,
     $  2.7040d0,  2.2900d0,  2.4280d0,  2.5660d0,  2.5660d0,  2.5660d0,
     &   2.5660d0,
     $  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,
     &   2.7040d0,
     $  2.7040d0,  2.7040d0,  2.7040d0,  2.7040d0,  2.7040d0,  2.7040d0,
     &   2.7040d0,
     $  2.2900d0,  2.4280d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,
     &   2.5660d0,
     $  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.5660d0,  2.7040d0,
     &   2.7040d0,
     $  2.7040d0,  2.7040d0,  2.7040d0,  2.7040d0,  2.7040d0,  2.7040d0,
     &   2.7040d0/
c...a 3
      data (aa(3,i),i=1,231)/
     $  0.0000d0,  2.8553d0,  2.9325d0,  2.9199d0,  2.9920d0,  3.0624d0,
     &   2.8983d0,
     $  3.0277d0,  3.1074d0,  3.1577d0,  2.8928d0,  3.0271d0,  3.1064d0,
     &   3.1568d0,
     $  3.1559d0,  2.9564d0,  3.1121d0,  3.2006d0,  3.2565d0,  3.2555d0,
     &   3.3619d0,
     $  3.0099d0,  3.1818d0,  3.2800d0,  3.3365d0,  3.3341d0,  3.4445d0,
     &   3.5348d0,
     $  3.0098d0,  3.1832d0,  3.2813d0,  3.3380d0,  3.3359d0,  3.4464d0,
     &   3.5363d0,
     $  3.5377d0,  3.0259d0,  3.1997d0,  3.2982d0,  3.3560d0,  3.3545d0,
     &   3.4659d0,
     $  3.5542d0,  3.5562d0,  3.5752d0,  3.0459d0,  3.2244d0,  3.3252d0,
     &   3.3831d0,
     $  3.3813d0,  3.4936d0,  3.5837d0,  3.5857d0,  3.6043d0,  3.6338d0,
     &   3.0562d0,
     $  3.2392d0,  3.3426d0,  3.3999d0,  3.3971d0,  3.5102d0,  3.6035d0,
     &   3.6049d0,
     $  3.6230d0,  3.6532d0,  3.6737d0,  3.3115d0,  3.5620d0,  3.7012d0,
     &   3.7626d0,
     $  3.7564d0,  3.8868d0,  4.0026d0,  4.0034d0,  4.0206d0,  4.0559d0,
     &   4.0830d0,
     $  4.5649d0,  2.9298d0,  3.0066d0,  3.1200d0,  3.1343d0,  3.1317d0,
     &   3.2291d0,
     $  3.3126d0,  3.3140d0,  3.3316d0,  3.3596d0,  3.3777d0,  3.7545d0,
     &   3.1943d0,
     $  2.8716d0,  2.9905d0,  3.0836d0,  3.1209d0,  3.1194d0,  3.2195d0,
     &   3.3008d0,
     $  3.3022d0,  3.3209d0,  3.3484d0,  3.3651d0,  3.7371d0,  1.6269d0,
     &   1.6142d0,
     $  2.8882d0,  3.0261d0,  3.1175d0,  3.1668d0,  3.1649d0,  3.2694d0,
     &   3.3561d0,
     $  3.3574d0,  3.3759d0,  3.4044d0,  3.4226d0,  3.8088d0,  3.1662d0,
     &   3.1389d0,
     $  3.1818d0,  3.0541d0,  3.2472d0,  3.3548d0,  3.4162d0,  3.4142d0,
     &   3.5313d0,
     $  3.6287d0,  3.6302d0,  3.6495d0,  3.6805d0,  3.7013d0,  4.1300d0,
     &   3.3913d0,
     $  3.3805d0,  3.4382d0,  3.7309d0,  3.0900d0,  3.2929d0,  3.4077d0,
     &   3.4682d0,
     $  3.4641d0,  3.5838d0,  3.6884d0,  3.6891d0,  3.7071d0,  3.7392d0,
     &   3.7629d0,
     $  4.2069d0,  3.4473d0,  3.4334d0,  3.4954d0,  3.7941d0,  3.8644d0,
     &   3.0788d0,
     $  3.2756d0,  3.3846d0,  3.4454d0,  3.4426d0,  3.5605d0,  3.6597d0,
     &   3.6611d0,
     $  3.6805d0,  3.7117d0,  3.7327d0,  4.1647d0,  3.4211d0,  3.4095d0,
     &   3.4682d0,
     $  3.7618d0,  3.8269d0,  3.7933d0,  3.0915d0,  3.2904d0,  3.4005d0,
     &   3.4628d0,
     $  3.4600d0,  3.5790d0,  3.6775d0,  3.6791d0,  3.6985d0,  3.7297d0,
     &   3.7508d0,
     $  4.1844d0,  3.4381d0,  3.4274d0,  3.4862d0,  3.7825d0,  3.8460d0,
     &   3.8134d0,
     $  3.8340d0,  3.1089d0,  3.3179d0,  3.4349d0,  3.4957d0,  3.4914d0,
     &   3.6124d0,
     $  3.7182d0,  3.7189d0,  3.7368d0,  3.7695d0,  3.7930d0,  4.2422d0,
     &   3.4753d0,
     $  3.4610d0,  3.5239d0,  3.8257d0,  3.8965d0,  3.8583d0,  3.8777d0,
     &   3.9291d0,
     $  3.1392d0,  3.3552d0,  3.4774d0,  3.5389d0,  3.5345d0,  3.6579d0,
     &   3.7653d0,
     $  3.7667d0,  3.7846d0,  3.8179d0,  3.8424d0,  4.3009d0,  3.5210d0,
     &   3.5064d0,
     $  3.5718d0,  3.8790d0,  3.9512d0,  3.9130d0,  3.9322d0,  3.9845d0,
     &   4.0404d0/
c...a 4
      data (aa(4,i),i=1,231)/
     $  0.0000d0,  0.1071d0,  0.2647d0,  0.1083d0,  0.2929d0,  0.3246d0,
     &   0.1923d0,
     $  0.3107d0,  0.3336d0,  0.3026d0,  0.1921d0,  0.3001d0,  0.3218d0,
     &   0.2913d0,
     $  0.2805d0,  0.2287d0,  0.3125d0,  0.3293d0,  0.2895d0,  0.2787d0,
     &   0.2729d0,
     $  0.2673d0,  0.3474d0,  0.3610d0,  0.3138d0,  0.3025d0,  0.2934d0,
     &   0.3121d0,
     $  0.2653d0,  0.3414d0,  0.3547d0,  0.3082d0,  0.2971d0,  0.2881d0,
     &   0.3066d0,
     $  0.3011d0,  0.2580d0,  0.3289d0,  0.3417d0,  0.2970d0,  0.2862d0,
     &   0.2777d0,
     $  0.2961d0,  0.2908d0,  0.2807d0,  0.2655d0,  0.3335d0,  0.3453d0,
     &   0.2998d0,
     $  0.2890d0,  0.2800d0,  0.2978d0,  0.2925d0,  0.2825d0,  0.2843d0,
     &   0.2811d0,
     $  0.3499d0,  0.3609d0,  0.3128d0,  0.3018d0,  0.2916d0,  0.3089d0,
     &   0.3035d0,
     $  0.2933d0,  0.2948d0,  0.3054d0,  0.3513d0,  0.3717d0,  0.3692d0,
     &   0.3193d0,
     $  0.3098d0,  0.2949d0,  0.3053d0,  0.3005d0,  0.2918d0,  0.2921d0,
     &   0.3001d0,
     $  0.2846d0,  0.1366d0,  0.2883d0,  0.3082d0,  0.3094d0,  0.2990d0,
     &   0.3035d0,
     $  0.3305d0,  0.3249d0,  0.3132d0,  0.3162d0,  0.3299d0,  0.3351d0,
     &   0.2814d0,
     $  0.1586d0,  0.2801d0,  0.2994d0,  0.2825d0,  0.2725d0,  0.2734d0,
     &   0.2975d0,
     $  0.2924d0,  0.2817d0,  0.2846d0,  0.2971d0,  0.3046d0,  2.9577d0,
     &   2.7639d0,
     $  0.1964d0,  0.3122d0,  0.3309d0,  0.2994d0,  0.2885d0,  0.2850d0,
     &   0.3065d0,
     $  0.3012d0,  0.2904d0,  0.2927d0,  0.3047d0,  0.3062d0,  0.3008d0,
     &   0.2782d0,
     $  0.2947d0,  0.2626d0,  0.3212d0,  0.3306d0,  0.2869d0,  0.2769d0,
     &   0.2676d0,
     $  0.2829d0,  0.2781d0,  0.2688d0,  0.2702d0,  0.2798d0,  0.2745d0,
     &   0.3031d0,
     $  0.2732d0,  0.2804d0,  0.2568d0,  0.2945d0,  0.3590d0,  0.3660d0,
     &   0.3171d0,
     $  0.3065d0,  0.2948d0,  0.3093d0,  0.3042d0,  0.2943d0,  0.2954d0,
     &   0.3050d0,
     $  0.2941d0,  0.3338d0,  0.3017d0,  0.3081d0,  0.2796d0,  0.3021d0,
     &   0.2678d0,
     $  0.3253d0,  0.3342d0,  0.2904d0,  0.2804d0,  0.2708d0,  0.2858d0,
     &   0.2810d0,
     $  0.2716d0,  0.2729d0,  0.2825d0,  0.2764d0,  0.3065d0,  0.2766d0,
     &   0.2835d0,
     $  0.2595d0,  0.2819d0,  0.2621d0,  0.2630d0,  0.3167d0,  0.3252d0,
     &   0.2823d0,
     $  0.2726d0,  0.2632d0,  0.2781d0,  0.2734d0,  0.2643d0,  0.2656d0,
     &   0.2749d0,
     $  0.2696d0,  0.2982d0,  0.2689d0,  0.2756d0,  0.2523d0,  0.2745d0,
     &   0.2550d0,
     $  0.2480d0,  0.2994d0,  0.3569d0,  0.3628d0,  0.3141d0,  0.3038d0,
     &   0.2918d0,
     $  0.3056d0,  0.3007d0,  0.2910d0,  0.2920d0,  0.3014d0,  0.2901d0,
     &   0.3310d0,
     $  0.2993d0,  0.3050d0,  0.2763d0,  0.2983d0,  0.2786d0,  0.2713d0,
     &   0.2946d0,
     $  0.3051d0,  0.3575d0,  0.3616d0,  0.3129d0,  0.3027d0,  0.2903d0,
     &   0.3036d0,
     $  0.2986d0,  0.2892d0,  0.2900d0,  0.2990d0,  0.2868d0,  0.3292d0,
     &   0.2979d0,
     $  0.3028d0,  0.2738d0,  0.2951d0,  0.2758d0,  0.2687d0,  0.2914d0,
     &   0.2882d0/
c...a 5
      data (aa(5,i),i=1,231)/
     $  0.0000d0,  2.0415d0,  2.1542d0,  2.1203d0,  2.2149d0,  2.2873d0,
     &   2.2964d0,
     $  2.5584d0,  2.6795d0,  2.9757d0,  2.3318d0,  2.6030d0,  2.7269d0,
     &   3.0084d0,
     $  3.0396d0,  2.4511d0,  2.7714d0,  2.9151d0,  3.1690d0,  3.1950d0,
     &   3.3441d0,
     $  2.4533d0,  2.8100d0,  2.9648d0,  3.2233d0,  3.2480d0,  3.4021d0,
     &   3.4672d0,
     $  2.4787d0,  2.8371d0,  2.9930d0,  3.2444d0,  3.2684d0,  3.4200d0,
     &   3.4848d0,
     $  3.5020d0,  2.5224d0,  2.8767d0,  3.0340d0,  3.2744d0,  3.2975d0,
     &   3.4449d0,
     $  3.5083d0,  3.5251d0,  3.5472d0,  2.5360d0,  2.8979d0,  3.0578d0,
     &   3.2963d0,
     $  3.3188d0,  3.4662d0,  3.5309d0,  3.5474d0,  3.5691d0,  3.5911d0,
     &   2.5132d0,
     $  2.8866d0,  3.0482d0,  3.2942d0,  3.3169d0,  3.4681d0,  3.5351d0,
     &   3.5518d0,
     $  3.5738d0,  3.5964d0,  3.6019d0,  2.7955d0,  3.2446d0,  3.4327d0,
     &   3.6339d0,
     $  3.6483d0,  3.7896d0,  3.8660d0,  3.8791d0,  3.8950d0,  3.9180d0,
     &   3.9289d0,
     $  4.2425d0,  2.3150d0,  2.4401d0,  2.6038d0,  2.9021d0,  2.9431d0,
     &   3.1228d0,
     $  3.1757d0,  3.2017d0,  3.2395d0,  3.2631d0,  3.2561d0,  3.6286d0,
     &   2.9149d0,
     $  2.3937d0,  2.6233d0,  2.7595d0,  3.0457d0,  3.0800d0,  3.2425d0,
     &   3.2946d0,
     $  3.3165d0,  3.3485d0,  3.3703d0,  3.3667d0,  3.7094d0,  1.8735d0,
     &   1.7146d0,
     $  2.4112d0,  2.6929d0,  2.8296d0,  3.1276d0,  3.1597d0,  3.3246d0,
     &   3.3832d0,
     $  3.4042d0,  3.4341d0,  3.4566d0,  3.4553d0,  3.7988d0,  3.0794d0,
     &   3.2059d0,
     $  3.2850d0,  2.6967d0,  3.0752d0,  3.2393d0,  3.4734d0,  3.4949d0,
     &   3.6411d0,
     $  3.7073d0,  3.7231d0,  3.7439d0,  3.7660d0,  3.7720d0,  4.0894d0,
     &   3.4426d0,
     $  3.5470d0,  3.6325d0,  3.9397d0,  2.6060d0,  3.0094d0,  3.1779d0,
     &   3.4351d0,
     $  3.4577d0,  3.6148d0,  3.6861d0,  3.7029d0,  3.7254d0,  3.7485d0,
     &   3.7546d0,
     $  4.0872d0,  3.3900d0,  3.5045d0,  3.5964d0,  3.9238d0,  3.9078d0,
     &   2.6986d0,
     $  3.0817d0,  3.2468d0,  3.4814d0,  3.5025d0,  3.6498d0,  3.7168d0,
     &   3.7327d0,
     $  3.7539d0,  3.7760d0,  3.7819d0,  4.1005d0,  3.4502d0,  3.5546d0,
     &   3.6407d0,
     $  3.9491d0,  3.9339d0,  3.9589d0,  2.7380d0,  3.1200d0,  3.2869d0,
     &   3.5125d0,
     $  3.5324d0,  3.6763d0,  3.7427d0,  3.7580d0,  3.7780d0,  3.7999d0,
     &   3.8066d0,
     $  4.1206d0,  3.4876d0,  3.5870d0,  3.6718d0,  3.9733d0,  3.9593d0,
     &   3.9829d0,
     $  4.0061d0,  2.6478d0,  3.0582d0,  3.2291d0,  3.4775d0,  3.4989d0,
     &   3.6533d0,
     $  3.7250d0,  3.7412d0,  3.7627d0,  3.7858d0,  3.7922d0,  4.1212d0,
     &   3.4382d0,
     $  3.5477d0,  3.6389d0,  3.9603d0,  3.9460d0,  3.9703d0,  3.9950d0,
     &   3.9835d0,
     $  2.6782d0,  3.0957d0,  3.2709d0,  3.5119d0,  3.5323d0,  3.6851d0,
     &   3.7572d0,
     $  3.7734d0,  3.7939d0,  3.8171d0,  3.8246d0,  4.1516d0,  3.4787d0,
     &   3.5834d0,
     $  3.6749d0,  3.9915d0,  3.9796d0,  4.0023d0,  4.0260d0,  4.0166d0,
     &   4.0489d0/
c...a 6
      data (aa(6,i),i=1,231)/
     $  0.0000d0, -0.5179d0, -1.0146d0, -0.8127d0, -1.3487d0, -1.6334d0,
     &  -0.5249d0,
     $ -1.0171d0, -1.3478d0, -1.0181d0, -0.3906d0, -0.8002d0, -1.0901d0,
     &  -0.8022d0,
     $ -0.6265d0, -0.3974d0, -0.8113d0, -1.1036d0, -0.8127d0, -0.6349d0,
     &  -0.6427d0,
     $ -0.6672d0, -1.2027d0, -1.5345d0, -1.2020d0, -0.9574d0, -0.9688d0,
     &  -1.3971d0,
     $ -0.6183d0, -1.1396d0, -1.4723d0, -1.1391d0, -0.9039d0, -0.9147d0,
     &  -1.3306d0,
     $ -1.2652d0, -0.5701d0, -1.0748d0, -1.4059d0, -1.0746d0, -0.8496d0,
     &  -0.8598d0,
     $ -1.2616d0, -1.1976d0, -1.1316d0, -0.6535d0, -1.1926d0, -1.5324d0,
     &  -1.1918d0,
     $ -0.9473d0, -0.9586d0, -1.3890d0, -1.3217d0, -1.2520d0, -1.3805d0,
     &  -0.7509d0,
     $ -1.3110d0, -1.6418d0, -1.3095d0, -1.0486d0, -1.0609d0, -1.5106d0,
     &  -1.4424d0,
     $ -1.3710d0, -1.5042d0, -1.6274d0, -1.0303d0, -1.6108d0, -1.8755d0,
     &  -1.6071d0,
     $ -1.3131d0, -1.3279d0, -1.8022d0, -1.7380d0, -1.6677d0, -1.8047d0,
     &  -1.9131d0,
     $ -2.1165d0, -0.8466d0, -1.3945d0, -1.6820d0, -1.3927d0, -1.1274d0,
     &  -1.1428d0,
     $ -1.5830d0, -1.5195d0, -1.4517d0, -1.5812d0, -1.6923d0, -1.9276d0,
     &  -1.7314d0,
     $ -0.5554d0, -1.0481d0, -1.3708d0, -1.0483d0, -0.8288d0, -0.8412d0,
     &  -1.2308d0,
     $ -1.1683d0, -1.1041d0, -1.2215d0, -1.3374d0, -1.6262d0, -1.4157d0,
     &  -1.0772d0,
     $ -0.5270d0, -0.9710d0, -1.2545d0, -0.9706d0, -0.7702d0, -0.7815d0,
     &  -1.1339d0,
     $ -1.0782d0, -1.0205d0, -1.1264d0, -1.2291d0, -1.4809d0, -1.2947d0,
     &  -0.9956d0,
     $ -0.9188d0, -0.6469d0, -1.1533d0, -1.4641d0, -1.1514d0, -0.9182d0,
     &  -0.9310d0,
     $ -1.3355d0, -1.2727d0, -1.2073d0, -1.3281d0, -1.4426d0, -1.7161d0,
     &  -1.5096d0,
     $ -1.1780d0, -1.0845d0, -1.2757d0, -0.9617d0, -1.5480d0, -1.8401d0,
     &  -1.5448d0,
     $ -1.2552d0, -1.2716d0, -1.7461d0, -1.6791d0, -1.6069d0, -1.7456d0,
     &  -1.8612d0,
     $ -2.0963d0, -1.8928d0, -1.5672d0, -1.4307d0, -1.6639d0, -2.0653d0,
     &  -0.7561d0,
     $ -1.2841d0, -1.5827d0, -1.2814d0, -1.0307d0, -1.0446d0, -1.4686d0,
     &  -1.4052d0,
     $ -1.3384d0, -1.4640d0, -1.5770d0, -1.8311d0, -1.6300d0, -1.3058d0,
     &  -1.1974d0,
     $ -1.4009d0, -1.7878d0, -1.5239d0, -0.7119d0, -1.2313d0, -1.5351d0,
     &  -1.2288d0,
     $ -0.9852d0, -0.9985d0, -1.4147d0, -1.3515d0, -1.2853d0, -1.4091d0,
     &  -1.5227d0,
     $ -1.7849d0, -1.5817d0, -1.2541d0, -1.1516d0, -1.3502d0, -1.7379d0,
     &  -1.4742d0,
     $ -1.4241d0, -0.9666d0, -1.5560d0, -1.8501d0, -1.5527d0, -1.2615d0,
     &  -1.2780d0,
     $ -1.7552d0, -1.6877d0, -1.6152d0, -1.7546d0, -1.8709d0, -2.1078d0,
     &  -1.9030d0,
     $ -1.5752d0, -1.4381d0, -1.6725d0, -2.0765d0, -1.7972d0, -1.7470d0,
     &  -2.0877d0,
     $ -0.9820d0, -1.5680d0, -1.8539d0, -1.5646d0, -1.2731d0, -1.2896d0,
     &  -1.7648d0,
     $ -1.6983d0, -1.6264d0, -1.7650d0, -1.8791d0, -2.1070d0, -1.9064d0,
     &  -1.5861d0,
     $ -1.4471d0, -1.6814d0, -2.0785d0, -1.8033d0, -1.7542d0, -2.0898d0,
     &  -2.0912d0/
c...a 7
      data (aa(7,i),i=1,231)/
     $  0.0000d0,  0.6001d0,  0.7120d0,  0.5636d0,  0.6696d0,  0.6294d0,
     &   0.5821d0,
     $  0.6911d0,  0.6505d0,  0.6708d0,  0.5002d0,  0.5937d0,  0.5589d0,
     &   0.5762d0,
     $  0.4949d0,  0.5009d0,  0.5947d0,  0.5600d0,  0.5772d0,  0.4958d0,
     &   0.4965d0,
     $  0.5964d0,  0.7087d0,  0.6672d0,  0.6881d0,  0.5909d0,  0.5920d0,
     &   0.7060d0,
     $  0.5903d0,  0.7014d0,  0.6604d0,  0.6809d0,  0.5848d0,  0.5858d0,
     &   0.6986d0,
     $  0.6913d0,  0.5810d0,  0.6902d0,  0.6499d0,  0.6700d0,  0.5754d0,
     &   0.5764d0,
     $  0.6874d0,  0.6802d0,  0.6692d0,  0.6071d0,  0.7214d0,  0.6791d0,
     &   0.7003d0,
     $  0.6015d0,  0.6025d0,  0.7185d0,  0.7110d0,  0.6995d0,  0.7312d0,
     &   0.6079d0,
     $  0.7225d0,  0.6801d0,  0.7016d0,  0.6026d0,  0.6037d0,  0.7199d0,
     &   0.7124d0,
     $  0.7010d0,  0.7327d0,  0.7341d0,  0.5455d0,  0.6493d0,  0.6111d0,
     &   0.6313d0,
     $  0.5422d0,  0.5434d0,  0.6481d0,  0.6414d0,  0.6312d0,  0.6596d0,
     &   0.6607d0,
     $  0.5944d0,  0.5675d0,  0.6744d0,  0.6342d0,  0.6552d0,  0.5629d0,
     &   0.5642d0,
     $  0.6722d0,  0.6653d0,  0.6547d0,  0.6842d0,  0.6853d0,  0.6159d0,
     &   0.6391d0,
     $  0.5695d0,  0.6765d0,  0.6369d0,  0.6567d0,  0.5640d0,  0.5652d0,
     &   0.6736d0,
     $  0.6666d0,  0.6559d0,  0.6856d0,  0.6869d0,  0.6183d0,  0.6416d0,
     &   0.6428d0,
     $  0.5022d0,  0.5967d0,  0.5619d0,  0.5793d0,  0.4975d0,  0.4986d0,
     &   0.5944d0,
     $  0.5881d0,  0.5787d0,  0.6049d0,  0.6061d0,  0.5458d0,  0.5661d0,
     &   0.5671d0,
     $  0.5004d0,  0.5521d0,  0.6564d0,  0.6185d0,  0.6373d0,  0.5473d0,
     &   0.5485d0,
     $  0.6541d0,  0.6472d0,  0.6368d0,  0.6657d0,  0.6671d0,  0.6011d0,
     &   0.6231d0,
     $  0.6241d0,  0.5507d0,  0.6061d0,  0.5758d0,  0.6851d0,  0.6449d0,
     &   0.6658d0,
     $  0.5718d0,  0.5733d0,  0.6835d0,  0.6764d0,  0.6655d0,  0.6956d0,
     &   0.6968d0,
     $  0.6271d0,  0.6499d0,  0.6521d0,  0.5755d0,  0.6337d0,  0.6615d0,
     &   0.5557d0,
     $  0.6607d0,  0.6221d0,  0.6417d0,  0.5511d0,  0.5524d0,  0.6586d0,
     &   0.6517d0,
     $  0.6412d0,  0.6702d0,  0.6715d0,  0.6046d0,  0.6269d0,  0.6284d0,
     &   0.5544d0,
     $  0.6104d0,  0.6376d0,  0.6144d0,  0.5553d0,  0.6602d0,  0.6218d0,
     &   0.6411d0,
     $  0.5506d0,  0.5518d0,  0.6580d0,  0.6511d0,  0.6406d0,  0.6696d0,
     &   0.6710d0,
     $  0.6043d0,  0.6265d0,  0.6278d0,  0.5539d0,  0.6097d0,  0.6372d0,
     &   0.6138d0,
     $  0.6133d0,  0.5783d0,  0.6881d0,  0.6478d0,  0.6687d0,  0.5743d0,
     &   0.5758d0,
     $  0.6865d0,  0.6794d0,  0.6684d0,  0.6986d0,  0.6999d0,  0.6300d0,
     &   0.6529d0,
     $  0.6549d0,  0.5780d0,  0.6365d0,  0.6645d0,  0.6404d0,  0.6400d0,
     &   0.6675d0,
     $  0.5675d0,  0.6754d0,  0.6358d0,  0.6564d0,  0.5638d0,  0.5652d0,
     &   0.6739d0,
     $  0.6669d0,  0.6562d0,  0.6858d0,  0.6871d0,  0.6184d0,  0.6408d0,
     &   0.6429d0,
     $  0.5674d0,  0.6249d0,  0.6523d0,  0.6287d0,  0.6283d0,  0.6552d0,
     &   0.6432d0/


      jd1=0
      jd2=0
      id1=kchg(kc1,5)
      id2=kchg(kc2,5)
      if(id1.eq.id_nucl)  jd1=1
      if(id1.eq.id_delt)  jd1=2
      if(id1.eq.id_nucls) jd1=3+(kc1-mstc(22))/2
      if(id1.eq.id_delts) jd1=4+(mstc(23)-mstc(22))/2+(kc1-mstc(24))/4
      if(id2.eq.id_nucl)  jd2=1
      if(id2.eq.id_delt)  jd2=2
      if(id2.eq.id_nucls) jd2=3+(kc2-mstc(22))/2
      if(id2.eq.id_delts) jd2=4+(mstc(23)-mstc(22))/2+(kc2-mstc(24))/4
      if(jd1.eq.0.or.jd2.eq.0) then
        write(check(1),'(''kc='',i4,1x,i4)')kc1,kc2
        write(check(2),'(''kf='',i9,1x,i9)')kchg(kc1,4),kchg(kc2,4)
        write(check(3),'(''id='',i4,1x,i4)')id1,id2
        call jamerrm(30,3,'(jambwtbl:)jd1,jd2=0')
      endif
      imin=min(jd1,jd2)
      imax=max(jd1,jd2)
      i=(imax*(imax-1))/2+imin
      if(i.le.1.or.i.gt.231) then
        write(check(1),'(''kc='',i4,1x,i4)')kc1,kc2
        write(check(2),'(''kf='',i9,1x,i9)')kchg(kc1,4),kchg(kc2,4)
        write(check(3),'(''id='',i4,1x,i4)')id1,id2
        write(check(4),'(''jd i='',i4,1x,i4)')jd1,jd2,i
        call jamerrm(30,4,'(jambwtbl:)invalid i')
      endif

c...f(x)=(x< 7.0)? a0*(x/a1-1)**a2*a3/((x/a4-1)**2+a3**2):b0+b1*x    
      jambwtbl=0.0d0
      if(srt.le.aa(2,i)) then
        jambwtbl=0.0d0
      else if(srt.le.7.0d0) then
        jambwtbl=aa(1,i)*(srt/aa(2,i)-1)**aa(3,i)*aa(4,i)
     $              /( (srt/aa(5,i)-1)**2 + aa(4,i)**2 )
      else
        jambwtbl=aa(6,i)+aa(7,i)*srt
      endif

      end

c***********************************************************************

      subroutine jambres2(srt,idd,emax,kfr,kcr)

c...Purpose: to calculate non-strange baryonic resonance prob.
c...idd     : 1:n* 2:p*  3:d*- 4:d*0 5:d*+ 6:d*++
c...in      : n* type

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      parameter(maxbr=70)
      dimension ibranch(maxbr),pwid(maxbr)
      dimension bw(20),kcb(20)

      if(idd.le.2) then         ! N*
        kcmin=mstc(22)+idd-1
        kcma=mstc(23)+idd-1
        isp=2
      else if(idd.le.6) then    ! D*
        kcmin=mstc(24)+idd-3
        kcma=mstc(25)+idd-3
        isp=4
      else
        write(check(1),'(''idd emax='',i4,g10.3)')idd,emax
        call jamerrm(1,1,'(jambres2:)invalid idd')
      endif

      kcmax=kcma
      do kc=kcmin,kcma,isp
       if(emax.lt.pmas(kc,1)-pmas(kc,3)) then
         kcmax=kc-1
         goto 100
       endif
      end do
100   continue
      if(kcmax.lt.kcmin) then
        call jamerrm(30,0,'(jambres2)invalid')
      endif
      kcmax=kcmax
      

c....Only spin dependent
      if(mstc(64).eq.1) then
        spin=0.0d0
        do i=kcmin,kcmax,isp
         kf=kchg(i,4)
         spin=spin+max(1,mod(kf,10))
        end do
        xran=rn(0)*spin
        do i=kcmin,kcmax,isp
          kf=kchg(i,4)
          xran=xran-max(1,mod(kf,10))
          if(xran.le.0.0d0) then
            kfr=kf
            kcr=i
            goto 200
          endif
        end do
        kfr=kchg(kcmax,4)
        kcr=kcmax
200     continue

c...Breit-Wigner
      else
        bwtot=0.0d0
        j=0
        do i=kcmin,kcmax,isp
         j=j+1
         em0=pmas(i,1)
         gam=pmas(i,2)
         if(mstc(64).eq.3)
     $     call jamwidm(i,1,0,0,0,srt,ibranch,pwid,gam,itag)
         kcb(j)=i
         spin=max(1,mod(kchg(i,4),10))
         bw(j)=spin*gam/((srt-em0)**2+0.25d0*gam**2)
         bwtot=bwtot+bw(j)
        end do
   
        xran=rn(0)*bwtot
        do i=1,j
          xran=xran-bw(i)
          if(xran.le.0.0d0) then
            kcr=kcb(i)
            goto 300
          endif
        end do 
        kcr=kcb(j)
 300    continue
        kfr=kchg(kcr,4)
      endif
 
      end

c***********************************************************************

      subroutine jamdetb1(msel,kc,srt,em1,em2,pr,fac,nhlf)

c...Purpose to calculate the correction factor of BR->BN cross section.
c...Integration by using Simpson (1/3) rule.
c...kc : KC particle code               (input)
c...srt: c.m. energy of the reaction.
c...em1: ingoing resonance mass.
c...em2: ingoing particle mass (not resonance assumed).
c---------------------------------------------------------------------
c...isw1=1: non-rel
c...isw1=2: rel
c...isw1=3: rel2
c...isw2=1: no momentum dependence for integrand.
c...isw2=2: momentum dependence
c...isw2=3: momentum dependence (P.Danielewicz,G.F.Bertch)
c...isw2=4: momentum-squeard dependence
c---------------------------------------------------------------------
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      real*8 jambwf1
      parameter(emnuc=0.939d0,empion=0.138d0)
      data mxhlf/30/     ! max. number of bin for integration
      data eps/ 1.0d-3/  ! accuracy

      isw=mstc(62)
      isw1=mod(isw/10,10)
      isw2=mod(isw,10)
      if((isw1.lt.1.or.isw1.gt.3).or.(isw2.lt.1.or.isw2.gt.4)) then
        write(check(1),'(''isw='',i3)')isw
        call jamerrm(30,1,'(jamdetb1:)invalid isw')
      endif

      wid=pmas(kc,2)     ! decay width
      if(wid.le.1d-5) then
        fac=1.0d0
        return
      endif
      emr0=pmas(kc,1)               ! pole mass
      emin=pmas(kc,1)-pmas(kc,3)    ! min. mass
      emax=srt-em2-0.001d0
      id=kchg(kc,5)
      if(id.eq.id_nucls) then
        emin=emnuc+2*empion+parc(41)
      else if(id.eq.id_delt) then
        emin=emnuc+empion+parc(41)
      else if(id.eq.id_delts) then
        emin=emnuc+3*empion+parc(41)
      endif

      if(emax.le.emin) then
        if(msel.eq.1) fac=1.0d0/(pr*pr)
        if(msel.eq.2) fac=0.0d0
        return
      endif


c...Option for constant width and no final momentum dependence.
      if(isw.le.1) then
        den1=atan((emax*emax-emr0*emr0)/(wid*emr0))
        den2=atan((emin*emin-emr0*emr0)/(wid*emr0))
        fac=paru(1)/max(den1-den2,1.d-4)/(pr*pr)
        return
      endif

c....Initialization for integration.
      a=emin
      b=emax
      jmax=0
      nhlf=0
      h=b-a
      s1=jambwf1(isw,kc,a,emr0,srt,em2)+jambwf1(isw,kc,b,emr0,srt,em2)
      s2=0.0d0
      s4=0.0d0
      sum=0.5d0*h*s1

c...Start the Simpson method.
100   continue
      jmax=2*jmax+1
      nhlf=nhlf+1
      h=h*0.5d0
      s2=s2+s4
      s4=0.0d0
      do 110 j=1,jmax,2
       s4=s4+jambwf1(isw,kc,a+j*h,emr0,srt,em2)
110   continue
      ss=sum
      sum=(h/3.0d0)*(s1+2.0d0*s2+4.0d0*s4)
      ds=abs(sum-ss)
      if(ds.le.eps*abs(sum)) goto 200
      if(nhlf.eq.mxhlf) then
        write(6,999)'simpson method does not converge. mxhlf=',mxhlf
        nhlf=-nhlf
        goto 200
      endif
      goto 100

200   continue
      fac=sum/paru(1)
      if(msel.eq.2) return

      fac=1.d0/max(sum/paru(1),1d-5)
c     fac=paru(1)/sum
      if(isw2.eq.1) then
        fac=fac/(pr*pr)
      else if(isw2.eq.2.or.isw2.eq.3) then
        fac=fac/pr
        if(isw2.eq.3) fac=fac*(em1/emr0)
      endif

      return
999   format(a,i3,2(a,1pe22.15))

      end

c***********************************************************************

      real*8 function jambwf1(idetsw,kc,em,emr0,srt,em2)

c...Provide B-W function for jamdetb1.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      parameter(maxbr=70)
      dimension ibranch(maxbr),pwid(maxbr)
c...Functions: momentum in two-particle cm.
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a) 
c...Options.
c...isw:
c...    =1  : no final-momentum dependence for integrand.
c...    =2,3: final-momentum dependence.
c...    =4  : final-momentum-squeard dependence.

      isw=mod(idetsw,10)
      if(idetsw.ge.10) then
        call jamwidm(kc,1,0,0,0,em,ibranch,pwid,gam,itag)
      else
        gam=pmas(kc,2)
      endif

c...Non-rel
      if(idetsw.le.20) then
        jambwf1=gam*0.5d0/((em-emr0)**2+0.25d0*gam**2)
c...Rel1
      else if(idetsw.le.30) then
        jambwf1=2*em*emr0*gam/((em**2-emr0**2)**2+(emr0*gam)**2)
c...Rel2
      else if(idetsw.le.40) then
        jambwf1=2*em**2*gam/((em*em-emr0**2)**2+(em*gam)**2)
      endif

      if(isw.eq.2.or.isw.eq.3) then
        if(srt.ge.em+em2) then
          pf=sqrt((srt**2-(em+em2)**2)*(srt**2-(em-em2)**2))/(2.d0*srt)
        else
          write(check(1),8000)srt,em,em2,chaf(kc,1)
          call jamerrm(1,1,'(jambwf1:)srt<em1+em2')
          pf=0.0d0
        endif
        jambwf1=jambwf1*pf
      else if(isw.eq.4) then
        jambwf1=jambwf1*pawt(srt,em,em2)**2
      endif
 
 8000 format('srt em em2=',3(g12.3,1x),a8)
      end

c***********************************************************************

      subroutine jamdetb2(msel,kc1,kc2,srt0,em1,em2,pr,fac,nhlf)

c...Purpose to calculate the correction factor of RR->NN cross section.
c...Integration by using Simpson (1/3) rule.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      common/dbdat/srt,em01,em02,idbsw
      real*8 jambwf2
      parameter(emnuc=0.939d0,empion=0.138d0)
      parameter(mxhlf=30,eps=1.0d-3)

c...Replace integral to constant width for high energy.
      if(msel.eq.1.and.srt.ge.3.5d0) then
        call jamdetb3(kc1,kc2,srt0,em1,em2,pr,fac,nhlf)
        fac=fac/(pr*pr)
        return
      endif

c...Check options.
      isw=mstc(62)
      isw1=mod(isw/10,10)
      isw2=mod(isw,10)
      if((isw1.lt.1.or.isw1.gt.3).or.(isw2.lt.1.or.isw2.gt.4)) then
        write(check(1),'(''isw='',i3)')isw
        call jamerrm(30,1,'(jamdetb2:)invalid isw')
      endif

      idbsw=isw
      srt=srt0
c...Find min. and max. mass for integration.
      emin1=pmas(kc1,1)-pmas(kc1,3)
      emin2=pmas(kc2,1)-pmas(kc2,3)
      em01=pmas(kc1,1)
      em02=pmas(kc2,1)
      id1=kchg(kc1,5)
      if(id1.eq.id_nucls) then
        emin1=emnuc+2*empion+parc(41)
      else if(id1.eq.id_delt) then
        emin1=emnuc+empion+parc(41)
      else if(id1.eq.id_delts) then
        emin1=emnuc+3*empion+parc(41)
      endif
      id2=kchg(kc2,5)
      if(id2.eq.id_nucls) then
        emin2=emnuc+2*empion+parc(41)
      else if(id2.eq.id_delt) then
        emin2=emnuc+empion+parc(41)
      else if(id2.eq.id_delts) then
        emin2=emnuc+3*empion+parc(41)
      endif

      emax1=srt-emin2-0.001d0
      emax2=srt-emin1-0.001d0
      if(emin1.gt.emax1.or.emin2.gt.emax2) then
        if(msel.eq.1) fac=1.d0/(pr*pr)
        if(msel.eq.2) fac=0.0d0
        return
      endif

c...Min. and max. for first integral.
      a=emin1
      b=emax1

c....Initialization
      jmax=0
      nhlf=0
      h=b-a
      s1=jambwf2(kc1,kc2,a,emin2)+jambwf2(kc1,kc2,b,emin2)
      s2=0.0d0
      s4=0.0d0
      sum=0.5d0*h*s1

c...Start the Simpson method
100   continue
      jmax=2*jmax+1
      nhlf=nhlf+1
      h=h*0.5d0
      s2=s2+s4
      s4=0.0d0
      do 110 j=1,jmax,2
       s4=s4+jambwf2(kc1,kc2,a+j*h,emin2)
110   continue
      ss=sum
      sum=(h/3.0d0)*(s1+2.0d0*s2+4.0d0*s4)
      ds=abs(sum-ss)
      if(ds.le.eps*abs(sum)) goto 200
      if(nhlf.eq.mxhlf) then
        write(6,999)'(jamdetb2:)simpson does not converge. mxhlf=',mxhlf
        nhlf=-nhlf
        goto 200
      endif
      goto 100

200   continue

      fac=sum/paru(1)**2
      if(msel.eq.2) return

      fac=1.d0/max(sum/paru(1)**2,1d-5)
c     fac=paru(1)**2/sum
      if(isw2.eq.1) then
        fac=fac/(pr*pr)
      else if(isw2.eq.2.or.isw2.eq.3) then
        fac=fac/pr
        if(isw2.eq.3) fac=fac*(em1/em01)*(em2/em02)
      endif

999   format(a,i3,2(a,1pe22.15))

      end

c***********************************************************************

      real*8 function jambwf2(kc1,kc2,em1,emin2)

c...Provide B-W function for jamdetb2.
      implicit double precision(a-h, o-z)
      real*8 jambwf3
      common/dbdat/srt,em01,em02,idbsw
      data eps/ 1.0d-3/
      data mxhlf/30/

      a=emin2
      b=srt-em1-0.001d0
      if(a.gt.b) then
        jambwf2=0.0d0
        return
      endif

c...Initialization.
      jmax=0
      nhlf=0
      h=b-a
      s1=jambwf3(kc1,kc2,em1,a)+jambwf3(kc1,kc2,em1,b)
      s2=0.0d0
      s4=0.0d0
      sum=0.5d0*h*s1

c...Start the Simpson method.
100   continue
      jmax=2*jmax+1
      nhlf=nhlf+1
      h=h*0.5d0
      s2=s2+s4
      s4=0.0d0
      do 110 j=1,jmax,2
       s4=s4+jambwf3(kc1,kc2,em1,a+j*h)
110   continue
      ss=sum
      sum=(h/3.0d0)*(s1+2.0d0*s2+4.0d0*s4)
      ds=abs(sum-ss)
      if(nhlf.eq.mxhlf) then
        write(6,999)'simpson method does not converge. mxhlf=',mxhlf
      endif
      if(ds.le.eps*abs(sum)) goto 200
      if(nhlf.ge.mxhlf)      goto 200
      goto 100
200   continue
      jambwf2=sum
      return

999   format(a,i3,2(a,1pe22.15))

      end

c***********************************************************************

      real*8 function jambwf3(kc1,kc2,em1,em2)

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      common/dbdat/srt,em01,em02,idbsw
      parameter(maxbr=70)
      dimension ibranch(maxbr),pwid(maxbr)
c...Functions: momentum in two-particle cm.
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a) 

      isw=mod(idbsw,10)
c...Calculate momentum dependence decay width.
      if(idbsw.ge.10) then
        call jamwidm(kc1,1,0,0,0,em1,ibranch,pwid,gam1,itag)
        call jamwidm(kc2,1,0,0,0,em2,ibranch,pwid,gam2,itag)
c...Constant decay width.
      else
        gam1=pmas(kc1,2)
        gam2=pmas(kc2,2)
      endif

c...Non-rel
      if(idbsw.le.20) then
        jambwf3=0.25d0*gam1*gam2/((em1-em01)**2+0.25d0*gam1**2)
     $                     /((em2-em02)**2+0.25d0*gam2**2)
c...Rel1
      else if(idbsw.le.30) then
        jambwf3=4*em1*em01*gam1/((em1**2-em01**2)**2+(em01*gam1)**2)
     $         *em2*em02*gam2/((em2**2-em02**2)**2+(em02*gam2)**2)
c...Rel2
      else if(idbsw.le.40) then
        jambwf3=4*em1**2*gam1/((em1**2-em01**2)**2+(em1*gam1)**2)
     $         *em2**2*gam2/((em2**2-em02**2)**2+(em2*gam2)**2)
      endif

      if(isw.eq.2.or.isw.eq.3) then
        if(srt.ge.em1+em2) then
          pf=sqrt((srt**2-(em1 
     & +em2)**2)*(srt**2-(em1-em2)**2))/(2.d0*srt)
        else
          write(check(1),8000)srt,em1,em2
          write(check(2),8100)chaf(kc1,1),chaf(kc2,1)
          call jamerrm(1,2,'(jambwf3:)srt<em1+em2')
          pf=0.0d0
        endif
        jambwf3=jambwf3*pf
      else if(isw.eq.4) then
        jambwf3=jambwf3*pawt(srt,em1,em2)**2
      endif
 8000 format('srt em1 em2',3(g12.3,1x))
 8100 format(a8,' + ',a8)

      end

c***********************************************************************

      subroutine jamdetb3(kc1,kc2,srt0,em1,em2,pr,sum,nhlf)

c...Calculate the correction factor for RR->NN with constant width.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      parameter(emnuc=0.939d0,empion=0.138d0)
      common/dbdat2/srt,em01,em02,gam1,gam2,emin2
      data eps/ 1.0d-3/
      data mxhlf/30/

      srt=srt0
      id1=kchg(kc1,5)
      id2=kchg(kc2,5)

c...Min. masses.
      emin1=pmas(kc1,1)-pmas(kc1,3)
      if(id1.eq.id_nucls) then
        emin1=emnuc+2*empion+parc(41)
      else if(id1.eq.id_delt) then
        emin1=emnuc+empion+parc(41)
      else if(id1.eq.id_delts) then
        emin1=emnuc+3*empion+parc(41)
      endif

      emin2=pmas(kc2,1)-pmas(kc2,3)
      if(id2.eq.id_nucls) then
        emin2=emnuc+2*empion+parc(41)
      else if(id2.eq.id_delt) then
        emin2=emnuc+empion+parc(41)
      else if(id2.eq.id_delts) then
        emin2=emnuc+3*empion+parc(41)
      endif

c...Pole masses.
      em01=pmas(kc1,1)
      em02=pmas(kc2,1)

c...Decay width.
      gam1=pmas(kc1,2)/2.d0
      gam2=pmas(kc2,2)/2.d0

c...Max. masses.
      emax1=srt-emin2-0.001d0
      emax2=srt-emin1-0.001d0

      if(emin1.gt.emax1.or.emin2.gt.emax2) then
        fac=1.d0
        return
      endif

c...Min. and max. for first integral.
      a=emin1
      b=emax1

c...Initialization.
      jmax=0
      nhlf=0
      h=b-a
      s1=gam1/((a-em01)**2+gam1**2)
     $    *( atan((srt-a-em02)/gam2)-atan((emin2-em02)/gam2) )
     $  +gam1/((b-em01)**2+gam1**2)
     $    *( atan((srt-b-em02)/gam2)-atan((emin2-em02)/gam2) )
      s2=0.0d0
      s4=0.0d0
      sum=0.5d0*h*s1

c...Start of the Simpson method.
100   continue
      jmax=2*jmax+1
      nhlf=nhlf+1
      h=h*0.5d0
      s2=s2+s4
      s4=0.0d0
      do 110 j=1,jmax,2
       x=a+j*h
       s4=s4+gam1/((x-em01)**2+gam1**2)
     $    *( atan((srt-x-em02)/gam2)-atan((emin2-em02)/gam2) )
110   continue
      ss=sum
      sum=(h/3.0d0)*(s1+2.0d0*s2+4.0d0*s4)
      ds=abs(sum-ss)
      if(ds.le.eps*abs(sum)) goto 200
      if(nhlf.eq.mxhlf) then
        write(6,999)'simpson method does not converge. mxhlf=',mxhlf
        goto 200
      endif
      goto 100

200   continue
      sum=max(sum/paru(1)**2,1.0d-5)
      return

999   format(a,i3,2(a,1pe22.15))

      end


c***********************************************************************

      function jamcpair(id1,id2)

c...Order collision pair.
      implicit double precision(a-h, o-z)

      idmin=min(id1,id2)
      idmax=max(id1,id2)
      jamcpair=(idmax*(idmax-1))/2+idmin

      end

c***********************************************************************
c...Cross section part
c***********************************************************************

      double precision function jamrgg92(srt,i)

c...Total cross section fit by Regge theory based formula.
c...Quark counting rule is used for the unknown reactions.
c...Fit were made using experimental data with srt>6GeV.
c...A. Donnachie and P.V.Landshoff, Phys. Lett. B296,277(1992)

      implicit double precision(a-h, o-z)
      dimension xpar(36),ypar(36)

C...=  1 : p + p;
C...=  2 : pbar + p;
C...=  3 : p + n;
C...=  4 : p~ + n;
C...=  5 : lambda + p;
C...=  6 : lambda~ + p;
C...=  7 : xi + p;
C...=  8 : xi~ + p;
C...=  9 : lambda_c + p;
C...= 10 : lambda_c~ + p;

C...= 11 : pi+ + p;
C...= 12 : pi- + p;
C...= 13 : pi0 + p;
C...= 14 : k+ + p;
C...= 15 : k- + p;
C...= 16 : phi + p;
C...= 17 : J/psi + p;

C...= 21 : rho + rho;
C...= 22 : rho + phi;
C...= 23 : rho + J/psi;
C...= 24 : phi + phi;
C...= 25 : phi + J/psi;
C...= 26 : J/psi + J/psi;

C...= 31 : gamma + p (DL);
C...= 32 : gamma + p (VDM).
C...= 33 : gamma + pi (DL);
C...= 34 : gamma + pi (VDM);
C...= 35 : gamma + gamma (DL);
C...= 36 : gamma + gamma (VDM).

C...X and Y parameters of sigmatot = X * s**epsilon + Y * s**(-eta).
      data xpar/4*21.70d0, 2*19.89d0,  4*0.0d0,
     $   3*13.63d0, 2*11.82d0, 10.01d0,0.970d0,3*0.d0,
     &   8.56d0,   6.29d0, 0.609d0,
     $   4.62d0, 0.447d0, 0.0434d0,4*0.d0,
     &   0.0677d0,0.0534d0,0.0425d0,0.0335d0,2.11d-4,1.31d-4/

      data ypar/56.08d0, 98.39d0, 54.77d0,92.72d0, 45.11d0,5*0.0d0,
     $  27.56d0, 36.02d0, 31.79d0, 8.15d0,26.36d0,     -1.51d0,-0.146d0,
     &  3*0.d0,
     &  13.08d0, -0.62d0, -0.060d0, 0.030d0, -0.0028d0, 0.00028d0, 
     & 4*0.d0,
     &   0.129d0,0.115d0,0.081d0,0.072d0,2.15d-4,1.70d-4/

      jamrgg92=xpar(i)*(srt*srt)**0.0808d0
     $        +ypar(i)*(srt*srt)**(-0.4525d0)

      end

c**********************************************************************

      double precision function jamrgg96(srt,i)

c...fit for high energy cross sections by regge theory
c...(particle data group  1996 Phys. Rev. D54)
c...sig_tot=X*s^eta + Y*s^eps
c   1: pp    total    (50GeV/c)
c   2: p~ p  total    (50GeV/c)
c   3: np    total    (50GeV/c)
c   4: p~ n  total    (50GeV/c)
c   5: pd    total    (50GeV/c)
c   6: p~ d total     (50GeV/c)
c   7: pi+ p total    (10GeV/c)
c   8: pi- p total    (10GeV/c)
c   9: pi+- d total   (10GeV/c)
c  10: k+ p  total    (10GeV/c)
c  11: k+ n  total    (10GeV/c)
c  12: k- p  total    (10GeV/c)
c  13: k- n  total    (10GeV/c)
c  14: k+ d  totol    (10GeV/c)
c  15: k- d  total    (10GeV/c)
c  16: gamm p total   (12GeV/c)

      implicit double precision(a-h, o-z)
      dimension x(16),y(16),eta(16),eps(16)

      data eta/4*0.46d0,2*0.45d0,2*0.45d0,0.43d0,4*0.5d0,2*0.47d0,
     $         0.46d0/
      data eps/4*0.079d0,2*0.090d0,2*0.079d0,0.088d0,4*0.079d0, 
     & 2*0.082d0,
     $         0.075d0/
      data x/2*22.d0,2*22.3d0, 2*35.7d0,
     $ 2*13.7d0,23.2d0,  4*12.2d0, 2*21.7d0, 0.071d0/
      data y/56.1d0,98.2d0,55.0d0,92.7d0,  179.0d0,270.6d0,
     $  27.8d0,35.9d0, 85.5d0,
     $  2*8.3d0,2*26.4d0,26.2d0,64.8d0,
     $  0.12d0/

      jamrgg96=x(i)*(srt*srt)**eps(i) + y(i)*(srt*srt)**(-eta(i))

      end

c**********************************************************************

      double precision function jamchc96(i,plab)

c...Fit formula for cross sections (particle data group  1996)
c...by the CERN-HERA and COMPAS groups.
c   1: gamm p total   (3.0-183)
c   2: gamm d total   (2.0-17.8)
c   3: pi+ p total    (4.0-340)
c   4: pi+ p elastic  (2.0-200)
c   5: pi- p total    (2.5-370)
c   6: pi- p elastic  (2.0-360)
c   7: pi+- d total   (2.5-370)

c   8: k+ p  total    (2.0-310)
c   9: k+ p  elastic  (2.0-175)
c  10: k+ n  total    (2.0-310)
c  11: k+ d  totol    (2.0-310)
c  12: k- p  total    (3.0-310)
c  13: k- p  elastic  (3.0-175)
c  14: k- n  total    (1.8-310)
c  15: k- d  total    (3.0-310)

c  16: pp    total    (3.0-2100)
c  17: pp    elastic  (2.0-2100)
c  18: np    total    (2.0-370)
c  19: pd    total    (3.0-370)
c  20: pd    elastic  (3.0-384)

c  21: p~ p total     (5.0-1730000)
c  22: p~ p elastic   (5.0-1.73e+6)
c  23: p~ n total     (1.1-280)
c  24: p~ n elastic   (1.1-5.55)
c  25: p~ d total     (2.0-280)
c
c  26: k+ n ==> k0 p charge exchange    (2.0-12.8)
c  27: k- p ==> k0 n charge exchange    (3.0-40)
c  28: lambda p total (0.12-21.0)
c  29: lambda p elastic (2.0-24.0)
c
c...sig(x) = A + B*x**n + C*(log(x))**2 + D*log(x)

      implicit double precision(a-h, o-z)
      dimension a(29),b(29),xn(29),c(29),d(29)

      data a/0.147d0, 0.3d0,
     $  16.4d0, 0.0d0, 33.0d0, 1.76d0, 56.8d0,
     $  18.1d0, 5.0d0, 18.7d0, 34.2d0, 32.1d0, 7.3d0, 25.2d0, 57.6d0,
     $  48.0d0, 11.9d0, 47.30d0, 91.3d0, 16.1d0,
     $  38.4d0, 10.2d0, 0.0d0, 36.5d0, 112,
     $  -1.439d0, -0.2173d0,  18.0d0, 3.49d0/

      data b/2*0.0d0,
     $   19.3d0, 11.4d0, 14.0d0, 11.2d0, 42.2d0,
     $   0.0d0, 8.1d0, 0.0d0, 7.9d0, 4*0.0d0,
     $   0.0d0, 26.9d0, 3*0.0d0,
     $   77.6d0, 52.7d0, 133.6d0, 0.0d0, 125.0d0,
     $   8.279d0, 3.641d0,0.121d0,26.2d0/

      data xn/2*0.0d0,
     $    -0.42d0, -0.4d0, -1.36d0, -0.64d0, -1.45d0,
     $     0.0d0, -1.8d0, 0.0d0, -2.1d0, 4*0.0d0,
     $     0.0d0,-1.21d0,3*0.0d0,
     $     -0.64d0, -1.16d0, -0.70d0, 0.0d0, -1.08d0,
     $     -1.618d0, -1.661d0,-3.92d0,-1.01d0/

      data c/0.0022d0,0.0095d0,
     $ 0.19d0, 0.079d0, 0.456d0, 0.043d0, 0.65d0,
     $ 0.26d0, 0.16d0,  0.21d0,  0.346d0, 0.66d0, 0.29d0, 0.38d0,1.17d0,
     $ 0.522d0, 0.169d0, 0.513d0, 1.05d0, 0.32d0,
     $ 0.26d0, 0.125d0, -1.22d0, 0.0d0, 1.14d0,
     $ -0.156d0, -0.0171d0, 6.38d0, 0.0d0/

      data d/-0.017d0,-0.057d0,
     $    2*0.0d0, -4.03d0,0.0d0,-5.39d0,
     $   -1.0d0,-1.3d0,-0.89d0,-0.99d0,-5.6d0,-2.40d0,-2.9d0,-9.5d0,
     $   -4.51d0,-1.85d0,-4.27d0,-8.8d0,-3.4d0,
     $   -1.2d0,-1.28d0,13.7d0,-11.9d0,-12.4d0,
     $    0.928d0, 0.1217d0,2*0.0d0/

c     data pmin/3.0,2.0,  4.0,2.0,2.5,2.0,2.5,
c    $ 2.0,2.0,2.0,2.0,3.0,3.0,1.8,3.0,
c    $ 3.0,2.0,3.0,3.0,2.0,
c    $ 5.0,5.0,1.1,1.1,2.0/

      if(i.ge.1.and.i.le.29) then
        jamchc96=a(i)+b(i)*plab**xn(i)+c(i)*(log(plab))**2
     $                                         +d(i)*log(plab)
      else
        print *,'(jamchc96:) invalid i=',i
        jamchc96=0.0d0
      endif

      end

c**********************************************************************

      double precision function jamchc88(plab,i)

c...Fit formula for cross sections (particle data group  1988)
c...high-energy parametrizations by CERN-HERA and COMPAS group
c
c     i=1: lambda p total (0.12-21.0)  2: lambda p elastic (2.0-24.0)
c       3: pi+ p total    (4.0-34.0)   4: pi+ p elastic    (2.0-250)
c       5: pi- p total    (2.5-370)    6: pi- p elastic    (2.0-360)
c       7: k+ p  total    (2.0-310)    8: k+ p  elastic    (1.5-250)
c       9: k+ n  total    (2.0-310)
c      10: k- p  total    (3.0-310)    11: k- p  elastic    (2.0-175)
c      12: k- n  total    (2.5-310)
c      13: pp    total    (3.0-2100)   14: pp    elastic    (2.0-2100)
c      15: np    total    (2.0-280)
c      16: pba p total    (5.0-432000) 17: pba p elastic   (2.0-159000)
c      18: pba n total    (1.13-280)   19: pba n elastic   (1.13-5.55)
c
c      20: k+ n ==> k0 p charge exchange    (2.0-12.8)
c      21: k- p ==> k0 n charge exchange    (3.0-40)
c
c k+ n charge exchange
c a = -1.43876
c b = 8.27893
c n = -1.61792
c c = -0.156198
c d = 0.9276
csig(x) = a + b*x**n + c*(log(x))**2 + d*log(x)
c k- p charge exchange
c a = -0.217313
c b = 3.64072
c n = -1.66109
c c = -0.0171123
c d = 0.121654

      implicit double precision(a-h, o-z)
      dimension a(21),b(21),xn(21),c(21),d(21)
      data a/18.0d0, 3.49d0, 32.10d0, 7.07d0, 33.1d0, 1.73d0, 17.1d0, 
     &  5.84d0, 18.4d0,
     a       -21.1d0, 7.24d0, -1040, 45.6d0, 11.2d0, 47.4d0, 41.1d0, 
     &  10.6d0,
     a        41.9d0, 37.5d0, -1.439d0, -0.2173d0/
      data b/0.121d0, 26.2d0, 48700, 11.3d0, 15.0d0, 11.2d0, 5.54d0, 
     &  17.2d0, 175,
     a       56.2d0,  46.0d0, 1060,  219,  25.5d0, -100, 77.2d0, 53.1d0,
     &  96.2d0,
     a       -2.63d0, 8.279d0, 3.641d0/
      data xn/-3.92d0, -1.01d0, -7.85d0, -1.6d0,  -1.41d0, -0.63d0, 
     &  -2.67d0, -3.06d0,
     a        -7.85d0, -0.27d0, -4.71d0, -0.03d0, -4.23d0, -1.12d0, 
     &  -4.57d0, -0.68d0,
     a        -1.19d0, -0.99d0, -2.58d0, -1.618d0, -1.661d0/
      data c/6.38d0,   0.0d0,  0.540d0, 0.16d0, 0.458d0, 0.0040d0, 
     &  0.139d0, 0.206d0,
     a       0.198d0, -0.155d0,0.279d0, 0.0d0,  0.41d0,  0.151d0, 
     &   0.512d0, 0.293d0,
     a       0.136d0, -0.154d0,-12.6d0, -0.156d0, -0.0171d0/
      data d/0.0d0, 0.0d0, -4.41d0, -1.56d0,-4.06d0,0.0d0,-0.27d0, 
     & -1.71d0,-0.753d0,6.24d0,
     a       -2.35d0, 27.8d0, -3.41d0, -1.62d0, -4.29d0, -1.82d0, 
     & -1.41d0, 0.0d0,0.0d0,
     a        0.928d0, 0.1217d0/

      if(i.le.0.or.i.gt.21.or.plab.le.1d-5) then
        print *,'jamchc88 i=',i,plab
        stop
      endif
      jamchc88=a(i)+b(i)*plab**xn(i)+c(i)*(log(plab))**2+d(i)*log(plab)

      end

c***********************************************************************

      subroutine jamsighh(sig,isig,s)

c...Give tableted low energy hh cross sections.

c...isig= 1: pp total
c         2: pp elastic
c         3: pn total
c         4: pn elastic
c         5: app total
c         6: app elastic
c         7: pion+ p total
c         8: pion+ p elsstic
c         9: pion- p total
c        10: pion- p elsstic
c        11: k- p total
c        12: k- p elastic
c        13: k- n total
c        14: k- n elastic
c        15: k+ p total
c        16: k+ p elastic
c        17: k+ n total
c        18: k+ n elastic
c        19: k- p -> k0 n charge exchange cross section
c        20: k+ n -> k0 p charge exchange cross section
c        21: k- p => lambda pi0                                        
c        22: k- p => sigma- pi+                                       
c        23: k- p => sigma+ pi-                                        
c        24: lambda p total
c        25: lambda p elastic
c
c        26: pp => pp pi0  sigma_{11}
c        27: pp => pp pi0  sigma_{10}
c        28: pp => pp pi0  sigma_{01}
c        29: factor for calculatin of dn->nn by delailed balance.
c
c....Table parameters
      implicit double precision(a-h, o-z)
      parameter (itblsz=100)
      parameter (tblhig=5.075668d0)
      parameter (tbllow=1.886d0,tlow=0.6344582d0) ! tlow=log(tbllow)
      parameter (tblstp=0.01d0)
      parameter (isgmax=29)
      common /sigma1/ sigfit(isgmax,itblsz)
 
      if(isig.ge.7.and.isig.le.10) then
        tblsft=0.9d0
      else if(isig.ge.11.and.isig.le.23) then
        tblsft=0.444d0
      else
        tblsft=0.0d0
      endif
      xpt=log(s+tblsft)
      if(xpt.lt.tlow) then
        xpt=tlow
      end if
      index=int((xpt-tlow)/tblstp)+1
      if(index.ge.itblsz) then
        write(6,*)'warning:index.ge.itblsz in sigma index=',index
        write(6,*)'srt isig',s,isig
        index=itblsz-1
      end if

c...Find slopes and cross sections
      x1=(index-1)*tblstp+tlow
      y1in=sigfit(isig,index)
      y2in=sigfit(isig,index+1)
      slin=(y2in-y1in)/tblstp
      sig=slin*(xpt-x1)+y1in

      end

c***********************************************************************

      subroutine jamxnnin(x,sigin,iporn)

c...Output NN inelastic resonance production cross section for given x.
c...(INPUT)
c     x: c.m.enrgy in GeV
c     iporn:
c        = 1: t=1 cross sections for pp
c        = 2: t=1 cross sections for pn:
c               only deuteron contribution is different
c        = 0: t=0 crosss sections
c   (OUTPUT)
c     sigin: inelastic cross sections
c
c...1) NN -> ND  (deuteron production is included)
c...2) NN -> NN*
c...3) NN -> DD
c...4) NN -> ND*
c...5) NN -> N*D
c...6) NN -> DD*
c...7) NN -> N*N*
c...8) NN -> N*D*
c...9) NN -> D*D*
c...10)NN-> s-wave pion

      implicit double precision(a-h, o-z)
      dimension sigin(10),a1(5,9),a0(5,9),ad1(5),ad0(5),as1(5),as0(5)

c...NN->ND(1232)
      data rppd0,gppd0,sppd0,thd0/
     $ 0.633719d0,  2.11477d0,0.0171405d0,2.0139999d0/
      data sgmppd0, sppd1, fppd1, gppd1, thd1/
     $ 48.64510d0,  2.06672d0, 0.480085d0, 0.576422d0, 2.124d0/   

c...T=1
      data a1/
     $  0.00000d0,  0.00000d0,  0.00000d0, 0.000000d0, 0.000d0,
     $ 24.94700d0,  2.48150d0,  2.63330d0, 0.425358d0, 2.162d0,
     $  7.63181d0,  1.41140d0,  2.67784d0, 0.311722d0, 2.252d0,
     $  8.01615d0,  2.74161d0,  3.34503d0, 0.259703d0, 2.340d0,
     $ 13.14580d0,  2.06775d0,  2.75682d0, 0.247810d0, 2.300d0,
     $ 19.63220d0,  2.01946d0,  2.80619d0, 0.297073d0, 2.528d0,
     $ 11.67320d0,  2.31682d0,  2.96359d0, 0.259223d0, 2.438d0,
     $  2.99086d0,  2.29380d0,  3.54392d0, 0.090438d0, 2.666d0,
     $ 35.13780d0,  2.25498d0,  3.14299d0, 0.215611d0, 2.804d0/
c...T=0
      data a0/
     $  0.00000d0,  0.00000d0,  0.00000d0, 0.000000d0, 0.000d0,
     $166.60600d0,  2.10128d0,  2.34635d0, 0.284955d0, 2.162d0,
     $ 39.99770d0,  1.83576d0,  2.40348d0, 0.288931d0, 2.252d0,
     $  0.00000d0,  0.00000d0,  0.00000d0, 0.000000d0, 0.000d0,
     $  0.00000d0,  0.00000d0,  0.00000d0, 0.000000d0, 0.000d0,
     $ 56.32490d0,  2.00679d0,  2.71312d0, 0.362132d0, 2.528d0,
     $  2.14575d0,  0.21662d0,  3.40108d0, 0.252889d0, 2.438d0,
     $  0.00000d0,  0.00000d0,  0.00000d0, 0.000000d0, 0.000d0,
     $  4.14197d0,  1.67026d0,  3.75133d0, 0.476595d0, 2.804d0/


c...deuteron 1
      data ad1/ .14648d0,   .20807d0,  2.13072d0,  .042475d0, 2.024d0/
c...deuteron 2
      data ad0/ .12892d0,   .08448d0,  2.18138d0,  .059207d0, 2.054d0/
c...s-wave pion production pp
      data as1/15.644100d0, 1.675220d0,  2.07706d0,  .658047d0, 2.014d0/
c...s-wave pion production pn
      data as0/78.868103d0,  .746742d0,  1.25223d0,  .404072d0, 2.014d0/

c...Definition of Fitting Function.
      bw(a,b)=b/((a**2-1)**2+b**2)

      do i=1,10
        sigin(i)=0.0d0
      end do

c...I=1 cross sections
      if(iporn.eq.1.or.iporn.eq.2) then

c...Low energy.
      if(x.le.thd0) return

c...NN->ND(1232)
      if(x.le.thd0) then
         sigd0 = 0.0d0
      else
         sigd0 = sgmppd0*rppd0*sqrt((x-thd0)/thd0)*gppd0
     $         /((x-sppd0)**2+gppd0**2)/100
      endif

      if(x.le.thd1) then
         sigd1=0.0d0
      else
       sigd1=sgmppd0*(x/thd1-1)**fppd1*bw(x/sppd1,gppd1)
      endif

      sigd=sigd1+sigd0  ! pp->ND  (-> NNpi)

c...Deuteron pi.
      if(x.le.ad1(5)) then
        sigdeut1=0.0d0
      else
        sigdeut1 = ad1(1)*(x/ad1(5)-1)**ad1(2)*bw(x/ad1(3),ad1(4))
      endif
      if(x.le.ad0(5)) then
         sigdeut2 = 0.0d0
      else
         sigdeut2 = ad0(1)*(x/ad0(5)-1)**ad0(2)*bw(x/ad0(3),ad0(4))
      endif

c...[1]NN -> ND,dpi
      if(iporn.eq.1) sigin(1)=sigd+sigdeut1
      if(iporn.eq.2) sigin(1)=sigd+sigdeut1+sigdeut2

      do i=2,9
       if(x.gt.a1(5,i))
     $  sigin(i)=a1(1,i)*(x/a1(5,i)-1)**a1(2,i)*bw(x/a1(3,i),a1(4,i))
      end do

c...[10]s-wave pion
      if(x.gt.as1(5))
     $  sigin(10)=as1(1)*(x/as1(5)-1)**as1(2)*bw(x/as1(3),as1(4))


c...T=0 cross sections.
      else

      do i=2,9
       if(a0(5,i).gt.2.0d0.and.x.gt.a0(5,i))
     $  sigin(i)=a0(1,i)*(x/a0(5,i)-1)**a0(2,i)*bw(x/a0(3,i),a0(4,i))
      end do

c...[10] s-wave pion
      if(x.gt.as0(5))
     $  sigin(10)=as0(1)*(x/as0(5)-1)**as0(2)*bw(x/as0(3),as0(4))


      endif

      end

c***********************************************************************

      double precision function jamsigkn(isig,s,plab)

c...Purpose: to give kaon-nucleon one/two-pion production x-sections.
c...1:K+ N => K delta  isospin=1
c...2:K+ p => K+*(890) p
c...3:K+ n => K+*(890) n
c...4:K+ n => K0*(890) p
c...5:K+ N => D(1232)+K(892) T=1

      implicit double precision(a-h, o-z)
      dimension a1(5),a2(5),a3(5),a4(5),b1(5),b2(5),sth(5),plth(5)
      data sth/4*1.5719d0,1.71231d0/
c     data a1/0.584417d0,0.519173d0,0.309147d0,0.465779d0,.262122d0/
      data a1/0.584417d0,0.519173d0,0.309147d0,0.465779d0,.349496d0/
      data a2/1.09807d0,1.50663d0,1.37519d0,0.720443d0,0.929574d0/
      data a3/1.83136d0,1.95501d0,1.97189d0,1.99835d0,2.50283d0/
      data a4/0.0365562d0,0.0531106d0,0.0972668d0,0.0681271d0
     $       ,0.143572d0/
      data plth/4*2.5d0,4.0d0/
c     data b1/5.44791d0,6.85573d0,2.40975d0,14.0587d0,11.1962d0/
      data b1/5.44791d0,6.85573d0,2.40975d0,14.0587d0,14.92827d0/
      data b2/1.66538d0,1.70302d0,1.18874d0,2.28464d0,1.79165d0/

      if(s.lt.sth(isig)) then
        jamsigkn=0.0d0
      else if(plab.le.plth(isig)) then
        jamsigkn=a1(isig)*(s-sth(isig))**a2(isig)
     $                 /((s-a3(isig))**2+a4(isig))
      else
        jamsigkn=b1(isig)*plab**(-b2(isig))
      endif

      end

c***********************************************************************

      block data jamsigda

c...Data for low-energy hh cross sections.

      implicit double precision(a-h, o-z)
      parameter(itblsz=100)
      parameter(isgmax=29)

      common /sigma1/sigfit(isgmax,itblsz)
c     common/sigma2/signd(3,itblsz),cofdel(itblsz)

c....pp total

      data (sigfit(1,i),i=1,itblsz)/
c    & 179.967, 42.241,30.861,26.501,24.533,23.915,23.756,23.945,
     & 120.967, 40.241,25.861,24.501,22.533,21.915,21.756,21.945,
     & 24.452, 25.224,26.807,29.646,32.835,36.025,39.013,41.575,
     & 43.593, 45.088,46.103,46.741,47.095,47.256,47.296,47.261,
     & 47.180, 47.066,46.923,46.762,46.585,46.395,46.198,45.994,
     & 45.788, 45.581,45.376,45.172,44.972,44.775,44.581,44.390,
     & 44.200, 44.013,43.828,43.648,43.474,43.307,43.148,42.998,
     & 42.859, 42.730,42.611,42.500,42.398,42.308,42.223,42.145,
     & 42.075, 42.010,41.948,41.892,41.838,41.785,41.735,41.685,
     & 41.635, 41.584,41.533,41.480,41.425,41.369,41.311,41.252,
     & 41.190, 41.128,41.064,40.999,40.934,40.869,40.804,40.740,
     & 40.677, 40.614,40.554,40.495,40.437,40.381,40.326,40.273,
     & 40.221, 40.170,40.120,40.071,40.022,39.974,39.926,39.879,
     & 39.832, 39.786, 39.740, 39.695
     &/

c...pp elastic

      data (sigfit(2,i),i=1,itblsz)/
     & 120.967, 40.241,25.861,24.501,22.533,21.915,21.756,21.165,
     & 23.494, 23.823,23.773,23.865,24.257,24.560,24.574,24.403,
     & 24.297, 24.271,24.248,24.150,23.947,23.640,23.276,22.888,
     & 22.496, 22.113,21.742,21.386,21.045,20.717,20.398,20.087,
     & 19.779, 19.473,19.164,18.850,18.533,18.212,17.893,17.577,
     & 17.268, 16.968,16.678,16.398,16.129,15.868,15.619,15.380,
     & 15.148, 14.924,14.710,14.501,14.299,14.105,13.915,13.731,
     & 13.552, 13.377,13.208,13.042,12.880,12.723,12.568,12.417,
     & 12.270, 12.125,11.984,11.847,11.713,11.583,11.458,11.336,
     & 11.218, 11.106,11.000,10.899,10.804,10.716,10.635,10.561,
     & 10.494, 10.435,10.382,10.336,10.294,10.258,10.225,10.196,
     & 10.169, 10.143,10.118,10.094,10.069,10.044,10.018, 9.991,
     &  9.962,  9.931,  9.899,  9.866
     &/

c...pn total

      data (sigfit(3,i),i=1,itblsz)/
c    & 485.849, 151.193,83.194,62.123,52.687,46.879,42.924,40.170,
     & 385.849, 118.193,73.194,62.123,48.687,46.879,40.924,38.170,
     & 37.249, 36.990,36.222,35.832,35.718,35.791,35.980,36.232,
     & 36.515, 36.811,37.115,37.425,37.743,38.068,38.402,38.740,
     & 39.078, 39.413,39.739,40.053,40.349,40.623,40.874,41.100,
     & 41.299, 41.469,41.610,41.720,41.801,41.856,41.887,41.898,
     & 41.893, 41.876,41.849,41.816,41.778,41.736,41.691,41.644,
     & 41.593, 41.540,41.484,41.426,41.365,41.301,41.235,41.166,
     & 41.097, 41.025,40.952,40.879,40.806,40.734,40.662,40.592,
     & 40.524, 40.458,40.395,40.334,40.277,40.223,40.172,40.124,
     & 40.081, 40.038,40.000,39.963,39.930,39.898,39.867,39.838,
     & 39.809, 39.781,39.753,39.725,39.697,39.668,39.639,39.609,
     & 39.576, 39.543,39.508,39.472,39.436,39.396,39.357,39.315,
     & 39.273, 39.230, 39.185, 39.141
     &/

c...pn elastic
      data (sigfit(4,i),i=1,itblsz)/
     & 385.849,118.193,73.194, 62.123, 48.687, 46.879, 40.924, 37.949,
     & 36.585, 35.783, 34.351, 32.354, 30.777, 28.331, 26.629, 25.329,
     & 24.297, 24.271,24.248,24.150,23.947,23.640,23.276,22.888,
     & 22.496, 22.113,21.742,21.386,21.045,20.717,20.398,20.087,
     & 19.779, 19.473,19.164,18.850,18.533,18.212,17.893,17.577,
     & 17.268, 16.968,16.678,16.398,16.129,15.868,15.619,15.380,
     & 15.148, 14.924,14.710,14.501,14.299,14.105,13.915,13.731,
     & 13.552, 13.377,13.208,13.042,12.880,12.723,12.568,12.417,
     & 12.270, 12.125,11.984,11.847,11.713,11.583,11.458,11.336,
     & 11.218, 11.106,11.000,10.899,10.804,10.716,10.635,10.561,
     & 10.494, 10.435,10.382,10.336,10.294,10.258,10.225,10.196,
     & 10.090, 10.132, 10.167, 10.199, 10.228, 10.251, 10.274, 10.295,
     & 10.316, 10.338, 10.359, 10.382
     &/

c...pp~ total
      data (sigfit(5,i),i=1,itblsz)/
     & 235.72, 227.5,185.12,165.21,151.52,143.79,138.05,133.22,
     & 128.94, 125.04,121.44,118.11,115.01,112.13,109.47,107.01,
     & 104.768,102.72,100.85,99.14,97.568,96.100,94.709,93.366,
     & 92.043, 90.725,89.410,88.105,86.823,85.580,84.391,83.270,
     & 82.229, 81.266,80.374,79.534,78.728,77.940,77.161,76.388,
     & 75.624, 74.870,74.130,73.405,72.696,72.003,71.329,70.673,
     & 70.036, 69.418,68.821,68.244,67.688,67.151,66.635,66.137,
     & 65.658, 65.199,64.757,64.331,63.921,63.527,63.148,62.781,
     & 62.427, 62.084,61.751,61.427,61.110,60.801,60.497,60.198,
     & 59.904, 59.613,59.324,59.037,58.751,58.466,58.182,57.898,
     & 57.614, 57.330,57.046,56.763,56.479,56.197,55.916,55.636,
     & 55.359, 55.083,54.811,54.543,54.278,54.018,53.763,53.513,
     & 53.269, 53.031, 52.799, 52.573
     &/

c...pp~ elastic
      data (sigfit(6,i),i=1,itblsz)/
     & 73.829, 75.373,66.216,59.243,56.159,54.083,52.219,50.518,
     & 48.973, 47.571,46.303,45.153,44.108,43.150,42.260,41.419,
     & 40.611, 39.825,39.053,38.287,37.521,36.746,35.959,35.157,
     & 34.344, 33.523,32.701,31.885,31.082,30.300,29.545,28.822,
     & 28.132, 27.478,26.859,26.274,25.720,25.195,24.698,24.225,
     & 23.774, 23.343,22.930,22.532,22.150,21.781,21.425,21.080,
     & 20.744, 20.420,20.103,19.795,19.495,19.203,18.918,18.640,
     & 18.368, 18.102,17.843,17.589,17.342,17.101,16.865,16.635,
     & 16.411, 16.192,15.979,15.772,15.570,15.373,15.182,14.997,
     & 14.817, 14.642,14.472,14.307,14.148,13.994,13.845,13.700,
     & 13.559, 13.422,13.290,13.162,13.037,12.916,12.798,12.684,
     & 12.573, 12.464,12.358,12.255,12.154,12.055,11.959,11.865,
     & 11.773, 11.682, 11.594, 11.506
     &/

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c...pion+ p total
      data (sigfit(7,i),i=1,itblsz)/
     &  0.000,  0.000, 0.000, 0.000, 0.000, 6.058, 6.147,15.689,
     & 33.768, 68.878,120.614,173.692,192.154,160.644,123.116,93.466,
     & 70.706, 53.969,42.287,34.251,28.386,23.954,20.633,18.191,
     & 16.612, 15.924,16.127,17.231,19.039,21.065,22.910,24.363,
     & 25.480, 26.576,27.954,29.801,32.141,34.786,37.266,38.929,
     & 39.316, 38.584,37.202,35.606,34.088,32.791,31.763,31.003,
     & 30.494, 30.173,29.993,29.921,29.902,29.910,29.919,29.914,
     & 29.881, 29.814,29.711,29.577,29.417,29.244,29.067,28.893,
     & 28.728, 28.570,28.417,28.265,28.109,27.946,27.778,27.606,
     & 27.435, 27.270,27.111,26.961,26.819,26.683,26.554,26.432,
     & 26.316, 26.209,26.110,26.018,25.931,25.856,25.783,25.717,
     & 25.654, 25.597,25.539,25.487,25.435,25.384,25.336,25.289,
     & 25.241, 25.194, 25.148, 25.103
     &/

c...pion+ p elastic
      data (sigfit(8,i),i=1,itblsz)/
     &  0.000,  0.000,  0.000,  0.000,  0.000,  6.058,  6.147, 15.689, 
     & 33.768, 68.878, 120.614,173.692,192.154,160.644,123.116,93.466, 
     & 70.706, 53.969, 42.287, 34.251, 28.386, 22.630, 18.205, 15.126, 
     & 12.819, 11.155, 10.137,  9.649,  9.595,  9.894, 10.464, 11.204, 
     & 12.014, 12.837, 13.663, 14.496, 15.324, 16.099, 16.722, 17.037, 
     & 16.896, 16.300, 15.412, 14.421, 13.459, 12.598, 11.870, 11.271, 
     & 10.783, 10.382, 10.047,  9.759,  9.505,  9.279,  9.070,  8.878, 
     &  8.697,  8.525,  8.362,  8.205,  8.054,  7.909,  7.768,  7.633, 
     &  7.502,  7.375,  7.253,  7.136,  7.022,  6.912,  6.806,  6.704, 
     &  6.606,  6.511,  6.420,  6.333,  6.248,  6.166,  6.087,  6.011, 
     &  5.937,  5.866,  5.797,  5.731,  5.666,  5.603,  5.543,  5.483, 
     &  5.426,  5.370,  5.316,  5.263,  5.212,  5.162,  5.113,  5.066, 
     &  5.020,  4.975,  4.931,  4.889
     &/

c...pion- p total
      data (sigfit(9,i),i=1,itblsz)/
     &  8.102,  8.230, 8.357, 8.485, 8.613, 8.740, 8.868,11.766,
     & 16.957, 26.810,44.361,61.961,67.514,57.631,44.819,35.576,
     & 30.090, 27.407,26.941,27.539,28.711,30.694,33.987,38.823,
     & 43.836, 44.236,41.079,38.625,39.577,44.125,51.113,56.464,
     & 54.333, 48.098,42.487,38.917,37.131,36.450,36.219,36.013,
     & 35.652, 35.212,34.871,34.715,34.726,34.843,34.998,35.136,
     & 35.218, 35.223,35.148,34.997,34.786,34.533,34.251,33.954,
     & 33.654, 33.359,33.073,32.801,32.543,32.296,32.062,31.837,
     & 31.620, 31.410,31.204,31.004,30.808,30.617,30.433,30.254,
     & 30.083, 29.919,29.762,29.610,29.464,29.324,29.189,29.061,
     & 28.938, 28.820,28.708,28.600,28.497,28.396,28.299,28.203,
     & 28.109, 28.017,27.925,27.835,27.745,27.656,27.567,27.480,
     & 27.393, 27.306, 27.221, 27.137
     &/

c...pion- p elastic
      data (sigfit(10,i),i=1,itblsz)/
     &  1.687,  1.714, 1.740, 1.767, 1.793, 1.820, 1.847, 3.126,
     &  5.355,  8.770,13.946,19.657,21.690,20.101,17.433,15.358,
     & 13.889, 12.908,12.404,12.384,12.863,13.809,15.114,16.548,
     & 17.656, 17.832,17.163,16.617,16.944,18.415,20.763,22.696,
     & 22.511, 20.511,18.126,16.099,14.536,13.344,12.426,11.722,
     & 11.188, 10.778,10.451,10.177, 9.938, 9.723, 9.525, 9.339,
     &  9.163,  8.995, 8.834, 8.679, 8.528, 8.382, 8.240, 8.102,
     &  7.968,  7.838, 7.712, 7.589, 7.470, 7.355, 7.244, 7.137,
     &  7.033,  6.934, 6.838, 6.746, 6.658, 6.574, 6.492, 6.413,
     &  6.337,  6.263, 6.190, 6.120, 6.050, 5.982, 5.915, 5.848,
     &  5.782,  5.718, 5.654, 5.591, 5.530, 5.470, 5.411, 5.354,
     &  5.298,  5.244, 5.193, 5.143, 5.095, 5.049, 5.005, 4.962,
     &  4.922,  4.884,  4.847,  4.812
     &/

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c...k- p total
      data (sigfit(11,i),i=1,itblsz)/
     & 111.295, 113.049,88.118,73.940,77.140,59.757,48.483,42.075,
     & 38.023, 35.747,35.030,35.806,37.788,39.555,40.985,42.810,
     & 44.977, 46.829,47.124,45.053,41.529,38.054,35.266,33.371,
     & 32.370, 32.078,32.222,32.539,32.822,32.935,32.811,32.462,
     & 31.961, 31.421,30.929,30.523,30.196,29.922,29.669,29.408,
     & 29.121, 28.814,28.509,28.228,27.977,27.753,27.551,27.371,
     & 27.210, 27.066,26.935,26.813,26.697,26.585,26.476,26.368,
     & 26.263, 26.159,26.056,25.955,25.855,25.757,25.661,25.566,
     & 25.472, 25.380,25.290,25.201,25.113,25.026,24.940,24.855,
     & 24.772, 24.689,24.607,24.526,24.445,24.365,24.285,24.206,
     & 24.127, 24.049,23.970,23.892,23.815,23.737,23.659,23.582,
     & 23.505, 23.427,23.350,23.273,23.197,23.120,23.044,22.968,
     & 22.892, 22.818, 22.743, 22.670
     &/

c...k- p elastic
      data (sigfit(12,i),i=1,itblsz)/
     & 84.158, 59.910,46.330,38.086,33.725,29.278,25.318,22.122,
     & 19.527, 17.513,16.215,15.940,17.039,18.872,20.100,20.634,
     & 20.823, 21.050,20.487,18.981,17.018,15.022,13.239,11.835,
     & 10.784,  9.995, 9.419, 9.031, 8.790, 8.643, 8.527, 8.393,
     &  8.223,  8.019, 7.800, 7.584, 7.381, 7.191, 7.010, 6.832,
     &  6.657,  6.485, 6.318, 6.157, 6.003, 5.857, 5.719, 5.589,
     &  5.467,  5.352, 5.246, 5.147, 5.055, 4.969, 4.890, 4.815,
     &  4.746,  4.681, 4.621, 4.564, 4.510, 4.459, 4.409, 4.364,
     &  4.319,  4.277, 4.235, 4.197, 4.158, 4.122, 4.087, 4.052,
     &  4.019,  3.986, 3.955, 3.924, 3.894, 3.865, 3.836, 3.808,
     &  3.781,  3.754, 3.728, 3.702, 3.677, 3.652, 3.628, 3.604,
     &  3.581,  3.558, 3.535, 3.513, 3.491, 3.470, 3.449, 3.428,
     &  3.408,  3.388,  3.368,  3.348
     &/

c...k- n total
      data (sigfit(13,i),i=1,itblsz)/
     & 22.954, 23.315,23.677,24.039,24.401,24.763,25.124,25.486,
     & 25.848, 26.248,27.828,28.982,29.933,30.831,31.681,32.386,
     & 32.850, 33.033,32.967,32.707,32.318,31.844,31.322,30.777,
     & 30.229, 29.686,29.160,28.654,28.170,27.713,27.282,26.876,
     & 26.496, 26.141,25.809,25.500,25.212,24.945,24.696,24.464,
     & 24.249, 24.049,23.862,23.688,23.526,23.375,23.234,23.102,
     & 22.979, 22.863,22.755,22.653,22.557,22.467,22.381,22.301,
     & 22.225, 22.153,22.084,22.019,21.957,21.898,21.841,21.787,
     & 21.736, 21.686,21.638,21.592,21.548,21.505,21.464,21.424,
     & 21.385, 21.347,21.310,21.274,21.240,21.206,21.172,21.140,
     & 21.108, 21.076,21.046,21.015,20.986,20.956,20.927,20.899,
     & 20.871, 20.843,20.815,20.788,20.761,20.735,20.708,20.682,
     & 20.656, 20.630, 20.605, 20.579
     &/

c...k- n elastic
      data (sigfit(14,i),i=1,itblsz)/
     &  4.216,  4.282, 4.349, 4.415, 4.482, 4.548, 4.614, 4.681,
     &  4.747,  6.379, 9.021,10.028,10.927,12.414,14.564,16.319,
     & 16.985, 16.860,16.419,15.873,15.305,14.744,14.201,13.682,
     & 13.185, 12.710,12.256,11.820,11.403,11.001,10.614,10.241,
     &  9.880,  9.531, 9.193, 8.865, 8.547, 8.237, 7.937, 7.645,
     &  7.360,  7.084, 6.815, 6.553, 6.298, 6.050, 5.809, 5.576,
     &  5.349,  5.130, 4.918, 4.714, 4.517, 4.329, 4.149, 3.979,
     &  3.818,  3.669, 3.532, 3.409, 3.301, 3.212, 3.145, 3.104,
     &  3.000,  3.000, 3.000, 3.000, 3.000, 3.000, 3.000, 3.000,
     &  3.000,  3.000, 3.000, 3.000, 3.000, 3.000, 3.000, 3.000,
     &  3.000,  3.000, 3.000, 3.000, 3.000, 3.000, 3.000, 3.000,
     &  3.000,  3.000, 3.000, 3.000, 3.000, 3.000, 3.000, 3.000,
     &  3.000,  3.000,  3.000,  3.000
     &/

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c...k+ p total
      data (sigfit(15,i),i=1,itblsz)/
     &  9.179, 11.411,11.716,12.237,12.640,12.840,12.879,12.835,
     & 12.751, 12.582,12.600,12.610,12.634,13.100,13.689,14.368,
     & 15.109, 15.854,16.535,17.107,17.545,17.845,18.022,18.103,
     & 18.113, 18.077,18.017,17.950,17.887,17.834,17.790,17.752,
     & 17.714, 17.676,17.636,17.597,17.562,17.531,17.504,17.481,
     & 17.458, 17.436,17.413,17.387,17.361,17.337,17.318,17.306,
     & 17.303, 17.306,17.315,17.329,17.345,17.363,17.383,17.403,
     & 17.424, 17.446,17.468,17.489,17.511,17.532,17.553,17.573,
     & 17.593, 17.613,17.632,17.650,17.668,17.685,17.702,17.718,
     & 17.733, 17.748,17.762,17.776,17.789,17.801,17.812,17.823,
     & 17.833, 17.843,17.852,17.860,17.867,17.874,17.879,17.885,
     & 17.889, 17.892,17.895,17.897,17.897,17.897,17.896,17.894,
     & 17.891, 17.887, 17.882, 17.876
     &/

c...k+ p elastic
      data (sigfit(16,i),i=1,itblsz)/
     &  9.179, 11.411, 11.716, 12.237, 12.640, 12.840, 12.879, 12.350, 
     & 12.387, 12.390, 12.392, 12.415, 12.446, 12.424, 12.301, 12.092, 
     & 11.847, 11.599, 11.367, 11.165, 10.998, 10.855, 10.694, 10.478, 
     & 10.218,  9.945,  9.673,  9.398,  9.109,  8.798,  8.471,  8.135, 
     &  7.803,  7.488,  7.197,  6.937,  6.708,  6.505,  6.325,  6.163, 
     &  6.014,  5.876,  5.745,  5.622,  5.504,  5.392,  5.285,  5.182, 
     &  5.083,  4.989,  4.900,  4.815,  4.734,  4.656,  4.583,  4.514, 
     &  4.448,  4.385,  4.326,  4.269,  4.216,  4.166,  4.118,  4.074, 
     &  4.031,  3.991,  3.954,  3.919,  3.886,  3.855,  3.826,  3.798, 
     &  3.772,  3.748,  3.726,  3.704,  3.684,  3.665,  3.647,  3.630, 
     &  3.613,  3.598,  3.583,  3.569,  3.555,  3.542,  3.529,  3.516, 
     &  3.504,  3.492,  3.480,  3.468,  3.457,  3.445,  3.434,  3.422, 
     &  3.411,  3.399,  3.388,  3.376
     &/

c...k+ n total
      data (sigfit(17,i),i=1,itblsz)/
     & 13.012, 13.217,13.422,13.627,13.832,14.037,14.242,14.448,
     & 14.653, 14.858,15.063,15.268,15.473,15.966,16.522,17.109,
     & 17.758, 18.446,19.052,19.478,19.697,19.724,19.613,19.431,
     & 19.240, 19.078,18.959,18.876,18.820,18.780,18.749,18.721,
     & 18.691, 18.659,18.623,18.577,18.516,18.443,18.363,18.284,
     & 18.209, 18.140,18.077,18.023,17.977,17.937,17.901,17.870,
     & 17.843, 17.820,17.801,17.784,17.770,17.757,17.746,17.735,
     & 17.726, 17.718,17.710,17.703,17.696,17.690,17.685,17.679,
     & 17.675, 17.670,17.666,17.662,17.659,17.656,17.653,17.650,
     & 17.648, 17.646,17.644,17.643,17.641,17.640,17.639,17.638,
     & 17.638, 17.637,17.637,17.637,17.637,17.638,17.638,17.639,
     & 17.640, 17.641,17.642,17.644,17.646,17.648,17.650,17.652,
     & 17.655, 17.658, 17.661, 17.664
     &/

c...k+ n elastic (blow pion production not appicable!)
      data (sigfit(18,i),i=1,itblsz)/
     & 13.012, 13.217, 13.422, 13.627, 13.832, 14.037, 14.242, 14.150, 
     & 13.847, 13.590, 12.892, 12.615, 12.446, 12.424, 12.301, 12.092, 
     & 11.847, 11.599, 11.367, 11.165, 10.998, 10.855, 10.694, 10.478, 
     & 10.218,  9.945,  9.673,  9.398,  9.109,  8.798,  8.471,  8.135, 
     &  7.803,  7.488,  7.197,  6.937,  6.708,  6.505,  6.325,  6.163, 
     &  6.014,  5.876,  5.745,  5.622,  5.504,  5.392,  5.285,  5.182, 
     &  5.083,  4.989,  4.900,  4.815,  4.734,  4.656,  4.583,  4.514, 
     &  4.448,  4.385,  4.326,  4.269,  4.216,  4.166,  4.118,  4.074, 
     &  4.031,  3.991,  3.954,  3.919,  3.886,  3.855,  3.826,  3.798, 
     &  3.772,  3.748,  3.726,  3.704,  3.684,  3.665,  3.647,  3.630, 
     &  3.613,  3.598,  3.583,  3.569,  3.555,  3.542,  3.529,  3.516, 
     &  3.504,  3.492,  3.480,  3.468,  3.457,  3.445,  3.434,  3.422, 
     &  3.411,  3.399,  3.388,  3.376
     &/

c...k- p -> k0 n charge exchange cross section
      data (sigfit(19,i),i=1,itblsz)/
     & 22.591, 12.394, 8.021, 7.641, 9.301, 6.995, 5.166, 4.362,
     &  3.969,  3.662, 3.300, 2.985, 3.234, 3.929, 4.448, 4.916,
     &  5.602,  6.522, 6.977, 6.299, 4.928, 3.610, 2.644, 2.081,
     &  1.867,  1.836, 1.852, 1.853, 1.822, 1.753, 1.651, 1.524,
     &  1.385,  1.250, 1.133, 1.042, 0.975, 0.926, 0.886, 0.846,
     &  0.802,  0.753, 0.702, 0.652, 0.609, 0.571, 0.538, 0.510,
     &  0.484,  0.461, 0.439, 0.419, 0.400, 0.381, 0.363, 0.346,
     &  0.329,  0.313, 0.298, 0.283, 0.268, 0.254, 0.241, 0.229,
     &  0.217,  0.206, 0.195, 0.186, 0.177, 0.169, 0.161, 0.154,
     &  0.147,  0.141, 0.136, 0.130, 0.125, 0.121, 0.116, 0.112,
     &  0.108,  0.104, 0.101, 0.097, 0.094, 0.091, 0.088, 0.085,
     &  0.082,  0.079, 0.077, 0.074, 0.072, 0.070, 0.067, 0.065,
     &  0.063,  0.061,  0.059,  0.057
     &/

c... k+ n -> k0 p charge exchange cross section
      data (sigfit(20,i),i=1,itblsz)/
     &  1.025,  1.041, 2.002, 2.868, 3.621, 4.347, 5.010, 5.569,
     &  6.011,  6.321, 6.536, 6.660, 6.720, 6.720, 6.674, 6.577,
     &  6.435,  6.248, 6.013, 5.739, 5.431, 5.100, 4.757, 4.418,
     &  4.094,  3.793, 3.518, 3.268, 3.040, 2.830, 2.636, 2.456,
     &  2.289,  2.134, 1.990, 1.857, 1.734, 1.620, 1.516, 1.420,
     &  1.332,  1.250, 1.174, 1.103, 1.037, 0.976, 0.919, 0.866,
     &  0.817,  0.770, 0.727, 0.687, 0.649, 0.613, 0.580, 0.548,
     &  0.518,  0.489, 0.463, 0.437, 0.414, 0.392, 0.372, 0.352,
     &  0.334,  0.317, 0.301, 0.285, 0.270, 0.257, 0.244, 0.231,
     &  0.220,  0.209, 0.199, 0.189, 0.180, 0.171, 0.163, 0.155,
     &  0.147,  0.140, 0.133, 0.127, 0.121, 0.115, 0.110, 0.105,
     &  0.100,  0.095, 0.090, 0.086, 0.082, 0.078, 0.075, 0.072,
     &  0.068,  0.065,  0.062,  0.059
     &/

c...k- p => lambda pi0                                        
      data (sigfit(21,i),i=1,itblsz)/
     & 10.005, 10.162, 6.746, 5.132, 4.229, 3.635, 3.190, 2.887,
     &  2.697,  2.617, 2.631, 2.689, 2.870, 2.812, 2.842, 3.130,
     &  3.232,  2.999, 2.563, 2.097, 1.722, 1.478, 1.312, 1.222,
     &  1.209,  1.207, 1.185, 1.141, 1.074, 0.987, 0.890, 0.797,
     &  0.720,  0.659, 0.609, 0.567, 0.530, 0.496, 0.466, 0.437,
     &  0.411,  0.386, 0.363, 0.340, 0.319, 0.299, 0.280, 0.261,
     &  0.244,  0.227, 0.211, 0.197, 0.183, 0.169, 0.157, 0.145,
     &  0.135,  0.125, 0.116, 0.107, 0.100, 0.093, 0.087, 0.081,
     &  0.076,  0.071, 0.067, 0.064, 0.060, 0.057, 0.054, 0.051,
     &  0.049,  0.047, 0.044, 0.042, 0.040, 0.039, 0.037, 0.035,
     &  0.034,  0.032, 0.031, 0.030, 0.028, 0.027, 0.026, 0.025,
     &  0.024,  0.023, 0.022, 0.021, 0.021, 0.020, 0.019, 0.018,
     &  0.018,  0.017,  0.017,  0.016
     &/

c...k- p => sigma- pi+                                       
      data (sigfit(22,i),i=1,itblsz)/
     &  8.576,  8.711, 9.213, 6.917, 6.911, 5.448, 4.134, 3.319,
     &  2.903,  2.839, 3.011, 3.089, 2.859, 2.445, 2.025, 1.690,
     &  1.450,  1.280, 1.146, 1.027, 0.907, 0.781, 0.651, 0.525,
     &  0.429,  0.387, 0.370, 0.355, 0.336, 0.313, 0.289, 0.264,
     &  0.239,  0.215, 0.192, 0.172, 0.154, 0.139, 0.125, 0.113,
     &  0.103,  0.094, 0.086, 0.079, 0.072, 0.066, 0.059, 0.054,
     &  0.048,  0.043, 0.038, 0.033, 0.029, 0.024, 0.021, 0.018,
     &  0.015,  0.013, 0.011, 0.010, 0.009, 0.008, 0.008, 0.007,
     &  0.007,  0.006, 0.005, 0.005, 0.004, 0.003, 0.003, 0.000,
     &  0.000,  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000,  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000,  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000,  0.000,  0.000,  0.000
     &/

c...k- p => sigma+ pi-                                        
      data (sigfit(23,i),i=1,itblsz)/
     & 36.056, 16.855,11.691,11.018,10.679, 8.768, 7.060, 5.816,
     &  4.843,  4.064, 3.443, 2.961, 2.585, 2.259, 1.962, 1.747,
     &  1.672,  1.681, 1.654, 1.550, 1.425, 1.374, 1.374, 1.350,
     &  1.282,  1.207, 1.152, 1.112, 1.072, 1.021, 0.960, 0.891,
     &  0.821,  0.755, 0.694, 0.640, 0.594, 0.554, 0.519, 0.489,
     &  0.463,  0.440, 0.418, 0.399, 0.380, 0.362, 0.345, 0.328,
     &  0.312,  0.296, 0.281, 0.266, 0.251, 0.238, 0.224, 0.212,
     &  0.200,  0.189, 0.179, 0.170, 0.161, 0.153, 0.146, 0.139,
     &  0.133,  0.128, 0.122, 0.118, 0.113, 0.109, 0.105, 0.101,
     &  0.098,  0.094, 0.091, 0.088, 0.085, 0.082, 0.079, 0.077,
     &  0.074,  0.072, 0.069, 0.067, 0.065, 0.062, 0.060, 0.058,
     &  0.056,  0.054, 0.052, 0.050, 0.048, 0.047, 0.045, 0.043,
     &  0.041,  0.040,  0.038,  0.037
     &/

c...Lambda p total
      data (sigfit(24,i),i=1,itblsz)/
     &  0.000,  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 114.450,26.369,12.922,17.513,23.050,24.958,25.415,
     & 24.205, 23.916,24.695,25.546,26.238,26.812,27.309,27.752,
     & 28.161, 28.547,28.916,29.274,29.623,29.958,30.288,30.606,
     & 30.917, 31.212,31.507,31.776,32.042,32.307,32.537,32.765,
     & 32.993, 33.202,33.385,33.569,33.753,33.930,34.063,34.197,
     & 34.330, 34.463,34.597,34.687,34.768,34.849,34.930,35.011,
     & 35.092, 35.159,35.190,35.222,35.253,35.285,35.316,35.347,
     & 35.379, 35.403,35.394,35.386,35.377,35.369,35.360,35.352,
     & 35.343, 35.335,35.326,35.313,35.283,35.253,35.223,35.193,
     & 35.162, 35.132,35.102,35.072,35.041,35.011,34.981,34.951,
     & 34.927, 34.904,34.881,34.858,34.836,34.813,34.790,34.767,
     & 34.744, 34.722, 34.699, 34.676
     &/

c...Lambda p elastic
      data (sigfit(25,i),i=1,itblsz)/
     &  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, 
     &  0.000, 114.450, 26.369, 12.922, 17.513, 23.050, 24.958, 25.415, 
     & 13.494, 13.452, 13.417, 13.423, 13.473, 13.554, 13.643, 13.718, 
     & 13.765, 13.782, 13.770, 13.737, 13.687, 13.625, 13.554, 13.476, 
     & 13.391, 13.300, 13.202, 13.097, 12.986, 12.867, 12.742, 12.610, 
     & 12.472, 12.327, 12.177, 12.023, 11.864, 11.702, 11.538, 11.373, 
     & 11.208, 11.043, 10.879, 10.717, 10.558, 10.401, 10.249, 10.100, 
     &  9.953,  9.812,  9.675,  9.540,  9.411,  9.283,  9.161,  9.040, 
     &  8.924,  8.810,  8.700,  8.590,  8.485,  8.380,  8.278,  8.178, 
     &  8.078,  7.982,  7.886,  7.792,  7.700,  7.607,  7.516,  7.427, 
     &  7.338,  7.251,  7.165,  7.079,  6.993,  6.911,  6.829,  6.748, 
     &  6.668,  6.591,  6.514,  6.437,  6.365,  6.296,  6.227,  6.157, 
     &  6.098,  6.040,  5.983,  5.926
     &/

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c... one pion production cross sections up to 2 gev
c ref.  g.shimizu et.al n.p.a386(1982)571
c       b.j.verwest and r.a.arndt, p.r.c25(1982)1979
c note: appricability of sigma11,sigma10,sigma01
c       is blow 2 gev incident energy
c...pp => pp pi0  sigma11
       data (sigfit(26,i),i=1,itblsz)/
     &     .0,     .0,     .0,     .0,     .0,     .0,     .0,     .0,
     &    .04,    .12,    .29,    .62,   1.22,   2.11,   2.97,   3.56,
     & 3.9033, 4.1065, 4.2257, 4.2933, 4.3269, 4.3370, 4.3304, 4.3116,
     & 4.2838, 4.2494, 4.2100, 4.1669, 4.1212, 4.0735, 4.0245, 3.9747,
     & 3.9244, 3.8738, 3.8233, 3.7729, 3.7229, 3.6733, 3.6243, 3.5758,
     & 3.5280, 3.4809, 3.4345, 3.3888, 3.3439, 3.2997, 3.2562, 3.2135,
     & 3.1716, 3.1303, 3.0899, 3.0501, 3.0110, 2.9726, 2.9350, 2.8979,
     & 2.8616, 2.8258, 2.7908, 2.7563, 2.7224, 2.6891, 2.6564, 2.6242,
     & 2.5926, 2.5616, 2.5310, 2.5010, 2.4714, 2.4424, 2.4138, 2.3857,
     & 2.3580, 2.3308, 2.3040, 2.2777, 2.2517, 2.2262, 2.2010, 2.1763,
     & 2.1519, 2.1279, 2.1042, 2.0809, 2.0580, 2.0354, 2.0131, 1.9911,
     & 1.9695, 1.9481, 1.9271, 1.9064, 1.8859, 1.8658, 1.8459, 1.8263,
     & 1.807,  1.788,  1.769,  1.751
     &/

c...sigma10
       data (sigfit(27,i),i=1,itblsz)/
     &  .000,   .000,   .000,   .000,   .000,   .000,   .001,   .105,
     &  .327,   .694,  1.283,  2.229,  3.702,  5.733,  7.990,  9.989,
     & 11.507, 12.565, 13.254, 13.668, 13.878, 13.940, 13.892, 13.764,
     & 13.577, 13.347, 13.087, 12.806, 12.511, 12.207, 11.898, 11.587,
     & 11.278, 10.971, 10.667, 10.369, 10.077,  9.792,  9.513,  9.241,
     &  8.976,  8.719,  8.469,  8.226,  7.991,  7.762,  7.541,  7.326,
     &  7.118,  6.917,  6.722,  6.533,  6.349,  6.172,  6.000,  5.834,
     &  5.672,  5.516,  5.365,  5.218,  5.076,  4.938,  4.804,  4.675,
     &  4.549,  4.427,  4.309,  4.195,  4.084,  3.976,  3.872,  3.770,
     &  3.672,  3.576,  3.483,  3.393,  3.306,  3.221,  3.138,  3.058,
     &  2.980,  2.905,  2.831,  2.760,  2.690,  2.623,  2.557,  2.493,
     &  2.431,  2.371,  2.312,  2.255,  2.199,  2.145,  2.092,  2.041,
     &  1.991,  1.943,  1.895,  1.849
     &/

c...sigma01
       data (sigfit(28,i),i=1,itblsz)/
     &   .000,   .000,   .000,   .000,   .000,   .000,   .000,   .002,
     &   .004,   .008,   .013,   .019,   .027,   .037,   .049,   .066,
     &   .088,   .118,   .161,   .222,   .310,   .437,   .616,   .852,
     &  1.141,  1.473,  1.837,  2.224,  2.630,  3.055,  3.499,  3.963,
     &  4.451,  4.965,  5.509,  6.088,  6.704,  7.365,  8.074,  8.838,
     &  9.663, 10.556, 11.525, 12.577, 13.722, 14.966, 16.317, 17.784,
     & 19.369, 21.077, 22.903, 24.840, 26.869, 28.961, 31.074, 33.149,
     & 35.114, 36.883, 38.365, 39.472, 40.130, 40.291, 39.942, 39.106,
     & 37.839, 36.221, 34.345, 32.302, 30.178, 28.043, 25.952, 23.945,
     & 22.048, 20.276, 18.634, 17.124, 15.740, 14.478, 13.329, 12.284,
     & 11.334, 10.471,  9.686,  8.973,  8.323,  7.732,  7.192,  6.699,
     &  6.247,  5.834,  5.455,  5.107,  4.787,  4.492,  4.219,  3.968,
     &  3.735,  3.520,  3.320,  3.134
     &/

c...ref. Gy.wolf et al., nucl.phys. A517(1990)615
c...     Gy.wolf et al., nucl.phys. A552(1993)549
c...     A.engel et al., nucl.phys. A572(1994)657
c... delta absorption coefficient

       data (sigfit(29,i),i=1,itblsz)/
     &    .0,     .0,     .0,     .0,     .0,     .0,     .0,  2365.0,
     & 296.99,  80.77,  30.27,  13.25,   6.51,   3.71,   2.54,   2.02,
     & 1.7485, 1.5944, 1.4965, 1.4296, 1.3814, 1.3453, 1.3174, 1.2953,
     & 1.2775, 1.2630, 1.2509, 1.2407, 1.2322, 1.2248, 1.2185, 1.2131,
     & 1.2083, 1.2042, 1.2005, 1.1973, 1.1945, 1.1919, 1.1897, 1.1876,
     & 1.1858, 1.1842, 1.1827, 1.1814, 1.1802, 1.1791, 1.1781, 1.1772,
     & 1.1764, 1.1756, 1.1749, 1.1743, 1.1737, 1.1732, 1.1727, 1.1722,
     & 1.1718, 1.1714, 1.1711, 1.1707, 1.1704, 1.1701, 1.1699, 1.1696,
     & 1.1694, 1.1692, 1.1690, 1.1688, 1.1686, 1.1685, 1.1683, 1.1682,
     & 1.1681, 1.1680, 1.1678, 1.1677, 1.1676, 1.1675, 1.1674, 1.1674,
     & 1.1673, 1.1672, 1.1672, 1.1671, 1.1670, 1.1670, 1.1669, 1.1669,
     & 1.1668, 1.1668, 1.1667, 1.1667, 1.1667, 1.1666, 1.1666, 1.1666,
     & 1.167,  1.166,  1.166,  1.164
     &/

      end

c***********************************************************************

      subroutine jamsigS1(sig,isig,s)

c...S=-1 low energy Baryon-Baryon cross sections below particle production
c...from Nijmegen model D

      implicit double precision(a-h, o-z)
c...Table parameters.
      parameter (itblsz=100)
      parameter (tblhig=2.51227903)
      parameter (tbllow=2.05588)    ! lambda p threthold + 1MeV
      parameter (tblstp=0.002)
      parameter (isgmax=26)
      dimension ems(isgmax)
      common /stbl1/sfits1(isgmax,itblsz)
      data ems/
     a  1.116,1.116,1.189,1.192,
     a  1.116,1.116,1.192,1.197,
     a  1.197,1.189,
     a  1.197,1.197,1.116,1.192,
     a  1.192,1.192,1.116,1.197,
     a  1.192,1.192,1.116,1.189,
     a  1.189,1.189,1.116,1.192/

      if(s.le.ems(isig)+0.938+0.001) then
        sig=0.0
        return
      endif

c...1.  lambda p total
c...2.  lambda p elastic
c...3.  lambda p => sigma+ n
c...4.  lambda p => sigma0 p
c...5.  lambda n total
c...6.  lambda n elastic
c...7.  lambda n => sigma0 n
c...8.  lambda n => sigma- p
c...9.  sigma- n total
c...10. sigma+ p total
c...11. sigma- p total
c...12. sigma- p elastic
c...13. sigma- p => lambda n
c...14. sigma- p => sigma0 n
c...15. sigma0 n total
c...16. sigma0 n elastic
c...17. sigma0 n => lambda n
c...18. sigma0 n => sigma- p
c...19. sigma0 p total
c...20. sigma0 p elastic
c...21. sigma0 p => lambda p
c...22. sigma0 p => sigma+ n
c...23. sigma+ n total
c...24. sigma+ n elastic
c...25. sigma+ n => lambda p
c...26. sigma+ n => sigma0 p

      xpt=log(s)
      tlow=log(tbllow)
      if(xpt.lt.tlow) xpt=tlow
      index=int((xpt-tlow)/tblstp)+1
      if(index.ge.itblsz) then
         sig=sfits1(isig,itblsz)
         return
      end if
c  find slopes and crossections
      x1=(index-1)*tblstp+tlow
      y1in=sfits1(isig,index)
      y2in=sfits1(isig,index+1)
      slin=(y2in-y1in)/tblstp
      sig=slin*(xpt-x1)+y1in
      if(sig.lt.0.0d0) then
         sig=0d0
         write(6,*)'(sigs1) sig index isig',sig,index,isig
         write(6,*)'y1in y2in',y1in,y2in
         write(6,*)'slin',slin
         write(6,*)'slin',slin
      endif

      end

c*************************************************************************

      block data jamss1da

c...S=-1 baryon - baryon cross sections extracted from njimegen model D.
      implicit double precision(a-h, o-z)
      parameter (itblsz=100)
      parameter (isgmax=26)

      common /stbl1/sfits1(isgmax,itblsz)

c...Lambda p total
      data (sfits1(1,i),i=1,itblsz)/
     $ 301.896,155.724,96.181,64.964,46.451,34.660,26.811,21.438,
     $ 17.727,15.138,13.353,12.168,11.451,11.131,11.192,11.721,
     $ 13.093,17.873,24.639,20.698,20.935,21.774,22.630,23.348,
     $ 23.901,24.298,24.566,24.730,24.814,24.837,24.815,24.762,
     $ 24.687,24.599,24.503,24.404,24.305,24.208,24.117,24.032,
     $ 23.954,23.886,23.827,23.778,23.741,23.716,23.704,23.704,
     $ 23.718,23.747,23.789,23.844,23.914,23.996,24.091,24.197,
     $ 24.313,24.437,24.569,24.706,24.846,24.988,25.131,25.272,
     $ 25.410,25.545,25.674,25.798,25.916,26.027,26.132,26.231,
     $ 26.323,26.409,26.490,26.566,26.637,26.705,26.769,26.831,
     $ 26.890,26.948,27.004,27.060,27.116,27.171,27.227,27.283,
     $ 27.341,27.399,27.459,27.521,27.585,27.650,27.717,27.787,
     $ 27.859, 27.933, 28.009, 28.088
     $/
c...Lambda p elastic
      data (sfits1(2,i),i=1,itblsz)/
     $ 301.896,155.724,96.181,64.964,46.451,34.660,26.811,21.438,
     $ 17.727,15.138,13.353,12.168,11.451,11.131,11.192,11.721,
     $ 13.093,17.873,17.195,12.808,12.263,12.256,12.393,12.556,
     $ 12.705,12.830,12.935,13.026,13.111,13.194,13.278,13.366,
     $ 13.459,13.558,13.663,13.773,13.889,14.010,14.135,14.265,
     $ 14.399,14.536,14.678,14.824,14.974,15.129,15.287,15.450,
     $ 15.617,15.789,15.966,16.146,16.330,16.518,16.709,16.902,
     $ 17.097,17.293,17.488,17.683,17.875,18.065,18.252,18.434,
     $ 18.612,18.785,18.952,19.114,19.270,19.421,19.567,19.708,
     $ 19.844,19.976,20.104,20.229,20.352,20.471,20.589,20.706,
     $ 20.821,20.935,21.049,21.162,21.276,21.390,21.505,21.620,
     $ 21.736,21.853,21.972,22.092,22.213,22.335,22.459,22.585,
     $ 22.712, 22.840, 22.971, 23.103
     $/
c...Lambda p => Sigma+ n
      data (sfits1(3,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 6.396, 5.597, 6.049, 6.580, 7.020, 7.347,
     $  7.576, 7.725, 7.812, 7.845, 7.834, 7.787, 7.712, 7.614,
     $  7.500, 7.374, 7.239, 7.099, 6.955, 6.810, 6.665, 6.521,
     $  6.380, 6.242, 6.107, 5.978, 5.852, 5.733, 5.618, 5.510,
     $  5.407, 5.312, 5.222, 5.139, 5.062, 4.992, 4.927, 4.869,
     $  4.816, 4.769, 4.726, 4.688, 4.653, 4.621, 4.591, 4.563,
     $  4.537, 4.511, 4.486, 4.461, 4.435, 4.408, 4.381, 4.353,
     $  4.323, 4.293, 4.261, 4.228, 4.194, 4.159, 4.123, 4.087,
     $  4.049, 4.012, 3.973, 3.935, 3.896, 3.856, 3.817, 3.778,
     $  3.739, 3.699, 3.660, 3.622, 3.583, 3.545, 3.507, 3.470,
     $  3.433,  3.396,  3.360,  3.325
     $/
c...Lambda p => sigma0 p
      data (sfits1(4,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 1.048, 2.292, 2.623, 2.938, 3.216, 3.445,
     $  3.620, 3.743, 3.820, 3.859, 3.869, 3.856, 3.825, 3.781,
     $  3.727, 3.666, 3.601, 3.532, 3.461, 3.389, 3.317, 3.246,
     $  3.176, 3.108, 3.041, 2.977, 2.914, 2.855, 2.798, 2.744,
     $  2.693, 2.646, 2.601, 2.560, 2.522, 2.487, 2.455, 2.426,
     $  2.400, 2.376, 2.355, 2.336, 2.318, 2.303, 2.288, 2.274,
     $  2.261, 2.249, 2.236, 2.224, 2.211, 2.198, 2.184, 2.170,
     $  2.156, 2.141, 2.125, 2.109, 2.092, 2.075, 2.057, 2.039,
     $  2.020, 2.001, 1.982, 1.963, 1.944, 1.924, 1.905, 1.885,
     $  1.866, 1.847, 1.827, 1.808, 1.789, 1.770, 1.751, 1.732,
     $  1.714,  1.696,  1.678,  1.660
     $/
c...Lambda n total
      data (sfits1(5,i),i=1,itblsz)/
     $ 303.343,183.072,108.489,71.523,50.253,36.953,28.200,22.250,
     $ 18.151,15.284,13.286,11.923,11.040,10.538,10.361,10.500,
     $ 11.024,12.251,16.408,23.052,20.266,20.308,21.124,22.001,
     $ 22.753,23.343,23.781,24.089,24.289,24.405,24.457,24.461,
     $ 24.431,24.378,24.309,24.231,24.148,24.064,23.983,23.905,
     $ 23.833,23.768,23.712,23.665,23.628,23.603,23.589,23.588,
     $ 23.600,23.626,23.666,23.719,23.786,23.866,23.959,24.063,
     $ 24.178,24.302,24.434,24.572,24.714,24.859,25.005,25.149,
     $ 25.292,25.431,25.566,25.695,25.818,25.934,26.044,26.147,
     $ 26.244,26.335,26.420,26.499,26.574,26.645,26.712,26.776,
     $ 26.837,26.897,26.955,27.012,27.069,27.125,27.182,27.239,
     $ 27.297,27.356,27.416,27.478,27.542,27.607,27.674,27.744,
     $ 27.816, 27.889, 27.966, 28.044
     $/
c...Lambda n elastic
      data (sfits1(6,i),i=1,itblsz)/
     $ 303.343,183.072,108.489,71.523,50.253,36.953,28.200,22.250,
     $ 18.151,15.284,13.286,11.923,11.040,10.538,10.361,10.500,
     $ 11.024,12.251,15.894,17.488,12.953,12.275,12.255,12.402,
     $ 12.581,12.746,12.886,13.003,13.103,13.195,13.282,13.369,
     $ 13.459,13.553,13.652,13.757,13.867,13.982,14.101,14.226,
     $ 14.355,14.488,14.626,14.768,14.914,15.065,15.220,15.379,
     $ 15.544,15.713,15.886,16.064,16.247,16.433,16.623,16.815,
     $ 17.010,17.206,17.402,17.598,17.792,17.984,18.173,18.358,
     $ 18.539,18.715,18.885,19.050,19.210,19.364,19.512,19.655,
     $ 19.794,19.928,20.059,20.185,20.309,20.430,20.549,20.666,
     $ 20.782,20.896,21.011,21.124,21.238,21.352,21.467,21.582,
     $ 21.698,21.815,21.933,22.053,22.174,22.296,22.419,22.545,
     $ 22.671, 22.800, 22.930, 23.061
     $/
c...Lambda n => sigma0 n
      data (sfits1(7,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.514, 4.448, 3.094, 3.191, 3.404, 3.571,
     $  3.678, 3.743, 3.784, 3.806, 3.810, 3.799, 3.775, 3.738,
     $  3.692, 3.639, 3.580, 3.517, 3.451, 3.383, 3.314, 3.246,
     $  3.178, 3.111, 3.046, 2.982, 2.921, 2.862, 2.805, 2.751,
     $  2.700, 2.652, 2.607, 2.565, 2.526, 2.490, 2.457, 2.428,
     $  2.401, 2.377, 2.355, 2.336, 2.318, 2.302, 2.287, 2.274,
     $  2.261, 2.248, 2.236, 2.224, 2.211, 2.199, 2.186, 2.172,
     $  2.158, 2.143, 2.128, 2.112, 2.095, 2.078, 2.061, 2.043,
     $  2.025, 2.006, 1.987, 1.968, 1.949, 1.929, 1.910, 1.890,
     $  1.871, 1.851, 1.832, 1.812, 1.793, 1.774, 1.755, 1.737,
     $  1.718,  1.700,  1.682,  1.664
     $/
c...Lambda n => sigma- p
      data (sfits1(8,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 1.116, 4.219, 4.841, 5.464, 6.028,
     $  6.495, 6.854, 7.112, 7.280, 7.376, 7.411, 7.400, 7.353,
     $  7.279, 7.186, 7.077, 6.958, 6.831, 6.700, 6.567, 6.433,
     $  6.300, 6.169, 6.040, 5.915, 5.793, 5.677, 5.565, 5.458,
     $  5.357, 5.262, 5.173, 5.090, 5.013, 4.943, 4.879, 4.820,
     $  4.767, 4.720, 4.677, 4.639, 4.604, 4.573, 4.544, 4.517,
     $  4.492, 4.468, 4.444, 4.421, 4.397, 4.372, 4.346, 4.320,
     $  4.292, 4.263, 4.233, 4.202, 4.170, 4.137, 4.102, 4.067,
     $  4.031, 3.994, 3.957, 3.920, 3.882, 3.843, 3.805, 3.766,
     $  3.728, 3.689, 3.651, 3.613, 3.575, 3.537, 3.500, 3.463,
     $  3.426,  3.390,  3.354,  3.319
     $/
c...Sigma- n total
      data (sfits1(9,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000,212.730,173.436,82.977,71.805,79.624,
     $ 91.695,100.172,103.073,101.693,97.959,93.284,88.467,83.900,
     $ 79.718,75.953,72.586,69.580,66.902,64.498,62.332,60.372,
     $ 58.593,56.972,55.485,54.113,52.845,51.670,50.573,49.546,
     $ 48.583,47.676,46.820,46.008,45.239,44.508,43.811,43.146,
     $ 42.511,41.903,41.320,40.762,40.226,39.711,39.217,38.741,
     $ 38.284,37.844,37.421,37.013,36.621,36.244,35.881,35.532,
     $ 35.197,34.874,34.563,34.265,33.979,33.704,33.440,33.188,
     $ 32.946,32.714,32.492,32.281,32.079,31.886,31.702,31.528,
     $ 31.362,31.205,31.057,30.917,30.784,30.660,30.544,30.435,
     $ 30.334, 30.240, 30.153, 30.073
     $/
c...Sigma+ p total
      data (sfits1(10,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,137.722,78.144,71.989,81.003,92.699,100.529,
     $ 102.902,101.196,97.387,92.762,88.038,83.555,79.453,75.765,
     $ 72.462,69.508,66.861,64.484,62.348,60.412,58.650,57.039,
     $ 55.563,54.204,52.945,51.774,50.682,49.661,48.702,47.798,
     $ 46.944,46.136,45.369,44.639,43.944,43.281,42.647,42.039,
     $ 41.458,40.900,40.365,39.850,39.356,38.880,38.423,37.983,
     $ 37.560,37.152,36.760,36.382,36.018,35.669,35.332,35.008,
     $ 34.697,34.398,34.111,33.834,33.570,33.316,33.072,32.839,
     $ 32.615,32.402,32.198,32.004,31.819,31.642,31.475,31.316,
     $ 31.165,31.023,30.889,30.762,30.644,30.533,30.429,30.333,
     $ 30.243, 30.161, 30.086, 30.017
     $/

c..Sigma- p total
      data (sfits1(11,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,824.565,826.744,364.573,218.894,174.939,156.028,
     $ 145.346,137.225,129.927,123.126,116.806,110.975,105.625,100.724,
     $ 96.220,92.072,88.242,84.701,81.425,78.380,75.544,72.897,
     $ 70.426,68.114,65.945,63.907,61.989,60.187,58.486,56.879,
     $ 55.362,53.930,52.574,51.291,50.076,48.927,47.836,46.802,
     $ 45.823,44.892,44.008,43.169,42.370,41.609,40.885,40.193,
     $ 39.531,38.899,38.292,37.710,37.152,36.615,36.098,35.600,
     $ 35.120,34.657,34.211,33.781,33.366,32.965,32.579,32.207,
     $ 31.849,31.504,31.172,30.853,30.546,30.251,29.969,29.698,
     $ 29.439,29.191,28.955,28.729,28.514,28.310,28.116,27.932,
     $ 27.758, 27.594, 27.439, 27.294
     $/
c...Sigma- p elastic
      data (sfits1(12,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,199.113,199.640,111.729,71.260,58.100,53.681,
     $ 52.434,52.006,51.543,50.863,49.984,48.962,47.847,46.679,
     $ 45.487,44.292,43.107,41.944,40.812,39.714,38.655,37.634,
     $ 36.655,35.716,34.816,33.954,33.130,32.342,31.589,30.868,
     $ 30.180,29.524,28.897,28.298,27.727,27.184,26.665,26.171,
     $ 25.702,25.254,24.829,24.424,24.040,23.674,23.326,22.995,
     $ 22.681,22.382,22.097,21.825,21.567,21.321,21.086,20.862,
     $ 20.649,20.446,20.252,20.068,19.893,19.727,19.570,19.421,
     $ 19.280,19.147,19.023,18.905,18.796,18.694,18.600,18.513,
     $ 18.433,18.360,18.294,18.234,18.182,18.136,18.096,18.064,
     $ 18.037, 18.016, 18.002, 17.993
     $/
c...Sigma- p => lambda n
      data (sfits1(13,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,358.270,359.217,134.564,73.007,54.705,45.855,
     $ 40.393,36.463,33.332,30.703,28.432,26.439,24.676,23.110,
     $ 21.707,20.446,19.309,18.282,17.352,16.507,15.735,15.031,
     $ 14.387,13.797,13.255,12.756,12.297,11.875,11.485,11.126,
     $ 10.794,10.488,10.206, 9.946, 9.705, 9.483, 9.277, 9.087,
     $  8.909, 8.744, 8.589, 8.442, 8.304, 8.172, 8.045, 7.922,
     $  7.803, 7.687, 7.572, 7.460, 7.348, 7.238, 7.128, 7.019,
     $  6.911, 6.803, 6.696, 6.590, 6.485, 6.382, 6.279, 6.178,
     $  6.078, 5.980, 5.883, 5.788, 5.694, 5.603, 5.513, 5.425,
     $  5.339, 5.254, 5.172, 5.091, 5.012, 4.935, 4.860, 4.787,
     $  4.715,  4.645,  4.576,  4.510
     $/
c...Sigma- p => sigma0 n
      data (sfits1(14,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,267.181,267.888,118.281,74.628,62.135,56.492,
     $ 52.519,48.756,45.052,41.559,38.389,35.575,33.102,30.936,
     $ 29.026,27.334,25.826,24.475,23.261,22.159,21.154,20.232,
     $ 19.384,18.601,17.874,17.196,16.563,15.970,15.412,14.885,
     $ 14.388,13.918,13.471,13.047,12.644,12.260,11.893,11.544,
     $ 11.212,10.894,10.591,10.302,10.027, 9.764, 9.514, 9.275,
     $  9.047, 8.830, 8.623, 8.426, 8.237, 8.057, 7.884, 7.719,
     $  7.561, 7.409, 7.263, 7.122, 6.987, 6.856, 6.730, 6.609,
     $  6.491, 6.377, 6.267, 6.159, 6.055, 5.954, 5.856, 5.761,
     $  5.668, 5.577, 5.489, 5.404, 5.320, 5.239, 5.159, 5.082,
     $  5.007,  4.933,  4.861,  4.791
     $/

c...Sigma0 n total
      data (sfits1(15,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,564.416,314.411,163.623,124.689,116.850,117.768,
     $ 118.963,117.528,113.797,108.839,103.507,98.283,93.393,88.917,
     $ 84.843,81.146,77.786,74.729,71.945,69.393,67.045,64.878,
     $ 62.874,61.017,59.286,57.669,56.155,54.737,53.403,52.145,
     $ 50.959,49.839,48.778,47.773,46.821,45.917,45.057,44.239,
     $ 43.462,42.722,42.015,41.342,40.699,40.084,39.497,38.934,
     $ 38.395,37.878,37.383,36.906,36.448,36.008,35.584,35.177,
     $ 34.784,34.406,34.042,33.692,33.354,33.030,32.717,32.417,
     $ 32.129,31.852,31.586,31.331,31.088,30.854,30.631,30.418,
     $ 30.215,30.022,29.838,29.664,29.498,29.342,29.195,29.056,
     $ 28.926, 28.804, 28.690, 28.585
     $/
c...Sigma0 n elastic
      data (sfits1(16,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,149.300,118.977,53.088,40.775,43.258,49.751,
     $ 55.263,57.973,58.309,57.203,55.407,53.370,51.324,49.379,
     $ 47.570,45.903,44.374,42.973,41.689,40.507,39.417,38.407,
     $ 37.471,36.601,35.787,35.026,34.311,33.639,33.005,32.406,
     $ 31.840,31.304,30.794,30.311,29.851,29.414,28.997,28.600,
     $ 28.221,27.860,27.514,27.185,26.869,26.567,26.278,26.002,
     $ 25.737,25.483,25.239,25.006,24.782,24.568,24.362,24.166,
     $ 23.977,23.797,23.625,23.461,23.305,23.157,23.015,22.882,
     $ 22.756,22.637,22.525,22.421,22.323,22.232,22.148,22.071,
     $ 22.000,21.936,21.878,21.826,21.781,21.742,21.708,21.681,
     $ 21.659, 21.644, 21.633, 21.629
     $/
c...Sigma0 n => lambda n
      data (sfits1(17,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,415.116,173.864,49.766,32.536,25.856,21.996,
     $ 19.323,17.326,15.755,14.475,13.399,12.472,11.663,10.947,
     $ 10.308, 9.733, 9.214, 8.743, 8.316, 7.926, 7.569, 7.242,
     $  6.942, 6.667, 6.414, 6.180, 5.964, 5.765, 5.582, 5.412,
     $  5.255, 5.110, 4.977, 4.853, 4.739, 4.634, 4.537, 4.446,
     $  4.362, 4.283, 4.210, 4.140, 4.074, 4.012, 3.951, 3.893,
     $  3.836, 3.780, 3.725, 3.671, 3.618, 3.564, 3.512, 3.459,
     $  3.407, 3.354, 3.303, 3.251, 3.200, 3.150, 3.100, 3.050,
     $  3.002, 2.954, 2.906, 2.860, 2.814, 2.769, 2.725, 2.682,
     $  2.640, 2.599, 2.558, 2.519, 2.480, 2.442, 2.405, 2.369,
     $  2.334,  2.299,  2.266,  2.233
     $/
c...Sigma0 n => sigma- p
      data (sfits1(18,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000,21.570,60.768,51.378,47.736,46.021,
     $ 44.377,42.230,39.733,37.161,34.701,32.441,30.406,28.590,
     $ 26.965,25.509,24.198,23.014,21.940,20.959,20.059,19.228,
     $ 18.461,17.749,17.085,16.463,15.880,15.333,14.816,14.327,
     $ 13.864,13.425,13.007,12.609,12.230,11.869,11.523,11.194,
     $ 10.879,10.578,10.291,10.017, 9.755, 9.505, 9.267, 9.040,
     $  8.823, 8.615, 8.418, 8.229, 8.048, 7.876, 7.710, 7.552,
     $  7.400, 7.254, 7.114, 6.979, 6.849, 6.723, 6.602, 6.485,
     $  6.371, 6.261, 6.155, 6.051, 5.950, 5.853, 5.758, 5.665,
     $  5.575, 5.488, 5.402, 5.319, 5.238, 5.159, 5.082, 5.006,
     $  4.933,  4.861,  4.791,  4.723
     $/

c...Sigma0 p total
      data (sfits1(19,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,1780.569,218.820,142.584,123.241,119.973,120.363,
     $ 119.411,116.195,111.539,106.301,101.038,96.031,91.399,87.172,
     $ 83.323,79.820,76.627,73.712,71.048,68.598,66.338,64.248,
     $ 62.310,60.510,58.829,57.256,55.782,54.398,53.095,51.866,
     $ 50.704,49.607,48.568,47.582,46.647,45.759,44.914,44.110,
     $ 43.345,42.616,41.920,41.256,40.622,40.015,39.435,38.879,
     $ 38.345,37.834,37.343,36.871,36.417,35.980,35.560,35.155,
     $ 34.765,34.390,34.028,33.680,33.345,33.022,32.712,32.414,
     $ 32.127,31.852,31.588,31.335,31.093,30.861,30.639,30.428,
     $ 30.226,30.035,29.852,29.679,29.515,29.360,29.214,29.077,
     $ 28.948, 28.827, 28.715, 28.610
     $/
c...Sigma0 p elastic
      data (sfits1(20,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,200.242,57.972,38.734,37.845,43.299,49.512,
     $ 53.682,55.359,55.291,54.230,52.697,50.992,49.273,47.619,
     $ 46.059,44.603,43.250,41.996,40.835,39.757,38.755,37.820,
     $ 36.948,36.133,35.368,34.648,33.970,33.331,32.726,32.153,
     $ 31.610,31.094,30.604,30.138,29.694,29.271,28.868,28.483,
     $ 28.115,27.764,27.428,27.107,26.799,26.505,26.222,25.952,
     $ 25.692,25.444,25.205,24.976,24.756,24.545,24.343,24.149,
     $ 23.964,23.786,23.617,23.455,23.301,23.155,23.016,22.884,
     $ 22.760,22.642,22.532,22.429,22.333,22.243,22.161,22.085,
     $ 22.015,21.952,21.895,21.845,21.801,21.762,21.730,21.704,
     $ 21.683, 21.668, 21.659, 21.655
     $/
c...Sigma0 p => lambda p
      data (sfits1(21,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,683.356,58.465,34.701,26.639,22.533,19.927,
     $ 18.025,16.506,15.222,14.107,13.127,12.257,11.482,10.790,
     $ 10.168, 9.606, 9.097, 8.636, 8.217, 7.835, 7.485, 7.165,
     $  6.872, 6.602, 6.354, 6.126, 5.915, 5.720, 5.541, 5.375,
     $  5.222, 5.081, 4.951, 4.831, 4.719, 4.617, 4.521, 4.433,
     $  4.350, 4.273, 4.201, 4.132, 4.067, 4.005, 3.944, 3.886,
     $  3.829, 3.773, 3.718, 3.664, 3.610, 3.556, 3.503, 3.450,
     $  3.398, 3.345, 3.293, 3.242, 3.191, 3.140, 3.090, 3.041,
     $  2.992, 2.944, 2.897, 2.851, 2.805, 2.761, 2.717, 2.674,
     $  2.632, 2.591, 2.551, 2.511, 2.473, 2.435, 2.399, 2.363,
     $  2.328,  2.294,  2.260,  2.228
     $/
c...Sigma0 p => sigma+ n
      data (sfits1(22,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,896.972,102.383,69.149,58.757,54.141,50.925,
     $ 47.704,44.330,41.026,37.964,35.215,32.783,30.644,28.763,
     $ 27.096,25.611,24.280,23.080,21.995,21.006,20.098,19.262,
     $ 18.490,17.774,17.107,16.483,15.897,15.347,14.828,14.337,
     $ 13.872,13.432,13.013,12.613,12.233,11.871,11.525,11.194,
     $ 10.879,10.579,10.291,10.017, 9.755, 9.506, 9.268, 9.041,
     $  8.824, 8.617, 8.420, 8.231, 8.051, 7.879, 7.714, 7.556,
     $  7.404, 7.258, 7.118, 6.983, 6.853, 6.727, 6.606, 6.489,
     $  6.375, 6.265, 6.158, 6.055, 5.954, 5.856, 5.761, 5.669,
     $  5.579, 5.491, 5.406, 5.323, 5.242, 5.163, 5.085, 5.010,
     $  4.937,  4.865,  4.795,  4.727
     $/
c...Sigmap n total
      data (sfits1(23,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,679.244,270.610,190.435,161.526,148.528,140.777,
     $ 134.226,127.748,121.354,115.240,109.525,104.247,99.401,94.958,
     $ 90.870,87.101,83.616,80.389,77.399,74.615,72.017,69.588,
     $ 67.317,65.188,63.187,61.303,59.528,57.857,56.278,54.784,
     $ 53.373,52.039,50.775,49.578,48.444,47.370,46.351,45.385,
     $ 44.469,43.598,42.771,41.985,41.237,40.524,39.845,39.196,
     $ 38.575,37.981,37.412,36.865,36.340,35.835,35.348,34.880,
     $ 34.428,33.992,33.573,33.167,32.777,32.400,32.037,31.687,
     $ 31.350,31.026,30.714,30.414,30.126,29.850,29.586,29.332,
     $ 29.090,28.859,28.638,28.428,28.228,28.038,27.858,27.689,
     $ 27.528, 27.377, 27.235, 27.102
     $/
c...Sigmap n elastic
      data (sfits1(24,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,256.945,97.451,68.045,58.512,55.709,55.126,
     $ 54.809,54.133,53.102,51.846,50.477,49.065,47.651,46.258,
     $ 44.901,43.586,42.319,41.101,39.935,38.821,37.755,36.738,
     $ 35.768,34.843,33.960,33.117,32.313,31.547,30.816,30.118,
     $ 29.452,28.819,28.214,27.637,27.089,26.566,26.069,25.595,
     $ 25.146,24.719,24.312,23.927,23.561,23.213,22.883,22.570,
     $ 22.272,21.989,21.720,21.464,21.221,20.989,20.769,20.559,
     $ 20.359,20.169,19.989,19.818,19.655,19.501,19.356,19.218,
     $ 19.089,18.967,18.853,18.747,18.648,18.556,18.472,18.395,
     $ 18.324,18.261,18.204,18.153,18.110,18.072,18.041,18.016,
     $ 17.998, 17.985, 17.978, 17.977
     $/
c...Sigmap n => lambda p
      data (sfits1(25,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,351.914,100.485,65.404,51.746,44.114,38.957,
     $ 35.106,32.052,29.515,27.347,25.460,23.797,22.320,21.002,
     $ 19.816,18.745,17.774,16.891,16.088,15.354,14.682,14.065,
     $ 13.499,12.979,12.499,12.056,11.647,11.271,10.923,10.601,
     $ 10.303,10.029, 9.776, 9.542, 9.325, 9.125, 8.940, 8.767,
     $  8.607, 8.457, 8.316, 8.182, 8.055, 7.934, 7.816, 7.702,
     $  7.591, 7.482, 7.374, 7.268, 7.162, 7.057, 6.953, 6.849,
     $  6.745, 6.643, 6.540, 6.439, 6.338, 6.238, 6.140, 6.042,
     $  5.946, 5.852, 5.759, 5.667, 5.577, 5.489, 5.402, 5.317,
     $  5.234, 5.152, 5.073, 4.995, 4.919, 4.844, 4.772, 4.701,
     $  4.631,  4.563,  4.497,  4.433
     $/
c...Sigmap n => sigma0 p
      data (sfits1(26,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,70.386,72.674,56.985,51.268,48.705,46.694,
     $ 44.311,41.562,38.738,36.047,33.588,31.385,29.430,27.697,
     $ 26.153,24.770,23.524,22.398,21.376,20.440,19.580,18.785,
     $ 18.050,17.366,16.728,16.130,15.568,15.039,14.539,14.066,
     $ 13.617,13.191,12.786,12.399,12.030,11.679,11.343,11.022,
     $ 10.715,10.423,10.143, 9.876, 9.621, 9.377, 9.145, 8.924,
     $  8.712, 8.510, 8.317, 8.133, 7.957, 7.788, 7.627, 7.472,
     $  7.323, 7.181, 7.043, 6.911, 6.783, 6.660, 6.541, 6.426,
     $  6.315, 6.207, 6.102, 6.000, 5.901, 5.805, 5.711, 5.620,
     $  5.532, 5.446, 5.362, 5.280, 5.200, 5.122, 5.046, 4.972,
     $  4.899,  4.829,  4.760,  4.693
     $/

      end

c***********************************************************************

      subroutine jamsigS2(sig,isig,s)

c...S=-2 low energy Baryon-Baryon cross sections below particle production
c...from Nijmegen model D with had core radius 0.5fm.
      implicit double precision(a-h, o-z)

c...Table parameters.
      parameter (itblsz=100)
      parameter (tblhig=2.72885799)
      parameter (tbllow=2.2332)    ! lambda lambda threthold + 1MeV
      parameter (tblstp=0.002)
      parameter (isgmax=65)

      dimension ems(2,65)
      common /tbls2/sfits2(isgmax,itblsz)
      data ems/
     1 1.32132, 0.938,  1.32132, 0.938, 1.116, 1.197,  1.192, 1.197,
     5 1.32132, 0.938,  1.32132, 0.938, 1.116, 1.116,  1.3149, 0.938,
     9 1.116,   1.192,  1.189, 1.197,   1.192, 1.192,  1.3149,0.938,
     3 1.3149,  0.938,  1.116,1.116,    1.32132,0.938, 1.116,1.192,
     7 1.189,   1.197,  1.192,1.192,    1.3149,0.938,  1.3149,0.938,
     1 1.116,   1.189,  1.189,1.192,    1.116,1.116,   1.116,1.116,
     5 1.3149,  0.938,  1.32132,0.938,  1.192,1.192,   1.189,1.197,
     9 1.116,   1.197,  1.116,1.197,    1.32132,0.938, 1.192,1.197,
     3 1.116,   1.192,  1.116, 1.192,   1.3149,0.938,  1.32132,0.938,
     3 1.189,   1.197,  1.116,1.189,    1.116,1.189,   1.3149,0.938,
     1 1.189,   1.192,  1.198,1.198,    1.192,1.197,   1.192,1.197,
     4 1.32132, 0.938,  1.116,1.197,    1.192,1.192,   1.192,1.192,
     4 1.116,   1.116,  1.3149,0.938,   1.32132,0.938, 1.116,1.192,
     5 1.189,   1.197,  1.189,1.197,    1.189,1.197,   1.116,1.116,
     7 1.3149,  0.938,  1.32132,0.938,  1.116,1.192,   1.192,1.192,
     6 1.189,   1.192,  1.189,1.192,    1.32132,0.938, 1.116,1.189,
     5 1.189,   1.189/

      if(s.le.ems(1,isig)+ems(2,isig)+0.001) then
        sig=0.0d0
        return
      endif

c...1. xi- n total(mb)                   5. xi- p total(mb)
c...2. xi- n elastic                     6. xi- p elastic
c...3. xi- n => lambda sigma-            7. xi- p => ll
c...4. xi- n => sigma0 sigma-            8. xi- p => g0n
c...                                     9. xi- p => ls0
c...                                    10. xi- p => s+s-
c...                                    11. xi- p => s0s0
c...12. g0n total(mb)                   19. g0p total(mb)
c...13. g0n elastic ll                  20. g0p elastic
c...14. g0n ll                          21. g0p lam s+
c...15. g0n g-p                         22. g0p s+s0
c...16. g0n ls0
c...17. g0n s+s-
c...18. g0n s0s0

c...23. lam lam total
c...24. lam lam elastic
c...25. lam lam => xi(0)n
c...26. lam lam => xi(-)p
c...27. lam lam => s(0) s(0)
c...28. lam lam => s(+) s(-)

c...29. lambda s- total(mb)
c...30. lambda s- elastic
c...31. lambda s- g-n
c...32. lambda s- s0s-

c...33. lambda sigma0 total(mb)
c...34. lambda sigma0 elastic
c...35. lambda sigma0 => xi0 n
c...36. lambda sigma0 => xi- p
c...37. lambda sigma0 => s+s-

c...38. lambda sigma+ total(mb)
c...39. lambda sigma+ elastic
c...40. lambda sigma+  => xi0 p
c...41. lambda sigma+  => s(+) s(0)

c...42. s(-) s(-) => s(-) s(-)
c...43. sigma0 sigma- total(mb)
c...44. sigma0 sigma- elastic
c...45. sigma0 sigma- => xi- n
c...46. sigma0 sigma- => lambda sigma-

c...47. sigma0 sigma0 total
c...48. sigma0 sigma0 elastic
c...49. sigma0 sigma0 => lam lam
c...50. sigma0 sigma0 => xi0 n
c...51. sigma0 sigma0 => xi- p
c...52. sigma0 sigma0 => lam s0
c...53. sigma0 sigma0 => s+ s-

c...54. sigma+ sigma- total(mb)
c...55. sigma+ sigma- elastic
c...56. sigma+ sigma- => ll
c...57. sigma+ sigma- => xi0 n
c...58. sigma+ sigma- => xi- p
c...59. sigma+ sigma- => l s0
c...60. sigma+ sigma- => s0s0

c...61. sigma+ sigma0 total(mb)
c...62. sigma+ sigma0 elastic
c...63. sigma+ sigma0 => xi0 p
c...64. sigma+ sigma0 => lam s+
c...65. s(+) s(+) => s(+) s(+)

      xpt=log(s)
      tlow=log(tbllow)
      if(xpt.lt.tlow) xpt=tlow
      index=int((xpt-tlow)/tblstp)+1
      if(index.ge.itblsz) then
c       write(6,*)'warning:index.ge.itblsz in tbls2 index=',index
c       write(6,*)'srt isig',s,isig
c       index=itblsz-1
        sig=sfits2(isig,itblsz)
        return
      end if

c...Find slopes and crossections.
      x1=(index-1)*tblstp+tlow
      y1in=sfits2(isig,index)
      y2in=sfits2(isig,index+1)
      slin=(y2in-y1in)/tblstp
      sig=slin*(xpt-x1)+y1in

      end

c***********************************************************************

      block data jams2da

c...S=-2 baryon - baryon cross sections extracted from Njimegen model D.
      implicit double precision(a-h, o-z)
      parameter (itblsz=100)
      parameter (isgmax=65)

      common /tbls2/sfits2(isgmax,itblsz)

c...1. xi- n total(mb)
      data (sfits2(1,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,58.823,42.293,
     $ 26.024,18.701,15.669,14.905,15.289,16.290,17.616,19.107,
     $ 20.740,22.788,31.425,34.616,36.447,37.577,38.136,38.239,
     $ 38.020,37.618,37.121,36.592,36.076,35.596,35.171,34.806,
     $ 34.516,34.311,34.217,33.919,33.582,33.199,32.810,32.444,
     $ 32.122,31.846,31.619,31.435,31.287,31.171,31.081,31.013,
     $ 30.963,30.928,30.905,30.894,30.892,30.899,30.913,30.936,
     $ 30.966,31.004,31.048,31.100,31.160,31.227,31.301,31.383,
     $ 31.472,31.569,31.672,31.782,31.898,32.020,32.146,32.275,
     $ 32.407,32.541,32.674,32.807,32.938,33.065,33.187,33.304,
     $ 33.414,33.516,33.611,33.697,33.774,33.841,33.900,33.949,
     $ 33.990,34.023,34.048,34.065,34.076,34.081,34.080,34.075,
     $ 34.065, 34.052, 34.036, 34.018
     $/
c...2. xi- n elastic
      data (sfits2(2,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,58.823,42.293,
     $ 26.024,18.701,15.669,14.905,15.289,16.290,17.616,19.107,
     $ 20.740,22.788,25.088,25.902,26.509,26.642,26.407,25.941,
     $ 25.367,24.777,24.222,23.728,23.305,22.956,22.679,22.468,
     $ 22.326,22.254,22.298,22.118,21.882,21.628,21.393,21.195,
     $ 21.044,20.934,20.863,20.822,20.806,20.807,20.823,20.849,
     $ 20.882,20.921,20.963,21.008,21.055,21.104,21.154,21.205,
     $ 21.257,21.309,21.363,21.418,21.474,21.531,21.589,21.648,
     $ 21.709,21.771,21.833,21.897,21.961,22.026,22.091,22.155,
     $ 22.219,22.283,22.344,22.404,22.461,22.516,22.568,22.617,
     $ 22.662,22.704,22.742,22.776,22.806,22.832,22.855,22.874,
     $ 22.890,22.903,22.913,22.920,22.926,22.929,22.931,22.931,
     $ 22.930, 22.929, 22.926, 22.923
     $/
c...3. xi- n => lambda sigma-
      data (sfits2(3,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 6.337, 8.714, 9.937,10.935,11.729,12.299,
     $ 12.654,12.841,12.899,12.864,12.771,12.641,12.492,12.338,
     $ 12.190,12.057,11.874,11.538,11.167,10.763,10.369,10.012,
     $  9.707, 9.452, 9.245, 9.076, 8.941, 8.833, 8.746, 8.677,
     $  8.621, 8.577, 8.543, 8.516, 8.496, 8.481, 8.472, 8.467,
     $  8.466, 8.469, 8.476, 8.486, 8.499, 8.516, 8.536, 8.560,
     $  8.586, 8.615, 8.646, 8.680, 8.716, 8.754, 8.793, 8.834,
     $  8.874, 8.915, 8.955, 8.994, 9.031, 9.067, 9.100, 9.130,
     $  9.156, 9.180, 9.199, 9.215, 9.227, 9.235, 9.240, 9.240,
     $  9.238, 9.232, 9.224, 9.213, 9.200, 9.185, 9.168, 9.150,
     $  9.131,  9.112,  9.091,  9.071
     $/
c...4. xi- n => sigma0 sigma-
      data (sfits2(4,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.045, 0.263, 0.533, 0.808, 1.048, 1.236,
     $  1.370, 1.460, 1.511, 1.536, 1.540, 1.531, 1.513, 1.488,
     $  1.460, 1.430, 1.399, 1.369, 1.340, 1.313, 1.288, 1.264,
     $  1.243, 1.225, 1.209, 1.197, 1.187, 1.180, 1.176, 1.175,
     $  1.178, 1.184, 1.193, 1.205, 1.221, 1.240, 1.261, 1.286,
     $  1.314, 1.343, 1.375, 1.409, 1.445, 1.482, 1.519, 1.557,
     $  1.595, 1.633, 1.670, 1.706, 1.741, 1.774, 1.806, 1.835,
     $  1.863, 1.888, 1.911, 1.932, 1.950, 1.967, 1.981, 1.993,
     $  2.003,  2.012,  2.018,  2.023
     $/
c...5. xi- p total(mb)
      data (sfits2(5,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,405.080,91.818,
     $ 50.229,34.833,28.198,25.486,24.592,24.625,25.134,25.893,
     $ 26.969,31.564,32.926,33.715,34.160,34.290,34.183,33.914,
     $ 33.549,33.144,32.733,32.341,31.980,31.658,31.379,31.143,
     $ 30.958,30.885,31.145,31.018,30.833,30.637,30.451,30.286,
     $ 30.149,30.039,29.953,29.888,29.841,29.808,29.787,29.774,
     $ 29.768,29.767,29.770,29.774,29.781,29.789,29.798,29.808,
     $ 29.819,29.831,29.844,29.859,29.875,29.894,29.915,29.939,
     $ 29.966,29.995,30.028,30.064,30.103,30.145,30.189,30.235,
     $ 30.283,30.333,30.383,30.433,30.483,30.532,30.579,30.624,
     $ 30.668,30.708,30.745,30.779,30.809,30.836,30.860,30.880,
     $ 30.897,30.911,30.922,30.930,30.937,30.941,30.943,30.944,
     $ 30.944, 30.943, 30.942, 30.940
     $/
c...6. xi- p elastic
      data (sfits2(6,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,177.437,63.493,
     $ 37.587,27.117,22.653,21.076,20.877,21.380,22.226,23.241,
     $ 24.472,25.467,25.828,26.034,25.966,25.671,25.248,24.772,
     $ 24.295,23.850,23.451,23.106,22.812,22.568,22.369,22.207,
     $ 22.079,21.982,22.028,21.910,21.776,21.659,21.574,21.521,
     $ 21.499,21.500,21.520,21.553,21.595,21.642,21.691,21.742,
     $ 21.791,21.838,21.883,21.924,21.963,21.998,22.031,22.060,
     $ 22.087,22.112,22.135,22.156,22.176,22.195,22.214,22.232,
     $ 22.250,22.268,22.285,22.303,22.320,22.338,22.355,22.371,
     $ 22.386,22.401,22.415,22.427,22.438,22.447,22.454,22.459,
     $ 22.462,22.463,22.462,22.459,22.454,22.447,22.438,22.428,
     $ 22.417,22.404,22.391,22.376,22.361,22.346,22.330,22.314,
     $ 22.298, 22.282, 22.266, 22.251
     $/
c...7. xi- p => ll
      data (sfits2(7,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,101.227,12.416,
     $  5.736, 3.724, 2.881, 2.468, 2.231, 2.081, 1.977, 1.897,
     $  1.831, 1.770, 1.717, 1.668, 1.623, 1.580, 1.540, 1.503,
     $  1.467, 1.434, 1.402, 1.372, 1.343, 1.315, 1.287, 1.260,
     $  1.232, 1.212, 1.223, 1.211, 1.194, 1.177, 1.159, 1.140,
     $  1.121, 1.100, 1.079, 1.057, 1.035, 1.013, 0.990, 0.969,
     $  0.947, 0.926, 0.906, 0.887, 0.869, 0.851, 0.835, 0.819,
     $  0.805, 0.792, 0.779, 0.767, 0.757, 0.747, 0.738, 0.729,
     $  0.721, 0.714, 0.708, 0.701, 0.696, 0.691, 0.686, 0.681,
     $  0.677, 0.673, 0.669, 0.666, 0.663, 0.660, 0.657, 0.654,
     $  0.652, 0.650, 0.647, 0.645, 0.644, 0.642, 0.640, 0.639,
     $  0.637, 0.636, 0.634, 0.633, 0.632, 0.631, 0.630, 0.629,
     $  0.628,  0.627,  0.626,  0.626
     $/
c...8. xi- p => g0n
      data (sfits2(8,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,126.417,15.909,
     $  6.906, 3.992, 2.664, 1.943, 1.484, 1.165, 0.930, 0.754,
     $  0.665, 0.773, 0.712, 0.729, 0.786, 0.860, 0.937, 1.010,
     $  1.075, 1.132, 1.181, 1.225, 1.265, 1.303, 1.342, 1.386,
     $  1.439, 1.532, 1.466, 1.480, 1.503, 1.518, 1.523, 1.520,
     $  1.511, 1.500, 1.487, 1.475, 1.463, 1.452, 1.441, 1.432,
     $  1.423, 1.414, 1.406, 1.399, 1.391, 1.383, 1.375, 1.367,
     $  1.359, 1.351, 1.343, 1.335, 1.327, 1.320, 1.312, 1.305,
     $  1.299, 1.292, 1.287, 1.282, 1.279, 1.276, 1.274, 1.274,
     $  1.275, 1.277, 1.280, 1.285, 1.291, 1.298, 1.307, 1.317,
     $  1.329, 1.341, 1.355, 1.370, 1.385, 1.401, 1.419, 1.436,
     $  1.454, 1.473, 1.492, 1.511, 1.530, 1.549, 1.568, 1.587,
     $  1.606,  1.625,  1.644,  1.662
     $/
c...9. xi- p => ls0
      data (sfits2(9,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 3.554, 4.669, 5.284, 5.786, 6.178, 6.457, 6.629,
     $  6.711, 6.728, 6.698, 6.638, 6.559, 6.471, 6.380, 6.291,
     $  6.207, 6.133, 5.986, 5.810, 5.609, 5.404, 5.214, 5.048,
     $  4.909, 4.794, 4.701, 4.626, 4.565, 4.517, 4.477, 4.446,
     $  4.421, 4.401, 4.385, 4.373, 4.363, 4.356, 4.352, 4.350,
     $  4.349, 4.351, 4.354, 4.358, 4.365, 4.373, 4.382, 4.393,
     $  4.406, 4.419, 4.434, 4.450, 4.467, 4.485, 4.503, 4.521,
     $  4.540, 4.558, 4.576, 4.593, 4.610, 4.625, 4.639, 4.651,
     $  4.662, 4.670, 4.677, 4.682, 4.686, 4.687, 4.687, 4.684,
     $  4.681, 4.676, 4.669, 4.662, 4.653, 4.644, 4.634, 4.624,
     $  4.613,  4.602,  4.591,  4.580
     $/
c...10. xi- p => s+s-
      data (sfits2(10,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.322, 0.476, 0.619, 0.746, 0.848, 0.925,
     $  0.977, 1.011, 1.030, 1.040, 1.043, 1.042, 1.038, 1.034,
     $  1.029, 1.025, 1.020, 1.017, 1.015, 1.014, 1.014, 1.014,
     $  1.015, 1.018, 1.021, 1.024, 1.029, 1.034, 1.041, 1.048,
     $  1.056, 1.065, 1.075, 1.086, 1.098, 1.112, 1.126, 1.141,
     $  1.158, 1.175, 1.193, 1.212, 1.231, 1.251, 1.270, 1.290,
     $  1.310, 1.330, 1.349, 1.367, 1.385, 1.403, 1.419, 1.434,
     $  1.449, 1.463, 1.475, 1.487, 1.497, 1.507, 1.516, 1.524,
     $  1.531,  1.538,  1.543,  1.549
     $/
c...11. xi- p => s0s0
      data (sfits2(11,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.026, 0.119, 0.130, 0.132, 0.132, 0.132, 0.132,
     $  0.133, 0.134, 0.135, 0.137, 0.140, 0.144, 0.148, 0.153,
     $  0.158, 0.163, 0.169, 0.175, 0.181, 0.186, 0.192, 0.198,
     $  0.203, 0.208, 0.213, 0.217, 0.221, 0.225, 0.229, 0.232,
     $  0.234, 0.237, 0.239, 0.241, 0.243, 0.244, 0.246, 0.247,
     $  0.248, 0.249, 0.250, 0.250, 0.251, 0.252, 0.252, 0.253,
     $  0.253, 0.254, 0.254, 0.255, 0.256, 0.256, 0.257, 0.258,
     $  0.259, 0.260, 0.261, 0.262, 0.263, 0.264, 0.265, 0.267,
     $  0.268,  0.270,  0.271,  0.273
     $/
c...12. g0n total(mb)
      data (sfits2(12,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000,252.422,97.383,48.561,
     $ 33.659,27.359,24.768,23.986,24.037,24.479,25.092,25.784,
     $ 26.648,30.726,32.077,32.825,33.232,33.348,33.253,33.017,
     $ 32.699,32.349,31.997,31.662,31.356,31.086,30.854,30.660,
     $ 30.510,30.457,30.739,30.651,30.502,30.337,30.178,30.036,
     $ 29.918,29.823,29.750,29.695,29.656,29.630,29.613,29.605,
     $ 29.603,29.605,29.610,29.617,29.626,29.636,29.648,29.660,
     $ 29.673,29.688,29.703,29.721,29.740,29.762,29.786,29.813,
     $ 29.842,29.875,29.911,29.950,29.992,30.036,30.083,30.132,
     $ 30.182,30.234,30.286,30.338,30.389,30.439,30.488,30.535,
     $ 30.578,30.619,30.657,30.691,30.722,30.749,30.773,30.793,
     $ 30.810,30.824,30.835,30.844,30.850,30.854,30.857,30.858,
     $ 30.858, 30.858, 30.857, 30.855
     $/
c...13. g0n elastic ll
      data (sfits2(13,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000,174.611,64.199,34.287,
     $ 25.250,21.485,20.201,20.175,20.722,21.519,22.402,23.306,
     $ 24.300,25.179,25.583,25.749,25.643,25.331,24.912,24.455,
     $ 24.008,23.596,23.232,22.921,22.659,22.444,22.270,22.129,
     $ 22.018,21.927,21.988,21.895,21.781,21.681,21.607,21.563,
     $ 21.546,21.550,21.571,21.603,21.643,21.687,21.734,21.781,
     $ 21.826,21.870,21.911,21.949,21.985,22.017,22.046,22.073,
     $ 22.098,22.120,22.141,22.161,22.180,22.197,22.215,22.232,
     $ 22.249,22.266,22.283,22.300,22.317,22.333,22.350,22.365,
     $ 22.380,22.394,22.407,22.418,22.428,22.436,22.442,22.446,
     $ 22.448,22.447,22.445,22.441,22.434,22.426,22.416,22.405,
     $ 22.393,22.379,22.364,22.348,22.332,22.316,22.299,22.282,
     $ 22.265, 22.249, 22.232, 22.217
     $/
c...14. g0n ll
      data (sfits2(14,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000,77.811,20.144, 6.485,
     $  3.926, 2.941, 2.470, 2.215, 2.058, 1.950, 1.870, 1.804,
     $  1.747, 1.698, 1.651, 1.607, 1.566, 1.527, 1.491, 1.457,
     $  1.425, 1.395, 1.367, 1.340, 1.314, 1.288, 1.262, 1.237,
     $  1.211, 1.193, 1.204, 1.192, 1.176, 1.160, 1.143, 1.125,
     $  1.106, 1.086, 1.065, 1.044, 1.022, 1.000, 0.979, 0.957,
     $  0.937, 0.916, 0.897, 0.878, 0.860, 0.843, 0.827, 0.812,
     $  0.798, 0.785, 0.773, 0.761, 0.751, 0.741, 0.732, 0.724,
     $  0.717, 0.710, 0.703, 0.697, 0.692, 0.687, 0.682, 0.678,
     $  0.674, 0.670, 0.666, 0.663, 0.660, 0.657, 0.655, 0.652,
     $  0.650, 0.648, 0.646, 0.644, 0.642, 0.640, 0.639, 0.637,
     $  0.636, 0.634, 0.633, 0.632, 0.631, 0.630, 0.629, 0.628,
     $  0.627,  0.627,  0.626,  0.625
     $/
c...15. g0n g-p
      data (sfits2(15,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,13.040, 7.789,
     $  4.483, 2.933, 2.098, 1.597, 1.258, 1.010, 0.820, 0.674,
     $  0.601, 0.705, 0.654, 0.674, 0.730, 0.803, 0.879, 0.951,
     $  1.015, 1.072, 1.122, 1.167, 1.208, 1.247, 1.288, 1.332,
     $  1.385, 1.477, 1.416, 1.432, 1.456, 1.473, 1.479, 1.478,
     $  1.471, 1.462, 1.451, 1.440, 1.430, 1.421, 1.412, 1.404,
     $  1.397, 1.390, 1.383, 1.377, 1.370, 1.363, 1.357, 1.350,
     $  1.343, 1.336, 1.329, 1.322, 1.315, 1.309, 1.302, 1.296,
     $  1.290, 1.285, 1.280, 1.276, 1.273, 1.271, 1.270, 1.270,
     $  1.272, 1.274, 1.278, 1.283, 1.290, 1.298, 1.307, 1.317,
     $  1.329, 1.342, 1.356, 1.371, 1.387, 1.403, 1.420, 1.438,
     $  1.457, 1.475, 1.494, 1.513, 1.533, 1.552, 1.571, 1.590,
     $  1.609,  1.628,  1.647,  1.666
     $/
c...16. g0n ls0
      data (sfits2(16,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 3.144, 4.189, 4.794, 5.293, 5.686, 5.971, 6.154,
     $  6.250, 6.285, 6.275, 6.234, 6.176, 6.107, 6.035, 5.963,
     $  5.896, 5.836, 5.710, 5.553, 5.369, 5.181, 5.007, 4.855,
     $  4.727, 4.621, 4.537, 4.468, 4.414, 4.370, 4.336, 4.308,
     $  4.287, 4.270, 4.258, 4.248, 4.242, 4.238, 4.236, 4.236,
     $  4.238, 4.242, 4.248, 4.255, 4.263, 4.274, 4.285, 4.299,
     $  4.313, 4.329, 4.346, 4.364, 4.383, 4.402, 4.423, 4.443,
     $  4.463, 4.483, 4.503, 4.522, 4.539, 4.556, 4.571, 4.585,
     $  4.597, 4.607, 4.615, 4.621, 4.625, 4.627, 4.628, 4.627,
     $  4.624, 4.620, 4.614, 4.608, 4.600, 4.592, 4.583, 4.573,
     $  4.563,  4.553,  4.542,  4.532
     $/
c...17. g0n s+s-
      data (sfits2(17,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.307, 0.455, 0.592, 0.715, 0.814, 0.888,
     $  0.940, 0.974, 0.994, 1.005, 1.009, 1.010, 1.007, 1.004,
     $  1.001, 0.998, 0.995, 0.993, 0.992, 0.992, 0.992, 0.994,
     $  0.996, 0.999, 1.002, 1.007, 1.012, 1.018, 1.025, 1.033,
     $  1.042, 1.051, 1.062, 1.074, 1.087, 1.101, 1.116, 1.132,
     $  1.148, 1.166, 1.185, 1.204, 1.223, 1.244, 1.264, 1.284,
     $  1.304, 1.324, 1.343, 1.362, 1.380, 1.398, 1.414, 1.430,
     $  1.444, 1.458, 1.471, 1.482, 1.493, 1.503, 1.512, 1.520,
     $  1.527,  1.533,  1.539,  1.544
     $/
c...18. g0n s0s0
      data (sfits2(18,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.025, 0.114, 0.125, 0.127, 0.128, 0.128, 0.128,
     $  0.129, 0.130, 0.132, 0.134, 0.137, 0.141, 0.145, 0.150,
     $  0.155, 0.161, 0.166, 0.172, 0.178, 0.184, 0.190, 0.195,
     $  0.200, 0.206, 0.210, 0.215, 0.219, 0.223, 0.226, 0.229,
     $  0.232, 0.234, 0.237, 0.238, 0.240, 0.242, 0.243, 0.244,
     $  0.245, 0.246, 0.247, 0.248, 0.249, 0.249, 0.250, 0.251,
     $  0.251, 0.252, 0.252, 0.253, 0.254, 0.255, 0.255, 0.256,
     $  0.257, 0.258, 0.259, 0.260, 0.261, 0.262, 0.264, 0.265,
     $  0.267,  0.268,  0.270,  0.271
     $/
c...19. g0p total(mb)
      data (sfits2(19,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000,49.867,29.026,20.091,
     $ 16.173,14.850,14.941,15.803,17.044,18.484,20.067,21.922,
     $ 28.981,33.776,35.834,37.188,37.939,38.179,38.070,37.733,
     $ 37.269,36.755,36.239,35.750,35.307,34.923,34.608,34.368,
     $ 34.245,34.019,33.704,33.333,32.943,32.567,32.228,31.935,
     $ 31.691,31.490,31.330,31.202,31.102,31.026,30.968,30.927,
     $ 30.899,30.884,30.878,30.881,30.892,30.911,30.938,30.972,
     $ 31.013,31.062,31.118,31.181,31.251,31.329,31.414,31.507,
     $ 31.607,31.714,31.826,31.945,32.068,32.196,32.327,32.460,
     $ 32.594,32.728,32.860,32.990,33.116,33.236,33.351,33.459,
     $ 33.559,33.651,33.734,33.809,33.874,33.929,33.977,34.015,
     $ 34.045,34.068,34.084,34.093,34.097,34.095,34.089,34.079,
     $ 34.065, 34.049, 34.030, 34.010
     $/
c...20. g0p elastic
      data (sfits2(20,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000,49.867,29.026,20.091,
     $ 16.173,14.850,14.941,15.803,17.044,18.484,20.067,21.922,
     $ 24.581,25.470,26.205,26.506,26.398,26.006,25.475,24.899,
     $ 24.341,23.834,23.392,23.024,22.725,22.495,22.331,22.232,
     $ 22.232,22.152,21.934,21.683,21.440,21.230,21.062,20.938,
     $ 20.855,20.804,20.779,20.775,20.786,20.809,20.840,20.877,
     $ 20.918,20.963,21.010,21.059,21.109,21.160,21.211,21.264,
     $ 21.318,21.372,21.428,21.485,21.543,21.602,21.662,21.723,
     $ 21.786,21.849,21.913,21.978,22.043,22.109,22.173,22.237,
     $ 22.300,22.361,22.421,22.478,22.532,22.583,22.631,22.675,
     $ 22.716,22.753,22.786,22.815,22.841,22.863,22.881,22.897,
     $ 22.909,22.919,22.926,22.931,22.934,22.936,22.937,22.936,
     $ 22.934, 22.932, 22.929, 22.926
     $/
c...21. g0p lam s+
      data (sfits2(21,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  4.400, 8.306, 9.629,10.682,11.541,12.173,12.595,12.834,
     $ 12.928,12.922,12.846,12.726,12.582,12.428,12.277,12.136,
     $ 12.002,11.686,11.336,10.938,10.537,10.166, 9.841, 9.567,
     $  9.342, 9.158, 9.010, 8.891, 8.795, 8.718, 8.657, 8.608,
     $  8.569, 8.538, 8.515, 8.497, 8.485, 8.478, 8.475, 8.476,
     $  8.480, 8.489, 8.500, 8.515, 8.533, 8.555, 8.579, 8.607,
     $  8.637, 8.669, 8.704, 8.741, 8.780, 8.819, 8.860, 8.900,
     $  8.941, 8.981, 9.019, 9.056, 9.091, 9.123, 9.152, 9.177,
     $  9.199, 9.218, 9.232, 9.243, 9.250, 9.253, 9.253, 9.249,
     $  9.243, 9.234, 9.222, 9.209, 9.193, 9.176, 9.158, 9.139,
     $  9.119,  9.099,  9.079,  9.058
     $/
c...22. g0p s+s0
      data (sfits2(22,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.010, 0.181, 0.434, 0.711, 0.965, 1.171, 1.325, 1.430,
     $  1.494, 1.529, 1.540, 1.536, 1.521, 1.499, 1.472, 1.443,
     $  1.412, 1.382, 1.353, 1.325, 1.298, 1.274, 1.252, 1.232,
     $  1.215, 1.201, 1.189, 1.181, 1.175, 1.173, 1.174, 1.178,
     $  1.185, 1.195, 1.209, 1.226, 1.245, 1.268, 1.294, 1.322,
     $  1.353, 1.386, 1.420, 1.456, 1.493, 1.531, 1.569, 1.607,
     $  1.644, 1.681, 1.716, 1.750, 1.783, 1.814, 1.842, 1.869,
     $  1.894, 1.916, 1.936, 1.954, 1.969, 1.983, 1.994, 2.004,
     $  2.012,  2.018,  2.023,  2.026
     $/
c...23. lam lam total
      data (sfits2(23,i),i=1,itblsz)/
     $ 480.293,118.594,63.477,42.569,32.653,27.600,19.857,14.985,
     $ 12.823,11.417,10.411, 9.654, 9.052, 8.562, 8.159, 7.826,
     $  7.553, 7.336, 7.170, 7.050, 6.973, 6.937, 6.936, 6.967,
     $  7.028, 7.113, 7.218, 7.341, 7.477, 7.621, 7.768, 7.912,
     $  8.044, 8.149, 8.888, 9.271, 9.585, 9.868,10.130,10.376,
     $ 10.610,10.835,11.052,11.263,11.468,11.668,11.862,12.051,
     $ 12.234,12.411,12.583,12.748,12.907,13.059,13.206,13.345,
     $ 13.479,13.606,13.728,13.843,13.953,14.058,14.158,14.254,
     $ 14.345,14.432,14.516,14.596,14.673,14.747,14.818,14.887,
     $ 14.954,15.018,15.081,15.142,15.201,15.259,15.315,15.370,
     $ 15.424,15.477,15.528,15.579,15.628,15.676,15.724,15.771,
     $ 15.816,15.862,15.906,15.949,15.992,16.034,16.075,16.116,
     $ 16.156, 16.195, 16.233, 16.271
     $/
c...24. lam lam elastic
      data (sfits2(24,i),i=1,itblsz)/
     $ 480.293,118.594,63.477,42.569,32.653,23.864,14.023,11.277,
     $  9.941, 8.909, 8.072, 7.378, 6.793, 6.297, 5.883, 5.541,
     $  5.263, 5.046, 4.885, 4.775, 4.711, 4.690, 4.707, 4.758,
     $  4.838, 4.944, 5.071, 5.216, 5.374, 5.542, 5.715, 5.886,
     $  6.047, 6.135, 6.378, 6.667, 6.938, 7.199, 7.453, 7.701,
     $  7.944, 8.182, 8.415, 8.643, 8.864, 9.079, 9.288, 9.489,
     $  9.683, 9.869,10.048,10.218,10.381,10.536,10.683,10.823,
     $ 10.957,11.084,11.205,11.320,11.429,11.534,11.634,11.731,
     $ 11.823,11.912,11.998,12.080,12.160,12.238,12.313,12.386,
     $ 12.458,12.527,12.594,12.660,12.724,12.787,12.848,12.908,
     $ 12.966,13.023,13.079,13.133,13.186,13.237,13.288,13.337,
     $ 13.385,13.431,13.477,13.521,13.564,13.606,13.646,13.686,
     $ 13.724, 13.761, 13.797, 13.832
     $/
c...25. lam lam => xi(0)n
      data (sfits2(25,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 3.735, 3.828, 1.913,
     $  1.480, 1.300, 1.221, 1.190, 1.180, 1.178, 1.180, 1.180,
     $  1.179, 1.176, 1.171, 1.163, 1.154, 1.144, 1.134, 1.123,
     $  1.111, 1.100, 1.089, 1.077, 1.065, 1.052, 1.039, 1.025,
     $  1.010, 1.001, 1.016, 1.011, 1.003, 0.994, 0.983, 0.972,
     $  0.959, 0.945, 0.931, 0.915, 0.899, 0.883, 0.866, 0.850,
     $  0.834, 0.818, 0.802, 0.787, 0.773, 0.759, 0.747, 0.735,
     $  0.723, 0.713, 0.703, 0.694, 0.686, 0.678, 0.671, 0.665,
     $  0.659, 0.653, 0.648, 0.643, 0.639, 0.635, 0.632, 0.628,
     $  0.625, 0.622, 0.620, 0.617, 0.615, 0.613, 0.611, 0.609,
     $  0.608, 0.606, 0.605, 0.603, 0.602, 0.601, 0.600, 0.599,
     $  0.598, 0.597, 0.596, 0.596, 0.595, 0.594, 0.594, 0.593,
     $  0.593,  0.592,  0.592,  0.592
     $/
c...26. lam lam => xi(-)p
      data (sfits2(26,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.006, 1.795,
     $  1.402, 1.208, 1.119, 1.086, 1.080, 1.086, 1.096, 1.105,
     $  1.112, 1.114, 1.115, 1.113, 1.109, 1.103, 1.095, 1.087,
     $  1.078, 1.069, 1.059, 1.048, 1.037, 1.026, 1.014, 1.001,
     $  0.987, 0.979, 0.994, 0.990, 0.983, 0.975, 0.965, 0.954,
     $  0.943, 0.930, 0.916, 0.901, 0.886, 0.870, 0.855, 0.839,
     $  0.823, 0.807, 0.792, 0.778, 0.764, 0.751, 0.738, 0.727,
     $  0.716, 0.705, 0.696, 0.687, 0.679, 0.671, 0.664, 0.658,
     $  0.652, 0.647, 0.642, 0.637, 0.633, 0.629, 0.626, 0.623,
     $  0.620, 0.617, 0.614, 0.612, 0.610, 0.608, 0.606, 0.604,
     $  0.602, 0.601, 0.599, 0.598, 0.597, 0.596, 0.595, 0.594,
     $  0.593, 0.592, 0.591, 0.591, 0.590, 0.590, 0.589, 0.588,
     $  0.588,  0.588,  0.587,  0.587
     $/
c...27. lam lam => s(0) s(0)
      data (sfits2(27,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.333, 0.402, 0.441, 0.468, 0.487, 0.501,
     $  0.512, 0.522, 0.531, 0.540, 0.549, 0.560, 0.572, 0.585,
     $  0.598, 0.613, 0.629, 0.645, 0.661, 0.677, 0.693, 0.708,
     $  0.723, 0.738, 0.751, 0.763, 0.775, 0.785, 0.794, 0.802,
     $  0.810, 0.816, 0.821, 0.826, 0.829, 0.832, 0.834, 0.836,
     $  0.837, 0.838, 0.838, 0.838, 0.838, 0.838, 0.837, 0.836,
     $  0.835, 0.835, 0.834, 0.833, 0.832, 0.832, 0.831, 0.831,
     $  0.831, 0.831, 0.831, 0.832, 0.832, 0.833, 0.835, 0.836,
     $  0.838,  0.840,  0.842,  0.844
     $/
c...28. lam lam => s(+) s(-)
      data (sfits2(28,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.035, 0.168, 0.200, 0.220, 0.233, 0.241, 0.247,
     $  0.252, 0.256, 0.260, 0.264, 0.269, 0.275, 0.281, 0.288,
     $  0.296, 0.303, 0.311, 0.320, 0.328, 0.336, 0.344, 0.352,
     $  0.360, 0.367, 0.373, 0.379, 0.385, 0.390, 0.394, 0.398,
     $  0.401, 0.404, 0.407, 0.409, 0.411, 0.412, 0.413, 0.414,
     $  0.414, 0.414, 0.414, 0.414, 0.414, 0.414, 0.413, 0.413,
     $  0.413, 0.412, 0.412, 0.411, 0.411, 0.411, 0.410, 0.410,
     $  0.410, 0.410, 0.410, 0.410, 0.411, 0.411, 0.412, 0.412,
     $  0.413,  0.414,  0.415,  0.417
     $/
c...29. ls- total(mb)
      data (sfits2(29,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.812, 0.597,76.443,64.513,58.536,55.061,
     $ 52.696,50.809,49.178,47.712,46.371,45.151,44.053,43.078,
     $ 42.242,41.569,42.255,42.368,41.807,41.126,40.428,39.764,
     $ 39.165,38.639,38.191,37.810,37.491,37.223,36.999,36.811,
     $ 36.652,36.518,36.404,36.309,36.227,36.159,36.102,36.057,
     $ 36.021,35.995,35.979,35.972,35.975,35.987,36.009,36.041,
     $ 36.081,36.132,36.191,36.259,36.336,36.420,36.511,36.607,
     $ 36.709,36.814,36.921,37.029,37.137,37.244,37.347,37.446,
     $ 37.541,37.629,37.711,37.785,37.851,37.910,37.960,38.002,
     $ 38.037,38.064,38.084,38.097,38.105,38.106,38.103,38.096,
     $ 38.084, 38.070, 38.053, 38.034
     $/
c...30. ls- elastic
      data (sfits2(30,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,59.332,32.114,22.668,19.372,18.492,18.632,
     $ 19.127,19.663,20.130,20.493,20.763,20.955,21.095,21.207,
     $ 21.314,21.454,22.163,22.303,22.306,22.228,22.114,21.998,
     $ 21.901,21.831,21.792,21.780,21.791,21.821,21.866,21.921,
     $ 21.985,22.054,22.126,22.200,22.276,22.352,22.427,22.502,
     $ 22.576,22.649,22.722,22.795,22.867,22.939,23.011,23.084,
     $ 23.157,23.231,23.306,23.382,23.459,23.537,23.615,23.694,
     $ 23.774,23.854,23.934,24.014,24.093,24.172,24.249,24.325,
     $ 24.399,24.471,24.541,24.608,24.673,24.736,24.797,24.854,
     $ 24.910,24.963,25.014,25.063,25.110,25.155,25.198,25.239,
     $ 25.279, 25.318, 25.356, 25.393
     $/
c...31. ls- g-n
      data (sfits2(31,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.310,74.546,53.775,45.142,40.044,36.429,
     $ 33.569,31.146,29.047,27.219,25.609,24.197,22.959,21.872,
     $ 20.928,20.115,19.306,18.313,17.343,16.387,15.503,14.720,
     $ 14.052,13.485,13.013,12.616,12.282,12.000,11.758,11.551,
     $ 11.371,11.214,11.076,10.955,10.848,10.753,10.668,10.594,
     $ 10.528,10.470,10.420,10.377,10.340,10.310,10.286,10.267,
     $ 10.254,10.246,10.242,10.243,10.247,10.255,10.265,10.278,
     $ 10.292,10.307,10.323,10.338,10.352,10.364,10.375,10.383,
     $ 10.388,10.389,10.387,10.382,10.373,10.360,10.344,10.324,
     $ 10.302,10.276,10.248,10.218,10.185,10.152,10.117,10.081,
     $ 10.044, 10.007,  9.971,  9.934
     $/
c...32. ls- s0s-
      data (sfits2(32,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.787, 1.751, 2.159, 2.512, 2.811, 3.046,
     $  3.212, 3.323, 3.386, 3.414, 3.418, 3.403, 3.375, 3.339,
     $  3.297, 3.250, 3.202, 3.153, 3.103, 3.054, 3.007, 2.961,
     $  2.917, 2.875, 2.837, 2.801, 2.768, 2.738, 2.712, 2.689,
     $  2.670, 2.655, 2.643, 2.635, 2.630, 2.628, 2.630, 2.635,
     $  2.642, 2.652, 2.664, 2.677, 2.692, 2.707, 2.723, 2.739,
     $  2.755, 2.769, 2.783, 2.795, 2.805, 2.813, 2.820, 2.824,
     $  2.825, 2.825, 2.822, 2.817, 2.810, 2.800, 2.789, 2.776,
     $  2.761,  2.744,  2.726,  2.707
     $/
c...33. lambda sigma0 total(mb)
      data (sfits2(33,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000,205.999,102.021,74.139,62.776,57.321,54.109,51.887,
     $ 50.123,48.577,47.171,45.878,44.688,43.602,42.623,41.745,
     $ 40.958,40.144,40.807,40.681,40.262,39.732,39.181,38.658,
     $ 38.189,37.780,37.432,37.138,36.893,36.689,36.519,36.377,
     $ 36.258,36.160,36.077,36.009,35.953,35.908,35.873,35.848,
     $ 35.831,35.824,35.825,35.835,35.854,35.882,35.919,35.965,
     $ 36.019,36.083,36.155,36.235,36.322,36.416,36.515,36.619,
     $ 36.727,36.837,36.948,37.059,37.169,37.276,37.379,37.476,
     $ 37.568,37.654,37.731,37.802,37.864,37.917,37.963,38.001,
     $ 38.030,38.053,38.069,38.078,38.082,38.081,38.075,38.065,
     $ 38.051, 38.035, 38.017, 37.997
     $/
c...34. lambda sigma0 elastic
      data (sfits2(34,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000,53.518,30.068,21.562,18.573,17.947,18.214,18.781,
     $ 19.368,19.871,20.262,20.546,20.747,20.886,20.982,21.050,
     $ 21.087,20.989,21.216,21.496,21.596,21.598,21.559,21.514,
     $ 21.481,21.467,21.475,21.502,21.545,21.600,21.665,21.737,
     $ 21.814,21.893,21.974,22.056,22.138,22.219,22.299,22.378,
     $ 22.457,22.534,22.611,22.687,22.764,22.840,22.916,22.993,
     $ 23.070,23.148,23.227,23.307,23.387,23.468,23.550,23.632,
     $ 23.715,23.797,23.879,23.961,24.041,24.121,24.199,24.275,
     $ 24.349,24.421,24.491,24.558,24.622,24.685,24.744,24.801,
     $ 24.856,24.908,24.958,25.006,25.052,25.096,25.138,25.179,
     $ 25.219, 25.257, 25.294, 25.331
     $/
c...35. lambda sigma0 => xi0 n
      data (sfits2(35,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000,75.229,35.616,26.088,21.961,19.573,17.847,16.461,
     $ 15.291,14.273,13.380,12.597,11.906,11.299,10.766,10.297,
     $  9.889, 9.535, 9.108, 8.667, 8.214, 7.783, 7.395, 7.058,
     $  6.774, 6.534, 6.333, 6.165, 6.023, 5.902, 5.798, 5.709,
     $  5.631, 5.563, 5.504, 5.451, 5.405, 5.363, 5.327, 5.295,
     $  5.267, 5.242, 5.221, 5.204, 5.189, 5.177, 5.168, 5.161,
     $  5.157, 5.156, 5.156, 5.158, 5.162, 5.167, 5.173, 5.180,
     $  5.187, 5.195, 5.203, 5.209, 5.216, 5.221, 5.225, 5.227,
     $  5.228, 5.227, 5.224, 5.220, 5.213, 5.205, 5.195, 5.183,
     $  5.170, 5.155, 5.140, 5.123, 5.106, 5.088, 5.070, 5.051,
     $  5.032,  5.013,  4.994,  4.975
     $/
c...36. lambda sigma0 => xi- p
      data (sfits2(36,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000,77.251,36.336,26.489,22.241,19.801,18.048,16.646,
     $ 15.463,14.433,13.529,12.736,12.034,11.417,10.875,10.397,
     $  9.982, 9.619, 9.180, 8.729, 8.268, 7.831, 7.437, 7.094,
     $  6.806, 6.562, 6.359, 6.188, 6.045, 5.922, 5.818, 5.728,
     $  5.649, 5.581, 5.520, 5.467, 5.420, 5.378, 5.341, 5.308,
     $  5.279, 5.254, 5.232, 5.213, 5.198, 5.185, 5.175, 5.167,
     $  5.162, 5.159, 5.158, 5.159, 5.162, 5.166, 5.171, 5.177,
     $  5.184, 5.190, 5.197, 5.203, 5.209, 5.213, 5.217, 5.219,
     $  5.219, 5.218, 5.215, 5.210, 5.203, 5.195, 5.185, 5.173,
     $  5.160, 5.146, 5.130, 5.114, 5.097, 5.079, 5.061, 5.042,
     $  5.023,  5.004,  4.985,  4.967
     $/
c...37. lambda sigma0 => s+s-
      data (sfits2(37,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 1.303, 1.790, 2.184, 2.521, 2.791, 2.992,
     $  3.129, 3.217, 3.264, 3.283, 3.282, 3.265, 3.238, 3.204,
     $  3.165, 3.123, 3.079, 3.035, 2.991, 2.948, 2.906, 2.866,
     $  2.829, 2.793, 2.761, 2.731, 2.704, 2.680, 2.660, 2.643,
     $  2.630, 2.620, 2.614, 2.611, 2.611, 2.615, 2.621, 2.630,
     $  2.641, 2.654, 2.670, 2.686, 2.703, 2.721, 2.738, 2.756,
     $  2.772, 2.788, 2.802, 2.814, 2.825, 2.833, 2.839, 2.843,
     $  2.845, 2.844, 2.840, 2.835, 2.827, 2.818, 2.806, 2.793,
     $  2.778,  2.761,  2.743,  2.724
     $/
c...38. lambda sigma+ total(mb)
      data (sfits2(38,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $ 302.839,122.652,82.533,67.002,59.742,55.811,53.212,51.233,
     $ 49.567,48.074,46.712,45.471,44.344,43.338,42.463,41.731,
     $ 41.541,42.458,42.004,41.356,40.664,39.989,39.368,38.816,
     $ 38.342,37.937,37.597,37.310,37.070,36.869,36.699,36.556,
     $ 36.435,36.332,36.245,36.172,36.111,36.061,36.021,35.992,
     $ 35.971,35.960,35.959,35.967,35.984,36.011,36.047,36.093,
     $ 36.148,36.212,36.284,36.365,36.453,36.547,36.646,36.750,
     $ 36.856,36.965,37.074,37.182,37.288,37.391,37.489,37.582,
     $ 37.669,37.748,37.821,37.885,37.941,37.989,38.030,38.062,
     $ 38.088,38.106,38.118,38.124,38.125,38.121,38.112,38.101,
     $ 38.086, 38.068, 38.049, 38.029
     $/
c...39. lambda sigma+ elastic
      data (sfits2(39,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $ 68.049,36.517,24.419,19.795,18.389,18.353,18.791,19.342,
     $ 19.847,20.256,20.564,20.787,20.951,21.077,21.189,21.315,
     $ 21.682,22.208,22.256,22.205,22.105,21.990,21.887,21.809,
     $ 21.760,21.739,21.742,21.766,21.806,21.858,21.919,21.986,
     $ 22.058,22.132,22.208,22.284,22.360,22.436,22.511,22.585,
     $ 22.659,22.732,22.804,22.877,22.949,23.022,23.096,23.170,
     $ 23.244,23.320,23.397,23.474,23.553,23.632,23.712,23.792,
     $ 23.873,23.953,24.033,24.113,24.191,24.268,24.344,24.417,
     $ 24.489,24.559,24.626,24.691,24.753,24.813,24.871,24.926,
     $ 24.979,25.029,25.078,25.124,25.169,25.212,25.254,25.294,
     $ 25.333, 25.370, 25.407, 25.443
     $/
c...40. lambda sigma+  => xi0 p
      data (sfits2(40,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $ 234.790,86.135,58.114,47.206,41.353,37.458,34.422,31.891,
     $ 29.720,27.818,26.148,24.684,23.393,22.262,21.274,20.415,
     $ 19.664,18.680,17.719,16.755,15.842,15.024,14.315,13.711,
     $ 13.206,12.781,12.424,12.122,11.865,11.644,11.453,11.287,
     $ 11.142,11.014,10.900,10.800,10.711,10.632,10.562,10.501,
     $ 10.447,10.400,10.361,10.327,10.300,10.279,10.263,10.252,
     $ 10.246,10.245,10.247,10.254,10.263,10.275,10.288,10.303,
     $ 10.319,10.335,10.350,10.364,10.376,10.385,10.393,10.397,
     $ 10.397,10.395,10.388,10.378,10.364,10.347,10.327,10.303,
     $ 10.277,10.249,10.218,10.185,10.151,10.116,10.081,10.044,
     $ 10.008,  9.971,  9.935,  9.899
     $/
c...41. lambda sigma+  => s(+) s(0)
      data (sfits2(41,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.196, 1.570, 2.029, 2.396, 2.717, 2.974, 3.165, 3.296,
     $  3.376, 3.418, 3.430, 3.422, 3.400, 3.367, 3.327, 3.283,
     $  3.235, 3.187, 3.137, 3.088, 3.040, 2.993, 2.948, 2.906,
     $  2.866, 2.828, 2.794, 2.763, 2.735, 2.710, 2.689, 2.671,
     $  2.657, 2.647, 2.640, 2.637, 2.637, 2.640, 2.646, 2.654,
     $  2.664, 2.677, 2.691, 2.706, 2.721, 2.737, 2.753, 2.768,
     $  2.782, 2.795, 2.806, 2.816, 2.824, 2.829, 2.832, 2.833,
     $  2.832, 2.828, 2.822, 2.814, 2.804, 2.792, 2.778, 2.762,
     $  2.745,  2.727,  2.708,  2.687
     $/
c...42. s(-) s(-) => s(-) s(-)
      data (sfits2(42,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000,319.483,175.764,103.954,70.316,51.063,
     $ 39.040,30.799,25.030,20.822,17.738,15.483,13.851,12.729,
     $ 11.995,11.603,11.470,11.563,11.831,12.241,12.757,13.350,
     $ 13.992,14.661,15.337,16.002,16.646,17.257,17.831,18.359,
     $ 18.844,19.280,19.673,20.022,20.331,20.602,20.840,21.046,
     $ 21.225,21.380,21.512,21.626,21.723,21.807,21.878,21.938,
     $ 21.990,22.034,22.072,22.104,22.132,22.157,22.178,22.197,
     $ 22.214,22.230,22.244,22.257,22.270,22.283,22.295,22.308,
     $ 22.321, 22.334, 22.348, 22.362
     $/
c...43. sigma0 sigma- total(mb)
      data (sfits2(43,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,855.254,273.486,163.394,116.470,90.946,75.161,
     $ 64.738,57.278,51.837,47.725,44.611,42.260,40.501,39.239,
     $ 38.362,37.827,37.554,37.509,37.638,37.909,38.281,38.728,
     $ 39.221,39.740,40.265,40.782,41.280,41.751,42.191,42.593,
     $ 42.960,43.289,43.585,43.846,44.078,44.282,44.461,44.618,
     $ 44.754,44.872,44.973,45.058,45.128,45.183,45.225,45.254,
     $ 45.269,45.271,45.259,45.235,45.197,45.147,45.085,45.010,
     $ 44.924,44.827,44.719,44.603,44.477,44.344,44.204,44.059,
     $ 43.908, 43.753, 43.594, 43.433
     $/
c...44. sigma0 sigma- elastic
      data (sfits2(44,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,472.366,200.539,117.257,79.877,59.496,47.197,
     $ 39.440,34.191,30.626,28.141,26.446,25.341,24.682,24.393,
     $ 24.386,24.624,25.047,25.626,26.320,27.101,27.938,28.805,
     $ 29.681,30.546,31.386,32.189,32.945,33.649,34.297,34.886,
     $ 35.418,35.892,36.315,36.685,37.009,37.290,37.532,37.738,
     $ 37.912,38.057,38.174,38.269,38.341,38.394,38.430,38.449,
     $ 38.453,38.444,38.422,38.389,38.346,38.293,38.231,38.161,
     $ 38.083,37.999,37.910,37.815,37.716,37.614,37.508,37.400,
     $ 37.291, 37.180, 37.069, 36.957
     $/
c...45. sigma0 sigma- => xi- n
      data (sfits2(45,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,21.425, 8.135, 8.058, 8.216, 8.164, 7.912,
     $  7.525, 7.073, 6.603, 6.144, 5.712, 5.313, 4.950, 4.623,
     $  4.328, 4.065, 3.829, 3.619, 3.430, 3.262, 3.112, 2.978,
     $  2.859, 2.754, 2.661, 2.580, 2.509, 2.448, 2.397, 2.354,
     $  2.319, 2.293, 2.274, 2.261, 2.256, 2.257, 2.263, 2.275,
     $  2.292, 2.312, 2.337, 2.365, 2.395, 2.427, 2.461, 2.495,
     $  2.529, 2.563, 2.596, 2.628, 2.658, 2.685, 2.711, 2.734,
     $  2.754, 2.771, 2.785, 2.797, 2.806, 2.812, 2.816, 2.817,
     $  2.816,  2.813,  2.808,  2.802
     $/
c...46. sigma0 sigma- => lambda sigma-
      data (sfits2(46,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,361.463,64.812,38.078,28.376,23.285,20.051,
     $ 17.774,16.014,14.608,13.440,12.453,11.607,10.868,10.222,
     $  9.648, 9.138, 8.678, 8.264, 7.888, 7.545, 7.232, 6.945,
     $  6.681, 6.440, 6.217, 6.013, 5.826, 5.655, 5.497, 5.354,
     $  5.223, 5.104, 4.996, 4.900, 4.813, 4.735, 4.666, 4.605,
     $  4.551, 4.503, 4.461, 4.424, 4.391, 4.362, 4.335, 4.310,
     $  4.286, 4.263, 4.240, 4.217, 4.194, 4.169, 4.143, 4.116,
     $  4.087, 4.057, 4.025, 3.991, 3.955, 3.919, 3.880, 3.841,
     $  3.801,  3.759,  3.717,  3.674
     $/
c...47. sigma0 sigma0 total
      data (sfits2(47,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000,249.633,133.727,87.052,65.163,52.596,44.083,37.802,
     $ 33.024,29.250,26.274,23.889,21.992,20.497,19.330,18.452,
     $ 17.804,17.364,17.085,16.947,16.918,16.977,17.100,17.271,
     $ 17.472,17.690,17.914,18.134,18.345,18.541,18.720,18.879,
     $ 19.019,19.139,19.241,19.326,19.395,19.451,19.494,19.527,
     $ 19.551,19.568,19.579,19.586,19.589,19.589,19.588,19.586,
     $ 19.583,19.581,19.579,19.579,19.579,19.581,19.585,19.590,
     $ 19.597,19.606,19.618,19.631,19.646,19.663,19.682,19.703,
     $ 19.726, 19.750, 19.776, 19.804
     $/
c...48. sigma0 sigma0 elastic
      data (sfits2(48,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 6.328, 4.863, 3.156, 2.409, 1.983, 1.698, 1.493,
     $  1.340, 1.224, 1.134, 1.063, 1.009, 0.966, 0.932, 0.907,
     $  0.886, 0.871, 0.859, 0.849, 0.842, 0.836, 0.831, 0.826,
     $  0.821, 0.817, 0.812, 0.807, 0.802, 0.797, 0.791, 0.785,
     $  0.778, 0.771, 0.764, 0.757, 0.749, 0.741, 0.734, 0.726,
     $  0.718, 0.710, 0.703, 0.695, 0.688, 0.680, 0.673, 0.666,
     $  0.660, 0.653, 0.647, 0.641, 0.636, 0.631, 0.626, 0.621,
     $  0.616, 0.612, 0.609, 0.605, 0.602, 0.599, 0.596, 0.594,
     $  0.592,  0.591,  0.589,  0.588
     $/
c...49. sigma0 sigma0 => lam lam
      data (sfits2(49,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 3.769, 2.759, 1.656, 1.179, 0.925, 0.768, 0.662,
     $  0.589, 0.536, 0.498, 0.469, 0.448, 0.433, 0.422, 0.415,
     $  0.410, 0.407, 0.406, 0.406, 0.406, 0.407, 0.408, 0.409,
     $  0.410, 0.411, 0.411, 0.411, 0.411, 0.411, 0.410, 0.408,
     $  0.407, 0.405, 0.403, 0.401, 0.399, 0.397, 0.394, 0.392,
     $  0.389, 0.386, 0.384, 0.381, 0.379, 0.376, 0.374, 0.372,
     $  0.370, 0.368, 0.366, 0.364, 0.363, 0.361, 0.360, 0.359,
     $  0.358, 0.357, 0.356, 0.356, 0.355, 0.355, 0.355, 0.355,
     $  0.355,  0.356,  0.356,  0.356
     $/
c...50. sigma0 sigma0 => xi0 n
      data (sfits2(50,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 3.798, 2.777, 1.664, 1.183, 0.927, 0.768, 0.662,
     $  0.588, 0.534, 0.495, 0.466, 0.445, 0.430, 0.419, 0.411,
     $  0.406, 0.404, 0.402, 0.402, 0.402, 0.403, 0.404, 0.406,
     $  0.407, 0.407, 0.408, 0.408, 0.408, 0.408, 0.407, 0.406,
     $  0.404, 0.403, 0.401, 0.399, 0.397, 0.394, 0.392, 0.389,
     $  0.387, 0.384, 0.382, 0.379, 0.377, 0.374, 0.372, 0.370,
     $  0.368, 0.366, 0.364, 0.362, 0.361, 0.359, 0.358, 0.357,
     $  0.356, 0.355, 0.354, 0.354, 0.353, 0.353, 0.353, 0.353,
     $  0.353,  0.353,  0.354,  0.354
     $/
c...51. sigma0 sigma0 => xi- p
      data (sfits2(51,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000,  0.000,  0.000,  0.000
     $/
c...52. sigma0 sigma0 => lam s0
      data (sfits2(52,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,28.496,21.114,18.547,17.070,15.759,14.510,
     $ 13.357,12.318,11.402,10.595, 9.889, 9.274, 8.737, 8.270,
     $  7.860, 7.504, 7.191, 6.916, 6.672, 6.456, 6.262, 6.087,
     $  5.928, 5.781, 5.644, 5.516, 5.393, 5.276, 5.162, 5.052,
     $  4.945, 4.841, 4.739, 4.639, 4.541, 4.446, 4.352, 4.261,
     $  4.172, 4.085, 4.000, 3.917, 3.837, 3.759, 3.683, 3.609,
     $  3.538, 3.468, 3.401, 3.336, 3.273, 3.212, 3.152, 3.095,
     $  3.040, 2.986, 2.934, 2.884, 2.835, 2.788, 2.743, 2.699,
     $  2.657,  2.616,  2.576,  2.538
     $/
c...53. sigma0 sigma0 => s+ s-
      data (sfits2(53,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000,235.738,94.831,59.461,41.845,31.691,25.090,20.475,
     $ 17.149,14.638,12.746,11.296,10.201, 9.394, 8.820, 8.450,
     $  8.241, 8.178, 8.227, 8.374, 8.595, 8.875, 9.195, 9.543,
     $  9.906,10.274,10.638,10.992,11.331,11.650,11.950,12.228,
     $ 12.484,12.719,12.934,13.130,13.309,13.473,13.622,13.760,
     $ 13.886,14.003,14.112,14.213,14.309,14.400,14.486,14.569,
     $ 14.649,14.726,14.801,14.875,14.947,15.019,15.089,15.158,
     $ 15.228,15.296,15.364,15.432,15.500,15.567,15.634,15.701,
     $ 15.768, 15.835, 15.902, 15.968
     $/
c...54. sigma+ sigma- total(mb)
      data (sfits2(54,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,172.640,112.283,93.250,84.536,78.924,74.501,
     $ 70.706,67.307,64.250,61.489,59.012,56.804,54.845,53.126,
     $ 51.615,50.307,49.167,48.186,47.335,46.602,45.966,45.414,
     $ 44.929,44.502,44.120,43.778,43.466,43.181,42.917,42.673,
     $ 42.445,42.232,42.033,41.847,41.674,41.513,41.363,41.224,
     $ 41.095,40.975,40.862,40.757,40.657,40.561,40.468,40.377,
     $ 40.286,40.195,40.102,40.006,39.908,39.805,39.698,39.587,
     $ 39.472,39.352,39.228,39.101,38.970,38.836,38.699,38.561,
     $ 38.421, 38.280, 38.138, 37.997
     $/
c...55. sigma+ sigma- elastic
      data (sfits2(55,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,65.049,46.170,39.282,36.574,35.283,34.485,
     $ 33.850,33.248,32.643,32.038,31.448,30.890,30.376,29.918,
     $ 29.518,29.180,28.900,28.675,28.498,28.363,28.264,28.193,
     $ 28.145,28.113,28.094,28.083,28.077,28.073,28.069,28.065,
     $ 28.059,28.051,28.040,28.028,28.014,27.998,27.980,27.962,
     $ 27.941,27.920,27.898,27.874,27.849,27.823,27.795,27.765,
     $ 27.733,27.699,27.662,27.623,27.582,27.538,27.492,27.443,
     $ 27.392,27.340,27.285,27.228,27.171,27.112,27.052,26.991,
     $ 26.931, 26.870, 26.809, 26.748
     $/
c...56. sigma+ sigma- => ll
      data (sfits2(56,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,14.227, 7.649, 5.457, 4.367, 3.688, 3.216,
     $  2.874, 2.612, 2.411, 2.253, 2.129, 2.031, 1.953, 1.893,
     $  1.845, 1.807, 1.778, 1.755, 1.737, 1.722, 1.709, 1.698,
     $  1.688, 1.678, 1.669, 1.659, 1.648, 1.637, 1.625, 1.613,
     $  1.599, 1.586, 1.571, 1.557, 1.542, 1.526, 1.511, 1.495,
     $  1.479, 1.464, 1.449, 1.433, 1.419, 1.404, 1.390, 1.376,
     $  1.363, 1.350, 1.338, 1.326, 1.315, 1.304, 1.294, 1.285,
     $  1.276, 1.268, 1.261, 1.254, 1.247, 1.241, 1.236, 1.232,
     $  1.227,  1.224,  1.221,  1.219
     $/
c...57. sigma+ sigma- => xi0 n
      data (sfits2(57,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,10.970, 7.253, 6.189, 5.656, 5.247, 4.875,
     $  4.526, 4.200, 3.904, 3.637, 3.400, 3.192, 3.009, 2.850,
     $  2.710, 2.589, 2.484, 2.391, 2.311, 2.240, 2.178, 2.124,
     $  2.076, 2.034, 1.997, 1.964, 1.936, 1.912, 1.891, 1.874,
     $  1.860, 1.850, 1.842, 1.837, 1.835, 1.835, 1.838, 1.843,
     $  1.850, 1.859, 1.869, 1.881, 1.894, 1.907, 1.922, 1.936,
     $  1.950, 1.965, 1.978, 1.992, 2.004, 2.015, 2.026, 2.035,
     $  2.043, 2.050, 2.055, 2.060, 2.063, 2.066, 2.067, 2.067,
     $  2.067,  2.066,  2.064,  2.062
     $/
c...58. sigma+ sigma- => xi- p
      data (sfits2(58,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,11.065, 7.316, 6.242, 5.703, 5.289, 4.912,
     $  4.559, 4.229, 3.928, 3.658, 3.418, 3.206, 3.021, 2.859,
     $  2.718, 2.595, 2.487, 2.394, 2.312, 2.240, 2.177, 2.122,
     $  2.073, 2.030, 1.992, 1.959, 1.930, 1.905, 1.884, 1.866,
     $  1.852, 1.840, 1.832, 1.826, 1.823, 1.823, 1.825, 1.829,
     $  1.836, 1.844, 1.854, 1.865, 1.877, 1.890, 1.904, 1.918,
     $  1.932, 1.946, 1.960, 1.972, 1.985, 1.996, 2.006, 2.015,
     $  2.024, 2.031, 2.036, 2.041, 2.045, 2.047, 2.049, 2.050,
     $  2.050,  2.049,  2.047,  2.045
     $/
c...59. sigma+ sigma- => l s0
      data (sfits2(59,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,29.277,18.326,14.937,13.294,12.187,11.291,
     $ 10.511, 9.809, 9.178, 8.610, 8.099, 7.641, 7.228, 6.857,
     $  6.522, 6.219, 5.944, 5.695, 5.468, 5.262, 5.073, 4.902,
     $  4.746, 4.604, 4.475, 4.359, 4.253, 4.159, 4.074, 4.000,
     $  3.934, 3.876, 3.827, 3.785, 3.750, 3.721, 3.699, 3.681,
     $  3.668, 3.659, 3.654, 3.651, 3.650, 3.650, 3.652, 3.653,
     $  3.655, 3.655, 3.654, 3.652, 3.647, 3.641, 3.632, 3.621,
     $  3.607, 3.591, 3.572, 3.552, 3.529, 3.504, 3.477, 3.449,
     $  3.419,  3.388,  3.355,  3.322
     $/
c...60. sigma+ sigma- => s0s0
      data (sfits2(60,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000,42.053,25.569,21.142,18.942,17.231,15.722,
     $ 14.387,13.209,12.186,11.293,10.518, 9.845, 9.258, 8.749,
     $  8.303, 7.916, 7.574, 7.276, 7.011, 6.776, 6.565, 6.375,
     $  6.202, 6.042, 5.893, 5.754, 5.622, 5.495, 5.373, 5.256,
     $  5.141, 5.030, 4.921, 4.815, 4.711, 4.610, 4.511, 4.415,
     $  4.321, 4.229, 4.140, 4.053, 3.968, 3.886, 3.806, 3.729,
     $  3.653, 3.581, 3.510, 3.442, 3.375, 3.311, 3.249, 3.188,
     $  3.130, 3.074, 3.019, 2.966, 2.915, 2.866, 2.818, 2.771,
     $  2.727,  2.684,  2.642,  2.601
     $/
c...61. sigma+ sigma0 total(mb)
      data (sfits2(61,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000,356.195,190.075,129.958,98.545,80.135,68.129,59.738,
     $ 53.689,49.126,45.695,43.078,41.113,39.672,38.649,37.985,
     $ 37.602,37.469,37.522,37.733,38.059,38.473,38.944,39.450,
     $ 39.971,40.490,40.997,41.481,41.935,42.356,42.741,43.090,
     $ 43.404,43.683,43.931,44.151,44.345,44.515,44.664,44.793,
     $ 44.905,45.001,45.080,45.146,45.197,45.235,45.260,45.272,
     $ 45.270,45.255,45.228,45.187,45.134,45.068,44.991,44.902,
     $ 44.803,44.694,44.575,44.448,44.314,44.173,44.027,43.876,
     $ 43.721, 43.563, 43.402, 43.240
     $/
c...62. sigma+ sigma0 elastic
      data (sfits2(62,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000,257.757,138.172,90.678,65.518,51.010,41.899,35.857,
     $ 31.783,28.929,26.981,25.668,24.851,24.426,24.309,24.451,
     $ 24.794,25.309,25.950,26.692,27.500,28.351,29.220,30.087,
     $ 30.936,31.754,32.531,33.259,33.933,34.550,35.111,35.614,
     $ 36.064,36.461,36.811,37.115,37.379,37.605,37.798,37.960,
     $ 38.094,38.204,38.289,38.355,38.402,38.432,38.447,38.447,
     $ 38.435,38.410,38.374,38.329,38.273,38.210,38.138,38.060,
     $ 37.975,37.885,37.790,37.691,37.588,37.483,37.376,37.267,
     $ 37.157, 37.046, 36.936, 36.825
     $/
c...63. sigma+ sigma0 => xi0 p
      data (sfits2(63,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 8.733, 7.949, 8.135, 8.171, 7.987, 7.647, 7.221,
     $  6.760, 6.299, 5.861, 5.452, 5.079, 4.741, 4.435, 4.162,
     $  3.916, 3.697, 3.501, 3.325, 3.168, 3.029, 2.904, 2.793,
     $  2.695, 2.609, 2.533, 2.468, 2.412, 2.365, 2.327, 2.297,
     $  2.274, 2.258, 2.249, 2.246, 2.250, 2.259, 2.273, 2.291,
     $  2.314, 2.340, 2.369, 2.400, 2.433, 2.466, 2.501, 2.535,
     $  2.569, 2.601, 2.632, 2.662, 2.689, 2.713, 2.735, 2.754,
     $  2.771, 2.784, 2.795, 2.803, 2.809, 2.812, 2.813, 2.811,
     $  2.808,  2.802,  2.795,  2.787
     $/
c...64. sigma+ sigma0 => lam s+
      data (sfits2(64,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000,89.705,43.954,31.145,24.856,21.138,18.583,16.660,
     $ 15.147,13.898,12.853,11.958,11.183,10.505, 9.904, 9.372,
     $  8.892, 8.462, 8.071, 7.716, 7.391, 7.093, 6.820, 6.570,
     $  6.339, 6.127, 5.933, 5.754, 5.590, 5.440, 5.303, 5.179,
     $  5.066, 4.964, 4.872, 4.789, 4.716, 4.651, 4.593, 4.542,
     $  4.497, 4.457, 4.422, 4.391, 4.362, 4.336, 4.312, 4.289,
     $  4.266, 4.244, 4.221, 4.197, 4.172, 4.146, 4.118, 4.088,
     $  4.057, 4.025, 3.990, 3.954, 3.917, 3.879, 3.839, 3.798,
     $  3.756,  3.714,  3.671,  3.628
     $/
c...65. s(+) s(+) => s(+) s(+)
      data (sfits2(65,i),i=1,itblsz)/
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     $ 308.657,148.960,92.007,64.086,47.336,36.569,29.154,23.833,
     $ 19.984,17.123,15.044,13.541,12.506,11.850,11.498,11.406,
     $ 11.517,11.803,12.218,12.738,13.327,13.965,14.628,15.297,
     $ 15.958,16.596,17.205,17.775,18.304,18.788,19.228,19.623,
     $ 19.977,20.290,20.567,20.809,21.021,21.205,21.365,21.503,
     $ 21.622,21.724,21.811,21.886,21.950,22.005,22.053,22.093,
     $ 22.128,22.158,22.184,22.207,22.228,22.246,22.262,22.277,
     $ 22.291,22.305,22.318,22.330,22.343,22.355,22.368,22.382,
     $ 22.396, 22.410, 22.425, 22.441
     $/

      end
