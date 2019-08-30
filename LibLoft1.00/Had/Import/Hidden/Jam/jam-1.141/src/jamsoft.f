c***********************************************************************
c***********************************************************************
c                                                                      *
c        PART 5: Soft interaction: String excitation and fragmentation *
c                                                                      *
c   List of subprograms in rough order of relevance with main purpose  *
c      (S = subroutine, F = function, B = block data, E = entry)       *
c                                                                      *
c  s  jamsoft  to give masses to the excited nucleons                  *
c  s  jamasft  to set particles after soft interaction                 *
c  s  jamsfex  to generate the soft interaction for each binary coll   *
c  s  jamdiffr to calculate diffractive cross section                  *
c  s  jamjetm1 to save parton properties for jet system                *
c  s  jamjetm2 to save parton properties after annihilation            *
c  s  jamdmass to give minimum and min. excited mass for the particle  *
c  s  jamjdec  to fragment jet system by Lund model                    *
c  s  jamfrg   to perform string fragmentation                         *
c                                                                      *
c  s  jamidres to determine resonance ID corresponding to mass         *
c  s  jamkfres to determine baryon resonance ID corresponding to mass  *
c  s  jamrobo  to perform rotations and boosts                         *
c  s  jamindb  to find index of newly produced baryon/antibaryon       *
c  f  xsamp1   to sample x of valence quarks for baryon, DPM type      *
c  f  xsamp2   to sample x of valence quarks for mseon DPM type        *
c  f  xsamp3   to sample x for diffractive scattering                  *
c  s  attrad   to perform soft radiations by Lund dipole approx.       *
c  s  ar3jet   to calculate scaled energy variables of gluon and quark *
c  s  arorie   to give pt to the soft radiative gluon                  *
c  s  atrobo   to make boost and rotation to entries                   *
c  s  str2had  to check mass of string system                          *
c  s  attflv   to give spin and quarkflavour to the leading particles  *
c  f  ifrkfcpj   to return the kf code for flavor having ia ib ic        *
c  f  kfprop   to give charge, baryon number and strangeness           *
c  s  kfcnst   to construct kf code from flavor pair                   *
c  s  attflv2  to give quark content                                   *
c  f  jamflav  to give baron number or flavour content of particles    *
c                                                                      *
c***********************************************************************
c***********************************************************************

      subroutine jamsoft(icltyp,jp,jt,icon)

c...Purpose: to collide two-particles to form strings.
      include 'jam1.inc'
      include 'jam2.inc'
      real*8 jamemjet
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      common/hijjet5/nfp(2,10)

      dimension p1(15),p2(15),pcm(5)

c...Initialize for a new event
      icon=0
      mste(4)=0  ! flag for diffractive scatt.

      kf1=k(2,jp)
      kf2=k(2,jt)
      kf01=kf1
      kf02=kf2
      if(kf1.ne.92) then
        do i=1,10
          vq(i,jp)=0.0d0
        end do
        do i=1,2
          kq(i,jp)=0
        enddo
      endif
      if(kf2.ne.92) then
        do i=1,10
          vq(i,jt)=0.0d0
        end do
        do i=1,2
          kq(i,jt)=0
        enddo
      endif

      do j=1,15
        p1(j)=0.0d0
        p2(j)=0.0d0
      end do
      do i=1,10
        nfp(1,i)=0
        nfp(2,i)=0
      end do

      jp0=jp
      jt0=jt
      if(k(10,jp).eq.0) then
        ipart1=0
        jpp=0
        ks1=k(1,jp)
        kf1=k(2,jp)
        if(kf1.ne.92) then
          kc1=jamcomp(kf1)
          call attflv(kf1,jq1,jqq1)
          p1(1)=p(1,jp)
          p1(2)=p(2,jp)
          p1(3)=p(3,jp)
          p1(4)=p(4,jp)
          p1(5)=p(5,jp)
        else
          jq1=kq(1,jp)
          jqq1=kq(2,jp)
          p1(1)=vq(1,jp)+vq(6,jp)
          p1(2)=vq(2,jp)+vq(7,jp)
          p1(3)=vq(3,jp)+vq(8,jp)
          p1(4)=vq(4,jp)+vq(9,jp)
          p1(5)=sqrt(p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2)
          call kfcnst(jq1,jqq1,kf1,0.0d0)
          if(kf1.eq.0)then
            call jamerrm(30,0,'(jamsoft:)kf1=0 at k1=4')
          endif
        endif
      else 
c....jp is a parton system
        ipart1=1
        kf1=0
        ks1=3
        kc1=0
        jp1=k(10,jp)
        jp2=k(11,jp)
        if(jp1.eq.jp) then
          jpp=jp2
        else if(jp2.eq.jp) then
          jpp=jp1
        else
          write(check(1),8000)jp1,jp2,jp,k(2,jp1),k(2,jp2)
 8000     format('jp1 jp2 jp q1 q2',3i10)
          call jamerrm(30,1,'(jamsoft:) ???? jp1 jp2 jp')
        endif
        if(abs(k(2,jp1)).ge.1.and.abs(k(2,jp1)).le.10) then
          jq1=k(2,jp1)
          jqq1=k(2,jp2)
        else
          jq1=k(2,jp2)
          jqq1=k(2,jp1)
        endif
        do j=1,4
        p1(j)=p(j,jp1)+p(j,jp2)
        end do
        p1(5)=sqrt(p1(4)**2-(p1(1)**2+p1(2)**2+p1(3)**2))
      endif

      if(k(10,jt).eq.0) then
        ipart2=0
        jtt=0
        ks2=k(1,jt)
        kf2=k(2,jt)
        if(kf2.ne.92)then
          kc2=jamcomp(kf2)
          call attflv(kf2,jq2,jqq2)
          p2(1)=p(1,jt)
          p2(2)=p(2,jt)
          p2(3)=p(3,jt)
          p2(4)=p(4,jt)
          p2(5)=p(5,jt)
        else
          jq2=kq(1,jt)
          jqq2=kq(2,jt)
          p2(1)=vq(1,jt)+vq(6,jt)
          p2(2)=vq(2,jt)+vq(7,jt)
          p2(3)=vq(3,jt)+vq(8,jt)
          p2(4)=vq(4,jt)+vq(9,jt)
          p2(5)=sqrt(p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2)
          call kfcnst(jq2,jqq2,kf2,0.0d0)
          if(kf2.eq.0)then
            call jamerrm(30,0,'(jamsoft:)kf2=0 at k1=4')
          endif
        endif
      else
        ipart2=1
        kf2=0
        ks2=3
        kc2=0
        jt1=k(10,jt)
        jt2=k(11,jt)
        if(jt1.eq.jt) then
          jtt=jt2
        else if(jt2.eq.jt) then
          jtt=jt1
        else
          write(check(1),8010)jt1,jt2,jt,k(2,jt1),k(2,jt2)
 8010     format('jt1 jt2 jt q1 q2',5i10)
          call jamerrm(30,1,'(jamsoft:) ???? jt1 jt2 jt')
        endif
        if(abs(k(2,jt1)).ge.1.and.abs(k(2,jt1)).le.10) then
          jq2=k(2,jt1)
          jqq2=k(2,jt2)
        else
          jq2=k(2,jt2)
          jqq2=k(2,jt1)
        endif
        do j=1,4
        p2(j)=p(j,jt1)+p(j,jt2)
        end do
        p2(5)=sqrt(p2(4)**2-(p2(1)**2+p2(2)**2+p2(3)**2))
      endif

      em01=p1(5)
      em02=p2(5)
      jq01=jq1
      jqq01=jqq1
      jq02=jq2
      jqq02=jqq2
      do j=1,4
        pcm(j)=p1(j)+p2(j)
      end do
      pcm(5)=sqrt(pcm(4)**2-(pcm(1)**2+pcm(2)**2+pcm(3)**2))

c...Lorentz-transformation in two-body c.m. system.

      bex=pcm(1)/pcm(4)
      bey=pcm(2)/pcm(4)
      bez=pcm(3)/pcm(4)
      gam=pcm(4)/pcm(5)
      call jamrobo(0.0d0,0.0d0,-bex,-bey,-bez,gam,p1(1),p1(2),p1(3), 
     & p1(4))
      call jamrobo(0.0d0,0.0d0,-bex,-bey,-bez,gam,p2(1),p2(2),p2(3), 
     & p2(4))

c...Rotation so that p_x=p_y=0.
      phi=pjangl(p1(1),p1(2))
      call jamrobo(0.0d0,-phi,0.0d0,0.0d0,0.0d0,1.0d0,p1(1),p1(2),p1(3),
     & p1(4))
      call jamrobo(0.0d0,-phi,0.0d0,0.0d0,0.0d0,1.0d0,p2(1),p2(2),p2(3),
     & p2(4))
      theta=pjangl(p1(3),p1(1))
      call jamrobo(-theta,0.d0,0.d0,0.d0,0.d0,1.0d0,p1(1),p1(2),p1(3), 
     & p1(4))
      call jamrobo(-theta,0.d0,0.d0,0.d0,0.d0,1.0d0,p2(1),p2(2),p2(3), 
     & p2(4))

      if(p1(3).lt.0.d0) then
        write(mstc(38),*)'p1',jp,k(1,jp),k(2,jp),(p1(j),j=1,5)
        write(mstc(38),*)'p2',jt,k(1,jt),k(2,jt),(p2(j),j=1,5)
        call jamerrm(1,0,'(jamsoft:) p1(3)<0????')
        icon=1
        return
      endif

c...Minimum jet masses.
      emjet1=jamemjet(jq1,jqq1)
      emjet2=jamemjet(jq2,jqq2)

c==============================================
c...Conduct soft scattering between JP and JT
c==============================================
 160  continue

      kf10=kf01
      kf20=kf02
      call jamsfex(kc1,kc2,kf10,kf20,kf1,kf2,jq1,jqq1,jq2,jqq2
     $                 ,ks1,ks2,p1,p2
     $                 ,emjet1,emjet2
     $                 ,ipath,ierror)


      if(ierror.eq.999) then
        if(mstc(8).ge.1) write(mstc(38),*)'error occured in jamsfex999'
        icon=2
        return
      endif
      if(ierror.ne.0) then
        icon=ierror
        return
      endif

300   continue
c------------------------------------------------------------
c...Store particle propaties in JAM arrary.

      mstd(29)=mstd(29)+1
      call jamasft(1,jp,jp1,jp2,ipart1,kf10,jq1,jqq1,emjet1,p1,
     $                    bex,bey,bez,gam,phi,theta,kf1,kf01)

      if(ipart1.eq.0) then
        kc=jamcomp(kf1)
        call jamkupda(1,jp,kf1,kc,ks1,k(6,jt),mstd(29),5)
      else
        k(5,jp)=mstd(29)
      endif

      call jamasft(-1,jt,jt1,jt2,ipart2,kf20,jq2,jqq2,emjet2,p2,
     $                     bex,bey,bez,gam,phi,theta,kf2,kf02)
      if(ipart2.eq.0) then
        kc=jamcomp(kf2)
        call jamkupda(1,jt,kf2,kc,ks2,k(6,jp),mstd(29),5)
      else
        k(5,jt)=mstd(29)
      endif

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      if(k(2,jp).eq.92.and.p(5,jp).lt.emjet1) ierror=10
      if(k(2,jt).eq.92.and.p(5,jt).lt.emjet2) ierror=10

c...Check energy momentum conservation.
      srt1=sqrt((p1(4)+p2(4))**2-(p1(1)+p2(1))**2
     a        -(p1(2)+p2(2))**2
     a        -(p1(3)+p2(3))**2)

      if(     (abs( pcm(1)-(p1(1)+p2(1) ) ).gt.0.1d0)
     $    .or.(abs( pcm(2)-(p1(2)+p2(2) ) ).gt.0.1d0)
     $    .or.(abs( pcm(3)-(p1(3)+p2(3) ) ).gt.0.1d0)
     $    .or.(abs( pcm(4)-(p1(4)+p2(4) ) ).gt.0.1d0)
     $    .or.(abs( pcm(5)-srt1).gt.0.1d0 ) ) then

        call jamerrm(11,0,'(jamsoft:) energy mom. not conserved')
        write(mstc(38),*)'ipath em1 em2',ipath,p1(5),p2(5)
        write(mstc(38),*)'pcm1=',pcm(1),p1(1)+p2(1)
        write(mstc(38),*)'pcm2=',pcm(2),p1(2)+p2(2)
        write(mstc(38),*)'pcm3=',pcm(3),p1(3)+p2(3)
        write(mstc(38),*)'pcm4=',pcm(4),p1(4)+p2(4)
        write(mstc(38),*)'pcm5=',pcm(5),srt1
        ierror=10
      endif

c...Print information after soft interaction.
      if(mstc(8).ge.4.or.ierror.eq.10) then

        ih=mstc(38)
        write(ih,*)' '
        write(ih,*)'=== after JAMSOFT === srt mste(4)',
     $      pcm(5),mste(4),emjet1,emjet2,ipart1,ipart2,ipath
        if(ipart1.eq.1)write(ih,*)'k10 k11 jp',k(10,jp),k(11,jp),jp
        if(ipart2.eq.1)write(ih,*)'k10 k11jt',k(10,jt),k(11,jt),jt

        write(ih,800)jp0,ks1,kf01,em01,jq01,jqq01
        write(ih,801)jt0,ks2,kf02,em02,jq02,jqq02
        write(ih,802)jp,k(1,jp),k(2,jp),p1(5),kf10,
     $   jamk(1,jp),kq(1,jp),kq(2,jp)
        write(ih,803)jt,k(1,jt),k(2,jt),p2(5),kf20,
     $   jamk(1,jt),kq(1,jt),kq(2,jt)

        write(ih,*)'r1',(r(j,jp),j=1,5)
        write(ih,*)'r2',(r(j,jt),j=1,5)

800   format('jp0=',i5,' ks1',i3,' kf01',i6,' em01',
     $    f8.3,' jq1 jq2',i4,1x,i4)
801   format('jt0=',i5,' ks2',i3,' kf02',i6,' em02',
     $    f8.3,' jq1 jq2',i4,1x,i4)
802   format('jp= ',i5,' ks1',i3,' kf1 ',i6,' em1'
     $   ,f10.3,' kf=',i6,' ch=',i3,' jq1 jq2',i4,1x,i4)
803   format('jt= ',i5,' ks2',i3,' kf2 ',i6,' em2'
     $   ,f10.3,' kf=',i6,' ch=',i3,' jq1 jq2',i4,1x,i4)

        if(ierror.ne.0) call jamerrm(30,0,'Error after soft')
      endif

c...Update the collision array.
      if(mstc(6).ge.0) then
        if(jpp.ne.0) call jamcupda(jpp,-1,1)
        if(jtt.ne.0) call jamcupda(jtt,-1,1)
      endif

      end

c***********************************************************************

      subroutine jamasft(is,ip,j1,j2,ipart,kfv,jq,jqq,emj,pj,
     $  bex,bey,bez,gam,phi,the,kf,kf0)

c...Set particle momentum and properties after soft interaction.
      include 'jam1.inc'
      include 'jam2.inc'
      dimension pj(15)
      parameter(gamm=0.1d0)
      data ismear/2/

c...Decide whether Jet excitation or resonance.
      ijet=1
      if(ipart.eq.0.and.kf0.ne.92)then
        if(pj(5).le.emj) then
          ijet=0
        else
          kfx=kfv
          call jamidres(kfx,pj(5),icon)
          if(icon.ne.10) then
            if(ismear.ge.1)then
              if(ismear.eq.1) probr=2.d0/(1d0+exp((pj(5)-emj)**2/0.2d0))
              if(ismear.eq.2) probr=gamm**2/4.0d0/((pj(5)-emj)**2
     $                           + gamm**2/4.0d0)
              if(rn(0).le.probr) ijet=0
            endif
          endif
        endif
      endif

c...Back to the original frame:string.
      if(ijet.eq.1) then
        call jamjetm1(is,ip,pj,jq,jqq,bex,bey,bez,gam,phi,the)
c...Recalculate energies in order to achieve numerical accuracy.
        pj(1)=vq(1,ip)+vq(6,ip)
        pj(2)=vq(2,ip)+vq(7,ip)
        pj(3)=vq(3,ip)+vq(8,ip)
        pj(4)=sqrt(pj(5)**2+pj(1)**2+pj(2)**2+pj(3)**2)
      else
       call jamrobo(the,phi,bex,bey,bez,gam,pj(1),pj(2),pj(3),pj(4))
       pj(4)=sqrt(pj(5)**2+pj(1)**2+pj(2)**2+pj(3)**2)
      endif

c...Store particle propaties in JAM arrary.

      if(ipart.eq.0) then

        if(ijet.eq.1) then
          p(1,ip)=bex
          p(2,ip)=bey
          p(3,ip)=bez
          p(4,ip)=pj(4)
          p(5,ip)=pj(5)
          v(1,ip)=phi
          v(2,ip)=the
          v(3,ip)=gam
        else

          do j=1,5
            p(j,ip)=pj(j)
          end do
          kq(1,ip)=999999
          kq(2,ip)=0
          vq(1,ip)=bex
          vq(2,ip)=bey
          vq(3,ip)=bez
          vq(4,ip)=gam
          vq(5,ip)=phi
          vq(6,ip)=the
c....Determine resonance kf code.
          if(kf.eq.92) then
            kf=kfv
            call jamidres(kf,pj(5),icon)
            if(kf.eq.0) then
              write(check(1),'(''kf pj5'',i9,g12.3)')kf,pj(5)
              call jamerrm (30,1,'(jamasft:)after jamidres kf')
            endif
          endif

        endif


      else
        do j=1,5
          p(j,j1)=vq(j,ip)
          p(j,j2)=vq(j+5,ip)
        end do
        vq(1,j1)=bex
        vq(2,j1)=bey
        vq(3,j1)=bez
        vq(4,j1)=gam
        vq(5,j1)=phi
        vq(6,j1)=the
      endif

      end

c***********************************************************************

      subroutine jamsfex(kc1,kc2,kf10,kf20,kf1,kf2,jq1,jqq1,jq2,jqq2
     $                 ,ks01,ks02,p1,p2
     $                 ,emjet1,emjet2
     $                 ,ipath,ierror)

c...Purpose: to generate the soft interaction for each binary h-h
c...collision (same as hijsft() of HIJING model).
c=======================================================================
c mste(4)->
c        1=double diffrac
c        2=single diffrac proj.,
c        3=single diffrac targ.,
c        4=non-single diffrac.
c
c (inputs)
c kf1 kf2 : flavor codes of proj. and targ. (Particle data code).
c jq1 jqq1: flavor codes of the valence quark in proj.
c jq2 jqq2: flavor codes of the valence quark in targ.
c ks1 ks2 : status code =1: stable 2:resonance 3:string
c icltyp  : collision type 1:BB 2:MB 3:MM 4:anti-BB
c p1(1)-p1(5)=(p_x,p_y,p_z,E,m) four momentum and the invariant
c            mass of projectile hadron.
c=======================================================================
      include 'jam1.inc'
      include 'jam2.inc'
      real*8 jamrnd2
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      common/hijdat/hidat0(10,10),hidat(10)
      common/hijjet5/nfp(2,10)
      save  /hiparnt/,/hijdat/,/hijjet5/
      dimension p1(15),p2(15),probd(4)
c...Functions: momentum in two-particle cm.
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a) 
      data ibb/0/  ! TEST

      ierror=0
      mste(4)=0

c...Initial light-cone momenta p+,p- of proj. and targ.
      epp=p1(4)+p1(3)
      epm=p1(4)-p1(3)
      etp=p2(4)+p2(3)
      etm=p2(4)-p2(3)
      pr=abs(p1(3))
      jq01=jq1
      jq02=jq2
      kfs1=kf1
      kfs2=kf2

c...Total W+,W- and center-of-mass energy
      wp=epp+etp
      wm=epm+etm
      sw=wp*wm
      srt0=sqrt(sw)

      if(wp.lt.0.0d0 .or. wm.lt.0.0d0) go to 1000

ccc   ihnt2(11)=jp
ccc   ihnt2(12)=jt

      ks1=mod(abs(ks01),10)
      ks2=mod(abs(ks02),10)
        pkc1=0.0d0
        pkc2=0.0d0
        pkc11=0.0d0
        pkc12=0.0d0
        pkc21=0.0d0
        pkc22=0.0d0
        dpkc11=0.0d0
        dpkc12=0.0d0
        dpkc21=0.0d0
        dpkc22=0.0d0

c...In the case of hard scattering
           if(nfp(1,10).eq.1) then
              phi1=pjangl(p1(10),p1(11))
              ppjet=sqrt(p1(10)**2+p1(11)**2)
              pkc1=ppjet
              pkc11=p1(10)
              pkc12=p1(11)
           endif
           if(nfp(2,10).eq.1) then
              phi2=pjangl(p2(10),p2(11))
              ptjet=sqrt(p2(10)**2+p2(11)**2)
              pkc2=ptjet
              pkc21=p2(10)
              pkc22=p2(11)
           endif

      miss=0
10    continue
      ptp02=p1(1)**2+p1(2)**2
      ptt02=p2(1)**2+p2(2)**2

c...In case of excited strings, no diffractive scatt.
      if(mstc(71).eq.0.or.(ks1.eq.3.or.ks2.eq.3)) then
        call jamdiffr(0,srt0,pr,kf1,kf2,ks1,ks2,probd,i_sng)
        i_sng=0

c...Decide whether to have single-diffractive.
      else
        call jamdiffr(1,srt0,pr,kf1,kf2,ks1,ks2,probd,i_sng)
      endif

      icltyp=mste(2)
      iflex=0

      if(kf1.eq.0.or.ks1.eq.3.or.kf1.eq.92) then
        emin1=emjet1+0.1d0
        edmin1=emjet1+0.1d0
        kfd1=92
        kf1=92
      else

c...Flavor exchange according to the SU(6).
        if(icltyp.eq.1) then
        if(kf2.ne.0.and.i_sng.eq.0) then
          call kfcnst(jq2,jqq1,kf1,0.0d0)
          iflex=1
          if(kf1.eq.0) then
            kf1=kfs1
            iflex=0
          endif
        endif
        endif
        call jamdmass(kf1,kfm1,kfd1,emin1,edmin1)
      endif

      if(kf2.eq.0.or.ks2.eq.3.or.kf2.eq.92) then
        emin2=emjet2+0.1d0
        edmin2=emjet2+0.1d0
        kfd2=92
        kf2=92
      else
c...Flavor exchange according to the SU(6).
        if(iflex.eq.1) then        
          call kfcnst(jq1,jqq2,kf2,0.0d0)
          if(kf1.eq.0) then
            kf1=kfs1
            kf2=kfs2
          endif
        endif
        call jamdmass(kf2,kfm2,kfd2,emin2,edmin2)
      endif

      if(kf1.eq.0.or.kf2.eq.0) then
        if(srt0.le.edmin1+edmin2+0.3d0) then
         ierror=999
         return
        endif
      endif

c..test:B*B* ->BBTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      if(ibb.eq.1.and.mste(1).ne.5
     $             .and.(ks1.le.2.and.ks2.le.2)) then
      if(ks1.eq.2.or.ks2.eq.2) then
        if(rn(0).lt.(edmin1+edmin2)**2/sw) then
          if(srt0.lt.emin1+emin2+0.01d0) then
            write(mstc(38),*)'(jamsfex:) srt0<emin1+emin2'
     $                             ,srt0,emin1,emin2,kf1,kf2
            goto 4000
          endif
          prel=pawt(srt0,emin1,emin2)
          swptx=4.0d0*(max(emin1,emin2)**2)
          if(sw.le.swptx) then
            pkcmx=0.0d0
          else if(sw.gt.swptx) then
            pkcmx=sqrt(sw/4.0d0-max(emin1,emin2)**2)
          endif
          widg=0.65d0
          pkc=widg*sqrt(-dlog(1.0d0-rn(0)*(1.0d0-exp( 
     & -pkcmx**2/widg**2))))
          phi0=2.0d0*paru(1)*rn(0)
          p1(1)=pkc*sin(phi0)
          p1(2)=pkc*cos(phi0)
          p2(1)=-p1(1)
          p2(2)=-p1(2)
          p1(3)=sqrt(prel**2-pkc**2)
          p2(3)=-p1(3)
          p1(4)=sqrt(emin1**2+prel**2)
          p2(4)=sqrt(emin2**2+prel**2)
          p1(5)=emin1
          p2(5)=emin2
          kf1=kfm1
          kf2=kfm2
          kf10=kfm1
          kf20=kfm1
          return
        endif
      endif
      endif
c...end test TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

      amp0=emin1
      amt0=emin2
      amx1=edmin1
      amx2=edmin2

      if(srt0.gt.emin1+emin2) then
        pl=pawt(srt0,emin1,emin2)
      else
        write(mstc(38),*)'(jamsfex:)srt<emin1+emin2',srt0,emin1,emin2
        write(mstc(38),*)'kf1 kf2 em1 em2',kf1,kf2,p1(5),p2(5)
        ierror=10
        return
      endif
c     epp=sqrt(emin1**2+pl**2)+pl
c     epm=sqrt(emin1**2+pl**2)-pl
c     etp=sqrt(emin2**2+pl**2)-pl
c     etm=sqrt(emin2**2+pl**2)+pl
c     wp=epp+etp
c     wm=epm+etm
c     sw=wp*wm

      if((i_sng.eq.0.or.i_sng.eq.3).and.mste(1).eq.5) then
c       if(srt0.le.emjet1+emjet2+0.001d0) then
        if(rn(0).le.0.5d0) then
          amp0=amx1
          amx1=emjet1
        else
          amt0=amx2
          amx2=emjet2
        endif
c       else
c        amp0=amx1
c        amx1=emjet1
c        amt0=amx2
c        amx2=emjet2
c       endif
      endif


      if(iflex.eq.1) then
        jq1=jq02
        jq2=jq01
        kf10=kf1
        kf20=kf2
      endif

c...test2:
      if(ibb.eq.2.and.mste(1).ne.5
     $             .and.(ks1.le.2.and.ks2.le.2)) then

        facsp1 =max(1,mod(abs(kfs1),10))*max(1,mod(abs(kfs2),10))
        facsp2 =max(1,mod(abs(kf1),10))*max(1,mod(abs(kf2),10))
        probnn=0.0d0
      if(ks1.eq.2.and.ks2.eq.2) then

c...B*B*->BB
        d1=edmin1**2/sw
        d2=edmin2**2/sw
        bb1=1.0d0+d1-d2
        bb2=1.0d0+d2-d1
        bbb1=bb1**2-4.0d0*d1
        if(bbb1.le.0.0d0) goto 11
        bbb2=bb2**2-4.0d0*d2
        if(bbb2.le.0.0d0) goto 11
        xmin1=(bb1-sqrt(bbb1))/2.0d0
        xmax1=(bb1+sqrt(bbb1))/2.0d0
        xmin2=(bb2-sqrt(bbb2))/2.0d0
        xmax2=(bb2+sqrt(bbb2))/2.0d0
        d1=p1(5)**2/sw
        d2=p2(5)**2/sw
        bb1=1.0d0+d1-d2
        bb2=1.0d0+d2-d1
        bbb1=bb1**2-4.0d0*d1
        if(bbb1.le.0.0d0) goto 11
        bbb2=bb2**2-4.0d0*d2
        if(bbb2.le.0.0d0) goto 11
        x1=(bb1-sqrt(bbb1))/2.0d0
        x2=(bb2-sqrt(bbb2))/2.0d0
        probnn=log(xmax1/xmin1)*log(xmax2/xmin2)/(x1*x2)
     $                               *(probd(1)+probd(4))
       
c...B1*B2->B1B2
      else if(ks1.eq.2.and.ks2.eq.1) then
        d1=edmin1**2/sw
        d2=emin2**2/sw
        bb1=1.0d0+d1-d2
        bbb=bb1**2-4.0d0*d1
        if(bbb.le.0.0d0) goto 11
        xmin1=(bb1-sqrt(bbb))/2.0d0
        xmax1=(bb1+sqrt(bbb))/2.0d0

        d1=p1(5)**2/sw
        d2=p2(5)**2/sw
        bb1=1.0d0+d1-d2
        bb2=1.0d0+d2-d1
        bbb=bb1**2-4.0d0*d1
        if(bbb.le.0.0d0) goto 11
        x1=(bb1-sqrt(bb1**2-4.0d0*d1))/2.0d0

        probnn=log(xmax1/xmin1)/x1*probd(2)

c...B1B2*->B1B2
      else if(ks1.eq.1.and.ks2.eq.2) then
        d1=emin1**2/sw
        d2=edmin2**2/sw
        bb2=1.0d0+d2-d1
        bbb=bb2**2-4.0d0*d1
        if(bbb.le.0.0d0) goto 11
        xmin2=(bb2-sqrt(bbb))/2.0d0
        xmax2=(bb2+sqrt(bbb))/2.0d0
        d1=p1(5)**2/sw
        d2=p2(5)**2/sw
        bb1=1.0d0+d1-d2
        bb2=1.0d0+d2-d1
        bbb=bb2**2-4.0d0*d2
        if(bbb.le.0.0d0) goto 11
        x2=(bb2-sqrt(bbb))/2.0d0
        probnn=log(xmax2/xmin2)/x2*probd(3)
      endif

      probnn=facsp2/facsp1*(pl/pr)**2*probnn
      if(kfs1.ne.kfs2.and.kfm1.eq.kfm2) probnn=0.5d0*probnn

      if(rn(0).le.probnn/(1+probnn)) then
          if(sw.le.swptx) then
            pkcmx=0.0d0
          else if(sw.gt.swptx) then
            pkcmx=sqrt(sw/4.0d0-max(emin1,emin2)**2)
          endif
          widg=0.65d0
          pkc=widg*sqrt(-dlog(1.0d0-rn(0)*(1.0d0-exp( 
     & -pkcmx**2/widg**2))))

          phi0=2.0d0*paru(1)*rn(0)
          p1(1)=pkc*sin(phi0)
          p1(2)=pkc*cos(phi0)
          p2(1)=-p1(1)
          p2(2)=-p1(2)
          p1(3)=sqrt(pl**2-pkc**2)
          p2(3)=-p1(3)
          p1(4)=sqrt(emin1**2+pl**2)
          p2(4)=sqrt(emin2**2+pl**2)
          p1(5)=emin1
          p2(5)=emin2
          kf1=kfm1
          kf2=kfm2
          kf10=kfm1
          kf20=kfm2
          return
        endif

      endif
c...end test2


11    continue
c...Scatter only if sw>snn
      if(sw.lt.(amp0+amt0)*2+0.001d0) go to 4000

c-----give some PT kick to the two exited strings-----

c...Find maximun PT kick.
20    swptn=4.0d0*(max(amp0,amt0)**2+max(ptp02,ptt02))
      swptx=4.0d0*(max(amx1,amx2)**2+max(ptp02,ptt02))

      if(sw.le.swptn) then
        pkcmx=0.0d0
      else if(sw.le.swptx) then
        pkcmx=sqrt(sw/4.0d0-max(amp0,amt0)**2)-sqrt(max(ptp02,ptt02))
      else
        pkcmx=sqrt(sw/4.0d0-max(amx1,amx2)**2)-sqrt(max(ptp02,ptt02))
      endif

c...If the valence quarks had a hard-collision
c...the pt kick is the pt from hard-collision.
        if(nfp(1,10).eq.1.or.nfp(2,10).eq.1) then
                if(pkc1.gt.pkcmx) then
                        pkc1=pkcmx
                        pkc11=pkc1*cos(phi1)
                        pkc12=pkc1*sin(phi1)
                        dpkc11=-(p1(10)-pkc11)/2.0d0
                        dpkc12=-(p1(11)-pkc12)/2.0d0
                endif
                if(pkc2.gt.pkcmx) then
                        pkc2=pkcmx
                        pkc21=pkc2*cos(phi2)
                        pkc22=pkc2*sin(phi2)
                        dpkc21=-(p2(10)-pkc21)/2.0d0
                        dpkc22=-(p2(11)-pkc22)/2.0d0
                endif
                dpkc1=dpkc11+dpkc21
                dpkc2=dpkc12+dpkc22
                nfp(1,10)=-nfp(1,10)
                nfp(2,10)=-nfp(2,10)
                go to 40
        endif


c---------Select PT kick--------
c...Option for Gaussian pt kick
      if(mstc(75).eq.0) then
        pkc=parc(66)*sqrt(-dlog(1.0d0-rn(0)
     &                  *(1.0d0-exp(-pkcmx**2/parc(66)**2))))
        widg=0.65d0
        if(i_sng.eq.1.or.i_sng.eq.2) pkc=widg*sqrt(
     &          -dlog(1.0d0-rn(0)*(1.0d0-exp(-pkcmx**2/widg**2))))
      else

c.....HIJING parametrization.
        pkc=sqrt(hirnd2(3,0.0d0,pkcmx**2))
c       if(pkcmx.gt.6.0d0) then
c          print *,'pkcmx=',pkcmx
c       endif
c       pkc=jamrnd2(3,0.0d0,min(6.0d0,pkcmx))
        if(pkc.gt.parc(68)) 
     &     pkc=parc(66)*sqrt(-dlog(exp(-parc(68)**2/parc(66)**2)
     &         -rn(0)*(exp(-parc(68)**2/parc(66)**2)-
     &         exp(-pkcmx**2/parc(66)**2))))

        widg=0.65d0
        if(i_sng.eq.1.or.i_sng.eq.2) pkc=widg*sqrt(
     &          -dlog(1.0d0-rn(0)*(1.0d0-exp(-pkcmx**2/widg**2))))

      endif

30    phi0=2.0d0*paru(1)*rn(0)
      pkc11=pkc*sin(phi0)
      pkc12=pkc*cos(phi0)
      pkc21=-pkc11
      pkc22=-pkc12
      dpkc1=0.0d0
      dpkc2=0.0d0

40    pp11=p1(1)+pkc11-dpkc1
      pp12=p1(2)+pkc12-dpkc2
      pt11=p2(1)+pkc21-dpkc1
      pt12=p2(2)+pkc22-dpkc2
      ptp2=pp11**2+pp12**2
      ptt2=pt11**2+pt12**2

      ampn=sqrt(amp0**2+ptp2)
      amtn=sqrt(amt0**2+ptt2)
      snn=(ampn+amtn)**2+0.001d0

      wp=epp+etp
      wm=epm+etm
      sw=wp*wm

      if(sw.lt.snn) then
        miss=miss+1
        if(miss.le.100) then
          pkc=0.0d0
          go to 30
        endif
        if(mstc(8).ne.0) 
     &  write(mstc(38),*) 'error occured in pt kick section of jamsfex'
        go to 4000
      endif

      ampx=sqrt(amx1**2+ptp2)
      amtx=sqrt(amx2**2+ptt2)

      dpn=ampn**2/sw
      dtn=amtn**2/sw
      dpx=ampx**2/sw
      dtx=amtx**2/sw

c...C.M. energy if proj=N,targ=N*
      spntx=(ampn+amtx)**2

c...C.M. energy if proj=N*,targ=N
      spxtn=(ampx+amtn)**2

c...C.M. energy if proj=delta, targ=delta
      sxx=(ampx+amtx)**2

      ipath=0 
 45   continue
      if(sw.gt.sxx+0.001d0) then

c....Non diffractive
        if(i_sng.eq.0) then
 46       d1=dpx
          d2=dtx
          kf1=92
          kf2=92
          ipath=1
          go to 400

c...Have single diffractive collision
        else if(i_sng.eq.1) then
          d1=dpn
          d2=dtx
          kf1=kfm1
          kf2=92
          ipath=2
          go to 220
        else if(i_sng.eq.2) then
          d1=dpx
          d2=dtn
          kf1=92
          kf2=kfm2
          ipath=3
          go to 240

c...Double diffractive
        else if(i_sng.eq.3) then
          d1=dpx
          d2=dtx
          kf1=92
          kf2=92
          ipath=4
          go to 400
        endif

      else if(sw.gt.max(spntx,spxtn)+0.001d0) then

           if(rn(0).gt.0.5d0) then
              d1=dpn
              d2=dtx
              kf1=kfm1
              kf2=92
              ipath=5
              go to 220
           else
              d1=dpx
              d2=dtn
              kf1=92
              kf2=kfm2
              ipath=6
              go to 240
           endif

      else if(sw.gt.min(spntx,spxtn)+0.001d0) then

           if(spntx.le.spxtn) then
              d1=dpn
              d2=dtx
              kf1=kfm1
              kf2=92
              ipath=7
              go to 220
           else
              d1=dpx
              d2=dtn
              kf1=92
              kf2=kfm2
              ipath=8
              go to 240
           endif

      endif

      ierror=999
      if(mstc(8).eq.0.or.mstc(13).eq.0) return
      if(mstc(13).eq.1.and.(mstd(25).gt.mstc(14))) return
      mstd(25)=mstd(25)+1
      write(mstc(38),*)' error in jamsfex: there is no path to here'
      write(mstc(38),*)'srt spntx spxtn sxx',
     $       srt0,sqrt(spntx),sqrt(spxtn),sqrt(sxx)
      write(mstc(38),*)'kf1 em1',kf1,p1(5),amp0,amx1
      write(mstc(38),*)'kf2 em2',kf2,p2(5),amt0,amx2
      return

c===============  elastic scattering ===============
c       this is like elastic, both proj and targ mass
c       must be fixed
c===================================================
100   continue

      bb1=1.0d0+d1-d2
      bb2=1.0d0+d2-d1
      if(bb1**2.lt.4.0d0*d1 .or. bb2**2.lt.4.0d0*d2) then
        miss=miss+1
        if(miss.gt.100.or.pkc.eq.0.0d0) go to 3000
        pkc=pkc*0.5d0
        go to 30
      endif
      if(rn(0).lt.0.5d0) then
        x1=(bb1-sqrt(bb1**2-4.0d0*d1))/2.0d0
        x2=(bb2-sqrt(bb2**2-4.0d0*d2))/2.0d0
      else
        x1=(bb1+sqrt(bb1**2-4.0d0*d1))/2.0d0
        x2=(bb2+sqrt(bb2**2-4.0d0*d2))/2.0d0
      endif
      mste(4)=1
      go to 600

c========== Single diffractive =======================
c either proj or targ's mass is fixed
c=====================================================
c...Fix proj mass
220   continue

      bb2=1.0d0+d2-d1
      if(bb2**2.lt.4.0d0*d2) then
        miss=miss+1
        if(miss.gt.100.or.pkc.eq.0.0d0) go to 3000
        pkc=pkc*0.5d0
        go to 30
      endif

      xmin=(bb2-sqrt(bb2**2-4.0d0*d2))/2.0d0
      xmax=(bb2+sqrt(bb2**2-4.0d0*d2))/2.0d0
      if(xmin.le.1.d-7) xmin=0.00001d0

      miss4=0
222   continue

      x2=xsamp3(xmin,xmax,srt0)
      x1=d1/(1.0d0-x2)
      if(x2*(1.0d0-x1).lt.(d2+1.d-4/sw)) then
        miss4=miss4+1
        if(miss4.le.1000) go to 222
        go to 5000
      endif

      mste(4)=2
      go to 600

c...Fix targ mass
240   continue

      bb1=1.0d0+d1-d2
      if(bb1**2.lt.4.0d0*d1) then
        miss=miss+1
        if(miss.gt.100.or.pkc.eq.0.0d0) go to 3000
        pkc=pkc*0.5d0
        go to 30
      endif

      xmin=(bb1-sqrt(bb1**2-4.0d0*d1))/2.0d0
      xmax=(bb1+sqrt(bb1**2-4.0d0*d1))/2.0d0

      miss4=0
242   continue

      x1=xsamp3(xmin,xmax,srt0)
      x2=d2/(1.0d0-x1)
      if(x1*(1.0d0-x2).lt.(d1+1.d-4/sw)) then
        miss4=miss4+1
        if(miss4.le.1000) go to 242
        go to 5000
      endif

      mste(4)=3
      go to 600

c============= non-single diffractive ====================
c       both proj and targ may not be fixed in mass 
c=========================================================
400   continue

      bb1=1.0d0+d1-d2
      bb2=1.0d0+d2-d1
      if(bb1**2.lt.4.0d0*d1 .or. bb2**2.lt.4.0d0*d2) then
        miss=miss+1
        if(miss.gt.100.or.pkc.eq.0.0d0) go to 3000
        pkc=pkc*0.5d0
        go to 30
      endif
      xmin1=(bb1-sqrt(bb1**2-4.0d0*d1))/2.0d0
      xmax1=(bb1+sqrt(bb1**2-4.0d0*d1))/2.0d0
      xmin2=(bb2-sqrt(bb2**2-4.0d0*d2))/2.0d0
      xmax2=(bb2+sqrt(bb2**2-4.0d0*d2))/2.0d0

      miss4=0 
410   continue

      if(i_sng.eq.0) then
      if(ks1.eq.1.and.ks2.eq.1) then
c....qq-q~q~ or q-q~
        if(abs(jq1*jqq1).ge.100.and.abs(jq1*jqq1).le.1000000) then
           x1=xsamp1(xmin1,xmax1,srt0)
        else
           x1=xsamp2(xmin1,xmax1,srt0)
        endif
        if(abs(jq2*jqq2).ge.100.and.abs(jq2*jqq2).le.1000000) then
          x2=xsamp1(xmin2,xmax2,srt0)
        else
          x2=xsamp2(xmin2,xmax2,srt0)
        endif
      else
        x1=xsamp3(xmin1,xmax1,srt0)
        x2=xsamp3(xmin2,xmax2,srt0)
      endif

      else
        x1=xsamp3(xmin1,xmax1,srt0)
        x2=xsamp3(xmin2,xmax2,srt0)
      endif

      if(abs(jq1*jqq1).gt.1000000) x1=1.0d0-x1
      xxp=x1*(1.0d0-x2)
      xxt=x2*(1.0d0-x1)
      if(xxp.lt.(d1+1.d-4/sw) .or. xxt.lt.(d2+1.d-4/sw)) then
        miss4=miss4+1
        if(miss4.le.1000) go to 410
        go to 5000
      endif
      mste(4)=4

c==========================================================
c...Now set masses and momentum of the outgoing particles
c==========================================================
600   continue
      if(x1*(1.0d0-x2).lt.(ampn**2-1.d-4)/sw.or.
     &                  x2*(1.0d0-x1).lt.(amtn**2-1.d-4)/sw) then
        miss=miss+1
        if(miss.gt.100.or.pkc.eq.0.0d0) go to 2000
        pkc=0.0d0
        go to 30
      endif

      epp1=(1.0d0-x2)*wp
      epm1=x1*wm
      etp1=x2*wp
      etm1=(1.0d0-x1)*wm

      if((epm1.gt.epp1 .or. etp1.gt.etm1)
     $   .or.(epp1/(epm1+0.01d0).lt.etp1/(etm1+0.01d0))) then
        goto 30
      endif

c...Set proj. energy, momentum and flavor.
      p1(3)=(epp1-epm1)/2.0d0
      p1(4)=(epp1+epm1)/2.0d0
      if(epp1*epm1-ptp2.lt.0.0d0) go to 6000
      p1(5)=sqrt(epp1*epm1-ptp2)

c...Set targ. energy, momentum and flavor.
      p2(3)=(etp1-etm1)/2.0d0
      p2(4)=(etp1+etm1)/2.0d0
      if(etp1*etm1-ptt2.lt.0.0d0) go to 6000
      p2(5)=sqrt(etp1*etm1-ptt2)

C...Recoil pt from hard-inter is shared by two end-partons
c...so that pt=p1+p2
      kickdip=1
      kickdit=1
      if(abs(jq1*jqq1).gt.1000000.or.abs(jq1*jqq1).lt.100) kickdip=0
      if(abs(jq2*jqq2).gt.1000000.or.abs(jq2*jqq2).lt.100) kickdit=0

c...pt kick for proj.
      p1(1)=pp11-pkc11
      p1(2)=pp12-pkc12
      if((kickdip.eq.0.and.rn(0).lt.0.5d0)
     &     .or.(kickdip.ne.0.and.rn(0)
     &     .lt.0.5d0/(1.0d0+(pkc11**2+pkc12**2)/hipr1(22)**2))) then
        p1(6)=(p1(1)-p1(6)-p1(8)-dpkc1)/2.0d0+p1(6)
        p1(7)=(p1(2)-p1(7)-p1(9)-dpkc2)/2.0d0+p1(7)
        p1(8)=(p1(1)-p1(6)-p1(8)-dpkc1)/2.0d0+p1(8)+pkc11
        p1(9)=(p1(2)-p1(7)-p1(9)-dpkc2)/2.0d0+p1(9)+pkc12
      else
        p1(8)=(p1(1)-p1(6)-p1(8)-dpkc1)/2.0d0+p1(8)
        p1(9)=(p1(2)-p1(7)-p1(9)-dpkc2)/2.0d0+p1(9)
        p1(6)=(p1(1)-p1(6)-p1(8)-dpkc1)/2.0d0+p1(6)+pkc11
        p1(7)=(p1(2)-p1(7)-p1(9)-dpkc2)/2.0d0+p1(7)+pkc12
      endif
      p1(1)=p1(6)+p1(8)
      p1(2)=p1(7)+p1(9)


      p2(1)=pt11-pkc21
      p2(2)=pt12-pkc22
c...pt kick for targ.
      if((kickdit.eq.0.and.rn(0).lt.0.5d0)
     &     .or.(kickdit.ne.0.and.rn(0)
     &     .lt.0.5d0/(1.0d0+(pkc21**2+pkc22**2)/hipr1(22)**2))) then
        p2(6)=(p2(1)-p2(6)-p2(8)-dpkc1)/2.0d0+p2(6)
        p2(7)=(p2(2)-p2(7)-p2(9)-dpkc2)/2.0d0+p2(7)
        p2(8)=(p2(1)-p2(6)-p2(8)-dpkc1)/2.0d0+p2(8)+pkc21
        p2(9)=(p2(2)-p2(7)-p2(9)-dpkc2)/2.0d0+p2(9)+pkc22
      else
        p2(8)=(p2(1)-p2(6)-p2(8)-dpkc1)/2.0d0+p2(8)
        p2(9)=(p2(2)-p2(7)-p2(9)-dpkc2)/2.0d0+p2(9)
        p2(6)=(p2(1)-p2(6)-p2(8)-dpkc1)/2.0d0+p2(6)+pkc21
        p2(7)=(p2(2)-p2(7)-p2(9)-dpkc2)/2.0d0+p2(7)+pkc22
      endif
      p2(1)=p2(6)+p2(8)
      p2(2)=p2(7)+p2(9)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      d1=sqrt(p1(5)**2+p1(1)**2+p1(2)**2+p1(3)**2)
      d2=sqrt(p2(5)**2+p2(1)**2+p2(2)**2+p2(3)**2)
      if((abs(d1-p1(4)).gt.0.001d0).or.(abs(d1-p1(4)).gt.0.001d0)) then
        ih=mstc(38)
        write(ih,*)'(jamsfex:) ipath i_sng kf1 kf2',ipath,i_sng,kf1,kf2
        write(ih,*)'D1 p4',d1,p1(4)
        write(ih,*)'D2 p4',d2,p2(4)
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      return
c=======================================================================
c...Print errors
c=======================================================================

1000  ierror=999
      if(mstc(8).eq.0.or.mstc(13).eq.0) return
      if(mstc(13).eq.1.and.(mstd(25).gt.mstc(14))) return
      mstd(25)=mstd(25)+1
      ih=mstc(38)
      write(ih,*)'    fatal jamsfex start error,abandon this event'
      write(ih,*)'    proj E+,E-,W+',kf1,epp,epm,wp
      write(ih,*)'    targ E+,E-,W-',kf2,etp,etm,wm
      write(ih,*)'    W+*W-, (apn+atn)^2',sw,snn
      return

2000  ierror=2000
      if(mstc(8).eq.0.or.mstc(13).eq.0) return
      if(mstc(13).eq.1.and.(mstd(25).gt.mstc(14))) return
      mstd(25)=mstd(25)+1
      ih=mstc(38)
      write(ih,*)'    (2)energy partition fail,',kf1,kf2
      write(ih,*)'    jamsfex not performed, but continue'
      write(ih,*)'    mp1,mpn',x1*(1.0d0-x2)*sw,ampn**2
      write(ih,*)'    mt2,mtn',x2*(1.0d0-x1)*sw,amtn**2
      return

3000  ierror=3000
      if(mstc(8).eq.0.or.mstc(13).eq.0) return
      if(mstc(13).eq.1.and.(mstd(25).gt.mstc(14))) return
      mstd(25)=mstd(25)+1
      ih=mstc(38)
      write(ih,*)'    (3)something is wrong with the pt kick, '
      write(ih,*)'    jamsfex not performed, but continue'
      write(ih,*)'    d1=',d1,' d2=',d2,' sw=',sw
      return

4000  ierror=4000
      if(mstc(8).eq.0.or.mstc(13).eq.0) return
      if(mstc(13).eq.1.and.(mstd(25).gt.mstc(14))) return
      mstd(25)=mstd(25)+1
      ih=mstc(38)
      write(ih,*)'***(4)unable to choose process, but not harmful'
      write(ih,*)' jamsfex not performed, but continue'
      write(ih,*)' ptp=',ptp2,' ptt=',ptt2
      write(ih,*)' sw=',sw,' snn=',snn
      write(ih,*)' amcut=',amx1,amx2
      write(ih,*)' kf1 em1',kf1,p1(5)
      write(ih,*)' kf2 em2',kf2,p2(5)
      return

5000  ierror=5000
      if(mstc(8).eq.0.or.mstc(13).eq.0) return
      if(mstc(13).eq.1.and.(mstd(25).gt.mstc(14))) return
      mstd(25)=mstd(25)+1
      ih=mstc(38)
      write(ih,*)'    energy partition failed(5),for limited try'
      write(ih,*)'    jamsfex not performed, but continue'
      write(ih,*)'    kf1=',kf1,' kf2=',kf2
      write(ih,*)'    d1',d1,' x1(1-x2)',x1*(1.0d0-x2)
      write(ih,*)'    d2',d2,' x2(1-x1)',x2*(1.0d0-x1)
      return

6000  pkc=0.0d0
      miss=miss+1
      if(miss.lt.100) go to 30
      if(mstc(8).eq.0.or.mstc(13).eq.0) return
      if(mstc(13).eq.1.and.(mstd(25).gt.mstc(14))) return
      mstd(25)=mstd(25)+1
      ih=mstc(38)
      ierror=1
      if(mstc(8).eq.0) return
      write(ih,*)'error occured, jamsfex not performed,abort this event'
      write(ih,*)'mtp,ptp2',epp*epm,ptp2,'  mtt,ptt2',etp*etm,ptt2 

      end

c***********************************************************************

      subroutine jamdiffr(msel,srt,pr,kf1,kf2,ks1,ks2,probd,idffra)

c...Purpose: to calculate diffractive cross sections.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      dimension sigt(0:6,0:6,0:5),probd(4)

c...sig = sig0 + sig_el + 2sig_sd + sig_dd
      idffra=0
      do i=1,4
       probd(i)=0.0d0
      end do
c...Do not allow excited strings to have single-diffr.
      if(ks1.eq.3.or.ks2.eq.3) return

c...Calculate cross sections.
      call jamxtot(kf1,kf2,srt,pr,sigt)
      sigd1=sigt(0,0,2)
      sigd2=sigt(0,0,3)
      sigdd=sigt(0,0,4)
      signo=sigt(0,0,5)
      sigint=sigd1+sigd2+sigdd+signo
      if(sigint.le.0.0d0) then
        return
      endif
      probd(1)=signo/sigint
      probd(2)=sigd1/sigint
      probd(3)=sigd2/sigint
      probd(4)=sigdd/sigint
      if(msel.eq.0) return

      if(mste(1).eq.4) then
        rans=rn(0)*sigint
c...Single diffractive dissociation.  a + b ==>  X + b
      if(rans.lt.sigd1) then
        idffra=2
      else if(rans.lt.sigd1+sigd2) then ! a + b ==> a + X
        idffra=1
c...Double diffractive dissociation.  a + b ==> X1 + X2
      else if(rans.lt.sigd1+sigd2+sigdd) then
        idffra=3
c...Non diffractive
      else
        idffra=0
      endif

c....In this case,contribution of resonance should be deleted due to
c....the explicit parametrization.
      else if(mste(1).eq.5) then

c....2-string not arrowed.
       if(srt.le.parc(51)+0.95d0) then
         idffra=2
         if(rn(0).le.0.5d0) idffra=1
         return
c....Somoth transition to high energy part of diffractive cross section.
       else if(srt.le.10.0d0) then
         sigd1=10.5372d0*(srt-3.0d0)**1.06245d0/(15.2027d0 
     &                               +(srt-9.43717d0)**2)/2.d0
         sigd2=sigd1
       endif

       pare(3)=pare(3)-sigd1
       if(pare(3).lt.0.0d0) then
         idffra=2
         return
       endif
       pare(3)=pare(3)-sigd2
       if(pare(3).lt.0.0d0) then
         idffra=1
         return
       endif
       pare(3)=pare(3)-sigdd
       if(pare(3).lt.0.0d0) then
         idffra=3
         return
       endif
       idffra=0
      else
        call jamerrm(30,0,'(jamdiffr:)???mste(1)')
      endif

      end
 
c***********************************************************************

      subroutine jamjetm1(is,ip,prel,jq1,jqq1,bex,bey,bez,gam,phi,the)

c...Purpose: to save parton properties for jet system.
      include 'jam1.inc'
      include 'jam2.inc'
      dimension prel(15)
      Dimension dps(5)

      qmass1=pjmass(jq1)
      qmass2=pjmass(jqq1)
      emt=prel(5)**2+prel(1)**2+prel(2)**2
      emt1=qmass1**2+prel(6)**2+prel(7)**2
      emt2=qmass2**2+prel(8)**2+prel(9)**2
 
      if(emt.gt.emt1+emt2) then
        pzcm=sqrt(abs(emt**2+emt1**2+emt2**2-2.0d0*emt*emt1
     &       -2.0d0*emt*emt2-2.0d0*emt1*emt2))/2.0d0/sqrt(emt)
      else
        write(check(1),8000)emt,emt1,emt2
        write(check(2),8010)'prel',(prel(l),l=1,5)
        call jamerrm(30,2,'(jamjetm1:)pzcm<0')
 8000   format('emt emt1 emt2=',3(g12.4,1x))
 8010   format('prel=',5(g12.4,1x))
      endif
 
      q11x=prel(6)
      q11y=prel(7)
      q11z=-pzcm*is
      q11e=sqrt(emt1+pzcm**2)
 
      q12x=prel(8)
      q12y=prel(9)
      q12z=pzcm*is
      q12e=sqrt(emt2+pzcm**2)
      emt=sqrt(emt)

c...Transform to two-body cm
      gg1=prel(4)/emt
      call jamrobo(0.0d0,0.0d0,0.0d0,0.0d0,prel(3)/prel(4),gg1,
     $              q11x,q11y,q11z,q11e)
      call jamrobo(0.0d0,0.0d0,0.0d0,0.0d0,prel(3)/prel(4),gg1,
     $              q12x,q12y,q12z,q12e)
 
c...Transform to comp. frame
      call jamrobo(the,phi,bex,bey,bez,gam,q11x,q11y,q11z,q11e)
      call jamrobo(the,phi,bex,bey,bez,gam,q12x,q12y,q12z,q12e)

      kq(1,ip)=jq1
      vq(1,ip)=q11x
      vq(2,ip)=q11y
      vq(3,ip)=q11z
      vq(4,ip)=sqrt(qmass1**2+q11x**2+q11y**2+q11z**2)
      vq(5,ip)=qmass1

      kq(2,ip)=jqq1
      vq(6,ip)=q12x
      vq(7,ip)=q12y
      vq(8,ip)=q12z
      vq(9,ip)=sqrt(qmass2**2+q12x**2+q12y**2+q12z**2)
      vq(10,ip)=qmass2

cRRRRRRRRRRRRRRRRRRRRRRRRRR
      do  j=1,5
       dps(j)=vq(j,ip)+vq(j+5,ip)
      end do
      emsq=dps(4)**2-dps(1)**2-dps(2)**2-dps(3)**2
      if(emsq.lt.(0.9d0*parj(32)+dps(5))**2) then
        Write(3,*)'(jamjetm1:)???? ip',ip
      endif
cRRRRRRRRRRRRRRRRRRRRRRRRRR

      end

c***********************************************************************

      subroutine jamjetm2(ip)

c...Purpose: to save parton properties of jet system after annihilation.
      include 'jam1.inc'
      include 'jam2.inc'
c...ip: line number of particle.
 
      kflv=k(2,ip)
      if(kflv.eq.0) then
        call jamerrm(30,0,'(jamjetm2:) kflv=0')
      endif

      call attflv(kflv,ifla,iflb)
      qmass1=pjmass(ifla)
      qmass2=pjmass(iflb)

      emj=p(5,ip)
      pcmq=(emj*emj-(qmass1+qmass2)**2)*(emj*emj-(qmass1-qmass2)**2)
      if(pcmq.gt.0.0d0) then
        pcmq=sqrt(pcmq)/(2.d0*emj)
      else
        call jamerrm(30,0,'(:) pcmq<0')
      endif
 
      q11x=0.0d0
      q11y=0.0d0
      q11z=-pcmq
      q11e=sqrt(qmass1**2+pcmq**2)
 
      q12x=0.0d0
      q12y=0.0d0
      q12z=pcmq
      q12e=sqrt(qmass2**2+pcmq**2)

c...Transform to comp. frame.
      the=pjangl(p(3,ip),sqrt(p(1,ip)**2+p(2,ip)**2)) 
      phi=pjangl(p(1,ip),p(2,ip))
      bex=p(1,ip)/p(4,ip)
      bey=p(2,ip)/p(4,ip)
      bez=p(3,ip)/p(4,ip)
      gg1=p(4,ip)/emj
      call jamrobo(the,phi,bex,bey,bez,gg1,q11x,q11y,q11z,q11e)
      call jamrobo(the,phi,bex,bey,bez,gg1,q12x,q12y,q12z,q12e)

      kq(1,ip)=ifla
      vq(1,ip)=q11x
      vq(2,ip)=q11y
      vq(3,ip)=q11z
      vq(4,ip)=sqrt(qmass1**2+q11x**2+q11y**2+q11z**2)
      vq(5,ip)=qmass1

      kq(2,ip)=iflb
      vq(6,ip)=q12x
      vq(7,ip)=q12y
      vq(8,ip)=q12z
      vq(9,ip)=sqrt(qmass2**2+q12x**2+q12y**2+q12z**2)
      vq(10,ip)=qmass2

      if(mstc(8).ge.1) then
      if( (abs( p(1,ip)-(vq(1,ip)+vq(6,ip)) ).gt.0.01d0)
     $ .or. (abs( p(2,ip)-(vq(2,ip)+vq(7,ip)) ).gt.0.01d0)
     $ .or. (abs( p(3,ip)-(vq(3,ip)+vq(8,ip)) ).gt.0.01d0)
     $ .or. (abs( p(4,ip)-(vq(4,ip)+vq(9,ip)) ).gt.0.01d0) ) then
         ih=mstc(38)
         kc1=jamcomp(kflv)
         write(ih,*)'after jamjetm2??? i1 k1 ',k(1,ip)
         write(ih,*)'kf',ip,kflv,p(5,ip)
     $        ,' ',chaf(kc1,(3-isign(1,kflv))/2),the,phi
         write(ih,*)'p1 vq1+vq6',p(1,ip),vq(1,ip)+vq(6,ip)
         write(ih,*)'p2 vq2+vq7',p(2,ip),vq(2,ip)+vq(7,ip)
         write(ih,*)'p3 vq3+vq8',p(3,ip),vq(3,ip)+vq(8,ip)
         write(ih,*)'p4 vq4+vq9',p(4,ip),vq(4,ip)+vq(9,ip)
      endif
      endif

      end

c***********************************************************************

      subroutine jamdmass(kf1,kfm,kfd,emin,emdn)

c...To give minimum and min. excited mass for the particle.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      parameter(dltm=0.3d0)

      kc1=jamcomp(kf1)
      if(kchg(kc1,3).eq.0) then
        kfsign=1
      else
        kfsign=isign(1,kf1)
      endif
      kfa=abs(kf1)
      id1=kchg(kc1,5)
      iz1a=kchg(kc1,1)/3
      ibar=kchg(kc1,6)

c...Baryons.
      if(ibar.eq.3) then
      if(id1.eq.id_nucl) then
         kfm=kf1
         kfd=(10000+abs(kf1))*kfsign
         emin=pmas(kc1,1)
         emdn=parc(64)
      else if(id1.eq.id_nucls) then
         if(iz1a.eq.1) then
           emin=parc(25)
           kfm=2212*kfsign
           kfd=12212*kfsign
         else if(iz1a.eq.0) then
           emin=parc(24)
           kfm=2112*kfsign
           kfd=12112*kfsign
         endif
         emdn=parc(64)
      else if(id1.eq.id_delt.or.id1.eq.id_delts) then
         if(iz1a.eq.1) then
           emin=parc(25)
           kfm=2212*kfsign
           kfd=12212*kfsign
         else if(iz1a.eq.0) then
           emin=parc(24)
           kfm=2112*kfsign
           kfd=12112*kfsign
         else if(iz1a.eq.-1) then
           kfm=1114*kfsign
           kfd=kfm
           emin=parc(64)
         else if(iz1a.eq.2) then
           kfm=2224*kfsign
           kfd=kfm
           emin=parc(64)
         endif
         emdn=parc(64)
      else if(id1.eq.id_lamb.or.id1.eq.id_lambs) then
        kfm=3122*kfsign
        kfd=13122*kfsign
        emin=pmas(jamcomp(3122),1)
        emdn=pmas(mstc(26),1)-pmas(mstc(26),3)
      else if(id1.eq.id_sigm.or.id1.eq.id_sigms) then
        if(iz1a.eq.-1) then
            kcp=jamcomp(3114)
            emdn=pmas(kcp,1)-pmas(kcp,3)
            kcp=jamcomp(3112)
            emin=pmas(kcp,1)
            kfd=3114*kfsign
            kfm=3112*kfsign
        else if(iz1a.eq. 0) then
            kcp=jamcomp(3214)
            emdn=pmas(kcp,1)-pmas(kcp,3)
            kcp=jamcomp(3212)
            emin=pmas(kcp,1)
            kfd=3214*kfsign
            kfm=3212*kfsign
        else if(iz1a.eq. 1) then
            kcp=jamcomp(3224)
            emdn=pmas(kcp,1)-pmas(kcp,3)
            kcp=jamcomp(3222)
            emin=pmas(kcp,1)
            kfd=3224*kfsign
            kfm=3222*kfsign
        endif
      else if(id1.eq.id_xi.or.id1.eq.id_xis) then
        if(iz1a.eq.-1) then
           kfd=3314*kfsign
           kfm=3312*kfsign
           kcp=jamcomp(3314)
           emdn=pmas(kcp,1)-pmas(kcp,3)
           kcp=jamcomp(3312)
           emin=pmas(kcp,1)
        else if(iz1a.eq. 0) then
           kcp=jamcomp(3324)
           emdn=pmas(kcp,1)-pmas(kcp,3)
           emin=pmas(jamcomp(3322),1)
           kfd=3324*kfsign
           kfm=3322*kfsign
        endif
      else if(id1.eq.id_omega) then
        kfd=3334*kfsign
        kfm=3334*kfsign
        emdn=1.8d0
        emin=pmas(jamcomp(3334),1)
      else
        kcp=jamcomp(kf1)
        emdn=pmas(kcp,1)-pmas(kcp,3)+dltm
        emin=pmas(kcp,1)
        kfm=kf1
        kfd=kf1
      endif

c...Mesons.
      else if(ibar.eq.0) then
        kflb=mod(kfa/100,10) 
        kflc=mod(kfa/10,10) 
        if(kflb.le.2.and.kflc.le.2) then
          emdn=0.37d0
        else if(kflb.eq.3.and.kflc.eq.3.and.kflb+kflc.lt.6) then
          emdn=0.8d0
        else if(kflb.eq.3.and.kflc.eq.3) then
          emdn=1.2d0
        else
          pmb=parf(100+kflb) 
          pmc=parf(100+kflc) 
          emdn=pmb+pmc+1.0d0
        endif
        kfd=kf1
        if(id1.eq.id_pi.or.id1.eq.id_str) kfd=(kfa+2)*kfsign
        kfmes=10*kflb+kflc
        kfm=(kfmes*10+1)*kfsign
        emin=pmas(jamcomp(kfm),1)
      else
        write(check(1),'(''ibar kf1'',i4,1x,i9)')ibar,kf1
        call jamerrm(30,1,'(jamdmass:) Unrecognized ibar code')
      endif

      end

c***********************************************************************
c
c                string fragmentation part
c
c***********************************************************************

      subroutine jamjdec(ip,indd,npar,icon)

c...Purpose: to fragment jet system by Lund model.
c...ip   : line number of the decaying particle.
c...indd : contains line number of produced particles.
c...npar : number of particle produced in this time. 
c...icon : condition code;
c...       =0: Ok.

      include 'jam1.inc'
      include 'jam2.inc'
c...PYTHIA commonblocks.
c... Note: meaning of vjet() is different from original code.
      common/jyjets/njet,npad,kjet(1000,5),pjet(1000,5),vjet(1000,5)
      common/jampos1/jqconst(2),kfcq(4),icq(4),icms
      save /jyjets/,/jampos1/

      real*8 jamdtim
      character*16 charn
      dimension xo(4),delt(4)
      dimension indd(100),idel(100)

c...Save some variables.
      icon=0
      nmeson0=nmeson
      nbary0=nbary
      nv0=nv
      k1old=k(1,ip)
      kf0=k(2,ip)
      k05=k(5,ip)
      k06=k(6,ip)

c...Parton system
      if(k1old.eq.4) then
 
        if(k(10,ip).eq.0.or.k(11,ip).eq.0) then
          write(check(1),8000)ip,k(2,ip),p(5,ip)
 8000     format('ip kf p5',i9,1x,i9,1x,g12.3)
          call jamerrm(30,1,'(jamjdec:) parton system k10 k11=0')
        endif

       ntp=3
       ibar0=0
       ibar0a=0
       px=0.0d0
       py=0.0d0
       pz=0.0d0
       e0=0.0d0
       em0=0.0d0
       nch0=0
       do j=1,4
        xo(j)=0.0d0
       end do
       mpa=k(11,ip)-k(10,ip)+1
c....Save total momentum and charge.
       do i=k(10,ip),k(11,ip) 
        nch0=nch0+jamk(1,i)
        px=px+p(1,i)
        py=py+p(2,i)
        pz=pz+p(3,i)
        e0=e0+p(4,i)
        do j=1,4
          xo(j)=xo(j)+r(j,i)/mpa
        end do
       end do
       em0=sqrt(e0**2-px**2-py**2-pz**2)
       bex=px/e0
       bey=py/e0
       bez=pz/e0
       gam=e0/em0
       ks0=0
       k10=k(10,ip)
       k11=k(11,ip)
       jqconst(1)=1
       jqconst(2)=1
       if(abs(k(2,k10)).gt.1000) jqconst(1)=2
       if(abs(k(2,k11)).gt.1000) jqconst(2)=2
       ks0=jqconst(1)+jqconst(2)
       bvx=vq(1,k10)
       bvy=vq(2,k10)
       bvz=vq(3,k10)
       gmv=vq(4,k10)
       phi=vq(5,k10)
       the=vq(6,k10)

c...Hadron
      else
        ibar0=k(9,ip)
        ibar0a=abs(ibar0)
        nch0=jamk(1,ip)
        ntp=1
        px=vq(1,ip)+vq(6,ip)
        py=vq(2,ip)+vq(7,ip)
        pz=vq(3,ip)+vq(8,ip)
        e0=vq(4,ip)+vq(9,ip)
        em0=sqrt(e0**2-px**2-py**2-pz**2)
        bex=px/e0
        bey=py/e0
        bez=pz/e0
        gam=e0/em0
        xo(1)=r(1,ip)
        xo(2)=r(2,ip)
        xo(3)=r(3,ip)
        xo(4)=r(4,ip)

        bvx=p(1,ip)
        bvy=p(2,ip)
        bvz=p(3,ip)
        phi=v(1,ip)
        the=v(2,ip)
        gmv=v(3,ip)

c...Count original const. quarks.
        if(k1old.ge.3) then
            if(ibar0a.eq.3) then
             ks0=3
             jqconst(1)=1
             jqconst(2)=2
            else if(ibar0.eq.0) then
             jqconst(1)=1
             jqconst(2)=1
             ks0=2
            else
             write(check(1),8100)ibar0,k1old,kf0
 8100        format('ibar0 k1old kf0',i4,1x,i4,1x,i9)
             call jamerrm(30,1,'(jamjdec:)invalid baryon #')
            endif

c....Some quarks are vertual.
        else
            ks0=mod(abs(k1old)/10,10)
            if(ks0.eq.1.or.ks0.eq.2) then
              jqconst(1)=0
              jqconst(2)=ks0
            else if(ks0.eq.3) then
              ks0=2
              jqconst(1)=1
              jqconst(2)=1
            else
             write(check(1),'(''k1old'',i4)')k1old
             call jamerrm(30,1,'(jamjdec:)funny ks')
            endif
        endif
      endif

c...Fragment the string systems.
      ntry=0
1000  continue
      ierror=0
      nmeson=nmeson0
      nbary=nbary0
      mdel=0
      nv=nv0
      ntry=ntry+1
      if(ntry.ge.3) then
c       call jamlist(1)
        write(mstc(38),*)'ip k',ip,(k(j,ip),j=1,11)
        call jamerrm(30,0,'string decay(jamjdec) ntry>3')
      endif
c----------------------------------------------------------------------

c...Fragment jet system.
      call jamfrg(kf0,ip,ntp,ierror)

c...Check errors.
      if(mstu(24).ne.0.or.ierror.gt.0) then
        ih=mstc(38)
        write(ih,*) 'Error occured(jdecy), repeat the event ibar',ibar0
        write(ih,*) 'mstu24 mstu28 ierror',mstu(24),mstu(28),ierror,ntry
        write(ih,*)'ind kf nch0 em0',ip,kf0,nch0,em0
        call pjlist(1)
        mstu(24)=0
        mstu(28)=0
        go to 1000
      endif

c...Remove already hadronized partons.
       if(ntp.eq.3) then
         do i=k(10,ip),k(11,ip) 
           if(k(2,i).lt.10) mstd(83)=mstd(83)-1
           if(abs(k(2,i)).gt.1000) mstd(83)=mstd(83)-2
           if(k(2,i).eq.21) mstd(84)=mstd(84)-1
           mstd(30)=mstd(30)+1
           call jamzero(i)
           k(1,i)=31
           if(mstc(6).ge.0) call jamcupda(i,-1,0)
         end do
       else if(ibar0a.eq.0) then
          mstd(30)=mstd(30)+1
          call jamzero(ip)
          k(1,ip)=42
          if(mstc(6).ge.0) call jamcupda(ip,-1,0)
       endif

c...Copy PYTHIA array to the JAM array.
      psum1=0.0d0
      psum2=0.0d0
      psum3=0.0d0
      psum4=0.0d0
      npar=0
      nch1=0
      n_str=0
      ibtag=0
      ikaon=0
      jcq=0

c...Loop over produced hadrons.
      do 390 i=1,njet
 
c......This is already removed particle.
        if(kjet(i,1).ge.10.or.kjet(i,1).eq.0) then
           n_str=n_str+1
           go to 390
        endif

        npar=npar+1

c...Convert Ks, KL into k0, k0bar
        if(kjet(i,2).eq.130.or.kjet(i,2).eq.310) then
             if(ikaon.eq.0) then
                kjet(i,2)=311
                ikaon=1
                if(rn(0).gt.0.5d0) then
                  kjet(i,2)=-311
                  ikaon=-1
                endif
             else if(ikaon.eq.1) then
               kjet(i,2)=-311
               ikaon=0
             else if(ikaon.eq.-1) then
               kjet(i,2)=311
               ikaon=0
             endif
        endif

c....Particle code
        kf1=kjet(i,2)
        kc1=jamcomp(kf1)
c.....Baryon number
        iba=abs(kchg(kc1,6)*isign(1,kf1))
        ibb=1
c...Mother is a baryonic string
        if(ibar0a.eq.3) then
          if(iba.ne.0.and.ibtag.eq.0) then
           indx=ip
           ibtag=1
           ibb=0
          else if(iba.eq.3) then
           ibb=2
          endif
        else if(iba.eq.3) then
          ibb=2
        endif

c....Meson daughter
        if(ibb.eq.0) then
        else if(ibb.eq.1.or.mstc(4).ge.10) then
c       if(ibb.eq.1) then
          nmeson=nmeson+1
          nv=nv+1
          if(nv.gt.mxv) then
            call jamerrm(30,0,'(jamjdec:)Particle too large [mxv]')
          endif
          indx=nv

c.....Baryon daughter
        else if(ibb.eq.2) then
          call jamindb(indx,indd,npar,mdel,idel)
        endif

        if(nv.ne.nbary+nmeson) then
          call jamlist(1)
          call jamerrm(30,0,'(jamjdec:)nv .ne. nbary+nmeson')
        endif

cxxxxxxxxxxxxxxxxxxxx
c       if(ibb.eq.2) then
c         print *,'indx=',indx,'npar=',npar,'kf=',kf1,'nv=',nv
c         pause
c       endif
cxxxxxxxxxxxxxxxxxxxx


c...Updat partilce array.
c---------------------------
        indd(npar)=indx
        call jamkupda(2,indx,kf1,kc1,1,k06,k05,0)
        k1=k(1,indx)

        p(1,indx)=pjet(i,1)
        p(2,indx)=pjet(i,2)
        p(3,indx)=pjet(i,3)
        p(4,indx)=sqrt(pjet(i,5)**2
     $            +pjet(i,1)**2+pjet(i,2)**2+pjet(i,3)**2)
        p(5,indx)=pjet(i,5)

c.....Gamma factor.
        if(p(5,indx).ge.0.00001d0) then
          gg=p(4,indx)/p(5,indx)
        else
          gg=1.0d0
        endif

c...Reset quark properties.
        kq(1,indx)=999999
        kq(2,indx)=0
        vq(1,indx)=bvx
        vq(2,indx)=bvy
        vq(3,indx)=bvz
        vq(4,indx)=gmv
        vq(5,indx)=phi
        vq(6,indx)=the
        do l=7,10
         vq(l,indx)=0.0d0
        end do



c....Newly produced hadron.
        if(kjet(i,4).eq.0) then 

c....Set formation point of hadron.
          if(mstc(72).le.3) then
            do j=1,4
              r(j,indx)=xo(j)+vjet(i,j)/parc(54)
              v(j,indx)=r(j,indx)
            end do
            r(5,indx)=r(4,indx)
            wtime=vjet(i,4)/parc(54)/gg

          else

c...Option for constant formation time.
            wtime=-parc(55)*dlog(rn(0))
            r(5,indx)=xo(4)+wtime*gg
            r(4,indx)=r(5,indx)
            v(4,indx)=r(5,indx)
            do j=1,3
              deltx=wtime*gg*p(j,indx)/p(4,indx)
              r(j,indx)=xo(j)+deltx
              v(j,indx)=r(j,indx)
            end do

c.....Add small distance from the mother.
            if(indx.ne.ip) then
              deltx=parc(42)*sqrt(rn(0))
              cos1=1.d0-2.d0*rn(0)
              sin1=sqrt(1.d0-cos1**2)
              phi1=2*paru(1)*rn(0)
              dxr=deltx*sin1*cos(phi1)
              dyr=deltx*sin1*sin(phi1)
              dzr=deltx*cos1
              dtr=0.0d0
              call jamrobo(0.0d0,0.0d0,bex,bey,bez,gam,dxr,dyr,dzr,dtr)
              r(1,indx)=r(1,indx)+dxr
              r(2,indx)=r(2,indx)+dyr
              r(3,indx)=r(3,indx)+dzr
              r(4,indx)=r(5,indx)+dtr
              r(5,indx)=r(4,indx)
              do j=1,4
                v(j,indx)=r(j,indx)
              end do
            endif

          endif

          pard(88)=pard(88)+wtime
          mstd(56)=mstd(56)+1

c......This hadron contains original const. quarks.
        else

            if(mstc(72).le.3) then
               delt(1)=vjet(i,1)/parc(54)
               delt(2)=vjet(i,2)/parc(54)
               delt(3)=vjet(i,3)/parc(54)
               delt(4)=vjet(i,4)/parc(54)
               wtime=delt(4)/gg
            else
              wtime=-parc(55)*dlog(rn(0))
              delt(1)=wtime*gg*p(1,indx)/p(4,indx)
              delt(2)=wtime*gg*p(2,indx)/p(4,indx)
              delt(3)=wtime*gg*p(3,indx)/p(4,indx)
              delt(4)=wtime*gg
c.....Add small distance from the mother.
              if(indx.ne.ip) then
                deltx=parc(42)*sqrt(rn(0))
                cos1=1.d0-2.d0*rn(0)
                sin1=sqrt(1.d0-cos1**2)
                phi1=2*paru(1)*rn(0)
                dxr=deltx*sin1*cos(phi1)
                dyr=deltx*sin1*sin(phi1)
                dzr=deltx*cos1
                dtr=0.0d0
                call jamrobo(0.0d0,0.0d0,bex,bey,bez,gam,dxr,dyr,dzr, 
     & dtr)
                delt(1)=delt(1)+dxr
                delt(2)=delt(2)+dyr
                delt(3)=delt(3)+dzr
                delt(4)=delt(4)+dtr
               endif
            endif

c.......Additive quark cross sectin within t_form
          if(mstc(53).eq.1) then
            do j=1,4
c             r(j,indx)=xo(j)+delt(j)-delt(4)*p(j,indx)/p(4,indx)
              r(j,indx)=xo(j)
              v(j,indx)=xo(j)+delt(j)
            end do
            v(4,indx)=r(4,indx)
            r(5,indx)=xo(4)+delt(4)
c.......k(i,4) is a original const. quark number.
            k(1,indx)=-(k1+kjet(i,4)*10)

c...Option for no collsion within a formation time.
          else if(mstc(53).eq.2) then
            do j=1,4
              r(j,indx)=xo(j)+delt(j)
              v(j,indx)=xo(j)+delt(j)
            end do
            r(5,indx)=r(4,indx)

c...Option for no formation time.
          else
            do j=1,4
              r(j,indx)=xo(j)
              v(j,indx)=xo(j)
            end do
            r(5,indx)=xo(4)
            wtime=0.0d0
          endif

          pard(89)=pard(89)+wtime
          mstd(57)=mstd(57)+1

          iqnum=kjet(i,4)
          if(kjet(i,4).eq.3) iqnum=2
          jcq=jcq+iqnum

        endif

c...Life time in the case of resonance.
        v(5,indx)=r(5,indx)+jamdtim(1,kf1,kc1,k1,p(5,indx),p(4,indx))

c...Option: print information.
        if(mstc(8).ge.3) then
          call pjname(kf1,charn)
          write(mstc(38),810)indx,k(1,indx),k(2,indx),
     $      kjet(i,4)
     $    ,kjet(i,5),p(5,indx),charn(1:8),ks0,r(4,indx),r(5,indx)
 810      format('(STR DECAY) (',i4,')',i3,1x,i6,
     $    1x,2i2,1x,f7.4,1x,a8,1x,'ks0=',i2,1x,e8.3,1x,e8.3)
        endif


c...Total momentum and charge.
        nch1=nch1+jamk(1,indx)
        psum1=psum1+pjet(i,1)
        psum2=psum2+pjet(i,2)
        psum3=psum3+pjet(i,3)
        psum4=psum4+pjet(i,4)
cTTTTTTTTTTTTTTTTTTTTTTTTTT
        if(kjet(i,4).eq.2) mstd(198)=mstd(198)+1

390   continue 

      if(ks0.ne.jcq) then
        write(check(1),8200)ks0,jcq,(jqconst(i),i=1,2),k1old,kf0,em0
 8200   format('ks0',i4,' jcq',i3,' jqcnst k1 kf em',
     $ 2i3,1x,i4,1x,i7,1x,g12.3)
        if(mstc(8).ge.2)then
          call pjlist(2)
          call jamerrm(1,1,'(jamjdec:) ks0 jcq jqconst')
        endif
      endif

c...Updat collision.
      if(mstc(6).ge.0.and.mstc(52).ne.1) then
      if(mdel.ge.1) then
       do i=1,mdel
         call jamcupda(idel(i),-1,1)
       end do
      endif
      endif

      if(mstc(8).le.1) return
c=======================================================================
c....Check possible errors after jet fragmentation.
      ierr=0
      if(nch1.ne.nch0) then
        ih=mstc(38)
        write(ih,*)'--------------------------------------------------'
        write(ih,*)'charge not conserved in jamjdec ind ntp=',ip,ntp
        write(ih,*)'k3 nch1 kf0 em0= ',nch0,nch1,kf0,em0,jamchge(kf0)
        write(ih,*)'--------------------------------------------------'
        ierr=1
      endif

      sf=sqrt(psum4**2-psum1**2-psum2**2-psum3**2)
      if(     (abs(psum1-px).gt.0.2d0)
     $        .or.(abs(psum2-py).gt.0.2d0)
     $           .or.(abs(psum3-pz).gt.0.2d0)
     $             .or.(abs(psum4-e0).gt.0.2d0) ) then
         ih=mstc(38)
         write(ih,*)'-------------------------------------'
         write(ih,*)'p not conserved after jet dec.'
         write(ih,*)'p   =',px,py,pz,e0,em0
         write(ih,*)'psum=',psum1,psum2,psum3,psum4,sf
         write(ih,*)'id em0 s=',kf0,em0,sf
         write(ih,*)'-------------------------------------'
         if(ierr.eq.1) then
           ierr=3
         else
           ierr=2
         endif
      endif

      if(ierr.ge.1) then
         call pjlist(1)
         ih=mstc(38)
         write(ih,*)'========================='
         write(ih,*)'njet=',njet
         write(ih,*)'========================='
         do i=1,njet
           kf=kjet(i,2)
           kc=jamcomp(kf)
           call pjname(kf,charn)
           write(ih,*)'-------------------------'
           write(ih,*)'kc itp nres e m',
     $                kc,pjet(i,4),pjet(i,5),charn
           write(ih,*)'pjet',(pjet(i,l),l=1,3)
           write(ih,*)'vjet',(vjet(i,l),l=1,4)
           write(ih,*)'-------------------------'
         enddo
         if(ierr.eq.1)
     $   write(ih,*)'after string decay charge not conserved'
         if(ierr.eq.2)
     $   write(ih,*)'after string decay energy not conserved'
         if(ierr.eq.3)
     $   write(ih,*)'after string decay p and charge not conserved'
         icon=999
         call jamerrm(30,0,'(jdec:) after string decay error')
      endif

      end

c***********************************************************************

      subroutine jamfrg(kf0,jtp,ntp,ierror)

c...Purpose: to fragment all leading strings
c... ntp=1: fragment proj string,
c... ntp=2: targ string, 
c... ntp=3: independent strings from jets.
c... jtp  : the line number of the string.

      include 'jam1.inc'
      include 'jam2.inc'

c...HIJING common block
      common/hijdat/hidat0(10,10),hidat(10)
      common/jyjets/njet,npad,kjet(1000,5),pjet(1000,5),vjet(1000,5)
      save  /hijdat/,/jyjets/
c...Local values.
      dimension pcm(5)
      double precision dps(5)
      common/jampos1/jqconst(2),kfcq(4),icq(4),icms
      save /jampos1/
        
c...Initialize the document lines.
      ierror=0
      njet=0
      icms=0

c...Fragment independent strings from jets.
      if(ntp.eq.3) then
 
        n1=k(10,jtp)
        n2=k(11,jtp) 
        idq1=0
        idq2=0
        if(abs(k(2,n1)).gt.1000) idq1=1
        if(abs(k(2,n2)).gt.1000) idq2=1
        if(idq1.eq.1.and.idq2.eq.0) then
          m1=n2
          m2=n1
          istep=-1
          itmp=jqconst(1) 
          jqconst(1)=jqconst(2)
          jqconst(2)=itmp
        else
          m1=n1
          m2=n2
          istep=1
        endif

        do 100 i=m1,m2,istep
          njet=njet+1
          kjet(njet,1)=2
          kjet(njet,2)=k(2,i)
          kjet(njet,3)=0
          kjet(njet,4)=0
          kjet(njet,5)=0
          pjet(njet,1)=p(1,i)
          pjet(njet,2)=p(2,i)
          pjet(njet,3)=p(3,i)
          pjet(njet,4)=p(4,i)
          pjet(njet,5)=p(5,i)
 100    continue
        kjet(njet,1)=1
        goto 5000
      endif

c=======================================================================
c....Copy string momentum into Pythia array.
 1000 do j=3,5
        kjet(1,j)=0
        kjet(2,j)=0
      end do
      njet=2
      kjet(1,1)=2
      kjet(1,2)=kq(1,jtp)
      pjet(1,1)=vq(1,jtp)
      pjet(1,2)=vq(2,jtp)
      pjet(1,3)=vq(3,jtp)
      pjet(1,4)=vq(4,jtp)
      pjet(1,5)=vq(5,jtp)
      kjet(2,1)=1
      kjet(2,2)=kq(2,jtp)
      pjet(2,1)=vq(6,jtp)
      pjet(2,2)=vq(7,jtp)
      pjet(2,3)=vq(8,jtp)
      pjet(2,4)=vq(9,jtp)
      pjet(2,5)=vq(10,jtp)

c------------------------------------------------------------------
      jetot=0
c...Conduct soft radiations.
      if(mstc(74).gt.0.and.
     $          (abs(kjet(1,2)).gt.10.or.abs(kjet(2,2)).gt.10)) then

        is=0
20      is=is+1
        if(is.eq.10) goto 30
        if(hidat0(10,is).le.pare(2)) goto 20
30      if(is.eq.1) is=2
        do j=2,3
           hidat(j)=hidat0(j,is-1)+(hidat0(j,is)-hidat0(j,is-1))
     &     *(pare(2)-hidat0(10,is-1))/(hidat0(10,is)-hidat0(10,is-1))
        end do

        if(rn(0).le.hidat(3)) then
c         hidat2=hidat(2)
          hidat2=2.0d0
          cutoff=parc(56)
c         if(ihpr2(8).eq.0.and.ihpr2(3).eq.0.and.ihpr2(9).eq.0)
c    &       hidat2=2.0d0
          if(pare(2).ge.1000.0d0.and.jetot.eq.0)then
            hidat2=3.0d0
            cutoff=5.0d0
          endif
          call attrad(ierror,hidat2,cutoff,parc(66))

        else if(jetot.eq.0.and.pare(2).ge.1000.0d0
     $                                      .and.rn(0).le.0.8d0) then
c         hidat2=3.0d0
          hidat2=2.0d0
          cutoff=5.0d0
c         if(ihpr2(8).eq.0.and.ihpr2(3).eq.0.and.ihpr2(9).eq.0)
c    &    hidat2=2.0d0
          call attrad(ierror,hidat2,cutoff,parc(66))
        endif
        if(ierror.ne.0) goto 1000

c.....Check energy-momentum conservation.
        pc1=0.0d0
        pc2=0.0d0
        pc3=0.0d0
        pc4=0.0d0
        do i=1,njet
          pc1=pc1+pjet(i,1)
          pc2=pc2+pjet(i,2)
          pc3=pc3+pjet(i,3)
          pjet(i,4)=
     $    sqrt(pjet(i,5)**2+pjet(i,1)**2+pjet(i,2)**2+pjet(i,3)**2)
          pc4=pc4+pjet(i,4)
        end do
        srt=sqrt(pc4**2-(pc1**2+pc2**2+pc3**2))
        if(abs(srt-p(5,jtp)).gt.0.01d0) then
         write(check(1),'(''srt p5'',g12.3,1x,g12.3)')srt,p(5,jtp)
         call jamerrm(1,1,'after attrad energy not conserved')
         ierror=1
         goto 1000
        endif

      endif
c=======================================================================
c....Boost string to its cms.
 5000 continue
      ierror=0
      icms=1
      nsave=njet
      do j=1,5
      dps(j)=0.0d0
      end do
      do i=1,njet
        do j=1,5
          if(j.ne.4) dps(j)=dps(j)+pjet(i,j)
        enddo
        dps(4)=dps(4)
     $    +sqrt(pjet(i,1)**2+pjet(i,2)**2+pjet(i,3)**2+pjet(i,5)**2)
      end do

c     mbst=0
      mstu(33)=1
      call pjrobo(1,njet,0.d0,0.d0
     $            ,-dps(1)/dps(4),-dps(2)/dps(4),-dps(3)/dps(4))

c...Rotation so that p_x=p_y=0.
      phi=pjangl(pjet(1,1),pjet(1,2))
      call pjrobo(1,njet,0.0d0,-phi,0d0,0d0,0d0)
      pxy=sqrt(pjet(1,1)**2+pjet(1,2)**2)
      the=pjangl(pjet(1,3),pxy)
      call pjrobo(1,njet,-the,0.d0,0d0,0d0,0d0)
        
c...Fragment jet system.
      call pjexec
      mstu(33)=0
      call pjrobo(nsave+1,njet,the,phi,dps(1)/dps(4),dps(2)/dps(4),
     &              dps(3)/dps(4))

c...Write particle list.
      if(mstc(8).ge.5) call pjlist(1)

      end

c***********************************************************************
c                                                                      *
c                String related utilites                               *
c                                                                      *
c***********************************************************************

      subroutine jamidres(kf0,emr,icon)

c...Purpose: to determine resonance ID corresponding to mass. 

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      dimension ncount(30)
      dimension kfmm(3)

      icon=0
      kf00=kf0
      kfsign=isign(1,kf00)
      kc0=jamcomp(kf0)
      if(kc0.ge.1) then
        id0=kchg(kc0,5)
        iz0=kchg(kc0,1)/3
        ibary=kchg(kc0,6)
      else
       write(check(1),'(''kf0='',i9)')kf0
       call jamerrm(30,1,'(jamidres:) Unrecognized KF code')
       icon=1
       return
      endif


      if(ibary.eq.3) then
      if(id0.eq.id_nucl.or.id0.eq.id_nucls) then
c....Find N* channel.
        if(emr.le.1.2d0) then
          if(iz0.eq. 0) kf0=2114*kfsign
          if(iz0.eq. 1) kf0=2214*kfsign
          return
        else if(emr.le.1.35d0) then
          if(iz0.eq. 0) kf0=2114*kfsign
          if(iz0.eq. 1) kf0=2214*kfsign
          return
        endif
        iof=iz0
        kcmin=mstc(22)+iof
        kcmax=mstc(23)+iof
        istep=2
        if(emr.ge.pmas(kcmax,1)+pmas(kcmax,3)) then
          kf0=kchg(kcmax,4)*isign(1,kf00)
          return
        endif
c...Delta*
      else if(id0.eq.id_delt.or.id0.eq.id_delts) then
        if(emr.le.1.5d0) then
          if(iz0.eq.-1) kf0=1114*kfsign
          if(iz0.eq. 0) kf0=2114*kfsign
          if(iz0.eq. 1) kf0=2214*kfsign
          if(iz0.eq. 2) kf0=2224*kfsign
          return
        endif
        iof=iz0+1
        kcmin=mstc(24)+iof
        kcmax=mstc(25)+iof
        istep=4
        if(emr.ge.pmas(kcmax,1)+pmas(kcmax,3)) then
          kf0=kchg(kcmax,4)*kfsign
          return
        endif
c...Lambda*
      else if(id0.eq.id_lamb.or.id0.eq.id_lambs) then
        kcmin=mstc(26)
        kcmax=mstc(27)
        istep=1
        if(emr.le.pmas(kcmin,1)) then
          kf0=kchg(kcmin,4)*kfsign
        endif
        if(emr.ge.pmas(kcmax,1)+pmas(kcmax,3)) then
          kf0=kchg(kcmax,4)*kfsign
          return
        endif
c...Sigma*
      else if(id0.eq.id_sigm.or.id0.eq.id_sigms) then
        iof=1+iz0
        kcmin=mstc(28)+iof
        kcmax=mstc(29)+iof
        istep=3
        if(emr.le.pmas(kcmin,1)-pmas(kcmin,3)) then
           if(iz0.eq.-1) kf0=3114*kfsign
           if(iz0.eq. 0) kf0=3214*kfsign
           if(iz0.eq. 1) kf0=3224*kfsign
           return
        endif
        if(emr.ge.pmas(kcmax,1)+pmas(kcmax,3)) then
          kf0=kchg(kcmax,4)*kfsign
          return
        endif
c...Xi*
      else if(id0.eq.id_xi.or.id0.eq.id_xis) then
        iof=1+iz0
        kcmin=mstc(30)+iof
        kcmax=mstc(31)+iof
        istep=2
        if(emr.le.pmas(kcmin,1)-pmas(kcmin,3)) then
           if(iz0.eq.-1) kf0=3314*kfsign
           if(iz0.eq. 0) kf0=3324*kfsign
           return
        endif
        if(emr.ge.pmas(kcmax,1)+pmas(kcmax,3)) then
          kf0=kchg(kcmax,4)*kfsign
          return
        endif
      else
        return
      endif

c...Mesons.
      else if(ibary.eq.0) then

        kfa=abs(kf0)
        kf1=mod(kfa/100,10)
        kf2=mod(kfa/10,10)
        kfm=100*kf1+10*kf2
        if(kchg(kc0,3).eq.0) then
          kfsign=1
        endif
        itry=0
1000    continue
        icount=0
        do ir=0,5
        do is=1,7,2
          kfmes=ir*10000+kfm+is
          kcm=jamcomp(kfmes)
          if(kcm.ne.0) then
            if( emr.ge.pmas(kcm,1)-pmas(kcm,3)
     $              .and.emr.lt.pmas(kcm,1)+pmas(kcm,3)) then
             icount=icount+1
             ncount(icount)=kfmes
             endif
           endif
        end do
        end do

        if(icount.ge.1) then
          ik=1+rn(0)*(icount-1)
          kf0=ncount(ik)*kfsign
        else
          itry=itry+1
          if(kfm.eq.220) then
             kfm=330
          else if(kfm.eq.110) then
             kfm=220
          else if(kfm.eq.330) then
             kfm=110
          endif
          if(itry.le.2) goto 1000

          icon=10
          return
        endif

        if(kchg(jamcomp(kf0),3).eq.0) kf0=abs(kf0)
        return

      else
        write(check(1),'(''kf0='',i9)')kf0
        call jamerrm(30,1,'(jamidres:) kf0=')
      endif

c....Find baryon resonance KF code.
      if(emr.lt.pmas(kcmin,1)-pmas(kcmin,3)) then
        write(check(1),'(''kf0 emr'',i9,1x,g13.3)')kf0,emr
        write(check(2),'(''kcmin kfmin'',i9,1x,i9)')kcmin,kchg(kcmin,4)
        write(check(3),'(''pmas1 pmas3'',2(g13.3,1x))')pmas(kcmin,1),
     $  pmas(kcmin,3)
        call jamerrm(1,3,'(jamidres:) Mass too small emr')
        icon=1
        return
      endif

      icount=0
      do i=kcmin,kcmax,istep
       if(emr.ge.pmas(i,1)-pmas(i,3)
     $           .and.emr.lt.pmas(i,1)+pmas(i,3)) then
         icount=icount+1
         ncount(icount)=i
       endif
      end do

      if(icount.ge.1) then
        ik=1+rn(0)*(icount-1)
        kf0=kchg(ncount(ik),4)*kfsign
      else
        kf0=kchg(kcmin,4)*kfsign
        icon=10
      endif

      end

c***********************************************************************

      subroutine jamkfres(kf0,kf1)

c...Purpose: to determine baryon resonance ID corresponding to mass. 
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

      kfsign=isign(1,kf0)
      kc0=jamcomp(kf0)
      if(kc0.ge.1) then
        id0=kchg(kc0,5)
        iz0=kchg(kc0,1)/3
        ibary=kchg(kc0,6)
      else
       write(check(1),'(''kf0='',i9)')kf0
       call jamerrm(1,1,'(jamidres:) Unrecognized KF code')
       icon=1
       return
      endif

      if(id0.eq.id_nucl.or.id0.eq.id_nucls) then
c....Find N* channel.
        kcmin=mstc(22)+iz0
        kcmax=mstc(23)+iz0
        istep=2
c          kf1=kchg(kcmin,4)*kfsign
c          return
c...Delta*
      else if(id0.eq.id_delt.or.id0.eq.id_delts) then
        kcmin=mstc(24)+iz0+1
        kcmax=mstc(25)+iz0+1
        istep=4
c          kf1=kchg(kcmin,4)*kfsign
c          return
c...Lambda*
      else if(id0.eq.id_lamb.or.id0.eq.id_lambs) then
        kcmin=mstc(26)
        kcmax=mstc(27)
        istep=1
c...Sigma*
      else if(id0.eq.id_sigm.or.id0.eq.id_sigms) then
        kcmin=mstc(28)+iz0+1
        kcmax=mstc(29)+iz0+1
        istep=3
c...Xi*
      else if(id0.eq.id_xi.or.id0.eq.id_xis) then
        kcmin=mstc(30)+iz0+1
        kcmax=mstc(31)+iz0+1
        istep=2
      else
        kf1=kf0
        return
      endif

c....Find baryon resonance KF code.
c     tot=(kcmax-kcmin)/dble(istep)+1
      tot=0.0d0
      do i=kcmin,kcmax,istep
        tot=tot+max(1,mod(kchg(i,4),10))
      end do
      xran=tot*rn(0)
      do i=kcmin,kcmax,istep
        xran=xran-max(1,mod(kchg(i,4),10))
        if(xran.le.0.0d0) then
           ik=i
           goto 100
        endif
      end do
      kf1=kchg(kcmax,4)*kfsign
      return
100   continue
      kf1=kchg(ik,4)*kfsign

      end

c***********************************************************************

      subroutine jamrobo(the,phi,bex,bey,bez,dga,p1,p2,p3,p4)    
    
c...Purpose: to perform rotations and boosts.   

      implicit double precision(a-h, o-z)
      dimension rot(3,3),pr(3),dp(4)
    
c...Rotate, typically from z axis to direction (theta,phi). 
      if(the**2+phi**2.gt.1d-20) then   
        rot(1,1)=cos(the)*cos(phi)  
        rot(1,2)=-sin(phi)  
        rot(1,3)=sin(the)*cos(phi)  
        rot(2,1)=cos(the)*sin(phi)  
        rot(2,2)=cos(phi)   
        rot(2,3)=sin(the)*sin(phi)  
        rot(3,1)=-sin(the)
        rot(3,2)=0.d0 
        rot(3,3)=cos(the)
        pr(1)=p1
        pr(2)=p2
        pr(3)=p3
        p1=rot(1,1)*pr(1)+rot(1,2)*pr(2)+rot(1,3)*pr(3) 
        p2=rot(2,1)*pr(1)+rot(2,2)*pr(2)+rot(2,3)*pr(3) 
        p3=rot(3,1)*pr(1)+rot(3,2)*pr(2)+rot(3,3)*pr(3) 
      endif 
    
c...Boost, typically from rest to momentum/energy=beta. 
      if(bex**2+bey**2+bez**2.gt.1d-20) then    
        dbx=bex
        dby=bey
        dbz=bez
        db=sqrt(dbx**2+dby**2+dbz**2)   
        eps1=1d0-1d-12
c...Rescale boost vector if too close to unity.
        if(db.gt.eps1) then
          call jamerrm(3,0,'(jamrobo:) boost vector too large')
          dbx=dbx*(eps1/db)
          dby=dby*(eps1/db)
          dbz=dbz*(eps1/db)
          db=eps1
        endif
c       dga=1d0/sqrt(1d0-db**2) 
        dp(1)=p1
        dp(2)=p2
        dp(3)=p3
        dp(4)=p4
        dbp=dbx*dp(1)+dby*dp(2)+dbz*dp(3)   
        dgabp=dga*(dga*dbp/(1d0+dga)+dp(4)) 
        p1=dp(1)+dgabp*dbx  
        p2=dp(2)+dgabp*dby  
        p3=dp(3)+dgabp*dbz  
        p4=dga*(dp(4)+dbp)  
      endif 
    
      end   

c***********************************************************************

      subroutine jamindb(indx,indd,npar,mdel,idel)

c...Find index of newly produced baryon/antibaryon.
      include 'jam1.inc'
      include 'jam2.inc'
      dimension indd(100),idel(100)

      if(nv+1.gt.mxv.or.nbary+1.gt.mxv) 
     $  call jamerrm(30,0,'(jamindb:)Particle too large [mxv]')

      nbary=nbary+1

      if(k(1,nbary).gt.10.or.k(1,nbary).eq.0) then
          indx=nbary
          if(nbary.gt.nv) then
             nv=nv+1
          else
             nmeson=nmeson-1
          endif
      else if(k(1,nbary).eq.4) then

          if(nbary.ne.k(10,nbary)) then
            write(check(1),'(''nbary k1 k10'',3(i9,1x))')
     $                  nbary,k(1,nbary),k(10,nbary)
            call jamerrm(30,1,'k(1,nbary).ne.k(10,nabry) nbary')
          endif
          k10=k(10,nbary)
          k11=k(11,nbary)
          iof=k11-k10+1
          if(nv+iof.gt.mxv) 
     $      call jamerrm(30,0,'(jamindb:)Particle too large [mxv]')
          k10new=nv+1
          k11new=nv+iof
          nmeson=nmeson+iof-1
          do j=k10,k11
            nv=nv+1
            call jamexch(nv,j)
            k(10,nv)=k10new 
            k(11,nv)=k11new 
            call jamzero(j)
            k(1,j)=40
            if(mstc(6).ge.0) call jamcupda(j,-1,0)
            mdel=mdel+1
            idel(mdel)=nv
          end do
          indx=nbary
      else
          nv=nv+1
          call jamexch(nv,nbary)
          indx=nbary
          mdel=mdel+1
          idel(mdel)=nv
c bugfix 02/23/2002
          do ii=1,npar
            if(indd(ii).eq.nbary) then
              indd(ii)=nv
            endif
          end do
c bugfix end
      endif

      end

c***********************************************************************

      function xsamp1(xmin,xmax,srt)

c...Sampe x of valence quarks for baryon, DPM type.
      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      save  /hiparnt/

c...Fritiof type  distribuiton: 1/x
c     if(srt.lt.10.0d0) then
c       xsamp1=xmin*(xmax/xmin)**rn(0)
c       return
c     endif

c...p(x)=(1-x)^d/(x^2+cutt^2)^b
c     xsamp1=hirnd2(4,xmin,xmax)

      b=hipr1(46)
      b1=2*b
      d=hipr1(44)
      cutt=hipr1(45)/srt
11    continue
      cut1=(xmin+cutt)**(1.d0-b1)
      cut2=(xmax+cutt)**(1.d0-b1)
      x1=( (cut2-cut1)*rn(0)+cut1 )**(1.d0/(1.d0-b1))-cutt
      if((1.d0-x1)**d*((x1+cutt)**2/(2.d0*(x1**2+cutt**2)))**b.lt.rn(0))
     $     goto 11
      xsamp1=x1

      end

c***********************************************************************

      function xsamp2(xmin,xmax,srt)

c...Sampe x of valence quarks for mseon DPM type.
      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      save  /hiparnt/

c...1/x
c     if(srt.lt.10.0d0) then
c       xsamp2=xmin*(xmax/xmin)**rn(0)
c       return
c     endif

c...p(x)=1./((x^2+cutt^2)^b*((1-x)^2+cutt^2)^b)
c     xsamp2=hirnd2(5,xmin,xmax)
      b=hipr1(46)
      cutt=hipr1(45)/srt
      cutt1=asin((2*xmin-1.d0)/(1.d0+2*cutt))
      cutt2=asin((2*xmax-1.d0)/(1.d0+2*cutt))
      cutt3=cutt2-cutt1
11    continue
      bb=sin(rn(0)*cutt3+cutt1)
      x1=0.5d0*(bb*(1+2*cutt)+1)
       if(x1.lt.xmin.or.x1.gt.xmax) go to 11
      fstrum=1.0d0/((1.0d0-x1)**2+cutt**2)**b/(x1**2+cutt**2)**b
      g=1.3d0/sqrt((x1+cutt)*(1-x1+cutt))
      if(fstrum.lt.rn(0)*g) goto 11
      xsamp2=x1

      end

c***********************************************************************

      function xsamp3(xmin,xmax,srt)

c....Sample x for diffractive scattering.
      implicit double precision(a-h, o-z)
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      save  /hiparnt/

c     if(srt.lt.10.0d0) then
c       xsamp3=xmin*(xmax/xmin)**rn(0)
c       return
c     endif

c...1/x
c     xsamp3=xmin*(xmax/xmin)**rn(0)

c...p(x)=1/sqrt(x^2+(c/srt)^2)
c     xsamp3=hirnd2(6,xmin,xmax)
      cut=(hipr1(45)/srt)**2
      cut1=xmax+sqrt(xmax**2+cut)
      cut2=xmin+sqrt(xmin**2+cut)
      ar=cut2*(cut1/cut2)**rn(0)
      xsamp3=0.5d0*(ar-cut/ar)


      end

c***********************************************************************

      subroutine attrad(ierror,hidat2,cutoff,widpt)

c...Purpose: to conduct soft radiation according to dipole approxiamtion
 
      implicit double precision(a-h, o-z)
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      save  /jyjets/

      ierror=0
c.....s invariant mass-squared between partons i and i+1......
c.....sm is the largest mass-squared....

40      sm=0.d0
        do 30 i=1,n-1
           s=2.d0*(p(i,4)*p(i+1,4)-p(i,1)*p(i+1,1)-p(i,2)*p(i+1,2)
     &          -p(i,3)*p(i+1,3))+p(i,5)**2+p(i+1,5)**2
           if(s.lt.0.d0) s=0.d0
           wp=sqrt(s)-1.5d0*(p(i,5)+p(i+1,5))
           if(wp.gt.sm) then
              pbt1=p(i,1)+p(i+1,1)
              pbt2=p(i,2)+p(i+1,2)
              pbt3=p(i,3)+p(i+1,3)
              pbt4=p(i,4)+p(i+1,4)
              btt=(pbt1**2+pbt2**2+pbt3**2)/pbt4**2
              if(btt.ge.1.0d0-1.0d-10) go to 30
              if((i.ne.1.or.i.ne.n-1).and.
     &             (k(i,2).ne.21.and.k(i+1,2).ne.21)) go to 30
              jl=i
              sm=wp
           endif
30      continue
        s=(sm+1.5d0*(p(jl,5)+p(jl+1,5)))**2
        if(sm.lt.cutoff) goto 2
     
C.....Make place for one gluon.....
        if(jl+1.eq.n) goto 190
        do 160 j=n,jl+2,-1
          k(j+1,1)=k(j,1)
          k(j+1,2)=k(j,2)
          do 150 m=1,5
150         p(j+1,m)=p(j,m)
160     continue
190     n=n+1
     
C.....Boost to rest system for particles jl and jl+1.....
        p1=p(jl,1)+p(jl+1,1)
        p2=p(jl,2)+p(jl+1,2)
        p3=p(jl,3)+p(jl+1,3)
        p4=p(jl,4)+p(jl+1,4)
        bex=-p1/p4
        bey=-p2/p4
        bez=-p3/p4
        imin=jl
        imax=jl+1
        call atrobo(0.d0,0.d0,bex,bey,bez,imin,imax,ierror)
        if(ierror.ne.0) return

C.....Rotate to z-axis....
        cth=p(jl,3)/sqrt(p(jl,4)**2-p(jl,5)**2)
        if(abs(cth).gt.1.0d0)  cth=max(-1.d0,min(1.d0,cth))
        theta=acos(cth)
        phi=pjangl(p(jl,1),p(jl,2))
        call atrobo(0.d0,-phi,0.d0,0.d0,0.d0,imin,imax,ierror)
        call atrobo(-theta,0.d0,0.d0,0.d0,0.d0,imin,imax,ierror)
     
C.....Create one gluon and orientate.....
1       call ar3jet(s,x1,x3,jl)
        call arorie(s,x1,x3,jl)         
        if(hidat2.gt.0.0d0) then
           ptg1=sqrt(p(jl,1)**2+p(jl,2)**2)
           ptg2=sqrt(p(jl+1,1)**2+p(jl+1,2)**2)
           ptg3=sqrt(p(jl+2,1)**2+p(jl+2,2)**2)
           ptg=max(ptg1,ptg2,ptg3)
           if(ptg.gt.hidat2) then
              fmfact=exp(-(ptg**2-hidat2**2)/widpt**2)
              if(rn(0).gt.fmfact) go to 1
           endif
        endif

C.....Rotate and boost back.....
        imin=jl
        imax=jl+2
        call atrobo(theta,phi,-bex,-bey,-bez,imin,imax,ierror)
        if(ierror.ne.0) return

C.....Enumerate the gluons.....
        k(jl+2,1)=k(jl+1,1)
        k(jl+2,2)=k(jl+1,2)
        k(jl+2,3)=k(jl+1,3)
        k(jl+2,4)=k(jl+1,4)
        k(jl+2,5)=k(jl+1,5)
        p(jl+2,5)=p(jl+1,5)
        k(jl+1,1)=2
        k(jl+1,2)=21
        k(jl+1,3)=0
        k(jl+1,4)=0
        k(jl+1,5)=0
        p(jl+1,5)=0.d0

c----Theta function damping of the emitted gluons. for hadron-hadron.
c----r0=vfr(2)
c       if(vfr(2).gt.0.) then
c       ptg=sqrt(p(jl+1,1)**2+p(jl+1,2)**2)
c       ptgmax=wstri/2.
c       dopt=sqrt((4.*par(71)*vfr(2))/wstri)
c       ptopt=(dopt*wstri)/(2.*vfr(2))
c       if(ptg.gt.ptopt) iorder=iorder-1
c       if(ptg.gt.ptopt) goto 1
c       endif
c-----
        if(sm.ge.cutoff) goto 40

2       k(1,1)=2
        k(1,3)=0
        k(1,4)=0
        k(1,5)=0
        k(n,1)=1
        k(n,3)=0
        k(n,4)=0
        k(n,5)=0

        end

c***********************************************************************

      subroutine ar3jet(s,x1,x3,jl)

c...Calculate scaled energy variables of gluon and quark.
      implicit double precision(a-h, o-z)
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      save  /jyjets/

        c=1.d0/3.d0
        if(k(jl,2).ne.21 .and. k(jl+1,2).ne.21) c=8.d0/27.d0
        exp1=3
        exp3=3
        if(k(jl,2).ne.21) exp1=2
        if(k(jl+1,2).ne.21) exp3=2
        a=0.24d0**2/s
        yma=dlog(.5d0/sqrt(a)+sqrt(.25d0/a-1))
        d=4.d0*c*yma
        sm1=p(jl,5)**2/s
        sm3=p(jl+1,5)**2/s
        xt2m=(1.d0-2.d0*sqrt(sm1)+sm1-sm3)*(1.d0-2.d0*sqrt(sm3)-sm1+sm3)
        xt2m=min(.25d0,xt2m)
        ntry=0
1       if(ntry.eq.5000) then
                x1=.5d0*(2.d0*sqrt(sm1)+1.d0+sm1-sm3)
                x3=.5d0*(2.d0*sqrt(sm3)+1.d0-sm1+sm3)
                return
        endif
        ntry=ntry+1
     
        xt2=a*(xt2m/a)**(rn(0)**(1.d0/d))
     
        ymax=dlog(.5d0/sqrt(xt2)+sqrt(.25d0/xt2-1.d0))
        y=(2.d0*rn(0)-1.d0)*ymax
        x1=1.d0-sqrt(xt2)*exp(y)
        x3=1.d0-sqrt(xt2)*exp(-y)
        x2=2.d0-x1-x3
        neg=0
        if(k(jl,2).ne.21 .or. k(jl+1,2).ne.21) then
        if((1.d0-x1)*(1.d0-x2)*(1.d0-x3)-x2*sm1*(1.d0-x1)-x2*sm3*(1.d0 
     & -x3).
     &  le.0.d0.or.x1.le.2.d0*sqrt(sm1)-sm1+sm3.or.x3.le.2.d0*sqrt(sm3)
     &  -sm3+sm1) neg=1
        x1=x1+sm1-sm3
        x3=x3-sm1+sm3
        endif
        if(neg.eq.1) goto 1
     
        fg=2.d0*ymax*c*(x1**exp1+x3**exp3)/d
        xt2m=xt2
        if(fg.lt.rn(0)) goto 1
     
        return
        end

c***********************************************************************

      subroutine arorie(s,x1,x3,jl)

c...Gives pt to the soft radiative gluon.
      implicit double precision(a-h, o-z)
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      save  /jyjets/

        w=sqrt(s)
        x2=2.d0-x1-x3
        d1=.5d0*x1*w
        d3=.5d0*x3*w
        p1=sqrt(d1**2-p(jl,5)**2)
        p3=sqrt(d3**2-p(jl+1,5)**2)
        cbet=1.d0
        if(p1.gt.0.d0.and.p3.gt.0.d0) cbet=(p(jl,5)**2
     &           +p(jl+1,5)**2+2.d0*d1*d3-s*(1.d0-x2))/(2.d0*p1*p3)
        if(abs(cbet).gt.1.0d0) cbet=max(-1.d0,min(1.d0,cbet))
        bet=acos(cbet)
     
C...Minimize pt1-squared plus pt3-squared.....
        if(p1.ge.p3) then
           psi=.5d0*pjangl(p1**2+p3**2*cos(2.d0*bet), 
     & -p3**2*sin(2.d0*bet))
           pt1=p1*sin(psi)
           pz1=p1*cos(psi)
           pt3=p3*sin(psi+bet)
           pz3=p3*cos(psi+bet)
        else if(p3.gt.p1) then
           psi=.5d0*pjangl(p3**2+p1**2*cos(2.d0*bet), 
     & -p1**2*sin(2.d0*bet))
           pt1=p1*sin(bet+psi)
           pz1=-p1*cos(bet+psi)
           pt3=p3*sin(psi)
           pz3=-p3*cos(psi)
        endif
     
        del=2.0d0*3.1415926d0*rn(0)
        p(jl,4)=d1
        p(jl,1)=pt1*sin(del)
        p(jl,2)=-pt1*cos(del)
        p(jl,3)=pz1
        p(jl+2,4)=d3
        p(jl+2,1)=pt3*sin(del)
        p(jl+2,2)=-pt3*cos(del)
        p(jl+2,3)=pz3
        p(jl+1,4)=w-d1-d3
        p(jl+1,1)=-p(jl,1)-p(jl+2,1)
        p(jl+1,2)=-p(jl,2)-p(jl+2,2)
        p(jl+1,3)=-p(jl,3)-p(jl+2,3)
        return
        end

c***********************************************************************

      subroutine atrobo(the,phi,bex,bey,bez,imin,imax,ierror)

c...Purpose: to make  boost and rotation to entries from imin to imax

      implicit double precision(a-h, o-z)
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      save  /jyjets/

      dimension rot(3,3),pv(3)
      double precision dp(4),dbex,dbey,dbez,dga,dga2,dbep,dgabep

      ierror=0
      if(imin.le.0 .or. imax.gt.n .or. imin.gt.imax) return

        if(the**2+phi**2.gt.1d-20) then
C...Rotate (typically from z axis to direction theta,phi)
           rot(1,1)=cos(the)*cos(phi)
           rot(1,2)=-sin(phi)
           rot(1,3)=sin(the)*cos(phi)
           rot(2,1)=cos(the)*sin(phi)
           rot(2,2)=cos(phi)
           rot(2,3)=sin(the)*sin(phi)
           rot(3,1)=-sin(the)
           rot(3,2)=0.d0
           rot(3,3)=cos(the)

           do 120 i=imin,imax
C**************    if(mod(k(i,1)/10000,10).ge.6) goto 120
              do 100 j=1,3
 100             pv(j)=p(i,j)
                 do 110 j=1,3
 110                p(i,j)=rot(j,1)*pv(1)+rot(j,2)*pv(2)
     &                     +rot(j,3)*pv(3)
 120       continue

        endif
     
        if(bex**2+bey**2+bez**2.gt.1d-20) then
C...Lorentz boost (typically from rest to momentum/energy=beta)
                dbex=bex
                dbey=bey
                dbez=bez
                dga2=1d0-dbex**2-dbey**2-dbez**2
                if(dga2.le.0d0) then
                   write(3,*)'(atrobo)bex bey bez',bex,bey,bez,dga2
                   ierror=1
                   return
                endif
                dga=1d0/dsqrt(dga2)
                do 140 i=imin,imax
c*************     if(mod(k(i,1)/10000,10).ge.6) goto 140
                   do 130 j=1,4
130                dp(j)=p(i,j)
                   dbep=dbex*dp(1)+dbey*dp(2)+dbez*dp(3)
                   dgabep=dga*(dga*dbep/(1d0+dga)+dp(4))
                   p(i,1)=dp(1)+dgabep*dbex
                   p(i,2)=dp(2)+dgabep*dbey
                   p(i,3)=dp(3)+dgabep*dbez
                   p(i,4)=dga*(dp(4)+dbep)
140             continue
        endif
     
        end

c***********************************************************************

      subroutine str2had(jmod,ind,kf,kc,emj,ic)

c...Purpose: to check mass of string system and if mass is
c...small, it is converted to hadron.
c..jmod=1: only check
c..ic=0: remain string system
c..ic=1: converted hadron
c..ic=2: mass too small

      include 'jam1.inc'
      include 'jam2.inc'
      real*8 jamdtim,jamemjet
      dimension idel(100)

      ic=0
      in1=k(10,ind)
      in2=k(11,ind)
ccc   if(abs(in2-in1).gt.1) return
      if(k(2,in1).eq.21.and.k(2,in2).eq.21) return
      ind1=in1 
      ind2=in2 
      px=0.0d0
      py=0.0d0
      pz=0.0d0
      e0=0.0d0
      do i=ind1,ind2
          px=px+p(1,i)
          py=py+p(2,i)
          pz=pz+p(3,i)
          e0=e0+p(4,i)
      end do
      emj=sqrt(e0**2-(px**2+py**2+pz**2))
      kf1=k(2,ind1)
      kf2=k(2,ind2)

      if(abs(kf1).le.10.and.abs(kf2).le.10) then
        iba=0
      else
        iba=1
      endif
      emmi=jamemjet(kf1,kf2)

c...Check mass
      if(emj.gt.emmi) return

      if(iba.eq.1) then
        if(emj.lt.1.08d0) then
          write(mstc(38),*)'(str2had:)1 ind emj=',ind,emj
          ic=2
          return
        endif
      else
        if(emj.lt.2*parc(27)+0.0001d0) then
          write(mstc(38),*)'(str2had:)2 ind emj=',ind,emj
          ic=2
          return
        endif
      endif

      ic=1
      if(jmod.eq.1) return
       iii=1
       if(iii.eq.1) return

c...Convert string system to a hadron.
         mstd(30)=mstd(30)+1
         p(1,ind1)=px
         p(2,ind1)=py
         p(3,ind1)=pz
         p(4,ind1)=e0
         p(5,ind1)=emj
         k8=k(8,ind1)
         if(abs(kf1).le.10) then
           irev=0
         else if(abs(kf1).gt.10) then
           irev=1
           kftmp=kf1
           kf1=kf2
           kf2=kftmp
           do j=1,5
            ptmp=p(j,ind1)
            p(j,ind1)=p(j,ind2)
            p(j,ind2)=ptmp
           end do
         endif
         call kfcnst(kf1,kf2,kf,emj)
         if(kf.eq.0)then
           call jamerrm(30,0,'(str2had:)kf=0')
         endif
         k(2,ind1)=kf
         kq(1,ind1)=kf1
         kq(2,ind1)=kf2
      kc=jamcomp(kf)
      k(3,ind1)=10
      k(8,ind1)=k8
      k(9,ind1)=kchg(kc,6)*isign(1,kf)
      k(4,ind1)=10
      if(pmas(kc,2).ge.1.d-5.and.mdcy(kc,1).ne.0) then
        k(1,ind1)=2
      else
        k(1,ind1)=1
      endif

      do j=1,5
       vq(j,ind1)=p(j,ind1)
       vq(j+5,ind1)=p(j,ind2)
      end do
      do j=1,4
      r(j,ind1)=(r(j,ind1)+r(j,ind2))/2.d0
      v(j,ind1)=r(j,ind1)
      end do
      r(5,ind1)=r(4,ind1)
      v(5,ind1)=r(5,ind1)+jamdtim(1,kf,kc,k(1,ind1),p(5,ind1),p(4,ind1))

      if(mstc(8).ge.1.or.mstc(13).ge.1) then
      if(mstc(13).ge.2.or.
     $                 (mstc(13).eq.1.and.(mstd(25).lt.mstc(14)))) then
        mstd(25)=mstd(25)+1
        ih=mstc(38)
        write(ih,*)'convert into hadron',ind1,ind2,r(4,ind1),r(5,ind1)
        write(ih,*)'nv nbary nmeson emj emmi',nv,nbary,nmeson,emj,emmi
        write(ih,*)'convert into hadron',ind1,ind2,r(4,ind1),r(5,ind1)
      endif
      endif

      k(1,ind2)=50
      call jamzero(ind2)

c...Change line number of hadron in case of baryon.
      if(iba.eq.1) then
         if(nv+1.gt.mxv) 
     $           call jamerrm(30,0,'(jdecay:)Particle too large [mxv]')
         nbary=nbary+1

         if(nbary.gt.mxv) 
     $         call jamerrm(30,0,'(jdecay:)Particle too large [mxv]')

        mdel=0
        if(k(1,nbary).gt.10.or.k(1,nbary).eq.0) then
                indx=nbary
                if(nbary.gt.nv) then
                  nv=nv+1
                else
                  nmeson=nmeson-1
                endif
        else if(k(1,nbary).eq.4) then

          if(nbary.ne.k(10,nbary)) then
            write(mstc(38),*)'k(1,nbary).ne.k(10,nabry) nbary',
     $                              nbary,k(1,nbary),k(10,nbary)
            call jamerrm(30,0,'(str2had:)k1 k10?')
          endif

                k10=k(10,nbary)
                k11=k(11,nbary)
                iof=k11-k10+1
                if(nv+iof.gt.mxv) 
     $            call jamerrm(30,0,'(jdecay:)Particle too large [mxv]')
                k10new=nv+1
                k11new=nv+iof
                nmeson=nmeson+iof-1
                do j=k10,k11
                  nv=nv+1
                  mdel=mdel+1
                  idel(mdel)=nv
                  call jamexch(nv,j)
                  k(10,nv)=k10new 
                  k(11,nv)=k11new 
                  call jamzero(j)
                  k(1,j)=40
                  if(mstc(6).ge.0) call jamcupda(j,-1,0)
                end do
                indx=nbary
        else
               nv=nv+1
               call jamexch(nv,nbary)
               indx=nbary
               mdel=mdel+1
               idel(mdel)=nv
        endif

        call jamexch(indx,ind1)
        call jamzero(ind1)
        k(1,ind1)=40
        if(mstc(6).ge.0) then
        if(mstc(52).ne.1) then
          call jamcupda(ind1,ind2,0)
          call jamcupda(indx,-1,1)
          if(mdel.ge.1) then
            do i=1,mdel
              call jamcupda(idel(i),-1,1)
            end do
          endif
          write(mstc(38),*)'(str2had:) indx k',indx,k(1,indx),k(2,indx)
     $                                      ,nv,nbary,nmeson
        else
          call jamcupda(ind1,ind2,1)
        endif
        endif

      endif

      end

c***********************************************************************

      subroutine attflv(kf,ifla,iflb)

c...Purpse: to give spin and quarkflavour to the ends of
c...the excited strings.
c...For mesons, the order of the end flavors is randomly given;
c...For baryons, where a quark-diquark combination,
c...the diquark is always assigned to iflb.
c...Use SU(6) weight proton=1/3d(uu)1 + 1/6u(ud)1 + 1/2u(ud)0
c...                 nurtron=1/3u(dd)1 + 1/6d(ud)1 + 1/2d(ud)0

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      integer ifrkfcpj
      data vfr13,vfr14,vfr15/0.167d0,  .333d0,  .5d0/
c...for proton
c...vfr13 probability of finding a diquark ud with spin 1 1/6
c...vfr14 probability of finding a diquark uu with spin 1 1/3
c...vfr15 probability of finding a diquark ud with spin 0 1/2
c...54321
c...    x quark
c... xx0x di-quark
c...   1x leptons
c...   xx gauge and Higgs bosons
c...  xxx meson
c...x0xxx meson
c... xxxx baryon
c...xxxxx baryon

      j = abs(kf)
      jc1=mod(j,10)      ! spin
      jc2=mod(j/10,10)
      jc3=mod(j/100,10)
      jc4=mod(j/1000,10)
      jc5=mod(j/10000,10)

      ifla=0
      iflb=0
c...Quark
      if(j.lt.10) then
         return
c...Lepton in idtrns
      else if(j.lt.20) then
         return
c...Gauge or higgs boson
      else if(j.le.100) then
         return
      endif

c...Di-quark
      if(jc2.eq.0) then
         return
      endif

      id1=isign(1,kf)*(jc4*100+jc3*10+jc2)
      spin=pjr(0)

c..Mesons
      if(jc4.eq.0) then
c...Identify the quark and antiquark in mesons:
c       j100=  j/100
c       j10 = (j-j100*100)/10
        j100=jc3
        j10=jc2
        isgn = (-1)**max(j100, j10)
        if(kf.lt.0) isgn = -isgn 
        if(isgn.gt.0) j10 = -j10
        if(isgn.lt.0) j100 = -j100

        if(id1.eq.11) then  ! pi0 eta rho0 omega...
          ranx=pjr(0)
          if(ranx.lt.0.5d0) then
             j100=1
             j10=-1
          else
             j100=2
             j10=-2
          endif
        endif

        if(j.eq.331) then  ! eta'  (1/6 1/6 2/3)
          ranx=pjr(0)
          if(ranx.lt.0.1667d0) then
             j100=1
             j10=-1
          else if(ranx.lt.0.3333d0) then
             j100=2
             j10=-2
          else
             j100=3
             j10=-3
          endif
        endif

        if(spin.lt..5d0) then
          ifla=j100
          iflb=j10
        else
          ifla=j10
          iflb=j100
        endif

c...Baryons
      elseif(jc4.ne.0) then
c       j1000=  j/1000
c       j100 = (j-j1000*1000)/100
c       j10 = (j-j1000*1000-j100*100)/10

        j1000=jc4
        j100=jc3
        j10=jc2
        if(kf.lt.0) then
          j1000=  -j1000
          j100 = -j100
          j10 = -j10
        endif

c... spin 1/2 baryons
        if(jc1.eq.2) then

        if(spin.lt.vfr13) then
          ifla=j1000
          iflb=ifrkfcpj(j100,j10,0,1.d0)
        elseif(spin.lt.vfr13+vfr14) then
          ifla=j10
          iflb=ifrkfcpj(j1000,j100,0,1.d0)
        elseif(spin.lt.vfr13+vfr14+vfr15) then
          ifla=j100
          s=0.0d0
          if(j1000.eq.j10) s=1.0d0
          iflb=ifrkfcpj(j1000,j10,0,s)
        endif

c...Certain Lambda-like hadrons have two lightest quarks in spin-0:
          if(abs(j100).lt.abs(j10)) then
        if(spin.lt.vfr13) then
          ifla=j1000
          s=0.0d0
          if(j100.eq.j10) s=1.0d0
          iflb=ifrkfcpj(j100,j10,0,s)
        elseif(spin.lt.vfr13+vfr14) then
          ifla=j10
          iflb=ifrkfcpj(j1000,j100,0,1.d0)
        elseif(spin.lt.vfr13+vfr14+vfr15) then
          ifla=j100
          iflb=ifrkfcpj(j1000,j10,0,1.d0)
        endif
          endif

c...spin 3/2 baryons
        else
      if(j1000.eq.j100 .and. j100.eq.j10) then
         ifla=j1000
         iflb=ifrkfcpj(j100,j10,0,1.d0)
      else if(j1000.eq.j100 .and. j100.ne.j10) then
         if(spin.lt.0.3333d0) then
           ifla=j10
           iflb=ifrkfcpj(j1000,j100,0,1.d0)
         else
           ifla=j100
           iflb=ifrkfcpj(j1000,j10,0,1.d0)
         endif
      else if(j1000.eq.j10 .and. j10.ne.j100) then
         if(spin.lt.0.3333d0) then
           ifla=j100
           iflb=ifrkfcpj(j1000,j10,0,1.d0)
         else
           ifla=j10
           iflb=ifrkfcpj(j1000,j100,0,1.d0)
         endif
      else
         if(spin.lt.0.3333d0) then
           ifla=j1000
           iflb=ifrkfcpj(j100,j10,0,1.d0)
         else if(spin.lt.0.6667d0) then
           ifla=j100
           iflb=ifrkfcpj(j1000,j10,0,1.d0)
         else
           ifla=j10
           iflb=ifrkfcpj(j1000,j100,0,1.d0)
         endif
      endif
      endif
        
      else
        call jamerrm(30,0,'(attflv:) Unrecognized particle code KF')
      endif

      if(ifla.eq.0.and.iflb.eq.0) then
        write(check(1),'(''kf='',i9)')kf
        write(check(2),'(''jc5='',i2)')jc5
        write(check(3),'(''jc4='',i2)')jc4
        write(check(4),'(''jc3='',i2)')jc3
        write(check(5),'(''jc2='',i2)')jc2
        write(check(6),'(''jc2='',i2)')jc1
        call jamerrm(30,6,'(attflv:) Something was wrong')
      endif

      if(jamcomp(ifla).eq.0.or.jamcomp(iflb).eq.0) then
        write(check(1),'(''kf ifla iflb='',3i10)')kf,ifla,iflb
        call jamerrm(30,1,'(attflv:)invalid kf code')
      endif

      end

c***********************************************************************

      function ifrkfcpj(ia,ib,ic,s)

c...Purpose: to return the kf code for flavor having ia ib ic.
c...The kf code for a 2- or 3-quark system of spin s composed by 
c...flavor ia, ib, ic: (the system must be qq or qqq, not qqbar, etc).
c...it corresponds to a diquark system if ic=0.........................
      implicit double precision(a-h, o-z)

      ia0 = max( iabs(ia), max(iabs(ib),iabs(ic)))
      ic0 = min( iabs(ia), min(iabs(ib),iabs(ic)))
      ib0 = iabs(ia+ib+ic)-ia0-ic0
      ifrkfcpj = 1000*ia0 + 100*ib0 + 10*ic0 + int(2.d0*(s+0.2d0))+ 1
      if(ia.ne.iabs(ia).or.ib.ne.iabs(ib)) ifrkfcpj = -ifrkfcpj

      end

c***********************************************************************
 
      subroutine kfcnst(kfl10,kfl20,kf,emf) 
 
c...Purpose: to construct kf code from flavour pair
c...kfl1: quark
c...kfl2: anti-/di-quark
c...emf : mass (GeV/c**2)
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
 
c...Default flavour values. input consistency checks. 
      if(kfl10*kfl20.lt.0) then
        if(kfl10.ge.1.and.kfl10.lt.10) then
          kfl1=kfl10
          kfl2=kfl20
        else
          kfl1=kfl20
          kfl2=kfl10
        endif
      else
        if(abs(kfl10).ge.1.and.abs(kfl10).lt.10) then
          kfl1=kfl10
          kfl2=kfl20
        else
          kfl1=kfl20
          kfl2=kfl10
        endif
      endif
      kf1a=iabs(kfl1)
      kf2a=iabs(kfl2)
      kf=0
      if(kf1a.eq.0) return


c...Flavour for meson, possibly with new flavour. 
      if((kf1a.ge.1.and.kf1a.le.10).and.(kf2a.ge.1.and.kf2a.le.10) )then

110     kfs=isign(1,kfl1) 
        kfla=max(kf1a,kf2a) 
        kflb=min(kf1a,kf2a) 
        if(kfla.ne.kf1a) kfs=-kfs 
 
c...Form meson, with spin and flavour mixing for diagonal states. 
        if(kfla.le.2) kmul=int(parj(11)+pjr(0)) 
        if(kfla.eq.3) kmul=int(parj(12)+pjr(0)) 
        if(kfla.ge.4) kmul=int(parj(13)+pjr(0)) 

        if(kmul.eq.0.and.parj(14).gt.0.d0) then 
          if(pjr(0).lt.parj(14)) kmul=2 
        elseif(kmul.eq.1.and.parj(15)+parj(16)+parj(17).gt.0.d0) then 
          rmul=pjr(0) 
          if(rmul.lt.parj(15)) kmul=3 
          if(kmul.eq.1.and.rmul.lt.parj(15)+parj(16)) kmul=4 
          if(kmul.eq.1.and.rmul.lt.parj(15)+parj(16)+parj(17)) kmul=5 
        endif 

        kfls=3 
        if(kmul.eq.0.or.kmul.eq.3) kfls=1 
        if(kmul.eq.5) kfls=5 

        if(kfla.ne.kflb) then 
          kf=(100*kfla+10*kflb+kfls)*kfs*(-1)**kfla 
        else 
          rmix=pjr(0) 
          imix=2*kfla+10*kmul 
          if(kfla.le.3) kf=110*(1+int(rmix+parf(imix-1))+ 
     &    int(rmix+parf(imix)))+kfls 
          if(kfla.ge.4) kf=110*kfla+kfls 
        endif 

        if(kmul.eq.2.or.kmul.eq.3) kf=kf+isign(10000,kf) 
        if(kmul.eq.4) kf=kf+isign(20000,kf) 
 
c...Optional extra suppression of eta and eta'. 
        if(kf.eq.221) then 
          if(pjr(0).gt.parj(25)) goto 110 
        elseif(kf.eq.331) then 
          if(pjr(0).gt.parj(26)) goto 110 
c...eta(1295)
c       else if(kf.eq.20221) then 
c         if(pjr(0).gt.parj(25)) goto 110 
        endif 
 
        if(kf.eq.-211.and.emf.gt.2*pjmass(kf)+0.001d0) kf=-213
        if(kf.eq.211.and.emf.gt.2*pjmass(kf)+0.001d0) kf=213
        if(kf.eq.111.and.emf.gt.2*pjmass(kf)+0.001d0) kf=113

c...Baryon
      else 
 
130       kfla=kf1a 
          kflb=mod(kf2a/1000,10) 
          kflc=mod(kf2a/100,10) 
          kflds=mod(kf2a,10) 
          if(kflb.eq.0.or.kflc.eq.0.or.kflds.eq.0) then
            write(check(1),8000)kfl1,kfl2,emf
 8000       format('kfl1 kfl2 emf',i9,1x,i9,1x,g12.3)
            call jamerrm(1,1,'Error(kfcnst) ikjfl2 should be di-quark')
            return
          endif
 
 
c...SU(6) factors for formation of baryon. try again if fails. 
        kbary=kflds 
        if(kflds.eq.3.and.kflb.ne.kflc) kbary=5 
        if(kfla.ne.kflb.and.kfla.ne.kflc) kbary=kbary+1 

c       wt=parf(60+kbary)+parj(18)*parf(70+kbary) 
c       if(mbary.eq.1.and.mstj(12).ge.2) then 
c         wtdq=pars0 
c         if(max(kflb,kflc).eq.3) wtdq=pars1 
c         if(min(kflb,kflc).eq.3) wtdq=pars2 
c         if(kflds.eq.1) wtdq=wtdq/(3.*par4m) 
c         if(kflds.eq.1) wt=wt*(1.+wtdq)/(1.+parsm/(3.*par4m)) 
c         if(kflds.eq.3) wt=wt*(1.+wtdq)/(1.+parsm) 
c       endif 
c       if(kf2a.eq.0.and.wt.lt.pjr(0)) goto 130 
 
c...Form baryon. distinguish lambda- and sigmalike baryons. 
        kfld=max(kfla,kflb,kflc) 
        kflf=min(kfla,kflb,kflc) 
        kfle=kfla+kflb+kflc-kfld-kflf 
        kfls=2 
        if((parf(60+kbary)+parj(18)*parf(70+kbary))*pjr(0).gt. 
     &  parf(60+kbary)) kfls=4 
        kfll=0 
        if(kfls.eq.2.and.kfld.gt.kfle.and.kfle.gt.kflf) then 
          if(kflds.eq.1.and.kfla.eq.kfld) kfll=1 
          if(kflds.eq.1.and.kfla.ne.kfld) kfll=int(0.25d0+pjr(0)) 
          if(kflds.eq.3.and.kfla.ne.kfld) kfll=int(0.75d0+pjr(0)) 
        endif 
        if(kfll.eq.0) kf=isign(1000*kfld+100*kfle+10*kflf+kfls,kfl1) 
        if(kfll.eq.1) kf=isign(1000*kfld+100*kflf+10*kfle+kfls,kfl1) 
        if(kf.eq.4322.and.emf.le.2.5d0) kf=4232

c...Check mass
        kfa=abs(kf)
        if(kfa.eq.2112.or.kfa.eq.2212) then

          if(emf.ge.parc(28)+parc(29).and.emf.le.1.5d0) then
c            kf=isign(10000+1000*kfld+100*kflf+10*kfle+kfls,kfl1) 
             kf=(kfa+10000)*isign(1,kf)
          else if(emf.le.1.6d0) then
             kf=(kfa+20000)*isign(1,kf)
          else if(emf.le.1.7d0) then
             kf=(kfa+30000)*isign(1,kf)
          else
             kf=(kfa+40000)*isign(1,kf)
          endif

        else if(kfa.eq.1114.or.kfa.eq.2114.or.kfa.eq.2214
     $          .or.kfa.eq.2224) then

           if(emf.le.parc(28)+parc(29)) then
           else if(emf.ge.1.5d0.and.emf.le.1.65d0) then
             kf=(kfa+30000)*isign(1,kf)
           else if(emf.le.1.85d0) then
             kf=(kfa+10000)*isign(1,kf)
           else
             kf=(kfa+20000)*isign(1,kf)
           endif

        else if(kfa.eq.3122) then
           if(emf.ge.
     $         pmas(jamcomp(3122),1)+pmas(jamcomp(211),1)+parc(41)) then
             kf=(kfa+20000)*isign(1,kf)
           endif
        else if(kfa.eq.3112.or.kfa.eq.3212.or.kfa.eq.3222) then
        endif

      endif 

      end 

c***********************************************************************

      function kfprop(kf,k)

c...Purpose: to give charge,baryon number and strangeness
c...of particle PDG(jetset) ID kf
c k=1: 3*charge
c k=2: 3*baryon number
c k=3: strangeness

      implicit double precision(a-h, o-z)
      dimension iq(0:100)
c... d u s c b t l h....
      data iq/0,
     &     -1,2,-1,2,-1,2,-1,2,2*0,
     &     -3,0,-3,0,-3,0,-3,3*0,
     &   3*0,3,9*0,3,2*0,3,63*0/

c...This is a quark
      if(abs(kf).le.10) then
        kfprop=0
        if(k.eq.1) kfprop=iq(abs(kf))*isign(1,kf)
        if(k.eq.2) kfprop=isign(1,kf)
        if(k.eq.3) then
          if(kf.eq.3)  kfprop=-1
          if(kf.eq.-3) kfprop=1
        endif
        return
c...This is a boson
      else if(abs(kf).gt.10.and.abs(kf).le.40) then
        kfprop=0
        return
      endif

      j=abs(kf)
c     jc1=mod(j,10)      ! spin
      jc2=mod(j/10,10)
      jc3=mod(j/100,10)
      jc4=mod(j/1000,10)
c     jc5=mod(j/10000,10)

c...This is a di-quark
      if(j.ge.1000.and.jc2.eq.0) then
        kf1=0
        kf2=kf
        ibary=1
c...This is a hadron
      else 
        call attflv(kf,kf1,kf2)
        if(jc4.eq.0) then
          ibary=0
        else
          ibary=1
        endif
      endif

      j1 = abs(kf1)
      j2 = abs(kf2)

c...Get three times charge
      if(k.eq.1) then
        ich1=iq(j1)*isign(1,kf1)
        if(j2.ge.1000) then
          jc3=mod(j2/100,10)
          jc4=mod(j2/1000,10)
         ich2=(iq(jc3)+iq(jc4))*isign(1,kf2)
        else
         ich2=iq(j2)*isign(1,kf2)
        endif
        kfprop=ich1+ich2

c...Get three times baryon number
      else if(k.eq.2) then
        jj1=0
        if(j1.ge.1) jj1=1
        jj2=1
        if(j2.ge.1000) jj2=2
        kfprop=jj1*isign(1,kf1)+jj2*isign(1,kf2)

c...Get strangeness number
      else if(k.eq.3) then

        kp1=0
        if(ibary.eq.1) then
          if(j1.eq.3) kp1=-1*isign(1,kf1)
          if(j2.ge.1000) then
            jc3=mod(j2/100,10)
            jc4=mod(j2/1000,10)
            kp2=0
            if(jc3.eq.3) kp2=kp2+1
            if(jc4.eq.3) kp2=kp2+1
            kp2=-kp2*isign(1,kf2)
          endif
          kfprop=kp1+kp2
        else
          kp1=0
          if(jc3.eq.5) then
               if(jc2.eq.3) kp1=-1
          elseif((jc2.eq.3.and.jc3.ne.3).or.(jc2.ne.3.and.jc3.eq.3))then
               kp1=1
          endif
          if(kf.le.0) kp1=-kp1
          kfprop=kp1
        endif

      endif

      end

c***********************************************************************

      subroutine attflv2(kf,ifla,iflb,iflc)

c...Purpse: to give quarkflavour contents.
c...For mesons, the order of the end flavors is randomly given;

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      integer ifrkfcpj
      data vfr13,vfr14,vfr15/0.167d0,  .333d0,  .5d0/
c...    x quark
c... xx0x di-quark
c...   1x leptons
c...   xx gauge and Higgs bosons
c...  xxx meson
c...x0xxx meson
c... xxxx baryon
c...xxxxx baryon

      j = abs(kf)
      jc1=mod(j,10)      ! spin
      jc2=mod(j/10,10)
      jc3=mod(j/100,10)
      jc4=mod(j/1000,10)
      jc5=mod(j/10000,10)
      isg=isign(1,kf)

      ifla=0
      iflb=0
      iflc=0
c...Quark
      if(j.lt.10) then
         ifla=kf
         return
c...Lepton in idtrns
      else if(j.lt.20) then
         return
c...Gauge or higgs boson
      else if(j.le.100) then
         return
      endif

c...Di-quark
      if(jc2.eq.0) then
         ifla=jc3*isg
         iflb=jc4*isg
         return
      endif

      id1=isign(1,kf)*(jc4*100+jc3*10+jc2)
      spin=pjr(0)

c..Mesons
      if(jc4.eq.0) then
c...Identify the quark and antiquark in mesons:
c       j100=  j/100
c       j10 = (j-j100*100)/10
        j100=jc3
        j10=jc2
        isgn = (-1)**max(j100, j10)
        if(kf.lt.0) isgn = -isgn 
        if(isgn.gt.0) j10 = -j10
        if(isgn.lt.0) j100 = -j100

        if(id1.eq.11) then  ! pi0 eta rho0 omega...
          ranx=pjr(0)
          if(ranx.lt.0.5d0) then
             j100=1
             j10=-1
          else
             j100=2
             j10=-2
          endif
        endif

        if(j.eq.331) then  ! eta'  (1/6 1/6 2/3)
          ranx=pjr(0)
          if(ranx.lt.0.1667d0) then
             j100=1
             j10=-1
          else if(ranx.lt.0.3333d0) then
             j100=2
             j10=-2
          else
             j100=3
             j10=-3
          endif
        endif

        if(spin.lt..5d0) then
          ifla=j100
          iflb=j10
        else
          ifla=j10
          iflb=j100
        endif

c...Baryons
      elseif(jc4.ne.0) then

        ifla=jc4*isg
        iflb=jc3*isg
        iflc=jc2*isg
        
      else
        call jamerrm(30,0,'(attflv:) Unrecognized particle code KF')
      endif

      if(ifla.eq.0.and.iflb.eq.0) then
        write(check(1),'(''kf='',i9)')kf
        write(check(2),'(''jc5='',i2)')jc5
        write(check(3),'(''jc4='',i2)')jc4
        write(check(4),'(''jc3='',i2)')jc3
        write(check(5),'(''jc2='',i2)')jc2
        write(check(6),'(''jc2='',i2)')jc1
        call jamerrm(30,6,'(attflv2:) Something was wrong')
      endif

c     if(jamcomp(ifla).eq.0.or.jamcomp(iflb).eq.0) then
c       write(check(1),'(''kf ifla iflb='',3i10)')kf,ifla,iflb
c       call jamerrm(30,1,'(attflv2:)invalid kf code')
c     endif

      end

C***********************************************************************

      function jamflav(kf,msel)

c...Give baron number or flavour content of particles.
c...msel=1: net d
c...msel=2: net u
c...msel=3: net s
c...msel=4: net c
c...msel=5: net b
c...msel=6: net t
c...msel=11: baryon number

      dimension kfl(3)

      jamflav=0
      kfs=isign(1,kf)
      kfa=iabs(kf)
      kfl(1)=mod(kfa/1000,10)
      kfl(2)=mod(kfa/100,10)
      kfl(3)=mod(kfa/10,10)

      if(kfa.le.10) then
        if(msel.le.10) then
           if(msel.eq.kfa) jamflav=kfs
        elseif(msel.eq.11) then
           jamflav=kfs
        endif
      else if(kfa.le.100) then

        jamflav=0

c....Diquark.
      else if(kfl(3).eq.0) then
        if(msel.le.10) then
          if(kfl(1).eq.msel) jamflav=kfs
          if(kfl(2).eq.msel) jamflav=jamflav+kfs
        else if(msel.eq.11) then
          jamflav=2*kfs
        endif

c....Mesons.
      else if(kfl(1).eq.0) then
c       kfma=kfl(2)*10+kfl(3)
c       if(kfma.eq.11.and.pjr(0).gt.0.5d0)then
c         kfl(2)=2
c         kfl(3)=2
c       else if(kfma.eq.22.and.pjr(0).gt.0.5)then
c         kfl(2)=1
c         kfl(3)=1
c       endif
        if(msel.le.10)then
          kfl(2)=(kfl(2)*(-1)**kfl(2))*kfs
          kfl(3)=(-kfl(3)*(-1)**iabs(kfl(2)))*kfs
          if(abs(kfl(2)).eq.msel) jamflav=isign(1,kfl(2))
          if(abs(kfl(3)).eq.msel) jamflav=jamflav+isign(1,kfl(3))
        else if(msel.eq.11) then
          jamflav=0
        endif

c....Baryons.
      else
        if(msel.le.10) then
          if(kfl(1).eq.msel) jamflav=kfs
          if(kfl(2).eq.msel) jamflav=jamflav+kfs
          if(kfl(3).eq.msel) jamflav=jamflav+kfs
        else if(msel.eq.11) then
          jamflav=3*kfs
        endif
      endif
 
      end

