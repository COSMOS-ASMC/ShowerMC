c***********************************************************************
c                                                                      *
c        PART 5: Particle decay part                                   *
c                                                                      *
c   List of subprograms in rough order of relevance with main purpose  *
c      (S = subroutine, F = function, B = block data, E = entry)       *
c                                                                      *
c                                                                      *
c s jamdec   to administrate the fragmentation of jet system or decay  *
c s jamrdec  to calculate decay of hadrons                             *
c s jamwidm  to calculate momentum dependent partical decay width      *
c f jamdtim  to calculate life time of the unstable particle           *
c s jambwmas to generate mass according to the B.W. distribution.      *
c f jamdlwid to calculate d(1232) momentum dependent total decay width *
c s ludecy   to do the decay of a particle                             *
c s jamfdec  to decay unstable particles at the end of simulation      *
c s jamsetd  to set decay switch of the particles.                     *
c                                                                      *
c***********************************************************************

      subroutine jamdec(ip,indd,nadd,icon)

c...Purpose: to administrate the fragmentation of jet system or decay
c...of unstable particle.

      include 'jam1.inc'
      include 'jam2.inc'
     
c...PYTHIA common block
      common/jyjets/njet,npad,kjet(1000,5),pjet(1000,5),vjet(1000,5)
      save  /jyjets/
      common/daxis/vbex,vbey,vbez,vgam,vphi,vthe
      real*8 jamdtim
      dimension pd(5),pd1(5),indd(100)
      character*16 charn


      icon=0
      nadd=0
      ks0=mod(abs(k(1,ip)),10)
      kc0=mste(22)
      nmes0=nmeson
      nv0=nv
      nbary0=nbary

c...Jet decay.
      if(ks0.eq.3.or.k(1,ip).eq.4) then
         call jamjdec(ip,indd,nadd,icon)
         jetd=2

c...Resonance decay.
      else if(k(1,ip).eq.2.or.mdcy(kc0,1).eq.1) then

          id0=kchg(kc0,5)
          jetd=3

          k01=k(1,ip)
          kf0=k(2,ip)
          k05=k(5,ip)
          k06=k(6,ip)

          v05=v(5,ip)
          x01=r(1,ip)
          x02=r(2,ip)
          x03=r(3,ip)
          emd=p(5,ip)

          do j=1,5
          pd(j)=p(j,ip)
          pd1(j)=p(j,ip)
          end do
          bex=pd(1)/pd(4)
          bey=pd(2)/pd(4)
          bez=pd(3)/pd(4)
          gam=pd(4)/pd(5)
          ires=0
          if(kq(1,ip).eq.999999) then
            ires=1
            vbex=vq(1,ip)
            vbey=vq(2,ip)
            vbez=vq(3,ip)
            vgam=vq(4,ip)
            vphi=vq(5,ip)
            vthe=vq(6,ip)
c...[ AO, 050614
          else
            vbex=0.0d0
            vbey=0.0d0
            vbez=0.0d0
            vgam=0.0d0
            vphi=0.0d0
            vthe=0.0d0
c...] AO, 050614
          endif

          iang=0
c.....2010/6/29 pjdec() does not work for X(1690) 
c... because invariant  mass of some decay chanel chennel may be
c...less than than the mass of daughter.
c         if(abs(kf0).eq.13322.or.abs(kf0).eq.13312) iang=1
          ibar=k(9,ip)
          if(abs(ibar).eq.3) iang=1
          if(id0.eq.id_charmb) iang=0
          if(id0.eq.id_bottb) iang=0
          if(id0.eq.id_omega) iang=0
          if(abs(kf0).eq.3122) iang=0 ! Lambada
          if(abs(kf0).eq.3112) iang=0 ! Sigma-
          if(abs(kf0).eq.3212) iang=0 ! Sigma0
          if(abs(kf0).eq.3222) iang=0 ! Sigma+
          if(abs(kf0).eq.3312) iang=0 ! Xi-
          if(abs(kf0).eq.3322) iang=0 ! Xi0
          if(abs(kf0).eq.3314) iang=0 ! Xi*-
          if(abs(kf0).eq.3324) iang=0 ! Xi*0

          iswang=mstc(61)
c         if(k(1,ip).eq.2.and.abs(k(9,ip)).eq.3) then
          if(k(1,ip).eq.2.and.kq(1,ip).eq.999999) then  ! resonance
            iang=1                                      ! isotropic
            if(emd.gt.1.1d0)  iang=2                    ! anisotropic
            if(id0.eq.id_delt.and.emd.gt.1.1d0) then
                if(mstc(61).ne.0) iswang=5
            endif
          endif


          if(mdcy(kc0,1).eq.0) then
            write(check(1),'(''kc0='',i3)')kc0
            write(check(2),'(3i10)')mste(21),mste(22),kcp(1,1)
            write(check(3),'(''ip k1 kf'',3i10)')ip,k01,kf0
            call jamerrm(1,3,'(jamdec:) invalid k(1,ip)???')
            k(1,ip)=1
            return
          endif

          njet=1
          kjet(1,1)=1
          kjet(1,2)=kf0
          kjet(1,3)=0
          kjet(1,4)=0
          kjet(1,5)=0
          pjet(1,1)=pd1(1)
          pjet(1,2)=pd1(2)
          pjet(1,3)=pd1(3)
          pjet(1,4)=pd1(4)
          pjet(1,5)=pd1(5)

          if(iang.ge.1) then
            call jamrdec(icon,iang,iswang,parc(43))
          else
            call pjdecy(1,icon)
          endif

c....All decay channel closed. This may be that mass is too small.
          if(icon.ne.0) then
            print *,'(decay:)icon.ne.0 after jamrdec',icon,kf0,emd
            k(1,ip)=1
            v(5,ip)=1d+25
            return
          endif

c......No decay???
          if(njet.eq.1) then
            call pjname(kf0,charn)
            write(check(1),8000)icon,ip,kf0,p(5,ip),charn
 8000       format('icon',i4,'ip',i9,'kf0=',i9,'p5',g12.3,a16)
            call pjlist(1)
            call jamerrm(3,1,'(jamdec:)after PYdecy njet=1')
            k(1,ip)=1
            return
          endif

c...Remove decayed particles.
          call pjedit(2)

c...Update particle after decay
          nadd=0
          ikaon=0
c----------------------------------------------------------------
          do i=1,njet  ! **** Loop over produced particles.
c----------------------------------------------------------------

           nadd=nadd+1
           if(i.eq.1) then
             indx=ip
             improd=0
           else
             improd=1
             nmeson=nmeson+1
             nv=nv+1
             if(nv.gt.mxv) then
               call jamerrm(30,0,'(jamdec:) particle too large'
     $                       //'mxv should be changed')
             endif
             indx=nv
           endif

c...Convert Ks, KL into k0, k0bar
           if(mste(40).eq.0) then
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
           endif
c...Update particle array.
           indd(nadd)=indx

           kf1=kjet(i,2)
           kc1=jamcomp(kf1)
           call jamkupda(2,indx,kf1,kc1,1,k06,k05,0)


           if(iang.eq.2.and.ires.eq.1) then
             kq(1,indx)=999999
             kq(2,indx)=0
             vq(1,indx)=vbex
             vq(2,indx)=vbey
             vq(3,indx)=vbez
             vq(4,indx)=vgam
             vq(5,indx)=vphi
             vq(6,indx)=vthe
             do l=7,10
               vq(l,indx)=0.0d0
             end do
           else
             kq(1,indx)=0
             kq(2,indx)=0
             do l=1,10
               vq(l,indx)=0.0d0
             end do
           endif

           p(1,indx)=pjet(i,1)
           p(2,indx)=pjet(i,2)
           p(3,indx)=pjet(i,3)
           p(4,indx)=pjet(i,4)
           p(5,indx)=pjet(i,5)

c.....Add small distance from the mother.
           if(improd.eq.1) then
             deltx=parc(42)*sqrt(rn(0))
             cos1=1.d0-2.d0*rn(0)
             sin1=sqrt(1.d0-cos1**2)
             phi1=2*paru(1)*rn(0)
             dxr=deltx*sin1*cos(phi1)
             dyr=deltx*sin1*sin(phi1)
             dzr=deltx*cos1
             dtr=0.0d0
             call jamrobo(0.0d0,0.0d0,bex,bey,bez,gam,dxr,dyr,dzr,dtr)
             r(1,indx)=x01+dxr
             r(2,indx)=x02+dyr
             r(3,indx)=x03+dzr
             r(4,indx)=v05+dtr
             r(5,indx)=r(4,indx)
           else
             r(1,indx)=x01
             r(2,indx)=x02
             r(3,indx)=x03
             r(4,indx)=v05
c...Formation time.
             r(5,indx)=r(4,indx)
           endif

c...Vertex point.
           v(1,indx)=r(1,indx)
           v(2,indx)=r(2,indx)
           v(3,indx)=r(3,indx)
           v(4,indx)=r(4,indx)

c...Life time.
           v(5,indx)=r(5,indx)
     $                +jamdtim(1,kf1,kc1,k01,p(5,indx),p(4,indx))

c----------------------------------------------------------------
          end do   ! **** End loop over produced particles.
c----------------------------------------------------------------

c...qmd:Recalculate momentum to recover total energy.
cq         if(mstc(57).ge.1) then
c......Sorry now only two-body decay can work.
cq         if(nadd.eq.2) then
cq           i1=indd(1)
cq           i2=indd(2)
cq           gamma=pd(4)/pd(5)
cq           pcs=pd(1)*p(1,i1)+pd(2)*p(2,i1)+pd(3)*p(3,i1)
cq           transf=(pcs/(pd(4)+pd(5))-p(4,i1))/pd(5)
cq           pxr=p(1,i1)+pd(1)*transf
cq           pyr=p(2,i1)+pd(2)*transf
cq           pzr=p(3,i1)+pd(3)*transf
cq           pr=sqrt(pxr**2+pyr**2+pzr**2)
cq           call setmom(2,pare(11),pd,gamma,p(5,i1),p(5,i2),
cq   $                                   pr,pxr,pyr,pzr,i1,i2,icon )
cq           if(icon.ne.0) then
cq             call jamsave(2,1,ip)
cq             nmeson=nmes0
cq             nv=nv0
cq             call caldis2(ip,ip)
cq             icon=3
cq             if(mstc(8).ge.1) 
cq   $             call jamerrm(1,0,'(jamdec:) decay can not recover'
cq   $                     //' energy')
cq             return
cq           end if
cq         endif
cq         endif

c...Check Pauli blocking.
       if(mstc(56).ge.1) then
         do j=1,nadd
           i=indd(j)
           if(k(2,i).eq.2212.or.k(2,i).eq.2112) then
             call jampauli(i,ntag,phase)
             if(ntag.eq.1) then
               call jamsave(2,1,ip)
               nmeson=nmes0
               nbary=nbary0
               nv=nv0
               icon=5
               mstd(51)=mstd(51)+1
               v(5,ip)=r(5,ip)
     $              +jamdtim(1,k(2,ip),mste(22),k(1,ip),p(5,ip),p(4,ip))
               return
             endif
           endif
         end do
       endif

      else

        write(check(1),8100)ip,k(1,ip),k(2,ip)
        write(check(2),8110)v(5,ip),p(5,ip)
 8100   format('ip k1 k2',i10,i4,i9)
 8110   format('v5 p5',2(g12.3,1x))
        call jamerrm(1,2,'(jamdec:) Invalid ip k(1,ip) k1 k2 v5 em')
        kc=jamcomp(k(2,ip))
        if(kc.eq.0) then
           v(5,ip)=1.d+27
           write(3,*)'(jamdec:) kc=0???'
        else
          if(pmas(kc,2).le.1d-7.or.mdcy(kc,1).eq.0
     $              .or.mdcy(kc,2).eq.0.or.mdcy(kc,3).eq.0)then
            v(5,ip)=r(5,ip)+jamdtim(1,k(2,ip),kc,2,p(5,ip),p(4,ip))
          else
            v(5,ip)=1.d+28
          endif
        endif

        return
      endif
c...end decay

c...Count number of decay and collision history.
      if(icon.eq.0) then
        mstd(50)=mstd(50)+1
        if(mstc(162).eq.1) call jamclhst(jetd,5)
      endif

      end
 
c***********************************************************************
 
      subroutine jamrdec(icon,iang,iswang1,ptx0) 
 
c...Purpose: to handle the decay of unstable particles. 
      implicit double precision(a-h, o-z)
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200) 
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/jydat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
      save /jyjets/,/jydat1/,/jydat2/,/jydat3/ 
      common/daxis/vbex,vbey,vbez,vgam,vphi,vthe

      dimension pd(5),pv(10,5),rord(10),ue(3),be(3),wtcor(10)
      parameter(maxbr=70)
      dimension ibranch(maxbr),pwid(maxbr)
      character chekc*70
      parameter( pi=3.141593d0 )
      parameter(utrat=0.00d0)
      data wtcor/2.d0,5.d0,15.d0,60.d0,250.d0,1500.d0,1.2d4,1.2d5, 
     & 150.d0,16.d0/ 
 
c...Functions: momentum in two-particle decays, four-product and 
c...Matrix element times phase space in weak decays. 
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a) 
      four(i,j)=p(i,4)*p(j,4)-p(i,1)*p(j,1)-p(i,2)*p(j,2)-p(i,3)*p(j,3) 
 
c...Initial values. 
      icon=0
      ntry=0 
      nsav=n 
      kfa=iabs(k(n,2)) 
      kfs=isign(1,k(n,2)) 
      kc=jamcomp(kfa) 
      emdec=p(n,5)
      do j=1,5
       pd(j)=p(n,j)
      end do
 
 
c...B-B~ mixing: flip sign of meson appropriately. 
      mmix=0 
      if((kfa.eq.511.or.kfa.eq.531).and.mstj(26).ge.1) then 
        xbbmix=parj(76) 
        if(kfa.eq.531) xbbmix=parj(77) 
        vip5=-pmas(kc,4)*log(pjr(0))
        if(sin(0.5d0*xbbmix*vip5/pmas(kc,4))**2.gt.pjr(0)) mmix=1 
        if(mmix.eq.1) kfs=-kfs 
      endif 
 
c...Check existence of decay channels. Particle/antiparticle rules. 
      kca=kc 
      if(mdcy(kc,2).gt.0) then 
        mdmdcy=mdme(mdcy(kc,2),2) 
        if(mdmdcy.gt.80.and.mdmdcy.le.90) kca=mdmdcy 
      endif 
      if(mdcy(kca,2).le.0.or.mdcy(kca,3).le.0) then 
        write(chekc,'(i9,1x,i3,1x,g12.3)')k(n,2),kca,emdec
        call jamerrm(9,0,'(rdecy:) no decay channel defined kf kc em='
     $              //chekc) 
        icon=1
        return 
      endif 

      if(mod(kfa/1000,10).eq.0.and.(kca.eq.85.or.kca.eq.87)) kfs=-kfs 

      if(kchg(kc,3).eq.0) then 
        kfsp=1 
        kfsn=0 
        if(pjr(0).gt.0.5d0) kfs=-kfs 
      elseif(kfs.gt.0) then 
        kfsp=1 
        kfsn=0 
      else 
        kfsp=0 
        kfsn=1 
      endif 
 
c...Sum branching ratios of allowed decay channels. 
      call jamwidm(kca,kfsp,kfsn,0,0,emdec,ibranch,pwid,brsu,itag)

  235  if(brsu.le.0.0d0)then
        write(mstu(11),*)'================================'
        call pjerrm(2,'(jamrdec:) All decay channels closed') 
        write(mstu(11),*)'totwid=0',k(n,2),emdec
        write(mstu(11),*)'================================'
        icon=120
        return
       endif
 
c...Select decay channel among allowed ones. 
  240 rbr=brsu*pjr(0) 
      idl=mdcy(kca,2)-1 
      maxb=mdcy(kca,3)
      do ibra=1,maxb
        idl=idl+1
        ibrac=ibra
        rbr=rbr-pwid(ibra)
        if(rbr.le.0.0d0) goto 255
      end do
255   continue
      idc=idl
 
c...Start readout of decay channel: matrix element, reset counters. 
      mmat=mdme(idc,2) 
      if(mmat.ge.11) then
        icon=11
        return
      endif

  260 ntry=ntry+1 
      if(ntry.gt.1000) then 
        call jamerrm(4,0,'(rdecy:) caught in infinite loop') 
        icon=140
        return
      endif 

      i=n 
      np=0 
      nq=0 
      mbst=0 

      do 270 j=1,5
       pv(1,j)=p(nsav,j) 
  270 continue 

      ps=0.d0 
c...Read out decay products. convert to standard flavour code. 
      jtmax=5 
      if(mdme(idc+1,2).eq.101) jtmax=10 

      do 280 jt=1,jtmax 

        if(jt.le.5) kp=kfdp(idc,jt) 
        if(jt.ge.6) kp=kfdp(idc+1,jt-5) 
        if(kp.eq.0) goto 280 
        kpa=iabs(kp) 
        kcp=jamcomp(kpa) 

c....This particle is its own antiparticle.
        if(kchg(kcp,3).eq.0.and.kpa.ne.81.and.kpa.ne.82) then 
          kfp=kp 
        elseif(kpa.ne.81.and.kpa.ne.82) then 
          kfp=kfs*kp 
        else
          icon=81
          return
        endif

c...Add decay product to event record or to quark flavour list. 
        i=i+1 
        np=np+1 

        k(i,1)=1+mod(nq,2) 
        k(i,2)=kfp 
        k(i,3)=nsav
        k(i,4)=0 
        k(i,5)=0 
        if(pmas(kcp,2).le.1.d-7) then
          p(i,5)=pjmass(kfp) 
        else
          emmin=pmas(kcp,1)-pmas(kcp,3)+parj(64)
          emmax=emdec-ps-parj(64)
          do j1=jt+1,jtmax 
             kp=kfdp(idc,j1) 
             if(kp.ne.0) then
               kpc=jamcomp(kp)
               emmax=emmax-pmas(kpc,1)-pmas(kpc,3)
             endif
          end do
          if(emmax.lt.emmin) then
            write(mstu(11),*)'(jamrdec:) emax<emin',emmin,emmax
            ibranch(ibrac)=0
            brsu=brsu-pwid(ibrac)
            pwid(ibrac)=0.0d0
            if(brsu.le.0.0d0) then
              icon=22
              return
            endif
            goto 235
          endif
          if(emmax.le.emmin) then
            p(i,5)=emmin
          else
            call jambwmas(emmin,emmax,pmas(kcp,1),pmas(kcp,2),p(i,5)
     $       ,icon)
          endif
        endif
        ps=ps+p(i,5) 
  280 continue 
 
c...Fully specified final state: check mass broadening effects. 
      if(np.ge.2.and.ps+parj(64).gt.pv(1,5)) then
        goto 260
      endif
      nd=np 
      if(np.le.1) then
        write(mstu(11),*)'nasv emd=',nsav,k(nsav,2),emdec
        icon=10
        return
      endif
 
 
c...Determine position of grandmother, number of sisters, q -> w sign. 
      nm=0 
      kfas=0 
      msgn=0 
      if(mmat.eq.3.or.mmat.eq.46) then 
        im=k(nsav,3) 
        if(im.lt.0.or.im.ge.nsav) im=0 
        if(mmat.eq.46.and.mstj(27).eq.1) then 
          im=0 
        elseif(mmat.eq.46.and.mstj(27).ge.2.and.im.ne.0) then 
          if(k(im,2).eq.94) then 
            im=k(k(im,3),3) 
            if(im.lt.0.or.im.ge.nsav) im=0 
          endif 
        endif 
        if(im.ne.0) kfam=iabs(k(im,2)) 
        if(im.ne.0.and.mmat.eq.3) then 
          do 390 il=max(nsav-2,im+1),min(nsav+2,n) 
          if(k(il,3).eq.im) nm=nm+1 
          if(k(il,3).eq.im.and.il.ne.nsav) isis=il 
  390     continue 
          if(nm.ne.2.or.kfam.le.100.or.mod(kfam,10).ne.1.or. 
     &    mod(kfam/1000,10).ne.0) nm=0 
          if(nm.eq.2) then 
            kfas=iabs(k(isis,2)) 
            if((kfas.le.100.or.mod(kfas,10).ne.1.or. 
     &      mod(kfas/1000,10).ne.0).and.kfas.ne.22) nm=0 
          endif 
        elseif(im.ne.0.and.mmat.eq.46) then 
          msgn=isign(1,k(im,2)*k(nsav,2)) 
          if(kfam.gt.100.and.mod(kfam/1000,10).eq.0) msgn= 
     &    msgn*(-1)**mod(kfam/100,10) 
        endif 
      endif 
 
c...Kinematics of one-particle decays. 
      if(nd.eq.1) then 
        do 400 j=1,4 
        p(n+1,j)=p(nsav,j) 
  400   continue 
        goto 660 
cvvvv
c...[
      else if(nd.eq.2.and.iswang1.ne.0) then

       i1=n+1
       i2=n+2
       em1=p(i1,5)
       em2=p(i2,5)
       prsq=(pd(5)*pd(5)-(em1+em2)**2)*(pd(5)*pd(5)-(em1-em2)**2)
       pr=sqrt(prsq)/(2*pd(5))

c...[
       if(iang.ne.2) then       ! isotropic decay
          cos1=1.0d0-2*rn(0)
          sin1=sqrt(1.0d0-cos1**2)
          phi1=paru(2)*rn(0)
          prx=pr*sin1*cos(phi1)
          pry=pr*sin1*sin(phi1)
          prz=pr*cos1
c...[
       else     ! iang=2, anisotropic decay

c...[ iswang1=5
        if(iswang1.eq.5) then   ! fixed L(=1) resonance decay
         itry=0
   30    cos1=-1d0+2d0*rn(0)
         itry=itry+1
         if(itry.le.200) then
           if(rn(0).gt.(1d0+cos1**2)/2d0) goto 30
         else
           call jamerrm(1,0,'(jamrdec:)delta(1232)itry>200')
         endif
         sin1=sqrt(1d0-cos1**2)
         pt=pr*sin1
         prz=pr*cos1

c...[ iswang1=1
c...Cut Gauss
       else if(iswang1.eq.1) then
         pt=ptx0*sqrt(-log(1.d0-rn(0)*(1.d0-exp(-pr*pr/ptx0/ptx0))))
         prz=sqrt(pr**2-pt**2)
c... u-t mixing ratio
         if(rn(0).lt.utrat) prz=-prz

c...[ iswang1=2
c...Gauss + Isotropic
       else if(iswang1.eq.2) then
         ptx=min(pr,ptx0)
         pt=ptx0*sqrt(-log(max(1.d-10,rn(0))))
         if(pt.gt.pr-0.001d0) then
           prz=pr*(1.0d0-2*rn(0))
           pt=sqrt(pr*pr-prz*prz)
         else
           prz=sqrt(pr**2-pt**2)
c... u-t mixing ratio
         if(rn(0).lt.utrat) prz=-prz
         endif

c...[ iswang1=3
c...Modified version 7: pt --> pr*theta1
c...See the comments in coll1.f (s angrr)
        else if(iswang1.eq.3) then
         ptx=pr*pi/2
         expf1=1.0d0-exp(-ptx*ptx/ptx0/ptx0)
         pt=ptx0*sqrt(-log(1.d0-rn(0)*expf1))
         theta1=pt/pr
         pt=pr*sin(theta1)
         prz=pr*cos(theta1)
c... u-t mixing ratio
         if(rn(0).lt.utrat) prz=-prz

c...[ iswang1=4
c...Modified version ?:
       else if(iswang1.eq.4) then
         ptx=min(pr,ptx0)
         pt=ptx*sqrt(-log(1.d0-rn(0)*(1.d0-exp(-pr*pr/ptx/ptx))))
         prz=sqrt(pr**2-pt**2)
c... u-t mixing ratio
         if(rn(0).lt.utrat) prz=-prz

c...[ iswang1=0, or >=5
c...Isotropic
        else
          cos1=1.0d0-2*rn(0)
          sin1=sqrt(1.0d0-cos1**2)
          pt=pr*sin1
          prz=pr*cos1
        endif
c...]]]]]]
c...
        phi1=paru(2)*rn(0)
        prx=pt*cos(phi1)
        pry=pt*sin(phi1)
        prs=prx*prx+pry*pry+prz*prz
        ee1=sqrt(em1*em1+prs)

c...[ Decay direction: AO, 050615
        pdd1=pd(1)
        pdd2=pd(2)
        pdd3=pd(3)
        pdd4=pd(4)
        isw_axis=1

c isw_axis=0: previous projectile direction is favored
        if(isw_axis.eq.0) then
c         call jamrobo(0d0,0d0,-vbex,-vbey,-vbez,vgam
c    &                                          ,pdd1,pdd2,pdd3,pdd4)
c         call jamrobo(0d0,-vphi,0d0,0d0,0d0,1d0,pdd1,pdd2,pdd3,pdd4)
c         call jamrobo(-vthe,0d0,0d0,0d0,0d0,1d0,pdd1,pdd2,pdd3,pdd4)
c         call jamrobo(-vthe,-vphi,0d0,-vbex,-vbey,-vbez,vgam
          call jamrobo(-vthe,-vphi,-vbex,-vbey,-vbez,vgam
     &                                          ,pdd1,pdd2,pdd3,pdd4)
          if(pdd3.lt.0) prz=-prz ! resonance is moving backward in the previous collision frame
          call jamrobo(vthe,0d0,0d0,0d0,0d0,1d0,prx,pry,prz,ee1)
          call jamrobo(0d0,vphi,0d0,0d0,0d0,1d0,prx,pry,prz,ee1)
c isw_axis=1: resonance direction in previous collision frame is favored
        else
c...[ AO, 050614 (comment):
c     These are angle of pdd in the previous
c     collision frame, where the previous leading particle have
c     momentum in the z-direction
          call jamrobo(0d0,0d0,-vbex,-vbey,-vbez,vgam
     &                ,pdd1,pdd2,pdd3,pdd4)
          phi=pjangl(pdd1,pdd2)
          the=pjangl(pdd3,sqrt(pdd1**2+pdd2**2))
c...]
          call jamrobo(the,phi,0d0,0d0,0d0,1d0,prx,pry,prz,ee1)
        endif
c...] Decay Direction: AO, 050614

        endif
c...(iang=2)]]

       pcs=pd(1)*prx+pd(2)*pry+pd(3)*prz
       ecm1=sqrt(em1**2+prx**2+pry**2+prz**2)
       transf=(pcs/(pd(4)+pd(5))+ecm1)/pd(5)
       p(i1,1)=prx+pd(1)*transf
       p(i1,2)=pry+pd(2)*transf
       p(i1,3)=prz+pd(3)*transf
       p(i1,4)=sqrt(em1**2+p(i1,1)**2+p(i1,2)**2+p(i1,3)**2)
       ecm2=sqrt(em2**2+prx**2+pry**2+prz**2)
       transf=(-pcs/(pd(4)+pd(5))+ecm2)/pd(5)
       p(i2,1)=-prx+pd(1)*transf
       p(i2,2)=-pry+pd(2)*transf
       p(i2,3)=-prz+pd(3)*transf
       p(i2,4)=sqrt(em2**2+p(i2,1)**2+p(i2,2)**2+p(i2,3)**2)

c       if(iang.eq.2) then
c        call jamrobo(the,phi,vbx,vby,vbz,vbg
c     &     ,p(i1,1),p(i1,2),p(i1,3),p(i1,4))
c        call jamrobo(the,phi,vbx,vby,vbz,vbg
c     &     ,p(i2,1),p(i2,2),p(i2,3),p(i2,4))
c       write(*,*) 'p(i1)',(p(i1,l)/p(i1,4),l=1,3)
c       endif

       n=n+2
       goto 800
      endif 
c...]
 
c...Calculate maximum weight nd-particle decay. 
      pv(nd,5)=p(n+nd,5) 
      if(nd.ge.3) then 
        wtmax=1.d0/wtcor(nd-2) 
        pmax=pv(1,5)-ps+p(n+nd,5) 
        pmin=0.d0 
        do 410 il=nd-1,1,-1 
        pmax=pmax+p(n+il,5) 
        pmin=pmin+p(n+il+1,5) 
        wtmax=wtmax*pawt(pmax,pmin,p(n+il,5)) 
  410   continue 
      endif 
 
c...Find virtual gamma mass in Dalitz decay. 
  420 if(nd.eq.2) then 
      elseif(mmat.eq.2) then 
        pmes=4.d0*pmas(11,1)**2 
        pmrho2=pmas(131,1)**2 
        pgrho2=pmas(131,2)**2 
  430   pmst=pmes*(p(nsav,5)**2/pmes)**pjr(0) 
        wt=(1+0.5d0*pmes/pmst)*sqrt(max(0.d0,1.d0-pmes/pmst))* 
     &  (1.d0-pmst/p(nsav,5)**2)**3*(1.d0+pgrho2/pmrho2)/ 
     &  ((1.d0-pmst/pmrho2)**2+pgrho2/pmrho2) 
        if(wt.lt.pjr(0)) goto 430 
        pv(2,5)=max(2.00001d0*pmas(11,1),sqrt(pmst)) 
 
c...M-generator gives weight. If rejected, try again. 
      else 
  440   rord(1)=1.d0 
        do 470 il1=2,nd-1 
        rsav=pjr(0) 
        do 450 il2=il1-1,1,-1 
        if(rsav.le.rord(il2)) goto 460 
        rord(il2+1)=rord(il2) 
  450   continue 
  460   rord(il2+1)=rsav 
  470   continue 
        rord(nd)=0.d0 
        wt=1.d0 
        do 480 il=nd-1,1,-1 
        pv(il,5)=pv(il+1,5)+p(n+il,5)+(rord(il)-rord(il+1))*(pv(1,5)-ps)
        wt=wt*pawt(pv(il,5),pv(il+1,5),p(n+il,5)) 
  480   continue 
        if(wt.lt.pjr(0)*wtmax) goto 440 
      endif 
 
c...Perform two-particle decays in respective CM frame. 
  490 do 510 il=1,nd-1 
      pa=pawt(pv(il,5),pv(il+1,5),p(n+il,5)) 
      ue(3)=2.d0*pjr(0)-1.d0 
      phi=paru(2)*pjr(0) 
      ue(1)=sqrt(1.d0-ue(3)**2)*cos(phi) 
      ue(2)=sqrt(1.d0-ue(3)**2)*sin(phi) 
      do 500 j=1,3 
      p(n+il,j)=pa*ue(j) 
      pv(il+1,j)=-pa*ue(j) 
  500 continue 
      p(n+il,4)=sqrt(pa**2+p(n+il,5)**2) 
      pv(il+1,4)=sqrt(pa**2+pv(il+1,5)**2) 
  510 continue 
 
c...Lorentz transform decay products to lab frame. 
      do 520 j=1,4 
      p(n+nd,j)=pv(nd,j) 
  520 continue 
      do 560 il=nd-1,1,-1 
      do 530 j=1,3 
      be(j)=pv(il,j)/pv(il,4) 
  530 continue 
      ga=pv(il,4)/pv(il,5) 
      do 550 i=n+il,n+nd 
      bep=be(1)*p(i,1)+be(2)*p(i,2)+be(3)*p(i,3) 
      do 540 j=1,3 
      p(i,j)=p(i,j)+ga*(ga*bep/(1.d0+ga)+p(i,4))*be(j) 
  540 continue 
      p(i,4)=ga*(p(i,4)+bep) 
  550 continue 
  560 continue 
 
c...Check that no infinite loop in matrix element weight. 
      ntry=ntry+1 
      if(ntry.gt.800) goto 590 
 
c...Matrix elements for omega and phi decays. 
      if(mmat.eq.1) then 
        wt=(p(n+1,5)*p(n+2,5)*p(n+3,5))**2-(p(n+1,5)*four(n+2,n+3))**2 
     &  -(p(n+2,5)*four(n+1,n+3))**2-(p(n+3,5)*four(n+1,n+2))**2 
     &  +2.d0*four(n+1,n+2)*four(n+1,n+3)*four(n+2,n+3) 
        if(max(wt*wtcor(9)/p(nsav,5)**6,0.001d0).lt.pjr(0)) goto 420 
 
c...Matrix elements for pi0 or eta Dalitz decay to gamma e+ e-. 
      elseif(mmat.eq.2) then 
        four12=four(n+1,n+2) 
        four13=four(n+1,n+3) 
        wt=(pmst-0.5d0*pmes)*(four12**2+four13**2)+ 
     &  pmes*(four12*four13+four12**2+four13**2) 
        if(wt.lt.pjr(0)*0.25d0*pmst*(p(nsav,5)**2-pmst)**2) goto 490 
 
c...Matrix element for S0 -> S1 + V1 -> S1 + S2 + S3 (S scalar, 
c...V vector), of form cos**2(theta02) in V1 rest frame, and for 
c...S0 -> gamma + V1 -> gamma + S2 + S3, of form sin**2(theta02). 
      elseif(mmat.eq.3.and.nm.eq.2) then 
        four10=four(nsav,im) 
        four12=four(nsav,n+1) 
        four02=four(im,n+1) 
        pms1=p(nsav,5)**2 
        pms0=p(im,5)**2 
        pms2=p(n+1,5)**2 
        if(kfas.ne.22) hnum=(four10*four12-pms1*four02)**2 
        if(kfas.eq.22) hnum=pms1*(2.d0*four10*four12*four02- 
     &  pms1*four02**2-pms0*four12**2-pms2*four10**2+pms1*pms0*pms2) 
        hnum=max(1d-6*pms1**2*pms0*pms2,hnum) 
        hden=(four10**2-pms1*pms0)*(four12**2-pms1*pms2) 
        if(hnum.lt.pjr(0)*hden) goto 490 
 
c...Matrix element for "onium" -> g + g + g or gamma + g + g. 
      elseif(mmat.eq.4) then 
        hx1=2.d0*four(nsav,n+1)/p(nsav,5)**2 
        hx2=2.d0*four(nsav,n+2)/p(nsav,5)**2 
        hx3=2.d0*four(nsav,n+3)/p(nsav,5)**2 
        wt=((1.d0-hx1)/(hx2*hx3))**2+((1.d0-hx2)/(hx1*hx3))**2+ 
     &  ((1.d0-hx3)/(hx1*hx2))**2 
        if(wt.lt.2.d0*pjr(0)) goto 420 
        if(k(nsav+1,2).eq.22
     $      .and.(1.d0-hx1)*p(nsav,5)**2.lt.4.d0*parj(32)**2) 
     &  goto 420 

      endif 
 
  590 continue

c...Scale back energy and reattach spectator. 
c 590 if(mrem.eq.1) then 
c       do 600 j=1,5 
c       pv(1,j)=pv(1,j)/(1.-pqt) 
c 600   continue 
c       nd=nd+1 
c       mrem=0 
c     endif 
 
c...Check invariant mass of w jets. may give one particle or start over.
 
 
c...Boost back for rapidly moving particle. 
  660 n=n+nd 
      if(mbst.eq.1) then 
        do 670 j=1,3 
        be(j)=p(nsav,j)/p(nsav,4) 
  670   continue 
        ga=p(nsav,4)/p(nsav,5) 
        do 690 i=nsav+1,n 
        bep=be(1)*p(i,1)+be(2)*p(i,2)+be(3)*p(i,3) 
        do 680 j=1,3 
        p(i,j)=p(i,j)+ga*(ga*bep/(1.d0+ga)+p(i,4))*be(j) 
  680   continue 
        p(i,4)=ga*(p(i,4)+bep) 
  690   continue 
      endif 
 
c...Fill in position of decay vertex. 
c     do 710 i=nsav+1,n 
c     do 700 j=1,4 
c     v(i,j)=vdcy(j) 
c 700 continue 
c     v(i,5)=0. 
c 710 continue 
 
 800  continue 
c...Mark decayed particle; special option for B-B~ mixing. 
      if(k(nsav,1).eq.5) k(nsav,1)=15 
      if(k(nsav,1).le.10) k(nsav,1)=11 
      if(mmix.eq.1.and.mstj(26).eq.2.and.k(nsav,1).eq.11) k(nsav,1)=12 
      k(nsav,4)=nsav+1 
      k(nsav,5)=n 
 
      return 
      end 
 

c***********************************************************************

      subroutine jamwidm(kc,kfsp,kfsn,kf1,kf2,emcm,ibranch,pwid,totwid,
     $                 itag)

c...Purpose: to calculate momentum dependent partical decay width
c...for resonances.
c======================================================================*
c     kc     : Compressed particle code                                *
c     kf1 kf2: ingoing particle KF code.                               *
c     emcm   :  Mass of the particle (GeV)                             *
c     pwid   :  partial decay width  (GeV) (output)                    *
c     totwid :  total width                                            *
c     itag   : branch which is identical to the kf1 and kf2.           *
c======================================================================*

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      real*8 jamdlwid
      parameter(maxbr=70)
      dimension ibranch(maxbr),pwid(maxbr)
      parameter(delta=0.09d0)
c...Functions: momentum in two-particle cm.
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a) 

      icon=0
      totwid=0.0d0
      itag=0

      minb=mdcy(kc,2)
      maxb=mdcy(kc,2)+mdcy(kc,3)-1

c....Special for Delta(1232).
      if(kchg(kc,5).eq.id_delt) then
         totwd=jamdlwid(mstc(63),emcm)
         i=0
         itag=1
         do idl=minb,maxb
           i=i+1
           ibranch(i)=0
           pwid(i)=0.0d0
           if(totwd.gt.1.d-5) then
             ibranch(i)=1
             pwid(i)=totwd*brat(idl)
             totwid=totwid+pwid(i)
           endif
         end do
         goto 1000
      endif

c...Resonance pole mass and width.
      emres=pmas(kc,1)
      widr=pmas(kc,2)
      ibra=0
      do 10 idl = minb,maxb
        ibra=ibra+1
        pwid(ibra)=0.0d0
        ibranch(ibra)=0
        if(mdme(idl,1).ne.1.and.kfsp*mdme(idl,1).ne.2.and. 
     &                         kfsn*mdme(idl,1).ne.3) goto 10
        if(mdme(idl,2).gt.100) goto 10
        widp=brat(idl)
        if(widp.le.1.d-5)  goto 10

        kfd3=kfdp(idl,3)
        kfd1=kfdp(idl,1)
        kfd2=kfdp(idl,2)
c.....Lepton, etc.
        if(abs(kfd1).le.100.or.abs(kfd2).le.100) goto 10

        kcd1=jamcomp(kfd1)
        kcd2=jamcomp(kfd2)
        if(kcd1.eq.0.or.kcd2.eq.0) goto 10
        em1=pmas(kcd1,1)-pmas(kcd1,3)
        em2=pmas(kcd2,1)-pmas(kcd2,3)

cbug fix 07/17/2003
c       if(emres.lt.em1+em2+parj(64)) goto 10
        if(emcm.lt.em1+em2+parj(64)) goto 10

        if( (kfd1.eq.kf2.and.kfd2.eq.kf1.and.kfd3.eq.0)
     $      .or.(kfd1.eq.kf1.and.kfd2.eq.kf2.and.kfd3.eq.0) ) then
            itag=ibra
        endif

c....Three body decay
        if(mstc(65).eq.0.or.kfdp(idl,3).ne.0) then
          pwid(ibra)=widp*widr

c...Mom. dep. width
cbug fix 07/17/2003
        else if(emres.gt.em1+em2+parj(64)) then
          prres=pawt(emres,em1,em2)
          ldec=mdme(idl,3)
          if(emcm.ge.em1+em2+parj(64)) then
            pr=pawt(emcm,em1,em2)
            if(mstc(63).eq.2) then 
              form=((prres**2+delta)/(pr**2+delta))**(ldec+1)
            else
              form=1.2d0/(1.d0+0.2d0*(pr/prres)**(2*ldec))
            endif
          else
            goto 10
          end if
          pwid(ibra)=widp*(pr/prres)**(2*ldec+1)*(emres/emcm)*form*widr
        else
          pwid(ibra)=widp*widr
        endif

        if(pwid(ibra).gt.1.d-4) ibranch(ibra)=1
        totwid=totwid+pwid(ibra)

 10   continue
1000  continue

      end
      
c*******************************************************************

      subroutine jambwmas(emmin,emmax,emr,wid,em,icon)

c...Purpose: to generate mass according to the B.W. distribution.
c==================================================================*
c  emmin : minimam mass           (input)
c  emmax : max. mass              (input)
c  emr   : resonance peak mass    (input)
c  wid   : resonance full width   (input)
c  em    : resonance mass         (output)
c==================================================================*
      implicit double precision(a-h, o-z)
      parameter( pi=3.141593d0 )
      icon=0

c...Check boundary.
      if(emmax.le.emmin) then
        em=emmax
        icon=999
        return
      endif

c...Breit Wigner distribution.
      const=2.d0*(emmin-emr)/wid
      const1=atan(const)
      const2=pi/2.d0-const1
      xmax=(atan(2.d0*(emmax-emr)/wid)-const1)/const2
      xmax=min(1.0d0,xmax)
      x=xmax*rn(0)
      t=tan(x*const2)
      em=emr+0.5d0*wid*(const+t)/(1.0d0-const*t)

      end

c***********************************************************************

      double precision function jamdtim(im,kf,kc,ks01,emd,ee)

c...Purpose: to calculate life time of the unstable particle
c======================================================================*
c variables:                                                           *
c im       -  switch for action (input)                                *
c             0: decay width (GeV)                                     *
c             1: predicted decay time in comp. frame (fm/c)            *
c kf       -  particle code                                            *
c kc       -  Compressed particle code                                 *
c ks01     -  status code, =3: jet systm                               *
c emd      -  particle mass                                            *
c ee       -  particle energy in the comp. frame                       *
c jamdtim  -  decay time in the comp. frame        (dble,input)        *
c======================================================================*

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      real*8 jamdlwid
      parameter(maxbr=70)
      dimension ibranch(maxbr),pwid(maxbr)

      jamdtim=1d+34
      ks1=mod(abs(ks01),10)
      if(kc.le.0.or.kc.gt.mstu(6)) then
       write(check(1),'(''kc kf='',i4,i9)')kc,kf
       call jamerrm(30,1,'(jamdtim:)invalid kc')
      endif

c....Jet system.
      if(ks1.eq.3.or.kf.eq.92) then

        if(mstc(76).eq.2) then ! not working now (2002.10)
          totwid=parc(57)*emd
          if(totwid.lt.1.d-7) then
            jamdtim=0.0d0
            if(im.eq.0) jamdtim=totwid
            return
          endif
        else
          jamdtim=0.0d0
          if(im.eq.0) jamdtim=totwid
          return
        endif

c...Unstable hadrons.
      else

c.....Calculate momentum dependent width.
          if(mstc(65).eq.1) then

c...........Special for Delta decay.
              if(kchg(kc,5).eq.id_delt) then
                totwid=jamdlwid(mstc(63),emd)
c...........effective s-wave N(1440)
              else if((abs(kf).eq.12212.or.abs(kf).eq.12112)
     $                                         .and.emd.le.1.2d0)then
                totwid=pmas(kc,2)
              else

                if(pmas(kc,2).le.1d-7.or.mdcy(kc,1).eq.0
     $              .or.mdcy(kc,2).eq.0.or.mdcy(kc,3).eq.0) then
                  jamdtim=1.d+30
                  if(im.eq.0) jamdtim=0.0d0
                  return
                endif

                kfa=abs(kf)
                kfs=isign(1,kf)
c...B-B~ mixing: flip sign of meson appropriately. 
                mmix=0 
                if((kfa.eq.511.or.kfa.eq.531).and.mstj(26).ge.1) then 
                  xbbmix=parj(76) 
                  if(kfa.eq.531) xbbmix=parj(77) 
                  vip5=-pmas(kc,4)*log(pjr(0))
                  if(sin(0.5d0*xbbmix*vip5/pmas(kc,4))**2.gt.pjr(0))
     $                mmix=1 
                  if(mmix.eq.1) kfs=-kfs 
                endif 

                if(kchg(kc,3).eq.0) then 
                  kfsp=1 
                  kfsn=0 
c                 if(rlu(0).gt.0.5) kfs=-kfs 
                elseif(kfs.gt.0) then 
                  kfsp=1 
                  kfsn=0 
                else 
                  kfsp=0 
                  kfsn=1 
                endif 
                call jamwidm(kc,kfsp,kfsn,0,0,
     $                 emd,ibranch,pwid,totwid,itag)

              endif

          else
              totwid=pmas(kc,2)
          endif

      endif

      if(im.eq.0) then
        jamdtim=totwid
        return
      endif

c...Convert width from GeV to (fm/c)^-1
      totwid   = totwid / paru(3)

c...Compulte rest frame decay time.
      if(totwid.gt.1.d-7) then
        tlife = - log( max( rn(0), 1.d-35 ) ) / totwid
      else
        jamdtim=1.d+22
        return
      endif

c...Gamma factor.
      if(emd.ge.1.d-4) then
        gg=ee/emd
      else
        gg=1.d0
      endif

c...Apply time dilation.
      jamdtim=tlife*gg

      return
      end

c***********************************************************************

      double precision function jamdlwid(iwidth,emd)

c....Purpose: to calculate momentum dependent total decay width
c...  of delta(1232).
c
c...iwidth = 0: no decay(frozen delta)
c...       1: Frankfurt
c...       2: Giessen
c...       3: Randrup
c...           Ref: Randrup, NP A314 (1979) 429.
c...                Rittenberg, REV.MOD.PHYS. 43 (1971) S1.
c...       4: Kitazoe
c...       5: Barz/Iwe  NP A453(1986)728

      implicit double precision(a-h, o-z)
      parameter(emnuc=0.9383d0,empion=0.138d0,ekinmi=0.0001d0)
      parameter(emdelt=1.232d0,widdlt=0.12d0)
      parameter(bet2=0.090d0, qqr2=0.051936d0, gamr=0.11d0)
      parameter(pscal1= 0.238d0, pscal2= 0.318d0, p0ref=0.227d0)
      character chekc*80
c...Functions: momentum in two-particle cm.
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a) 

      if(emd.lt.emnuc+empion+ekinmi) then
c        call jamlist(1)
         write(chekc,'(g14.6)')emd
         call jamerrm(30,0,'(jamdlwid:)invalid d(1232) mass'//chekc)
       else
        pp=pawt(emd,emnuc,empion)
       endif

      pp2=pp*pp
      if(pp.le.0.0d0) then
          jamdlwid=0.0d0
      else

        if(iwidth.eq.1) then
           jamdlwid=0.12d0*emdelt/emd*sqrt(pp2/qqr2)**3*1.2d0/(1 
     & +0.2d0*pp2/qqr2)
        else if(iwidth.eq.2) then
            form= (1.d0+qqr2/bet2)/(1.d0+pp2/bet2)
            jamdlwid=sqrt(pp2/qqr2)**3*emdelt/emd*gamr*form**2
        else if(iwidth.eq.3) then
            jamdlwid=widdlt*(pp**3/(1.d0+(pp/pscal1)**2+(pp/pscal2)**4))
     a     /(p0ref**3/(1.d0+(p0ref/pscal1)**2+(p0ref/pscal2)**4))
        else if(iwidth.eq.4) then
           jamdlwid=0.47d0/(1.0d0+0.6d0*pp2/empion**2)*pp2/empion**2*pp
        else if(iwidth.eq.5) then
         jamdlwid=29.d0*pp**3/(1.d0+40.d0*pp2)
        else
           jamdlwid=widdlt
        endif

      endif

      end

 
c***********************************************************************

      subroutine jamfdec

c...Decay unstable particles at the end of simulation.
      include 'jam1.inc'
      include 'jam2.inc'
      dimension indd(100)

      mste(40)=1

c....Option for weak decay.
      if(mstc(42).eq.0) call jamsetd(1)

      itry=0
2200  continue
      itry=itry+1
      if(itry.ge.1000) then
        call jamerrm(3,0,
     $   '(jamfdec:) Infinit loop??? 1000 decay after simulation')
        call jamlist(1)
        return
      endif

      do 200 ip=1,nv

        k1=k(1,ip)
        kf=k(2,ip)

c....Dead particle.
        if(k1.gt.10) goto 200
c.....This particle is still within a formation time.
        if(k1.lt.0)  k(1,ip)=mod(abs(k1),10)

        kc=jamcomp(kf)
        if(kf.ne.92.and.k1.ne.4.and.mdcy(kc,1).eq.0) goto 200
        if(mstc(42).ne.0.and.k(1,ip).le.1) goto 200

c...Forbid unstable particle decay.
        if(mstc(41).eq.0.and.(kf.ne.92.and.k1.ne.4)) goto 200

c...Avoid numerical error for long lived particles.
c       if(v(5,ip).lt.1000.0d0) then
          dect=v(5,ip)-r(4,ip)
c       else
c         dect=1000.0d0
c       endif

cxxxxxxxxxxxxxxxxxxx
        if(mstc(42).ne.0.and.v(5,ip).gt.1e+7) then
          ih=mstc(38)
          write(ih,*)'v5? in final decay v5=',v(5,ip)
          write(ih,*)'kf=',k(2,ip),' k3=',k(3,ip),' k4=',k(4,ip)
          write(ih,*)'k=',(k(j,ip),j=1,7)
          write(ih,*)'r=',(r(j,ip),j=1,5)
          write(ih,*)'v=',(v(j,ip),j=1,5)
          write(ih,*)'p=',(p(j,ip),j=1,5)
c         stop
        endif
cxxxxxxxxxxxxxxxxxxx

        do j=1,3
         r(j,ip)=r(j,ip)+dect*p(j,ip)/p(4,ip)
        end do
        r(4,ip)=v(5,ip)

c.....qmd:Save initital energy (resonance decay only)
cq      if(mstc(57).ge.1.and.k(1,ip).eq.2) then
cq        call epotall(epot,epotpa)
cq        pare(11) = p(4,ip) + epot
cq      endif

c...Set collision type.
        mste(2)=-2
        call jamsave(1,1,ip)
        call jamdec(ip,indd,nadd,icon)

c...Update number of decay.
        mstd(53)=mstd(53)+1

cTABARA
c       call ttchk(indd,nadd)
c...Print information after decay.
        if(mstc(8).ge.2) call jamprcl(indd,nadd)

  200 continue

c...Find next decay if possible.
      do 300 ip=1,nv
          k1=k(1,ip)
          if(k1.gt.10) goto 300
          kf=k(2,ip)
c...Forbid unstable particle decay.
          if(mstc(41).eq.0.and.(kf.ne.92.and.k1.ne.4)) goto 300
          kc=jamcomp(kf)
          if(kc.le.0.or.kc.gt.mstu(6)) then
            write(6,*)'(jamfdec:)kf??',ip,k1,kf,kc
            goto 300
          endif
c...[AOtest
c       if(kf.eq.413) then
c         write(99,*) '413(D*):k1,kc=',k1,kc,(mdcy(kc,ik),ik=1,3)
c       endif
        if((k1.ge.2.and.k1.le.10).or.mdcy(kc,1).eq.1) then
        if(mstc(42).ne.0.and.pmas(kc,2).eq.0) then
c         write(99,'(1x,i7,1x,i3,1x,i3,3(1x,i4),2(1x,1pe10.3))')
c    &      kf,kc,k1,(mdcy(kc,ik),ik=1,3),(pmas(kc,ik),ik=1,2)
        else
          goto 2200
        endif
        endif
c...]AOtest
c         if(k1.ge.2.and.k1.le.10) goto 2200
c         if(mdcy(kc,1).eq.1) goto 2200
 300  continue

c....Reset decay swicth.
      if(mstc(42).eq.0) call jamsetd(0)

      mste(40)=0

      end

c***********************************************************************

      subroutine jamsetd(i)

c...Decay unstable particles at the end of simulation.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

c       mdcy(jamcomp(111),1)=i    ! pi0
        mdcy(jamcomp(-311),1)=i   ! ak0
        mdcy(jamcomp(311),1)=i    ! k0
        mdcy(jamcomp(310),1)=i    ! k0_S
c       mdcy(jamcomp(130),1)=i    ! k0_L
        mdcy(jamcomp(411),1)=i    ! D+
        mdcy(jamcomp(421),1)=i    ! D0
        mdcy(jamcomp(221),1)=i    ! eta
c       mdcy(jamcomp(331),1)=i    ! eta'
        mdcy(jamcomp(441),1)=i    ! eta_c
        mdcy(jamcomp(310),1)=i
        mdcy(jamcomp(431),1)=i
        mdcy(jamcomp(511),1)=i
        mdcy(jamcomp(521),1)=i
        mdcy(jamcomp(531),1)=i
        mdcy(jamcomp(3122),1)=i
        mdcy(jamcomp(3112),1)=i
        mdcy(jamcomp(3212),1)=i
        mdcy(jamcomp(3222),1)=i
        mdcy(jamcomp(3312),1)=i
        mdcy(jamcomp(3322),1)=i
        mdcy(jamcomp(3334),1)=i
c       call pjgive('mdcy(c111,1)=1')
c       call pjgive('mdcy(c3122,1)=1;mdcy(c-3122,1)=1')

      end