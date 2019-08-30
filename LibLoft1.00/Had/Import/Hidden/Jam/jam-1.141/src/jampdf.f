C*********************************************************************
C*********************************************************************
C                                                                    *
C  List of subprograms in order of appearance, with main purpose     *
C  (S = subroutine, F = function, B = block data)                    *
C                                                                    *
C                                                                    *
C  S   PYPDFU   to evaluate parton distributions                     *
C  S   PYPDFL   to evaluate parton distributions at low x and Q^2    *
C  S   PYPDEL   to evaluate electron parton distributions            *
C  S   PYPDGA   to evaluate photon parton distributions (generic)    *
C  S   PYGGAM   to evaluate photon parton distributions (SaS sets)   *
C  S   PYGVMD   to evaluate VMD part of photon parton distributions  *
C  S   PYGANO   to evaluate anomalous part of photon pdf's           *
C  S   PYGBEH   to evaluate Bethe-Heitler part of photon pdf's       *
C  S   PYGDIR   to evaluate direct contribution to photon pdf's      *
c                                                                    *
C  S   PYPDPI   to evaluate pion parton distributions                *
C  S   PYPDPR   to evaluate proton parton distributions              *
C  F   PYCTEQ   to evaluate the CTEQ 3 proton parton distributions   *
C  S   PYGRVL   to evaluate the GRV 94L pronton parton distributions *
C  S   PYGRVM   to evaluate the GRV 94M pronton parton distributions *
C  S   PYGRVD   to evaluate the GRV 94D pronton parton distributions *
C  F   PYGRVV   auxiliary to the PYGRV* routines                     *
C  F   PYGRVW   auxiliary to the PYGRV* routines                     *
C  F   PYGRVS   auxiliary to the PYGRV* routines                     *
c  f   pyctq2   to give the CTEQ 2 parton distribution function sets *
c  s   ehlq     to give proton structure functions from EHLQ         *
c  s   dukeowen to give Proton structure functions from Duke, Owens  *
C  S   PDFSET   dummy routine to be removed when using PDFLIB        *
C  S   STRUCTM  dummy routine to be removed when using PDFLIB        *
C                                                                    *
C*********************************************************************
 
C...PYPDFU
C...Gives electron, photon, pi+, neutron, proton and hyperon
C...parton distributions according to a few different parametrizations.
C...Note that what is coded is x times the probability distribution,
C...i.e. xq(x,Q2) etc.
 
      subroutine pjpdfu(kf,x,q2,xpq)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      include 'jam2.inc'
c     common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
c     common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint8/xpvmd(-6:6),xpanl(-6:6),xpanh(-6:6),xpbeh(-6:6),
     &xpdir(-6:6)
      save /pjpars/,/pjint1/,/pjint8/
C...Local arrays.
      dimension xpq(-25:25),xpel(-25:25),xpga(-6:6),vxpga(-6:6),
     &xppi(-6:6),xppr(-6:6)
      character check1*120
 
C...Interface to PDFLIB.
      common/w50513/xmin,xmax,q2min,q2max
      save /w50513/
      double precision xx,qq,upv,dnv,usea,dsea,str,chm,bot,top,glu,
     &value(20),xmin,xmax,q2min,q2max
      character*20 parm(20)
      data value/20*0d0/,parm/20*' '/
 
C...Data related to Schuler-Sjostrand photon distributions.
      data alamga/0.2d0/, pmcga/1.3d0/, pmbga/4.6d0/
 
C...Reset parton distributions.
      mint(92)=0
      do 100 kfl=-25,25
        xpq(kfl)=0d0
  100 continue
 
C...Check x and particle species.
      if(x.le.0d0.or.x.ge.1d0) then
        write(mstu(11),5000) x
        return
      endif
      kfa=iabs(kf)

c     if(kfa.ne.11.and.kfa.ne.22.and.kfa.ne.211.and.kfa.ne.2112.and.
c    &kfa.ne.2212.and.kfa.ne.3122.and.kfa.ne.3112.and.kfa.ne.3212
c    &.and.kfa.ne.3222.and.kfa.ne.3312.and.kfa.ne.3322.and.
c    &kfa.ne.3334.and.kfa.ne.111) then
c       write(mstu(11),5100) kf
c       return
c     endif

      kc=jamcomp(kf)
      id=kchg(kc,5)
      ibar=kchg(kc,6)
      iz=kchg(kc,1)
      imeso=0
      kfma=0
      ipion=999
c...Mesons.
      if(ibar.eq.0.and.kc.gt.100) then
        imeso=1
        ipion=0
        kfm1=mod(kfa/100,10)
        kfm2=mod(kfa/10,10)
        kfma=10*kfm1+kfm2
        if(kfm1.eq.2.and.kfm2.eq.1) ipion=1
        if(id.eq.id_str.and.kf.gt.0) then
           kfm1=-kfm1
        else if((kfma.eq.41.or.kfma.eq.42).and.kf.lt.0) then
           kfm1=-kfm1
        else if(kfma.eq.43.and.kf.lt.0) then
           kfm1=-kfm1
        else if((kfma.eq.51.or.kfma.eq.52).and.kf.gt.0) then
           kfm1=-kfm1
        else if(kfma.eq.53.and.kf.gt.0) then
           kfm1=-kfm1
        else
          kfm2=-kfm2
        endif
      endif
      if(kc.eq.0) then
        write(mstu(11),5100) kf
        return
      endif
 
C...Electron parton distribution call.
      if(kfa.eq.11) then
        call pjpdel(x,q2,xpel)
        do 110 kfl=-25,25
          xpq(kfl)=xpel(kfl)
  110   continue
 
C...Photon parton distribution call (VDM+anomalous).
      elseif(kfa.eq.22.and.mint(109).le.1) then
        if(mstp(56).eq.1.and.mstp(55).eq.1) then
          call pjpdga(x,q2,xpga)
          do 120 kfl=-6,6
            xpq(kfl)=xpga(kfl)
  120     continue
        elseif(mstp(56).eq.1.and.mstp(55).ge.5.and.mstp(55).le.8) then
          q2mx=q2
          p2mx=0.36d0
          if(mstp(55).ge.7) p2mx=4.0d0
          if(mstp(57).eq.0) q2mx=p2mx
          call pjggam(mstp(55)-4,x,q2mx,0d0,mstp(60),f2gam,xpga)
          do 130 kfl=-6,6
            xpq(kfl)=xpga(kfl)
  130     continue
          vint(231)=p2mx
        elseif(mstp(56).eq.1.and.mstp(55).ge.9.and.mstp(55).le.12) then
          q2mx=q2
          p2mx=0.36d0
          if(mstp(55).ge.11) p2mx=4.0d0
          if(mstp(57).eq.0) q2mx=p2mx
          call pjggam(mstp(55)-8,x,q2mx,0d0,mstp(60),f2gam,xpga)
          do 140 kfl=-6,6
            xpq(kfl)=xpvmd(kfl)+xpanl(kfl)+xpbeh(kfl)+xpdir(kfl)
  140     continue
          vint(231)=p2mx
        elseif(mstp(56).eq.2) then
C...Call PDFLIB parton distributions.
          parm(1)='NPTYPE'
          value(1)=3
          parm(2)='NGROUP'
          value(2)=mstp(55)/1000
          parm(3)='NSET'
          value(3)=mod(mstp(55),1000)
          if(mint(93).ne.3000000+mstp(55)) then
            call pjpdfset(parm,value)
            mint(93)=3000000+mstp(55)
          endif
          xx=x
          qq=sqrt(max(0d0,q2min,q2))
          if(mstp(57).eq.0) qq=sqrt(q2min)
          call pjstructm(xx,qq,upv,dnv,usea,dsea,str,chm,bot,top,glu)
          vint(231)=q2min
          xpq(0)=glu
          xpq(1)=dnv
          xpq(-1)=dnv
          xpq(2)=upv
          xpq(-2)=upv
          xpq(3)=str
          xpq(-3)=str
          xpq(4)=chm
          xpq(-4)=chm
          xpq(5)=bot
          xpq(-5)=bot
          xpq(6)=top
          xpq(-6)=top
        else
          write(mstu(11),5200) kf,mstp(56),mstp(55)
        endif
 
C...Pion/gammaVDM parton distribution call.
      elseif(kfa.eq.211.or.kfa.eq.111.or.(kfa.eq.22.and.
     &  mint(109).eq.2)) then
        if(kfa.eq.22.and.mstp(56).eq.1.and.mstp(55).ge.5.and.
     &  mstp(55).le.12) then
          iset=1+mod(mstp(55)-1,4)
          q2mx=q2
          p2mx=0.36d0
          if(iset.ge.3) p2mx=4.0d0
          if(mstp(57).eq.0) q2mx=p2mx
          call pjgvmd(iset,2,x,q2mx,p2mx,alamga,xpga,vxpga)
          do 150 kfl=-6,6
            xpq(kfl)=xpga(kfl)
  150     continue
          vint(231)=p2mx
        elseif(mstp(54).eq.1.and.mstp(53).ge.1.and.mstp(53).le.3) then
          call pjpdpi(x,q2,xppi)
          do 160 kfl=-6,6
            xpq(kfl)=xppi(kfl)
  160     continue
        elseif(mstp(54).eq.2) then
C...Call PDFLIB parton distributions.
          parm(1)='NPTYPE'
          value(1)=2
          parm(2)='NGROUP'
          value(2)=mstp(53)/1000
          parm(3)='NSET'
          value(3)=mod(mstp(53),1000)
          if(mint(93).ne.2000000+mstp(53)) then
            call pjpdfset(parm,value)
            mint(93)=2000000+mstp(53)
          endif
          xx=x
          qq=sqrt(max(0d0,q2min,q2))
          if(mstp(57).eq.0) qq=sqrt(q2min)
          call pjstructm(xx,qq,upv,dnv,usea,dsea,str,chm,bot,top,glu)
          vint(231)=q2min
          xpq(0)=glu
          xpq(1)=dsea
          xpq(-1)=upv+dsea
          xpq(2)=upv+usea
          xpq(-2)=usea
          xpq(3)=str
          xpq(-3)=str
          xpq(4)=chm
          xpq(-4)=chm
          xpq(5)=bot
          xpq(-5)=bot
          xpq(6)=top
          xpq(-6)=top
        else
          write(mstu(11),5200) kf,mstp(54),mstp(53)
        endif
 
C...Anomalous photon parton distribution call.
      elseif(kfa.eq.22.and.mint(109).eq.3) then
        q2mx=q2
        p2mx=parp(15)**2
        if(mstp(56).eq.1.and.mstp(55).le.8) then
          if(mstp(55).eq.5.or.mstp(55).eq.6) p2mx=0.36d0
          if(mstp(55).eq.7.or.mstp(55).eq.8) p2mx=4.0d0
          if(mstp(57).eq.0) q2mx=p2mx
          call pjgano(0,x,q2mx,p2mx,alamga,xpga,vxpga)
          do 170 kfl=-6,6
            xpq(kfl)=xpga(kfl)
  170     continue
          vint(231)=p2mx
        elseif(mstp(56).eq.1) then
          if(mstp(55).eq.9.or.mstp(55).eq.10) p2mx=0.36d0
          if(mstp(55).eq.11.or.mstp(55).eq.12) p2mx=4.0d0
          if(mstp(57).eq.0) q2mx=p2mx
          call pjggam(mstp(55)-8,x,q2mx,0d0,mstp(60),f2gm,xpga)
          do 180 kfl=-6,6
            xpq(kfl)=max(0d0,xpanl(kfl)+xpbeh(kfl)+xpdir(kfl))
  180     continue
          vint(231)=p2mx
        elseif(mstp(56).eq.2) then
          if(mstp(57).eq.0) q2mx=p2mx
          call pjgano(0,x,q2mx,p2mx,alamga,xpga,vxpga)
          do 190 kfl=-6,6
            xpq(kfl)=xpga(kfl)
  190     continue
          vint(231)=p2mx
        elseif(mstp(55).ge.1.and.mstp(55).le.5) then
          if(mstp(57).eq.0) q2mx=p2mx
          call pjgvmd(0,mstp(55),x,q2mx,p2mx,parp(1),xpga,vxpga)
          do 200 kfl=-6,6
            xpq(kfl)=xpga(kfl)
  200     continue
          vint(231)=p2mx
        else
  210     rkf=11d0*pjr(0)
          kfr=1
          if(rkf.gt.1d0) kfr=2
          if(rkf.gt.5d0) kfr=3
          if(rkf.gt.6d0) kfr=4
          if(rkf.gt.10d0) kfr=5
          if(kfr.eq.4.and.q2.lt.pmcga**2) goto 210
          if(kfr.eq.5.and.q2.lt.pmbga**2) goto 210
          if(mstp(57).eq.0) q2mx=p2mx
          call pjgvmd(0,kfr,x,q2mx,p2mx,parp(1),xpga,vxpga)
          do 220 kfl=-6,6
            xpq(kfl)=xpga(kfl)
  220     continue
          vint(231)=p2mx
        endif
 
C...jam:Proton parton distribution call.
      else if(iabs(ibar).eq.3) then

        if(mstp(52).eq.1.and.mstp(51).ge.1.and.mstp(51).le.17) then
          call pjpdpr(x,q2,xppr)
          do 230 kfl=-6,6
            xpq(kfl)=xppr(kfl)
  230     continue
        elseif(mstp(52).eq.2) then
C...Call PDFLIB parton distributions.
          parm(1)='NPTYPE'
          value(1)=1
          parm(2)='NGROUP'
          value(2)=mstp(51)/1000
          parm(3)='NSET'
          value(3)=mod(mstp(51),1000)
          if(mint(93).ne.1000000+mstp(51)) then
            call pjpdfset(parm,value)
            mint(93)=1000000+mstp(51)
          endif
          xx=x
          qq=sqrt(max(0d0,q2min,q2))
          if(mstp(57).eq.0) qq=sqrt(q2min)
          call pjstructm(xx,qq,upv,dnv,usea,dsea,str,chm,bot,top,glu)
          vint(231)=q2min
          xpq(0)=glu
          xpq(1)=dnv+dsea
          xpq(-1)=dsea
          xpq(2)=upv+usea
          xpq(-2)=usea
          xpq(3)=str
          xpq(-3)=str
          xpq(4)=chm
          xpq(-4)=chm
          xpq(5)=bot
          xpq(-5)=bot
          xpq(6)=top
          xpq(-6)=top
        else
          write(mstu(11),5200) kf,mstp(52),mstp(51)
        endif
      else
        write(mstu(11),5200) kf,mstp(52),mstp(51)
        call pjerrm(30,'(pjpdfu:)unknown kf')
      endif

cjam++
      if(iabs(ibar).eq.3) then

C...Isospin conjugation for neutron.
        if(id.eq.id_nucl.or.id.eq.id_nucls
     $       .or.id.eq.id_delt.or.id.eq.id_delts) then
          if(iz.eq.0) then
            xps=xpq(1)
            xpq(1)=xpq(2)
            xpq(2)=xps
            xps=xpq(-1)
            xpq(-1)=xpq(-2)
            xpq(-2)=xps
            goto 1000
          endif
        endif
 
C...Simple recipes for hyperon (average valence parton distribution).
c     elseif(kfa.eq.3122.or.kfa.eq.3112.or.kfa.eq.3212.or.kfa.eq.3222
c    &  .or.kfa.eq.3312.or.kfa.eq.3322.or.kfa.eq.3334) then

        xpval=(xpq(1)+xpq(2)-xpq(-1)-xpq(-2))/3d0
        xpsea=0.5d0*(xpq(-1)+xpq(-2))
        xpq(1)=xpsea
        xpq(2)=xpsea
        xpq(-1)=xpsea
        xpq(-2)=xpsea
        xpq(mod(kfa/1000,10))=xpq(mod(kfa/1000,10))+xpval   ! JAM
        xpq(mod(kfa/100,10))=xpq(mod(kfa/100,10))+xpval
        xpq(mod(kfa/10,10))=xpq(mod(kfa/10,10))+xpval

c...Mesons.
      else if(imeso.eq.1) then
 
        if(ipion.eq.0)then
C...Isospin average for pi0
          if(kfma.eq.11) then
            xps=0.5d0*(xpq(1)+xpq(-2))
            xpv=0.5d0*(xpq(2)+xpq(-1))-xps
            xpq(2)=xps
            xpq(-1)=xps
            xpq(1)=xpq(1)+0.5d0*xpv
            xpq(-1)=xpq(-1)+0.5d0*xpv
            xpq(2)=xpq(2)+0.5d0*xpv
            xpq(-2)=xpq(-2)+0.5d0*xpv
          else
            xpval=0.5d0*(xpq(2)+xpq(-1)-xpq(-2)-xpq(1))
            spsea=0.5d0*(xpq(1)+xpq(-2))
            xpq(1)=xpsea
            xpq(2)=xpsea
            xpq(-1)=xpsea
            xpq(-2)=xpsea
            xpq(kfm1)=xpq(kfm1)+xpval
            xpq(kfm2)=xpq(kfm2)+xpval
          endif
        endif

C...Isospin average for pi0/gammaVDM.
      else if(kfa.eq.22.and.mint(109).eq.2) then !(

        if(mstp(55).ge.5.and.mstp(55).le.12) then
          xpv=xpq(2)-xpq(1)
          xpq(2)=xpq(1)
          xpq(-2)=xpq(-1)
        else
          xps=0.5d0*(xpq(1)+xpq(-2))
          xpv=0.5d0*(xpq(2)+xpq(-1))-xps
          xpq(2)=xps
          xpq(-1)=xps
        endif

        if(mint(105).le.223) then
          xpq(1)=xpq(1)+0.2d0*xpv
          xpq(-1)=xpq(-1)+0.2d0*xpv
          xpq(2)=xpq(2)+0.8d0*xpv
          xpq(-2)=xpq(-2)+0.8d0*xpv
        elseif(mint(105).eq.333) then
          xpq(3)=xpq(3)+xpv
          xpq(-3)=xpq(-3)+xpv
        elseif(mint(105).eq.443) then
          xpq(4)=xpq(4)+xpv
          xpq(-4)=xpq(-4)+xpv
          if(mstp(55).ge.9) then
            do 240 kfl=-6,6
              xpq(kfl)=0d0
  240       continue
          endif
        else
          xpq(1)=xpq(1)+0.5d0*xpv
          xpq(-1)=xpq(-1)+0.5d0*xpv
          xpq(2)=xpq(2)+0.5d0*xpv
          xpq(-2)=xpq(-2)+0.5d0*xpv
        endif
 
C...Rescale for gammaVDM by effective gamma -> rho coupling.
        if(mint(109).eq.2) then
          do 250 kfl=-6,6
            xpq(kfl)=vint(281)*xpq(kfl)
  250     continue
          vint(232)=vint(281)*xpv
        endif
 
      else !)
        write(check1,'(i9)')kf
        call pjerrm(30,'(pjpdfu:)unknown kf'//check1)
      endif !)
cjam--
 
C...Charge conjugation for antiparticle.
 1000 if(kf.lt.0) then
        do 260 kfl=1,25
          if(kfl.eq.21.or.kfl.eq.22.or.kfl.eq.23.or.kfl.eq.25) goto 260
          xps=xpq(kfl)
          xpq(kfl)=xpq(-kfl)
          xpq(-kfl)=xps
  260   continue
      endif
 
C...Allow gluon also in position 21.
      xpq(21)=xpq(0)
 
C...Check positivity and reset above maximum allowed flavour.
      do 270 kfl=-25,25
        xpq(kfl)=max(0d0,xpq(kfl))
        if(iabs(kfl).gt.mstp(58).and.iabs(kfl).le.8) xpq(kfl)=0d0
  270 continue
 
C...Formats for error printouts.
 5000 format(' Error: x value outside physical range; x =',1p,d12.3)
 5100 format(' Error: illegal particle code for parton distribution;',
     &' KF =',i5)
 5200 format(' Error: unknown parton distribution; KF, library, set =',
     &3i5)
 
      return
      end
 
C*********************************************************************
 
C...PYPDFL
C...Gives proton parton distribution at small x and/or Q^2 according to
C...correct limiting behaviour.
 
      subroutine pjpdfl(kf,x,q2,xpq)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      save /jydat1/,/jydat2/,/pjpars/,/pjint1/
C...Local arrays.
      dimension xpq(-25:25),xpa(-25:25),xpb(-25:25),wtsb(-3:3)
      data rmr/0.92d0/,rmp/0.38d0/,wtsb/0.5d0,1d0,1d0,5d0,1d0,1d0,0.5d0/
 
C...Send everything but protons/neutrons/VMD pions directly to PYPDFU.
      mint(92)=0
      kfa=iabs(kf)
      iacc=0
      if((kfa.eq.2212.or.kfa.eq.2112).and.mstp(57).ge.2) iacc=1
      if(kfa.eq.211.and.mstp(57).ge.3) iacc=1
      if(kfa.eq.22.and.mint(109).eq.2.and.mstp(57).ge.3) iacc=1
      if(iacc.eq.0) then
        call pjpdfu(kf,x,q2,xpq)
        return
      endif
 
C...Reset. Check x.
      do 100 kfl=-25,25
        xpq(kfl)=0d0
  100 continue
      if(x.le.0d0.or.x.ge.1d0) then
        write(mstu(11),5000) x
        return
      endif
 
C...Define valence content.
      kfc=kf
      nv1=2
      nv2=1
      if(kf.eq.2212) then
        kfv1=2
        kfv2=1
      elseif(kf.eq.-2212) then
        kfv1=-2
        kfv2=-1
      elseif(kf.eq.2112) then
        kfv1=1
        kfv2=2
      elseif(kf.eq.-2112) then
        kfv1=-1
        kfv2=-2
      elseif(kf.eq.211) then
        nv1=1
        kfv1=2
        kfv2=-1
      elseif(kf.eq.-211) then
        nv1=1
        kfv1=-2
        kfv2=1
      elseif(mint(105).le.223) then
        kfv1=1
        wtv1=0.2d0
        kfv2=2
        wtv2=0.8d0
      elseif(mint(105).eq.333) then
        kfv1=3
        wtv1=1.0d0
        kfv2=1
        wtv2=0.0d0
      elseif(mint(105).eq.443) then
        kfv1=4
        wtv1=1.0d0
        kfv2=1
        wtv2=0.0d0
      endif
 
C...Do naive evaluation and find min Q^2, boundary Q^2 and x_0.
      call pjpdfu(kfc,x,q2,xpa)
      q2mn=max(3d0,vint(231))
      q2b=2d0+0.052d0**2*exp(3.56d0*sqrt(max(0d0,-log(3d0*x))))
      xmn=exp(-(log((q2mn-2d0)/0.052d0**2)/3.56d0)**2)/3d0
 
C...Large Q2 and large x: naive call is enough.
      if(q2.gt.q2mn.and.q2.gt.q2b) then
        do 110 kfl=-25,25
          xpq(kfl)=xpa(kfl)
  110   continue
        mint(92)=1
 
C...Small Q2 and large x: dampen boundary value.
      elseif(x.gt.xmn) then
 
C...Evaluate at boundary and define dampening factors.
        call pjpdfu(kfc,x,q2mn,xpa)
        fv=(q2*(q2mn+rmr)/(q2mn*(q2+rmr)))**(0.55d0*(1d0-x)/(1d0-xmn))
        fs=(q2*(q2mn+rmp)/(q2mn*(q2+rmp)))**1.08d0
 
C...Separate valence and sea parts of parton distribution.
        if(kfa.ne.22) then
          xfv1=xpa(kfv1)-xpa(-kfv1)
          xpa(kfv1)=xpa(-kfv1)
          xfv2=xpa(kfv2)-xpa(-kfv2)
          xpa(kfv2)=xpa(-kfv2)
        else
          xpa(kfv1)=xpa(kfv1)-wtv1*vint(232)
          xpa(-kfv1)=xpa(-kfv1)-wtv1*vint(232)
          xpa(kfv2)=xpa(kfv2)-wtv2*vint(232)
          xpa(-kfv2)=xpa(-kfv2)-wtv2*vint(232)
        endif
 
C...Dampen valence and sea separately. Put back together.
        do 120 kfl=-25,25
          xpq(kfl)=fs*xpa(kfl)
  120   continue
        if(kfa.ne.22) then
          xpq(kfv1)=xpq(kfv1)+fv*xfv1
          xpq(kfv2)=xpq(kfv2)+fv*xfv2
        else
          xpq(kfv1)=xpq(kfv1)+fv*wtv1*vint(232)
          xpq(-kfv1)=xpq(-kfv1)+fv*wtv1*vint(232)
          xpq(kfv2)=xpq(kfv2)+fv*wtv2*vint(232)
          xpq(-kfv2)=xpq(-kfv2)+fv*wtv2*vint(232)
        endif
        mint(92)=2
 
C...Large Q2 and small x: interpolate behaviour.
      elseif(q2.gt.q2mn) then
 
C...Evaluate at extremes and define coefficients for interpolation.
        call pjpdfu(kfc,xmn,q2mn,xpa)
        vi232a=vint(232)
        call pjpdfu(kfc,x,q2b,xpb)
        vi232b=vint(232)
        fla=log(q2b/q2)/log(q2b/q2mn)
        fva=(x/xmn)**0.45d0*fla
        fsa=(x/xmn)**(-0.08d0)*fla
        fb=1d0-fla
 
C...Separate valence and sea parts of parton distribution.
        if(kfa.ne.22) then
          xfva1=xpa(kfv1)-xpa(-kfv1)
          xpa(kfv1)=xpa(-kfv1)
          xfva2=xpa(kfv2)-xpa(-kfv2)
          xpa(kfv2)=xpa(-kfv2)
          xfvb1=xpb(kfv1)-xpb(-kfv1)
          xpb(kfv1)=xpb(-kfv1)
          xfvb2=xpb(kfv2)-xpb(-kfv2)
          xpb(kfv2)=xpb(-kfv2)
        else
          xpa(kfv1)=xpa(kfv1)-wtv1*vi232a
          xpa(-kfv1)=xpa(-kfv1)-wtv1*vi232a
          xpa(kfv2)=xpa(kfv2)-wtv2*vi232a
          xpa(-kfv2)=xpa(-kfv2)-wtv2*vi232a
          xpb(kfv1)=xpb(kfv1)-wtv1*vi232b
          xpb(-kfv1)=xpb(-kfv1)-wtv1*vi232b
          xpb(kfv2)=xpb(kfv2)-wtv2*vi232b
          xpb(-kfv2)=xpb(-kfv2)-wtv2*vi232b
        endif
 
C...Interpolate for valence and sea. Put back together.
        do 130 kfl=-25,25
          xpq(kfl)=fsa*xpa(kfl)+fb*xpb(kfl)
  130   continue
        if(kfa.ne.22) then
          xpq(kfv1)=xpq(kfv1)+(fva*xfva1+fb*xfvb1)
          xpq(kfv2)=xpq(kfv2)+(fva*xfva2+fb*xfvb2)
        else
          xpq(kfv1)=xpq(kfv1)+wtv1*(fva*vi232a+fb*vi232b)
          xpq(-kfv1)=xpq(-kfv1)+wtv1*(fva*vi232a+fb*vi232b)
          xpq(kfv2)=xpq(kfv2)+wtv2*(fva*vi232a+fb*vi232b)
          xpq(-kfv2)=xpq(-kfv2)+wtv2*(fva*vi232a+fb*vi232b)
        endif
        mint(92)=3
 
C...Small Q2 and small x: dampen boundary value and add term.
      else
 
C...Evaluate at boundary and define dampening factors.
        call pjpdfu(kfc,xmn,q2mn,xpa)
        fb=(xmn-x)*(q2mn-q2)/(xmn*q2mn)
        fa=1d0-fb
        fvc=(x/xmn)**0.45d0*(q2/(q2+rmr))**0.55d0
        fva=fvc*fa*((q2mn+rmr)/q2mn)**0.55d0
        fvb=fvc*fb*1.10d0*xmn**0.45d0*0.11d0
        fsc=(x/xmn)**(-0.08d0)*(q2/(q2+rmp))**1.08d0
        fsa=fsc*fa*((q2mn+rmp)/q2mn)**1.08d0
        fsb=fsc*fb*0.21d0*xmn**(-0.08d0)*0.21d0
 
C...Separate valence and sea parts of parton distribution.
        if(kfa.ne.22) then
          xfv1=xpa(kfv1)-xpa(-kfv1)
          xpa(kfv1)=xpa(-kfv1)
          xfv2=xpa(kfv2)-xpa(-kfv2)
          xpa(kfv2)=xpa(-kfv2)
        else
          xpa(kfv1)=xpa(kfv1)-wtv1*vint(232)
          xpa(-kfv1)=xpa(-kfv1)-wtv1*vint(232)
          xpa(kfv2)=xpa(kfv2)-wtv2*vint(232)
          xpa(-kfv2)=xpa(-kfv2)-wtv2*vint(232)
        endif
 
C...Dampen valence and sea separately. Add constant terms.
C...Put back together.
        do 140 kfl=-25,25
          xpq(kfl)=fsa*xpa(kfl)
  140   continue
        if(kfa.ne.22) then
          do 150 kfl=-3,3
            xpq(kfl)=xpq(kfl)+fsb*wtsb(kfl)
  150     continue
          xpq(kfv1)=xpq(kfv1)+(fva*xfv1+fvb*nv1)
          xpq(kfv2)=xpq(kfv2)+(fva*xfv2+fvb*nv2)
        else
          do 160 kfl=-3,3
            xpq(kfl)=xpq(kfl)+vint(281)*fsb*wtsb(kfl)
  160     continue
          xpq(kfv1)=xpq(kfv1)+wtv1*(fva*vint(232)+fvb*vint(281))
          xpq(-kfv1)=xpq(-kfv1)+wtv1*(fva*vint(232)+fvb*vint(281))
          xpq(kfv2)=xpq(kfv2)+wtv2*(fva*vint(232)+fvb*vint(281))
          xpq(-kfv2)=xpq(-kfv2)+wtv2*(fva*vint(232)+fvb*vint(281))
        endif
        xpq(21)=xpq(0)
        mint(92)=4
      endif
 
C...Format for error printout.
 5000 format(' Error: x value outside physical range; x =',1p,d12.3)
 
      return
      end
 
C*********************************************************************
 
C...PYPDEL
C...Gives electron parton distribution.
 
      subroutine pjpdel(x,q2,xpel)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      save /jydat1/,/jydat2/,/pjpars/,/pjint1/
C...Local arrays.
      dimension xpel(-25:25),xpga(-6:6),sxp(0:6)
 
C...Interface to PDFLIB.
      common/w50513/xmin,xmax,q2min,q2max
      save /w50513/
      double precision xx,qq,upv,dnv,usea,dsea,str,chm,bot,top,glu,
     &value(20),xmin,xmax,q2min,q2max
      character*20 parm(20)
      data value/20*0d0/,parm/20*' '/
 
C...Some common constants.
      do 100 kfl=-25,25
        xpel(kfl)=0d0
  100 continue
      aem=paru(101)
      pme=pmas(11,1)
      xl=log(max(1d-10,x))
      x1l=log(max(1d-10,1d0-x))
      hle=log(max(3d0,q2/pme**2))
      hbe2=(aem/paru(1))*(hle-1d0)
 
C...Electron inside electron, see R. Kleiss et al., in Z physics at
C...LEP 1, CERN 89-08, p. 34
      if(mstp(59).le.1) then
        hde=1d0+(aem/paru(1))*(1.5d0*hle+1.289868d0)+(aem/paru(1))**2*
     &  (-2.164868d0*hle**2+9.840808d0*hle-10.130464d0)
        hee=hbe2*(1d0-x)**(hbe2-1d0)*sqrt(max(0d0,hde))-
     &  0.5d0*hbe2*(1d0+x)+hbe2**2/8d0*((1d0+x)*(-4d0*x1l+3d0*xl)-
     &  4d0*xl/(1d0-x)-5d0-x)
      else
        hee=hbe2*(1d0-x)**(hbe2-1d0)*exp(0.172784d0*hbe2)/
     &  pjgamm(1d0+hbe2)-0.5d0*hbe2*(1d0+x)+hbe2**2/8d0*((1d0+x)*
     &  (-4d0*x1l+3d0*xl)-4d0*xl/(1d0-x)-5d0-x)
      endif
      if(x.gt.0.9999d0.and.x.le.0.999999d0) then
        hee=hee*100d0**hbe2/(100d0**hbe2-1d0)
      elseif(x.gt.0.999999d0) then
        hee=0d0
      endif
      xpel(11)=x*hee
 
C...Photon and (transverse) W- inside electron.
      aemp=pjalem(pme*sqrt(max(0d0,q2)))/paru(2)
      if(mstp(13).le.1) then
        hlg=hle
      else
        hlg=log(max(1d0,(parp(13)/pme**2)*(1d0-x)/x**2))
      endif
      xpel(22)=aemp*hlg*(1d0+(1d0-x)**2)
      hlw=log(1d0+q2/pmas(24,1)**2)/(4d0*paru(102))
      xpel(-24)=aemp*hlw*(1d0+(1d0-x)**2)
 
C...Electron or positron inside photon inside electron.
      if(mstp(12).eq.1) then
        xfsea=0.5d0*(aemp*(hle-1d0))**2*(4d0/3d0+x-x**2-4d0*x**3/3d0+
     &  2d0*x*(1d0+x)*xl)
        xpel(11)=xpel(11)+xfsea
        xpel(-11)=xfsea
 
C...Initialize PDFLIB photon parton distributions.
        if(mstp(56).eq.2) then
          parm(1)='NPTYPE'
          value(1)=3
          parm(2)='NGROUP'
          value(2)=mstp(55)/1000
          parm(3)='NSET'
          value(3)=mod(mstp(55),1000)
          if(mint(93).ne.3000000+mstp(55)) then
            call pjpdfset(parm,value)
            mint(93)=3000000+mstp(55)
          endif
        endif
 
C...Quarks and gluons inside photon inside electron:
C...numerical convolution required.
        do 110 kfl=0,6
          sxp(kfl)=0d0
  110   continue
        sumxpp=0d0
        iter=-1
  120   iter=iter+1
        sumxp=sumxpp
        nstp=2**(iter-1)
        if(iter.eq.0) nstp=2
        do 130 kfl=0,6
          sxp(kfl)=0.5d0*sxp(kfl)
  130   continue
        wtstp=0.5d0/nstp
        if(iter.eq.0) wtstp=0.5d0
C...Pick grid of x_{gamma} values logarithmically even.
        do 150 istp=1,nstp
          if(iter.eq.0) then
            xle=xl*(istp-1)
          else
            xle=xl*(istp-0.5d0)/nstp
          endif
          xe=min(0.999999d0,exp(xle))
          xg=min(0.999999d0,x/xe)
C...Evaluate photon inside electron parton distribution for convolution.
          xpgp=1d0+(1d0-xe)**2
          if(mstp(13).le.1) then
            xpgp=xpgp*hle
          else
            xpgp=xpgp*log(max(1d0,(parp(13)/pme**2)*(1d0-xe)/xe**2))
          endif
C...Evaluate photon parton distributions for convolution.
          if(mstp(56).eq.1) then
            call pjpdga(xg,q2,xpga)
            do 140 kfl=0,5
              sxp(kfl)=sxp(kfl)+wtstp*xpgp*xpga(kfl)
  140       continue
          elseif(mstp(56).eq.2) then
C...Call PDFLIB parton distributions.
            xx=xg
            qq=sqrt(max(0d0,q2min,q2))
            if(mstp(57).eq.0) qq=sqrt(q2min)
            call pjstructm(xx,qq,upv,dnv,usea,dsea,str,chm,bot,top,glu)
            sxp(0)=sxp(0)+wtstp*xpgp*glu
            sxp(1)=sxp(1)+wtstp*xpgp*dnv
            sxp(2)=sxp(2)+wtstp*xpgp*upv
            sxp(3)=sxp(3)+wtstp*xpgp*str
            sxp(4)=sxp(4)+wtstp*xpgp*chm
            sxp(5)=sxp(5)+wtstp*xpgp*bot
            sxp(6)=sxp(6)+wtstp*xpgp*top
          endif
  150   continue
        sumxpp=sxp(0)+2d0*sxp(1)+2d0*sxp(2)
        if(iter.le.2.or.(iter.le.7.and.abs(sumxpp-sumxp).gt.
     &  parp(14)*(sumxpp+sumxp))) goto 120
 
C...Put convolution into output arrays.
        fconv=aemp*(-xl)
        xpel(0)=fconv*sxp(0)
        do 160 kfl=1,6
          xpel(kfl)=fconv*sxp(kfl)
          xpel(-kfl)=xpel(kfl)
  160   continue
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYPDGA
C...Gives photon parton distribution.
 
      subroutine pjpdga(x,q2,xpga)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      save /jydat1/,/pjpars/,/pjint1/
C...Local arrays.
      dimension xpga(-6:6),dgag(4,3),dgbg(4,3),dgcg(4,3),dgan(4,3),
     &dgbn(4,3),dgcn(4,3),dgdn(4,3),dgen(4,3),dgas(4,3),dgbs(4,3),
     &dgcs(4,3),dgds(4,3),dges(4,3)
 
C...The following data lines are coefficients needed in the
C...Drees and Grassie photon parton distribution parametrization.
      data dgag/-.207d0,.6158d0,1.074d0,0.d0,.8926d-2,.6594d0,
     &.4766d0,.1975d-1,.03197d0,1.018d0,.2461d0,.2707d-1/
      data dgbg/-.1987d0,.6257d0,8.352d0,5.024d0,.5085d-1,.2774d0,
     &-.3906d0,-.3212d0,-.618d-2,.9476d0,-.6094d0,-.1067d-1/
      data dgcg/5.119d0,-.2752d0,-6.993d0,2.298d0,-.2313d0,.1382d0,
     &6.542d0,.5162d0,-.1216d0,.9047d0,2.653d0,.2003d-2/
      data dgan/2.285d0,-.1526d-1,1330.d0,4.219d0,-.3711d0,1.061d0,
     &4.758d0,-.1503d-1,15.8d0,-.9464d0,-.5d0,-.2118d0/
      data dgbn/6.073d0,-.8132d0,-41.31d0,3.165d0,-.1717d0,.7815d0,
     &1.535d0,.7067d-2,2.742d0,-.7332d0,.7148d0,3.287d0/
      data dgcn/-.4202d0,.1778d-1,.9216d0,.18d0,.8766d-1,.2197d-1,
     &.1096d0,.204d0,.2917d-1,.4657d-1,.1785d0,.4811d-1/
      data dgdn/-.8083d-1,.6346d0,1.208d0,.203d0,-.8915d0,.2857d0,
     &2.973d0,.1185d0,-.342d-1,.7196d0,.7338d0,.8139d-1/
      data dgen/.5526d-1,1.136d0,.9512d0,.1163d-1,-.1816d0,.5866d0,
     &2.421d0,.4059d0,-.2302d-1,.9229d0,.5873d0,-.79d-4/
      data dgas/16.69d0,-.7916d0,1099.d0,4.428d0,-.1207d0,1.071d0,
     &1.977d0,-.8625d-2,6.734d0,-1.008d0,-.8594d-1,.7625d-1/
      data dgbs/.176d0,.4794d-1,1.047d0,.25d-1,25.d0,-1.648d0,
     &-.1563d-1,6.438d0,59.88d0,-2.983d0,4.48d0,.9686d0/
      data dgcs/-.208d-1,.3386d-2,4.853d0,.8404d0,-.123d-1,1.162d0,
     &.4824d0,-.11d-1,-.3226d-2,.8432d0,.3616d0,.1383d-2/
      data dgds/-.1685d-1,1.353d0,1.426d0,1.239d0,-.9194d-1,.7912d0,
     &.6397d0,2.327d0,-.3321d-1,.9475d0,-.3198d0,.2132d-1/
      data dges/-.1986d0,1.1d0,1.136d0,-.2779d0,.2015d-1,.9869d0,
     &-.7036d-1,.1694d-1,.1059d0,.6954d0,-.6663d0,.3683d0/
 
C...Photon parton distribution from Drees and Grassie.
C...Allowed variable range: 1 GeV^2 < Q^2 < 10000 GeV^2.
      do 100 kfl=-6,6
        xpga(kfl)=0d0
  100 continue
      vint(231)=1d0
      if(mstp(57).le.0) then
        t=log(1d0/0.16d0)
      else
        t=log(min(1d4,max(1d0,q2))/0.16d0)
      endif
      x1=1d0-x
      nf=3
      if(q2.gt.25d0) nf=4
      if(q2.gt.300d0) nf=5
      nfe=nf-2
      aem=paru(101)
 
C...Evaluate gluon content.
      dga=dgag(1,nfe)*t**dgag(2,nfe)+dgag(3,nfe)*t**(-dgag(4,nfe))
      dgb=dgbg(1,nfe)*t**dgbg(2,nfe)+dgbg(3,nfe)*t**(-dgbg(4,nfe))
      dgc=dgcg(1,nfe)*t**dgcg(2,nfe)+dgcg(3,nfe)*t**(-dgcg(4,nfe))
      xpgl=dga*x**dgb*x1**dgc
 
C...Evaluate up- and down-type quark content.
      dga=dgan(1,nfe)*t**dgan(2,nfe)+dgan(3,nfe)*t**(-dgan(4,nfe))
      dgb=dgbn(1,nfe)*t**dgbn(2,nfe)+dgbn(3,nfe)*t**(-dgbn(4,nfe))
      dgc=dgcn(1,nfe)*t**dgcn(2,nfe)+dgcn(3,nfe)*t**(-dgcn(4,nfe))
      dgd=dgdn(1,nfe)*t**dgdn(2,nfe)+dgdn(3,nfe)*t**(-dgdn(4,nfe))
      dge=dgen(1,nfe)*t**dgen(2,nfe)+dgen(3,nfe)*t**(-dgen(4,nfe))
      xpqn=x*(x**2+x1**2)/(dga-dgb*log(x1))+dgc*x**dgd*x1**dge
      dga=dgas(1,nfe)*t**dgas(2,nfe)+dgas(3,nfe)*t**(-dgas(4,nfe))
      dgb=dgbs(1,nfe)*t**dgbs(2,nfe)+dgbs(3,nfe)*t**(-dgbs(4,nfe))
      dgc=dgcs(1,nfe)*t**dgcs(2,nfe)+dgcs(3,nfe)*t**(-dgcs(4,nfe))
      dgd=dgds(1,nfe)*t**dgds(2,nfe)+dgds(3,nfe)*t**(-dgds(4,nfe))
      dge=dges(1,nfe)*t**dges(2,nfe)+dges(3,nfe)*t**(-dges(4,nfe))
      dgf=9d0
      if(nf.eq.4) dgf=10d0
      if(nf.eq.5) dgf=55d0/6d0
      xpqs=dgf*x*(x**2+x1**2)/(dga-dgb*log(x1))+dgc*x**dgd*x1**dge
      if(nf.le.3) then
        xpqu=(xpqs+9d0*xpqn)/6d0
        xpqd=(xpqs-4.5d0*xpqn)/6d0
      elseif(nf.eq.4) then
        xpqu=(xpqs+6d0*xpqn)/8d0
        xpqd=(xpqs-6d0*xpqn)/8d0
      else
        xpqu=(xpqs+7.5d0*xpqn)/10d0
        xpqd=(xpqs-5d0*xpqn)/10d0
      endif
 
C...Put into output arrays.
      xpga(0)=aem*xpgl
      xpga(1)=aem*xpqd
      xpga(2)=aem*xpqu
      xpga(3)=aem*xpqd
      if(nf.ge.4) xpga(4)=aem*xpqu
      if(nf.ge.5) xpga(5)=aem*xpqd
      do 110 kfl=1,6
        xpga(-kfl)=xpga(kfl)
  110 continue
 
      return
      end
 
C*********************************************************************
 
C...PYGGAM
C...Constructs the F2 and parton distributions of the photon
C...by summing homogeneous (VMD) and inhomogeneous (anomalous) terms.
C...For F2, c and b are included by the Bethe-Heitler formula;
C...in the 'MSbar' scheme additionally a Cgamma term is added.
C...Contains the SaS sets 1D, 1M, 2D and 2M.
C...Adapted from SaSgam library, authors G.A. Schuler and T. Sjostrand.
 
      subroutine pjggam(iset,x,q2,p2,ip2,f2gm,xpdfgm)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/pjint8/xpvmd(-6:6),xpanl(-6:6),xpanh(-6:6),xpbeh(-6:6),
     &xpdir(-6:6)
      common/pjint9/vxpvmd(-6:6),vxpanl(-6:6),vxpanh(-6:6),vxpdgm(-6:6)
      save /pjint8/,/pjint9/
C...Local arrays.
      dimension xpdfgm(-6:6),xpga(-6:6), vxpga(-6:6)
C...Charm and bottom masses (low to compensate for J/psi etc.).
      data pmc/1.3d0/, pmb/4.6d0/
C...alpha_em and alpha_em/(2*pi).
      data aem/0.007297d0/, aem2pi/0.0011614d0/
C...Lambda value for 4 flavours.
      data alam/0.20d0/
C...Mixture u/(u+d), = 0.5 for incoherent and = 0.8 for coherent sum.
      data fracu/0.8d0/
C...VMD couplings f_V**2/(4*pi).
      data frho/2.20d0/, fomega/23.6d0/, fphi/18.4d0/
C...Masses for rho (=omega) and phi.
      data pmrho/0.770d0/, pmphi/1.020d0/
C...Number of points in integration for IP2=1.
      data nstep/100/
 
C...Reset output.
      f2gm=0d0
      do 100 kfl=-6,6
        xpdfgm(kfl)=0d0
        xpvmd(kfl)=0d0
        xpanl(kfl)=0d0
        xpanh(kfl)=0d0
        xpbeh(kfl)=0d0
        xpdir(kfl)=0d0
        vxpvmd(kfl)=0d0
        vxpanl(kfl)=0d0
        vxpanh(kfl)=0d0
        vxpdgm(kfl)=0d0
  100 continue
 
C...Set Q0 cut-off parameter as function of set used.
      if(iset.le.2) then
        q0=0.6d0
      else
        q0=2d0
      endif
      q02=q0**2
 
C...Scale choice for off-shell photon; common factors.
      q2a=q2
      facnor=1d0
      if(ip2.eq.1) then
        p2mx=p2+q02
        q2a=q2+p2*q02/max(q02,q2)
        facnor=log(q2/q02)/nstep
      elseif(ip2.eq.2) then
        p2mx=max(p2,q02)
      elseif(ip2.eq.3) then
        p2mx=p2+q02
        q2a=q2+p2*q02/max(q02,q2)
      elseif(ip2.eq.4) then
        p2mx=q2*(q02+p2)/(q2+p2)*exp(p2*(q2-q02)/
     &  ((q2+p2)*(q02+p2)))
      elseif(ip2.eq.5) then
        p2mxa=q2*(q02+p2)/(q2+p2)*exp(p2*(q2-q02)/
     &  ((q2+p2)*(q02+p2)))
        p2mx=q0*sqrt(p2mxa)
        facnor=log(q2/p2mxa)/log(q2/p2mx)
      elseif(ip2.eq.6) then
        p2mx=q2*(q02+p2)/(q2+p2)*exp(p2*(q2-q02)/
     &  ((q2+p2)*(q02+p2)))
        p2mx=max(0d0,1d0-p2/q2)*p2mx+min(1d0,p2/q2)*max(p2,q02)
      else
        p2mxa=q2*(q02+p2)/(q2+p2)*exp(p2*(q2-q02)/
     &  ((q2+p2)*(q02+p2)))
        p2mx=q0*sqrt(p2mxa)
        p2mxb=p2mx
        p2mx=max(0d0,1d0-p2/q2)*p2mx+min(1d0,p2/q2)*max(p2,q02)
        p2mxb=max(0d0,1d0-p2/q2)*p2mxb+min(1d0,p2/q2)*p2mxa
        facnor=log(q2/p2mxa)/log(q2/p2mxb)
      endif
 
C...Call VMD parametrization for d quark and use to give rho, omega,
C...phi. Note dipole dampening for off-shell photon.
      call pjgvmd(iset,1,x,q2a,p2mx,alam,xpga,vxpga)
      xfval=vxpga(1)
      xpga(1)=xpga(2)
      xpga(-1)=xpga(-2)
      facud=aem*(1d0/frho+1d0/fomega)*(pmrho**2/(pmrho**2+p2))**2
      facs=aem*(1d0/fphi)*(pmphi**2/(pmphi**2+p2))**2
      do 110 kfl=-5,5
        xpvmd(kfl)=(facud+facs)*xpga(kfl)
  110 continue
      xpvmd(1)=xpvmd(1)+(1d0-fracu)*facud*xfval
      xpvmd(2)=xpvmd(2)+fracu*facud*xfval
      xpvmd(3)=xpvmd(3)+facs*xfval
      xpvmd(-1)=xpvmd(-1)+(1d0-fracu)*facud*xfval
      xpvmd(-2)=xpvmd(-2)+fracu*facud*xfval
      xpvmd(-3)=xpvmd(-3)+facs*xfval
      vxpvmd(1)=(1d0-fracu)*facud*xfval
      vxpvmd(2)=fracu*facud*xfval
      vxpvmd(3)=facs*xfval
      vxpvmd(-1)=(1d0-fracu)*facud*xfval
      vxpvmd(-2)=fracu*facud*xfval
      vxpvmd(-3)=facs*xfval
 
      if(ip2.ne.1) then
C...Anomalous parametrizations for different strategies
C...for off-shell photons; except full integration.
 
C...Call anomalous parametrization for d + u + s.
        call pjgano(-3,x,q2a,p2mx,alam,xpga,vxpga)
        do 120 kfl=-5,5
          xpanl(kfl)=facnor*xpga(kfl)
          vxpanl(kfl)=facnor*vxpga(kfl)
  120   continue
 
C...Call anomalous parametrization for c and b.
        call pjgano(4,x,q2a,p2mx,alam,xpga,vxpga)
        do 130 kfl=-5,5
          xpanh(kfl)=facnor*xpga(kfl)
          vxpanh(kfl)=facnor*vxpga(kfl)
  130   continue
        call pjgano(5,x,q2a,p2mx,alam,xpga,vxpga)
        do 140 kfl=-5,5
          xpanh(kfl)=xpanh(kfl)+facnor*xpga(kfl)
          vxpanh(kfl)=vxpanh(kfl)+facnor*vxpga(kfl)
  140   continue
 
      else
C...Special option: loop over flavours and integrate over k2.
        do 170 kf=1,5
          do 160 istep=1,nstep
            q2step=q02*(q2/q02)**((istep-0.5d0)/nstep)
            if((kf.eq.4.and.q2step.lt.pmc**2).or.
     &      (kf.eq.5.and.q2step.lt.pmb**2)) goto 160
            call pjgvmd(0,kf,x,q2,q2step,alam,xpga,vxpga)
            facq=aem2pi*(q2step/(q2step+p2))**2*facnor
            if(mod(kf,2).eq.0) facq=facq*(8d0/9d0)
            if(mod(kf,2).eq.1) facq=facq*(2d0/9d0)
            do 150 kfl=-5,5
              if(kf.le.3) xpanl(kfl)=xpanl(kfl)+facq*xpga(kfl)
              if(kf.ge.4) xpanh(kfl)=xpanh(kfl)+facq*xpga(kfl)
              if(kf.le.3) vxpanl(kfl)=vxpanl(kfl)+facq*vxpga(kfl)
              if(kf.ge.4) vxpanh(kfl)=vxpanh(kfl)+facq*vxpga(kfl)
  150       continue
  160     continue
  170   continue
      endif
 
C...Call Bethe-Heitler term expression for charm and bottom.
      call pjgbeh(4,x,q2,p2,pmc**2,xpbh)
      xpbeh(4)=xpbh
      xpbeh(-4)=xpbh
      call pjgbeh(5,x,q2,p2,pmb**2,xpbh)
      xpbeh(5)=xpbh
      xpbeh(-5)=xpbh
 
C...For MSbar subtraction call C^gamma term expression for d, u, s.
      if(iset.eq.2.or.iset.eq.4) then
        call pjgdir(x,q2,p2,q02,xpga)
        do 180 kfl=-5,5
          xpdir(kfl)=xpga(kfl)
  180   continue
      endif
 
C...Store result in output array.
      do 190 kfl=-5,5
        chsq=1d0/9d0
        if(iabs(kfl).eq.2.or.iabs(kfl).eq.4) chsq=4d0/9d0
        xpf2=xpvmd(kfl)+xpanl(kfl)+xpbeh(kfl)+xpdir(kfl)
        if(kfl.ne.0) f2gm=f2gm+chsq*xpf2
        xpdfgm(kfl)=xpvmd(kfl)+xpanl(kfl)+xpanh(kfl)
        vxpdgm(kfl)=vxpvmd(kfl)+vxpanl(kfl)+vxpanh(kfl)
  190 continue
 
      return
      end
 
C*********************************************************************
 
C...PYGVMD
C...Evaluates the VMD parton distributions of a photon,
C...evolved homogeneously from an initial scale P2 to Q2.
C...Does not include dipole suppression factor.
C...ISET is parton distribution set, see above;
C...additionally ISET=0 is used for the evolution of an anomalous photon
C...which branched at a scale P2 and then evolved homogeneously to Q2.
C...ALAM is the 4-flavour Lambda, which is automatically converted
C...to 3- and 5-flavour equivalents as needed.
C...Adapted from SaSgam library, authors G.A. Schuler and T. Sjostrand.
 
      subroutine pjgvmd(iset,kf,x,q2,p2,alam,xpga,vxpga)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Local arrays and data.
      dimension xpga(-6:6), vxpga(-6:6)
      data pmc/1.3d0/, pmb/4.6d0/
c    $  ,aem/0.007297d0/, aem2pi/0.0011614d0/
 
C...Reset output.
      do 100 kfl=-6,6
        xpga(kfl)=0d0
        vxpga(kfl)=0d0
  100 continue
      kfa=iabs(kf)
 
C...Calculate Lambda; protect against unphysical Q2 and P2 input.
      alam3=alam*(pmc/alam)**(2d0/27d0)
      alam5=alam*(alam/pmb)**(2d0/23d0)
      p2eff=max(p2,1.2d0*alam3**2)
      if(kfa.eq.4) p2eff=max(p2eff,pmc**2)
      if(kfa.eq.5) p2eff=max(p2eff,pmb**2)
      q2eff=max(q2,p2eff)
 
C...Find number of flavours at lower and upper scale.
      nfp=4
      if(p2eff.lt.pmc**2) nfp=3
      if(p2eff.gt.pmb**2) nfp=5
      nfq=4
      if(q2eff.lt.pmc**2) nfq=3
      if(q2eff.gt.pmb**2) nfq=5
 
C...Find s as sum of 3-, 4- and 5-flavour parts.
      s=0d0
      if(nfp.eq.3) then
        q2div=pmc**2
        if(nfq.eq.3) q2div=q2eff
        s=s+(6d0/27d0)*log(log(q2div/alam3**2)/log(p2eff/alam3**2))
      endif
      if(nfp.le.4.and.nfq.ge.4) then
        p2div=p2eff
        if(nfp.eq.3) p2div=pmc**2
        q2div=q2eff
        if(nfq.eq.5) q2div=pmb**2
        s=s+(6d0/25d0)*log(log(q2div/alam**2)/log(p2div/alam**2))
      endif
      if(nfq.eq.5) then
        p2div=pmb**2
        if(nfp.eq.5) p2div=p2eff
        s=s+(6d0/23d0)*log(log(q2eff/alam5**2)/log(p2div/alam5**2))
      endif
 
C...Calculate frequent combinations of x and s.
      x1=1d0-x
      xl=-log(x)
      s2=s**2
      s3=s**3
      s4=s**4
 
C...Evaluate homogeneous anomalous parton distributions below or
C...above threshold.
      if(iset.eq.0) then
        if(q2.le.p2.or.(kfa.eq.4.and.q2.lt.pmc**2).or.
     &  (kfa.eq.5.and.q2.lt.pmb**2)) then
          xval = x * 1.5d0 * (x**2+x1**2)
          xglu = 0d0
          xsea = 0d0
        else
          xval = (1.5d0/(1d0-0.197d0*s+4.33d0*s2)*x**2 +
     &    (1.5d0+2.10d0*s)/(1d0+3.29d0*s)*x1**2 +
     &    5.23d0*s/(1d0+1.17d0*s+19.9d0*s3)*x*x1) *
     &    x**(1d0/(1d0+1.5d0*s)) * (1d0-x**2)**(2.667d0*s)
          xglu = 4d0*s/(1d0+4.76d0*s+15.2d0*s2+29.3d0*s4) *
     &    x**(-2.03d0*s/(1d0+2.44d0*s)) * (x1*xl)**(1.333d0*s) *
     &    ((4d0*x**2+7d0*x+4d0)*x1/3d0 - 2d0*x*(1d0+x)*xl)
          xsea = s2/(1d0+4.54d0*s+8.19d0*s2+8.05d0*s3) *
     &    x**(-1.54d0*s/(1d0+1.29d0*s)) * x1**(2.667d0*s) *
     &    ((8d0-73d0*x+62d0*x**2)*x1/9d0 + (3d0-8d0*x**2/3d0)*x*xl +
     &    (2d0*x-1d0)*x*xl**2)
        endif
 
C...Evaluate set 1D parton distributions below or above threshold.
      elseif(iset.eq.1) then
        if(q2.le.p2.or.(kfa.eq.4.and.q2.lt.pmc**2).or.
     &  (kfa.eq.5.and.q2.lt.pmb**2)) then
          xval = 1.294d0 * x**0.80d0 * x1**0.76d0
          xglu = 1.273d0 * x**0.40d0 * x1**1.76d0
          xsea = 0.100d0 * x1**3.76d0
        else
          xval = 1.294d0/(1d0+0.252d0*s+3.079d0*s2) *
     &    x**(0.80d0-0.13d0*s) * x1**(0.76d0+0.667d0*s) * xl**(2d0*s)
          xglu = 7.90d0*s/(1d0+5.50d0*s) * exp(-5.16d0*s) *
     &    x**(-1.90d0*s/(1d0+3.60d0*s)) * x1**1.30d0 *
     &    xl**(0.50d0+3d0*s) + 1.273d0 * exp(-10d0*s) *
     &    x**0.40d0 * x1**(1.76d0+3d0*s)
          xsea = (0.1d0-0.397d0*s2+1.121d0*s3)/
     &    (1d0+5.61d0*s2+5.26d0*s3) * x**(-7.32d0*s2/(1d0+10.3d0*s2)) *
     &    x1**((3.76d0+15d0*s+12d0*s2)/(1d0+4d0*s))
          xsea0 = 0.100d0 * x1**3.76d0
        endif
 
C...Evaluate set 1M parton distributions below or above threshold.
      elseif(iset.eq.2) then
        if(q2.le.p2.or.(kfa.eq.4.and.q2.lt.pmc**2).or.
     &  (kfa.eq.5.and.q2.lt.pmb**2)) then
          xval = 0.8477d0 * x**0.51d0 * x1**1.37d0
          xglu = 3.42d0 * x**0.255d0 * x1**2.37d0
          xsea = 0d0
        else
          xval = 0.8477d0/(1d0+1.37d0*s+2.18d0*s2+3.73d0*s3) *
     &    x**(0.51d0+0.21d0*s) * x1**1.37d0 * xl**(2.667d0*s)
          xglu = 24d0*s/(1d0+9.6d0*s+0.92d0*s2+14.34d0*s3) *
     &    exp(-5.94d0*s) * x**((-0.013d0-1.80d0*s)/(1d0+3.14d0*s)) *
     &    x1**(2.37d0+0.4d0*s) * xl**(0.32d0+3.6d0*s) + 3.42d0 *
     &    exp(-12d0*s) * x**0.255d0 * x1**(2.37d0+3d0*s)
          xsea = 0.842d0*s/(1d0+21.3d0*s-33.2d0*s2+229d0*s3) *
     &    x**((0.13d0-2.90d0*s)/(1d0+5.44d0*s)) * x1**(3.45d0+0.5d0*s) *
     &    xl**(2.8d0*s)
          xsea0 = 0d0
        endif
 
C...Evaluate set 2D parton distributions below or above threshold.
      elseif(iset.eq.3) then
        if(q2.le.p2.or.(kfa.eq.4.and.q2.lt.pmc**2).or.
     &  (kfa.eq.5.and.q2.lt.pmb**2)) then
          xval = x**0.46d0 * x1**0.64d0 + 0.76d0 * x
          xglu = 1.925d0 * x1**2
          xsea = 0.242d0 * x1**4
        else
          xval = (1d0+0.186d0*s)/(1d0-0.209d0*s+1.495d0*s2) *
     &    x**(0.46d0+0.25d0*s) *
     &    x1**((0.64d0+0.14d0*s+5d0*s2)/(1d0+s)) * xl**(1.9d0*s) +
     &    (0.76d0+0.4d0*s) * x * x1**(2.667d0*s)
          xglu = (1.925d0+5.55d0*s+147d0*s2)/(1d0-3.59d0*s+3.32d0*s2) *
     &    exp(-18.67d0*s) *
     &    x**((-5.81d0*s-5.34d0*s2)/(1d0+29d0*s-4.26d0*s2))
     &    * x1**((2d0-5.9d0*s)/(1d0+1.7d0*s)) *
     &    xl**(9.3d0*s/(1d0+1.7d0*s))
          xsea = (0.242d0-0.252d0*s+1.19d0*s2)/
     &    (1d0-0.607d0*s+21.95d0*s2) *
     &    x**(-12.1d0*s2/(1d0+2.62d0*s+16.7d0*s2)) * x1**4 * xl**s
          xsea0 = 0.242d0 * x1**4
        endif
 
C...Evaluate set 2M parton distributions below or above threshold.
      elseif(iset.eq.4) then
        if(q2.le.p2.or.(kfa.eq.4.and.q2.lt.pmc**2).or.
     &  (kfa.eq.5.and.q2.lt.pmb**2)) then
          xval = 1.168d0 * x**0.50d0 * x1**2.60d0 + 0.965d0 * x
          xglu = 1.808d0 * x1**2
          xsea = 0.209d0 * x1**4
        else
          xval = (1.168d0+1.771d0*s+29.35d0*s2) * exp(-5.776d0*s) *
     &    x**((0.5d0+0.208d0*s)/(1d0-0.794d0*s+1.516d0*s2)) *
     &    x1**((2.6d0+7.6d0*s)/(1d0+5d0*s)) *
     &    xl**(5.15d0*s/(1d0+2d0*s)) +
     &    (0.965d0+22.35d0*s)/(1d0+18.4d0*s) * x * x1**(2.667d0*s)
          xglu = (1.808d0+29.9d0*s)/(1d0+26.4d0*s) * exp(-5.28d0*s) *
     &    x**((-5.35d0*s-10.11d0*s2)/(1d0+31.71d0*s)) *
     &    x1**((2d0-7.3d0*s+4d0*s2)/(1d0+2.5d0*s)) *
     &    xl**(10.9d0*s/(1d0+2.5d0*s))
          xsea = (0.209d0+0.644d0*s2)/(1d0+0.319d0*s+17.6d0*s2) *
     &    x**((-0.373d0*s-7.71d0*s2)/(1d0+0.815d0*s+11.0d0*s2)) *
     &    x1**(4d0+s) * xl**(0.45d0*s)
          xsea0 = 0.209d0 * x1**4
        endif
      endif
 
C...Threshold factors for c and b sea.
      sll=log(log(q2eff/alam**2)/log(p2eff/alam**2))
      xchm=0d0
      if(q2.gt.pmc**2.and.q2.gt.1.001d0*p2eff) then
        sch=max(0d0,log(log(pmc**2/alam**2)/log(p2eff/alam**2)))
        if(iset.eq.0) then
          xchm=xsea*(1d0-(sch/sll)**2)
        else
          xchm=max(0d0,xsea-xsea0*x1**(2.667d0*s))*(1d0-sch/sll)
        endif
      endif
      xbot=0d0
      if(q2.gt.pmb**2.and.q2.gt.1.001d0*p2eff) then
        sbt=max(0d0,log(log(pmb**2/alam**2)/log(p2eff/alam**2)))
        if(iset.eq.0) then
          xbot=xsea*(1d0-(sbt/sll)**2)
        else
          xbot=max(0d0,xsea-xsea0*x1**(2.667d0*s))*(1d0-sbt/sll)
        endif
      endif
 
C...Fill parton distributions.
      xpga(0)=xglu
      xpga(1)=xsea
      xpga(2)=xsea
      xpga(3)=xsea
      xpga(4)=xchm
      xpga(5)=xbot
      xpga(kfa)=xpga(kfa)+xval
      do 110 kfl=1,5
        xpga(-kfl)=xpga(kfl)
  110 continue
      vxpga(kfa)=xval
      vxpga(-kfa)=xval
 
      return
      end
 
C*********************************************************************
 
C...PYGANO
C...Evaluates the parton distributions of the anomalous photon,
C...inhomogeneously evolved from a scale P2 (where it vanishes) to Q2.
C...KF=0 gives the sum over (up to) 5 flavours,
C...KF<0 limits to flavours up to abs(KF),
C...KF>0 is for flavour KF only.
C...ALAM is the 4-flavour Lambda, which is automatically converted
C...to 3- and 5-flavour equivalents as needed.
C...Adapted from SaSgam library, authors G.A. Schuler and T. Sjostrand.
 
      subroutine pjgano(kf,x,q2,p2,alam,xpga,vxpga)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Local arrays and data.
      dimension xpga(-6:6), vxpga(-6:6), alamsq(3:5)
      data pmc/1.3d0/, pmb/4.6d0/, aem2pi/0.0011614d0/
c, aem/0.007297d0/
 
C...Reset output.
      do 100 kfl=-6,6
        xpga(kfl)=0d0
        vxpga(kfl)=0d0
  100 continue
      if(q2.le.p2) return
      kfa=iabs(kf)
 
C...Calculate Lambda; protect against unphysical Q2 and P2 input.
      alamsq(3)=(alam*(pmc/alam)**(2d0/27d0))**2
      alamsq(4)=alam**2
      alamsq(5)=(alam*(alam/pmb)**(2d0/23d0))**2
      p2eff=max(p2,1.2d0*alamsq(3))
      if(kf.eq.4) p2eff=max(p2eff,pmc**2)
      if(kf.eq.5) p2eff=max(p2eff,pmb**2)
      q2eff=max(q2,p2eff)
      xl=-log(x)
 
C...Find number of flavours at lower and upper scale.
      nfp=4
      if(p2eff.lt.pmc**2) nfp=3
      if(p2eff.gt.pmb**2) nfp=5
      nfq=4
      if(q2eff.lt.pmc**2) nfq=3
      if(q2eff.gt.pmb**2) nfq=5
 
C...Define range of flavour loop.
      if(kf.eq.0) then
        kflmn=1
        kflmx=5
      elseif(kf.lt.0) then
        kflmn=1
        kflmx=kfa
      else
        kflmn=kfa
        kflmx=kfa
      endif
 
C...Loop over flavours the photon can branch into.
      do 110 kfl=kflmn,kflmx
 
C...Light flavours: calculate t range and (approximate) s range.
        if(kfl.le.3.and.(kfl.eq.1.or.kfl.eq.kf)) then
          tdiff=log(q2eff/p2eff)
          s=(6d0/(33d0-2d0*nfq))*log(log(q2eff/alamsq(nfq))/
     &    log(p2eff/alamsq(nfq)))
          if(nfq.gt.nfp) then
            q2div=pmb**2
            if(nfq.eq.4) q2div=pmc**2
            snfq=(6d0/(33d0-2d0*nfq))*log(log(q2div/alamsq(nfq))/
     &      log(p2eff/alamsq(nfq)))
            snfp=(6d0/(33d0-2d0*(nfq-1)))*log(log(q2div/alamsq(nfq-1))/
     &      log(p2eff/alamsq(nfq-1)))
            s=s+(log(q2div/p2eff)/log(q2eff/p2eff))*(snfp-snfq)
          endif
          if(nfq.eq.5.and.nfp.eq.3) then
            q2div=pmc**2
            snf4=(6d0/(33d0-2d0*4))*log(log(q2div/alamsq(4))/
     &      log(p2eff/alamsq(4)))
            snf3=(6d0/(33d0-2d0*3))*log(log(q2div/alamsq(3))/
     &      log(p2eff/alamsq(3)))
            s=s+(log(q2div/p2eff)/log(q2eff/p2eff))*(snf3-snf4)
          endif
 
C...u and s quark do not need a separate treatment when d has been done.
        elseif(kfl.eq.2.or.kfl.eq.3) then
 
C...Charm: as above, but only include range above c threshold.
        elseif(kfl.eq.4) then
          if(q2.le.pmc**2) goto 110
          p2eff=max(p2eff,pmc**2)
          q2eff=max(q2eff,p2eff)
          tdiff=log(q2eff/p2eff)
          s=(6d0/(33d0-2d0*nfq))*log(log(q2eff/alamsq(nfq))/
     &    log(p2eff/alamsq(nfq)))
          if(nfq.eq.5.and.nfp.eq.4) then
            q2div=pmb**2
            snfq=(6d0/(33d0-2d0*nfq))*log(log(q2div/alamsq(nfq))/
     &      log(p2eff/alamsq(nfq)))
            snfp=(6d0/(33d0-2d0*(nfq-1)))*log(log(q2div/alamsq(nfq-1))/
     &      log(p2eff/alamsq(nfq-1)))
            s=s+(log(q2div/p2eff)/log(q2eff/p2eff))*(snfp-snfq)
          endif
 
C...Bottom: as above, but only include range above b threshold.
        elseif(kfl.eq.5) then
          if(q2.le.pmb**2) goto 110
          p2eff=max(p2eff,pmb**2)
          q2eff=max(q2,p2eff)
          tdiff=log(q2eff/p2eff)
          s=(6d0/(33d0-2d0*nfq))*log(log(q2eff/alamsq(nfq))/
     &    log(p2eff/alamsq(nfq)))
        endif
 
C...Evaluate flavour-dependent prefactor (charge^2 etc.).
        chsq=1d0/9d0
        if(kfl.eq.2.or.kfl.eq.4) chsq=4d0/9d0
        fac=aem2pi*2d0*chsq*tdiff
 
C...Evaluate parton distributions (normalized to unit momentum sum).
        if(kfl.eq.1.or.kfl.eq.4.or.kfl.eq.5.or.kfl.eq.kf) then
          xval= ((1.5d0+2.49d0*s+26.9d0*s**2)/(1d0+32.3d0*s**2)*x**2 +
     &    (1.5d0-0.49d0*s+7.83d0*s**2)/(1d0+7.68d0*s**2)*(1d0-x)**2 +
     &    1.5d0*s/(1d0-3.2d0*s+7d0*s**2)*x*(1d0-x)) *
     &    x**(1d0/(1d0+0.58d0*s)) * (1d0-x**2)**(2.5d0*s/(1d0+10d0*s))
          xglu= 2d0*s/(1d0+4d0*s+7d0*s**2) *
     &    x**(-1.67d0*s/(1d0+2d0*s)) * (1d0-x**2)**(1.2d0*s) *
     &    ((4d0*x**2+7d0*x+4d0)*(1d0-x)/3d0 - 2d0*x*(1d0+x)*xl)
          xsea= 0.333d0*s**2/(1d0+4.90d0*s+4.69d0*s**2+21.4d0*s**3) *
     &    x**(-1.18d0*s/(1d0+1.22d0*s)) * (1d0-x)**(1.2d0*s) *
     &    ((8d0-73d0*x+62d0*x**2)*(1d0-x)/9d0 +
     &    (3d0-8d0*x**2/3d0)*x*xl + (2d0*x-1d0)*x*xl**2)
 
C...Threshold factors for c and b sea.
          sll=log(log(q2eff/alam**2)/log(p2eff/alam**2))
          xchm=0d0
          if(q2.gt.pmc**2.and.q2.gt.1.001d0*p2eff) then
            sch=max(0d0,log(log(pmc**2/alam**2)/log(p2eff/alam**2)))
            xchm=xsea*(1d0-(sch/sll)**3)
          endif
          xbot=0d0
          if(q2.gt.pmb**2.and.q2.gt.1.001d0*p2eff) then
            sbt=max(0d0,log(log(pmb**2/alam**2)/log(p2eff/alam**2)))
            xbot=xsea*(1d0-(sbt/sll)**3)
          endif
        endif
 
C...Add contribution of each valence flavour.
        xpga(0)=xpga(0)+fac*xglu
        xpga(1)=xpga(1)+fac*xsea
        xpga(2)=xpga(2)+fac*xsea
        xpga(3)=xpga(3)+fac*xsea
        xpga(4)=xpga(4)+fac*xchm
        xpga(5)=xpga(5)+fac*xbot
        xpga(kfl)=xpga(kfl)+fac*xval
        vxpga(kfl)=vxpga(kfl)+fac*xval
  110 continue
      do 120 kfl=1,5
        xpga(-kfl)=xpga(kfl)
        vxpga(-kfl)=vxpga(kfl)
  120 continue
 
      return
      end
 
C*********************************************************************
 
C...PYGBEH
C...Evaluates the Bethe-Heitler cross section for heavy flavour
C...production.
C...Adapted from SaSgam library, authors G.A. Schuler and T. Sjostrand.
 
      subroutine pjgbeh(kf,x,q2,p2,pm2,xpbh)
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
 
C...Local data.
      data aem2pi/0.0011614d0/
 
C...Reset output.
      xpbh=0d0
      sigbh=0d0
 
C...Check kinematics limits.
      if(x.ge.q2/(4d0*pm2+q2+p2)) return
      w2=q2*(1d0-x)/x-p2
      beta2=1d0-4d0*pm2/w2
      if(beta2.lt.1d-10) return
      beta=sqrt(beta2)
      rmq=4d0*pm2/q2
 
C...Simple case: P2 = 0.
      if(p2.lt.1d-4) then
        if(beta.lt.0.99d0) then
          xbl=log((1d0+beta)/(1d0-beta))
        else
          xbl=log((1d0+beta)**2*w2/(4d0*pm2))
        endif
        sigbh=beta*(8d0*x*(1d0-x)-1d0-rmq*x*(1d0-x))+
     &  xbl*(x**2+(1d0-x)**2+rmq*x*(1d0-3d0*x)-0.5d0*rmq**2*x**2)
 
C...Complicated case: P2 > 0, based on approximation of
C...C.T. Hill and G.G. Ross, Nucl. Phys. B148 (1979) 373
      else
        rpq=1d0-4d0*x**2*p2/q2
        if(rpq.gt.1d-10) then
          rpbe=sqrt(rpq*beta2)
          if(rpbe.lt.0.99d0) then
            xbl=log((1d0+rpbe)/(1d0-rpbe))
            xbi=2d0*rpbe/(1d0-rpbe**2)
          else
            rpbesn=4d0*pm2/w2+(4d0*x**2*p2/q2)*beta2
            xbl=log((1d0+rpbe)**2/rpbesn)
            xbi=2d0*rpbe/rpbesn
          endif
          sigbh=beta*(6d0*x*(1d0-x)-1d0)+
     &    xbl*(x**2+(1d0-x)**2+rmq*x*(1d0-3d0*x)-0.5d0*rmq**2*x**2)+
     &    xbi*(2d0*x/q2)*(pm2*x*(2d0-rmq)-p2*x)
        endif
      endif
 
C...Multiply by charge-squared etc. to get parton distribution.
      chsq=1d0/9d0
      if(iabs(kf).eq.2.or.iabs(kf).eq.4) chsq=4d0/9d0
      xpbh=3d0*chsq*aem2pi*x*sigbh
 
      return
      end
 
C*********************************************************************
 
C...PYGDIR
C...Evaluates the direct contribution, i.e. the C^gamma term,
C...as needed in MSbar parametrizations.
C...Adapted from SaSgam library, authors G.A. Schuler and T. Sjostrand.
 
      subroutine pjgdir(x,q2,p2,q02,xpga)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Local array and data.
      dimension xpga(-6:6)
c     data pmc/1.3d0/, pmb/4.6d0/, aem2pi/0.0011614d0/
      data aem2pi/0.0011614d0/
 
C...Reset output.
      do 100 kfl=-6,6
        xpga(kfl)=0d0
  100 continue
 
C...Evaluate common x-dependent expression.
      xtmp = (x**2+(1d0-x)**2) * (-log(x)) - 1d0
      cgam = 3d0*aem2pi*x * (xtmp*(1d0+p2/(p2+q02)) + 6d0*x*(1d0-x))
 
C...d, u, s part by simple charge factor.
      xpga(1)=(1d0/9d0)*cgam
      xpga(2)=(4d0/9d0)*cgam
      xpga(3)=(1d0/9d0)*cgam
 
C...Also fill for antiquarks.
      do 110 kf=1,5
        xpga(-kf)=xpga(kf)
  110 continue
 
      return
      end
 
C*********************************************************************
 
C...PYPDPI
C...Gives pi+ parton distribution according to two different
C...parametrizations.
 
      subroutine pjpdpi(x,q2,xppi)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      save /jydat1/,/pjpars/,/pjint1/
C...Local arrays.
      dimension xppi(-6:6),cow(3,5,4,2),xq(9),ts(6)
 
C...The following data lines are coefficients needed in the
C...Owens pion parton distribution parametrizations, see below.
C...Expansion coefficients for up and down valence quark distributions.
      data ((cow(ip,is,1,1),is=1,5),ip=1,3)/
     &4.0000d-01,  7.0000d-01,  0.0000d+00,  0.0000d+00,  0.0000d+00,
     &-6.2120d-02,  6.4780d-01,  0.0000d+00,  0.0000d+00,  0.0000d+00,
     &-7.1090d-03,  1.3350d-02,  0.0000d+00,  0.0000d+00,  0.0000d+00/
      data ((cow(ip,is,1,2),is=1,5),ip=1,3)/
     &4.0000d-01,  6.2800d-01,  0.0000d+00,  0.0000d+00,  0.0000d+00,
     &-5.9090d-02,  6.4360d-01,  0.0000d+00,  0.0000d+00,  0.0000d+00,
     &-6.5240d-03,  1.4510d-02,  0.0000d+00,  0.0000d+00,  0.0000d+00/
C...Expansion coefficients for gluon distribution.
      data ((cow(ip,is,2,1),is=1,5),ip=1,3)/
     &8.8800d-01,  0.0000d+00,  3.1100d+00,  6.0000d+00,  0.0000d+00,
     &-1.8020d+00, -1.5760d+00, -1.3170d-01,  2.8010d+00, -1.7280d+01,
     &1.8120d+00,  1.2000d+00,  5.0680d-01, -1.2160d+01,  2.0490d+01/
      data ((cow(ip,is,2,2),is=1,5),ip=1,3)/
     &7.9400d-01,  0.0000d+00,  2.8900d+00,  6.0000d+00,  0.0000d+00,
     &-9.1440d-01, -1.2370d+00,  5.9660d-01, -3.6710d+00, -8.1910d+00,
     &5.9660d-01,  6.5820d-01, -2.5500d-01, -2.3040d+00,  7.7580d+00/
C...Expansion coefficients for (up+down+strange) quark sea distribution.
      data ((cow(ip,is,3,1),is=1,5),ip=1,3)/
     &9.0000d-01,  0.0000d+00,  5.0000d+00,  0.0000d+00,  0.0000d+00,
     &-2.4280d-01, -2.1200d-01,  8.6730d-01,  1.2660d+00,  2.3820d+00,
     &1.3860d-01,  3.6710d-03,  4.7470d-02, -2.2150d+00,  3.4820d-01/
      data ((cow(ip,is,3,2),is=1,5),ip=1,3)/
     &9.0000d-01,  0.0000d+00,  5.0000d+00,  0.0000d+00,  0.0000d+00,
     &-1.4170d-01, -1.6970d-01, -2.4740d+00, -2.5340d+00,  5.6210d-01,
     &-1.7400d-01, -9.6230d-02,  1.5750d+00,  1.3780d+00, -2.7010d-01/
C...Expansion coefficients for charm quark sea distribution.
      data ((cow(ip,is,4,1),is=1,5),ip=1,3)/
     &0.0000d+00, -2.2120d-02,  2.8940d+00,  0.0000d+00,  0.0000d+00,
     &7.9280d-02, -3.7850d-01,  9.4330d+00,  5.2480d+00,  8.3880d+00,
     &-6.1340d-02, -1.0880d-01, -1.0852d+01, -7.1870d+00, -1.1610d+01/
      data ((cow(ip,is,4,2),is=1,5),ip=1,3)/
     &0.0000d+00, -8.8200d-02,  1.9240d+00,  0.0000d+00,  0.0000d+00,
     &6.2290d-02, -2.8920d-01,  2.4240d-01, -4.4630d+00, -8.3670d-01,
     &-4.0990d-02, -1.0820d-01,  2.0360d+00,  5.2090d+00, -4.8400d-02/
 
C...Euler's beta function, requires ordinary Gamma function
      eulbet(x,y)=pjgamm(x)*pjgamm(y)/pjgamm(x+y)
 
C...Reset output array.
      do 100 kfl=-6,6
        xppi(kfl)=0d0
  100 continue
 
      if(mstp(53).le.2) then
C...Pion parton distributions from Owens.
C...Allowed variable range: 4 GeV^2 < Q^2 < approx 2000 GeV^2.
 
C...Determine set, Lambda and s expansion variable.
        nset=mstp(53)
        if(nset.eq.1) alam=0.2d0
        if(nset.eq.2) alam=0.4d0
        vint(231)=4d0
        if(mstp(57).le.0) then
          sd=0d0
        else
          q2in=min(2d3,max(4d0,q2))
          sd=log(log(q2in/alam**2)/log(4d0/alam**2))
        endif
 
C...Calculate parton distributions.
        do 120 kfl=1,4
          do 110 is=1,5
            ts(is)=cow(1,is,kfl,nset)+cow(2,is,kfl,nset)*sd+
     &      cow(3,is,kfl,nset)*sd**2
  110     continue
          if(kfl.eq.1) then
            xq(kfl)=x**ts(1)*(1d0-x)**ts(2)/eulbet(ts(1),ts(2)+1d0)
          else
            xq(kfl)=ts(1)*x**ts(2)*(1d0-x)**ts(3)*(1d0+ts(4)*x+
     &      ts(5)*x**2)
          endif
  120   continue
 
C...Put into output array.
        xppi(0)=xq(2)
        xppi(1)=xq(3)/6d0
        xppi(2)=xq(1)+xq(3)/6d0
        xppi(3)=xq(3)/6d0
        xppi(4)=xq(4)
        xppi(-1)=xq(1)+xq(3)/6d0
        xppi(-2)=xq(3)/6d0
        xppi(-3)=xq(3)/6d0
        xppi(-4)=xq(4)
 
C...Leading order pion parton distributions from Gluck, Reya and Vogt.
C...Allowed variable range: 0.25 GeV^2 < Q^2 < 10^8 GeV^2 and
C...10^-5 < x < 1.
      else
 
C...Determine s expansion variable and some x expressions.
        vint(231)=0.25d0
        if(mstp(57).le.0) then
          sd=0d0
        else
          q2in=min(1d8,max(0.25d0,q2))
          sd=log(log(q2in/0.232d0**2)/log(0.25d0/0.232d0**2))
        endif
        sd2=sd**2
        xl=-log(x)
        xs=sqrt(x)
 
C...Evaluate valence, gluon and sea distributions.
        xfval=(0.519d0+0.180d0*sd-0.011d0*sd2)*x**(0.499d0-0.027d0*sd)*
     &  (1d0+(0.381d0-0.419d0*sd)*xs)*(1d0-x)**(0.367d0+0.563d0*sd)
        xfglu=(x**(0.482d0+0.341d0*sqrt(sd))*((0.678d0+0.877d0*
     &  sd-0.175d0*sd2)+
     &  (0.338d0-1.597d0*sd)*xs+(-0.233d0*sd+0.406d0*sd2)*x)+
     &  sd**0.599d0*exp(-(0.618d0+2.070d0*sd)+sqrt(3.676d0*sd**1.263d0*
     &  xl)))*
     &  (1d0-x)**(0.390d0+1.053d0*sd)
        xfsea=sd**0.55d0*(1d0-0.748d0*xs+(0.313d0+0.935d0*sd)*x)*(1d0-
     &  x)**3.359d0*
     &  exp(-(4.433d0+1.301d0*sd)+sqrt((9.30d0-0.887d0*sd)*sd**0.56d0*
     &  xl))/
     &  xl**(2.538d0-0.763d0*sd)
        if(sd.le.0.888d0) then
          xfchm=0d0
        else
          xfchm=(sd-0.888d0)**1.02d0*(1d0+1.008d0*x)*(1d0-x)**(1.208d0+
     &    0.771d0*sd)*
     &    exp(-(4.40d0+1.493d0*sd)+sqrt((2.032d0+1.901d0*sd)*sd**0.39d0*
     &    xl))
        endif
        if(sd.le.1.351d0) then
          xfbot=0d0
        else
          xfbot=(sd-1.351d0)**1.03d0*(1d0-x)**(0.697d0+0.855d0*sd)*
     &    exp(-(4.51d0+1.490d0*sd)+sqrt((3.056d0+1.694d0*sd)*sd**0.39d0*
     &    xl))
        endif
 
C...Put into output array.
        xppi(0)=xfglu
        xppi(1)=xfsea
        xppi(2)=xfsea
        xppi(3)=xfsea
        xppi(4)=xfchm
        xppi(5)=xfbot
        do 130 kfl=1,5
          xppi(-kfl)=xppi(kfl)
  130   continue
        xppi(2)=xppi(2)+xfval
        xppi(-1)=xppi(-1)+xfval
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYPDPR
C...Gives proton parton distributions according to a few different
C...parametrizations.
c...mstp(51)=1 EHLQ set 1 (1986 updated version).
c...        =2 EHLQ set 2 (1986 updated version).
c...        =3 Duke-Owens set 1.
c...        =4 Duke-Owens set 2.
c...        =5 CTEQ2M (best MSbar fit)
c...        =6 CTEQ2MS (singuar at small x)
c...        =7 CTEQ2MF (flat at small x)
c...        =8 CTEQ2ML (large $\Lambda$)
c...        =9 CTEQ2L (best leading order fit)
c...        =10 CTEQ2D (best DIS fit)
c...        =11 CRV LO (1992 updated version)
c...        =12 CTEQ 3L (leading order). 
c...        =13 CTEQ 3M (MSbar).
c...        =14 CTEQ 3D (DIS).
c...        =15 GRV 94L (leading order).  % *
c...        =16 GRV 94M (MSbar).
c...        =17 GRV 94D (DIS).
c...References:
c...CTEQ Collaboration, H.L. Lai et al., Phys. Rev. D51 (1995) 4763; 
c...M. Gluck, E. Reya and A. Vogt, Z. Phys. C67 (1995) 433 

      subroutine pjpdpr(x,q2,xppr)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      save /jydat1/,/jydat2/,/pjpars/,/pjint1/
C...Arrays and data.
      dimension xppr(-6:6),q2min(6)
      data q2min/ 2.56d0, 2.56d0, 2.56d0, 0.4d0, 0.4d0, 0.4d0/
 
C...Reset output array.
      do 100 kfl=-6,6
        xppr(kfl)=0d0
  100 continue
 
cjam
      if(mstp(51).ge.12) then

C...Common preliminaries.
      nset=max(1,min(6,mstp(51)-11))
      vint(231)=q2min(nset)
      if(mstp(57).eq.0) then
        q2l=q2min(nset)
      else
        q2l=max(q2min(nset),q2)
      endif
 
      if(nset.ge.1.and.nset.le.3) then
C...Interface to the CTEQ 3 parton distributions.
        qrt=sqrt(max(1d0,q2l))
 
C...Loop over flavours.
        do 110 i=-6,6
          if(i.le.0) then
            xppr(i)=pjcteq(nset,i,x,qrt)
          elseif(i.le.2) then
            xppr(i)=pjcteq(nset,i,x,qrt)+xppr(-i)
          else
            xppr(i)=xppr(-i)
          endif
  110   continue
 
      elseif(nset.ge.4.and.nset.le.6) then
C...Interface to the GRV 94 distributions.
        if(nset.eq.4) then
          call pjgrvl (x, q2l, uv, dv, del, udb, sb, chm, bot, gl)
        elseif(nset.eq.5) then
          call pjgrvm (x, q2l, uv, dv, del, udb, sb, chm, bot, gl)
        else
          call pjgrvd (x, q2l, uv, dv, del, udb, sb, chm, bot, gl)
        endif
 
C...Put into output array.
        xppr(0)=gl
        xppr(-1)=0.5d0*(udb+del)
        xppr(-2)=0.5d0*(udb-del)
        xppr(-3)=sb
        xppr(-4)=chm
        xppr(-5)=bot
        xppr(1)=dv+xppr(-1)
        xppr(2)=uv+xppr(-2)
        xppr(3)=sb
        xppr(4)=chm
        xppr(5)=bot
 
      endif
cjam++
      else if(mstp(51).eq.1.or.mstp(51).eq.2) then
C...Eichten, Hinchliffe, Lane, Quigg.
        call ehlq(x,q2,xppr)
      elseif(mstp(51).eq.3.or.mstp(51).eq.4) then
C...Duke, Owens.
      call dukeowen(x,q2,xppr)
      elseif(mstp(51).ge.5.and.mstp(51).le.10) then
C...Interface to the CTEQ 2 structure functions.
        nset=mstp(51)-4
        qrt=sqrt(max(1.d0,q2))
 
C...Loop over flavours; put u and d in right order.
        do 200 i=-6,2
        kfl=i
        if(i.eq.1) kfl=2
        if(i.eq.2) kfl=1
        if(i.eq.-1) kfl=-2
        if(i.eq.-2) kfl=-1
        if(i.le.0) then
          xppr(kfl)=pjctq2(nset,i,x,qrt)
          xppr(-kfl)=xppr(kfl)
        else
          xppr(kfl)=pjctq2(nset,i,x,qrt)+xppr(-kfl)
        endif
  200   continue
 
C...Leading order proton structure functions from Gluck, Reya and Vogt.
C...Allowed variable range: 0.25 GeV^2 < Q^2 < 10^8 GeV^2 and
C...10^-5 < x < 1.
      else if(mstp(51).eq.11) then
C...Determine s expansion variable and some x expressions.
        vint(231)=0.25d0
        if(mstp(57).le.0) then
          sd=0.d0
        else
          q2in=min(1d8,max(0.25d0,q2))
          sd=log(log(q2in/0.232d0**2)/log(0.25d0/0.232d0**2))
        endif
        sd2=sd**2
        xl=-log(x)
        xs=sqrt(x)
 
C...Evaluate valence, gluon and sea distributions.
        xfvud=(0.663d0+0.191d0*sd-0.041d0*sd2+0.031d0*sd**3)*x**0.326d0*
     &  (1.d0+(-1.97d0+6.74d0*sd-1.96d0*sd2)*xs
     $                 +(24.4d0-20.7d0*sd+4.08d0*sd2)*x)*
     &  (1.d0-x)**(2.86d0+0.70d0*sd-0.02d0*sd2)
        xfvdd=(0.579d0+0.283d0*sd+0.047d0*sd2)*x**(0.523d0-0.015d0*sd)*
     &  (1d0+(2.22d0-0.59d0*sd-0.27d0*sd2)*xs
     $                    +(5.95d0-6.19d0*sd+1.55d0*sd2)*x)*
     &  (1d0-x)**(3.57d0+0.94d0*sd-0.16d0*sd2)
        xfglu=(x**(1d0-0.17d0*sd)*((4.879d0*sd-1.383d0*sd2)+
     &  (25.92d0-28.97d0*sd+5.596d0*sd2)*x
     $                +(-25.69d0+23.68d0*sd-1.975d0*sd2)*x**2)+
     &  sd**0.558d0*exp(-(0.595d0+2.138d0*sd)
     $                         +sqrt(4.066d0*sd**1.218d0*xl)))*
     &  (1.d0-x)**(2.537d0+1.718d0*sd+0.353d0*sd2)
        xfsea=(x**(0.412d0-0.171d0*sd)*(0.363d0-1.196d0*x+
     &  (1.029d0+1.785d0*sd-0.459d0*sd2)*x**2)*xl**(0.566d0-0.496d0*sd)+
     &  sd**1.396d0*exp(-(3.838d0+1.944d0*sd)
     $                              +sqrt(2.845d0*sd**1.331d0*xl)))*
     &  (1.d0-x)**(4.696d0+2.109d0*sd)
        xfstr=sd**0.803d0*(1.d0+(-3.055d0+1.024d0*sd**0.67d0)*xs+
     &  (27.4d0-20.0d0*sd**0.154d0)*x)*(1.d0-x)**6.22d0*
     &  exp(-(4.33d0+1.408d0*sd)
     $      +sqrt((8.27d0-0.437d0*sd)*sd**0.563d0*xl))
     $    /xl**(2.082d0-0.577d0*sd)
        if(sd.le.0.888d0) then
          xfchm=0.d0
        else
          xfchm=(sd-0.888d0)**1.01d0*(1.d0+(4.24d0-0.804d0*sd)*x)*
     &    (1.d0-x)**(3.46d0+1.076d0*sd)*exp(-(4.61d0+1.49d0*sd)+
     &    sqrt((2.555d0+1.961d0*sd)*sd**0.37d0*xl))
        endif
        if(sd.le.1.351d0) then
          xfbot=0.d0
        else
          xfbot=(sd-1.351d0)*(1.d0+1.848d0*x)
     $         *(1.d0-x)**(2.929d0+1.396d0*sd)*
     &         exp(-(4.71d0+1.514d0*sd)
     $         +sqrt((4.02d0+1.239d0*sd)*sd**0.51d0*xl))
        endif
 
C...Put into output array.
        xppr(0)=xfglu
        xppr(1)=xfvdd+xfsea
        xppr(2)=xfvud-xfvdd+xfsea
        xppr(3)=xfstr
        xppr(4)=xfchm
        xppr(5)=xfbot
        xppr(-1)=xfsea
        xppr(-2)=xfsea
        xppr(-3)=xfstr
        xppr(-4)=xfchm
        xppr(-5)=xfbot
 
      else
        call pjerrm(40,'(pjpdpr:) invalid mstp(51)')
      endif
cjam--
 
      return
      end
 
C*********************************************************************
 
C...PYCTEQ
C...Gives the CTEQ 3 parton distribution function sets in
C...parametrized form, of October 24, 1994.
C...Authors: H.L. Lai, J. Botts, J. Huston, J.G. Morfin, J.F. Owens,
C...J. Qiu, W.K. Tung and H. Weerts.
 
      function pjcteq (iset, iprt, x, q)
 
C...Double precision declaration.
      implicit double precision(a-h, o-z)
 
C...Data on Lambda values of fits, minimum Q and quark masses.
      dimension alm(3), qms(4:6)
      data alm / 0.177d0, 0.239d0, 0.247d0 /
      data qmn / 1.60d0 /, (qms(i), i=4,6) / 1.60d0, 5.00d0, 180.0d0 /
 
C....Check flavour thresholds. Set up QI for SB.
      ip = iabs(iprt)
      if(ip .ge. 4) then
        if(q .le. qms(ip)) then
          pjcteq = 0d0
          return
        endif
        qi = qms(ip)
      else
        qi = qmn
      endif
 
C...Use "standard lambda" of parametrization program for expansion.
      alam = alm (iset)
      sbl = log(q/alam) / log(qi/alam)
      sb = log (sbl)
      sb2 = sb*sb
      sb3 = sb2*sb
 
C...Expansion for CTEQ3L.
      if(iset .eq. 1) then
        if(iprt .eq. 2) then
          a0=exp( 0.1907d+00+0.4205d-01*sb +0.2752d+00*sb2-
     &    0.3171d+00*sb3)
          a1= 0.4611d+00+0.2331d-01*sb -0.3403d-01*sb2+0.3174d-01*sb3
          a2= 0.3504d+01+0.5739d+00*sb +0.2676d+00*sb2-0.1553d+00*sb3
          a3= 0.7452d+01-0.6742d+01*sb +0.2849d+01*sb2-0.1964d+00*sb3
          a4= 0.1116d+01-0.3435d+00*sb +0.2865d+00*sb2-0.1288d+00*sb3
          a5= 0.6659d-01+0.2714d+00*sb -0.2688d+00*sb2+0.2763d+00*sb3
        elseif(iprt .eq. 1) then
          a0=exp( 0.1141d+00+0.4764d+00*sb -0.1745d+01*sb2+
     &    0.7728d+00*sb3)
          a1= 0.4275d+00-0.1290d+00*sb +0.3609d+00*sb2-0.1689d+00*sb3
          a2= 0.3000d+01+0.2946d+01*sb -0.4117d+01*sb2+0.1989d+01*sb3
          a3=-0.1302d+01+0.2322d+01*sb -0.4258d+01*sb2+0.2109d+01*sb3
          a4= 0.2586d+01-0.1920d+00*sb -0.3754d+00*sb2+0.2731d+00*sb3
          a5=-0.2251d+00-0.5374d+00*sb +0.2245d+01*sb2-0.1034d+01*sb3
        elseif(iprt .eq. 0) then
          a0=exp(-0.7631d+00-0.7241d+00*sb -0.1170d+01*sb2+
     &    0.5343d+00*sb3)
          a1=-0.3573d+00+0.3469d+00*sb -0.3396d+00*sb2+0.9188d-01*sb3
          a2= 0.5604d+01+0.7458d+00*sb -0.5082d+00*sb2+0.1844d+00*sb3
          a3= 0.1549d+02-0.1809d+02*sb +0.1162d+02*sb2-0.3483d+01*sb3
          a4= 0.9881d+00+0.1364d+00*sb -0.4421d+00*sb2+0.2051d+00*sb3
          a5=-0.9505d-01+0.3259d+01*sb -0.1547d+01*sb2+0.2918d+00*sb3
        elseif(iprt .eq. -1) then
          a0=exp(-0.2449d+01-0.3513d+01*sb +0.4529d+01*sb2-
     &    0.2031d+01*sb3)
          a1=-0.4050d+00+0.3411d+00*sb -0.3669d+00*sb2+0.1109d+00*sb3
          a2= 0.7470d+01-0.2982d+01*sb +0.5503d+01*sb2-0.2419d+01*sb3
          a3= 0.1503d+02+0.1638d+01*sb -0.8772d+01*sb2+0.3852d+01*sb3
          a4= 0.1137d+01-0.1006d+01*sb +0.1485d+01*sb2-0.6389d+00*sb3
          a5=-0.5299d+00+0.3160d+01*sb -0.3104d+01*sb2+0.1219d+01*sb3
        elseif(iprt .eq. -2) then
          a0=exp(-0.2740d+01-0.7987d-01*sb -0.9015d+00*sb2-
     &    0.9872d-01*sb3)
          a1=-0.3909d+00+0.1244d+00*sb -0.4487d-01*sb2+0.1277d-01*sb3
          a2= 0.9163d+01+0.2823d+00*sb -0.7720d+00*sb2-0.9360d-02*sb3
          a3= 0.1080d+02-0.3915d+01*sb -0.1153d+01*sb2+0.2649d+01*sb3
          a4= 0.9894d+00-0.1647d+00*sb -0.9426d-02*sb2+0.2945d-02*sb3
          a5=-0.3395d+00+0.6998d+00*sb +0.7000d+00*sb2-0.6730d-01*sb3
        elseif(iprt .eq. -3) then
          a0=exp(-0.3640d+01+0.1250d+01*sb -0.2914d+01*sb2+
     &    0.8390d+00*sb3)
          a1=-0.3595d+00-0.5259d-01*sb +0.3122d+00*sb2-0.1642d+00*sb3
          a2= 0.7305d+01+0.9727d+00*sb -0.9788d+00*sb2-0.5193d-01*sb3
          a3= 0.1198d+02-0.1799d+02*sb +0.2614d+02*sb2-0.1091d+02*sb3
          a4= 0.9882d+00-0.6101d+00*sb +0.9737d+00*sb2-0.4935d+00*sb3
          a5=-0.1186d+00-0.3231d+00*sb +0.3074d+01*sb2-0.1274d+01*sb3
        elseif(iprt .eq. -4) then
          a0=sb** 0.1122d+01*exp(-0.3718d+01-0.1335d+01*sb +
     &    0.1651d-01*sb2)
          a1=-0.4719d+00+0.7509d+00*sb -0.8420d+00*sb2+0.2901d+00*sb3
          a2= 0.6194d+01-0.1641d+01*sb +0.4907d+01*sb2-0.2523d+01*sb3
          a3= 0.4426d+01-0.4270d+01*sb +0.6581d+01*sb2-0.3474d+01*sb3
          a4= 0.2683d+00+0.9876d+00*sb -0.7612d+00*sb2+0.1780d+00*sb3
          a5=-0.4547d+00+0.4410d+01*sb -0.3712d+01*sb2+0.1245d+01*sb3
        elseif(iprt .eq. -5) then
          a0=sb** 0.9838d+00*exp(-0.2548d+01-0.7660d+01*sb +
     &    0.3702d+01*sb2)
          a1=-0.3122d+00-0.2120d+00*sb +0.5716d+00*sb2-0.3773d+00*sb3
          a2= 0.6257d+01-0.8214d-01*sb -0.2537d+01*sb2+0.2981d+01*sb3
          a3=-0.6723d+00+0.2131d+01*sb +0.9599d+01*sb2-0.7910d+01*sb3
          a4= 0.9169d-01+0.4295d-01*sb -0.5017d+00*sb2+0.3811d+00*sb3
          a5= 0.2402d+00+0.2656d+01*sb -0.1586d+01*sb2+0.2880d+00*sb3
        elseif(iprt .eq. -6) then
          a0=sb** 0.1001d+01*exp(-0.6934d+01+0.3050d+01*sb -
     &    0.6943d+00*sb2)
          a1=-0.1713d+00-0.5167d+00*sb +0.1241d+01*sb2-0.1703d+01*sb3
          a2= 0.6169d+01+0.3023d+01*sb -0.1972d+02*sb2+0.1069d+02*sb3
          a3= 0.4439d+01-0.1746d+02*sb +0.1225d+02*sb2+0.8350d+00*sb3
          a4= 0.5458d+00-0.4586d+00*sb +0.9089d+00*sb2-0.4049d+00*sb3
          a5= 0.3207d+01-0.3362d+01*sb +0.5877d+01*sb2-0.7659d+01*sb3
        endif
 
C...Expansion for CTEQ3M.
      elseif(iset .eq. 2) then
        if(iprt .eq. 2) then
          a0=exp( 0.2259d+00+0.1237d+00*sb +0.3035d+00*sb2-
     &    0.2935d+00*sb3)
          a1= 0.5085d+00+0.1651d-01*sb -0.3592d-01*sb2+0.2782d-01*sb3
          a2= 0.3732d+01+0.4901d+00*sb +0.2218d+00*sb2-0.1116d+00*sb3
          a3= 0.7011d+01-0.6620d+01*sb +0.2557d+01*sb2-0.1360d+00*sb3
          a4= 0.8969d+00-0.2429d+00*sb +0.1811d+00*sb2-0.6888d-01*sb3
          a5= 0.8636d-01+0.2558d+00*sb -0.3082d+00*sb2+0.2535d+00*sb3
        elseif(iprt .eq. 1) then
          a0=exp(-0.7266d+00-0.1584d+01*sb +0.1259d+01*sb2-
     &    0.4305d-01*sb3)
          a1= 0.5285d+00-0.3721d+00*sb +0.5150d+00*sb2-0.1697d+00*sb3
          a2= 0.4075d+01+0.8282d+00*sb -0.4496d+00*sb2+0.2107d+00*sb3
          a3= 0.3279d+01+0.5066d+01*sb -0.9134d+01*sb2+0.2897d+01*sb3
          a4= 0.4399d+00-0.5888d+00*sb +0.4802d+00*sb2-0.1664d+00*sb3
          a5= 0.3678d+00-0.8929d+00*sb +0.1592d+01*sb2-0.5713d+00*sb3
        elseif(iprt .eq. 0) then
          a0=exp(-0.2318d+00-0.9779d+00*sb -0.3783d+00*sb2+
     &    0.1037d-01*sb3)
          a1=-0.2916d+00+0.1754d+00*sb -0.1884d+00*sb2+0.6116d-01*sb3
          a2= 0.5349d+01+0.7460d+00*sb +0.2319d+00*sb2-0.2622d+00*sb3
          a3= 0.6920d+01-0.3454d+01*sb +0.2027d+01*sb2-0.7626d+00*sb3
          a4= 0.1013d+01+0.1423d+00*sb -0.1798d+00*sb2+0.1872d-01*sb3
          a5=-0.5465d-01+0.2303d+01*sb -0.9584d+00*sb2+0.3098d+00*sb3
        elseif(iprt .eq. -1) then
          a0=exp(-0.2328d+01-0.3061d+01*sb +0.3620d+01*sb2-
     &    0.1602d+01*sb3)
          a1=-0.3358d+00+0.3198d+00*sb -0.4210d+00*sb2+0.1571d+00*sb3
          a2= 0.8478d+01-0.3112d+01*sb +0.5243d+01*sb2-0.2255d+01*sb3
          a3= 0.1971d+02+0.3389d+00*sb -0.5268d+01*sb2+0.2099d+01*sb3
          a4= 0.1128d+01-0.4701d+00*sb +0.7779d+00*sb2-0.3506d+00*sb3
          a5=-0.4708d+00+0.3341d+01*sb -0.3375d+01*sb2+0.1353d+01*sb3
        elseif(iprt .eq. -2) then
          a0=exp(-0.2906d+01-0.1069d+00*sb -0.1055d+01*sb2+
     &    0.2496d+00*sb3)
          a1=-0.2875d+00+0.6571d-01*sb -0.1987d-01*sb2-0.1800d-02*sb3
          a2= 0.9854d+01-0.2715d+00*sb -0.7407d+00*sb2+0.2888d+00*sb3
          a3= 0.1583d+02-0.7687d+01*sb +0.3428d+01*sb2-0.3327d+00*sb3
          a4= 0.9763d+00+0.7599d-01*sb -0.2128d+00*sb2+0.6852d-01*sb3
          a5=-0.8444d-02+0.9434d+00*sb +0.4152d+00*sb2-0.1481d+00*sb3
        elseif(iprt .eq. -3) then
          a0=exp(-0.3780d+01+0.2499d+01*sb -0.4962d+01*sb2+
     &    0.1936d+01*sb3)
          a1=-0.2639d+00-0.1575d+00*sb +0.3584d+00*sb2-0.1646d+00*sb3
          a2= 0.8082d+01+0.2794d+01*sb -0.5438d+01*sb2+0.2321d+01*sb3
          a3= 0.1811d+02-0.2000d+02*sb +0.1951d+02*sb2-0.6904d+01*sb3
          a4= 0.9822d+00+0.4972d+00*sb -0.8690d+00*sb2+0.3415d+00*sb3
          a5= 0.1772d+00-0.6078d+00*sb +0.3341d+01*sb2-0.1473d+01*sb3
        elseif(iprt .eq. -4) then
          a0=sb** 0.1122d+01*exp(-0.4232d+01-0.1808d+01*sb +
     &    0.5348d+00*sb2)
          a1=-0.2824d+00+0.5846d+00*sb -0.7230d+00*sb2+0.2419d+00*sb3
          a2= 0.5683d+01-0.2948d+01*sb +0.5916d+01*sb2-0.2560d+01*sb3
          a3= 0.2051d+01+0.4795d+01*sb -0.4271d+01*sb2+0.4174d+00*sb3
          a4= 0.1737d+00+0.1717d+01*sb -0.1978d+01*sb2+0.6643d+00*sb3
          a5= 0.8689d+00+0.3500d+01*sb -0.3283d+01*sb2+0.1026d+01*sb3
        elseif(iprt .eq. -5) then
          a0=sb** 0.9906d+00*exp(-0.1496d+01-0.6576d+01*sb +
     &    0.1569d+01*sb2)
          a1=-0.2140d+00-0.6419d-01*sb -0.2741d-02*sb2+0.3185d-02*sb3
          a2= 0.5781d+01+0.1049d+00*sb -0.3930d+00*sb2+0.5174d+00*sb3
          a3=-0.9420d+00+0.5511d+00*sb +0.8817d+00*sb2+0.1903d+01*sb3
          a4= 0.2418d-01+0.4232d-01*sb -0.1244d-01*sb2-0.2365d-01*sb3
          a5= 0.7664d+00+0.1794d+01*sb -0.4917d+00*sb2-0.1284d+00*sb3
        elseif(iprt .eq. -6) then
          a0=sb** 0.1000d+01*exp(-0.8460d+01+0.1154d+01*sb +
     &    0.8838d+01*sb2)
          a1=-0.4316d-01-0.2976d+00*sb +0.3174d+00*sb2-0.1429d+01*sb3
          a2= 0.4910d+01+0.2273d+01*sb +0.5631d+01*sb2-0.1994d+02*sb3
          a3= 0.1190d+02-0.2000d+02*sb -0.2000d+02*sb2+0.1292d+02*sb3
          a4= 0.5771d+00-0.2552d+00*sb +0.7510d+00*sb2+0.6923d+00*sb3
          a5= 0.4402d+01-0.1627d+01*sb -0.2085d+01*sb2-0.6737d+01*sb3
        endif
 
C...Expansion for CTEQ3D.
      elseif(iset .eq. 3) then
        if(iprt .eq. 2) then
          a0=exp( 0.2148d+00+0.5814d-01*sb +0.2734d+00*sb2-
     &    0.2902d+00*sb3)
          a1= 0.4810d+00+0.1657d-01*sb -0.3800d-01*sb2+0.3125d-01*sb3
          a2= 0.3509d+01+0.3923d+00*sb +0.4010d+00*sb2-0.1932d+00*sb3
          a3= 0.7055d+01-0.6552d+01*sb +0.3466d+01*sb2-0.5657d+00*sb3
          a4= 0.1061d+01-0.3453d+00*sb +0.4089d+00*sb2-0.1817d+00*sb3
          a5= 0.8687d-01+0.2548d+00*sb -0.2967d+00*sb2+0.2647d+00*sb3
        elseif(iprt .eq. 1) then
          a0=exp( 0.3961d+00+0.4914d+00*sb -0.1728d+01*sb2+
     &    0.7257d+00*sb3)
          a1= 0.4162d+00-0.1419d+00*sb +0.3680d+00*sb2-0.1618d+00*sb3
          a2= 0.3248d+01+0.3028d+01*sb -0.4307d+01*sb2+0.1920d+01*sb3
          a3=-0.1100d+01+0.2184d+01*sb -0.3820d+01*sb2+0.1717d+01*sb3
          a4= 0.2082d+01-0.2756d+00*sb +0.3043d+00*sb2-0.1260d+00*sb3
          a5=-0.4822d+00-0.5706d+00*sb +0.2243d+01*sb2-0.9760d+00*sb3
        elseif(iprt .eq. 0) then
          a0=exp(-0.4665d+00-0.7554d+00*sb -0.3323d+00*sb2-
     &    0.2734d-04*sb3)
          a1=-0.3359d+00+0.2395d+00*sb -0.2377d+00*sb2+0.7059d-01*sb3
          a2= 0.5451d+01+0.6086d+00*sb +0.8606d-01*sb2-0.1425d+00*sb3
          a3= 0.1026d+02-0.9352d+01*sb +0.4879d+01*sb2-0.1150d+01*sb3
          a4= 0.9935d+00-0.5017d-01*sb -0.1707d-01*sb2-0.1464d-02*sb3
          a5=-0.4160d-01+0.2305d+01*sb -0.1063d+01*sb2+0.3211d+00*sb3
        elseif(iprt .eq. -1) then
          a0=exp(-0.2714d+01-0.2868d+01*sb +0.3700d+01*sb2-
     &    0.1671d+01*sb3)
          a1=-0.3893d+00+0.3341d+00*sb -0.3897d+00*sb2+0.1420d+00*sb3
          a2= 0.8359d+01-0.3267d+01*sb +0.5327d+01*sb2-0.2245d+01*sb3
          a3= 0.2359d+02-0.5669d+01*sb -0.4602d+01*sb2+0.3153d+01*sb3
          a4= 0.1106d+01-0.4745d+00*sb +0.7739d+00*sb2-0.3417d+00*sb3
          a5=-0.5557d+00+0.3433d+01*sb -0.3390d+01*sb2+0.1354d+01*sb3
        elseif(iprt .eq. -2) then
          a0=exp(-0.3323d+01+0.2296d+00*sb -0.1109d+01*sb2+
     &    0.2223d+00*sb3)
          a1=-0.3410d+00+0.8847d-01*sb -0.1111d-01*sb2-0.5927d-02*sb3
          a2= 0.9753d+01-0.5182d+00*sb -0.4670d+00*sb2+0.1921d+00*sb3
          a3= 0.1977d+02-0.1600d+02*sb +0.9481d+01*sb2-0.1864d+01*sb3
          a4= 0.9818d+00+0.2839d-02*sb -0.1188d+00*sb2+0.3584d-01*sb3
          a5=-0.7934d-01+0.1004d+01*sb +0.3704d+00*sb2-0.1220d+00*sb3
        elseif(iprt .eq. -3) then
          a0=exp(-0.3985d+01+0.2855d+01*sb -0.5208d+01*sb2+
     &    0.1937d+01*sb3)
          a1=-0.3337d+00-0.1150d+00*sb +0.3691d+00*sb2-0.1709d+00*sb3
          a2= 0.7968d+01+0.3641d+01*sb -0.6599d+01*sb2+0.2642d+01*sb3
          a3= 0.1873d+02-0.1999d+02*sb +0.1734d+02*sb2-0.5813d+01*sb3
          a4= 0.9731d+00+0.5082d+00*sb -0.8780d+00*sb2+0.3231d+00*sb3
          a5=-0.5542d-01-0.4189d+00*sb +0.3309d+01*sb2-0.1439d+01*sb3
        elseif(iprt .eq. -4) then
          a0=sb** 0.1105d+01*exp(-0.3952d+01-0.1901d+01*sb +
     &    0.5137d+00*sb2)
          a1=-0.3543d+00+0.6055d+00*sb -0.6941d+00*sb2+0.2278d+00*sb3
          a2= 0.5955d+01-0.2629d+01*sb +0.5337d+01*sb2-0.2300d+01*sb3
          a3= 0.1933d+01+0.4882d+01*sb -0.3810d+01*sb2+0.2290d+00*sb3
          a4= 0.1806d+00+0.1655d+01*sb -0.1893d+01*sb2+0.6395d+00*sb3
          a5= 0.4790d+00+0.3612d+01*sb -0.3152d+01*sb2+0.9684d+00*sb3
        elseif(iprt .eq. -5) then
          a0=sb** 0.9818d+00*exp(-0.1825d+01-0.7464d+01*sb +
     &    0.2143d+01*sb2)
          a1=-0.2604d+00-0.1400d+00*sb +0.1702d+00*sb2-0.8476d-01*sb3
          a2= 0.6005d+01+0.6275d+00*sb -0.2535d+01*sb2+0.2219d+01*sb3
          a3=-0.9067d+00+0.1149d+01*sb +0.1974d+01*sb2+0.4716d+01*sb3
          a4= 0.3915d-01+0.5945d-01*sb -0.9844d-01*sb2+0.2783d-01*sb3
          a5= 0.5500d+00+0.1994d+01*sb -0.6727d+00*sb2-0.1510d+00*sb3
        elseif(iprt .eq. -6) then
          a0=sb** 0.1002d+01*exp(-0.8553d+01+0.3793d+00*sb +
     &    0.9998d+01*sb2)
          a1=-0.5870d-01-0.2792d+00*sb +0.6526d+00*sb2-0.1984d+01*sb3
          a2= 0.4716d+01+0.4473d+00*sb +0.1128d+02*sb2-0.1937d+02*sb3
          a3= 0.1289d+02-0.1742d+02*sb -0.1983d+02*sb2-0.9274d+00*sb3
          a4= 0.5647d+00-0.2732d+00*sb +0.1074d+01*sb2+0.5981d+00*sb3
          a5= 0.4390d+01-0.1262d+01*sb -0.9026d+00*sb2-0.9394d+01*sb3
        endif
      endif
 
C...Calculation of x * f(x, Q).
      pjcteq = max(0d0, a0 *(x**a1) *((1d0-x)**a2) *(1d0+a3*(x**a4))
     &   *(log(1d0+1d0/x))**a5 )
 
      return
      end
 
C*********************************************************************
 
C...PYGRVL
C...Gives the GRV 94 L (leading order) parton distribution function set
C...in parametrized form.
C...Authors: M. Glueck, E. Reya and A. Vogt.
 
      subroutine pjgrvl (x, q2, uv, dv, del, udb, sb, chm, bot, gl)
 
C...Double precision declaration.
      implicit double precision (a - z)
 
C...Common expressions.
      mu2  = 0.23d0
      lam2 = 0.2322d0 * 0.2322d0
      s  = log (log(q2/lam2) / log(mu2/lam2))
      ds = sqrt (s)
      s2 = s * s
      s3 = s2 * s
 
C...uv :
      nu  =  2.284d0 + 0.802d0 * s + 0.055d0 * s2
      aku =  0.590d0 - 0.024d0 * s
      bku =  0.131d0 + 0.063d0 * s
      au  = -0.449d0 - 0.138d0 * s - 0.076d0 * s2
      bu  =  0.213d0 + 2.669d0 * s - 0.728d0 * s2
      cu  =  8.854d0 - 9.135d0 * s + 1.979d0 * s2
      du  =  2.997d0 + 0.753d0 * s - 0.076d0 * s2
      uv  = pjgrvv (x, nu, aku, bku, au, bu, cu, du)
 
C...dv :
      nd  =  0.371d0 + 0.083d0 * s + 0.039d0 * s2
      akd =  0.376d0
      bkd =  0.486d0 + 0.062d0 * s
      ad  = -0.509d0 + 3.310d0 * s - 1.248d0 * s2
      bd  =  12.41d0 - 10.52d0 * s + 2.267d0 * s2
      cd  =  6.373d0 - 6.208d0 * s + 1.418d0 * s2
      dd  =  3.691d0 + 0.799d0 * s - 0.071d0 * s2
      dv  = pjgrvv (x, nd, akd, bkd, ad, bd, cd, dd)
 
C...del :
      ne  =  0.082d0 + 0.014d0 * s + 0.008d0 * s2
      ake =  0.409d0 - 0.005d0 * s
      bke =  0.799d0 + 0.071d0 * s
      ae  = -38.07d0 + 36.13d0 * s - 0.656d0 * s2
      be  =  90.31d0 - 74.15d0 * s + 7.645d0 * s2
      ce  =  0.0d0
      de  =  7.486d0 + 1.217d0 * s - 0.159d0 * s2
      del = pjgrvv (x, ne, ake, bke, ae, be, ce, de)
 
C...udb :
      alx =  1.451d0
      bex =  0.271d0
      akx =  0.410d0 - 0.232d0 * s
      bkx =  0.534d0 - 0.457d0 * s
      agx =  0.890d0 - 0.140d0 * s
      bgx = -0.981d0
      cx  =  0.320d0 + 0.683d0 * s
      dx  =  4.752d0 + 1.164d0 * s + 0.286d0 * s2
      ex  =  4.119d0 + 1.713d0 * s
      esx =  0.682d0 + 2.978d0 * s
      udb = pjgrvw (x, s, alx, bex, akx, bkx, agx, bgx, cx,
     & dx, ex, esx)
 
C...sb :
      sts =  0d0
      als =  0.914d0
      bes =  0.577d0
      aks =  1.798d0 - 0.596d0 * s
      as  = -5.548d0 + 3.669d0 * ds - 0.616d0 * s
      bs  =  18.92d0 - 16.73d0 * ds + 5.168d0 * s
      dst =  6.379d0 - 0.350d0 * s  + 0.142d0 * s2
      est =  3.981d0 + 1.638d0 * s
      ess =  6.402d0
      sb  = pjgrvs (x, s, sts, als, bes, aks, as, bs, dst, est, ess)
 
C...cb :
      stc =  0.888d0
      alc =  1.01d0
      bec =  0.37d0
      akc =  0d0
      ac  =  0d0
      bc  =  4.24d0  - 0.804d0 * s
      dct =  3.46d0  - 1.076d0 * s
      ect =  4.61d0  + 1.49d0  * s
      esc =  2.555d0 + 1.961d0 * s
      chm = pjgrvs (x, s, stc, alc, bec, akc, ac, bc, dct, ect, esc)
 
C...bb :
      stb =  1.351d0
      alb =  1.00d0
      beb =  0.51d0
      akb =  0d0
      ab  =  0d0
      bb  =  1.848d0
      dbt =  2.929d0 + 1.396d0 * s
      ebt =  4.71d0  + 1.514d0 * s
      esb =  4.02d0  + 1.239d0 * s
      bot = pjgrvs (x, s, stb, alb, beb, akb, ab, bb, dbt, ebt, esb)
 
C...gl :
      alg =  0.524d0
      beg =  1.088d0
      akg =  1.742d0 - 0.930d0 * s
      bkg =                         - 0.399d0 * s2
      ag  =  7.486d0 - 2.185d0 * s
      bg  =  16.69d0 - 22.74d0 * s  + 5.779d0 * s2
      cg  = -25.59d0 + 29.71d0 * s  - 7.296d0 * s2
      dg  =  2.792d0 + 2.215d0 * s  + 0.422d0 * s2 - 0.104d0 * s3
      eg  =  0.807d0 + 2.005d0 * s
      esg =  3.841d0 + 0.316d0 * s
      gl  = pjgrvw (x, s, alg, beg, akg, bkg, ag, bg, cg,
     & dg, eg, esg)
 
      return
      end
 
C*********************************************************************
 
C...PYGRVM
C...Gives the GRV 94 M (MSbar) parton distribution function set
C...in parametrized form.
C...Authors: M. Glueck, E. Reya and A. Vogt.
 
      subroutine pjgrvm (x, q2, uv, dv, del, udb, sb, chm, bot, gl)
 
C...Double precision declaration.
      implicit double precision (a - z)
 
C...Common expressions.
      mu2  = 0.34d0
      lam2 = 0.248d0 * 0.248d0
      s  = log (log(q2/lam2) / log(mu2/lam2))
      ds = sqrt (s)
      s2 = s * s
      s3 = s2 * s
 
C...uv :
      nu  =  1.304d0 + 0.863d0 * s
      aku =  0.558d0 - 0.020d0 * s
      bku =          0.183d0 * s
      au  = -0.113d0 + 0.283d0 * s - 0.321d0 * s2
      bu  =  6.843d0 - 5.089d0 * s + 2.647d0 * s2 - 0.527d0 * s3
      cu  =  7.771d0 - 10.09d0 * s + 2.630d0 * s2
      du  =  3.315d0 + 1.145d0 * s - 0.583d0 * s2 + 0.154d0 * s3
      uv  = pjgrvv (x, nu, aku, bku, au, bu, cu, du)
 
C...dv :
      nd  =  0.102d0 - 0.017d0 * s + 0.005d0 * s2
      akd =  0.270d0 - 0.019d0 * s
      bkd =  0.260d0
      ad  =  2.393d0 + 6.228d0 * s - 0.881d0 * s2
      bd  =  46.06d0 + 4.673d0 * s - 14.98d0 * s2 + 1.331d0 * s3
      cd  =  17.83d0 - 53.47d0 * s + 21.24d0 * s2
      dd  =  4.081d0 + 0.976d0 * s - 0.485d0 * s2 + 0.152d0 * s3
      dv  = pjgrvv (x, nd, akd, bkd, ad, bd, cd, dd)
 
C...del :
      ne  =  0.070d0 + 0.042d0 * s - 0.011d0 * s2 + 0.004d0 * s3
      ake =  0.409d0 - 0.007d0 * s
      bke =  0.782d0 + 0.082d0 * s
      ae  = -29.65d0 + 26.49d0 * s + 5.429d0 * s2
      be  =  90.20d0 - 74.97d0 * s + 4.526d0 * s2
      ce  =  0.0d0
      de  =  8.122d0 + 2.120d0 * s - 1.088d0 * s2 + 0.231d0 * s3
      del = pjgrvv (x, ne, ake, bke, ae, be, ce, de)
 
C...udb :
      alx =  0.877d0
      bex =  0.561d0
      akx =  0.275d0
      bkx =  0.0d0
      agx =  0.997d0
      bgx =  3.210d0 - 1.866d0 * s
      cx  =  7.300d0
      dx  =  9.010d0 + 0.896d0 * ds + 0.222d0 * s2
      ex  =  3.077d0 + 1.446d0 * s
      esx =  3.173d0 - 2.445d0 * ds + 2.207d0 * s
      udb = pjgrvw (x, s, alx, bex, akx, bkx, agx, bgx, cx,
     & dx, ex, esx)
 
C...sb :
      sts =  0d0
      als =  0.756d0
      bes =  0.216d0
      aks =  1.690d0 + 0.650d0 * ds - 0.922d0 * s
      as  = -4.329d0 + 1.131d0 * s
      bs  =  9.568d0 - 1.744d0 * s
      dst =  9.377d0 + 1.088d0 * ds - 1.320d0 * s + 0.130d0 * s2
      est =  3.031d0 + 1.639d0 * s
      ess =  5.837d0 + 0.815d0 * s
      sb  = pjgrvs (x, s, sts, als, bes, aks, as, bs, dst, est, ess)
 
C...cb :
      stc =  0.820d0
      alc =  0.98d0
      bec =  0d0
      akc = -0.625d0 - 0.523d0 * s
      ac  =  0d0
      bc  =  1.896d0 + 1.616d0 * s
      dct =  4.12d0  + 0.683d0 * s
      ect =  4.36d0  + 1.328d0 * s
      esc =  0.677d0 + 0.679d0 * s
      chm = pjgrvs (x, s, stc, alc, bec, akc, ac, bc, dct, ect, esc)
 
C...bb :
      stb =  1.297d0
      alb =  0.99d0
      beb =  0d0
      akb =          - 0.193d0 * s
      ab  =  0d0
      bb  =  0d0
      dbt =  3.447d0 + 0.927d0 * s
      ebt =  4.68d0  + 1.259d0 * s
      esb =  1.892d0 + 2.199d0 * s
      bot = pjgrvs (x, s, stb, alb, beb, akb, ab, bb, dbt, ebt, esb)
 
C...gl :
       alg =  1.014d0
       beg =  1.738d0
       akg =  1.724d0 + 0.157d0 * s
       bkg =  0.800d0 + 1.016d0 * s
       ag  =  7.517d0 - 2.547d0 * s
       bg  =  34.09d0 - 52.21d0 * ds + 17.47d0 * s
       cg  =  4.039d0 + 1.491d0 * s
       dg  =  3.404d0 + 0.830d0 * s
       eg  = -1.112d0 + 3.438d0 * s  - 0.302d0 * s2
       esg =  3.256d0 - 0.436d0 * s
       gl  = pjgrvw (x, s, alg, beg, akg, bkg, ag, bg, cg, dg, eg, esg)
 
       return
       end
 
C*********************************************************************
 
C...PYGRVD
C...Gives the GRV 94 D (DIS) parton distribution function set
C...in parametrized form.
C...Authors: M. Glueck, E. Reya and A. Vogt.
 
      subroutine pjgrvd (x, q2, uv, dv, del, udb, sb, chm, bot, gl)
 
C...Double precision declaration.
      implicit double precision (a - z)
 
C...Common expressions.
      mu2  = 0.34d0
      lam2 = 0.248d0 * 0.248d0
      s  = log (log(q2/lam2) / log(mu2/lam2))
      ds = sqrt (s)
      s2 = s * s
      s3 = s2 * s
 
C...uv :
      nu  =  2.484d0 + 0.116d0 * s + 0.093d0 * s2
      aku =  0.563d0 - 0.025d0 * s
      bku =  0.054d0 + 0.154d0 * s
      au  = -0.326d0 - 0.058d0 * s - 0.135d0 * s2
      bu  = -3.322d0 + 8.259d0 * s - 3.119d0 * s2 + 0.291d0 * s3
      cu  =  11.52d0 - 12.99d0 * s + 3.161d0 * s2
      du  =  2.808d0 + 1.400d0 * s - 0.557d0 * s2 + 0.119d0 * s3
      uv  = pjgrvv (x, nu, aku, bku, au, bu, cu, du)
 
C...dv :
      nd  =  0.156d0 - 0.017d0 * s
      akd =  0.299d0 - 0.022d0 * s
      bkd =  0.259d0 - 0.015d0 * s
      ad  =  3.445d0 + 1.278d0 * s + 0.326d0 * s2
      bd  = -6.934d0 + 37.45d0 * s - 18.95d0 * s2 + 1.463d0 * s3
      cd  =  55.45d0 - 69.92d0 * s + 20.78d0 * s2
      dd  =  3.577d0 + 1.441d0 * s - 0.683d0 * s2 + 0.179d0 * s3
      dv  = pjgrvv (x, nd, akd, bkd, ad, bd, cd, dd)
 
C...del :
      ne  =  0.099d0 + 0.019d0 * s + 0.002d0 * s2
      ake =  0.419d0 - 0.013d0 * s
      bke =  1.064d0 - 0.038d0 * s
      ae  = -44.00d0 + 98.70d0 * s - 14.79d0 * s2
      be  =  28.59d0 - 40.94d0 * s - 13.66d0 * s2 + 2.523d0 * s3
      ce  =  84.57d0 - 108.8d0 * s + 31.52d0 * s2
      de  =  7.469d0 + 2.480d0 * s - 0.866d0 * s2
      del = pjgrvv (x, ne, ake, bke, ae, be, ce, de)
 
C...udb :
      alx =  1.215d0
      bex =  0.466d0
      akx =  0.326d0 + 0.150d0 * s
      bkx =  0.956d0 + 0.405d0 * s
      agx =  0.272d0
      bgx =  3.794d0 - 2.359d0 * ds
      cx  =  2.014d0
      dx  =  7.941d0 + 0.534d0 * ds - 0.940d0 * s + 0.410d0 * s2
      ex  =  3.049d0 + 1.597d0 * s
      esx =  4.396d0 - 4.594d0 * ds + 3.268d0 * s
      udb = pjgrvw (x, s, alx, bex, akx, bkx, agx, bgx, cx,
     & dx, ex, esx)
 
C...sb :
      sts =  0d0
      als =  0.175d0
      bes =  0.344d0
      aks =  1.415d0 - 0.641d0 * ds
      as  =  0.580d0 - 9.763d0 * ds + 6.795d0 * s  - 0.558d0 * s2
      bs  =  5.617d0 + 5.709d0 * ds - 3.972d0 * s
      dst =  13.78d0 - 9.581d0 * s  + 5.370d0 * s2 - 0.996d0 * s3
      est =  4.546d0 + 0.372d0 * s2
      ess =  5.053d0 - 1.070d0 * s  + 0.805d0 * s2
      sb  = pjgrvs (x, s, sts, als, bes, aks, as, bs, dst, est, ess)
 
C...cb :
      stc =  0.820d0
      alc =  0.98d0
      bec =  0d0
      akc = -0.625d0 - 0.523d0 * s
      ac  =  0d0
      bc  =  1.896d0 + 1.616d0 * s
      dct =  4.12d0  + 0.683d0 * s
      ect =  4.36d0  + 1.328d0 * s
      esc =  0.677d0 + 0.679d0 * s
      chm = pjgrvs (x, s, stc, alc, bec, akc, ac, bc, dct, ect, esc)
 
C...bb :
      stb =  1.297d0
      alb =  0.99d0
      beb =  0d0
      akb =          - 0.193d0 * s
      ab  =  0d0
      bb  =  0d0
      dbt =  3.447d0 + 0.927d0 * s
      ebt =  4.68d0  + 1.259d0 * s
      esb =  1.892d0 + 2.199d0 * s
      bot = pjgrvs (x, s, stb, alb, beb, akb, ab, bb, dbt, ebt, esb)
 
C...gl :
      alg =  1.258d0
      beg =  1.846d0
      akg =  2.423d0
      bkg =  2.427d0 + 1.311d0 * s  - 0.153d0 * s2
      ag  =  25.09d0 - 7.935d0 * s
      bg  = -14.84d0 - 124.3d0 * ds + 72.18d0 * s
      cg  =  590.3d0 - 173.8d0 * s
      dg  =  5.196d0 + 1.857d0 * s
      eg  = -1.648d0 + 3.988d0 * s  - 0.432d0 * s2
      esg =  3.232d0 - 0.542d0 * s
      gl  = pjgrvw (x, s, alg, beg, akg, bkg, ag, bg, cg, dg, eg, esg)
 
      return
      end
 
C*********************************************************************
 
C...PYGRVV
C...Auxiliary for the GRV 94 parton distribution functions
C...for u and d valence and d-u sea.
C...Authors: M. Glueck, E. Reya and A. Vogt.
 
      function pjgrvv (x, n, ak, bk, a, b, c, d)
 
C...Double precision declaration.
      implicit double precision (a - z)
 
C...Evaluation.
      dx = sqrt (x)
      pjgrvv = n * x**ak * (1d0+ a*x**bk + x * (b + c*dx)) *
     & (1d0- x)**d
 
      return
      end
 
C*********************************************************************
 
C...PYGRVW
C...Auxiliary for the GRV 94 parton distribution functions
C...for d+u sea and gluon.
C...Authors: M. Glueck, E. Reya and A. Vogt.
 
      function pjgrvw (x, s, al, be, ak, bk, a, b, c, d, e, es)
 
C...Double precision declaration.
      implicit double precision (a - z)
 
C...Evaluation.
      lx = log (1d0/x)
      pjgrvw = (x**ak * (a + x * (b + x*c)) * lx**bk + s**al
     &     * exp (-e + sqrt (es * s**be * lx))) * (1d0- x)**d
 
      return
      end
 
C*********************************************************************
 
C...PYGRVS
C...Auxiliary for the GRV 94 parton distribution functions
C...for s, c and b sea.
C...Authors: M. Glueck, E. Reya and A. Vogt.
 
      function pjgrvs (x, s, sth, al, be, ak, ag, b, d, e, es)
 
C...Double precision declaration.
      implicit double precision (a - z)
 
C...Evaluation.
      if(s.le.sth) then
        pjgrvs = 0d0
      else
        dx = sqrt (x)
        lx = log (1d0/x)
        pjgrvs = (s - sth)**al / lx**ak * (1d0+ ag*dx + b*x) *
     &     (1d0- x)**d * exp (-e + sqrt (es * s**be * lx))
      endif
 
      return
      end
 
C*********************************************************************
 
      function pjctq2 (iset, iprt, x, q)
 
C...This routine gives the CTEQ 2 parton distribution function sets in
C...parametrized form. It is adapted from the revised parametrization
C...with extended range of November 12, 1993.
C...Authors: J. Botts, H.L. Lai, J.G. Morfin, J.F. Owens, J. Qiu,
C...W.K. Tung and H. Weerts.
      implicit double precision(a-h, o-z)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jydat2/
 
C...Data on Lambda values of fits, minimum Q and quark masses.
      dimension alm(6), qms(4:6)
      data alm / 0.213d0, 0.208d0, 0.208d0, 0.322d0, 0.190d0, 0.235d0 /
      data qmn / 1.60d0 /, (qms(i), i=4,6) / 1.60d0, 5.00d0, 180.0d0 /
      qms(6) = pmas(6,1)
 
C....Check flavour thresholds. Set up Qi for SB.
      ip = iabs(iprt)
      if (ip .ge. 4) then
        if (q .le. qms(ip)) then
          pjctq2 = 0.0d0
          return
        endif
        qi = qms(ip)
      else
        qi = qmn
      endif
 
C...Use "standard lambda" of parametrization program for expansion.
      alam = alm (iset)
      sbl = log(q/alam) / log(qi/alam)
      sb = log (sbl)
      sb2 = sb*sb
      sb3 = sb2*sb
 
C...Expansion for run le26 - CTEQ2M
      if (iset .eq. 1) then
      if (iprt .eq. 2) then
      a0=exp( 0.2143d+00+0.8417d+00*sb -0.2451d+01*sb2+0.9875d+00*sb3)
      a1= 0.5209d+00-0.2384d+00*sb +0.5086d+00*sb2-0.2123d+00*sb3
      a2= 0.3178d+01+0.5258d+01*sb -0.8102d+01*sb2+0.3334d+01*sb3
      a3=-0.8537d+00+0.5921d+01*sb -0.1007d+02*sb2+0.4146d+01*sb3
      a4= 0.1821d+01+0.2822d-01*sb +0.1662d+00*sb2-0.1058d+00*sb3
      a5= 0.0000d+00-0.1090d+01*sb +0.3136d+01*sb2-0.1301d+01*sb3
      elseif (iprt .eq. 1) then
      a0=exp(-0.1314d+01-0.1342d-01*sb +0.1136d+00*sb2-0.1557d+00*sb3)
      a1= 0.2780d+00+0.2558d-01*sb +0.4467d-02*sb2-0.2472d-02*sb3
      a2= 0.3672d+01+0.5324d+00*sb +0.3531d-01*sb2+0.7928d-03*sb3
      a3= 0.2957d+02-0.2000d+02*sb +0.5929d+01*sb2+0.3390d+00*sb3
      a4= 0.8069d+00-0.2877d+00*sb +0.3574d-01*sb2+0.5622d-02*sb3
      a5= 0.0000d+00+0.2287d+00*sb -0.4052d-01*sb2+0.5589d-01*sb3
      elseif (iprt .eq. 0) then
      a0=exp(-0.1059d+00-0.1461d+01*sb -0.2544d+00*sb2+0.4526d-01*sb3)
      a1=-0.2578d+00+0.1385d+00*sb -0.1383d+00*sb2+0.3811d-01*sb3
      a2= 0.5195d+01+0.9648d+00*sb -0.2103d+00*sb2-0.6701d-01*sb3
      a3= 0.5131d+01+0.2151d+01*sb -0.2880d+01*sb2+0.6608d+00*sb3
      a4= 0.1118d+01+0.2636d+00*sb -0.5140d+00*sb2+0.1613d+00*sb3
      a5= 0.0000d+00+0.2456d+01*sb -0.8741d+00*sb2+0.2136d+00*sb3
      elseif (iprt .eq. -1) then
      a0=exp(-0.2732d+00-0.3523d+01*sb +0.3657d+01*sb2-0.1415d+01*sb3)
      a1=-0.3807d+00+0.1211d+00*sb -0.1231d+00*sb2+0.3753d-01*sb3
      a2= 0.9698d+01-0.2596d+01*sb +0.2412d+01*sb2-0.9257d+00*sb3
      a3=-0.6165d+00+0.1120d+01*sb -0.1708d+01*sb2+0.6383d+00*sb3
      a4= 0.7292d-01-0.1339d+00*sb +0.2104d+00*sb2-0.7987d-01*sb3
      a5=-0.1370d+01+0.2452d+01*sb -0.1804d+01*sb2+0.6459d+00*sb3
      elseif (iprt .eq. -2) then
      a0=exp(-0.2319d+01-0.3182d+01*sb +0.3572d+01*sb2-0.1431d+01*sb3)
      a1=-0.2622d+00+0.3085d+00*sb -0.4394d+00*sb2+0.1496d+00*sb3
      a2= 0.9481d+01-0.3627d+01*sb +0.5640d+01*sb2-0.2265d+01*sb3
      a3= 0.5000d+02-0.1851d+02*sb +0.2640d+01*sb2-0.6001d+00*sb3
      a4= 0.1566d+01-0.7375d+00*sb +0.8736d+00*sb2-0.3449d+00*sb3
      a5=-0.7983d-01+0.3236d+01*sb -0.3373d+01*sb2+0.1236d+01*sb3
      elseif (iprt .eq. -3) then
      a0=exp(-0.1855d+01-0.5302d+01*sb +0.8433d+00*sb2-0.1236d+00*sb3)
      a1=-0.4000d-02-0.1345d+01*sb +0.1192d+01*sb2-0.3039d+00*sb3
      a2= 0.6870d+01+0.1246d+01*sb -0.8968d+00*sb2-0.9791d-01*sb3
      a3= 0.0000d+00+0.4616d+01*sb +0.1026d+02*sb2+0.2844d+02*sb3
      a4= 0.1000d-02+0.4098d+00*sb -0.4250d+00*sb2+0.1100d+00*sb3
      a5= 0.0000d+00-0.2151d+01*sb +0.2991d+01*sb2-0.7717d+00*sb3
      elseif (iprt .eq. -4) then
      a0=sb** 0.7722d+00*exp(-0.7241d+01-0.7885d-01*sb -0.1124d+01*sb2)
      a1=-0.3971d+00+0.9132d+00*sb -0.1175d+01*sb2+0.3573d+00*sb3
      a2= 0.6367d+01-0.6565d+01*sb +0.8114d+01*sb2-0.2666d+01*sb3
      a3= 0.2878d+02-0.2000d+02*sb +0.7000d+00*sb2+0.3000d+02*sb3
      a4= 0.1010d+00-0.4592d+00*sb +0.5877d+00*sb2-0.1472d+00*sb3
      a5= 0.1749d+00+0.3875d+01*sb -0.3768d+01*sb2+0.1316d+01*sb3
      elseif (iprt .eq. -5) then
      a0=sb** 0.1299d+00*exp(-0.4868d+01-0.4339d+01*sb +0.7080d+00*sb2)
      a1=-0.1705d+00-0.3381d+00*sb +0.5287d+00*sb2-0.2644d+00*sb3
      a2= 0.5610d+01-0.1365d+01*sb +0.1835d+01*sb2-0.5655d+00*sb3
      a3=-0.1001d+01+0.3044d+01*sb +0.2680d+01*sb2+0.1426d+02*sb3
      a4= 0.3814d-02+0.3430d+00*sb -0.6926d+00*sb2+0.3486d+00*sb3
      a5= 0.1156d+01+0.2016d+01*sb -0.1674d+01*sb2+0.5981d+00*sb3
      elseif (iprt .eq. -6) then
      a0=sb** 0.9819d+00*exp(-0.7859d+01+0.6819d+00*sb -0.3386d+01*sb2)
      a1=-0.1055d+00-0.1413d+01*sb +0.3451d+01*sb2-0.2466d+01*sb3
      a2= 0.4055d+01+0.8107d+01*sb -0.1576d+02*sb2+0.8094d+01*sb3
      a3= 0.3799d+01+0.9616d+01*sb -0.1984d+02*sb2+0.2641d+02*sb3
      a4= 0.3619d+00-0.8627d+00*sb -0.9390d-01*sb2+0.9196d+00*sb3
      a5= 0.3779d+01-0.6073d+01*sb +0.9999d+01*sb2-0.4304d+01*sb3
      endif
 
C...Expansion for run sa17 - CTEQ2MS
      elseif (iset .eq. 2) then
      if (iprt .eq. 2) then
      a0=exp( 0.2790d+00+0.7294d+00*sb -0.2202d+01*sb2+0.8599d+00*sb3)
      a1= 0.5380d+00-0.2261d+00*sb +0.4636d+00*sb2-0.1871d+00*sb3
      a2= 0.3259d+01+0.2141d+01*sb -0.2947d+01*sb2+0.1245d+01*sb3
      a3=-0.8390d+00+0.1448d+01*sb -0.2331d+01*sb2+0.8658d+00*sb3
      a4= 0.1847d+01-0.3943d+01*sb +0.5998d+01*sb2-0.2191d+01*sb3
      a5= 0.0000d+00-0.9719d+00*sb +0.2830d+01*sb2-0.1137d+01*sb3
      elseif (iprt .eq. 1) then
      a0=exp(-0.1318d+01+0.2328d-01*sb +0.5179d-01*sb2-0.1305d+00*sb3)
      a1= 0.2760d+00+0.4429d-01*sb -0.2626d-01*sb2+0.7143d-02*sb3
      a2= 0.3660d+01+0.5232d+00*sb +0.5491d-01*sb2-0.4115d-02*sb3
      a3= 0.2910d+02-0.2000d+02*sb +0.6631d+01*sb2-0.3050d-01*sb3
      a4= 0.8010d+00-0.2688d+00*sb +0.1051d-01*sb2+0.1195d-01*sb3
      a5= 0.0000d+00+0.2887d+00*sb -0.1398d+00*sb2+0.8194d-01*sb3
      elseif (iprt .eq. 0) then
      a0=exp(-0.1623d+01-0.7232d+00*sb +0.1889d+00*sb2+0.1140d+00*sb3)
      a1=-0.5000d+00+0.8611d-01*sb +0.2203d-01*sb2-0.1401d-01*sb3
      a2= 0.3821d+01+0.8976d+00*sb +0.1400d+00*sb2-0.9163d-01*sb3
      a3= 0.5809d+01-0.5060d+01*sb +0.3808d+00*sb2+0.2519d+00*sb3
      a4= 0.4500d+00-0.5121d+00*sb +0.1979d+00*sb2-0.2705d-01*sb3
      a5= 0.0000d+00+0.1210d+01*sb -0.2921d+00*sb2+0.1240d+00*sb3
      elseif (iprt .eq. -1) then
      a0=exp(-0.6986d-01-0.5954d+00*sb -0.1582d+01*sb2+0.5104d+00*sb3)
      a1=-0.8461d+00+0.2127d+00*sb +0.9425d-01*sb2-0.5264d-01*sb3
      a2= 0.1200d+02+0.1659d+01*sb -0.5354d+01*sb2+0.1795d+01*sb3
      a3= 0.2958d+02+0.3000d+02*sb +0.3000d+02*sb2-0.1965d+02*sb3
      a4= 0.4000d+01-0.4865d+00*sb +0.9460d+00*sb2+0.3432d+00*sb3
      a5=-0.3378d+01+0.1656d+01*sb +0.1123d+01*sb2-0.4667d+00*sb3
      elseif (iprt .eq. -2) then
      a0=exp(-0.1929d+01-0.2626d+01*sb +0.2926d+01*sb2-0.1297d+01*sb3)
      a1=-0.6627d+00+0.4561d+00*sb -0.3818d+00*sb2+0.1239d+00*sb3
      a2= 0.9506d+01-0.2724d+01*sb +0.4283d+01*sb2-0.1804d+01*sb3
      a3= 0.1897d+02+0.1642d+01*sb -0.8390d+01*sb2+0.3894d+01*sb3
      a4= 0.1024d+01-0.1786d+00*sb +0.4535d+00*sb2-0.2075d+00*sb3
      a5=-0.1746d+01+0.3572d+01*sb -0.2908d+01*sb2+0.1093d+01*sb3
      elseif (iprt .eq. -3) then
      a0=exp(-0.4913d+00-0.6866d+01*sb +0.1432d+01*sb2-0.1749d+00*sb3)
      a1=-0.1157d+00-0.1567d+01*sb +0.1439d+01*sb2-0.3724d+00*sb3
      a2= 0.7730d+01+0.9748d+00*sb -0.1157d+01*sb2-0.8358d-02*sb3
      a3=-0.6050d+00+0.1835d+01*sb +0.3788d+01*sb2+0.3000d+02*sb3
      a4= 0.1620d-08+0.4590d+00*sb -0.4070d+00*sb2+0.8900d-01*sb3
      a5=-0.7048d+00-0.2505d+01*sb +0.4000d+01*sb2-0.1161d+01*sb3
      elseif (iprt .eq. -4) then
      a0=sb** 0.7393d+00*exp(-0.6518d+01-0.3998d+00*sb -0.1111d+01*sb2)
      a1=-0.6482d+00+0.1125d+01*sb -0.1290d+01*sb2+0.3940d+00*sb3
      a2= 0.8487d+01-0.9235d+01*sb +0.9353d+01*sb2-0.2913d+01*sb3
      a3= 0.2265d+02-0.1999d+02*sb +0.4105d+01*sb2+0.2144d+02*sb3
      a4= 0.8990d-01-0.4372d+00*sb +0.5941d+00*sb2-0.1469d+00*sb3
      a5=-0.9690d+00+0.5068d+01*sb -0.4368d+01*sb2+0.1503d+01*sb3
      elseif (iprt .eq. -5) then
      a0=sb** 0.9880d+00*exp(-0.7180d+01-0.2494d+01*sb +0.3561d-01*sb2)
      a1=-0.4301d+00-0.2611d+00*sb +0.3914d+00*sb2-0.1638d+00*sb3
      a2= 0.5137d+01+0.1506d+01*sb -0.9588d+00*sb2-0.1596d+00*sb3
      a3= 0.1483d+02+0.2998d+02*sb +0.2357d+02*sb2-0.9353d+01*sb3
      a4= 0.2426d+00+0.1371d+00*sb -0.3791d+00*sb2+0.1948d+00*sb3
      a5= 0.1463d+01+0.1907d+00*sb +0.3557d+00*sb2+0.2097d-01*sb3
      elseif (iprt .eq. -6) then
      a0=sb** 0.1005d+01*exp(-0.5255d+01-0.9866d-01*sb -0.2737d+01*sb2)
      a1=-0.3140d+00-0.2055d+00*sb +0.5594d+00*sb2-0.2960d+00*sb3
      a2= 0.9227d+01-0.4569d+01*sb -0.9724d+01*sb2+0.1026d+02*sb3
      a3= 0.1131d+02-0.1972d+02*sb -0.1107d+02*sb2+0.2311d+02*sb3
      a4= 0.1488d+01+0.1737d+01*sb +0.4323d+01*sb2-0.9925d+01*sb3
      a5= 0.1895d+01-0.7350d+00*sb +0.3780d+01*sb2-0.1408d+01*sb3
      endif
 
      elseif (iset .eq. 3) then
C...Expansion for run fa06 - CTEQ2MF
      if (iprt .eq. 2) then
      a0=exp(-0.7913d+00-0.2789d+01*sb -0.7289d-01*sb2+0.1770d+00*sb3)
      a1= 0.4942d+00-0.7886d-01*sb +0.9057d-01*sb2-0.5259d-01*sb3
      a2= 0.3727d+01+0.1089d+01*sb -0.1004d+01*sb2+0.4345d+00*sb3
      a3= 0.1944d+01+0.7846d+01*sb +0.7984d+01*sb2+0.5548d+01*sb3
      a4= 0.2940d-02+0.8428d-04*sb +0.1266d+00*sb2-0.3517d-01*sb3
      a5=-0.1060d+00-0.1192d-01*sb +0.1130d+01*sb2-0.4527d+00*sb3
      elseif (iprt .eq. 1) then
      a0=exp(-0.1344d+01+0.7859d-02*sb +0.4623d-01*sb2-0.1273d+00*sb3)
      a1= 0.2760d+00+0.4201d-01*sb -0.1795d-01*sb2+0.3212d-02*sb3
      a2= 0.3660d+01+0.5247d+00*sb +0.4405d-01*sb2+0.1391d-02*sb3
      a3= 0.2981d+02-0.2000d+02*sb +0.6566d+01*sb2+0.2479d-01*sb3
      a4= 0.7950d+00-0.2732d+00*sb +0.2470d-01*sb2+0.6157d-02*sb3
      a5= 0.0000d+00+0.2793d+00*sb -0.9197d-01*sb2+0.5953d-01*sb3
      elseif (iprt .eq. 0) then
      a0=exp( 0.9746d+00-0.3252d+01*sb +0.1664d+01*sb2-0.6410d+00*sb3)
      a1=-0.5271d-02-0.3198d+00*sb +0.1279d+00*sb2-0.1256d-02*sb3
      a2= 0.5740d+01-0.3139d+01*sb +0.3841d+01*sb2-0.1415d+01*sb3
      a3= 0.7161d-01-0.4363d+01*sb +0.4925d+01*sb2-0.1614d+01*sb3
      a4= 0.1860d+01+0.1342d+01*sb -0.2234d+01*sb2+0.1047d+01*sb3
      a5= 0.7409d-01+0.2390d+01*sb -0.1457d+01*sb2+0.5853d+00*sb3
      elseif (iprt .eq. -1) then
      a0=exp(-0.8454d+00-0.3334d+01*sb +0.3591d+01*sb2-0.1485d+01*sb3)
      a1=-0.2826d-02-0.2810d+00*sb -0.3809d-01*sb2+0.6585d-01*sb3
      a2= 0.9139d+01-0.2811d+01*sb +0.4730d+01*sb2-0.2157d+01*sb3
      a3=-0.3120d+00+0.1217d+01*sb -0.1726d+01*sb2+0.6220d+00*sb3
      a4= 0.1793d-01-0.4608d-01*sb +0.5294d-01*sb2-0.1709d-01*sb3
      a5=-0.1471d+00+0.1104d+01*sb -0.1358d+01*sb2+0.7200d+00*sb3
      elseif (iprt .eq. -2) then
      a0=exp(-0.1398d+01-0.3536d+01*sb +0.3849d+01*sb2-0.1549d+01*sb3)
      a1=-0.1332d-01-0.2155d-01*sb -0.3404d+00*sb2+0.1569d+00*sb3
      a2= 0.9981d+01-0.3499d+01*sb +0.5448d+01*sb2-0.2198d+01*sb3
      a3= 0.3736d+02-0.2000d+02*sb +0.6675d+01*sb2-0.7276d+00*sb3
      a4= 0.1705d+01-0.1013d+01*sb +0.1122d+01*sb2-0.4057d+00*sb3
      a5=-0.1189d-01+0.2698d+01*sb -0.3429d+01*sb2+0.1389d+01*sb3
      elseif (iprt .eq. -3) then
      a0=exp(-0.2979d+01-0.6085d+01*sb +0.2428d+01*sb2-0.6482d+00*sb3)
      a1=-0.1372d+00-0.1281d+00*sb +0.1587d+00*sb2-0.9637d-01*sb3
      a2= 0.7009d+01-0.1609d+01*sb +0.2765d+01*sb2-0.1177d+01*sb3
      a3= 0.1308d+01+0.9583d+01*sb +0.2360d+02*sb2+0.2999d+02*sb3
      a4= 0.2509d-01+0.2106d+00*sb -0.4405d+00*sb2+0.2075d+00*sb3
      a5=-0.2069d-01+0.1971d+01*sb -0.1615d+01*sb2+0.6039d+00*sb3
      elseif (iprt .eq. -4) then
      a0=sb** 0.8072d+00*exp(-0.6920d+01-0.5031d+00*sb -0.9965d+00*sb2)
      a1=-0.2118d+00+0.7930d+00*sb -0.1101d+01*sb2+0.3302d+00*sb3
      a2= 0.8039d+01-0.7170d+01*sb +0.8657d+01*sb2-0.2893d+01*sb3
      a3= 0.2926d+02-0.1993d+02*sb +0.1841d+01*sb2+0.2996d+02*sb3
      a4= 0.1339d+00-0.5531d+00*sb +0.6505d+00*sb2-0.1595d+00*sb3
      a5= 0.7439d+00+0.3307d+01*sb -0.3284d+01*sb2+0.1152d+01*sb3
      elseif (iprt .eq. -5) then
      a0=sb** 0.9925d+00*exp(-0.2190d+01-0.3393d+01*sb -0.8631d+00*sb2)
      a1=-0.1261d+00-0.2368d+00*sb +0.4143d+00*sb2-0.1577d+00*sb3
      a2= 0.4585d+01+0.5227d+01*sb -0.3248d+01*sb2-0.2599d+00*sb3
      a3=-0.1094d+01+0.4927d+00*sb -0.9921d+00*sb2+0.3138d+01*sb3
      a4= 0.1396d+00+0.2562d+00*sb +0.1844d+00*sb2-0.1599d+00*sb3
      a5= 0.8621d+00+0.4715d+00*sb +0.2547d+01*sb2-0.8429d+00*sb3
      elseif (iprt .eq. -6) then
      a0=sb** 0.1016d+01*exp(-0.5397d+01-0.1979d+01*sb -0.2441d+00*sb2)
      a1=-0.1426d+00-0.2861d+00*sb +0.7434d+00*sb2-0.5214d+00*sb3
      a2= 0.6363d+01+0.4028d+00*sb -0.8356d+01*sb2+0.6814d+01*sb3
      a3=-0.2526d+00+0.2425d+01*sb -0.1407d+02*sb2+0.3000d+02*sb3
      a4= 0.1125d+00-0.1089d+01*sb +0.9977d+01*sb2+0.1000d+02*sb3
      a5= 0.2669d+01-0.6366d+00*sb +0.4355d+01*sb2-0.2919d+01*sb3
      endif
 
      elseif (iset .eq. 4) then
C...Expansion for run ll25 - CTEQ2ML
      if (iprt .eq. 2) then
      a0=exp( 0.3760d+00+0.5491d+00*sb -0.1845d+01*sb2+0.6803d+00*sb3)
      a1= 0.5650d+00-0.1953d+00*sb +0.3761d+00*sb2-0.1419d+00*sb3
      a2= 0.3464d+01+0.3817d+01*sb -0.5384d+01*sb2+0.2057d+01*sb3
      a3=-0.5850d+00+0.5566d+01*sb -0.9000d+01*sb2+0.3433d+01*sb3
      a4= 0.2322d+01-0.1431d+00*sb +0.3901d+00*sb2-0.1678d+00*sb3
      a5= 0.0000d+00-0.7370d+00*sb +0.2310d+01*sb2-0.8743d+00*sb3
      elseif (iprt .eq. 1) then
      a0=exp(-0.1324d+01+0.1169d-01*sb +0.1969d-01*sb2-0.7583d-01*sb3)
      a1= 0.2890d+00+0.5832d-01*sb -0.2921d-01*sb2+0.4701d-02*sb3
      a2= 0.3580d+01+0.5291d+00*sb -0.5662d-02*sb2+0.2746d-01*sb3
      a3= 0.3021d+02-0.1999d+02*sb +0.6250d+01*sb2-0.3035d+00*sb3
      a4= 0.7990d+00-0.2531d+00*sb +0.5556d-02*sb2+0.8272d-02*sb3
      a5= 0.0000d+00+0.3674d+00*sb -0.1383d+00*sb2+0.4665d-01*sb3
      elseif (iprt .eq. 0) then
      a0=exp(-0.1920d+00-0.7015d+00*sb -0.9113d+00*sb2+0.2352d+00*sb3)
      a1=-0.2120d+00+0.1133d-01*sb -0.1553d-01*sb2+0.2822d-02*sb3
      a2= 0.4549d+01+0.1250d+01*sb -0.4647d+00*sb2+0.9617d-01*sb3
      a3= 0.1197d+02-0.4156d+01*sb +0.1413d+00*sb2+0.1607d+00*sb3
      a4= 0.1616d+01+0.1082d+00*sb -0.6651d+00*sb2+0.2356d+00*sb3
      a5= 0.0000d+00+0.1824d+01*sb -0.2063d+00*sb2+0.1148d-01*sb3
      elseif (iprt .eq. -1) then
      a0=exp(-0.1388d+01-0.7408d+00*sb -0.6454d+00*sb2+0.2373d+00*sb3)
      a1=-0.2928d+00-0.1726d-01*sb +0.4033d-01*sb2-0.2514d-01*sb3
      a2= 0.9975d+01-0.2048d+01*sb -0.6060d+00*sb2+0.5225d+00*sb3
      a3= 0.2687d+02-0.4683d+01*sb -0.1999d+02*sb2+0.1188d+02*sb3
      a4= 0.4000d+01-0.6773d+00*sb +0.4301d+00*sb2+0.4524d+00*sb3
      a5=-0.7164d+00+0.7488d+00*sb +0.5766d+00*sb2-0.2609d+00*sb3
      elseif (iprt .eq. -2) then
      a0=exp(-0.2272d+01-0.2998d+01*sb +0.3282d+01*sb2-0.1203d+01*sb3)
      a1=-0.2062d+00+0.3320d+00*sb -0.5074d+00*sb2+0.1655d+00*sb3
      a2= 0.9667d+01-0.3497d+01*sb +0.5271d+01*sb2-0.1984d+01*sb3
      a3= 0.4996d+02-0.3241d+01*sb -0.1425d+02*sb2+0.3849d+01*sb3
      a4= 0.1619d+01-0.5354d+00*sb +0.5753d+00*sb2-0.2238d+00*sb3
      a5= 0.8755d-01+0.3195d+01*sb -0.3496d+01*sb2+0.1197d+01*sb3
      elseif (iprt .eq. -3) then
      a0=exp(-0.1864d+01-0.5258d+01*sb +0.1034d+01*sb2-0.1550d+00*sb3)
      a1= 0.1000d-02-0.1090d+01*sb +0.8345d+00*sb2-0.1887d+00*sb3
      a2= 0.6898d+01-0.4951d+00*sb +0.4279d+00*sb2-0.2727d+00*sb3
      a3= 0.0000d+00+0.4322d+01*sb +0.8181d+01*sb2+0.2309d+02*sb3
      a4= 0.1000d-02+0.3550d+00*sb -0.3220d+00*sb2+0.7294d-01*sb3
      a5= 0.0000d+00-0.1347d+01*sb +0.1896d+01*sb2-0.4491d+00*sb3
      elseif (iprt .eq. -4) then
      a0=sb** 0.7528d+00*exp(-0.7684d+01+0.6791d-01*sb -0.9094d+00*sb2)
      a1=-0.3732d+00+0.8408d+00*sb -0.1020d+01*sb2+0.3046d+00*sb3
      a2= 0.4984d+01-0.5534d+01*sb +0.6418d+01*sb2-0.1856d+01*sb3
      a3= 0.3761d+02-0.1999d+02*sb -0.3358d+01*sb2+0.2999d+02*sb3
      a4= 0.1161d+00-0.4680d+00*sb +0.5567d+00*sb2-0.1633d+00*sb3
      a5= 0.3028d+00+0.3339d+01*sb -0.3004d+01*sb2+0.9160d+00*sb3
      elseif (iprt .eq. -5) then
      a0=sb** 0.1011d+01*exp(-0.7217d+01-0.2288d+01*sb +0.3450d+00*sb2)
      a1=-0.1955d+00-0.3371d+00*sb +0.5111d+00*sb2-0.2210d+00*sb3
      a2= 0.4302d+01-0.1214d+01*sb +0.3104d+01*sb2-0.1408d+01*sb3
      a3= 0.1487d+02+0.1549d+02*sb +0.2875d+02*sb2-0.1922d+02*sb3
      a4= 0.8935d-02+0.3571d+00*sb -0.6668d+00*sb2+0.3037d+00*sb3
      a5= 0.1570d+01+0.7105d+00*sb -0.6070d+00*sb2+0.3796d+00*sb3
      elseif (iprt .eq. -6) then
      a0=sb** 0.9986d+00*exp(-0.5847d+01-0.2798d+00*sb -0.9882d+00*sb2)
      a1=-0.2154d+00-0.8282d-01*sb +0.3611d-01*sb2+0.2623d-01*sb3
      a2= 0.3250d+01+0.9635d+01*sb -0.1274d+02*sb2+0.4453d+01*sb3
      a3=-0.2594d+01+0.9097d+01*sb +0.1581d+02*sb2-0.9123d+01*sb3
      a4= 0.1768d+01-0.2749d+01*sb +0.9999d+01*sb2+0.9995d+01*sb3
      a5= 0.2521d+01-0.1802d-01*sb +0.4820d+00*sb2+0.2004d+00*sb3
      endif
 
      elseif (iset .eq. 5) then
C...Expansion for run lo24 - CTEQ2L
      if (iprt .eq. 2) then
      a0=exp( 0.7248d-01+0.3941d+00*sb -0.1772d+01*sb2+0.7629d+00*sb3)
      a1= 0.4964d+00-0.1224d+00*sb +0.3646d+00*sb2-0.1685d+00*sb3
      a2= 0.3000d+01+0.2780d+01*sb -0.4028d+01*sb2+0.1816d+01*sb3
      a3=-0.1064d+01+0.3062d+01*sb -0.5927d+01*sb2+0.2785d+01*sb3
      a4= 0.3193d+01+0.1499d+01*sb -0.2765d+01*sb2+0.1019d+01*sb3
      a5= 0.1524d-01-0.4541d+00*sb +0.2281d+01*sb2-0.1033d+01*sb3
      elseif (iprt .eq. 1) then
      a0=exp(-0.1794d+01-0.2055d+00*sb -0.3350d-01*sb2-0.5084d-01*sb3)
      a1= 0.1748d+00+0.4637d-01*sb -0.2048d-01*sb2+0.2596d-02*sb3
      a2= 0.3321d+01+0.6253d+00*sb +0.2148d-01*sb2+0.1288d-01*sb3
      a3= 0.4355d+02-0.2000d+02*sb +0.5486d+01*sb2+0.1536d+00*sb3
      a4= 0.9586d+00-0.3217d+00*sb +0.4458d-01*sb2-0.1404d-03*sb3
      a5=-0.6595d-02+0.3499d+00*sb -0.7048d-01*sb2+0.2619d-01*sb3
      elseif (iprt .eq. 0) then
      a0=exp(-0.6194d+00-0.2643d+00*sb -0.1875d+01*sb2+0.6011d+00*sb3)
      a1=-0.2600d+00+0.8704d-01*sb -0.7375d-01*sb2+0.1876d-01*sb3
      a2= 0.4620d+01+0.1578d+01*sb -0.8411d+00*sb2+0.1527d+00*sb3
      a3= 0.1604d+02-0.1230d+02*sb +0.6939d+01*sb2-0.2012d+01*sb3
      a4= 0.1255d+01+0.4769d+00*sb -0.9915d+00*sb2+0.3439d+00*sb3
      a5= 0.1116d-02+0.2409d+01*sb -0.4442d+00*sb2+0.3431d-01*sb3
      elseif (iprt .eq. -1) then
      a0=exp(-0.1571d+01-0.1905d+00*sb -0.8672d+00*sb2+0.2070d+00*sb3)
      a1=-0.3266d+00+0.6428d-01*sb -0.8694d-01*sb2+0.1778d-01*sb3
      a2= 0.8921d+01-0.5010d+00*sb -0.9658d+00*sb2+0.3893d+00*sb3
      a3= 0.1329d+02+0.4652d+01*sb -0.2000d+02*sb2+0.1001d+02*sb3
      a4= 0.3283d+01-0.3400d+00*sb -0.1957d+00*sb2+0.8063d+00*sb3
      a5=-0.5701d+00+0.4042d+00*sb +0.5239d+00*sb2-0.1665d+00*sb3
      elseif (iprt .eq. -2) then
      a0=exp(-0.2281d+01-0.2768d+01*sb +0.3137d+01*sb2-0.1278d+01*sb3)
      a1=-0.2624d+00+0.4142d+00*sb -0.5936d+00*sb2+0.1937d+00*sb3
      a2= 0.9438d+01-0.3179d+01*sb +0.5107d+01*sb2-0.2179d+01*sb3
      a3= 0.5000d+02-0.1802d+02*sb -0.7515d+01*sb2+0.2991d+01*sb3
      a4= 0.1809d+01-0.9121d+00*sb +0.8854d+00*sb2-0.3582d+00*sb3
      a5= 0.4056d-01+0.3033d+01*sb -0.3431d+01*sb2+0.1253d+01*sb3
      elseif (iprt .eq. -3) then
      a0=exp(-0.2318d+01-0.4104d+01*sb -0.1502d+00*sb2+0.1693d+00*sb3)
      a1=-0.2251d-01-0.1101d+01*sb +0.1037d+01*sb2-0.3290d+00*sb3
      a2= 0.6989d+01+0.1794d+01*sb -0.1811d+01*sb2+0.3061d+00*sb3
      a3= 0.7972d+00+0.7806d+01*sb +0.1869d+02*sb2+0.2999d+02*sb3
      a4= 0.4795d-01+0.1622d+00*sb -0.3977d+00*sb2+0.1920d+00*sb3
      a5=-0.5275d-01-0.2616d+01*sb +0.3076d+01*sb2-0.7425d+00*sb3
      elseif (iprt .eq. -4) then
      a0=sb** 0.8431d+00*exp(-0.6539d+01-0.1875d+00*sb -0.1346d+01*sb2)
      a1=-0.4970d+00+0.9062d+00*sb -0.1169d+01*sb2+0.3703d+00*sb3
      a2= 0.4939d+01-0.2995d+01*sb +0.4483d+01*sb2-0.1704d+01*sb3
      a3= 0.3113d+02-0.1997d+02*sb +0.1540d+01*sb2+0.3000d+02*sb3
      a4= 0.1349d+00-0.5418d+00*sb +0.6142d+00*sb2-0.1360d+00*sb3
      a5=-0.8590d+00+0.3956d+01*sb -0.3612d+01*sb2+0.1401d+01*sb3
      elseif (iprt .eq. -5) then
      a0=sb** 0.2639d-01*exp(-0.2099d+01-0.2681d+01*sb +0.2925d+00*sb2)
      a1=-0.2243d+00-0.5343d-01*sb -0.1953d-01*sb2+0.1586d-01*sb3
      a2= 0.4294d+01+0.1102d+01*sb -0.1822d+00*sb2-0.2481d+00*sb3
      a3=-0.9998d+00+0.8275d-01*sb +0.5494d+00*sb2-0.1982d+00*sb3
      a4= 0.5904d-04+0.9222d-01*sb -0.9293d-01*sb2+0.9159d-01*sb3
      a5= 0.2657d+00+0.1770d+01*sb -0.7111d+00*sb2+0.2525d+00*sb3
      elseif (iprt .eq. -6) then
      a0=sb** 0.1009d+01*exp(-0.7032d+01+0.4562d+01*sb -0.9081d+01*sb2)
      a1=-0.1412d+00-0.5076d+00*sb +0.9513d+00*sb2-0.4326d+00*sb3
      a2= 0.5385d+01+0.3023d+01*sb -0.1162d+02*sb2+0.7006d+01*sb3
      a3= 0.4997d+01-0.1600d+02*sb +0.1342d+02*sb2+0.1197d+02*sb3
      a4= 0.5825d+00+0.3994d+00*sb -0.1255d+01*sb2+0.6486d+00*sb3
      a5= 0.3365d+01-0.4026d+01*sb +0.8385d+01*sb2-0.2260d+01*sb3
      endif
 
      elseif (iset .eq. 6) then
C...Expansion for run da06 - CTEQ2D
      if (iprt .eq. 2) then
      a0=exp( 0.1590d+00+0.5580d+00*sb -0.1838d+01*sb2+0.7018d+00*sb3)
      a1= 0.5110d+00-0.1625d+00*sb +0.3547d+00*sb2-0.1412d+00*sb3
      a2= 0.3158d+01+0.3962d+01*sb -0.5866d+01*sb2+0.2375d+01*sb3
      a3=-0.6000d+00+0.6144d+01*sb -0.1056d+02*sb2+0.4345d+01*sb3
      a4= 0.2306d+01-0.4669d-01*sb +0.2711d+00*sb2-0.1640d+00*sb3
      a5= 0.0000d+00-0.6638d+00*sb +0.2239d+01*sb2-0.8843d+00*sb3
      elseif (iprt .eq. 1) then
      a0=exp(-0.1182d+01+0.1449d+00*sb +0.2753d-01*sb2-0.1009d+00*sb3)
      a1= 0.2540d+00+0.2686d-01*sb -0.1546d-01*sb2+0.5396d-02*sb3
      a2= 0.3442d+01+0.5576d+00*sb +0.1937d-01*sb2+0.6696d-02*sb3
      a3= 0.2545d+02-0.2000d+02*sb +0.7355d+01*sb2-0.7058d+00*sb3
      a4= 0.9170d+00-0.3090d+00*sb +0.1705d-01*sb2+0.8534d-02*sb3
      a5= 0.0000d+00+0.1449d+00*sb -0.7821d-01*sb2+0.6405d-01*sb3
      elseif (iprt .eq. 0) then
      a0=exp(-0.3410d+00-0.9613d+00*sb -0.4969d+00*sb2+0.9360d-01*sb3)
      a1=-0.2400d+00+0.1473d+00*sb -0.1593d+00*sb2+0.4538d-01*sb3
      a2= 0.4841d+01+0.9311d+00*sb +0.1601d-03*sb2-0.1331d+00*sb3
      a3= 0.7427d+01-0.1397d+01*sb +0.1489d+00*sb2-0.2848d+00*sb3
      a4= 0.9600d+00+0.3697d+00*sb -0.4246d+00*sb2+0.1032d+00*sb3
      a5= 0.0000d+00+0.2484d+01*sb -0.9908d+00*sb2+0.2568d+00*sb3
      elseif (iprt .eq. -1) then
      a0=exp( 0.1176d+00-0.3418d+01*sb +0.3529d+01*sb2-0.1367d+01*sb3)
      a1=-0.3654d+00+0.1914d+00*sb -0.2192d+00*sb2+0.6933d-01*sb3
      a2= 0.1099d+02-0.4281d+01*sb +0.3729d+01*sb2-0.1254d+01*sb3
      a3=-0.7514d+00+0.7696d+00*sb -0.1134d+01*sb2+0.4245d+00*sb3
      a4= 0.7690d-01-0.6558d-01*sb +0.8726d-01*sb2-0.3345d-01*sb3
      a5=-0.1447d+01+0.2617d+01*sb -0.2094d+01*sb2+0.7536d+00*sb3
      elseif (iprt .eq. -2) then
      a0=exp(-0.2412d+01-0.2522d+01*sb +0.3126d+01*sb2-0.1305d+01*sb3)
      a1=-0.2353d+00+0.3118d+00*sb -0.4864d+00*sb2+0.1689d+00*sb3
      a2= 0.9017d+01-0.2437d+01*sb +0.4659d+01*sb2-0.2044d+01*sb3
      a3= 0.5000d+02-0.1158d+02*sb -0.9260d+01*sb2+0.2847d+01*sb3
      a4= 0.1726d+01-0.6849d+00*sb +0.7864d+00*sb2-0.3300d+00*sb3
      a5= 0.5080d-01+0.2858d+01*sb -0.3297d+01*sb2+0.1246d+01*sb3
      elseif (iprt .eq. -3) then
      a0=exp(-0.1966d+01-0.4405d+01*sb +0.2436d+00*sb2+0.4576d-01*sb3)
      a1=-0.4000d-02-0.1229d+01*sb +0.1118d+01*sb2-0.2988d+00*sb3
      a2= 0.6902d+01+0.1266d+01*sb -0.1068d+01*sb2+0.3062d-01*sb3
      a3= 0.0000d+00+0.3987d+01*sb +0.9389d+01*sb2+0.1881d+02*sb3
      a4= 0.1000d-02+0.3528d+00*sb -0.4201d+00*sb2+0.1248d+00*sb3
      a5= 0.0000d+00-0.2149d+01*sb +0.2925d+01*sb2-0.7609d+00*sb3
      elseif (iprt .eq. -4) then
      a0=sb** 0.7561d+00*exp(-0.6960d+01+0.5634d-01*sb -0.1170d+01*sb2)
      a1=-0.4232d+00+0.9269d+00*sb -0.1161d+01*sb2+0.3470d+00*sb3
      a2= 0.6057d+01-0.5790d+01*sb +0.7352d+01*sb2-0.2435d+01*sb3
      a3= 0.2941d+02-0.1999d+02*sb -0.8345d+00*sb2+0.3000d+02*sb3
      a4= 0.1069d+00-0.4620d+00*sb +0.5614d+00*sb2-0.1336d+00*sb3
      a5=-0.1865d+00+0.3953d+01*sb -0.3791d+01*sb2+0.1315d+01*sb3
      elseif (iprt .eq. -5) then
      a0=sb** 0.5661d-02*exp(-0.2123d+01-0.3026d+01*sb +0.1912d+00*sb2)
      a1=-0.2011d+00-0.1338d-01*sb -0.3974d-01*sb2+0.1948d-01*sb3
      a2= 0.4906d+01+0.1740d+01*sb -0.1387d+01*sb2+0.1263d+00*sb3
      a3=-0.1000d+01+0.5767d-01*sb +0.6377d+00*sb2+0.4736d-01*sb3
      a4= 0.5927d-04+0.1039d+00*sb -0.9797d-01*sb2+0.6881d-01*sb3
      a5= 0.4017d+00+0.1981d+01*sb -0.7758d+00*sb2+0.2916d+00*sb3
      elseif (iprt .eq. -6) then
      a0=sb** 0.1008d+01*exp(-0.7211d+01+0.3273d+01*sb -0.6979d+01*sb2)
      a1=-0.1026d+00-0.4948d+00*sb +0.1188d+01*sb2-0.8016d+00*sb3
      a2= 0.5397d+01+0.2135d+01*sb -0.9531d+01*sb2+0.6115d+01*sb3
      a3= 0.4966d+01-0.1111d+02*sb +0.4732d+01*sb2+0.1568d+02*sb3
      a4= 0.5345d+00-0.1935d+00*sb +0.5816d+00*sb2-0.6794d+00*sb3
      a5= 0.3569d+01-0.3477d+01*sb +0.8756d+01*sb2-0.4139d+01*sb3
      endif
      endif
 
C...Calculation of x * f(x, Q).
      pjctq2 = max(0.d0, a0 *(x**a1) *((1.d0-x)**a2) *(1.d0+a3*(x**a4))
     &                 *(log(1.d0+1.d0/x))**a5 )
 
      return
      end
 
c***********************************************************************

      subroutine ehlq(x,q2,xppr)

C...Proton structure functions from Eichten, Hinchliffe, Lane, Quigg.
C...Allowed variable range: 5 GeV^2 < Q^2 < 1E8 GeV^2; 1E-4 < x < 1

      implicit double precision(a-h, o-z)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      save /jydat2/
      save /pjpars/,/pjint1/
      dimension xppr(-6:6),xq(9),tx(6),tt(6),nehlq(8,2),cehlq(6,6,2,8,2)

C...The following data lines are coefficients needed in the
C...Eichten, Hinchliffe, Lane, Quigg proton structure function
C...parametrizations, see below.
C...Powers of 1-x in different cases.
      data nehlq/3,4,7,5,7,7,7,7,3,4,7,6,7,7,7,7/
C...Expansion coefficients for up valence quark distribution.
      data (((cehlq(ix,it,nx,1,1),ix=1,6),it=1,6),nx=1,2)/
     1 7.677d-01,-2.087d-01,-3.303d-01,-2.517d-02,-1.570d-02,-1.000d-04,
     2-5.326d-01,-2.661d-01, 3.201d-01, 1.192d-01, 2.434d-02, 7.620d-03,
     3 2.162d-01, 1.881d-01,-8.375d-02,-6.515d-02,-1.743d-02,-5.040d-03,
     4-9.211d-02,-9.952d-02, 1.373d-02, 2.506d-02, 8.770d-03, 2.550d-03,
     5 3.670d-02, 4.409d-02, 9.600d-04,-7.960d-03,-3.420d-03,-1.050d-03,
     6-1.549d-02,-2.026d-02,-3.060d-03, 2.220d-03, 1.240d-03, 4.100d-04,
     1 2.395d-01, 2.905d-01, 9.778d-02, 2.149d-02, 3.440d-03, 5.000d-04,
     2 1.751d-02,-6.090d-03,-2.687d-02,-1.916d-02,-7.970d-03,-2.750d-03,
     3-5.760d-03,-5.040d-03, 1.080d-03, 2.490d-03, 1.530d-03, 7.500d-04,
     4 1.740d-03, 1.960d-03, 3.000d-04,-3.400d-04,-2.900d-04,-1.800d-04,
     5-5.300d-04,-6.400d-04,-1.700d-04, 4.000d-05, 6.000d-05, 4.000d-05,
     6 1.700d-04, 2.200d-04, 8.000d-05, 1.000d-05,-1.000d-05,-1.000d-05/
      data (((cehlq(ix,it,nx,1,2),ix=1,6),it=1,6),nx=1,2)/
     1 7.237d-01,-2.189d-01,-2.995d-01,-1.909d-02,-1.477d-02, 2.500d-04,
     2-5.314d-01,-2.425d-01, 3.283d-01, 1.119d-01, 2.223d-02, 7.070d-03,
     3 2.289d-01, 1.890d-01,-9.859d-02,-6.900d-02,-1.747d-02,-5.080d-03,
     4-1.041d-01,-1.084d-01, 2.108d-02, 2.975d-02, 9.830d-03, 2.830d-03,
     5 4.394d-02, 5.116d-02,-1.410d-03,-1.055d-02,-4.230d-03,-1.270d-03,
     6-1.991d-02,-2.539d-02,-2.780d-03, 3.430d-03, 1.720d-03, 5.500d-04,
     1 2.410d-01, 2.884d-01, 9.369d-02, 1.900d-02, 2.530d-03, 2.400d-04,
     2 1.765d-02,-9.220d-03,-3.037d-02,-2.085d-02,-8.440d-03,-2.810d-03,
     3-6.450d-03,-5.260d-03, 1.720d-03, 3.110d-03, 1.830d-03, 8.700d-04,
     4 2.120d-03, 2.320d-03, 2.600d-04,-4.900d-04,-3.900d-04,-2.300d-04,
     5-6.900d-04,-8.200d-04,-2.000d-04, 7.000d-05, 9.000d-05, 6.000d-05,
     6 2.400d-04, 3.100d-04, 1.100d-04, 0.000d+00,-2.000d-05,-2.000d-05/
C...Expansion coefficients for down valence quark distribution.
      data (((cehlq(ix,it,nx,2,1),ix=1,6),it=1,6),nx=1,2)/
     1 3.813d-01,-8.090d-02,-1.634d-01,-2.185d-02,-8.430d-03,-6.200d-04,
     2-2.948d-01,-1.435d-01, 1.665d-01, 6.638d-02, 1.473d-02, 4.080d-03,
     3 1.252d-01, 1.042d-01,-4.722d-02,-3.683d-02,-1.038d-02,-2.860d-03,
     4-5.478d-02,-5.678d-02, 8.900d-03, 1.484d-02, 5.340d-03, 1.520d-03,
     5 2.220d-02, 2.567d-02,-3.000d-05,-4.970d-03,-2.160d-03,-6.500d-04,
     6-9.530d-03,-1.204d-02,-1.510d-03, 1.510d-03, 8.300d-04, 2.700d-04,
     1 1.261d-01, 1.354d-01, 3.958d-02, 8.240d-03, 1.660d-03, 4.500d-04,
     2 3.890d-03,-1.159d-02,-1.625d-02,-9.610d-03,-3.710d-03,-1.260d-03,
     3-1.910d-03,-5.600d-04, 1.590d-03, 1.590d-03, 8.400d-04, 3.900d-04,
     4 6.400d-04, 4.900d-04,-1.500d-04,-2.900d-04,-1.800d-04,-1.000d-04,
     5-2.000d-04,-1.900d-04, 0.000d+00, 6.000d-05, 4.000d-05, 3.000d-05,
     6 7.000d-05, 8.000d-05, 2.000d-05,-1.000d-05,-1.000d-05,-1.000d-05/
      data (((cehlq(ix,it,nx,2,2),ix=1,6),it=1,6),nx=1,2)/
     1 3.578d-01,-8.622d-02,-1.480d-01,-1.840d-02,-7.820d-03,-4.500d-04,
     2-2.925d-01,-1.304d-01, 1.696d-01, 6.243d-02, 1.353d-02, 3.750d-03,
     3 1.318d-01, 1.041d-01,-5.486d-02,-3.872d-02,-1.038d-02,-2.850d-03,
     4-6.162d-02,-6.143d-02, 1.303d-02, 1.740d-02, 5.940d-03, 1.670d-03,
     5 2.643d-02, 2.957d-02,-1.490d-03,-6.450d-03,-2.630d-03,-7.700d-04,
     6-1.218d-02,-1.497d-02,-1.260d-03, 2.240d-03, 1.120d-03, 3.500d-04,
     1 1.263d-01, 1.334d-01, 3.732d-02, 7.070d-03, 1.260d-03, 3.400d-04,
     2 3.660d-03,-1.357d-02,-1.795d-02,-1.031d-02,-3.880d-03,-1.280d-03,
     3-2.100d-03,-3.600d-04, 2.050d-03, 1.920d-03, 9.800d-04, 4.400d-04,
     4 7.700d-04, 5.400d-04,-2.400d-04,-3.900d-04,-2.400d-04,-1.300d-04,
     5-2.600d-04,-2.300d-04, 2.000d-05, 9.000d-05, 6.000d-05, 4.000d-05,
     6 9.000d-05, 1.000d-04, 2.000d-05,-2.000d-05,-2.000d-05,-1.000d-05/
C...Expansion coefficients for up and down sea quark distributions.
      data (((cehlq(ix,it,nx,3,1),ix=1,6),it=1,6),nx=1,2)/
     1 6.870d-02,-6.861d-02, 2.973d-02,-5.400d-03, 3.780d-03,-9.700d-04,
     2-1.802d-02, 1.400d-04, 6.490d-03,-8.540d-03, 1.220d-03,-1.750d-03,
     3-4.650d-03, 1.480d-03,-5.930d-03, 6.000d-04,-1.030d-03,-8.000d-05,
     4 6.440d-03, 2.570d-03, 2.830d-03, 1.150d-03, 7.100d-04, 3.300d-04,
     5-3.930d-03,-2.540d-03,-1.160d-03,-7.700d-04,-3.600d-04,-1.900d-04,
     6 2.340d-03, 1.930d-03, 5.300d-04, 3.700d-04, 1.600d-04, 9.000d-05,
     1 1.014d+00,-1.106d+00, 3.374d-01,-7.444d-02, 8.850d-03,-8.700d-04,
     2 9.233d-01,-1.285d+00, 4.475d-01,-9.786d-02, 1.419d-02,-1.120d-03,
     3 4.888d-02,-1.271d-01, 8.606d-02,-2.608d-02, 4.780d-03,-6.000d-04,
     4-2.691d-02, 4.887d-02,-1.771d-02, 1.620d-03, 2.500d-04,-6.000d-05,
     5 7.040d-03,-1.113d-02, 1.590d-03, 7.000d-04,-2.000d-04, 0.000d+00,
     6-1.710d-03, 2.290d-03, 3.800d-04,-3.500d-04, 4.000d-05, 1.000d-05/
      data (((cehlq(ix,it,nx,3,2),ix=1,6),it=1,6),nx=1,2)/
     1 1.008d-01,-7.100d-02, 1.973d-02,-5.710d-03, 2.930d-03,-9.900d-04,
     2-5.271d-02,-1.823d-02, 1.792d-02,-6.580d-03, 1.750d-03,-1.550d-03,
     3 1.220d-02, 1.763d-02,-8.690d-03,-8.800d-04,-1.160d-03,-2.100d-04,
     4-1.190d-03,-7.180d-03, 2.360d-03, 1.890d-03, 7.700d-04, 4.100d-04,
     5-9.100d-04, 2.040d-03,-3.100d-04,-1.050d-03,-4.000d-04,-2.400d-04,
     6 1.190d-03,-1.700d-04,-2.000d-04, 4.200d-04, 1.700d-04, 1.000d-04,
     1 1.081d+00,-1.189d+00, 3.868d-01,-8.617d-02, 1.115d-02,-1.180d-03,
     2 9.917d-01,-1.396d+00, 4.998d-01,-1.159d-01, 1.674d-02,-1.720d-03,
     3 5.099d-02,-1.338d-01, 9.173d-02,-2.885d-02, 5.890d-03,-6.500d-04,
     4-3.178d-02, 5.703d-02,-2.070d-02, 2.440d-03, 1.100d-04,-9.000d-05,
     5 8.970d-03,-1.392d-02, 2.050d-03, 6.500d-04,-2.300d-04, 2.000d-05,
     6-2.340d-03, 3.010d-03, 5.000d-04,-3.900d-04, 6.000d-05, 1.000d-05/
C...Expansion coefficients for gluon distribution.
      data (((cehlq(ix,it,nx,4,1),ix=1,6),it=1,6),nx=1,2)/
     1 9.482d-01,-9.578d-01, 1.009d-01,-1.051d-01, 3.456d-02,-3.054d-02,
     2-9.627d-01, 5.379d-01, 3.368d-01,-9.525d-02, 1.488d-02,-2.051d-02,
     3 4.300d-01,-8.306d-02,-3.372d-01, 4.902d-02,-9.160d-03, 1.041d-02,
     4-1.925d-01,-1.790d-02, 2.183d-01, 7.490d-03, 4.140d-03,-1.860d-03,
     5 8.183d-02, 1.926d-02,-1.072d-01,-1.944d-02,-2.770d-03,-5.200d-04,
     6-3.884d-02,-1.234d-02, 5.410d-02, 1.879d-02, 3.350d-03, 1.040d-03,
     1 2.948d+01,-3.902d+01, 1.464d+01,-3.335d+00, 5.054d-01,-5.915d-02,
     2 2.559d+01,-3.955d+01, 1.661d+01,-4.299d+00, 6.904d-01,-8.243d-02,
     3-1.663d+00, 1.176d+00, 1.118d+00,-7.099d-01, 1.948d-01,-2.404d-02,
     4-2.168d-01, 8.170d-01,-7.169d-01, 1.851d-01,-1.924d-02,-3.250d-03,
     5 2.088d-01,-4.355d-01, 2.239d-01,-2.446d-02,-3.620d-03, 1.910d-03,
     6-9.097d-02, 1.601d-01,-5.681d-02,-2.500d-03, 2.580d-03,-4.700d-04/
      data (((cehlq(ix,it,nx,4,2),ix=1,6),it=1,6),nx=1,2)/
     1 2.367d+00, 4.453d-01, 3.660d-01, 9.467d-02, 1.341d-01, 1.661d-02,
     2-3.170d+00,-1.795d+00, 3.313d-02,-2.874d-01,-9.827d-02,-7.119d-02,
     3 1.823d+00, 1.457d+00,-2.465d-01, 3.739d-02, 6.090d-03, 1.814d-02,
     4-1.033d+00,-9.827d-01, 2.136d-01, 1.169d-01, 5.001d-02, 1.684d-02,
     5 5.133d-01, 5.259d-01,-1.173d-01,-1.139d-01,-4.988d-02,-2.021d-02,
     6-2.881d-01,-3.145d-01, 5.667d-02, 9.161d-02, 4.568d-02, 1.951d-02,
     1 3.036d+01,-4.062d+01, 1.578d+01,-3.699d+00, 6.020d-01,-7.031d-02,
     2 2.700d+01,-4.167d+01, 1.770d+01,-4.804d+00, 7.862d-01,-1.060d-01,
     3-1.909d+00, 1.357d+00, 1.127d+00,-7.181d-01, 2.232d-01,-2.481d-02,
     4-2.488d-01, 9.781d-01,-8.127d-01, 2.094d-01,-2.997d-02,-4.710d-03,
     5 2.506d-01,-5.427d-01, 2.672d-01,-3.103d-02,-1.800d-03, 2.870d-03,
     6-1.128d-01, 2.087d-01,-6.972d-02,-2.480d-03, 2.630d-03,-8.400d-04/
C...Expansion coefficients for strange sea quark distribution.
      data (((cehlq(ix,it,nx,5,1),ix=1,6),it=1,6),nx=1,2)/
     1 4.968d-02,-4.173d-02, 2.102d-02,-3.270d-03, 3.240d-03,-6.700d-04,
     2-6.150d-03,-1.294d-02, 6.740d-03,-6.890d-03, 9.000d-04,-1.510d-03,
     3-8.580d-03, 5.050d-03,-4.900d-03,-1.600d-04,-9.400d-04,-1.500d-04,
     4 7.840d-03, 1.510d-03, 2.220d-03, 1.400d-03, 7.000d-04, 3.500d-04,
     5-4.410d-03,-2.220d-03,-8.900d-04,-8.500d-04,-3.600d-04,-2.000d-04,
     6 2.520d-03, 1.840d-03, 4.100d-04, 3.900d-04, 1.600d-04, 9.000d-05,
     1 9.235d-01,-1.085d+00, 3.464d-01,-7.210d-02, 9.140d-03,-9.100d-04,
     2 9.315d-01,-1.274d+00, 4.512d-01,-9.775d-02, 1.380d-02,-1.310d-03,
     3 4.739d-02,-1.296d-01, 8.482d-02,-2.642d-02, 4.760d-03,-5.700d-04,
     4-2.653d-02, 4.953d-02,-1.735d-02, 1.750d-03, 2.800d-04,-6.000d-05,
     5 6.940d-03,-1.132d-02, 1.480d-03, 6.500d-04,-2.100d-04, 0.000d+00,
     6-1.680d-03, 2.340d-03, 4.200d-04,-3.400d-04, 5.000d-05, 1.000d-05/
      data (((cehlq(ix,it,nx,5,2),ix=1,6),it=1,6),nx=1,2)/
     1 6.478d-02,-4.537d-02, 1.643d-02,-3.490d-03, 2.710d-03,-6.700d-04,
     2-2.223d-02,-2.126d-02, 1.247d-02,-6.290d-03, 1.120d-03,-1.440d-03,
     3-1.340d-03, 1.362d-02,-6.130d-03,-7.900d-04,-9.000d-04,-2.000d-04,
     4 5.080d-03,-3.610d-03, 1.700d-03, 1.830d-03, 6.800d-04, 4.000d-04,
     5-3.580d-03, 6.000d-05,-2.600d-04,-1.050d-03,-3.800d-04,-2.300d-04,
     6 2.420d-03, 9.300d-04,-1.000d-04, 4.500d-04, 1.700d-04, 1.100d-04,
     1 9.868d-01,-1.171d+00, 3.940d-01,-8.459d-02, 1.124d-02,-1.250d-03,
     2 1.001d+00,-1.383d+00, 5.044d-01,-1.152d-01, 1.658d-02,-1.830d-03,
     3 4.928d-02,-1.368d-01, 9.021d-02,-2.935d-02, 5.800d-03,-6.600d-04,
     4-3.133d-02, 5.785d-02,-2.023d-02, 2.630d-03, 1.600d-04,-8.000d-05,
     5 8.840d-03,-1.416d-02, 1.900d-03, 5.800d-04,-2.500d-04, 1.000d-05,
     6-2.300d-03, 3.080d-03, 5.500d-04,-3.700d-04, 7.000d-05, 1.000d-05/
C...Expansion coefficients for charm sea quark distribution.
      data (((cehlq(ix,it,nx,6,1),ix=1,6),it=1,6),nx=1,2)/
     1 9.270d-03,-1.817d-02, 9.590d-03,-6.390d-03, 1.690d-03,-1.540d-03,
     2 5.710d-03,-1.188d-02, 6.090d-03,-4.650d-03, 1.240d-03,-1.310d-03,
     3-3.960d-03, 7.100d-03,-3.590d-03, 1.840d-03,-3.900d-04, 3.400d-04,
     4 1.120d-03,-1.960d-03, 1.120d-03,-4.800d-04, 1.000d-04,-4.000d-05,
     5 4.000d-05,-3.000d-05,-1.800d-04, 9.000d-05,-5.000d-05,-2.000d-05,
     6-4.200d-04, 7.300d-04,-1.600d-04, 5.000d-05, 5.000d-05, 5.000d-05,
     1 8.098d-01,-1.042d+00, 3.398d-01,-6.824d-02, 8.760d-03,-9.000d-04,
     2 8.961d-01,-1.217d+00, 4.339d-01,-9.287d-02, 1.304d-02,-1.290d-03,
     3 3.058d-02,-1.040d-01, 7.604d-02,-2.415d-02, 4.600d-03,-5.000d-04,
     4-2.451d-02, 4.432d-02,-1.651d-02, 1.430d-03, 1.200d-04,-1.000d-04,
     5 1.122d-02,-1.457d-02, 2.680d-03, 5.800d-04,-1.200d-04, 3.000d-05,
     6-7.730d-03, 7.330d-03,-7.600d-04,-2.400d-04, 1.000d-05, 0.000d+00/
      data (((cehlq(ix,it,nx,6,2),ix=1,6),it=1,6),nx=1,2)/
     1 9.980d-03,-1.945d-02, 1.055d-02,-6.870d-03, 1.860d-03,-1.560d-03,
     2 5.700d-03,-1.203d-02, 6.250d-03,-4.860d-03, 1.310d-03,-1.370d-03,
     3-4.490d-03, 7.990d-03,-4.170d-03, 2.050d-03,-4.400d-04, 3.300d-04,
     4 1.470d-03,-2.480d-03, 1.460d-03,-5.700d-04, 1.200d-04,-1.000d-05,
     5-9.000d-05, 1.500d-04,-3.200d-04, 1.200d-04,-6.000d-05,-4.000d-05,
     6-4.200d-04, 7.600d-04,-1.400d-04, 4.000d-05, 7.000d-05, 5.000d-05,
     1 8.698d-01,-1.131d+00, 3.836d-01,-8.111d-02, 1.048d-02,-1.300d-03,
     2 9.626d-01,-1.321d+00, 4.854d-01,-1.091d-01, 1.583d-02,-1.700d-03,
     3 3.057d-02,-1.088d-01, 8.022d-02,-2.676d-02, 5.590d-03,-5.600d-04,
     4-2.845d-02, 5.164d-02,-1.918d-02, 2.210d-03,-4.000d-05,-1.500d-04,
     5 1.311d-02,-1.751d-02, 3.310d-03, 5.100d-04,-1.200d-04, 5.000d-05,
     6-8.590d-03, 8.380d-03,-9.200d-04,-2.600d-04, 1.000d-05,-1.000d-05/
C...Expansion coefficients for bottom sea quark distribution.
      data (((cehlq(ix,it,nx,7,1),ix=1,6),it=1,6),nx=1,2)/
     1 9.010d-03,-1.401d-02, 7.150d-03,-4.130d-03, 1.260d-03,-1.040d-03,
     2 6.280d-03,-9.320d-03, 4.780d-03,-2.890d-03, 9.100d-04,-8.200d-04,
     3-2.930d-03, 4.090d-03,-1.890d-03, 7.600d-04,-2.300d-04, 1.400d-04,
     4 3.900d-04,-1.200d-03, 4.400d-04,-2.500d-04, 2.000d-05,-2.000d-05,
     5 2.600d-04, 1.400d-04,-8.000d-05, 1.000d-04, 1.000d-05, 1.000d-05,
     6-2.600d-04, 3.200d-04, 1.000d-05,-1.000d-05, 1.000d-05,-1.000d-05,
     1 8.029d-01,-1.075d+00, 3.792d-01,-7.843d-02, 1.007d-02,-1.090d-03,
     2 7.903d-01,-1.099d+00, 4.153d-01,-9.301d-02, 1.317d-02,-1.410d-03,
     3-1.704d-02,-1.130d-02, 2.882d-02,-1.341d-02, 3.040d-03,-3.600d-04,
     4-7.200d-04, 7.230d-03,-5.160d-03, 1.080d-03,-5.000d-05,-4.000d-05,
     5 3.050d-03,-4.610d-03, 1.660d-03,-1.300d-04,-1.000d-05, 1.000d-05,
     6-4.360d-03, 5.230d-03,-1.610d-03, 2.000d-04,-2.000d-05, 0.000d+00/
      data (((cehlq(ix,it,nx,7,2),ix=1,6),it=1,6),nx=1,2)/
     1 8.980d-03,-1.459d-02, 7.510d-03,-4.410d-03, 1.310d-03,-1.070d-03,
     2 5.970d-03,-9.440d-03, 4.800d-03,-3.020d-03, 9.100d-04,-8.500d-04,
     3-3.050d-03, 4.440d-03,-2.100d-03, 8.500d-04,-2.400d-04, 1.400d-04,
     4 5.300d-04,-1.300d-03, 5.600d-04,-2.700d-04, 3.000d-05,-2.000d-05,
     5 2.000d-04, 1.400d-04,-1.100d-04, 1.000d-04, 0.000d+00, 0.000d+00,
     6-2.600d-04, 3.200d-04, 0.000d+00,-3.000d-05, 1.000d-05,-1.000d-05,
     1 8.672d-01,-1.174d+00, 4.265d-01,-9.252d-02, 1.244d-02,-1.460d-03,
     2 8.500d-01,-1.194d+00, 4.630d-01,-1.083d-01, 1.614d-02,-1.830d-03,
     3-2.241d-02,-5.630d-03, 2.815d-02,-1.425d-02, 3.520d-03,-4.300d-04,
     4-7.300d-04, 8.030d-03,-5.780d-03, 1.380d-03,-1.300d-04,-4.000d-05,
     5 3.460d-03,-5.380d-03, 1.960d-03,-2.100d-04, 1.000d-05, 1.000d-05,
     6-4.850d-03, 5.950d-03,-1.890d-03, 2.600d-04,-3.000d-05, 0.000d+00/
C...Expansion coefficients for top sea quark distribution.
      data (((cehlq(ix,it,nx,8,1),ix=1,6),it=1,6),nx=1,2)/
     1 4.410d-03,-7.480d-03, 3.770d-03,-2.580d-03, 7.300d-04,-7.100d-04,
     2 3.840d-03,-6.050d-03, 3.030d-03,-2.030d-03, 5.800d-04,-5.900d-04,
     3-8.800d-04, 1.660d-03,-7.500d-04, 4.700d-04,-1.000d-04, 1.000d-04,
     4-8.000d-05,-1.500d-04, 1.200d-04,-9.000d-05, 3.000d-05, 0.000d+00,
     5 1.300d-04,-2.200d-04,-2.000d-05,-2.000d-05,-2.000d-05,-2.000d-05,
     6-7.000d-05, 1.900d-04,-4.000d-05, 2.000d-05, 0.000d+00, 0.000d+00,
     1 6.623d-01,-9.248d-01, 3.519d-01,-7.930d-02, 1.110d-02,-1.180d-03,
     2 6.380d-01,-9.062d-01, 3.582d-01,-8.479d-02, 1.265d-02,-1.390d-03,
     3-2.581d-02, 2.125d-02, 4.190d-03,-4.980d-03, 1.490d-03,-2.100d-04,
     4 7.100d-04, 5.300d-04,-1.270d-03, 3.900d-04,-5.000d-05,-1.000d-05,
     5 3.850d-03,-5.060d-03, 1.860d-03,-3.500d-04, 4.000d-05, 0.000d+00,
     6-3.530d-03, 4.460d-03,-1.500d-03, 2.700d-04,-3.000d-05, 0.000d+00/
      data (((cehlq(ix,it,nx,8,2),ix=1,6),it=1,6),nx=1,2)/
     1 4.260d-03,-7.530d-03, 3.830d-03,-2.680d-03, 7.600d-04,-7.300d-04,
     2 3.640d-03,-6.050d-03, 3.030d-03,-2.090d-03, 5.900d-04,-6.000d-04,
     3-9.200d-04, 1.710d-03,-8.200d-04, 5.000d-04,-1.200d-04, 1.000d-04,
     4-5.000d-05,-1.600d-04, 1.300d-04,-9.000d-05, 3.000d-05, 0.000d+00,
     5 1.300d-04,-2.100d-04,-1.000d-05,-2.000d-05,-2.000d-05,-1.000d-05,
     6-8.000d-05, 1.800d-04,-5.000d-05, 2.000d-05, 0.000d+00, 0.000d+00,
     1 7.146d-01,-1.007d+00, 3.932d-01,-9.246d-02, 1.366d-02,-1.540d-03,
     2 6.856d-01,-9.828d-01, 3.977d-01,-9.795d-02, 1.540d-02,-1.790d-03,
     3-3.053d-02, 2.758d-02, 2.150d-03,-4.880d-03, 1.640d-03,-2.500d-04,
     4 9.200d-04, 4.200d-04,-1.340d-03, 4.600d-04,-8.000d-05,-1.000d-05,
     5 4.230d-03,-5.660d-03, 2.140d-03,-4.300d-04, 6.000d-05, 0.000d+00,
     6-3.890d-03, 5.000d-03,-1.740d-03, 3.300d-04,-4.000d-05, 0.000d+00/

C...Determine set, Lambda and x and t expansion variables.
        nset=mstp(51)
        if(nset.eq.1) alam=0.2d0
        if(nset.eq.2) alam=0.29d0
        vint(231)=5.d0
        tmin=log(5.d0/alam**2)
        tmax=log(1d8/alam**2)
        if(mstp(57).eq.0) then
          t=tmin
        else
          t=log(max(1.d0,q2/alam**2))
        endif
        vt=max(-1.d0,min(1.d0,(2.d0*t-tmax-tmin)/(tmax-tmin)))
        nx=1
        if(x.le.0.1d0) nx=2
        if(nx.eq.1) vx=(2.d0*x-1.1d0)/0.9d0
        if(nx.eq.2) vx=max(-1.d0,(2.d0*log(x)+11.51293d0)/6.90776d0)
        cxs=1.d0
        if(x.lt.1d-4.and.abs(parp(51)-1.d0).gt.0.01d0) cxs=
     &  (1d-4/x)**(parp(51)-1.d0)
 
C...Chebyshev polynomials for x and t expansion.
        tx(1)=1.d0
        tx(2)=vx
        tx(3)=2.d0*vx**2-1.d0
        tx(4)=4.d0*vx**3-3.d0*vx
        tx(5)=8.d0*vx**4-8.d0*vx**2+1.d0
        tx(6)=16.d0*vx**5-20.d0*vx**3+5.d0*vx
        tt(1)=1.d0
        tt(2)=vt
        tt(3)=2.d0*vt**2-1.d0
        tt(4)=4.d0*vt**3-3.d0*vt
        tt(5)=8.d0*vt**4-8.d0*vt**2+1.d0
        tt(6)=16.d0*vt**5-20.d0*vt**3+5.d0*vt
 
C...Calculate structure functions.
        do 130 kfl=1,6
        xqsum=0.d0
        do 120 it=1,6
        do 110 ix=1,6
        xqsum=xqsum+cehlq(ix,it,nx,kfl,nset)*tx(ix)*tt(it)
  110   continue
  120   continue
       xq(kfl)=xqsum*(1.d0-x)**nehlq(kfl,nset)*cxs
  130   continue
 
C...Put into output array.
        xppr(0)=xq(4)
        xppr(1)=xq(2)+xq(3)
        xppr(2)=xq(1)+xq(3)
        xppr(3)=xq(5)
        xppr(4)=xq(6)
        xppr(-1)=xq(3)
        xppr(-2)=xq(3)
        xppr(-3)=xq(5)
        xppr(-4)=xq(6)
 
C...Special expansion for bottom (threshold effects).
        if(mstp(58).ge.5) then
          if(nset.eq.1) tmin=8.1905d0
          if(nset.eq.2) tmin=7.4474d0
          if(t.gt.tmin) then
            vt=max(-1.d0,min(1.d0,(2.d0*t-tmax-tmin)/(tmax-tmin)))
            tt(1)=1.d0
            tt(2)=vt
            tt(3)=2.d0*vt**2-1.d0
            tt(4)=4.d0*vt**3-3.d0*vt
            tt(5)=8.d0*vt**4-8.d0*vt**2+1.d0
            tt(6)=16.d0*vt**5-20.d0*vt**3+5.d0*vt
            xqsum=0.d0
            do 150 it=1,6
            do 140 ix=1,6
            xqsum=xqsum+cehlq(ix,it,nx,7,nset)*tx(ix)*tt(it)
  140       continue
  150       continue
            xppr(5)=xqsum*(1.d0-x)**nehlq(7,nset)*cxs
            xppr(-5)=xppr(5)
          endif
        endif
 
C...Special expansion for top (threshold effects).
        if(mstp(58).ge.6) then
          if(nset.eq.1) tmin=11.5528d0
          if(nset.eq.2) tmin=10.8097d0
          tmin=tmin+2.d0*log(pmas(6,1)/30.d0)
          tmax=tmax+2.d0*log(pmas(6,1)/30.d0)
          if(t.gt.tmin) then
            vt=max(-1.d0,min(1.d0,(2.d0*t-tmax-tmin)/(tmax-tmin)))
            tt(1)=1.d0
            tt(2)=vt
            tt(3)=2.d0*vt**2-1.d0
            tt(4)=4.d0*vt**3-3.d0*vt
            tt(5)=8.d0*vt**4-8.d0*vt**2+1.d0
            tt(6)=16.d0*vt**5-20.d0*vt**3+5.d0*vt
            xqsum=0.d0
            do 170 it=1,6
            do 160 ix=1,6
            xqsum=xqsum+cehlq(ix,it,nx,8,nset)*tx(ix)*tt(it)
  160       continue
  170       continue
            xppr(6)=xqsum*(1.d0-x)**nehlq(8,nset)*cxs
            xppr(-6)=xppr(6)
          endif
        endif

      end

c***********************************************************************

      subroutine dukeowen(x,q2,xppr)

C...Proton structure functions from Duke, Owens.
C...Allowed variable range: 4 GeV^2 < Q^2 < approx 1E6 GeV^2.

      implicit double precision(a-h, o-z)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      save /pjpars/,/pjint1/

      dimension xppr(-6:6),xq(9),ts(6),cdo(3,6,5,2)

C...The following data lines are coefficients needed in the
C...Duke, Owens proton structure function parametrizations, see below.
C...Expansion coefficients for (up+down) valence quark distribution.
      data ((cdo(ip,is,1,1),is=1,6),ip=1,3)/
     1 4.190d-01, 3.460d+00, 4.400d+00, 0.000d+00, 0.000d+00, 0.000d+00,
     2 4.000d-03, 7.240d-01,-4.860d+00, 0.000d+00, 0.000d+00, 0.000d+00,
     3-7.000d-03,-6.600d-02, 1.330d+00, 0.000d+00, 0.000d+00, 0.000d+00/
      data ((cdo(ip,is,1,2),is=1,6),ip=1,3)/
     1 3.740d-01, 3.330d+00, 6.030d+00, 0.000d+00, 0.000d+00, 0.000d+00,
     2 1.400d-02, 7.530d-01,-6.220d+00, 0.000d+00, 0.000d+00, 0.000d+00,
     3 0.000d+00,-7.600d-02, 1.560d+00, 0.000d+00, 0.000d+00, 0.000d+00/
C...Expansion coefficients for down valence quark distribution.
      data ((cdo(ip,is,2,1),is=1,6),ip=1,3)/
     1 7.630d-01, 4.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00,
     2-2.370d-01, 6.270d-01,-4.210d-01, 0.000d+00, 0.000d+00, 0.000d+00,
     3 2.600d-02,-1.900d-02, 3.300d-02, 0.000d+00, 0.000d+00, 0.000d+00/
      data ((cdo(ip,is,2,2),is=1,6),ip=1,3)/
     1 7.610d-01, 3.830d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00,
     2-2.320d-01, 6.270d-01,-4.180d-01, 0.000d+00, 0.000d+00, 0.000d+00,
     3 2.300d-02,-1.900d-02, 3.600d-02, 0.000d+00, 0.000d+00, 0.000d+00/
C...Expansion coefficients for (up+down+strange) sea quark distribution.
      data ((cdo(ip,is,3,1),is=1,6),ip=1,3)/
     1 1.265d+00, 0.000d+00, 8.050d+00, 0.000d+00, 0.000d+00, 0.000d+00,
     2-1.132d+00,-3.720d-01, 1.590d+00, 6.310d+00,-1.050d+01, 1.470d+01,
     3 2.930d-01,-2.900d-02,-1.530d-01,-2.730d-01,-3.170d+00, 9.800d+00/
      data ((cdo(ip,is,3,2),is=1,6),ip=1,3)/
     1 1.670d+00, 0.000d+00, 9.150d+00, 0.000d+00, 0.000d+00, 0.000d+00,
     2-1.920d+00,-2.730d-01, 5.300d-01, 1.570d+01,-1.010d+02, 2.230d+02,
     3 5.820d-01,-1.640d-01,-7.630d-01,-2.830d+00, 4.470d+01,-1.170d+02/
C...Expansion coefficients for charm sea quark distribution.
      data ((cdo(ip,is,4,1),is=1,6),ip=1,3)/
     1 0.000d+00,-3.600d-02, 6.350d+00, 0.000d+00, 0.000d+00, 0.000d+00,
     2 1.350d-01,-2.220d-01, 3.260d+00,-3.030d+00, 1.740d+01,-1.790d+01,
     3-7.500d-02,-5.800d-02,-9.090d-01, 1.500d+00,-1.130d+01, 1.560d+01/
       data ((cdo(ip,is,4,2),is=1,6),ip=1,3)/
     1 0.000d+00,-1.200d-01, 3.510d+00, 0.000d+00, 0.000d+00, 0.000d+00,
     2 6.700d-02,-2.330d-01, 3.660d+00,-4.740d-01, 9.500d+00,-1.660d+01,
     3-3.100d-02,-2.300d-02,-4.530d-01, 3.580d-01,-5.430d+00, 1.550d+01/
C...Expansion coefficients for gluon distribution.
      data ((cdo(ip,is,5,1),is=1,6),ip=1,3)/
     1 1.560d+00, 0.000d+00, 6.000d+00, 9.000d+00, 0.000d+00, 0.000d+00,
     2-1.710d+00,-9.490d-01, 1.440d+00,-7.190d+00,-1.650d+01, 1.530d+01,
     3 6.380d-01, 3.250d-01,-1.050d+00, 2.550d-01, 1.090d+01,-1.010d+01/
      data ((cdo(ip,is,5,2),is=1,6),ip=1,3)/
     1 8.790d-01, 0.000d+00, 4.000d+00, 9.000d+00, 0.000d+00, 0.000d+00,
     2-9.710d-01,-1.160d+00, 1.230d+00,-5.640d+00,-7.540d+00,-5.960d-01,
     3 4.340d-01, 4.760d-01,-2.540d-01,-8.170d-01, 5.500d+00, 1.260d-01/
 
C...Euler's beta function, requires ordinary Gamma function
      eulbet(x,y)=pjgamm(x)*pjgamm(y)/pjgamm(x+y)

C...Determine set, Lambda and s expansion parameter.
        nset=mstp(51)-2
        if(nset.eq.1) alam=0.2d0
        if(nset.eq.2) alam=0.4d0
        vint(231)=4.d0
        if(mstp(57).le.0) then
          sd=0.d0
        else
          q2in=min(1d6,max(4.d0,q2))
          sd=log(log(q2in/alam**2)/log(4.d0/alam**2))
        endif
 
C...Calculate structure functions.
        do 190 kfl=1,5
        do 180 is=1,6
        ts(is)=cdo(1,is,kfl,nset)+cdo(2,is,kfl,nset)*sd+
     &  cdo(3,is,kfl,nset)*sd**2
  180   continue
        if(kfl.le.2) then
          xq(kfl)=x**ts(1)*(1.d0-x)**ts(2)*(1.d0+ts(3)*x)/(eulbet(ts(1),
     &    ts(2)+1.d0)*(1.d0+ts(3)*ts(1)/(ts(1)+ts(2)+1.d0)))
        else
          xq(kfl)=ts(1)*x**ts(2)*(1.d0-x)**ts(3)*(1.d0+ts(4)*x 
     & +ts(5)*x**2+
     &    ts(6)*x**3)
        endif
  190   continue
 
C...Put into output arrays.
        xppr(0)=xq(5)
        xppr(1)=xq(2)+xq(3)/6.d0
        xppr(2)=3.d0*xq(1)-xq(2)+xq(3)/6.d0
        xppr(3)=xq(3)/6.d0
        xppr(4)=xq(4)
        xppr(-1)=xq(3)/6.d0
        xppr(-2)=xq(3)/6.d0
        xppr(-3)=xq(3)/6.d0
        xppr(-4)=xq(4)

c     write(70,*)'pjstpr'
c     write(70,'(7(f9.6,1x))')(xppr(i)-xppr(-i),xppr(-i),i=1,3),xppr(0)

      end

C*********************************************************************
 
C...PDFSET
C...Dummy routine, to be removed when PDFLIB is to be linked.
 
      subroutine pjpdfset(parm,value)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      save /jydat1/
C...Local arrays and character variables.
      character*20 parm(20)
      double precision value(20)
 
C...Stop program if this routine is ever called.
      write(mstu(11),5000)
      if(pjr(0).lt.10d0) stop
      parm(20)=parm(1)
      value(20)=value(1)
 
C...Format for error printout.
 5000 format(1x,'Error: you did not link PDFLIB correctly.'/
     &1x,'Dummy routine PDFSET in PYTHIA file called instead.'/
     &1x,'Execution stopped!')
 
      return
      end
 
C*********************************************************************
 
C...STRUCTM
C...Dummy routine, to be removed when PDFLIB is to be linked.
 
      subroutine pjstructm(xx,qq,upv,dnv,usea,dsea,str,chm,bot,top,glu)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      save /jydat1/
C...Local variables
      double precision xx,qq,upv,dnv,usea,dsea,str,chm,bot,top,glu
 
C...Stop program if this routine is ever called.
      write(mstu(11),5000)
      if(pjr(0).lt.10d0) stop
      upv=xx+qq
      dnv=xx+2d0*qq
      usea=xx+3d0*qq
      dsea=xx+4d0*qq
      str=xx+5d0*qq
      chm=xx+6d0*qq
      bot=xx+7d0*qq
      top=xx+8d0*qq
      glu=xx+9d0*qq
 
C...Format for error printout.
 5000 format(1x,'Error: you did not link PDFLIB correctly.'/
     &1x,'Dummy routine STRUCTM in PYTHIA file called instead.'/
     &1x,'Execution stopped!')
 
      return
      end
