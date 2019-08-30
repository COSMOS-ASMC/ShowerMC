C***********************************************************************
 
C...PYSIGH
C...Differential matrix elements for all included subprocesses
C...Note that what is coded is (disregarding the COMFAC factor)
C...1) for 2 -> 1 processes: s-hat/pi*d(sigma-hat), where,
C...when d(sigma-hat) is given in the zero-width limit, the delta
C...function in tau is replaced by a (modified) Breit-Wigner:
C...1/pi*s*H_res/((s*tau-m_res^2)^2+H_res^2),
C...where H_res = s-hat/m_res*Gamma_res(s-hat);
C...2) for 2 -> 2 processes: (s-hat)**2/pi*d(sigma-hat)/d(t-hat);
C...i.e., dimensionless quantities
C...3) for 2 -> 3 processes: abs(M)^2, where the total cross-section is
C...Integral abs(M)^2/(2shat') * (prod_(i=1)^3 d^3p_i/((2pi)^3*2E_i)) *
C...(2pi)^4 delta^4(P - sum p_i)
C...COMFAC contains the factor pi/s (or equivalent) and
C...the conversion factor from GeV^-2 to mb
 
      subroutine pjsigh(nchn,sigs)
 
C...Double precision declarations
      implicit double precision(a-h, o-z)
C...Parameter statement to help give large particle numbers.
      parameter (ksusy1=1000000,ksusy2=2000000,kexcit=4000000)
C...Commonblocks
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/jydat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      common/pjint3/xsfx(2,-40:40),isig(1000,3),sigh(1000)
      common/pjint4/mwid(500),wids(500,5)
      common/pjint5/ngenpd,ngen(0:500,3),xsec(0:500,3)
      common/pjint7/sigt(0:6,0:6,0:5)
      common/pjssmt/zmix(4,4),umix(2,2),vmix(2,2),smz(4),smw(2),
     &sfmix(16,4)
      save /jyjets/,/jydat1/,/jydat2/,/jydat3/,/pjsubs/,/pjpars/,
     &/pjint1/,/pjint2/,/pjint3/,/pjint4/,/pjint5/,/pjint7/,
     &/pjssmt/
C...Local arrays and complex variables
      dimension x(2),xpq(-25:25),kfac(2,-40:40),wdtp(0:200),
     &wdte(0:200,0:5),hgz(6,3),hl3(3),hr3(3),hl4(3),hr4(3)
      complex a004,a204,a114,a00u,a20u,a11u
      complex cigtot,ciztot,f0alp,f1alp,f2alp,f0bet,f1bet,f2bet,fif,
     &coulck,coulcp,coulcd,coulcr,coulcs
      real a00l,a11l,a20l,coulxx
 
C...Reset number of channels and cross-section
      nchn=0
      sigs=0d0
 
C...Convert H or A process into equivalent h one
      isub=mint(1)
      isubsv=isub
      ihigg=1
      kfhigg=25
      if((isub.ge.151.and.isub.le.160).or.(isub.ge.171.and.
     &isub.le.190)) then
        ihigg=2
        if(mod(isub-1,10).ge.5) ihigg=3
        kfhigg=33+ihigg
        if(isub.eq.151.or.isub.eq.156) isub=3
        if(isub.eq.152.or.isub.eq.157) isub=102
        if(isub.eq.153.or.isub.eq.158) isub=103
        if(isub.eq.171.or.isub.eq.176) isub=24
        if(isub.eq.172.or.isub.eq.177) isub=26
        if(isub.eq.173.or.isub.eq.178) isub=123
        if(isub.eq.174.or.isub.eq.179) isub=124
        if(isub.eq.181.or.isub.eq.186) isub=121
        if(isub.eq.182.or.isub.eq.187) isub=122
      endif
 
CMRENNA++
C...Convert almost equivalent SUSY processes into each other
C...Extract differences in flavours and couplings
      if(isub.ge.200.and.isub.le.280) then
        call pjerrm(30,'(pjsigh:) invalid isub') 
      endif
CMRENNA--
 
C...Read kinematical variables and limits
      istsb=iset(isubsv)
      taumin=vint(11)
      ystmin=vint(12)
      ctnmin=vint(13)
      ctpmin=vint(14)
      taupmn=vint(16)
      tau=vint(21)
      yst=vint(22)
      cth=vint(23)
      xt2=vint(25)
      taup=vint(26)
      taumax=vint(31)
      ystmax=vint(32)
      ctnmax=vint(33)
      ctpmax=vint(34)
      taupmx=vint(36)
 
C...Derive kinematical quantities
      taue=tau
      if(istsb.ge.3.and.istsb.le.5) taue=taup
      x(1)=sqrt(taue)*exp(yst)
      x(2)=sqrt(taue)*exp(-yst)
      if(mint(45).eq.2.and.istsb.ge.1) then
        if(x(1).gt.0.9999d0) return
      elseif(mint(45).eq.3) then
        x(1)=min(0.9999989d0,x(1))
      endif
      if(mint(46).eq.2.and.istsb.ge.1) then
        if(x(2).gt.0.9999d0) return
      elseif(mint(46).eq.3) then
        x(2)=min(0.9999989d0,x(2))
      endif
      sh=tau*vint(2)
      sqm3=vint(63)
      sqm4=vint(64)
      rm3=sqm3/sh
      rm4=sqm4/sh
      be34=sqrt(max(0d0,(1d0-rm3-rm4)**2-4d0*rm3*rm4))
      rpts=4d0*vint(71)**2/sh
      be34l=sqrt(max(0d0,(1d0-rm3-rm4)**2-4d0*rm3*rm4-rpts))
      rm34=max(1d-20,2d0*rm3*rm4)
      rsqm=1d0+rm34
      if(2d0*vint(71)**2/(vint(21)*vint(2)).lt.0.0001d0) rm34=max(rm34,
     &2d0*vint(71)**2/(vint(21)*vint(2)))
      rthm=(4d0*rm3*rm4+rpts)/(1d0-rm3-rm4+be34l)
      if(istsb.eq.0) then
        th=vint(45)
        uh=-0.5d0*sh*max(rthm,1d0-rm3-rm4+be34*cth)
        sqpth=max(vint(71)**2,0.25d0*sh*be34**2*vint(59)**2)
      else
        th=-0.5d0*sh*max(rthm,1d0-rm3-rm4-be34*cth)
        uh=-0.5d0*sh*max(rthm,1d0-rm3-rm4+be34*cth)
        sqpth=max(vint(71)**2,0.25d0*sh*be34**2*(1d0-cth**2))
      endif
      shr=sqrt(sh)
      sh2=sh**2
      th2=th**2
      uh2=uh**2
 
C...Choice of Q2 scale: hard, parton distributions, parton showers
      if(istsb.eq.1.or.istsb.eq.3.or.istsb.eq.5) then
        q2=sh
      elseif(mod(istsb,2).eq.0.or.istsb.eq.9) then
        if(mstp(32).eq.1) then
          q2=2d0*sh*th*uh/(sh**2+th**2+uh**2)
        elseif(mstp(32).eq.2) then
          q2=sqpth+0.5d0*(sqm3+sqm4)
        elseif(mstp(32).eq.3) then
          q2=min(-th,-uh)
        elseif(mstp(32).eq.4) then
          q2=sh
        elseif(mstp(32).eq.5) then
          q2=-th
        endif
        if(istsb.eq.9) q2=sqpth
        if((istsb.eq.9.and.mstp(82).ge.2).or.(istsb.ne.9.and.
     &  mstp(85).eq.1)) q2=q2+parp(82)**2
      endif
      q2sf=q2
      if(istsb.ge.3.and.istsb.le.5) then
        q2sf=pmas(23,1)**2
        if(isub.eq.8.or.isub.eq.76.or.isub.eq.77.or.isub.eq.124)
     &  q2sf=pmas(24,1)**2
        if(isub.eq.121.or.isub.eq.122) then
          q2sf=pmas(jamcomp(kfpr(isubsv,2)),1)**2
          if(mstp(39).eq.2) q2sf=q2sf+max(vint(202),vint(207))
          if(mstp(39).eq.3) q2sf=sh
          if(mstp(39).eq.4) q2sf=vint(26)*vint(2)
        endif
      endif
      q2ps=q2sf
      q2sf=q2sf*parp(34)
      if(mstp(68).ge.2.and.mint(47).eq.5) q2sf=vint(2)
      if(mstp(22).ge.1.and.(isub.eq.10.or.isub.eq.83).and.
     &(mint(43).eq.2.or.mint(43).eq.3)) then
        xbj=x(2)
        if(mint(43).eq.3) xbj=x(1)
        if(mstp(22).eq.1) then
          q2ps=-th
        elseif(mstp(22).eq.2) then
          q2ps=((1d0-xbj)/xbj)*(-th)
        elseif(mstp(22).eq.3) then
          q2ps=sqrt((1d0-xbj)/xbj)*(-th)
        else
          q2ps=(1d0-xbj)*max(1d0,-log(xbj))*(-th)
        endif
      endif
      if(mstp(68).ge.1.and.mint(47).eq.5) q2ps=vint(2)
 
C...Store derived kinematical quantities
      vint(41)=x(1)
      vint(42)=x(2)
      vint(44)=sh
      vint(43)=sqrt(sh)
      vint(45)=th
      vint(46)=uh
      vint(48)=sqpth
      vint(47)=sqrt(sqpth)
      vint(50)=taup*vint(2)
      vint(49)=sqrt(max(0d0,vint(50)))
      vint(52)=q2
      vint(51)=sqrt(q2)
      vint(54)=q2sf
      vint(53)=sqrt(q2sf)
      vint(56)=q2ps
      vint(55)=sqrt(q2ps)
 
C...Calculate parton distributions
      if(istsb.le.0) goto 170
      if(mint(47).ge.2) then
        do 110 i=3-min(2,mint(45)),min(2,mint(46))
          xsf=x(i)
          if(istsb.eq.9) xsf=x(i)/vint(142+i)
          mint(105)=mint(102+i)
          mint(109)=mint(106+i)
          if(mstp(57).le.1) then
            call pjpdfu(mint(10+i),xsf,q2sf,xpq)
          else
            call pjpdfl(mint(10+i),xsf,q2sf,xpq)
          endif
          do 100 kfl=-25,25
            xsfx(i,kfl)=xpq(kfl)
  100     continue
  110   continue
      endif
 
C...Calculate alpha_em, alpha_strong and K-factor
      xw=paru(102)
      xwv=xw
      if(mstp(8).ge.2.or.(isub.ge.71.and.isub.le.77)) xw=
     &1d0-(pmas(24,1)/pmas(23,1))**2
      xw1=1d0-xw
      xwc=1d0/(16d0*xw*xw1)
      aem=pjalem(q2)
      if(mstp(8).ge.1) aem=sqrt(2d0)*paru(105)*pmas(24,1)**2*xw/paru(1)
      if(mstp(33).ne.3) as=pjalps(parp(34)*q2)
      fack=1d0
      faca=1d0
      if(mstp(33).eq.1) then
        fack=parp(31)
      elseif(mstp(33).eq.2) then
        fack=parp(31)
        faca=parp(32)/parp(31)
      elseif(mstp(33).eq.3) then
        q2as=parp(33)*q2
        if(istsb.eq.9.and.mstp(82).ge.2) q2as=q2as+
     &  paru(112)*parp(82)
        as=pjalps(q2as)
      endif
      vint(138)=1d0
      vint(57)=aem
      vint(58)=as
 
C...Set flags for allowed reacting partons/leptons
      do 140 i=1,2
        do 120 j=-25,25
          kfac(i,j)=0
  120   continue
        if(mint(44+i).eq.1) then
          kfac(i,mint(10+i))=1
        elseif(mint(40+i).eq.1.and.mstp(12).eq.0) then
          kfac(i,mint(10+i))=1
          kfac(i,22)=1
          kfac(i,24)=1
          kfac(i,-24)=1
        else
          do 130 j=-25,25
            kfac(i,j)=kfin(i,j)
            if(iabs(j).gt.mstp(58).and.iabs(j).le.10) kfac(i,j)=0
            if(xsfx(i,j).lt.1d-10) kfac(i,j)=0
  130     continue
        endif
  140 continue
 
C...Lower and upper limit for fermion flavour loops
      mmin1=0
      mmax1=0
      mmin2=0
      mmax2=0
      do 150 j=-20,20
        if(kfac(1,-j).eq.1) mmin1=-j
        if(kfac(1,j).eq.1) mmax1=j
        if(kfac(2,-j).eq.1) mmin2=-j
        if(kfac(2,j).eq.1) mmax2=j
  150 continue
      mmina=min(mmin1,mmin2)
      mmaxa=max(mmax1,mmax2)
 
C...Common resonance mass and width combinations
      sqmz=pmas(23,1)**2
      sqmw=pmas(24,1)**2
      sqmh=pmas(kfhigg,1)**2
      gmmz=pmas(23,1)*pmas(23,2)
      gmmw=pmas(24,1)*pmas(24,2)
      gmmh=pmas(kfhigg,1)*pmas(kfhigg,2)
C...MRENNA+++
      zwid=pmas(23,2)
      wwid=pmas(24,2)
      tanw=sqrt(xw/xw1)
C...MRENNA---
 
C...Phase space integral in tau
      comfac=paru(1)*paru(5)/vint(2)
      if(mint(41).eq.2.and.mint(42).eq.2) comfac=comfac*fack
      if((mint(47).ge.2.or.(istsb.ge.3.and.istsb.le.5)).and.
     &istsb.ne.9) then
        atau1=log(taumax/taumin)
        atau2=(taumax-taumin)/(taumax*taumin)
        h1=coef(isubsv,1)+(atau1/atau2)*coef(isubsv,2)/tau
        if(mint(72).ge.1) then
          taur1=vint(73)
          gamr1=vint(74)
          ataud=log(taumax/taumin*(taumin+taur1)/(taumax+taur1))
          atau3=ataud/taur1
          if(ataud.gt.1d-6) h1=h1+
     &    (atau1/atau3)*coef(isubsv,3)/(tau+taur1)
          ataud=atan((taumax-taur1)/gamr1)-atan((taumin-taur1)/gamr1)
          atau4=ataud/gamr1
          if(ataud.gt.1d-6) h1=h1+
     &    (atau1/atau4)*coef(isubsv,4)*tau/((tau-taur1)**2+gamr1**2)
        endif
        if(mint(72).eq.2) then
          taur2=vint(75)
          gamr2=vint(76)
          ataud=log(taumax/taumin*(taumin+taur2)/(taumax+taur2))
          atau5=ataud/taur2
          if(ataud.gt.1d-6) h1=h1+
     &    (atau1/atau5)*coef(isubsv,5)/(tau+taur2)
          ataud=atan((taumax-taur2)/gamr2)-atan((taumin-taur2)/gamr2)
          atau6=ataud/gamr2
          if(ataud.gt.1d-6) h1=h1+
     &    (atau1/atau6)*coef(isubsv,6)*tau/((tau-taur2)**2+gamr2**2)
        endif
        if(mint(47).eq.5.and.(istsb.le.2.or.istsb.ge.5)) then
          atau7=log(max(2d-6,1d0-taumin)/max(2d-6,1d0-taumax))
          if(atau7.gt.1d-6) h1=h1+(atau1/atau7)*coef(isubsv,7)*tau/
     &    max(2d-6,1d0-tau)
        endif
        comfac=comfac*atau1/(tau*h1)
      endif
 
C...Phase space integral in y*
      if(mint(47).ge.4.and.istsb.ne.9) then
        ayst0=ystmax-ystmin
        if(ayst0.lt.1d-6) then
          comfac=0d0
        else
          ayst1=0.5d0*(ystmax-ystmin)**2
          ayst2=ayst1
          ayst3=2d0*(atan(exp(ystmax))-atan(exp(ystmin)))
          h2=(ayst0/ayst1)*coef(isubsv,8)*(yst-ystmin)+
     &    (ayst0/ayst2)*coef(isubsv,9)*(ystmax-yst)+
     &    (ayst0/ayst3)*coef(isubsv,10)/cosh(yst)
          if(mint(45).eq.3) then
            yst0=-0.5d0*log(taue)
            ayst4=log(max(1d-6,exp(yst0-ystmin)-1d0)/
     &      max(1d-6,exp(yst0-ystmax)-1d0))
            if(ayst4.gt.1d-6) h2=h2+(ayst0/ayst4)*coef(isubsv,11)/
     &      max(1d-6,1d0-exp(yst-yst0))
          endif
          if(mint(46).eq.3) then
            yst0=-0.5d0*log(taue)
            ayst5=log(max(1d-6,exp(yst0+ystmax)-1d0)/
     &      max(1d-6,exp(yst0+ystmin)-1d0))
            if(ayst5.gt.1d-6) h2=h2+(ayst0/ayst5)*coef(isubsv,12)/
     &      max(1d-6,1d0-exp(-yst-yst0))
          endif
          comfac=comfac*ayst0/h2
        endif
      endif
 
C...2 -> 1 processes: reduction in angular part of phase space integral
C...for case of decaying resonance
      acth0=ctnmax-ctnmin+ctpmax-ctpmin
      if((istsb.eq.1.or.istsb.eq.3.or.istsb.eq.5)) then
        if(mdcy(jamcomp(kfpr(isubsv,1)),1).eq.1) then
          if(kfpr(isub,1).eq.25.or.kfpr(isub,1).eq.37.or.
     &    kfpr(isub,1).eq.39) then
            comfac=comfac*0.5d0*acth0
          else
            comfac=comfac*0.125d0*(3d0*acth0+ctnmax**3-ctnmin**3+
     &      ctpmax**3-ctpmin**3)
          endif
        endif
 
C...2 -> 2 processes: angular part of phase space integral
      elseif(istsb.eq.2.or.istsb.eq.4) then
        acth1=log((max(rm34,rsqm-ctnmin)*max(rm34,rsqm-ctpmin))/
     &  (max(rm34,rsqm-ctnmax)*max(rm34,rsqm-ctpmax)))
        acth2=log((max(rm34,rsqm+ctnmax)*max(rm34,rsqm+ctpmax))/
     &  (max(rm34,rsqm+ctnmin)*max(rm34,rsqm+ctpmin)))
        acth3=1d0/max(rm34,rsqm-ctnmax)-1d0/max(rm34,rsqm-ctnmin)+
     &  1d0/max(rm34,rsqm-ctpmax)-1d0/max(rm34,rsqm-ctpmin)
        acth4=1d0/max(rm34,rsqm+ctnmin)-1d0/max(rm34,rsqm+ctnmax)+
     &  1d0/max(rm34,rsqm+ctpmin)-1d0/max(rm34,rsqm+ctpmax)
        h3=coef(isubsv,13)+
     &  (acth0/acth1)*coef(isubsv,14)/max(rm34,rsqm-cth)+
     &  (acth0/acth2)*coef(isubsv,15)/max(rm34,rsqm+cth)+
     &  (acth0/acth3)*coef(isubsv,16)/max(rm34,rsqm-cth)**2+
     &  (acth0/acth4)*coef(isubsv,17)/max(rm34,rsqm+cth)**2
        comfac=comfac*acth0*0.5d0*be34/h3
 
C...2 -> 2 processes: take into account final state Breit-Wigners
        comfac=comfac*vint(80)
      endif
 
C...2 -> 3, 4 processes: phace space integral in tau'
      if(mint(47).ge.2.and.istsb.ge.3.and.istsb.le.5) then
        ataup1=log(taupmx/taupmn)
        ataup2=((1d0-tau/taupmx)**4-(1d0-tau/taupmn)**4)/(4d0*tau)
        h4=coef(isubsv,18)+
     &  (ataup1/ataup2)*coef(isubsv,19)*(1d0-tau/taup)**3/taup
        if(mint(47).eq.5) then
          ataup3=log(max(2d-6,1d0-taupmn)/max(2d-6,1d0-taupmx))
          h4=h4+(ataup1/ataup3)*coef(isubsv,20)*taup/max(2d-6,1d0-taup)
        endif
        comfac=comfac*ataup1/h4
      endif
 
C...2 -> 3, 4 processes: effective W/Z parton distributions
      if(istsb.eq.3.or.istsb.eq.4) then
        if(1d0-tau/taup.gt.1.d-4) then
          fzw=(1d0+tau/taup)*log(taup/tau)-2d0*(1d0-tau/taup)
        else
          fzw=1d0/6d0*(1d0-tau/taup)**3*tau/taup
        endif
        comfac=comfac*fzw
      endif
 
C...2 -> 3 processes: phase space integrals for pT1, pT2, y3, mirror
      if(istsb.eq.5) then
        comfac=comfac*vint(205)*vint(210)*vint(212)*vint(214)/
     &  (128d0*paru(1)**4*vint(220))*(tau**2/taup)
      endif
 
C...2 -> 2 processes: optional dampening by pT^4/(pT0^2+pT^2)^2
      if(mstp(85).eq.1.and.mod(istsb,2).eq.0) comfac=comfac*
     &sqpth**2/(parp(82)**2+sqpth)**2
 
C...gamma + gamma: include factor 2 when different nature
      if(mint(11).eq.22.and.mint(12).eq.22.and.mint(123).ge.4)
     &comfac=2d0*comfac
 
C...Phase space integral for low-pT and multiple interactions
      if(istsb.eq.9) then
        comfac=paru(1)*paru(5)*fack*0.5d0*vint(2)/sh2
        atau1=log(2d0*(1d0+sqrt(1d0-xt2))/xt2-1d0)
        atau2=2d0*atan(1d0/xt2-1d0)/sqrt(xt2)
        h1=coef(isubsv,1)+(atau1/atau2)*coef(isubsv,2)/sqrt(tau)
        comfac=comfac*atau1/h1
        ayst0=ystmax-ystmin
        ayst1=0.5d0*(ystmax-ystmin)**2
        ayst3=2d0*(atan(exp(ystmax))-atan(exp(ystmin)))
        h2=(ayst0/ayst1)*coef(isubsv,8)*(yst-ystmin)+
     &  (ayst0/ayst1)*coef(isubsv,9)*(ystmax-yst)+
     &  (ayst0/ayst3)*coef(isubsv,10)/cosh(yst)
        comfac=comfac*ayst0/h2
        if(mstp(82).le.1) comfac=comfac*xt2**2*(1d0/vint(149)-1d0)
C...For MSTP(82)>=2 an additional factor (xT2/(xT2+VINT(149))**2 is
C...introduced to make cross-section finite for xT2 -> 0
        if(mstp(82).ge.2) comfac=comfac*xt2**2/(vint(149)*
     &  (1d0+vint(149)))
      endif
 
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron
      if((mstp(46).ge.3.and.mstp(46).le.6).and.(isub.eq.71.or.isub.eq.
     &72.or.isub.eq.73.or.isub.eq.76.or.isub.eq.77)) then
C...Calculate M_R and N_R functions for Higgs-like and QCD-like models
        if(mstp(46).le.4) then
          hdtlh=log(pmas(25,1)/parp(44))
          hdtmr=(4.5d0*paru(1)/sqrt(3d0)-74d0/9d0)/8d0+hdtlh/12d0
          hdtnr=-1d0/18d0+hdtlh/6d0
        else
          hdtnm=0.125d0*(1d0/(288d0*paru(1)**2)+(parp(47)/parp(45))**2)
          hdtlq=log(parp(45)/parp(44))
          hdtmr=-(4d0*paru(1))**2*0.5d0*hdtnm+hdtlq/12d0
          hdtnr=(4d0*paru(1))**2*hdtnm+hdtlq/6d0
        endif
 
C...Calculate lowest and next-to-lowest order partial wave amplitudes
        hdtv=1d0/(16d0*paru(1)*parp(47)**2)
        a00l=sngl(hdtv*sh)
        a20l=-0.5*a00l
        a11l=a00l/6.
        hdtls=log(sh/parp(44)**2)
        a004=sngl((hdtv*sh)**2/(4d0*paru(1)))*
     &  cmplx(sngl((176d0*hdtmr+112d0*hdtnr)/3d0+11d0/27d0-
     &  (50d0/9d0)*hdtls),sngl(4d0*paru(1)))
        a204=sngl((hdtv*sh)**2/(4d0*paru(1)))*
     &  cmplx(sngl(32d0*(hdtmr+2d0*hdtnr)/3d0+25d0/54d0-
     &  (20d0/9d0)*hdtls),sngl(paru(1)))
        a114=sngl((hdtv*sh)**2/(6d0*paru(1)))*
     &  cmplx(sngl(4d0*(-2d0*hdtmr+hdtnr)-1d0/18d0),sngl(paru(1)/6d0))
 
C...Unitarize partial wave amplitudes with Pade or K-matrix method
        if(mstp(46).eq.3.or.mstp(46).eq.5) then
          a00u=a00l/(1.-a004/a00l)
          a20u=a20l/(1.-a204/a20l)
          a11u=a11l/(1.-a114/a11l)
        else
          a00u=(a00l+real(a004))/(1.-cmplx(0.,a00l+real(a004)))
          a20u=(a20l+real(a204))/(1.-cmplx(0.,a20l+real(a204)))
          a11u=(a11l+real(a114))/(1.-cmplx(0.,a11l+real(a114)))
        endif
      endif
 
C...Supersymmetric processes - all of type 2 -> 2 :
C...correct final-state Breit-Wigners from fixed to running width.
      if(isub.ge.200.and.isub.le.280.and.mstp(42).gt.0) then
        do 160 i=1,2
        kflw=kfpr(isubsv,i)
        kcw=jamcomp(kflw)
        if(pmas(kcw,2).lt.parp(41)) goto 160
        if(i.eq.1) sqmi=sqm3
        if(i.eq.2) sqmi=sqm4
        sqms=pmas(kcw,1)**2
        gmms=pmas(kcw,1)*pmas(kcw,2)
        hbws=gmms/((sqmi-sqms)**2+gmms**2)
        call pjwidt(kflw,sqmi,wdtp,wdte)
        gmmi=sqrt(sqmi)*wdtp(0)
        hbwi=gmmi/((sqmi-sqms)**2+gmmi**2)
        comfac=comfac*(hbwi/hbws)
  160   continue
      endif
 
C...A: 2 -> 1, tree diagrams
 
  170 if(isub.le.10) then
        if(isub.eq.1) then
C...f + fbar -> gamma*/Z0
          mint(61)=2
          call pjwidt(23,sh,wdtp,wdte)
          hs=shr*wdtp(0)
          facz=4d0*comfac*3d0
          hp0=aem/3d0*sh
          hp1=aem/3d0*xwc*sh
          do 180 i=mmina,mmaxa
            if(i.eq.0.or.kfac(1,i)*kfac(2,-i).eq.0) goto 180
            ei=kchg(iabs(i),1)/3d0
            ai=sign(1d0,ei)
            vi=ai-4d0*ei*xwv
            hi0=hp0
            if(iabs(i).le.10) hi0=hi0*faca/3d0
            hi1=hp1
            if(iabs(i).le.10) hi1=hi1*faca/3d0
            nchn=nchn+1
            isig(nchn,1)=i
            isig(nchn,2)=-i
            isig(nchn,3)=1
            sigh(nchn)=facz*(ei**2/sh2*hi0*hp0*vint(111)+
     &      ei*vi*(1d0-sqmz/sh)/((sh-sqmz)**2+hs**2)*
     &      (hi0*hp1+hi1*hp0)*vint(112)+(vi**2+ai**2)/
     &      ((sh-sqmz)**2+hs**2)*hi1*hp1*vint(114))
  180     continue
 
        elseif(isub.eq.2) then
C...f + fbar' -> W+/-
          call pjwidt(24,sh,wdtp,wdte)
          hs=shr*wdtp(0)
          facbw=4d0*comfac/((sh-sqmw)**2+hs**2)*3d0
          hp=aem/(24d0*xw)*sh
          do 200 i=mmin1,mmax1
            if(i.eq.0.or.kfac(1,i).eq.0) goto 200
            ia=iabs(i)
            do 190 j=mmin2,mmax2
              if(j.eq.0.or.kfac(2,j).eq.0) goto 190
              ja=iabs(j)
              if(i*j.gt.0.or.mod(ia+ja,2).eq.0) goto 190
              if((ia.le.10.and.ja.gt.10).or.(ia.gt.10.and.ja.le.10))
     &        goto 190
              kchw=(kchg(ia,1)*isign(1,i)+kchg(ja,1)*isign(1,j))/3
              hi=hp*2d0
              if(ia.le.10) hi=hi*vckm((ia+1)/2,(ja+1)/2)*faca/3d0
              nchn=nchn+1
              isig(nchn,1)=i
              isig(nchn,2)=j
              isig(nchn,3)=1
              hf=shr*(wdte(0,1)+wdte(0,(5-kchw)/2)+wdte(0,4))
              sigh(nchn)=hi*facbw*hf
  190       continue
  200     continue
 
        elseif(isub.eq.3) then
C...f + fbar -> h0 (or H0, or A0)
          call pjwidt(kfhigg,sh,wdtp,wdte)
          hs=shr*wdtp(0)
          facbw=4d0*comfac/((sh-sqmh)**2+hs**2)
          if(abs(shr-pmas(kfhigg,1)).gt.parp(48)*pmas(kfhigg,2))
     &    facbw=0d0
          hp=aem/(8d0*xw)*sh/sqmw*sh
          hf=shr*(wdte(0,1)+wdte(0,2)+wdte(0,4))
          do 210 i=mmina,mmaxa
            if(i.eq.0.or.kfac(1,i)*kfac(2,-i).eq.0) goto 210
            ia=iabs(i)
            rmq=pmas(ia,1)**2/sh
            hi=hp*rmq
            if(ia.le.10) hi=hp*rmq*faca/3d0
            if(ia.le.10.and.mstp(37).eq.1.and.mstp(2).ge.1) hi=hi*
     &      (log(max(4d0,parp(37)**2*rmq*sh/paru(117)**2))/
     &      log(max(4d0,sh/paru(117)**2)))**(24d0/(33d0-2d0*mstu(118)))
            if(mstp(4).ge.1.or.ihigg.ge.2) then
              ikfi=1
              if(ia.le.10.and.mod(ia,2).eq.0) ikfi=2
              if(ia.gt.10) ikfi=3
              hi=hi*paru(150+10*ihigg+ikfi)**2
            endif
            nchn=nchn+1
            isig(nchn,1)=i
            isig(nchn,2)=-i
            isig(nchn,3)=1
            sigh(nchn)=hi*facbw*hf
  210     continue
 
        elseif(isub.eq.4) then
C...gamma + W+/- -> W+/-
 
        elseif(isub.eq.5) then
C...Z0 + Z0 -> h0
          call pjwidt(25,sh,wdtp,wdte)
          hs=shr*wdtp(0)
          facbw=4d0*comfac/((sh-sqmh)**2+hs**2)
          if(abs(shr-pmas(25,1)).gt.parp(48)*pmas(25,2)) facbw=0d0
          hp=aem/(8d0*xw)*sh/sqmw*sh
          hf=shr*(wdte(0,1)+wdte(0,2)+wdte(0,4))
          hi=hp/4d0
          faci=8d0/(paru(1)**2*xw1)*(aem*xwc)**2
          do 230 i=mmin1,mmax1
            if(i.eq.0.or.kfac(1,i).eq.0) goto 230
            do 220 j=mmin2,mmax2
              if(j.eq.0.or.kfac(2,j).eq.0) goto 220
              ei=kchg(iabs(i),1)/3d0
              ai=sign(1d0,ei)
              vi=ai-4d0*ei*xwv
              ej=kchg(iabs(j),1)/3d0
              aj=sign(1d0,ej)
              vj=aj-4d0*ej*xwv
              nchn=nchn+1
              isig(nchn,1)=i
              isig(nchn,2)=j
              isig(nchn,3)=1
              sigh(nchn)=faci*(vi**2+ai**2)*(vj**2+aj**2)*hi*facbw*hf
  220       continue
  230     continue
 
        elseif(isub.eq.6) then
C...Z0 + W+/- -> W+/-
 
        elseif(isub.eq.7) then
C...W+ + W- -> Z0
 
        elseif(isub.eq.8) then
C...W+ + W- -> h0
          call pjwidt(25,sh,wdtp,wdte)
          hs=shr*wdtp(0)
          facbw=4d0*comfac/((sh-sqmh)**2+hs**2)
          if(abs(shr-pmas(25,1)).gt.parp(48)*pmas(25,2)) facbw=0d0
          hp=aem/(8d0*xw)*sh/sqmw*sh
          hf=shr*(wdte(0,1)+wdte(0,2)+wdte(0,4))
          hi=hp/2d0
          faci=1d0/(4d0*paru(1)**2)*(aem/xw)**2
          do 250 i=mmin1,mmax1
            if(i.eq.0.or.kfac(1,i).eq.0) goto 250
            ei=sign(1d0,dble(i))*kchg(iabs(i),1)
            do 240 j=mmin2,mmax2
              if(j.eq.0.or.kfac(2,j).eq.0) goto 240
              ej=sign(1d0,dble(j))*kchg(iabs(j),1)
              if(ei*ej.gt.0d0) goto 240
              nchn=nchn+1
              isig(nchn,1)=i
              isig(nchn,2)=j
              isig(nchn,3)=1
              sigh(nchn)=faci*vint(180+i)*vint(180+j)*hi*facbw*hf
  240       continue
  250     continue
 
C...B: 2 -> 2, tree diagrams
 
        elseif(isub.eq.10) then
C...f + f' -> f + f' (gamma/Z/W exchange)
          facggf=comfac*aem**2*2d0*(sh2+uh2)/th2
          facgzf=comfac*aem**2*xwc*4d0*sh2/(th*(th-sqmz))
          faczzf=comfac*(aem*xwc)**2*2d0*sh2/(th-sqmz)**2
          facwwf=comfac*(0.5d0*aem/xw)**2*sh2/(th-sqmw)**2
          do 270 i=mmin1,mmax1
            if(i.eq.0.or.kfac(1,i).eq.0) goto 270
            ia=iabs(i)
            do 260 j=mmin2,mmax2
              if(j.eq.0.or.kfac(2,j).eq.0) goto 260
              ja=iabs(j)
C...Electroweak couplings
              ei=kchg(ia,1)*isign(1,i)/3d0
              ai=sign(1d0,kchg(ia,1)+0.5d0)*isign(1,i)
              vi=ai-4d0*ei*xwv
              ej=kchg(ja,1)*isign(1,j)/3d0
              aj=sign(1d0,kchg(ja,1)+0.5d0)*isign(1,j)
              vj=aj-4d0*ej*xwv
              epsij=isign(1,i*j)
C...gamma/Z exchange, only gamma exchange, or only Z exchange
              if(mstp(21).ge.1.and.mstp(21).le.4) then
                if(mstp(21).eq.1.or.mstp(21).eq.4) then
                  facncf=facggf*ei**2*ej**2+facgzf*ei*ej*
     &            (vi*vj*(1d0+uh2/sh2)+ai*aj*epsij*(1d0-uh2/sh2))+
     &            faczzf*((vi**2+ai**2)*(vj**2+aj**2)*(1d0+uh2/sh2)+
     &            4d0*vi*vj*ai*aj*epsij*(1d0-uh2/sh2))
                elseif(mstp(21).eq.2) then
                  facncf=facggf*ei**2*ej**2
                else
                  facncf=faczzf*((vi**2+ai**2)*(vj**2+aj**2)*
     &            (1d0+uh2/sh2)+4d0*vi*vj*ai*aj*epsij*(1d0-uh2/sh2))
                endif
                nchn=nchn+1
                isig(nchn,1)=i
                isig(nchn,2)=j
                isig(nchn,3)=1
                sigh(nchn)=facncf
              endif
C...W exchange
              if((mstp(21).eq.1.or.mstp(21).eq.5).and.ai*aj.lt.0d0) then
                facccf=facwwf*vint(180+i)*vint(180+j)
                if(epsij.lt.0d0) facccf=facccf*uh2/sh2
                if(ia.gt.10.and.mod(ia,2).eq.0) facccf=2d0*facccf
                if(ja.gt.10.and.mod(ja,2).eq.0) facccf=2d0*facccf
                nchn=nchn+1
                isig(nchn,1)=i
                isig(nchn,2)=j
                isig(nchn,3)=2
                sigh(nchn)=facccf
              endif
  260       continue
  270     continue
        endif
 
      elseif(isub.le.20) then
        if(isub.eq.11) then
C...f + f' -> f + f' (g exchange)
          facqq1=comfac*as**2*4d0/9d0*(sh2+uh2)/th2
          facqqb=comfac*as**2*4d0/9d0*((sh2+uh2)/th2*faca-
     &    mstp(34)*2d0/3d0*uh2/(sh*th))
          facqq2=comfac*as**2*4d0/9d0*((sh2+th2)/uh2-
     &    mstp(34)*2d0/3d0*sh2/(th*uh))
          if(mstp(5).ge.1) then
C...Modifications from contact interactions (compositeness)
            facci1=facqq1+comfac*(sh2/paru(155)**4)
            faccib=facqqb+comfac*(8d0/9d0)*(as*paru(156)/paru(155)**2)*
     &      (uh2/th+uh2/sh)+comfac*(5d0/3d0)*(uh2/paru(155)**4)
            facci2=facqq2+comfac*(8d0/9d0)*(as*paru(156)/paru(155)**2)*
     &      (sh2/th+sh2/uh)+comfac*(5d0/3d0)*(sh2/paru(155)**4)
            facci3=facqq1+comfac*(uh2/paru(155)**4)
          endif
          do 290 i=mmin1,mmax1
            ia=iabs(i)
            if(i.eq.0.or.ia.gt.mstp(58).or.kfac(1,i).eq.0) goto 290
            do 280 j=mmin2,mmax2
              ja=iabs(j)
              if(j.eq.0.or.ja.gt.mstp(58).or.kfac(2,j).eq.0) goto 280
              nchn=nchn+1
              isig(nchn,1)=i
              isig(nchn,2)=j
              isig(nchn,3)=1
              if(mstp(5).le.0.or.(mstp(5).eq.1.and.(ia.ge.3.or.
     &        ja.ge.3))) then
                sigh(nchn)=facqq1
                if(i.eq.-j) sigh(nchn)=facqqb
              else
                sigh(nchn)=facci1
                if(i*j.lt.0) sigh(nchn)=facci3
                if(i.eq.-j) sigh(nchn)=faccib
              endif
              if(i.eq.j) then
                sigh(nchn)=0.5d0*sigh(nchn)
                nchn=nchn+1
                isig(nchn,1)=i
                isig(nchn,2)=j
                isig(nchn,3)=2
                if(mstp(5).le.0.or.(mstp(5).eq.1.and.ia.ge.3)) then
                  sigh(nchn)=0.5d0*facqq2
                else
                  sigh(nchn)=0.5d0*facci2
                endif
              endif
  280       continue
  290     continue
 
        elseif(isub.eq.12) then
C...f + fbar -> f' + fbar' (q + qbar -> q' + qbar' only)
          call pjwidt(21,sh,wdtp,wdte)
          facqqb=comfac*as**2*4d0/9d0*(th2+uh2)/sh2*
     &    (wdte(0,1)+wdte(0,2)+wdte(0,4))
          if(mstp(5).eq.1) then
C...Modifications from contact interactions (compositeness)
            faccib=facqqb
            do 300 i=1,2
              faccib=faccib+comfac*(uh2/paru(155)**4)*(wdte(i,1)+
     &        wdte(i,2)+wdte(i,4))
  300       continue
          elseif(mstp(5).ge.2) then
            faccib=facqqb+comfac*(uh2/paru(155)**4)*
     &      (wdte(0,1)+wdte(0,2)+wdte(0,4))
          endif
          do 310 i=mmina,mmaxa
            if(i.eq.0.or.iabs(i).gt.mstp(58).or.
     &      kfac(1,i)*kfac(2,-i).eq.0) goto 310
            nchn=nchn+1
            isig(nchn,1)=i
            isig(nchn,2)=-i
            isig(nchn,3)=1
            if(mstp(5).le.0.or.(mstp(5).eq.1.and.iabs(i).ge.3)) then
              sigh(nchn)=facqqb
            else
              sigh(nchn)=faccib
            endif
  310     continue
 
        elseif(isub.eq.13) then
C...f + fbar -> g + g (q + qbar -> g + g only)
          facgg1=comfac*as**2*32d0/27d0*(uh/th-(2d0+mstp(34)*1d0/4d0)*
     &    uh2/sh2)
          facgg2=comfac*as**2*32d0/27d0*(th/uh-(2d0+mstp(34)*1d0/4d0)*
     &    th2/sh2)
          do 320 i=mmina,mmaxa
            if(i.eq.0.or.iabs(i).gt.mstp(58).or.
     &      kfac(1,i)*kfac(2,-i).eq.0) goto 320
            nchn=nchn+1
            isig(nchn,1)=i
            isig(nchn,2)=-i
            isig(nchn,3)=1
            sigh(nchn)=0.5d0*facgg1
            nchn=nchn+1
            isig(nchn,1)=i
            isig(nchn,2)=-i
            isig(nchn,3)=2
            sigh(nchn)=0.5d0*facgg2
  320     continue
 
        elseif(isub.eq.14) then
C...f + fbar -> g + gamma (q + qbar -> g + gamma only)
          facgg=comfac*as*aem*8d0/9d0*(th2+uh2)/(th*uh)
          do 330 i=mmina,mmaxa
            if(i.eq.0.or.iabs(i).gt.mstp(58).or.
     &      kfac(1,i)*kfac(2,-i).eq.0) goto 330
            ei=kchg(iabs(i),1)/3d0
            nchn=nchn+1
            isig(nchn,1)=i
            isig(nchn,2)=-i
            isig(nchn,3)=1
            sigh(nchn)=facgg*ei**2
  330     continue
 
        elseif(isub.eq.15) then
C...f + fbar -> g + (gamma*/Z0) (q + qbar -> g + (gamma*/Z0) only)
          faczg=comfac*as*aem*(8d0/9d0)*(th2+uh2+2d0*sqm4*sh)/(th*uh)
C...gamma, gamma/Z interference and Z couplings to final fermion pairs
          hfgg=0d0
          hfgz=0d0
          hfzz=0d0
          radc4=1d0+pjalps(sqm4)/paru(1)
          do 340 i=1,min(16,mdcy(23,3))
            idc=i+mdcy(23,2)-1
            if(mdme(idc,1).lt.0) goto 340
            imdm=0
            if(mdme(idc,1).eq.1.or.mdme(idc,1).eq.2.or.mdme(idc,1).eq.4)
     &      imdm=1
            if(i.le.8) then
              ef=kchg(i,1)/3d0
              af=sign(1d0,ef+0.1d0)
              vf=af-4d0*ef*xwv
            elseif(i.le.16) then
              ef=kchg(i+2,1)/3d0
              af=sign(1d0,ef+0.1d0)
              vf=af-4d0*ef*xwv
            endif
            rm1=pmas(iabs(kfdp(idc,1)),1)**2/sqm4
            if(4d0*rm1.lt.1d0) then
              fcof=1d0
              if(i.le.8) fcof=3d0*radc4
              be34=sqrt(max(0d0,1d0-4d0*rm1))
              if(imdm.eq.1) then
                hfgg=hfgg+fcof*ef**2*(1d0+2d0*rm1)*be34
                hfgz=hfgz+fcof*ef*vf*(1d0+2d0*rm1)*be34
                hfzz=hfzz+fcof*(vf**2*(1d0+2d0*rm1)+
     &          af**2*(1d0-4d0*rm1))*be34
              endif
            endif
  340     continue
C...Propagators: as simulated in PYOFSH and as desired
          hbw4=(1d0/paru(1))*gmmz/((sqm4-sqmz)**2+gmmz**2)
          mint(15)=1
          mint(61)=1
          call pjwidt(23,sqm4,wdtp,wdte)
          hfaem=(paru(108)/paru(2))*(2d0/3d0)
          hfgg=hfgg*hfaem*vint(111)/sqm4
          hfgz=hfgz*hfaem*vint(112)/sqm4
          hfzz=hfzz*hfaem*vint(114)/sqm4
C...Loop over flavours; consider full gamma/Z structure
          do 350 i=mmina,mmaxa
            if(i.eq.0.or.iabs(i).gt.mstp(58).or.
     &      kfac(1,i)*kfac(2,-i).eq.0) goto 350
            ei=kchg(iabs(i),1)/3d0
            ai=sign(1d0,ei)
            vi=ai-4d0*ei*xwv
            nchn=nchn+1
            isig(nchn,1)=i
            isig(nchn,2)=-i
            isig(nchn,3)=1
            sigh(nchn)=faczg*(ei**2*hfgg+ei*vi*hfgz+
     &      (vi**2+ai**2)*hfzz)/hbw4
  350     continue
 
        elseif(isub.eq.16) then
C...f + fbar' -> g + W+/- (q + qbar' -> g + W+/- only)
          facwg=comfac*as*aem/xw*2d0/9d0*(th2+uh2+2d0*sqm4*sh)/(th*uh)
C...Propagators: as simulated in PYOFSH and as desired
          hbw4=gmmw/((sqm4-sqmw)**2+gmmw**2)
          call pjwidt(24,sqm4,wdtp,wdte)
          gmmwc=sqrt(sqm4)*wdtp(0)
          hbw4c=gmmwc/((sqm4-sqmw)**2+gmmwc**2)
          facwg=facwg*hbw4c/hbw4
          do 370 i=mmin1,mmax1
            ia=iabs(i)
            if(i.eq.0.or.ia.gt.10.or.kfac(1,i).eq.0) goto 370
            do 360 j=mmin2,mmax2
              ja=iabs(j)
              if(j.eq.0.or.ja.gt.10.or.kfac(2,j).eq.0) goto 360
              if(i*j.gt.0.or.mod(ia+ja,2).eq.0) goto 360
              kchw=(kchg(ia,1)*isign(1,i)+kchg(ja,1)*isign(1,j))/3
              widsc=(wdte(0,1)+wdte(0,(5-kchw)/2)+wdte(0,4))/wdtp(0)
              fckm=vckm((ia+1)/2,(ja+1)/2)
              nchn=nchn+1
              isig(nchn,1)=i
              isig(nchn,2)=j
              isig(nchn,3)=1
              sigh(nchn)=facwg*fckm*widsc
  360       continue
  370     continue
 
        elseif(isub.eq.17) then
C...f + fbar -> g + h0 (q + qbar -> g + h0 only)
 
        elseif(isub.eq.18) then
C...f + fbar -> gamma + gamma
          facgg=comfac*aem**2*2d0*(th2+uh2)/(th*uh)
          do 380 i=mmina,mmaxa
            if(i.eq.0.or.kfac(1,i)*kfac(2,-i).eq.0) goto 380
            ei=kchg(iabs(i),1)/3d0
            fcoi=1d0
            if(iabs(i).le.10) fcoi=faca/3d0
            nchn=nchn+1
            isig(nchn,1)=i
            isig(nchn,2)=-i
            isig(nchn,3)=1
            sigh(nchn)=0.5d0*facgg*fcoi*ei**4
  380     continue
 
        elseif(isub.eq.19) then
C...f + fbar -> gamma + (gamma*/Z0)
          facgz=comfac*2d0*aem**2*(th2+uh2+2d0*sqm4*sh)/(th*uh)
C...gamma, gamma/Z interference and Z couplings to final fermion pairs
          hfgg=0d0
          hfgz=0d0
          hfzz=0d0
          radc4=1d0+pjalps(sqm4)/paru(1)
          do 390 i=1,min(16,mdcy(23,3))
            idc=i+mdcy(23,2)-1
            if(mdme(idc,1).lt.0) goto 390
            imdm=0
            if(mdme(idc,1).eq.1.or.mdme(idc,1).eq.2.or.mdme(idc,1).eq.4)
     &      imdm=1
            if(i.le.8) then
              ef=kchg(i,1)/3d0
              af=sign(1d0,ef+0.1d0)
              vf=af-4d0*ef*xwv
            elseif(i.le.16) then
              ef=kchg(i+2,1)/3d0
              af=sign(1d0,ef+0.1d0)
              vf=af-4d0*ef*xwv
            endif
            rm1=pmas(iabs(kfdp(idc,1)),1)**2/sqm4
            if(4d0*rm1.lt.1d0) then
              fcof=1d0
              if(i.le.8) fcof=3d0*radc4
              be34=sqrt(max(0d0,1d0-4d0*rm1))
              if(imdm.eq.1) then
                hfgg=hfgg+fcof*ef**2*(1d0+2d0*rm1)*be34
                hfgz=hfgz+fcof*ef*vf*(1d0+2d0*rm1)*be34
                hfzz=hfzz+fcof*(vf**2*(1d0+2d0*rm1)+
     &          af**2*(1d0-4d0*rm1))*be34
              endif
            endif
  390     continue
C...Propagators: as simulated in PYOFSH and as desired
          hbw4=(1d0/paru(1))*gmmz/((sqm4-sqmz)**2+gmmz**2)
          mint(15)=1
          mint(61)=1
          call pjwidt(23,sqm4,wdtp,wdte)
          hfaem=(paru(108)/paru(2))*(2d0/3d0)
          hfgg=hfgg*hfaem*vint(111)/sqm4
          hfgz=hfgz*hfaem*vint(112)/sqm4
          hfzz=hfzz*hfaem*vint(114)/sqm4
C...Loop over flavours; consider full gamma/Z structure
          do 400 i=mmina,mmaxa
            if(i.eq.0.or.kfac(1,i)*kfac(2,-i).eq.0) goto 400
            ei=kchg(iabs(i),1)/3d0
            ai=sign(1d0,ei)
            vi=ai-4d0*ei*xwv
            fcoi=1d0
            if(iabs(i).le.10) fcoi=faca/3d0
            nchn=nchn+1
            isig(nchn,1)=i
            isig(nchn,2)=-i
            isig(nchn,3)=1
            sigh(nchn)=facgz*fcoi*ei**2*(ei**2*hfgg+ei*vi*hfgz+
     &      (vi**2+ai**2)*hfzz)/hbw4
  400     continue
 
        elseif(isub.eq.20) then
C...f + fbar' -> gamma + W+/-
          facgw=comfac*0.5d0*aem**2/xw
C...Propagators: as simulated in PYOFSH and as desired
          hbw4=gmmw/((sqm4-sqmw)**2+gmmw**2)
          call pjwidt(24,sqm4,wdtp,wdte)
          gmmwc=sqrt(sqm4)*wdtp(0)
          hbw4c=gmmwc/((sqm4-sqmw)**2+gmmwc**2)
          facgw=facgw*hbw4c/hbw4
C...Anomalous couplings
          term1=(th2+uh2+2d0*sqm4*sh)/(th*uh)
          term2=0d0
          term3=0d0
          if(mstp(5).ge.1) then
            term2=paru(153)*(th-uh)/(th+uh)
            term3=0.5d0*paru(153)**2*(th*uh+(th2+uh2)*sh/
     &      (4d0*sqmw))/(th+uh)**2
          endif
          do 420 i=mmin1,mmax1
            ia=iabs(i)
            if(i.eq.0.or.ia.gt.20.or.kfac(1,i).eq.0) goto 420
            do 410 j=mmin2,mmax2
              ja=iabs(j)
              if(j.eq.0.or.ja.gt.20.or.kfac(2,j).eq.0) goto 410
              if(i*j.gt.0.or.mod(ia+ja,2).eq.0) goto 410
              if((ia.le.10.and.ja.gt.10).or.(ia.gt.10.and.ja.le.10))
     &        goto 410
              kchw=(kchg(ia,1)*isign(1,i)+kchg(ja,1)*isign(1,j))/3
              widsc=(wdte(0,1)+wdte(0,(5-kchw)/2)+wdte(0,4))/wdtp(0)
              if(ia.le.10) then
                facwr=uh/(th+uh)-1d0/3d0
                fckm=vckm((ia+1)/2,(ja+1)/2)
                fcoi=faca/3d0
              else
                facwr=-th/(th+uh)
                fckm=1d0
                fcoi=1d0
              endif
              facwk=term1*facwr**2+term2*facwr+term3
              nchn=nchn+1
              isig(nchn,1)=i
              isig(nchn,2)=j
              isig(nchn,3)=1
              sigh(nchn)=facgw*facwk*fcoi*fckm*widsc
  410       continue
  420     continue
        endif
 
      elseif(isub.le.30) then
        if(isub.eq.21) then
C...f + fbar -> gamma + h0
 
        elseif(isub.eq.22) then
C...f + fbar -> (gamma*/Z0) + (gamma*/Z0)
C...Kinematics dependence
          faczz=comfac*aem**2*((th2+uh2+2d0*(sqm3+sqm4)*sh)/(th*uh)-
     &    sqm3*sqm4*(1d0/th2+1d0/uh2))
C...gamma, gamma/Z interference and Z couplings to final fermion pairs
          do 440 i=1,6
            do 430 j=1,3
              hgz(i,j)=0d0
  430       continue
  440     continue
          radc3=1d0+pjalps(sqm3)/paru(1)
          radc4=1d0+pjalps(sqm4)/paru(1)
          do 450 i=1,min(16,mdcy(23,3))
            idc=i+mdcy(23,2)-1
            if(mdme(idc,1).lt.0) goto 450
            imdm=0
            if(mdme(idc,1).eq.1.or.mdme(idc,1).eq.2) imdm=1
            if(mdme(idc,1).eq.4.or.mdme(idc,1).eq.5) imdm=mdme(idc,1)-2
            if(i.le.8) then
              ef=kchg(i,1)/3d0
              af=sign(1d0,ef+0.1d0)
              vf=af-4d0*ef*xwv
            elseif(i.le.16) then
              ef=kchg(i+2,1)/3d0
              af=sign(1d0,ef+0.1d0)
              vf=af-4d0*ef*xwv
            endif
            rm1=pmas(iabs(kfdp(idc,1)),1)**2/sqm3
            if(4d0*rm1.lt.1d0) then
              fcof=1d0
              if(i.le.8) fcof=3d0*radc3
              be34=sqrt(max(0d0,1d0-4d0*rm1))
              if(imdm.ge.1) then
                hgz(1,imdm)=hgz(1,imdm)+fcof*ef**2*(1d0+2d0*rm1)*be34
                hgz(2,imdm)=hgz(2,imdm)+fcof*ef*vf*(1d0+2d0*rm1)*be34
                hgz(3,imdm)=hgz(3,imdm)+fcof*(vf**2*(1d0+2d0*rm1)+
     &          af**2*(1d0-4d0*rm1))*be34
              endif
            endif
            rm1=pmas(iabs(kfdp(idc,1)),1)**2/sqm4
            if(4d0*rm1.lt.1d0) then
              fcof=1d0
              if(i.le.8) fcof=3d0*radc4
              be34=sqrt(max(0d0,1d0-4d0*rm1))
              if(imdm.ge.1) then
                hgz(4,imdm)=hgz(4,imdm)+fcof*ef**2*(1d0+2d0*rm1)*be34
                hgz(5,imdm)=hgz(5,imdm)+fcof*ef*vf*(1d0+2d0*rm1)*be34
                hgz(6,imdm)=hgz(6,imdm)+fcof*(vf**2*(1d0+2d0*rm1)+
     &          af**2*(1d0-4d0*rm1))*be34
              endif
            endif
  450     continue
C...Propagators: as simulated in PYOFSH and as desired
          hbw3=(1d0/paru(1))*gmmz/((sqm3-sqmz)**2+gmmz**2)
          hbw4=(1d0/paru(1))*gmmz/((sqm4-sqmz)**2+gmmz**2)
          mint(15)=1
          mint(61)=1
          call pjwidt(23,sqm3,wdtp,wdte)
          hfaem=(paru(108)/paru(2))*(2d0/3d0)
          do 460 j=1,3
            hgz(1,j)=hgz(1,j)*hfaem*vint(111)/sqm3
            hgz(2,j)=hgz(2,j)*hfaem*vint(112)/sqm3
            hgz(3,j)=hgz(3,j)*hfaem*vint(114)/sqm3
  460     continue
          mint(61)=1
          call pjwidt(23,sqm4,wdtp,wdte)
          hfaem=(paru(108)/paru(2))*(2d0/3d0)
          do 470 j=1,3
            hgz(4,j)=hgz(4,j)*hfaem*vint(111)/sqm4
            hgz(5,j)=hgz(5,j)*hfaem*vint(112)/sqm4
            hgz(6,j)=hgz(6,j)*hfaem*vint(114)/sqm4
  470     continue
C...Loop over flavours; separate left- and right-handed couplings
          do 490 i=mmina,mmaxa
            if(i.eq.0.or.kfac(1,i)*kfac(2,-i).eq.0) goto 490
            ei=kchg(iabs(i),1)/3d0
            ai=sign(1d0,ei)
            vi=ai-4d0*ei*xwv
            vali=vi-ai
            vari=vi+ai
            fcoi=1d0
            if(iabs(i).le.10) fcoi=faca/3d0
            do 480 j=1,3
              hl3(j)=ei**2*hgz(1,j)+ei*vali*hgz(2,j)+vali**2*hgz(3,j)
              hr3(j)=ei**2*hgz(1,j)+ei*vari*hgz(2,j)+vari**2*hgz(3,j)
              hl4(j)=ei**2*hgz(4,j)+ei*vali*hgz(5,j)+vali**2*hgz(6,j)
              hr4(j)=ei**2*hgz(4,j)+ei*vari*hgz(5,j)+vari**2*hgz(6,j)
  480       continue
            faclr=hl3(1)*hl4(1)+hl3(1)*(hl4(2)+hl4(3))+
     &      hl4(1)*(hl3(2)+hl3(3))+hl3(2)*hl4(3)+hl4(2)*hl3(3)+
     &      hr3(1)*hr4(1)+hr3(1)*(hr4(2)+hr4(3))+
     &      hr4(1)*(hr3(2)+hr3(3))+hr3(2)*hr4(3)+hr4(2)*hr3(3)
            nchn=nchn+1
            isig(nchn,1)=i
            isig(nchn,2)=-i
            isig(nchn,3)=1
            sigh(nchn)=0.5d0*faczz*fcoi*faclr/(hbw3*hbw4)
  490     continue
 
        elseif(isub.eq.23) then
C...f + fbar' -> Z0 + W+/-
          faczw=comfac*0.5d0*(aem/xw)**2
          faczw=faczw*wids(23,2)
          thuh=max(th*uh-sqm3*sqm4,sh*ckin(3)**2)
          facbw=1d0/((sh-sqmw)**2+gmmw**2)
          do 510 i=mmin1,mmax1
            ia=iabs(i)
            if(i.eq.0.or.ia.gt.20.or.kfac(1,i).eq.0) goto 510
            do 500 j=mmin2,mmax2
              ja=iabs(j)
              if(j.eq.0.or.ja.gt.20.or.kfac(2,j).eq.0) goto 500
              if(i*j.gt.0.or.mod(ia+ja,2).eq.0) goto 500
              if((ia.le.10.and.ja.gt.10).or.(ia.gt.10.and.ja.le.10))
     &        goto 500
              kchw=(kchg(ia,1)*isign(1,i)+kchg(ja,1)*isign(1,j))/3
              ei=kchg(ia,1)/3d0
              ai=sign(1d0,ei+0.1d0)
              vi=ai-4d0*ei*xwv
              ej=kchg(ja,1)/3d0
              aj=sign(1d0,ej+0.1d0)
              vj=aj-4d0*ej*xwv
              if(vi+ai.gt.0) then
                visav=vi
                aisav=ai
                vi=vj
                ai=aj
                vj=visav
                aj=aisav
              endif
              fckm=1d0
              if(ia.le.10) fckm=vckm((ia+1)/2,(ja+1)/2)
              fcoi=1d0
              if(ia.le.10) fcoi=faca/3d0
              nchn=nchn+1
              isig(nchn,1)=i
              isig(nchn,2)=j
              isig(nchn,3)=1
              sigh(nchn)=faczw*fcoi*fckm*(facbw*((9d0-8d0*xw)/4d0*thuh+
     &        (8d0*xw-6d0)/4d0*sh*(sqm3+sqm4))+(thuh-sh*(sqm3+sqm4))*
     &        (sh-sqmw)*facbw*0.5d0*((vj+aj)/th-(vi+ai)/uh)+
     &        thuh/(16d0*xw1)*((vj+aj)**2/th2+(vi+ai)**2/uh2)+
     &        sh*(sqm3+sqm4)/(8d0*xw1)*(vi+ai)*(vj+aj)/(th*uh))*
     &        wids(24,(5-kchw)/2)
  500       continue
  510     continue
 
        elseif(isub.eq.24) then
C...f + fbar -> Z0 + h0 (or H0, or A0)
          thuh=max(th*uh-sqm3*sqm4,sh*ckin(3)**2)
          fachz=comfac*8d0*(aem*xwc)**2*
     &    (thuh+2d0*sh*sqm3)/((sh-sqmz)**2+gmmz**2)
          fachz=fachz*wids(23,2)*wids(kfhigg,2)
          if(mstp(4).ge.1.or.ihigg.ge.2) fachz=fachz*
     &    paru(154+10*ihigg)**2
          do 520 i=mmina,mmaxa
            if(i.eq.0.or.kfac(1,i)*kfac(2,-i).eq.0) goto 520
            ei=kchg(iabs(i),1)/3d0
            ai=sign(1d0,ei)
            vi=ai-4d0*ei*xwv
            fcoi=1d0
            if(iabs(i).le.10) fcoi=faca/3d0
            nchn=nchn+1
            isig(nchn,1)=i
            isig(nchn,2)=-i
            isig(nchn,3)=1
            sigh(nchn)=fachz*fcoi*(vi**2+ai**2)
  520     continue
 
        elseif(isub.eq.25) then
C...f + fbar -> W+ + W-
C...Propagators: Z0, W+- as simulated in PYOFSH and as desired
          call pjwidt(23,sh,wdtp,wdte)
          gmmzc=shr*wdtp(0)
          hbwzc=sh**2/((sh-sqmz)**2+gmmzc**2)
          hbw3=gmmw/((sqm3-sqmw)**2+gmmw**2)
          call pjwidt(24,sqm3,wdtp,wdte)
          gmmw3=sqrt(sqm3)*wdtp(0)
          hbw3c=gmmw3/((sqm3-sqmw)**2+gmmw3**2)
          hbw4=gmmw/((sqm4-sqmw)**2+gmmw**2)
          call pjwidt(24,sqm4,wdtp,wdte)
          gmmw4=sqrt(sqm4)*wdtp(0)
          hbw4c=gmmw4/((sqm4-sqmw)**2+gmmw4**2)
C...Kinematical functions
          thuh=max(th*uh-sqm3*sqm4,sh*ckin(3)**2)
          thuh34=(2d0*sh*(sqm3+sqm4)+thuh)/(sqm3*sqm4)
          gs=(((sh-sqm3-sqm4)**2-4d0*sqm3*sqm4)*thuh34+12d0*thuh)/sh2
          gt=thuh34+4d0*thuh/th2
          gst=((sh-sqm3-sqm4)*thuh34+4d0*(sh*(sqm3+sqm4)-thuh)/th)/sh
          gu=thuh34+4d0*thuh/uh2
          gsu=((sh-sqm3-sqm4)*thuh34+4d0*(sh*(sqm3+sqm4)-thuh)/uh)/sh
C...Common factors and couplings
          facww=comfac*(hbw3c/hbw3)*(hbw4c/hbw4)
          facww=facww*wids(24,1)
          cgg=aem**2/2d0
          cgz=aem**2/(4d0*xw)*hbwzc*(1d0-sqmz/sh)
          czz=aem**2/(32d0*xw**2)*hbwzc
          cng=aem**2/(4d0*xw)
          cnz=aem**2/(16d0*xw**2)*hbwzc*(1d0-sqmz/sh)
          cnn=aem**2/(16d0*xw**2)
C...Coulomb factor for W+W- pair
          if(mstp(40).ge.1.and.mstp(40).le.3) then
            coule=(sh-4d0*sqmw)/(4d0*pmas(24,1))
            coulp=max(1d-10,0.5d0*be34*sqrt(sh))
            if(coule.lt.100d0*pmas(24,2)) then
              coulp1=sqrt(0.5d0*pmas(24,1)*(sqrt(coule**2+
     &        pmas(24,2)**2)-coule))
            else
              coulp1=sqrt(0.5d0*pmas(24,1)*(0.5d0*pmas(24,2)**2/coule))
            endif
            if(coule.gt.-100d0*pmas(24,2)) then
              coulp2=sqrt(0.5d0*pmas(24,1)*(sqrt(coule**2+
     &        pmas(24,2)**2)+coule))
            else
              coulp2=sqrt(0.5d0*pmas(24,1)*(0.5d0*pmas(24,2)**2/
     &        abs(coule)))
            endif
            if(mstp(40).eq.1) then
              couldc=paru(1)-2d0*atan((coulp1**2+coulp2**2-coulp**2)/
     &        max(1d-10,2d0*coulp*coulp1))
              faccou=1d0+0.5d0*paru(101)*couldc/max(1d-5,be34)
            elseif(mstp(40).eq.2) then
              coulck=cmplx(sngl(coulp1),sngl(coulp2))
              coulcp=cmplx(0.,sngl(coulp))
              coulcd=(coulck+coulcp)/(coulck-coulcp)
              coulcr=1.+sngl(paru(101)*sqrt(sh))/(4.*coulcp)*log(coulcd)
              coulcs=cmplx(0.,0.)
              nstp=100
              do 530 istp=1,nstp
                coulxx=(istp-0.5)/nstp
                coulcs=coulcs+(1./coulxx)*log((1.+coulxx*coulcd)/
     &          (1.+coulxx/coulcd))
  530         continue
              coulcr=coulcr+sngl(paru(101)**2*sh)/(16.*coulcp*coulck)*
     &        (coulcs/nstp)
              faccou=abs(coulcr)**2
            elseif(mstp(40).eq.3) then
              couldc=paru(1)-2d0*(1d0-be34)**2*atan((coulp1**2+
     &        coulp2**2-coulp**2)/max(1d-10,2d0*coulp*coulp1))
              faccou=1d0+0.5d0*paru(101)*couldc/max(1d-5,be34)
            endif
          elseif(mstp(40).eq.4) then
            faccou=1d0+0.5d0*paru(101)*paru(1)/max(1d-5,be34)
          else
            faccou=1d0
          endif
          vint(95)=faccou
          facww=facww*faccou
C...Loop over allowed flavours
          do 540 i=mmina,mmaxa
            if(i.eq.0.or.kfac(1,i)*kfac(2,-i).eq.0) goto 540
            ei=kchg(iabs(i),1)/3d0
            ai=sign(1d0,ei+0.1d0)
            vi=ai-4d0*ei*xwv
            fcoi=1d0
            if(iabs(i).le.10) fcoi=faca/3d0
            if(ai.lt.0d0) then
              dsigww=(cgg*ei**2+cgz*vi*ei+czz*(vi**2+ai**2))*gs+
     &        (cng*ei+cnz*(vi+ai))*gst+cnn*gt
            else
              dsigww=(cgg*ei**2+cgz*vi*ei+czz*(vi**2+ai**2))*gs-
     &        (cng*ei+cnz*(vi+ai))*gsu+cnn*gu
            endif
            nchn=nchn+1
            isig(nchn,1)=i
            isig(nchn,2)=-i
            isig(nchn,3)=1
            sigh(nchn)=facww*fcoi*dsigww
  540     continue
 
        elseif(isub.eq.26) then
C...f + fbar' -> W+/- + h0 (or H0, or A0)
          thuh=max(th*uh-sqm3*sqm4,sh*ckin(3)**2)
          fachw=comfac*0.125d0*(aem/xw)**2*(thuh+2d0*sh*sqm3)/
     &    ((sh-sqmw)**2+gmmw**2)
          fachw=fachw*wids(kfhigg,2)
          if(mstp(4).ge.1.or.ihigg.ge.2) fachw=fachw*
     &    paru(155+10*ihigg)**2
          do 560 i=mmin1,mmax1
            ia=iabs(i)
            if(i.eq.0.or.ia.gt.20.or.kfac(1,i).eq.0) goto 560
            do 550 j=mmin2,mmax2
              ja=iabs(j)
              if(j.eq.0.or.ja.gt.20.or.kfac(1,j).eq.0) goto 550
              if(i*j.gt.0.or.mod(ia+ja,2).eq.0) goto 550
              if((ia.le.10.and.ja.gt.10).or.(ia.gt.10.and.ja.le.10))
     &        goto 550
              kchw=(kchg(ia,1)*isign(1,i)+kchg(ja,1)*isign(1,j))/3
              fckm=1d0
              if(ia.le.10) fckm=vckm((ia+1)/2,(ja+1)/2)
              fcoi=1d0
              if(ia.le.10) fcoi=faca/3d0
              nchn=nchn+1
              isig(nchn,1)=i
              isig(nchn,2)=j
              isig(nchn,3)=1
              sigh(nchn)=fachw*fcoi*fckm*wids(24,(5-kchw)/2)
  550       continue
  560     continue
 
        elseif(isub.eq.27) then
C...f + fbar -> h0 + h0
 
        elseif(isub.eq.28) then
C...f + g -> f + g (q + g -> q + g only)
          facqg1=comfac*as**2*4d0/9d0*((2d0+mstp(34)*1d0/4d0)*uh2/th2-
     &    uh/sh)*faca
          facqg2=comfac*as**2*4d0/9d0*((2d0+mstp(34)*1d0/4d0)*sh2/th2-
     &    sh/uh)
          do 580 i=mmina,mmaxa
            if(i.eq.0.or.iabs(i).gt.10) goto 580
            do 570 isde=1,2
              if(isde.eq.1.and.kfac(1,i)*kfac(2,21).eq.0) goto 570
              if(isde.eq.2.and.kfac(1,21)*kfac(2,i).eq.0) goto 570
              nchn=nchn+1
              isig(nchn,isde)=i
              isig(nchn,3-isde)=21
              isig(nchn,3)=1
              sigh(nchn)=facqg1
              nchn=nchn+1
              isig(nchn,isde)=i
              isig(nchn,3-isde)=21
              isig(nchn,3)=2
              sigh(nchn)=facqg2
  570       continue
  580     continue
 
        elseif(isub.eq.29) then
C...f + g -> f + gamma (q + g -> q + gamma only)
          fgq=comfac*faca*as*aem*1d0/3d0*(sh2+uh2)/(-sh*uh)
          do 600 i=mmina,mmaxa
            if(i.eq.0.or.iabs(i).gt.mstp(58)) goto 600
            ei=kchg(iabs(i),1)/3d0
            facgq=fgq*ei**2
            do 590 isde=1,2
              if(isde.eq.1.and.kfac(1,i)*kfac(2,21).eq.0) goto 590
              if(isde.eq.2.and.kfac(1,21)*kfac(2,i).eq.0) goto 590
              nchn=nchn+1
              isig(nchn,isde)=i
              isig(nchn,3-isde)=21
              isig(nchn,3)=1
              sigh(nchn)=facgq
  590       continue
  600     continue
 
        elseif(isub.eq.30) then
C...f + g -> f + (gamma*/Z0) (q + g -> q + (gamma*/Z0) only)
          fzq=comfac*faca*as*aem*(1d0/3d0)*(sh2+uh2+2d0*sqm4*th)/
     &    (-sh*uh)
C...gamma, gamma/Z interference and Z couplings to final fermion pairs
          hfgg=0d0
          hfgz=0d0
          hfzz=0d0
          radc4=1d0+pjalps(sqm4)/paru(1)
          do 610 i=1,min(16,mdcy(23,3))
            idc=i+mdcy(23,2)-1
            if(mdme(idc,1).lt.0) goto 610
            imdm=0
            if(mdme(idc,1).eq.1.or.mdme(idc,1).eq.2.or.mdme(idc,1).eq.4)
     &      imdm=1
            if(i.le.8) then
              ef=kchg(i,1)/3d0
              af=sign(1d0,ef+0.1d0)
              vf=af-4d0*ef*xwv
            elseif(i.le.16) then
              ef=kchg(i+2,1)/3d0
              af=sign(1d0,ef+0.1d0)
              vf=af-4d0*ef*xwv
            endif
            rm1=pmas(iabs(kfdp(idc,1)),1)**2/sqm4
            if(4d0*rm1.lt.1d0) then
              fcof=1d0
              if(i.le.8) fcof=3d0*radc4
              be34=sqrt(max(0d0,1d0-4d0*rm1))
              if(imdm.eq.1) then
                hfgg=hfgg+fcof*ef**2*(1d0+2d0*rm1)*be34
                hfgz=hfgz+fcof*ef*vf*(1d0+2d0*rm1)*be34
                hfzz=hfzz+fcof*(vf**2*(1d0+2d0*rm1)+
     &          af**2*(1d0-4d0*rm1))*be34
              endif
            endif
  610     continue
C...Propagators: as simulated in PYOFSH and as desired
          hbw4=(1d0/paru(1))*gmmz/((sqm4-sqmz)**2+gmmz**2)
          mint(15)=1
          mint(61)=1
          call pjwidt(23,sqm4,wdtp,wdte)
          hfaem=(paru(108)/paru(2))*(2d0/3d0)
          hfgg=hfgg*hfaem*vint(111)/sqm4
          hfgz=hfgz*hfaem*vint(112)/sqm4
          hfzz=hfzz*hfaem*vint(114)/sqm4
C...Loop over flavours; consider full gamma/Z structure
          do 630 i=mmina,mmaxa
            if(i.eq.0.or.iabs(i).gt.mstp(58)) goto 630
            ei=kchg(iabs(i),1)/3d0
            ai=sign(1d0,ei)
            vi=ai-4d0*ei*xwv
            faczq=fzq*(ei**2*hfgg+ei*vi*hfgz+
     &      (vi**2+ai**2)*hfzz)/hbw4
            do 620 isde=1,2
              if(isde.eq.1.and.kfac(1,i)*kfac(2,21).eq.0) goto 620
              if(isde.eq.2.and.kfac(1,21)*kfac(2,i).eq.0) goto 620
              nchn=nchn+1
              isig(nchn,isde)=i
              isig(nchn,3-isde)=21
              isig(nchn,3)=1
              sigh(nchn)=faczq
  620       continue
  630     continue
        endif
 
      elseif(isub.le.40) then

        if(isub.eq.31) then
C...f + g -> f' + W+/- (q + g -> q' + W+/- only)
          facwq=comfac*faca*as*aem/xw*1d0/12d0*
     &    (sh2+uh2+2d0*sqm4*th)/(-sh*uh)
C...Propagators: as simulated in PYOFSH and as desired
          hbw4=gmmw/((sqm4-sqmw)**2+gmmw**2)
          call pjwidt(24,sqm4,wdtp,wdte)
          gmmwc=sqrt(sqm4)*wdtp(0)
          hbw4c=gmmwc/((sqm4-sqmw)**2+gmmwc**2)
          facwq=facwq*hbw4c/hbw4
          do 650 i=mmina,mmaxa
            if(i.eq.0.or.iabs(i).gt.mstp(58)) goto 650
            ia=iabs(i)
            kchw=isign(1,kchg(ia,1)*isign(1,i))
            widsc=(wdte(0,1)+wdte(0,(5-kchw)/2)+wdte(0,4))/wdtp(0)
            do 640 isde=1,2
              if(isde.eq.1.and.kfac(1,i)*kfac(2,21).eq.0) goto 640
              if(isde.eq.2.and.kfac(1,21)*kfac(2,i).eq.0) goto 640
              nchn=nchn+1
              isig(nchn,isde)=i
              isig(nchn,3-isde)=21
              isig(nchn,3)=1
              sigh(nchn)=facwq*vint(180+i)*widsc
  640       continue
  650     continue
 
        elseif(isub.eq.32) then
C...f + g -> f + h0 (q + g -> q + h0 only)
 
        elseif(isub.eq.33) then
C...f + gamma -> f + g (q + gamma -> q + g only)
          fgq=comfac*as*aem*8d0/3d0*(sh2+uh2)/(-sh*uh)
          do 670 i=mmina,mmaxa
            if(i.eq.0.or.iabs(i).gt.mstp(58)) goto 670
            ei=kchg(iabs(i),1)/3d0
            facgq=fgq*ei**2
            do 660 isde=1,2
              if(isde.eq.1.and.kfac(1,i)*kfac(2,22).eq.0) goto 660
              if(isde.eq.2.and.kfac(1,22)*kfac(2,i).eq.0) goto 660
              nchn=nchn+1
              isig(nchn,isde)=i
              isig(nchn,3-isde)=22
              isig(nchn,3)=1
              sigh(nchn)=facgq
  660       continue
  670     continue
 
        elseif(isub.eq.34) then
C...f + gamma -> f + gamma
          fgq=comfac*aem**2*2d0*(sh2+uh2)/(-sh*uh)
          do 690 i=mmina,mmaxa
            if(i.eq.0) goto 690
            ei=kchg(iabs(i),1)/3d0
            facgq=fgq*ei**4
            do 680 isde=1,2
              if(isde.eq.1.and.kfac(1,i)*kfac(2,22).eq.0) goto 680
              if(isde.eq.2.and.kfac(1,22)*kfac(2,i).eq.0) goto 680
              nchn=nchn+1
              isig(nchn,isde)=i
              isig(nchn,3-isde)=22
              isig(nchn,3)=1
              sigh(nchn)=facgq
  680       continue
  690     continue
 
        elseif(isub.eq.35) then
C...f + gamma -> f + (gamma*/Z0)
          fzqn=comfac*2d0*aem**2*(sh2+uh2+2d0*sqm4*th)
          fzqd=sqpth*sqm4-sh*uh
C...gamma, gamma/Z interference and Z couplings to final fermion pairs
          hfgg=0d0
          hfgz=0d0
          hfzz=0d0
          radc4=1d0+pjalps(sqm4)/paru(1)
          do 700 i=1,min(16,mdcy(23,3))
            idc=i+mdcy(23,2)-1
            if(mdme(idc,1).lt.0) goto 700
            imdm=0
            if(mdme(idc,1).eq.1.or.mdme(idc,1).eq.2.or.mdme(idc,1).eq.4)
     &      imdm=1
            if(i.le.8) then
              ef=kchg(i,1)/3d0
              af=sign(1d0,ef+0.1d0)
              vf=af-4d0*ef*xwv
            elseif(i.le.16) then
              ef=kchg(i+2,1)/3d0
              af=sign(1d0,ef+0.1d0)
              vf=af-4d0*ef*xwv
            endif
            rm1=pmas(iabs(kfdp(idc,1)),1)**2/sqm4
            if(4d0*rm1.lt.1d0) then
              fcof=1d0
              if(i.le.8) fcof=3d0*radc4
              be34=sqrt(max(0d0,1d0-4d0*rm1))
              if(imdm.eq.1) then
                hfgg=hfgg+fcof*ef**2*(1d0+2d0*rm1)*be34
                hfgz=hfgz+fcof*ef*vf*(1d0+2d0*rm1)*be34
                hfzz=hfzz+fcof*(vf**2*(1d0+2d0*rm1)+
     &          af**2*(1d0-4d0*rm1))*be34
              endif
            endif
  700     continue
C...Propagators: as simulated in PYOFSH and as desired
          hbw4=(1d0/paru(1))*gmmz/((sqm4-sqmz)**2+gmmz**2)
          mint(15)=1
          mint(61)=1
          call pjwidt(23,sqm4,wdtp,wdte)
          hfaem=(paru(108)/paru(2))*(2d0/3d0)
          hfgg=hfgg*hfaem*vint(111)/sqm4
          hfgz=hfgz*hfaem*vint(112)/sqm4
          hfzz=hfzz*hfaem*vint(114)/sqm4
C...Loop over flavours; consider full gamma/Z structure
          do 720 i=mmina,mmaxa
            if(i.eq.0) goto 720
            ei=kchg(iabs(i),1)/3d0
            ai=sign(1d0,ei)
            vi=ai-4d0*ei*xwv
            faczq=ei**2*(ei**2*hfgg+ei*vi*hfgz+
     &      (vi**2+ai**2)*hfzz)/hbw4
            do 710 isde=1,2
              if(isde.eq.1.and.kfac(1,i)*kfac(2,22).eq.0) goto 710
              if(isde.eq.2.and.kfac(1,22)*kfac(2,i).eq.0) goto 710
              nchn=nchn+1
              isig(nchn,isde)=i
              isig(nchn,3-isde)=22
              isig(nchn,3)=1
              sigh(nchn)=faczq*fzqn/max(pmas(iabs(i),1)**2*sqm4,fzqd)
  710       continue
  720     continue
 
        elseif(isub.eq.36) then
C...f + gamma -> f' + W+/-
          fwq=comfac*aem**2/(2d0*xw)*
     &    (sh2+uh2+2d0*sqm4*th)/(sqpth*sqm4-sh*uh)
C...Propagators: as simulated in PYOFSH and as desired
          hbw4=gmmw/((sqm4-sqmw)**2+gmmw**2)
          call pjwidt(24,sqm4,wdtp,wdte)
          gmmwc=sqrt(sqm4)*wdtp(0)
          hbw4c=gmmwc/((sqm4-sqmw)**2+gmmwc**2)
          fwq=fwq*hbw4c/hbw4
          do 740 i=mmina,mmaxa
            if(i.eq.0) goto 740
            ia=iabs(i)
            eia=abs(kchg(iabs(i),1)/3d0)
            facwq=fwq*(eia-sh/(sh+uh))**2
            kchw=isign(1,kchg(ia,1)*isign(1,i))
            widsc=(wdte(0,1)+wdte(0,(5-kchw)/2)+wdte(0,4))/wdtp(0)
            do 730 isde=1,2
              if(isde.eq.1.and.kfac(1,i)*kfac(2,22).eq.0) goto 730
              if(isde.eq.2.and.kfac(1,22)*kfac(2,i).eq.0) goto 730
              nchn=nchn+1
              isig(nchn,isde)=i
              isig(nchn,3-isde)=22
              isig(nchn,3)=1
              sigh(nchn)=facwq*vint(180+i)*widsc
  730       continue
  740     continue
 
        elseif(isub.eq.37) then
C...f + gamma -> f + h0
 
        elseif(isub.eq.38) then
C...f + Z0 -> f + g (q + Z0 -> q + g only)
 
        elseif(isub.eq.39) then
C...f + Z0 -> f + gamma
 
        elseif(isub.eq.40) then
C...f + Z0 -> f + Z0
        endif
 
      elseif(isub.le.50) then
        if(isub.eq.41) then
C...f + Z0 -> f' + W+/-
 
        elseif(isub.eq.42) then
C...f + Z0 -> f + h0
 
        elseif(isub.eq.43) then
C...f + W+/- -> f' + g (q + W+/- -> q' + g only)
 
        elseif(isub.eq.44) then
C...f + W+/- -> f' + gamma
 
        elseif(isub.eq.45) then
C...f + W+/- -> f' + Z0
 
        elseif(isub.eq.46) then
C...f + W+/- -> f' + W+/-
 
        elseif(isub.eq.47) then
C...f + W+/- -> f' + h0
 
        elseif(isub.eq.48) then
C...f + h0 -> f + g (q + h0 -> q + g only)
 
        elseif(isub.eq.49) then
C...f + h0 -> f + gamma
 
        elseif(isub.eq.50) then
C...f + h0 -> f + Z0
        endif
 
      elseif(isub.le.60) then
        if(isub.eq.51) then
C...f + h0 -> f' + W+/-
 
        elseif(isub.eq.52) then
C...f + h0 -> f + h0
 
        elseif(isub.eq.53) then
C...g + g -> f + fbar (g + g -> q + qbar only)
          call pjwidt(21,sh,wdtp,wdte)
          facqq1=comfac*as**2*1d0/6d0*(uh/th-(2d0+mstp(34)*1d0/4d0)*
     &    uh2/sh2)*(wdte(0,1)+wdte(0,2)+wdte(0,3)+wdte(0,4))*faca
          facqq2=comfac*as**2*1d0/6d0*(th/uh-(2d0+mstp(34)*1d0/4d0)*
     &    th2/sh2)*(wdte(0,1)+wdte(0,2)+wdte(0,3)+wdte(0,4))*faca
          if(kfac(1,21)*kfac(2,21).eq.0) goto 750
          nchn=nchn+1
          isig(nchn,1)=21
          isig(nchn,2)=21
          isig(nchn,3)=1
          sigh(nchn)=facqq1
          nchn=nchn+1
          isig(nchn,1)=21
          isig(nchn,2)=21
          isig(nchn,3)=2
          sigh(nchn)=facqq2
  750     continue
 
        elseif(isub.eq.54) then
C...g + gamma -> f + fbar (g + gamma -> q + qbar only)
          call pjwidt(21,sh,wdtp,wdte)
          wdtesu=0d0
          do 760 i=1,min(8,mdcy(21,3))
            ef=kchg(i,1)/3d0
            wdtesu=wdtesu+ef**2*(wdte(i,1)+wdte(i,2)+wdte(i,3)+
     &      wdte(i,4))
  760     continue
          facqq=comfac*aem*as*wdtesu*(th2+uh2)/(th*uh)
          if(kfac(1,21)*kfac(2,22).ne.0) then
            nchn=nchn+1
            isig(nchn,1)=21
            isig(nchn,2)=22
            isig(nchn,3)=1
            sigh(nchn)=facqq
          endif
          if(kfac(1,22)*kfac(2,21).ne.0) then
            nchn=nchn+1
            isig(nchn,1)=22
            isig(nchn,2)=21
            isig(nchn,3)=1
            sigh(nchn)=facqq
          endif
 
        elseif(isub.eq.55) then
C...g + Z -> f + fbar (g + Z -> q + qbar only)
 
        elseif(isub.eq.56) then
C...g + W -> f + f'bar (g + W -> q + q'bar only)
 
        elseif(isub.eq.57) then
C...g + h0 -> f + fbar (g + h0 -> q + qbar only)
 
        elseif(isub.eq.58) then
C...gamma + gamma -> f + fbar
          call pjwidt(22,sh,wdtp,wdte)
          wdtesu=0d0
          do 770 i=1,min(12,mdcy(22,3))
            if(i.le.8) ef= kchg(i,1)/3d0
            if(i.ge.9) ef= kchg(9+2*(i-8),1)/3d0
            wdtesu=wdtesu+ef**2*(wdte(i,1)+wdte(i,2)+wdte(i,3)+
     &      wdte(i,4))
  770     continue
          facff=comfac*aem**2*wdtesu*2d0*(th2+uh2)/(th*uh)
          if(kfac(1,22)*kfac(2,22).ne.0) then
            nchn=nchn+1
            isig(nchn,1)=22
            isig(nchn,2)=22
            isig(nchn,3)=1
            sigh(nchn)=facff
          endif
 
        elseif(isub.eq.59) then
C...gamma + Z0 -> f + fbar
 
        elseif(isub.eq.60) then
C...gamma + W+/- -> f + fbar'
        endif
 
      elseif(isub.le.70) then
        if(isub.eq.61) then
C...gamma + h0 -> f + fbar
 
        elseif(isub.eq.62) then
C...Z0 + Z0 -> f + fbar
 
        elseif(isub.eq.63) then
C...Z0 + W+/- -> f + fbar'
 
        elseif(isub.eq.64) then
C...Z0 + h0 -> f + fbar
 
        elseif(isub.eq.65) then
C...W+ + W- -> f + fbar
 
        elseif(isub.eq.66) then
C...W+/- + h0 -> f + fbar'
 
        elseif(isub.eq.67) then
C...h0 + h0 -> f + fbar
 
        elseif(isub.eq.68) then
C...g + g -> g + g
          facgg1=comfac*as**2*9d0/4d0*(sh2/th2+2d0*sh/th+3d0+2d0*th/sh+
     &    th2/sh2)*faca
          facgg2=comfac*as**2*9d0/4d0*(uh2/sh2+2d0*uh/sh+3d0+2d0*sh/uh+
     &    sh2/uh2)*faca
          facgg3=comfac*as**2*9d0/4d0*(th2/uh2+2d0*th/uh+3d0+2d0*uh/th+
     &    uh2/th2)
          if(kfac(1,21)*kfac(2,21).eq.0) goto 780
          nchn=nchn+1
          isig(nchn,1)=21
          isig(nchn,2)=21
          isig(nchn,3)=1
          sigh(nchn)=0.5d0*facgg1
          nchn=nchn+1
          isig(nchn,1)=21
          isig(nchn,2)=21
          isig(nchn,3)=2
          sigh(nchn)=0.5d0*facgg2
          nchn=nchn+1
          isig(nchn,1)=21
          isig(nchn,2)=21
          isig(nchn,3)=3
          sigh(nchn)=0.5d0*facgg3
  780     continue
 
        elseif(isub.eq.69) then
C...gamma + gamma -> W+ + W-
          sqmwe=max(0.5d0*sqmw,sqrt(sqm3*sqm4))
          fprop=sh2/((sqmwe-th)*(sqmwe-uh))
          facww=comfac*6d0*aem**2*(1d0-fprop*(4d0/3d0+2d0*sqmwe/sh)+
     &    fprop**2*(2d0/3d0+2d0*(sqmwe/sh)**2))*wids(24,1)
          if(kfac(1,22)*kfac(2,22).eq.0) goto 790
          nchn=nchn+1
          isig(nchn,1)=22
          isig(nchn,2)=22
          isig(nchn,3)=1
          sigh(nchn)=facww
  790     continue
 
        elseif(isub.eq.70) then
C...gamma + W+/- -> Z0 + W+/-
          sqmwe=max(0.5d0*sqmw,sqrt(sqm3*sqm4))
          fprop=(th-sqmwe)**2/(-sh*(sqmwe-uh))
          faczw=comfac*6d0*aem**2*(xw1/xw)*
     &    (1d0-fprop*(4d0/3d0+2d0*sqmwe/(th-sqmwe))+
     &    fprop**2*(2d0/3d0+2d0*(sqmwe/(th-sqmwe))**2))*wids(23,2)
          do 810 kchw=1,-1,-2
            do 800 isde=1,2
              if(kfac(isde,22)*kfac(3-isde,24*kchw).eq.0) goto 800
              nchn=nchn+1
              isig(nchn,isde)=22
              isig(nchn,3-isde)=24*kchw
              isig(nchn,3)=1
              sigh(nchn)=faczw*wids(24,(5-kchw)/2)
  800       continue
  810     continue
        endif
 
      elseif(isub.le.80) then
        if(isub.eq.71) then
C...Z0 + Z0 -> Z0 + Z0
          if(sh.le.4.01d0*sqmz) goto 840
 
          if(mstp(46).le.2) then
C...Exact scattering ME:s for on-mass-shell gauge bosons
            be2=1d0-4d0*sqmz/sh
            th=-0.5d0*sh*be2*(1d0-cth)
            uh=-0.5d0*sh*be2*(1d0+cth)
            if(max(th,uh).gt.-1d0) goto 840
            shang=1d0/xw1*sqmw/sqmz*(1d0+be2)**2
            ashre=(sh-sqmh)/((sh-sqmh)**2+gmmh**2)*shang
            ashim=-gmmh/((sh-sqmh)**2+gmmh**2)*shang
            thang=1d0/xw1*sqmw/sqmz*(be2-cth)**2
            athre=(th-sqmh)/((th-sqmh)**2+gmmh**2)*thang
            athim=-gmmh/((th-sqmh)**2+gmmh**2)*thang
            uhang=1d0/xw1*sqmw/sqmz*(be2+cth)**2
            auhre=(uh-sqmh)/((uh-sqmh)**2+gmmh**2)*uhang
            auhim=-gmmh/((uh-sqmh)**2+gmmh**2)*uhang
            faczz=comfac*1d0/(4096d0*paru(1)**2*16d0*xw1**2)*
     &      (aem/xw)**4*(sh/sqmw)**2*(sqmz/sqmw)*sh2
            if(mstp(46).le.0) faczz=faczz*(ashre**2+ashim**2)
            if(mstp(46).eq.1) faczz=faczz*((ashre+athre+auhre)**2+
     &      (ashim+athim+auhim)**2)
            if(mstp(46).eq.2) faczz=0d0
 
          else
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron
            faczz=comfac*(aem/(16d0*paru(1)*xw*xw1))**2*(64d0/9d0)*
     &      abs(a00u+2.*a20u)**2
          endif
          faczz=faczz*wids(23,1)
 
          do 830 i=mmin1,mmax1
            if(i.eq.0.or.kfac(1,i).eq.0) goto 830
            ei=kchg(iabs(i),1)/3d0
            ai=sign(1d0,ei)
            vi=ai-4d0*ei*xwv
            avi=ai**2+vi**2
            do 820 j=mmin2,mmax2
              if(j.eq.0.or.kfac(2,j).eq.0) goto 820
              ej=kchg(iabs(j),1)/3d0
              aj=sign(1d0,ej)
              vj=aj-4d0*ej*xwv
              avj=aj**2+vj**2
              nchn=nchn+1
              isig(nchn,1)=i
              isig(nchn,2)=j
              isig(nchn,3)=1
              sigh(nchn)=0.5d0*faczz*avi*avj
  820       continue
  830     continue
  840     continue
 
        elseif(isub.eq.72) then
C...Z0 + Z0 -> W+ + W-
          if(sh.le.4.01d0*sqmz) goto 870
 
          if(mstp(46).le.2) then
C...Exact scattering ME:s for on-mass-shell gauge bosons
            be2=sqrt((1d0-4d0*sqmw/sh)*(1d0-4d0*sqmz/sh))
            cth2=cth**2
            th=-0.5d0*sh*(1d0-2d0*(sqmw+sqmz)/sh-be2*cth)
            uh=-0.5d0*sh*(1d0-2d0*(sqmw+sqmz)/sh+be2*cth)
            if(max(th,uh).gt.-1d0) goto 870
            shang=4d0*sqrt(sqmw/(sqmz*xw1))*(1d0-2d0*sqmw/sh)*
     &      (1d0-2d0*sqmz/sh)
            ashre=(sh-sqmh)/((sh-sqmh)**2+gmmh**2)*shang
            ashim=-gmmh/((sh-sqmh)**2+gmmh**2)*shang
            atwre=xw1/sqmz*sh/(th-sqmw)*((cth-be2)**2*(3d0/2d0+be2/2d0*
     &      cth-(sqmw+sqmz)/sh+(sqmw-sqmz)**2/(sh*sqmw))+4d0*
     &      ((sqmw+sqmz)/sh*(1d0-3d0*cth2)+8d0*sqmw*sqmz/sh2*
     &      (2d0*cth2-1d0)+4d0*(sqmw**2+sqmz**2)/sh2*cth2+
     &      2d0*(sqmw+sqmz)/sh*be2*cth))
            atwim=0d0
            auwre=xw1/sqmz*sh/(uh-sqmw)*((cth+be2)**2*(3d0/2d0-be2/2d0*
     &      cth-(sqmw+sqmz)/sh+(sqmw-sqmz)**2/(sh*sqmw))+4d0*
     &      ((sqmw+sqmz)/sh*(1d0-3d0*cth2)+8d0*sqmw*sqmz/sh2*
     &      (2d0*cth2-1d0)+4d0*(sqmw**2+sqmz**2)/sh2*cth2-
     &      2d0*(sqmw+sqmz)/sh*be2*cth))
            auwim=0d0
            a4re=2d0*xw1/sqmz*(3d0-cth2-4d0*(sqmw+sqmz)/sh)
            a4im=0d0
            facww=comfac*1d0/(4096d0*paru(1)**2*16d0*xw1**2)*
     &      (aem/xw)**4*(sh/sqmw)**2*(sqmz/sqmw)*sh2
            if(mstp(46).le.0) facww=facww*(ashre**2+ashim**2)
            if(mstp(46).eq.1) facww=facww*((ashre+atwre+auwre+a4re)**2+
     &      (ashim+atwim+auwim+a4im)**2)
            if(mstp(46).eq.2) facww=facww*((atwre+auwre+a4re)**2+
     &      (atwim+auwim+a4im)**2)
 
          else
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron
            facww=comfac*(aem/(16d0*paru(1)*xw*xw1))**2*(64d0/9d0)*
     &      abs(a00u-a20u)**2
          endif
          facww=facww*wids(24,1)
 
          do 860 i=mmin1,mmax1
            if(i.eq.0.or.kfac(1,i).eq.0) goto 860
            ei=kchg(iabs(i),1)/3d0
            ai=sign(1d0,ei)
            vi=ai-4d0*ei*xwv
            avi=ai**2+vi**2
            do 850 j=mmin2,mmax2
              if(j.eq.0.or.kfac(2,j).eq.0) goto 850
              ej=kchg(iabs(j),1)/3d0
              aj=sign(1d0,ej)
              vj=aj-4d0*ej*xwv
              avj=aj**2+vj**2
              nchn=nchn+1
              isig(nchn,1)=i
              isig(nchn,2)=j
              isig(nchn,3)=1
              sigh(nchn)=facww*avi*avj
  850       continue
  860     continue
  870     continue
 
        elseif(isub.eq.73) then
C...Z0 + W+/- -> Z0 + W+/-
          if(sh.le.2d0*sqmz+2d0*sqmw) goto 900
 
          if(mstp(46).le.2) then
C...Exact scattering ME:s for on-mass-shell gauge bosons
            be2=1d0-2d0*(sqmz+sqmw)/sh+((sqmz-sqmw)/sh)**2
            ep1=1d0-(sqmz-sqmw)/sh
            ep2=1d0+(sqmz-sqmw)/sh
            th=-0.5d0*sh*be2*(1d0-cth)
            uh=(sqmz-sqmw)**2/sh-0.5d0*sh*be2*(1d0+cth)
            if(max(th,uh).gt.-1d0) goto 900
            thang=(be2-ep1*cth)*(be2-ep2*cth)
            athre=(th-sqmh)/((th-sqmh)**2+gmmh**2)*thang
            athim=-gmmh/((th-sqmh)**2+gmmh**2)*thang
            aswre=-xw1/sqmz*sh/(sh-sqmw)*(-be2*(ep1+ep2)**4*cth+
     &      1d0/4d0*(be2+ep1*ep2)**2*((ep1-ep2)**2-4d0*be2*cth)+
     &      2d0*be2*(be2+ep1*ep2)*(ep1+ep2)**2*cth-
     &      1d0/16d0*sh/sqmw*(ep1**2-ep2**2)**2*(be2+ep1*ep2)**2)
            aswim=0d0
            auwre=xw1/sqmz*sh/(uh-sqmw)*(-be2*(ep2+ep1*cth)*
     &      (ep1+ep2*cth)*(be2+ep1*ep2)+be2*(ep2+ep1*cth)*
     &      (be2+ep1*ep2*cth)*(2d0*ep2-ep2*cth+ep1)-
     &      be2*(ep2+ep1*cth)**2*(be2-ep2**2*cth)-1d0/8d0*
     &      (be2+ep1*ep2*cth)**2*((ep1+ep2)**2+2d0*be2*(1d0-cth))+
     &      1d0/32d0*sh/sqmw*(be2+ep1*ep2*cth)**2*
     &      (ep1**2-ep2**2)**2-be2*(ep1+ep2*cth)*(ep2+ep1*cth)*
     &      (be2+ep1*ep2)+be2*(ep1+ep2*cth)*(be2+ep1*ep2*cth)*
     &      (2d0*ep1-ep1*cth+ep2)-be2*(ep1+ep2*cth)**2*
     &      (be2-ep1**2*cth)-1d0/8d0*(be2+ep1*ep2*cth)**2*
     &      ((ep1+ep2)**2+2d0*be2*(1d0-cth))+1d0/32d0*sh/sqmw*
     &      (be2+ep1*ep2*cth)**2*(ep1**2-ep2**2)**2)
            auwim=0d0
            a4re=xw1/sqmz*(ep1**2*ep2**2*(cth**2-1d0)-
     &      2d0*be2*(ep1**2+ep2**2+ep1*ep2)*cth-2d0*be2*ep1*ep2)
            a4im=0d0
            faczw=comfac*1d0/(4096d0*paru(1)**2*4d0*xw1)*(aem/xw)**4*
     &      (sh/sqmw)**2*sqrt(sqmz/sqmw)*sh2
            if(mstp(46).le.0) faczw=0d0
            if(mstp(46).eq.1) faczw=faczw*((athre+aswre+auwre+a4re)**2+
     &      (athim+aswim+auwim+a4im)**2)
            if(mstp(46).eq.2) faczw=faczw*((aswre+auwre+a4re)**2+
     &      (aswim+auwim+a4im)**2)
 
          else
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron
            faczw=comfac*aem**2/(64d0*paru(1)**2*xw**2*xw1)*16d0*
     &      abs(a20u+3.*a11u*sngl(cth))**2
          endif
          faczw=faczw*wids(23,2)
 
          do 890 i=mmin1,mmax1
            if(i.eq.0.or.kfac(1,i).eq.0) goto 890
            ei=kchg(iabs(i),1)/3d0
            ai=sign(1d0,ei)
            vi=ai-4d0*ei*xwv
            avi=ai**2+vi**2
            kchwi=isign(1,kchg(iabs(i),1)*isign(1,i))
            do 880 j=mmin2,mmax2
              if(j.eq.0.or.kfac(2,j).eq.0) goto 880
              ej=kchg(iabs(j),1)/3d0
              aj=sign(1d0,ej)
              vj=ai-4d0*ej*xwv
              avj=aj**2+vj**2
              kchwj=isign(1,kchg(iabs(j),1)*isign(1,j))
              nchn=nchn+1
              isig(nchn,1)=i
              isig(nchn,2)=j
              isig(nchn,3)=1
              sigh(nchn)=faczw*avi*vint(180+j)*wids(24,(5-kchwj)/2)
              nchn=nchn+1
              isig(nchn,1)=i
              isig(nchn,2)=j
              isig(nchn,3)=2
              sigh(nchn)=faczw*vint(180+i)*wids(24,(5-kchwi)/2)*avj
  880       continue
  890     continue
  900     continue
 
        elseif(isub.eq.75) then
C...W+ + W- -> gamma + gamma
 
        elseif(isub.eq.76) then
C...W+ + W- -> Z0 + Z0
          if(sh.le.4.01d0*sqmz) goto 930
 
          if(mstp(46).le.2) then
C...Exact scattering ME:s for on-mass-shell gauge bosons
            be2=sqrt((1d0-4d0*sqmw/sh)*(1d0-4d0*sqmz/sh))
            cth2=cth**2
            th=-0.5d0*sh*(1d0-2d0*(sqmw+sqmz)/sh-be2*cth)
            uh=-0.5d0*sh*(1d0-2d0*(sqmw+sqmz)/sh+be2*cth)
            if(max(th,uh).gt.-1d0) goto 930
            shang=4d0*sqrt(sqmw/(sqmz*xw1))*(1d0-2d0*sqmw/sh)*
     &      (1d0-2d0*sqmz/sh)
            ashre=(sh-sqmh)/((sh-sqmh)**2+gmmh**2)*shang
            ashim=-gmmh/((sh-sqmh)**2+gmmh**2)*shang
            atwre=xw1/sqmz*sh/(th-sqmw)*((cth-be2)**2*(3d0/2d0+be2/2d0*
     &      cth-(sqmw+sqmz)/sh+(sqmw-sqmz)**2/(sh*sqmw))+4d0*
     &      ((sqmw+sqmz)/sh*(1d0-3d0*cth2)+8d0*sqmw*sqmz/sh2*
     &      (2d0*cth2-1d0)+4d0*(sqmw**2+sqmz**2)/sh2*cth2+
     &      2d0*(sqmw+sqmz)/sh*be2*cth))
            atwim=0d0
            auwre=xw1/sqmz*sh/(uh-sqmw)*((cth+be2)**2*(3d0/2d0-be2/2d0*
     &      cth-(sqmw+sqmz)/sh+(sqmw-sqmz)**2/(sh*sqmw))+4d0*
     &      ((sqmw+sqmz)/sh*(1d0-3d0*cth2)+8d0*sqmw*sqmz/sh2*
     &      (2d0*cth2-1d0)+4d0*(sqmw**2+sqmz**2)/sh2*cth2-
     &      2d0*(sqmw+sqmz)/sh*be2*cth))
            auwim=0d0
            a4re=2d0*xw1/sqmz*(3d0-cth2-4d0*(sqmw+sqmz)/sh)
            a4im=0d0
            faczz=comfac*1d0/(4096d0*paru(1)**2)*(aem/xw)**4*
     &      (sh/sqmw)**2*sh2
            if(mstp(46).le.0) faczz=faczz*(ashre**2+ashim**2)
            if(mstp(46).eq.1) faczz=faczz*((ashre+atwre+auwre+a4re)**2+
     &      (ashim+atwim+auwim+a4im)**2)
            if(mstp(46).eq.2) faczz=faczz*((atwre+auwre+a4re)**2+
     &      (atwim+auwim+a4im)**2)
 
          else
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron
            faczz=comfac*(aem/(4d0*paru(1)*xw))**2*(64d0/9d0)*
     &      abs(a00u-a20u)**2
          endif
          faczz=faczz*wids(23,1)
 
          do 920 i=mmin1,mmax1
            if(i.eq.0.or.kfac(1,i).eq.0) goto 920
            ei=sign(1d0,dble(i))*kchg(iabs(i),1)
            do 910 j=mmin2,mmax2
              if(j.eq.0.or.kfac(2,j).eq.0) goto 910
              ej=sign(1d0,dble(j))*kchg(iabs(j),1)
              if(ei*ej.gt.0d0) goto 910
              nchn=nchn+1
              isig(nchn,1)=i
              isig(nchn,2)=j
              isig(nchn,3)=1
              sigh(nchn)=0.5d0*faczz*vint(180+i)*vint(180+j)
  910       continue
  920     continue
  930     continue
 
        elseif(isub.eq.77) then
C...W+/- + W+/- -> W+/- + W+/-
          if(sh.le.4.01d0*sqmw) goto 960
 
          if(mstp(46).le.2) then
C...Exact scattering ME:s for on-mass-shell gauge bosons
            be2=1d0-4d0*sqmw/sh
            be4=be2**2
            cth2=cth**2
            cth3=cth**3
            th=-0.5d0*sh*be2*(1d0-cth)
            uh=-0.5d0*sh*be2*(1d0+cth)
            if(max(th,uh).gt.-1d0) goto 960
            shang=(1d0+be2)**2
            ashre=(sh-sqmh)/((sh-sqmh)**2+gmmh**2)*shang
            ashim=-gmmh/((sh-sqmh)**2+gmmh**2)*shang
            thang=(be2-cth)**2
            athre=(th-sqmh)/((th-sqmh)**2+gmmh**2)*thang
            athim=-gmmh/((th-sqmh)**2+gmmh**2)*thang
            uhang=(be2+cth)**2
            auhre=(uh-sqmh)/((uh-sqmh)**2+gmmh**2)*uhang
            auhim=-gmmh/((uh-sqmh)**2+gmmh**2)*uhang
            sgzang=1d0/sqmw*be2*(3d0-be2)**2*cth
            asgre=xw*sgzang
            asgim=0d0
            aszre=xw1*sh/(sh-sqmz)*sgzang
            aszim=0d0
            tgzang=1d0/sqmw*(be2*(4d0-2d0*be2+be4)+be2*(4d0-10d0*be2+
     &      be4)*cth+(2d0-11d0*be2+10d0*be4)*cth2+be2*cth3)
            atgre=0.5d0*xw*sh/th*tgzang
            atgim=0d0
            atzre=0.5d0*xw1*sh/(th-sqmz)*tgzang
            atzim=0d0
            ugzang=1d0/sqmw*(be2*(4d0-2d0*be2+be4)-be2*(4d0-10d0*be2+
     &      be4)*cth+(2d0-11d0*be2+10d0*be4)*cth2-be2*cth3)
            augre=0.5d0*xw*sh/uh*ugzang
            augim=0d0
            auzre=0.5d0*xw1*sh/(uh-sqmz)*ugzang
            auzim=0d0
            a4are=1d0/sqmw*(1d0+2d0*be2-6d0*be2*cth-cth2)
            a4aim=0d0
            a4sre=2d0/sqmw*(1d0+2d0*be2-cth2)
            a4sim=0d0
            fww=comfac*1d0/(4096d0*paru(1)**2)*(aem/xw)**4*
     &      (sh/sqmw)**2*sh2
            if(mstp(46).le.0) then
              awware=ashre
              awwaim=ashim
              awwsre=0d0
              awwsim=0d0
            elseif(mstp(46).eq.1) then
              awware=ashre+athre+asgre+aszre+atgre+atzre+a4are
              awwaim=ashim+athim+asgim+aszim+atgim+atzim+a4aim
              awwsre=-athre-auhre+atgre+atzre+augre+auzre+a4sre
              awwsim=-athim-auhim+atgim+atzim+augim+auzim+a4sim
            else
              awware=asgre+aszre+atgre+atzre+a4are
              awwaim=asgim+aszim+atgim+atzim+a4aim
              awwsre=atgre+atzre+augre+auzre+a4sre
              awwsim=atgim+atzim+augim+auzim+a4sim
            endif
            awwa2=awware**2+awwaim**2
            awws2=awwsre**2+awwsim**2
 
          else
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron
            fwwa=comfac*(aem/(4d0*paru(1)*xw))**2*(64d0/9d0)*
     &      abs(a00u+0.5*a20u+4.5*a11u*sngl(cth))**2
            fwws=comfac*(aem/(4d0*paru(1)*xw))**2*64d0*abs(a20u)**2
          endif
 
          do 950 i=mmin1,mmax1
            if(i.eq.0.or.kfac(1,i).eq.0) goto 950
            ei=sign(1d0,dble(i))*kchg(iabs(i),1)
            do 940 j=mmin2,mmax2
              if(j.eq.0.or.kfac(2,j).eq.0) goto 940
              ej=sign(1d0,dble(j))*kchg(iabs(j),1)
              if(ei*ej.lt.0d0) then
C...W+W-
                if(mstp(45).eq.1) goto 940
                if(mstp(46).le.2) facww=fww*awwa2*wids(24,1)
                if(mstp(46).ge.3) facww=fwwa*wids(24,1)
              else
C...W+W+/W-W-
                if(mstp(45).eq.2) goto 940
                if(mstp(46).le.2) facww=fww*awws2
                if(mstp(46).ge.3) facww=fwws
                if(ei.gt.0d0) facww=facww*wids(24,4)
                if(ei.lt.0d0) facww=facww*wids(24,5)
              endif
              nchn=nchn+1
              isig(nchn,1)=i
              isig(nchn,2)=j
              isig(nchn,3)=1
              sigh(nchn)=facww*vint(180+i)*vint(180+j)
              if(ei*ej.gt.0d0) sigh(nchn)=0.5d0*sigh(nchn)
  940       continue
  950     continue
  960     continue
 
        elseif(isub.eq.78) then
C...W+/- + h0 -> W+/- + h0
 
        elseif(isub.eq.79) then
C...h0 + h0 -> h0 + h0
 
        elseif(isub.eq.80) then
C...q + gamma -> q' + pi+/-
          fqpi=comfac*(2d0*aem/9d0)*(-sh/th)*(1d0/sh2+1d0/th2)
          assh=pjalps(max(0.5d0,0.5d0*sh))
          q2fpsh=0.55d0/log(max(2d0,2d0*sh))
          delsh=uh*sqrt(assh*q2fpsh)
          asuh=pjalps(max(0.5d0,-0.5d0*uh))
          q2fpuh=0.55d0/log(max(2d0,-2d0*uh))
          deluh=sh*sqrt(asuh*q2fpuh)
          do 980 i=max(-2,mmina),min(2,mmaxa)
            if(i.eq.0) goto 980
            ei=kchg(iabs(i),1)/3d0
            ej=sign(1d0-abs(ei),ei)
            do 970 isde=1,2
              if(isde.eq.1.and.kfac(1,i)*kfac(2,22).eq.0) goto 970
              if(isde.eq.2.and.kfac(1,22)*kfac(2,i).eq.0) goto 970
              nchn=nchn+1
              isig(nchn,isde)=i
              isig(nchn,3-isde)=22
              isig(nchn,3)=1
              sigh(nchn)=fqpi*(ei*delsh+ej*deluh)**2
  970       continue
  980     continue
 
        endif
 
C...C: 2 -> 2, tree diagrams with masses
 
      elseif(isub.le.90) then
        if(isub.eq.81) then
C...q + qbar -> Q + Qbar
          facqqb=comfac*as**2*4d0/9d0*(((th-sqm3)**2+
     &    (uh-sqm3)**2)/sh2+2d0*sqm3/sh)
          if(mstp(35).ge.1) facqqb=facqqb*pjhfth(sh,sqm3,0d0)
          wid2=1d0
          if(mint(55).eq.6) wid2=wids(6,1)
          if(mint(55).eq.7.or.mint(55).eq.8) wid2=wids(mint(55),1)
          facqqb=facqqb*wid2
          do 990 i=mmina,mmaxa
            if(i.eq.0.or.iabs(i).gt.mstp(58).or.
     &      kfac(1,i)*kfac(2,-i).eq.0) goto 990
            nchn=nchn+1
            isig(nchn,1)=i
            isig(nchn,2)=-i
            isig(nchn,3)=1
            sigh(nchn)=facqqb
  990     continue
 
        elseif(isub.eq.82) then
C...g + g -> Q + Qbar
          if(mstp(34).eq.0) then
            facqq1=comfac*faca*as**2*(1d0/6d0)*((uh-sqm3)/(th-sqm3)-
     &      2d0*(uh-sqm3)**2/sh2+4d0*(sqm3/sh)*(th*uh-sqm3**2)/
     &      (th-sqm3)**2)
            facqq2=comfac*faca*as**2*(1d0/6d0)*((th-sqm3)/(uh-sqm3)-
     &      2d0*(th-sqm3)**2/sh2+4d0*(sqm3/sh)*(th*uh-sqm3**2)/
     &      (uh-sqm3)**2)
          else
            facqq1=comfac*faca*as**2*(1d0/6d0)*((uh-sqm3)/(th-sqm3)-
     &      2.25d0*(uh-sqm3)**2/sh2+4.5d0*(sqm3/sh)*(th*uh-sqm3**2)/
     &      (th-sqm3)**2+0.5d0*sqm3*th/(th-sqm3)**2-sqm3**2/
     &      (sh*(th-sqm3)))
            facqq2=comfac*faca*as**2*(1d0/6d0)*((th-sqm3)/(uh-sqm3)-
     &      2.25d0*(th-sqm3)**2/sh2+4.5d0*(sqm3/sh)*(th*uh-sqm3**2)/
     &      (uh-sqm3)**2+0.5d0*sqm3*uh/(uh-sqm3)**2-sqm3**2/
     &      (sh*(uh-sqm3)))
          endif
          if(mstp(35).ge.1) then
            fatre=pjhfth(sh,sqm3,2d0/7d0)
            facqq1=facqq1*fatre
            facqq2=facqq2*fatre
          endif
          wid2=1d0
          if(mint(55).eq.6) wid2=wids(6,1)
          if(mint(55).eq.7.or.mint(55).eq.8) wid2=wids(mint(55),1)
          facqq1=facqq1*wid2
          facqq2=facqq2*wid2
          if(kfac(1,21)*kfac(2,21).eq.0) goto 1000
          nchn=nchn+1
          isig(nchn,1)=21
          isig(nchn,2)=21
          isig(nchn,3)=1
          sigh(nchn)=facqq1
          nchn=nchn+1
          isig(nchn,1)=21
          isig(nchn,2)=21
          isig(nchn,3)=2
          sigh(nchn)=facqq2
 1000     continue
 
        elseif(isub.eq.83) then
C...f + q -> f' + Q
          facqqs=comfac*(0.5d0*aem/xw)**2*sh*(sh-sqm3)/(sqmw-th)**2
          facqqu=comfac*(0.5d0*aem/xw)**2*uh*(uh-sqm3)/(sqmw-th)**2
          do 1020 i=mmin1,mmax1
            if(i.eq.0.or.kfac(1,i).eq.0) goto 1020
            do 1010 j=mmin2,mmax2
              if(j.eq.0.or.kfac(2,j).eq.0) goto 1010
              if(i*j.gt.0.and.mod(iabs(i+j),2).eq.0) goto 1010
              if(i*j.lt.0.and.mod(iabs(i+j),2).eq.1) goto 1010
              if(iabs(i).lt.mint(55).and.mod(iabs(i+mint(55)),2).eq.1)
     &        then
                nchn=nchn+1
                isig(nchn,1)=i
                isig(nchn,2)=j
                isig(nchn,3)=1
                if(mod(mint(55),2).eq.0) facckm=vckm(mint(55)/2,
     &          (iabs(i)+1)/2)*vint(180+j)
                if(mod(mint(55),2).eq.1) facckm=vckm(iabs(i)/2,
     &          (mint(55)+1)/2)*vint(180+j)
                wid2=1d0
                if(i.gt.0) then
                  if(mint(55).eq.6) wid2=wids(6,2)
                  if(mint(55).eq.7.or.mint(55).eq.8) wid2=
     &            wids(mint(55),2)
                else
                  if(mint(55).eq.6) wid2=wids(6,3)
                  if(mint(55).eq.7.or.mint(55).eq.8) wid2=
     &            wids(mint(55),3)
                endif
                if(i*j.gt.0) sigh(nchn)=facqqs*facckm*wid2
                if(i*j.lt.0) sigh(nchn)=facqqu*facckm*wid2
              endif
              if(iabs(j).lt.mint(55).and.mod(iabs(j+mint(55)),2).eq.1)
     &        then
                nchn=nchn+1
                isig(nchn,1)=i
                isig(nchn,2)=j
                isig(nchn,3)=2
                if(mod(mint(55),2).eq.0) facckm=vckm(mint(55)/2,
     &          (iabs(j)+1)/2)*vint(180+i)
                if(mod(mint(55),2).eq.1) facckm=vckm(iabs(j)/2,
     &          (mint(55)+1)/2)*vint(180+i)
                if(j.gt.0) then
                  if(mint(55).eq.6) wid2=wids(6,2)
                  if(mint(55).eq.7.or.mint(55).eq.8) wid2=
     &            wids(mint(55),2)
                else
                  if(mint(55).eq.6) wid2=wids(6,3)
                  if(mint(55).eq.7.or.mint(55).eq.8) wid2=
     &            wids(mint(55),3)
                endif
                if(i*j.gt.0) sigh(nchn)=facqqs*facckm*wid2
                if(i*j.lt.0) sigh(nchn)=facqqu*facckm*wid2
              endif
 1010       continue
 1020     continue
 
        elseif(isub.eq.84) then
C...g + gamma -> Q + Qbar
          fmtu=sqm3/(sqm3-th)+sqm3/(sqm3-uh)
          facqq=comfac*as*aem*(kchg(iabs(mint(55)),1)/3d0)**2*
     &    ((sqm3-th)/(sqm3-uh)+(sqm3-uh)/(sqm3-th)+4d0*fmtu*(1d0-fmtu))
          if(mstp(35).ge.1) facqq=facqq*pjhfth(sh,sqm3,0d0)
          wid2=1d0
          if(mint(55).eq.6) wid2=wids(6,1)
          if(mint(55).eq.7.or.mint(55).eq.8) wid2=wids(mint(55),1)
          facqq=facqq*wid2
          if(kfac(1,21)*kfac(2,22).ne.0) then
            nchn=nchn+1
            isig(nchn,1)=21
            isig(nchn,2)=22
            isig(nchn,3)=1
            sigh(nchn)=facqq
          endif
          if(kfac(1,22)*kfac(2,21).ne.0) then
            nchn=nchn+1
            isig(nchn,1)=22
            isig(nchn,2)=21
            isig(nchn,3)=1
            sigh(nchn)=facqq
          endif
 
        elseif(isub.eq.85) then
C...gamma + gamma -> F + Fbar (heavy fermion, quark or lepton)
          fmtu=sqm3/(sqm3-th)+sqm3/(sqm3-uh)
          facff=comfac*aem**2*(kchg(iabs(mint(56)),1)/3d0)**4*2d0*
     &    ((sqm3-th)/(sqm3-uh)+(sqm3-uh)/(sqm3-th)+4d0*fmtu*(1d0-fmtu))
          if(iabs(mint(56)).lt.10) facff=3d0*facff
          if(iabs(mint(56)).lt.10.and.mstp(35).ge.1)
     &    facff=facff*pjhfth(sh,sqm3,1d0)
          wid2=1d0
          if(mint(56).eq.6) wid2=wids(6,1)
          if(mint(56).eq.7.or.mint(56).eq.8) wid2=wids(mint(56),1)
          if(mint(56).eq.17) wid2=wids(17,1)
          facff=facff*wid2
          if(kfac(1,22)*kfac(2,22).ne.0) then
            nchn=nchn+1
            isig(nchn,1)=22
            isig(nchn,2)=22
            isig(nchn,3)=1
            sigh(nchn)=facff
          endif
 
        elseif(isub.eq.86) then
C...g + g -> J/Psi + g
          facqqg=comfac*as**3*(5d0/9d0)*parp(38)*sqrt(sqm3)*
     &    (((sh*(sh-sqm3))**2+(th*(th-sqm3))**2+(uh*(uh-sqm3))**2)/
     &    ((th-sqm3)*(uh-sqm3))**2)/(sh-sqm3)**2
          if(kfac(1,21)*kfac(2,21).ne.0) then
            nchn=nchn+1
            isig(nchn,1)=21
            isig(nchn,2)=21
            isig(nchn,3)=1
            sigh(nchn)=facqqg
          endif
 
        elseif(isub.eq.87) then
C...g + g -> chi_0c + g
          pgtw=(sh*th+th*uh+uh*sh)/sh2
          qgtw=(sh*th*uh)/sh**3
          rgtw=sqm3/sh
          facqqg=comfac*as**3*4d0*(parp(39)/sqrt(sqm3))*(1d0/sh)*
     &    (9d0*rgtw**2*pgtw**4*(rgtw**4-2d0*rgtw**2*pgtw+pgtw**2)-
     &    6d0*rgtw*pgtw**3*qgtw*(2d0*rgtw**4-5d0*rgtw**2*pgtw+pgtw**2)-
     &    pgtw**2*qgtw**2*(rgtw**4+2d0*rgtw**2*pgtw-pgtw**2)+
     &    2d0*rgtw*pgtw*qgtw**3*(rgtw**2-pgtw)+6d0*rgtw**2*qgtw**4)/
     &    (qgtw*(qgtw-rgtw*pgtw)**4)
          if(kfac(1,21)*kfac(2,21).ne.0) then
            nchn=nchn+1
            isig(nchn,1)=21
            isig(nchn,2)=21
            isig(nchn,3)=1
            sigh(nchn)=facqqg
          endif
 
        elseif(isub.eq.88) then
C...g + g -> chi_1c + g
          pgtw=(sh*th+th*uh+uh*sh)/sh2
          qgtw=(sh*th*uh)/sh**3
          rgtw=sqm3/sh
          facqqg=comfac*as**3*12d0*(parp(39)/sqrt(sqm3))*(1d0/sh)*
     &    pgtw**2*(rgtw*pgtw**2*(rgtw**2-4d0*pgtw)+2d0*qgtw*(-rgtw**4+
     &    5d0*rgtw**2*pgtw+pgtw**2)-15d0*rgtw*qgtw**2)/
     &    (qgtw-rgtw*pgtw)**4
          if(kfac(1,21)*kfac(2,21).ne.0) then
            nchn=nchn+1
            isig(nchn,1)=21
            isig(nchn,2)=21
            isig(nchn,3)=1
            sigh(nchn)=facqqg
          endif
 
        elseif(isub.eq.89) then
C...g + g -> chi_2c + g
          pgtw=(sh*th+th*uh+uh*sh)/sh2
          qgtw=(sh*th*uh)/sh**3
          rgtw=sqm3/sh
          facqqg=comfac*as**3*4d0*(parp(39)/sqrt(sqm3))*(1d0/sh)*
     &    (12d0*rgtw**2*pgtw**4*(rgtw**4-2d0*rgtw**2*pgtw+pgtw**2)-
     &    3d0*rgtw*pgtw**3*qgtw*(8d0*rgtw**4-rgtw**2*pgtw+4d0*pgtw**2)+
     &    2d0*pgtw**2*qgtw**2*(-7d0*rgtw**4+43d0*rgtw**2*pgtw+pgtw**2)+
     &    rgtw*pgtw*qgtw**3*(16d0*rgtw**2-61d0*pgtw)+12d0*rgtw**2*
     &    qgtw**4)/(qgtw*(qgtw-rgtw*pgtw)**4)
          if(kfac(1,21)*kfac(2,21).ne.0) then
            nchn=nchn+1
            isig(nchn,1)=21
            isig(nchn,2)=21
            isig(nchn,3)=1
            sigh(nchn)=facqqg
          endif
        endif
 
C...D: Mimimum bias processes
 
      elseif(isub.le.100) then
        if(isub.eq.91) then
C...Elastic scattering
          sigs=sigt(0,0,1)
 
        elseif(isub.eq.92) then
C...Single diffractive scattering (first side, i.e. XB)
          sigs=sigt(0,0,2)
 
        elseif(isub.eq.93) then
C...Single diffractive scattering (second side, i.e. AX)
          sigs=sigt(0,0,3)
 
        elseif(isub.eq.94) then
C...Double diffractive scattering
          sigs=sigt(0,0,4)
 
        elseif(isub.eq.95) then
C...Low-pT scattering
          sigs=sigt(0,0,5)
 
        elseif(isub.eq.96) then
C...Multiple interactions: sum of QCD processes
          call pjwidt(21,sh,wdtp,wdte)
 
C...q + q' -> q + q'
          facqq1=comfac*as**2*4d0/9d0*(sh2+uh2)/th2
          facqqb=comfac*as**2*4d0/9d0*((sh2+uh2)/th2*faca-
     &    mstp(34)*2d0/3d0*uh2/(sh*th))
          facqq2=comfac*as**2*4d0/9d0*((sh2+th2)/uh2-
     &    mstp(34)*2d0/3d0*sh2/(th*uh))
          do 1040 i=-3,3
            if(i.eq.0) goto 1040
            do 1030 j=-3,3
              if(j.eq.0) goto 1030
              nchn=nchn+1
              isig(nchn,1)=i
              isig(nchn,2)=j
              isig(nchn,3)=111
              sigh(nchn)=facqq1
              if(i.eq.-j) sigh(nchn)=facqqb
              if(i.eq.j) then
                sigh(nchn)=0.5d0*sigh(nchn)
                nchn=nchn+1
                isig(nchn,1)=i
                isig(nchn,2)=j
                isig(nchn,3)=112
                sigh(nchn)=0.5d0*facqq2
              endif
 1030       continue
 1040     continue
 
C...q + qbar -> q' + qbar' or g + g
          facqqb=comfac*as**2*4d0/9d0*(th2+uh2)/sh2*
     &    (wdte(0,1)+wdte(0,2)+wdte(0,3)+wdte(0,4))
          facgg1=comfac*as**2*32d0/27d0*(uh/th-(2d0+mstp(34)*1d0/4d0)*
     &    uh2/sh2)
          facgg2=comfac*as**2*32d0/27d0*(th/uh-(2d0+mstp(34)*1d0/4d0)*
     &    th2/sh2)
          do 1050 i=-3,3
            if(i.eq.0) goto 1050
            nchn=nchn+1
            isig(nchn,1)=i
            isig(nchn,2)=-i
            isig(nchn,3)=121
            sigh(nchn)=facqqb
            nchn=nchn+1
            isig(nchn,1)=i
            isig(nchn,2)=-i
            isig(nchn,3)=131
            sigh(nchn)=0.5d0*facgg1
            nchn=nchn+1
            isig(nchn,1)=i
            isig(nchn,2)=-i
            isig(nchn,3)=132
            sigh(nchn)=0.5d0*facgg2
 1050     continue
 
C...q + g -> q + g
          facqg1=comfac*as**2*4d0/9d0*((2d0+mstp(34)*1d0/4d0)*uh2/th2-
     &    uh/sh)*faca
          facqg2=comfac*as**2*4d0/9d0*((2d0+mstp(34)*1d0/4d0)*sh2/th2-
     &    sh/uh)
          do 1070 i=-3,3
            if(i.eq.0) goto 1070
            do 1060 isde=1,2
              nchn=nchn+1
              isig(nchn,isde)=i
              isig(nchn,3-isde)=21
              isig(nchn,3)=281
              sigh(nchn)=facqg1
              nchn=nchn+1
              isig(nchn,isde)=i
              isig(nchn,3-isde)=21
              isig(nchn,3)=282
              sigh(nchn)=facqg2
 1060       continue
 1070     continue
 
C...g + g -> q + qbar or g + g
          facqq1=comfac*as**2*1d0/6d0*(uh/th-(2d0+mstp(34)*1d0/4d0)*
     &    uh2/sh2)*(wdte(0,1)+wdte(0,2)+wdte(0,3)+wdte(0,4))*faca
          facqq2=comfac*as**2*1d0/6d0*(th/uh-(2d0+mstp(34)*1d0/4d0)*
     &    th2/sh2)*(wdte(0,1)+wdte(0,2)+wdte(0,3)+wdte(0,4))*faca
          facgg1=comfac*as**2*9d0/4d0*(sh2/th2+2d0*sh/th+3d0+
     &    2d0*th/sh+th2/sh2)*faca
          facgg2=comfac*as**2*9d0/4d0*(uh2/sh2+2d0*uh/sh+3d0+
     &    2d0*sh/uh+sh2/uh2)*faca
          facgg3=comfac*as**2*9d0/4d0*(th2/uh2+2d0*th/uh+3+
     &    2d0*uh/th+uh2/th2)
          nchn=nchn+1
          isig(nchn,1)=21
          isig(nchn,2)=21
          isig(nchn,3)=531
          sigh(nchn)=facqq1
          nchn=nchn+1
          isig(nchn,1)=21
          isig(nchn,2)=21
          isig(nchn,3)=532
          sigh(nchn)=facqq2
          nchn=nchn+1
          isig(nchn,1)=21
          isig(nchn,2)=21
          isig(nchn,3)=681
          sigh(nchn)=0.5d0*facgg1
          nchn=nchn+1
          isig(nchn,1)=21
          isig(nchn,2)=21
          isig(nchn,3)=682
          sigh(nchn)=0.5d0*facgg2
          nchn=nchn+1
          isig(nchn,1)=21
          isig(nchn,2)=21
          isig(nchn,3)=683
          sigh(nchn)=0.5d0*facgg3
        endif
 
C...E: 2 -> 1, loop diagrams
      elseif(isub.le.110) then
 
      elseif(isub.le.120) then
 
C...G: 2 -> 3, tree diagrams
 
      elseif(isub.le.140) then
 
C...H: 2 -> 1, tree diagrams, non-standard model processes
 
      elseif(isub.le.160) then
 
C...I: 2 -> 2, tree diagrams, non-standard model processes
 
      elseif(isub.le.200) then
 
CMRENNA++
C...J: 2 -> 2, tree diagrams, SUSY processes
      elseif(isub.le.210) then

      elseif(isub.le.220) then
 
      elseif(isub.le.230) then
 
      elseif(isub.le.240) then
 
      elseif(isub.le.250) then
 
      elseif(isub.le.260) then
 
      elseif(isub.le.270) then
 
      elseif(isub.le.280) then
CMRENNA--
      endif
 
C...Multiply with parton distributions
      if(isub.le.90.or.isub.ge.96) then
        do 1880 ichn=1,nchn
          if(mint(45).ge.2) then
            kfl1=isig(ichn,1)
            sigh(ichn)=sigh(ichn)*xsfx(1,kfl1)
          endif
          if(mint(46).ge.2) then
            kfl2=isig(ichn,2)
            sigh(ichn)=sigh(ichn)*xsfx(2,kfl2)
          endif
          sigs=sigs+sigh(ichn)
 1880   continue
      endif
 
      return
      end
