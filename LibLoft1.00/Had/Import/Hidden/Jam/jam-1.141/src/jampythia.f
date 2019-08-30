C*********************************************************************
C*********************************************************************
C*                                                                  **
C*                                                    March 1997    **
C*                                                                  **
C*           The Lund Monte Carlo for Hadronic Processes            **
C*                                                                  **
C*                        PYTHIA version 6.1                        **
C*                                                                  **
C*                        Torbjorn Sjostrand                        **
C*                Department of Theoretical Physics 2               **
C*                         Lund University                          **
C*               Solvegatan 14A, S-223 62 Lund, Sweden              **
C*                    phone +46 - 46 - 222 48 16                    **
C*                    E-mail torbjorn@thep.lu.se                    **
C*                                                                  **
C*                          SUSY parts by                           **
C*                         Stephen Mrenna                           **
C*                    Argonne National Laboratory                   **
C*          9700 South Cass Avenue, Argonne, IL 60439, USA          **
C*                   phone + 1 - 630 - 252 - 7615                   **
C*                    E-mail mrenna@hep.anl.gov                     **
C*                                                                  **
C*         Several parts are written by Hans-Uno Bengtsson          **
C*          PYSHOW is written together with Mats Bengtsson          **
C*     advanced popcorn baryon production written by Patrik Eden    **
C*     CTEQ 3 parton distributions are by the CTEQ collaboration    **
C*      GRV 94 parton distributions are by Glueck, Reya and Vogt    **
C*   SaS photon parton distributions together with Gerhard Schuler  **
C*     g + g and q + qbar -> t + tbar + H code by Zoltan Kunszt     **
C*         MSSM Higgs mass calculation code by M. Carena,           **
C*           J.R. Espinosa, M. Quiros and C.E.M. Wagner             **
C*         PYGAUS adapted from CERN library (K.S. Kolbig)           **
C*                                                                  **
C*   The latest program version and documentation is found on WWW   **
C*       http://www.thep.lu.se/tf2/staff/torbjorn/Pythia.html       **
C*                                                                  **
C*              Copjright Torbjorn Sjostrand, Lund 1997             **
C*                                                                  **
C*********************************************************************
C*********************************************************************
C                                                                    *
C  List of subprograms in order of appearance, with main purpose     *
C  (S = subroutine, F = function, B = block data)                    *
C                                                                    *
C                                                                    *
C  S   PYINIT   to administer the initialization procedure           *
C  S   PYEVNT   to administer the generation of an event             *
C  S   PYSTAT   to print cross-section and other information         *
C  S   PYINRE   to initialize treatment of resonances                *
C  S   PYINBM   to read in beam, target and frame choices            *
C  S   PYINKI   to initialize kinematics of incoming particles       *
C  S   PYINPR   to set up the selection of included processes        *
C  S   PYXTOT   to give total, elastic and diffractive cross-sect.   *
C  S   PYMAXI   to find differential cross-section maxima            *
C  S   PYPILE   to select multiplicity of pileup events              *
C  S   PYSAVE   to save alternatives for gamma-p and gamma-gamma     *
C  S   PYRAND   to select subprocess and kinematics for event        *
C  S   PYSCAT   to set up kinematics and colour flow of event        *
C  S   PYSSPA   to simulate initial state spacelike showers          *
C  S   PYRESD   to perform resonance decays                          *
C  S   PYMULT   to generate multiple interactions                    *
C  S   PYREMN   to add on target remnants                            *
C  S   PYDIFF   to set up kinematics for diffractive events          *
C  S   PYDOCU   to compute cross-sections and handle documentation   *
C  S   PYFRAM   to perform boosts between different frames           *
C  S   PYWIDT   to calculate full and partial widths of resonances   *
C  S   PYOFSH   to calculate partial width into off-shell channels   *
C  S   PYRECO   to handle colour reconnection in W+W- events         *
C  S   PYKLIM   to calculate borders of allowed kinematical region   *
C  S   PYKMAP   to construct value of kinematical variable           *

C  S   PYSIGH   to calculate differential cross-sections             *

C  F   PYHFTH   to evaluate threshold factor for heavy flavour       *
C  S   PYSPLI   to find flavours left in hadron when one removed     *
C  F   PYGAMM   to evaluate ordinary Gamma function Gamma(x)         *
C  S   PYWAUX   to evaluate auxiliary functions W1(s) and W2(s)      *
C  S   PYI3AU   to evaluate auxiliary function I3(s,t,u,v)           *
C  F   PYSPEN   to evaluate Spence (dilogarithm) function Sp(x)      *
C  S   PYQQBH   to evaluate matrix element for g + g -> Q + Qbar + H *
C                                                                    *
C  S   PYKCUT   dummy routine for user kinematical cuts              *
C  S   PYEVWT   dummy routine for weighting events                   *
C  S   PYUPIN   dummy routine to initialize a user process           *
C  S   PYUPEV   dummy routine to generate a user process event       *
C  S   PYTAUD   dummy routine for interface to tau decay libraries   *
C                                                                    *
C*********************************************************************
 
 
C*********************************************************************
 
C...PYINIT
C...Initializes the generation procedure; finds maxima of the
C...differential cross-sections to be used for weighting.
 
      subroutine pjinit(frame,kfp,kft,win,icon)
 
C...Double precision.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/jydat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
      common/jydat4/chaf(500,2)
      character chaf*16
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      common/pjint5/ngenpd,ngen(0:500,3),xsec(0:500,3)
      save /jydat1/,/jydat2/,/jydat3/,/jydat4/,/pjsubs/,/pjpars/,
     &/pjint1/,/pjint2/,/pjint5/
C...Local arrays and character variables.
      dimension alamin(20),nfin(20)
      character*(*) frame
      character chfram*8,chlh(2)*6
 
C...Interface to PDFLIB.
      common/w50512/qcdl4,qcdl5
      save /w50512/
      double precision value(20),qcdl4,qcdl5
      character*20 parm(20)
      data value/20*0d0/,parm/20*' '/
 
C...Data:Lambda and n_f values for parton distributions; months.
      data alamin/0.20d0,0.29d0,0.20d0,0.40d0,0.213d0,0.208d0,
     $ 0.208d0,0.322d0,0.190d0,0.235d0,0.2d0,
     $ 0.177d0,0.239d0,0.247d0,0.2322d0,0.248d0,0.248d0,
     $ 3*0.2d0/,nfin/20*4/

cc    data alamin/0.177d0,0.239d0,0.247d0,0.2322d0,0.248d0,0.248d0,
cc   &14*0.2d0/,nfin/20*4/
      data chlh/'lepton','hadron'/
 
C...Reset MINT and VINT arrays. Write headers.
      do 100 j=1,400
        mint(j)=0
        vint(j)=0d0
  100 continue
      if(mstu(12).ge.1) call pjlist(0)
      if(mstp(122).ge.1) write(mstu(11),5100)
      icon=0
 
C...Maximum 4 generations; set maximum number of allowed flavours.
      mstp(1)=min(4,mstp(1))
      mstu(114)=min(mstu(114),2*mstp(1))
      mstp(58)=min(mstp(58),2*mstp(1))
 
C...Sum up Cabibbo-Kobayashi-Maskawa factors for each quark/lepton.
      do 120 i=-20,20
        vint(180+i)=0d0
        ia=iabs(i)
        if(ia.ge.1.and.ia.le.2*mstp(1)) then
          do 110 j=1,mstp(1)
            ib=2*j-1+mod(ia,2)
            if(ib.ge.6.and.mstp(9).eq.0) goto 110
            ipm=(5-isign(1,i))/2
            idc=j+mdcy(ia,2)+2
            if(mdme(idc,1).eq.1.or.mdme(idc,1).eq.ipm) vint(180+i)=
     &      vint(180+i)+vckm((ia+1)/2,(ib+1)/2)
  110     continue
        elseif(ia.ge.11.and.ia.le.10+2*mstp(1)) then
          vint(180+i)=1d0
        endif
  120 continue
 
C...Initialize parton distributions: PDFLIB.
      if(mstp(52).eq.2) then
        parm(1)='NPTYPE'
        value(1)=1
        parm(2)='NGROUP'
        value(2)=mstp(51)/1000
        parm(3)='NSET'
        value(3)=mod(mstp(51),1000)
        parm(4)='TMAS'
        value(4)=pmas(6,1)
        call pjpdfset(parm,value)
        mint(93)=1000000+mstp(51)
      endif
 
C...Choose Lambda value to use in alpha-strong.
      mstu(111)=mstp(2)
      if(mstp(3).ge.2) then
        alam=0.2d0
        nf=4
        if(mstp(52).eq.1.and.mstp(51).ge.1.and.mstp(51).le.17) then
          alam=alamin(mstp(51))
          nf=nfin(mstp(51))
        elseif(mstp(52).eq.2) then
          alam=qcdl4
          nf=4
        endif
        parp(1)=alam
        parp(61)=alam
        parp(72)=alam
        paru(112)=alam
        mstu(112)=nf
        if(mstp(3).eq.3) parj(81)=alam
      endif
 
C...Initialize the SUSY generation: couplings, masses,
C...decay modes, branching ratios, and so on.
cjam  call pjmsin
 
C...Initialize widths and partial widths for resonances.
      call pjinre
C...Set Z0 mass and width for e+e- routines.
      parj(123)=pmas(23,1)
      parj(124)=pmas(23,2)
 
C...Identify beam and target particles and frame of process.
cjam++
      chfram=frame//' '
      call pjinbm(chfram,kfp,kft,win)
      if(mint(65).eq.1) goto 170
cjam--
 
C...For gamma-p or gamma-gamma allow many (3 or 6) alternatives.
C...For e-gamma allow 2 alternatives.
      mint(121)=1
      mint(123)=mstp(14)
      if(mstp(14).eq.10.and.(msel.eq.1.or.msel.eq.2)) then
        if((mint(11).eq.22.or.mint(12).eq.22).and.
     &  (iabs(mint(11)).ge.28.or.iabs(mint(12)).ge.28)) mint(121)=3
        if(mint(11).eq.22.and.mint(12).eq.22) mint(121)=6
        if((mint(11).eq.22.or.mint(12).eq.22).and.
     &  (iabs(mint(11)).eq.11.or.iabs(mint(12)).eq.11)) mint(121)=2
      endif
 
C...Set up kinematics of process.
      call pjinki(0)
 
C...Precalculate flavour selection weights
cjam
      if(mstu(123).eq.0.and.mstj(12).gt.0)  call pjkfin
c     call pjkfin
 
C...Loop over gamma-p or gamma-gamma alternatives.
      do 160 iga=1,mint(121)
        mint(122)=iga
 
C...Select partonic subprocesses to be included in the simulation.
        call pjinpr
 
C...Count number of subprocesses on.
        mint(48)=0
        do 130 isub=1,500
          if(mint(50).eq.0.and.isub.ge.91.and.isub.le.96.and.
     &    msub(isub).eq.1) then
            write(mstu(11),5200) isub,chlh(mint(41)),chlh(mint(42))
            stop
          elseif(msub(isub).eq.1.and.iset(isub).eq.-1) then
            write(mstu(11),5300) isub
            stop
          elseif(msub(isub).eq.1.and.iset(isub).le.-2) then
            write(mstu(11),5400) isub
            stop
          elseif(msub(isub).eq.1) then
            mint(48)=mint(48)+1
          endif
  130   continue
        if(mint(48).eq.0) then
          write(mstu(11),5500)
          stop
        endif
        mint(49)=mint(48)-msub(91)-msub(92)-msub(93)-msub(94)
 
C...Reset variables for cross-section calculation.
        do 150 i=0,500
          do 140 j=1,3
            ngen(i,j)=0
            xsec(i,j)=0d0
  140     continue
  150   continue
 
C...Find parametrized total cross-sections.
        call pjxtot
 
C...Maxima of differential cross-sections.
        if(mstp(121).le.1) call pjmaxi
 
C...Initialize possibility of pileup events.
        if(mint(121).gt.1) mstp(131)=0
        if(mstp(131).ne.0) call pjpile(1)
 
C...Initialize multiple interactions with variable impact parameter.
        if(mint(50).eq.1.and.(mint(49).ne.0.or.mstp(131).ne.0).and.
     &  mstp(82).ge.2) call pjmult(1)
 
C...Save results for gamma-p and gamma-gamma alternatives.
        if(mint(121).gt.1) call pjsave(1,iga)
  160 continue
 
C...Initialization finished.
  170 if(mstp(122).ge.1) write(mstu(11),5600)
 
C...Formats for initialization information.
 5100 format('1',18('*'),1x,'PYINIT: initialization of PYTHIA ',
     &'routines',1x,17('*'))
 5200 format(1x,'Error: process number ',i3,' not meaningful for ',a6,
     &'-',a6,' interactions.'/1x,'Execution stopped!')
 5300 format(1x,'Error: requested subprocess',i4,' not implemented.'/
     &1x,'Execution stopped!')
 5400 format(1x,'Error: requested subprocess',i4,' not existing.'/
     &1x,'Execution stopped!')
 5500 format(1x,'Error: no subprocess switched on.'/
     &1x,'Execution stopped.')
 5600 format(/1x,22('*'),1x,'PYINIT: initialization completed',1x,
     &22('*'))
 
      return
      end
 
C*********************************************************************
 
C...PYEVNT
C...Administers the generation of a high-pT event via calls to
C...a number of subroutines.
 
      subroutine pjevnt
 
C...Double precision.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      common/pjint4/mwid(500),wids(500,5)
      common/pjint5/ngenpd,ngen(0:500,3),xsec(0:500,3)
      common/pjuppr/nup,kup(20,7),nfup,ifup(10,2),pup(20,5),q2up(0:10)
      save /jyjets/,/jydat1/,/jydat2/,/pjpars/,/pjint1/,/pjint2/,
     &/pjint4/,/pjint5/,/pjuppr/
C...Local array.
      dimension vtx(4)
 
C...Initial values for some counters.
      n=0
      mint(5)=mint(5)+1
      mint(7)=0           ! line number of outgoing scattered parton 1.
      mint(8)=0           ! line number of outgoing scattered parton 2.
      mint(83)=0
      mint(84)=mstp(126)
      mstu(24)=0          ! type of latest error
      mstu70=0
      mstj14=mstj(14)
 
C...If variable energies: redo incoming kinematics and cross-section.
      msti(61)=0 ! flag =0: o.k. 1:event was not generated.
      if(mstp(171).eq.1) then
        call pjinki(1)
        if(msti(61).eq.1) then
          mint(5)=mint(5)-1
          return
        endif
        if(mint(121).gt.1) call pjsave(3,1)
        call pjxtot
      endif
 
C...Loop over number of pileup events; check space left.
      if(mstp(131).le.0) then
        npile=1
      else
        call pjpile(2)
        npile=mint(81)
      endif
      do 260 ipile=1,npile
        if(mint(84)+100.ge.mstu(4)) then
          call pjerrm(11,
     &    '(PYEVNT:) no more space in PYJETS for pileup events')
          if(mstu(21).ge.1) goto 270
        endif
        mint(82)=ipile
 
C...Generate variables of hard scattering.
        mint(51)=0
        msti(52)=0
  100   continue
        if(mint(51).ne.0.or.mstu(24).ne.0) msti(52)=msti(52)+1
        mint(31)=0
        mint(51)=0
        mint(57)=0

        call pjrand

        if(msti(61).eq.1) then
          mint(5)=mint(5)-1
          return
        endif
        if(mint(51).eq.2) return
        isub=mint(1)
        if(mstp(111).eq.-1) goto 250
 
        if(isub.le.90.or.isub.ge.95) then
C...Hard scattering (including low-pT):
C...reconstruct kinematics and colour flow of hard scattering.
  110     mint(51)=0
          call pjscat
          if(mint(51).eq.1) goto 100
          ipu1=mint(84)+1
          ipu2=mint(84)+2
          if(isub.eq.95) goto 130
 
C...Showering of initial state partons (optional).
          alamsv=parj(81)
          parj(81)=parp(72)
          if(mstp(61).ge.1.and.mint(47).ge.2) call pjsspa(ipu1,ipu2)
          parj(81)=alamsv
          if(mint(51).eq.1) goto 100
 
C...Showering of final state partons (optional).
          alamsv=parj(81)
          parj(81)=parp(72)
          if(mstp(71).ge.1.and.iset(isub).ge.2.and.iset(isub).le.10)
     &    then
            ipu3=mint(84)+3
            ipu4=mint(84)+4
            if(iset(isub).eq.5) ipu4=-3
            qmax=vint(55)
            if(iset(isub).eq.2) qmax=sqrt(parp(71))*vint(55)
            call pjshow(ipu3,ipu4,qmax)
          elseif(mstp(71).ge.1.and.iset(isub).eq.11.and.nfup.ge.1) then
            do 120 iup=1,nfup
              ipu3=ifup(iup,1)+mint(84)
              ipu4=ifup(iup,2)+mint(84)
              qmax=sqrt(max(0d0,q2up(iup)))
              call pjshow(ipu3,ipu4,qmax)
  120       continue
          endif
          parj(81)=alamsv
 
C...Decay of final state resonances.
          mint(32)=0
          if(mstp(41).ge.1.and.iset(isub).le.10) call pjresd(0)
          if(mint(51).eq.1) goto 100
          mint(52)=n
 
C...Multiple interactions.
          if(mstp(81).ge.1.and.mint(50).eq.1) call pjmult(6)
          mint(53)=n
 
C...Hadron remnants and primordial kT.
  130     call pjremn(ipu1,ipu2)
          if(mint(51).eq.1.and.mint(57).ge.1.and.mint(57).le.5) goto 110
          if(mint(51).eq.1) goto 100
 
        else
C...Diffractive and elastic scattering.
          call pjdiff
        endif
 
C...Check that no odd resonance left undecayed.
        if(mstp(111).ge.1) then
          nfix=n
          do 140 i=mint(84)+1,nfix
            if(k(i,1).ge.1.and.k(i,1).le.10.and.k(i,2).ne.21.and.
     &      k(i,2).ne.22) then
              if(mwid(jamcomp(k(i,2))).ne.0) then
                call pjresd(i)
                if(mint(51).eq.1) goto 100
              endif
            endif
  140     continue
        endif
 
C...Recalculate energies from momenta and masses (if desired).
        if(mstp(113).ge.1) then
          do 150 i=mint(83)+1,n
            if(k(i,1).gt.0.and.k(i,1).le.10) p(i,4)=sqrt(p(i,1)**2+
     &      p(i,2)**2+p(i,3)**2+p(i,5)**2)
  150     continue
          nrecal=n
        endif
 
C...Rearrange partons along strings, check invariant mass cuts.
        mstu(28)=0
        if(mstp(111).le.0) mstj(14)=-1
        call pjprep(mint(84)+1)
        mstj(14)=mstj14
        if(mstp(112).eq.1.and.mstu(28).eq.3) goto 100

cjam++
C...Go back to lab frame (needed for vertices, also in fragmentation).
        call pjfram(1)
c...Calculate parton position.
        call jamhrdv
cjam--

        if(mstp(125).eq.0.or.mstp(125).eq.1) then
          do 180 i=mint(84)+1,n
            if(k(i,2).eq.94) then
              do 170 i1=i+1,min(n,i+3)
                if(k(i1,3).eq.i) then
                  k(i1,3)=mod(k(i1,4)/mstu(5),mstu(5))
                  if(k(i1,3).eq.0) then
                    do 160 ii=mint(84)+1,i-1
                        if(k(ii,2).eq.k(i1,2)) then
                          if(mod(k(ii,4),mstu(5)).eq.i1.or.
     &                    mod(k(ii,5),mstu(5)).eq.i1) k(i1,3)=ii
                        endif
  160               continue
                    if(k(i+1,3).eq.0) k(i+1,3)=k(i,3)
                  endif
                endif
  170         continue
            endif
  180     continue
          call pjedit(12)
          call pjedit(14)
          if(mstp(125).eq.0) call pjedit(15)
          if(mstp(125).eq.0) mint(4)=0
          do 200 i=mint(83)+1,n
            if(k(i,1).eq.11.and.k(i,4).eq.0.and.k(i,5).eq.0) then
              do 190 i1=i+1,n
                if(k(i1,3).eq.i.and.k(i,4).eq.0) k(i,4)=i1
                if(k(i1,3).eq.i) k(i,5)=i1
  190         continue
            endif
  200     continue
        endif
 
C...Introduce separators between sections in PYLIST event listing.
        if(ipile.eq.1.and.mstp(125).le.0) then
          mstu70=1
          mstu(71)=n
        elseif(ipile.eq.1) then
          mstu70=3
          mstu(71)=2
          mstu(72)=mint(4)
          mstu(73)=n
        endif
 
cjam
C...Go back to lab frame (needed for vertices, also in fragmentation).
c       call pjfram(1)
 
C...Set nonvanishing production vertex (optional).
        if(mstp(151).eq.1) then
          do 210 j=1,4
            vtx(j)=parp(150+j)*sqrt(-2d0*log(max(1d-10,pjr(0))))*
     &      sin(paru(2)*pjr(0))
  210     continue
          do 230 i=mint(83)+1,n
            do 220 j=1,4
              v(i,j)=v(i,j)+vtx(j)
  220       continue
  230     continue
        endif
 
C...Perform hadronization (if desired).
        if(mstp(111).ge.1) then
          call pjexec
          if(mstu(24).ne.0) goto 100
        endif
        if(mstp(113).ge.1) then
          do 240 i=nrecal,n
            if(p(i,5).gt.0d0) p(i,4)=sqrt(p(i,1)**2+
     &      p(i,2)**2+p(i,3)**2+p(i,5)**2)
  240     continue
        endif
        if(mstp(125).eq.0.or.mstp(125).eq.1) call pjedit(14)
 
C...Store event information and calculate Monte Carlo estimates of
C...subprocess cross-sections.
  250   if(ipile.eq.1) call pjdocu
 
C...Set counters for current pileup event and loop to next one.
        msti(41)=ipile
        if(ipile.ge.2.and.ipile.le.10) msti(40+ipile)=isub
        if(mstu70.lt.10) then
          mstu70=mstu70+1
          mstu(70+mstu70)=n
        endif
        mint(83)=n
        mint(84)=n+mstp(126)
        if(ipile.lt.npile) call pjfram(2)
  260 continue
 
C...Generic information on pileup events. Reconstruct missing history.
      if(mstp(131).eq.1.and.mstp(133).ge.1) then
        pari(91)=vint(132)
        pari(92)=vint(133)
        pari(93)=vint(134)
        if(mstp(133).ge.2) pari(93)=pari(93)*xsec(0,3)/vint(131)
      endif
      call pjedit(16)
 
C...Transform to the desired coordinate frame.
  270 call pjfram(mstp(124))
      mstu(70)=mstu70
      paru(21)=vint(1)
 
      return
      end
 
C***********************************************************************
 
C...PYSTAT
C...Prints out information about cross-sections, decay widths, branching
C...ratios, kinematical limits, status codes and parameter values.
 
      subroutine pjstat(mstat)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Parameter statement to help give large particle numbers.
      parameter (ksusy1=1000000,ksusy2=2000000,kexcit=4000000)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/jydat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      common/pjint4/mwid(500),wids(500,5)
      common/pjint5/ngenpd,ngen(0:500,3),xsec(0:500,3)
      common/pjint6/proc(0:500)
      character proc*28
      common/pjmssm/imss(0:99),rmss(0:99)
      save /jydat1/,/jydat2/,/jydat3/,/pjsubs/,/pjpars/,/pjint1/,
     &/pjint2/,/pjint4/,/pjint5/,/pjint6/,/pjmssm/
C...Local arrays, character variables and data.
      dimension wdtp(0:200),wdte(0:200,0:5)
      character proga(6)*28,chau*16,chkf*16,chd1*16,chd2*16,chd3*16,
     &chin(2)*12,state(-1:5)*4,chkin(21)*18,disga(2)*28
      data proga/
     &'VMD/hadron * VMD            ','VMD/hadron * direct         ',
     &'VMD/hadron * anomalous      ','direct * direct             ',
     &'direct * anomalous          ','anomalous * anomalous       '/
      data disga/'e * VMD','e * anomalous'/
      data state/'----','off ','on  ','on/+','on/-','on/1','on/2'/,
     &chkin/' m_hard (GeV/c^2) ',' p_T_hard (GeV/c) ',
     &'m_finite (GeV/c^2)','   y*_subsystem   ','     y*_large     ',
     &'     y*_small     ','    eta*_large    ','    eta*_small    ',
     &'cos(theta*)_large ','cos(theta*)_small ','       x_1        ',
     &'       x_2        ','       x_F        ',' cos(theta_hard)  ',
     &'m''_hard (GeV/c^2) ','       tau        ','        y*        ',
     &'cos(theta_hard^-) ','cos(theta_hard^+) ','      x_T^2       ',
     &'       tau''       '/
 
C...Cross-sections.
      if(mstat.le.1) then
        if(mint(121).gt.1) call pjsave(5,0)
        write(mstu(11),5000)
        write(mstu(11),5100)
        write(mstu(11),5200) 0,proc(0),ngen(0,3),ngen(0,1),xsec(0,3)
        do 100 i=1,500
          if(msub(i).ne.1) goto 100
          write(mstu(11),5200) i,proc(i),ngen(i,3),ngen(i,1),xsec(i,3)
  100   continue
        if(mint(121).gt.1) then
          write(mstu(11),5300)
          do 110 iga=1,mint(121)
            call pjsave(3,iga)
            if(mint(121).eq.2) then
              write(mstu(11),5200) iga,disga(iga),ngen(0,3),ngen(0,1),
     &        xsec(0,3)
            else
              write(mstu(11),5200) iga,proga(iga),ngen(0,3),ngen(0,1),
     &        xsec(0,3)
            endif
  110     continue
          call pjsave(5,0)
        endif
        write(mstu(11),5400) 1d0-dble(ngen(0,3))/
     &  max(1d0,dble(ngen(0,2)))
 
C...Decay widths and branching ratios.
      elseif(mstat.eq.2) then
        write(mstu(11),5500)
        write(mstu(11),5600)
        do 140 kc=1,500
          kf=kchg(kc,4)
          call pjname(kf,chkf)
          ioff=0
          if(kc.le.22) then
            if(kc.gt.2*mstp(1).and.kc.le.10) goto 140
            if(kc.gt.10+2*mstp(1).and.kc.le.20) goto 140
            if(kc.le.5.or.(kc.ge.11.and.kc.le.16)) ioff=1
            if(kc.eq.18.and.pmas(18,1).lt.1d0) ioff=1
            if(kc.eq.21.or.kc.eq.22) ioff=1
          else
            if(mwid(kc).le.0) goto 140
            if(imss(1).le.0.and.(kf/ksusy1.eq.1.or.
     &      kf/ksusy1.eq.2)) goto 140
          endif
C...Off-shell branchings.
          if(ioff.eq.1) then
            ngp=0
            if(kc.le.20) ngp=(mod(kc,10)+1)/2
            if(ngp.le.mstp(1)) write(mstu(11),5700) kf,chkf(1:10),
     &      pmas(kc,1),0d0,0d0,state(mdcy(kc,1)),0d0
            do 120 j=1,mdcy(kc,3)
              idc=j+mdcy(kc,2)-1
              ngp1=0
              if(iabs(kfdp(idc,1)).le.20) ngp1=
     &        (mod(iabs(kfdp(idc,1)),10)+1)/2
              ngp2=0
              if(iabs(kfdp(idc,2)).le.20) ngp2=
     &        (mod(iabs(kfdp(idc,2)),10)+1)/2
              call pjname(kfdp(idc,1),chd1)
              call pjname(kfdp(idc,2),chd2)
              if(kfdp(idc,3).eq.0) then
                if(mdme(idc,2).eq.102.and.ngp1.le.mstp(1).and.
     &          ngp2.le.mstp(1)) write(mstu(11),5800) idc,chd1(1:10),
     &          chd2(1:10),0d0,0d0,state(mdme(idc,1)),0d0
              else
                call pjname(kfdp(idc,3),chd3)
                if(mdme(idc,2).eq.102.and.ngp1.le.mstp(1).and.
     &          ngp2.le.mstp(1)) write(mstu(11),5900) idc,chd1(1:10),
     &          chd2(1:10),chd3(1:10),0d0,0d0,state(mdme(idc,1)),0d0
              endif
  120       continue
C...On-shell decays.
          else
            call pjwidt(kf,pmas(kc,1)**2,wdtp,wdte)
            brfin=1d0
            if(wdte(0,0).le.0d0) brfin=0d0
            write(mstu(11),5700) kf,chkf(1:10),pmas(kc,1),wdtp(0),1d0,
     &      state(mdcy(kc,1)),brfin
            do 130 j=1,mdcy(kc,3)
              idc=j+mdcy(kc,2)-1
              ngp1=0
              if(iabs(kfdp(idc,1)).le.20) ngp1=
     &        (mod(iabs(kfdp(idc,1)),10)+1)/2
              ngp2=0
              if(iabs(kfdp(idc,2)).le.20) ngp2=
     &        (mod(iabs(kfdp(idc,2)),10)+1)/2
              brfin=0d0
              if(wdte(0,0).gt.0d0) brfin=wdte(j,0)/wdte(0,0)
              call pjname(kfdp(idc,1),chd1)
              call pjname(kfdp(idc,2),chd2)
              if(kfdp(idc,3).eq.0) then
                if(ngp1.le.mstp(1).and.ngp2.le.mstp(1))
     &          write(mstu(11),5800) idc,chd1(1:10),
     &          chd2(1:10),wdtp(j),wdtp(j)/wdtp(0),
     &          state(mdme(idc,1)),brfin
              else
                call pjname(kfdp(idc,3),chd3)
                if(ngp1.le.mstp(1).and.ngp2.le.mstp(1))
     &          write(mstu(11),5900) idc,chd1(1:10),
     &          chd2(1:10),chd3(1:10),wdtp(j),wdtp(j)/wdtp(0),
     &          state(mdme(idc,1)),brfin
              endif
  130       continue
          endif
  140   continue
        write(mstu(11),6000)
 
C...Allowed incoming partons/particles at hard interaction.
      elseif(mstat.eq.3) then
        write(mstu(11),6100)
        call pjname(mint(11),chau)
        chin(1)=chau(1:12)
        call pjname(mint(12),chau)
        chin(2)=chau(1:12)
        write(mstu(11),6200) chin(1),chin(2)
        do 150 i=-20,22
          if(i.eq.0) goto 150
          ia=iabs(i)
          if(ia.gt.mstp(58).and.ia.le.10) goto 150
          if(ia.gt.10+2*mstp(1).and.ia.le.20) goto 150
          call pjname(i,chau)
          write(mstu(11),6300) chau,state(kfin(1,i)),chau,
     &    state(kfin(2,i))
  150   continue
        write(mstu(11),6400)
 
C...User-defined limits on kinematical variables.
      elseif(mstat.eq.4) then
        write(mstu(11),6500)
        write(mstu(11),6600)
        shrmax=ckin(2)
        if(shrmax.lt.0d0) shrmax=vint(1)
        write(mstu(11),6700) ckin(1),chkin(1),shrmax
        pthmin=max(ckin(3),ckin(5))
        pthmax=ckin(4)
        if(pthmax.lt.0d0) pthmax=0.5d0*shrmax
        write(mstu(11),6800) ckin(3),pthmin,chkin(2),pthmax
        write(mstu(11),6900) chkin(3),ckin(6)
        do 160 i=4,14
          write(mstu(11),6700) ckin(2*i-1),chkin(i),ckin(2*i)
  160   continue
        sprmax=ckin(32)
        if(sprmax.lt.0d0) sprmax=vint(1)
        write(mstu(11),6700) ckin(31),chkin(15),sprmax
        write(mstu(11),7000)
 
C...Status codes and parameter values.
      elseif(mstat.eq.5) then
        write(mstu(11),7100)
        write(mstu(11),7200)
        do 170 i=1,100
          write(mstu(11),7300) i,mstp(i),parp(i),100+i,mstp(100+i),
     &    parp(100+i)
  170   continue
 
C...List of all processes implemented in the program.
      elseif(mstat.eq.6) then
        write(mstu(11),7400)
        write(mstu(11),7500)
        do 180 i=1,500
          if(iset(i).lt.0) goto 180
          write(mstu(11),7600) i,proc(i),iset(i),kfpr(i,1),kfpr(i,2)
  180   continue
        write(mstu(11),7700)
      endif
 
C...Formats for printouts.
 5000 format('1',9('*'),1x,'PYSTAT:  Statistics on Number of ',
     &'Events and Cross-sections',1x,9('*'))
 5100 format(/1x,78('=')/1x,'I',34x,'I',28x,'I',12x,'I'/1x,'I',12x,
     &'Subprocess',12x,'I',6x,'Number of points',6x,'I',4x,'Sigma',3x,
     &'I'/1x,'I',34x,'I',28x,'I',12x,'I'/1x,'I',34('-'),'I',28('-'),
     &'I',4x,'(mb)',4x,'I'/1x,'I',34x,'I',28x,'I',12x,'I'/1x,'I',1x,
     &'N:o',1x,'Type',25x,'I',4x,'Generated',9x,'Tried',1x,'I',12x,
     &'I'/1x,'I',34x,'I',28x,'I',12x,'I'/1x,78('=')/1x,'I',34x,'I',28x,
     &'I',12x,'I')
 5200 format(1x,'I',1x,i3,1x,a28,1x,'I',1x,i12,1x,i13,1x,'I',1x,1p,
     &d10.3,1x,'I')
 5300 format(1x,'I',34x,'I',28x,'I',12x,'I'/1x,78('=')/
     &1x,'I',34x,'I',28x,'I',12x,'I')
 5400 format(1x,'I',34x,'I',28x,'I',12x,'I'/1x,78('=')//
     &1x,'********* Fraction of events that fail fragmentation ',
     &'cuts =',1x,f8.5,' *********'/)
 5500 format('1',27('*'),1x,'PYSTAT:  Decay Widths and Branching ',
     &'Ratios',1x,27('*'))
 5600 format(/1x,98('=')/1x,'I',49x,'I',13x,'I',12x,'I',6x,'I',12x,'I'/
     &1x,'I',5x,'Mother  -->  Branching/Decay Channel',8x,'I',1x,
     &'Width (GeV)',1x,'I',7x,'B.R.',1x,'I',1x,'Stat',1x,'I',2x,
     &'Eff. B.R.',1x,'I'/1x,'I',49x,'I',13x,'I',12x,'I',6x,'I',12x,'I'/
     &1x,98('='))
 5700 format(1x,'I',49x,'I',13x,'I',12x,'I',6x,'I',12x,'I'/1x,'I',1x,
     &i8,2x,a10,3x,'(m =',f10.3,')',2x,'-->',5x,'I',2x,1p,d10.3,0p,1x,
     &'I',1x,1p,d10.3,0p,1x,'I',1x,a4,1x,'I',1x,1p,d10.3,0p,1x,'I')
 5800 format(1x,'I',1x,i8,2x,a10,1x,'+',1x,a10,15x,'I',2x,
     &1p,d10.3,0p,1x,'I',1x,1p,d10.3,0p,1x,'I',1x,a4,1x,'I',1x,
     &1p,d10.3,0p,1x,'I')
 5900 format(1x,'I',1x,i8,2x,a10,1x,'+',1x,a10,1x,'+',1x,a10,2x,'I',2x,
     &1p,d10.3,0p,1x,'I',1x,1p,d10.3,0p,1x,'I',1x,a4,1x,'I',1x,
     &1p,d10.3,0p,1x,'I')
 6000 format(1x,'I',49x,'I',13x,'I',12x,'I',6x,'I',12x,'I'/1x,98('='))
 6100 format('1',7('*'),1x,'PYSTAT: Allowed Incoming Partons/',
     &'Particles at Hard Interaction',1x,7('*'))
 6200 format(/1x,78('=')/1x,'I',38x,'I',37x,'I'/1x,'I',1x,
     &'Beam particle:',1x,a12,10x,'I',1x,'Target particle:',1x,a12,7x,
     &'I'/1x,'I',38x,'I',37x,'I'/1x,'I',1x,'Content',6x,'State',19x,
     &'I',1x,'Content',6x,'State',18x,'I'/1x,'I',38x,'I',37x,'I'/1x,
     &78('=')/1x,'I',38x,'I',37x,'I')
 6300 format(1x,'I',1x,a9,5x,a4,19x,'I',1x,a9,5x,a4,18x,'I')
 6400 format(1x,'I',38x,'I',37x,'I'/1x,78('='))
 6500 format('1',12('*'),1x,'PYSTAT: User-Defined Limits on ',
     &'Kinematical Variables',1x,12('*'))
 6600 format(/1x,78('=')/1x,'I',76x,'I')
 6700 format(1x,'I',16x,1p,d10.3,0p,1x,'<',1x,a,1x,'<',1x,1p,d10.3,0p,
     &16x,'I')
 6800 format(1x,'I',3x,1p,d10.3,0p,1x,'(',1p,d10.3,0p,')',1x,'<',1x,a,
     &1x,'<',1x,1p,d10.3,0p,16x,'I')
 6900 format(1x,'I',29x,a,1x,'=',1x,1p,d10.3,0p,16x,'I')
 7000 format(1x,'I',76x,'I'/1x,78('='))
 7100 format('1',12('*'),1x,'PYSTAT: Summary of Status Codes and ',
     &'Parameter Values',1x,12('*'))
 7200 format(/3x,'I',4x,'MSTP(I)',9x,'PARP(I)',20x,'I',4x,'MSTP(I)',9x,
     &'PARP(I)'/)
 7300 format(1x,i3,5x,i6,6x,1p,d10.3,0p,18x,i3,5x,i6,6x,1p,d10.3)
 7400 format('1',13('*'),1x,'PYSTAT: List of implemented processes',
     &1x,13('*'))
 7500 format(/1x,65('=')/1x,'I',34x,'I',28x,'I'/1x,'I',12x,
     &'Subprocess',12x,'I',1x,'ISET',2x,'KFPR(I,1)',2x,'KFPR(I,2)',1x,
     &'I'/1x,'I',34x,'I',28x,'I'/1x,65('=')/1x,'I',34x,'I',28x,'I')
 7600 format(1x,'I',1x,i3,1x,a28,1x,'I',1x,i4,1x,i10,1x,i10,1x,'I')
 7700 format(1x,'I',34x,'I',28x,'I'/1x,65('='))
 
      return
      end
 
C*********************************************************************
 
C...PYINRE
C...Calculates full and effective widths of gauge bosons, stores
C...masses and widths, rescales coefficients to be used for
C...resonance production generation.
 
      subroutine pjinre
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Parameter statement to help give large particle numbers.
      parameter (ksusy1=1000000,ksusy2=2000000,kexcit=4000000)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/jydat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
      common/jydat4/chaf(500,2)
      character chaf*16
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      common/pjint4/mwid(500),wids(500,5)
      common/pjint6/proc(0:500)
      character proc*28
      common/pjmssm/imss(0:99),rmss(0:99)
      save /jydat1/,/jydat2/,/jydat3/,/jydat4/,/pjsubs/,/pjpars/,
     &/pjint1/,/pjint2/,/pjint4/,/pjint6/,/pjmssm/
C...Local arrays and data.
      dimension wdtp(0:200),wdte(0:200,0:5),wdtpm(0:200),
     &wdtem(0:200,0:5),kcord(500),pmord(500)
 
C...Born level couplings in MSSM Higgs doublet sector.
      xw=paru(102)
      xwv=xw
      if(mstp(8).ge.2) xw=1d0-(pmas(24,1)/pmas(23,1))**2
      xw1=1d0-xw
      if(mstp(4).eq.2) then
        tanbe=paru(141)
        ratbe=((1d0-tanbe**2)/(1d0+tanbe**2))**2
        sqmz=pmas(23,1)**2
        sqmw=pmas(24,1)**2
        sqmh=pmas(25,1)**2
        sqma=sqmh*(sqmz-sqmh)/(sqmz*ratbe-sqmh)
        sqmhp=0.5d0*(sqma+sqmz+sqrt((sqma+sqmz)**2-4d0*sqma*sqmz*ratbe))
        sqmhc=sqma+sqmw
        if(sqmh.ge.sqmz.or.min(sqma,sqmhp,sqmhc).le.0d0) then
          write(mstu(11),5000)
          stop
        endif
        pmas(35,1)=sqrt(sqmhp)
        pmas(36,1)=sqrt(sqma)
        pmas(37,1)=sqrt(sqmhc)
        alsu=0.5d0*atan(2d0*tanbe*(sqma+sqmz)/((1d0-tanbe**2)*
     &  (sqma-sqmz)))
        besu=atan(tanbe)
        paru(142)=1d0
        paru(143)=1d0
        paru(161)=-sin(alsu)/cos(besu)
        paru(162)=cos(alsu)/sin(besu)
        paru(163)=paru(161)
        paru(164)=sin(besu-alsu)
        paru(165)=paru(164)
        paru(168)=sin(besu-alsu)+0.5d0*cos(2d0*besu)*sin(besu+alsu)/xw
        paru(171)=cos(alsu)/cos(besu)
        paru(172)=sin(alsu)/sin(besu)
        paru(173)=paru(171)
        paru(174)=cos(besu-alsu)
        paru(175)=paru(174)
        paru(176)=cos(2d0*alsu)*cos(besu+alsu)-2d0*sin(2d0*alsu)*
     &  sin(besu+alsu)
        paru(177)=cos(2d0*besu)*cos(besu+alsu)
        paru(178)=cos(besu-alsu)-0.5d0*cos(2d0*besu)*cos(besu+alsu)/xw
        paru(181)=tanbe
        paru(182)=1d0/tanbe
        paru(183)=paru(181)
        paru(184)=0d0
        paru(185)=paru(184)
        paru(186)=cos(besu-alsu)
        paru(187)=sin(besu-alsu)
        paru(188)=paru(186)
        paru(189)=paru(187)
        paru(190)=0d0
        paru(195)=cos(besu-alsu)
      endif
 
C...Reset effective widths of gauge bosons.
      do 110 i=1,500
        do 100 j=1,5
          wids(i,j)=1d0
  100   continue
  110 continue
 
C...Order resonances by increasing mass (except Z0 and W+/-).
      nres=0
      do 140 kc=1,500
        kf=kchg(kc,4)
        if(kf.eq.0) goto 140
        if(mwid(kc).eq.0) goto 140
        if(kc.eq.7.or.kc.eq.8.or.kc.eq.17.or.kc.eq.18) then
          if(mstp(1).le.3) goto 140
        endif
        if(kf/ksusy1.eq.1.or.kf/ksusy1.eq.2) then
          if(imss(1).le.0) goto 140
        endif
        nres=nres+1
        pmres=pmas(kc,1)
        if(kc.eq.23.or.kc.eq.24) pmres=0d0
        do 120 i1=nres-1,1,-1
          if(pmres.ge.pmord(i1)) goto 130
          kcord(i1+1)=kcord(i1)
          pmord(i1+1)=pmord(i1)
  120   continue
  130   kcord(i1+1)=kc
        pmord(i1+1)=pmres
  140 continue
 
C...Loop over possible resonances.
      do 180 i=1,nres
        kc=kcord(i)
        kf=kchg(kc,4)
 
C...Check that no fourth generation channels on by mistake.
        if(mstp(1).le.3) then
          do 150 j=1,mdcy(kc,3)
            idc=j+mdcy(kc,2)-1
            kfa1=iabs(kfdp(idc,1))
            kfa2=iabs(kfdp(idc,2))
            if(kfa1.eq.7.or.kfa1.eq.8.or.kfa1.eq.17.or.kfa1.eq.18.or.
     &      kfa2.eq.7.or.kfa2.eq.8.or.kfa2.eq.17.or.kfa2.eq.18)
     &      mdme(idc,1)=-1
  150     continue
        endif
 
C...Check that no supersymmetric channels on by mistake.
        if(imss(1).le.0) then
          do 160 j=1,mdcy(kc,3)
            idc=j+mdcy(kc,2)-1
            kfa1s=iabs(kfdp(idc,1))/ksusy1
            kfa2s=iabs(kfdp(idc,2))/ksusy1
            if(kfa1s.eq.1.or.kfa1s.eq.2.or.kfa2s.eq.1.or.kfa2s.eq.2)
     &      mdme(idc,1)=-1
  160     continue
        endif
 
C...Find mass and evaluate width.
        pmr=pmas(kc,1)
        if(kf.eq.25.or.kf.eq.35.or.kf.eq.36) mint(62)=1
        if(mwid(kc).eq.3) mint(63)=1
        call pjwidt(kf,pmr**2,wdtp,wdte)
        mint(51)=0
 
C...Evaluate suppression factors due to non-simulated channels.
        if(kchg(kc,3).eq.0) then
          wids(kc,1)=((wdte(0,1)+wdte(0,2))**2+
     &    2d0*(wdte(0,1)+wdte(0,2))*(wdte(0,4)+wdte(0,5))+
     &    2d0*wdte(0,4)*wdte(0,5))/wdtp(0)**2
          wids(kc,2)=(wdte(0,1)+wdte(0,2)+wdte(0,4))/wdtp(0)
          wids(kc,3)=0d0
          wids(kc,4)=0d0
          wids(kc,5)=0d0
        else
          if(mwid(kc).eq.3) mint(63)=1
          call pjwidt(-kf,pmr**2,wdtpm,wdtem)
          mint(51)=0
          wids(kc,1)=((wdte(0,1)+wdte(0,2))*(wdtem(0,1)+wdtem(0,3))+
     &    (wdte(0,1)+wdte(0,2))*(wdtem(0,4)+wdtem(0,5))+
     &    (wdte(0,4)+wdte(0,5))*(wdtem(0,1)+wdtem(0,3))+
     &    wdte(0,4)*wdtem(0,5)+wdte(0,5)*wdtem(0,4))/wdtp(0)**2
          wids(kc,2)=(wdte(0,1)+wdte(0,2)+wdte(0,4))/wdtp(0)
          wids(kc,3)=(wdtem(0,1)+wdtem(0,3)+wdtem(0,4))/wdtp(0)
          wids(kc,4)=((wdte(0,1)+wdte(0,2))**2+
     &    2d0*(wdte(0,1)+wdte(0,2))*(wdte(0,4)+wdte(0,5))+
     &    2d0*wdte(0,4)*wdte(0,5))/wdtp(0)**2
          wids(kc,5)=((wdtem(0,1)+wdtem(0,3))**2+
     &    2d0*(wdtem(0,1)+wdtem(0,3))*(wdtem(0,4)+wdtem(0,5))+
     &    2d0*wdtem(0,4)*wdtem(0,5))/wdtp(0)**2
        endif
 
C...Set resonance widths and branching ratios;
C...also on/off switch for decays.
        if(mwid(kc).eq.1.or.mwid(kc).eq.3) then
          pmas(kc,2)=wdtp(0)
          pmas(kc,3)=min(0.9d0*pmas(kc,1),10d0*pmas(kc,2))
          mdcy(kc,1)=mstp(41)
          do 170 j=1,mdcy(kc,3)
            idc=j+mdcy(kc,2)-1
            brat(idc)=0d0
            if(wdtp(0).gt.0d0) brat(idc)=wdtp(j)/wdtp(0)
  170     continue
        endif
  180 continue
 
C...Flavours of leptoquark: redefine charge and name.
      kflqq=kfdp(mdcy(39,2),1)
      kflql=kfdp(mdcy(39,2),2)
      kchg(39,1)=kchg(jamcomp(kflqq),1)*isign(1,kflqq)+
     &kchg(jamcomp(kflql),1)*isign(1,kflql)
      ll=1
      if(iabs(kflql).eq.13) ll=2
      if(iabs(kflql).eq.15) ll=3
      chaf(39,1)='LQ_'//chaf(iabs(kflqq),1)(1:1)//
     &chaf(iabs(kflql),1)(1:ll)//' '
      chaf(39,2)=chaf(39,2)(1:4+ll)//'bar '
 
C...Special cases in treatment of gamma*/Z0: redefine process name.
      if(mstp(43).eq.1) then
        proc(1)='f + fbar -> gamma*'
        proc(15)='f + fbar -> g + gamma*'
        proc(19)='f + fbar -> gamma + gamma*'
        proc(30)='f + g -> f + gamma*'
        proc(35)='f + gamma -> f + gamma*'
      elseif(mstp(43).eq.2) then
        proc(1)='f + fbar -> Z0'
        proc(15)='f + fbar -> g + Z0'
        proc(19)='f + fbar -> gamma + Z0'
        proc(30)='f + g -> f + Z0'
        proc(35)='f + gamma -> f + Z0'
      elseif(mstp(43).eq.3) then
        proc(1)='f + fbar -> gamma*/Z0'
        proc(15)='f + fbar -> g + gamma*/Z0'
        proc(19)='f + fbar -> gamma + gamma*/Z0'
        proc(30)='f + g -> f + gamma*/Z0'
        proc(35)='f + gamma -> f + gamma*/Z0'
      endif
 
C...Special cases in treatment of gamma*/Z0/Z'0: redefine process name.
      if(mstp(44).eq.1) then
        proc(141)='f + fbar -> gamma*'
      elseif(mstp(44).eq.2) then
        proc(141)='f + fbar -> Z0'
      elseif(mstp(44).eq.3) then
        proc(141)='f + fbar -> Z''0'
      elseif(mstp(44).eq.4) then
        proc(141)='f + fbar -> gamma*/Z0'
      elseif(mstp(44).eq.5) then
        proc(141)='f + fbar -> gamma*/Z''0'
      elseif(mstp(44).eq.6) then
        proc(141)='f + fbar -> Z0/Z''0'
      elseif(mstp(44).eq.7) then
        proc(141)='f + fbar -> gamma*/Z0/Z''0'
      endif
 
C...Special cases in treatment of WW -> WW: redefine process name.
      if(mstp(45).eq.1) then
        proc(77)='W+ + W+ -> W+ + W+'
      elseif(mstp(45).eq.2) then
        proc(77)='W+ + W- -> W+ + W-'
      elseif(mstp(45).eq.3) then
        proc(77)='W+/- + W+/- -> W+/- + W+/-'
      endif
 
C...Format for error information.
 5000 format(1x,'Error: unphysical input tan^2(beta) and m_H ',
     &'combination'/1x,'Execution stopped!')
 
      return
      end
 
C*********************************************************************
 
C...PYINBM
C...Identifies the two incoming particles and the choice of frame.
 
       subroutine pjinbm(chfram,kfp,kft,win)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      save /jyjets/,/jydat1/,/jydat2/,/pjsubs/,/pjpars/,/pjint1/
C...Local arrays, character variables and data.
      character chfram*8,chbeam*16,chtarg*16,chcom(3)*16,chalp(2)*26,
     &chidnt(3)*16,chtemp*8,chcde(29)*8,chinit*76
      dimension len(3),kcde(29),pm(2)
      data chalp/'abcdefghijklmnopqrstuvwxyz',
     &'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data chcde/'e-      ','e+      ','nu_e    ','nu_ebar ',
     &'mu-     ','mu+     ','nu_mu   ','nu_mubar','tau-    ',
     &'tau+    ','nu_tau  ','nu_tauba','pi+     ','pi-     ',
     &'n0      ','nbar0   ','p+      ','pbar-   ','gamma   ',
     &'lambda0 ','sigma-  ','sigma0  ','sigma+  ','xi-     ',
     &'xi0     ','omega-  ','pi0     ','reggeon ','pomeron '/
      data kcde/11,-11,12,-12,13,-13,14,-14,15,-15,16,-16,
     &211,-211,2112,-2112,2212,-2212,22,3122,3112,3212,3222,
     &3312,3322,3334,111,28,29/
 
C...Store initial energy. Default frame.
      vint(290)=win
      mint(111)=0
      mint(11)=kfp
      pm(1)=pjmass(mint(11))
      vint(3)=pm(1)
      mint(12)=kft
      pm(2)=pjmass(mint(12))
      vint(4)=pm(2)
      call pjname(mint(11),chbeam)
      call pjname(mint(12),chtarg)
      if(mint(11).eq.0) write(mstu(11),5000) chbeam
      if(mint(12).eq.0) write(mstu(11),5100) chtarg
      if(mint(11).eq.0.or.mint(12).eq.0) stop
      chcom(1)(1:8)=chfram
      chcom(2)=chbeam
      chcom(3)=chtarg

C...Convert character variables to lowercase and find their length.
      do i=1,3
        if(i.eq.1) then
          len(i)=8
          lm=8
        else
          len(i)=16
          lm=16
        endif
        do 110 ll=lm,1,-1
          if(len(i).eq.ll.and.chcom(i)(ll:ll).eq.' ') len(i)=ll-1
          do 100 la=1,26
            if(chcom(i)(ll:ll).eq.chalp(2)(la:la))
     $                            chcom(i)(ll:ll)=chalp(1)(la:la)
  100     continue
  110   continue
        chidnt(i)=chcom(i)
      end do

C...Fix up bar, underscore and charge in particle name (if needed).
c       do 120 ll=1,6
c         if(chidnt(i)(ll:ll).eq.'~') then
c           chtemp=chidnt(1)
c           chidnt(1)=chtemp(1:ll-1)//'bar'//chtemp(ll+1:6)//'  '
c         endif
c 120   continue
c       if(chidnt(i)(7:7).eq.'~') chidnt(i)(7:8)='ba'
c       if(chidnt(i)(1:2).eq.'nu'.and.chidnt(i)(3:3).ne.'_') then
c         chtemp=chidnt(i)
c         chidnt(i)='nu_'//chtemp(3:7)
c       elseif(chidnt(i)(1:2).eq.'n ') then
c         chidnt(i)(1:3)='n0 '
c       elseif(chidnt(i)(1:4).eq.'nbar') then
c         chidnt(i)(1:5)='nbar0'
c       elseif(chidnt(i)(1:2).eq.'p ') then
c         chidnt(i)(1:3)='p+ '
c       elseif(chidnt(i)(1:4).eq.'pbar'.or.
c    &    chidnt(i)(1:2).eq.'p-') then
c         chidnt(i)(1:5)='pbar-'
c       elseif(chidnt(i)(1:6).eq.'lambda') then
c         chidnt(i)(7:7)='0'
c       elseif(chidnt(i)(1:3).eq.'reg') then
c         chidnt(i)(1:7)='reggeon'
c       elseif(chidnt(i)(1:3).eq.'pom') then
c         chidnt(i)(1:7)='pomeron'
c       endif
c 130 continue
 
C...Identify free initialization.
      if(chcom(1)(1:2).eq.'no') then
        mint(65)=1
        return
      endif
 
C...Identify incoming beam and target particles.
c     do 150 i=1,2
c       do 140 j=1,29
c         if(chidnt(i+1).eq.chcde(j)) mint(10+i)=kcde(j)
c 140   continue
c       pm(i)=pjmass(mint(10+i))
c       vint(2+i)=pm(i)
c 150 continue

 
C...Identify choice of frame and input energies.
      chinit=' '
 
C...Events defined in the CM frame.
      if(chcom(1)(1:2).eq.'cm') then
        mint(111)=1
        s=win**2
        if(mstp(122).ge.1) then
          if(chcom(2)(1:1).ne.'e') then
            loffs=(31-(len(2)+len(3)))/2
            chinit(loffs+1:76)='PYTHIA will be initialized for a '//
     &      chcom(2)(1:len(2))//' on '//chcom(3)(1:len(3))//
     &      ' collider'//' '
          else
            loffs=(30-(len(2)+len(3)))/2
            chinit(loffs+1:76)='PYTHIA will be initialized for an '//
     &      chcom(2)(1:len(2))//' on '//chcom(3)(1:len(3))//
     &      ' collider'//' '
          endif
          write(mstu(11),5200) chinit
          write(mstu(11),5300) win
        endif
 
C...Events defined in fixed target frame.
      elseif(chcom(1)(1:3).eq.'fix') then
        mint(111)=2
        s=pm(1)**2+pm(2)**2+2d0*pm(2)*sqrt(pm(1)**2+win**2)
        if(mstp(122).ge.1) then
          loffs=(29-(len(2)+len(3)))/2
          chinit(loffs+1:76)='PYTHIA will be initialized for '//
     &    chcom(2)(1:len(2))//' on '//chcom(3)(1:len(3))//
     &    ' fixed target'//' '
          write(mstu(11),5200) chinit
          write(mstu(11),5400) win
          write(mstu(11),5500) sqrt(s)
        endif
 
C...Frame defined by user three-vectors.
      elseif(chcom(1)(1:3).eq.'use') then
        mint(111)=3
        p(1,5)=pm(1)
        p(2,5)=pm(2)
        p(1,4)=sqrt(p(1,1)**2+p(1,2)**2+p(1,3)**2+p(1,5)**2)
        p(2,4)=sqrt(p(2,1)**2+p(2,2)**2+p(2,3)**2+p(2,5)**2)
        s=(p(1,4)+p(2,4))**2-(p(1,1)+p(2,1))**2-(p(1,2)+p(2,2))**2-
     &  (p(1,3)+p(2,3))**2
        if(mstp(122).ge.1) then
          loffs=(12-(len(2)+len(3)))/2
          chinit(loffs+1:76)='PYTHIA will be initialized for '//
     &    chcom(2)(1:len(2))//' on '//chcom(3)(1:len(3))//
     &    ' user-specified configuration'//' '
          write(mstu(11),5200) chinit
          write(mstu(11),5600)
          write(mstu(11),5700) chcom(2),p(1,1),p(1,2),p(1,3),p(1,4)
          write(mstu(11),5700) chcom(3),p(2,1),p(2,2),p(2,3),p(2,4)
          write(mstu(11),5500) sqrt(max(0d0,s))
        endif
 
C...Frame defined by user four-vectors.
      elseif(chcom(1)(1:4).eq.'four') then
        mint(111)=4
        pms1=p(1,4)**2-p(1,1)**2-p(1,2)**2-p(1,3)**2
        p(1,5)=sign(sqrt(abs(pms1)),pms1)
        pms2=p(2,4)**2-p(2,1)**2-p(2,2)**2-p(2,3)**2
        p(2,5)=sign(sqrt(abs(pms2)),pms2)
        s=(p(1,4)+p(2,4))**2-(p(1,1)+p(2,1))**2-(p(1,2)+p(2,2))**2-
     &  (p(1,3)+p(2,3))**2
        if(mstp(122).ge.1) then
          loffs=(12-(len(2)+len(3)))/2
          chinit(loffs+1:76)='PYTHIA will be initialized for '//
     &    chcom(2)(1:len(2))//' on '//chcom(3)(1:len(3))//
     &    ' user-specified configuration'//' '
          write(mstu(11),5200) chinit
          write(mstu(11),5600)
          write(mstu(11),5700) chcom(2),p(1,1),p(1,2),p(1,3),p(1,4)
          write(mstu(11),5700) chcom(3),p(2,1),p(2,2),p(2,3),p(2,4)
          write(mstu(11),5500) sqrt(max(0d0,s))
        endif
 
C...Frame defined by user five-vectors.
      elseif(chcom(1)(1:4).eq.'five') then
        mint(111)=5
        s=(p(1,4)+p(2,4))**2-(p(1,1)+p(2,1))**2-(p(1,2)+p(2,2))**2-
     &  (p(1,3)+p(2,3))**2
        if(mstp(122).ge.1) then
ch        loffs=(12-(len(2)+len(3)))/2
          loffs=(30-(len(2)+len(3)))/2
          chinit(loffs+1:76)='PYTHIA will be initialized for '//
     &    chcom(2)(1:len(2))//' on '//chcom(3)(1:len(3))//
     &    ' user-specified configuration'//' '
          write(mstu(11),5200) chinit
          write(mstu(11),5600)
          write(mstu(11),5700) chcom(2),p(1,1),p(1,2),p(1,3),p(1,4)
          write(mstu(11),5700) chcom(3),p(2,1),p(2,2),p(2,3),p(2,4)
          write(mstu(11),5500) sqrt(max(0d0,s))
        endif
 
C...Unknown frame.
      else
        write(mstu(11),5800) chfram(1:len(1))
        stop
      endif

c...Error for too low CM energy.
      if(s.lt.parp(2)**2) then
        write(mstu(11),5900) sqrt(s)
        stop
      endif
 
C...Formats for initialization and error information.
 5000 format(1x,'Error: unrecognized beam particle ''',a,'''D0'/
     &1x,'Execution stopped!')
 5100 format(1x,'Error: unrecognized target particle ''',a,'''D0'/
     &1x,'Execution stopped!')

 5200 format(/1x,78('=')/1x,'I',76x,'I'/1x,'I',a76,'I')
 5300 format(1x,'I',18x,'at',1x,f10.3,1x,'GeV center-of-mass energy',
     &19x,'I'/1x,'I',76x,'I'/1x,78('='))
 5400 format(1x,'I',22x,'at',1x,f10.3,1x,'GeV/c lab-momentum',22x,'I')
 5500 format(1x,'I',76x,'I'/1x,'I',11x,'corresponding to',1x,f10.3,1x,
     &'GeV center-of-mass energy',12x,'I'/1x,'I',76x,'I'/1x,78('='))
 5600 format(1x,'I',76x,'I'/1x,'I',18x,'px (GeV/c)',3x,'py (GeV/c)',3x,
     &'pz (GeV/c)',6x,'E (GeV)',9x,'I')
 5700 format(1x,'I',8x,a8,4(2x,f10.3,1x),8x,'I')
 5800 format(1x,'Error: unrecognized coordinate frame ''',a,'''D0'/
     &1x,'Execution stopped!')
 5900 format(1x,'Error: too low CM energy,',f8.3,' GeV for event ',
     &'generation.'/1x,'Execution stopped!')
 
      return
      end
 
C*********************************************************************
 
C...PYINKI
C...Sets up kinematics, including rotations and boosts to/from CM frame.
 
      subroutine pjinki(modki)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      save /jyjets/,/jydat1/,/jydat2/,/pjsubs/,/pjpars/,/pjint1/
 
C...Set initial flavour state.
      n=2
      do 100 i=1,2
        k(i,1)=1
        k(i,2)=mint(10+i)
  100 continue
 
C...Reset boost. Do kinematics for various cases.
      do 110 j=6,10
        vint(j)=0d0
  110 continue
 
C...Set up kinematics for events defined in CM frame.
      if(mint(111).eq.1) then
        win=vint(290)
        if(modki.eq.1) win=parp(171)*vint(290)
        s=win**2
        p(1,5)=vint(3)
        p(2,5)=vint(4)
        p(1,1)=0d0
        p(1,2)=0d0
        p(2,1)=0d0
        p(2,2)=0d0
        p(1,3)=sqrt(((s-p(1,5)**2-p(2,5)**2)**2-(2d0*p(1,5)*p(2,5))**2)/
     &  (4d0*s))
        p(2,3)=-p(1,3)
        p(1,4)=sqrt(p(1,3)**2+p(1,5)**2)
        p(2,4)=sqrt(p(2,3)**2+p(2,5)**2)
 
C...Set up kinematics for fixed target events.
      elseif(mint(111).eq.2) then
        win=vint(290)
        if(modki.eq.1) win=parp(171)*vint(290)
        p(1,5)=vint(3)
        p(2,5)=vint(4)
        p(1,1)=0d0
        p(1,2)=0d0
        p(2,1)=0d0
        p(2,2)=0d0
        p(1,3)=win
        p(1,4)=sqrt(p(1,3)**2+p(1,5)**2)
        p(2,3)=0d0
        p(2,4)=p(2,5)
        s=p(1,5)**2+p(2,5)**2+2d0*p(2,4)*p(1,4)
        vint(10)=p(1,3)/(p(1,4)+p(2,4))
        call pjrobo(0,0,0d0,0d0,0d0,0d0,-vint(10))
 
C...Set up kinematics for events in user-defined frame.
      elseif(mint(111).eq.3) then
        p(1,5)=vint(3)
        p(2,5)=vint(4)
        p(1,4)=sqrt(p(1,1)**2+p(1,2)**2+p(1,3)**2+p(1,5)**2)
        p(2,4)=sqrt(p(2,1)**2+p(2,2)**2+p(2,3)**2+p(2,5)**2)
        do 120 j=1,3
          vint(7+j)=(p(1,j)+p(2,j))/(p(1,4)+p(2,4))
  120   continue
        call pjrobo(0,0,0d0,0d0,-vint(8),-vint(9),-vint(10))
        vint(7)=pjangl(p(1,1),p(1,2))
        call pjrobo(0,0,0d0,-vint(7),0d0,0d0,0d0)
        vint(6)=pjangl(p(1,3),p(1,1))
        call pjrobo(0,0,-vint(6),0d0,0d0,0d0,0d0)
        s=p(1,5)**2+p(2,5)**2+2d0*(p(1,4)*p(2,4)-p(1,3)*p(2,3))
 
C...Set up kinematics for events with user-defined four-vectors.
      elseif(mint(111).eq.4) then
        pms1=p(1,4)**2-p(1,1)**2-p(1,2)**2-p(1,3)**2
        p(1,5)=sign(sqrt(abs(pms1)),pms1)
        pms2=p(2,4)**2-p(2,1)**2-p(2,2)**2-p(2,3)**2
        p(2,5)=sign(sqrt(abs(pms2)),pms2)
        do 130 j=1,3
          vint(7+j)=(p(1,j)+p(2,j))/(p(1,4)+p(2,4))
  130   continue
        call pjrobo(0,0,0d0,0d0,-vint(8),-vint(9),-vint(10))
        vint(7)=pjangl(p(1,1),p(1,2))
        call pjrobo(0,0,0d0,-vint(7),0d0,0d0,0d0)
        vint(6)=pjangl(p(1,3),p(1,1))
        call pjrobo(0,0,-vint(6),0d0,0d0,0d0,0d0)
        s=(p(1,4)+p(2,4))**2
 
C...Set up kinematics for events with user-defined five-vectors.
      elseif(mint(111).eq.5) then
        do 140 j=1,3
          vint(7+j)=(p(1,j)+p(2,j))/(p(1,4)+p(2,4))
  140   continue
        call pjrobo(0,0,0d0,0d0,-vint(8),-vint(9),-vint(10))
        vint(7)=pjangl(p(1,1),p(1,2))
        call pjrobo(0,0,0d0,-vint(7),0d0,0d0,0d0)
        vint(6)=pjangl(p(1,3),p(1,1))
        call pjrobo(0,0,-vint(6),0d0,0d0,0d0,0d0)
        s=(p(1,4)+p(2,4))**2
      endif
 
C...Return or error for too low CM energy.
      if(modki.eq.1.and.s.lt.parp(2)**2) then
        if(mstp(172).le.1) then
          call pjerrm(23,
     &    '(PYINKI:) too low invariant mass in this event')
        else
          msti(61)=1
          return
        endif
      endif
 
C...Save information on incoming particles.
      vint(1)=sqrt(s)
      vint(2)=s
      if(mint(111).ge.4) vint(3)=p(1,5)
      if(mint(111).ge.4) vint(4)=p(2,5)
      vint(5)=p(1,3)
      if(modki.eq.0) vint(289)=s
      do 150 j=1,5
        v(1,j)=0d0
        v(2,j)=0d0
        vint(290+j)=p(1,j)
        vint(295+j)=p(2,j)
  150 continue
 
C...Store pT cut-off and related constants to be used in generation.
      if(modki.eq.0) vint(285)=ckin(3)
      if(mstp(82).le.1) then
        if(mint(121).gt.1) parp(81)=1.30d0+0.15d0*log(vint(1)/200d0)/
     &  log(900d0/200d0)
        ptmn=parp(81)
      else
        if(mint(121).gt.1) parp(82)=1.25d0+0.15d0*log(vint(1)/200d0)/
     &  log(900d0/200d0)
        ptmn=parp(82)
      endif
      vint(149)=4d0*ptmn**2/s
 
      return
      end
 
C*********************************************************************
 
C...PYINPR
C...Selects partonic subprocesses to be included in the simulation.
 
      subroutine pjinpr
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      save /jydat1/,/jydat3/,/pjsubs/,/pjpars/,/pjint1/,/pjint2/
 
C...Reset processes to be included.
      if(msel.ne.0) then
        do 100 i=1,500
          msub(i)=0
  100   continue
      endif
 
C...For e-gamma witn MSTP(14)=10 allow mixture of VMD and anomalous.
      if(mint(121).eq.2) then
        msub(10)=1
        mint(123)=mint(122)+1
 
C...For gamma-p or gamma-gamma with MSTP(14)=10 allow mixture.
C...Here also set a few parameters otherwise normally not touched.
      elseif(mint(121).gt.1) then
 
C...Parton distributions dampened at small Q2; go to low energies,
C...alpha_s <1; no minimum pT cut-off a priori.
        mstp(57)=3
        mstp(85)=0
        parp(2)=2d0
        paru(115)=1d0
        ckin(5)=0.2d0
        ckin(6)=0.2d0
 
C...Define pT cut-off parameters and whether run involves low-pT.
        if(mstp(82).le.1) then
          ptmvmd=1.30d0+0.15d0*log(vint(1)/200d0)/log(900d0/200d0)
        else
          ptmvmd=1.25d0+0.15d0*log(vint(1)/200d0)/log(900d0/200d0)
        endif
        ptmdir=parp(15)
        ptmano=ptmvmd
        if(mstp(15).eq.5) ptmano=0.60d0+
     &  0.125d0*log(1d0+0.10d0*vint(1))**2
        iptl=1
        if(vint(285).gt.max(ptmvmd,ptmdir,ptmano)) iptl=0
        if(msel.eq.2) iptl=1
 
C...Set up for p/VMD * VMD.
        if(mint(122).eq.1) then
          mint(123)=2
          msub(11)=1
          msub(12)=1
          msub(13)=1
          msub(28)=1
          msub(53)=1
          msub(68)=1
          if(iptl.eq.1) msub(95)=1
          if(msel.eq.2) then
            msub(91)=1
            msub(92)=1
            msub(93)=1
            msub(94)=1
          endif
          parp(81)=ptmvmd
          parp(82)=ptmvmd
          if(iptl.eq.1) ckin(3)=0d0
 
C...Set up for p/VMD * direct gamma.
        elseif(mint(122).eq.2) then
          mint(123)=0
          if(mint(121).eq.6) mint(123)=5
          msub(33)=1
          msub(54)=1
          if(iptl.eq.1) ckin(3)=ptmdir
 
C...Set up for p/VMD * anomalous gamma.
        elseif(mint(122).eq.3) then
          mint(123)=3
          if(mint(121).eq.6) mint(123)=7
          msub(11)=1
          msub(12)=1
          msub(13)=1
          msub(28)=1
          msub(53)=1
          msub(68)=1
          if(mstp(82).ge.2) mstp(85)=1
          if(iptl.eq.1) ckin(3)=ptmano
 
C...Set up for direct * direct gamma (switch off leptons).
        elseif(mint(122).eq.4) then
          mint(123)=0
          msub(58)=1
          do 110 ii=mdcy(22,2),mdcy(22,2)+mdcy(22,3)-1
            if(iabs(kfdp(ii,1)).ge.10) mdme(ii,1)=min(0,mdme(ii,1))
  110     continue
          if(iptl.eq.1) ckin(3)=ptmdir
 
C...Set up for direct * anomalous gamma.
        elseif(mint(122).eq.5) then
          mint(123)=6
          msub(33)=1
          msub(54)=1
          if(iptl.eq.1) ckin(3)=ptmano
 
C...Set up for anomalous * anomalous gamma.
        elseif(mint(122).eq.6) then
          mint(123)=3
          msub(11)=1
          msub(12)=1
          msub(13)=1
          msub(28)=1
          msub(53)=1
          msub(68)=1
          if(mstp(82).ge.2) mstp(85)=1
          if(iptl.eq.1) ckin(3)=ptmano
        endif
 
C...End of special set up for gamma-p and gamma-gamma.
        ckin(1)=2d0*ckin(3)
      endif
 
C...Flavour information for individual beams.
      do 120 i=1,2
        mint(40+i)=1
        if(mint(123).ge.1.and.mint(10+i).eq.22) mint(40+i)=2
        if(iabs(mint(10+i)).gt.100) mint(40+i)=2
        if(mint(10+i).eq.28.or.mint(10+i).eq.29) mint(40+i)=2
        mint(44+i)=mint(40+i)
        if(mstp(11).ge.1.and.iabs(mint(10+i)).eq.11) mint(44+i)=3
  120 continue
 
C...If two gammas, whereof one direct, pick the first.
      if(mint(11).eq.22.and.mint(12).eq.22) then
        if(mint(123).ge.4.and.mint(123).le.6) then
          mint(41)=1
          mint(45)=1
        endif
      elseif(mint(11).eq.22.or.mint(12).eq.22) then
        if(mint(123).ge.4) call pjerrm(26,
     &  '(PYINPR:) unallowed MSTP(14) code for single photon')
      endif
 
C...Flavour information on combination of incoming particles.
      mint(43)=2*mint(41)+mint(42)-2
      mint(44)=mint(43)
      if(mint(123).le.0) then
        if(mint(11).eq.22) mint(43)=mint(43)+2
        if(mint(12).eq.22) mint(43)=mint(43)+1
      elseif(mint(123).le.3) then
        if(mint(11).eq.22) mint(44)=mint(44)-2
        if(mint(12).eq.22) mint(44)=mint(44)-1
      elseif(mint(11).eq.22.and.mint(12).eq.22) then
        mint(43)=4
        mint(44)=1
      endif
      mint(47)=2*min(2,mint(45))+min(2,mint(46))-2
      if(min(mint(45),mint(46)).eq.3) mint(47)=5
      mint(50)=0
      if(mint(41).eq.2.and.mint(42).eq.2) mint(50)=1
      if((mint(11).eq.22.or.mint(12).eq.22).and.mint(123).ge.3)
     &mint(50)=0
      mint(107)=0
      if(mint(11).eq.22) then
        mint(107)=mint(123)
        if(mint(123).ge.4) mint(107)=0
        if(mint(123).eq.7) mint(107)=2
      endif
      mint(108)=0
      if(mint(12).eq.22) then
        mint(108)=mint(123)
        if(mint(123).ge.4) mint(108)=mint(123)-3
        if(mint(123).eq.7) mint(108)=3
      endif
 
C...Select default processes according to incoming beams
C...(already done for gamma-p and gamma-gamma with MSTP(14)=10).
      if(mint(121).gt.1) then
      elseif(msel.eq.1.or.msel.eq.2) then
 
        if(mint(43).eq.1) then
C...Lepton + lepton -> gamma/Z0 or W.
          if(mint(11)+mint(12).eq.0) msub(1)=1
          if(mint(11)+mint(12).ne.0) msub(2)=1
 
        elseif(mint(43).le.3.and.mint(123).eq.0.and.
     &    (mint(11).eq.22.or.mint(12).eq.22)) then
C...Unresolved photon + lepton: Compton scattering.
          msub(34)=1
 
        elseif(mint(43).le.3) then
C...Lepton + hadron: deep inelastic scattering.
          msub(10)=1
 
        elseif(mint(123).eq.0.and.mint(11).eq.22.and.
     &    mint(12).eq.22) then
C...Two unresolved photons: fermion pair production.
          msub(58)=1
 
        elseif((mint(123).eq.0.and.(mint(11).eq.22.or.mint(12).eq.22))
     &    .or.(mint(123).ge.4.and.mint(123).le.6.and.mint(11).eq.22.and.
     &    mint(12).eq.22)) then
C...Unresolved photon + hadron: photon-parton scattering.
          msub(33)=1
          msub(34)=1
          msub(54)=1
 
        elseif(msel.eq.1) then
C...High-pT QCD processes:
          msub(11)=1
          msub(12)=1
          msub(13)=1
          msub(28)=1
          msub(53)=1
          msub(68)=1
          if(mstp(82).le.1.and.ckin(3).lt.parp(81)) msub(95)=1
          if(mstp(82).ge.2.and.ckin(3).lt.parp(82)) msub(95)=1
          if(msub(95).eq.1.and.mint(50).eq.0) msub(95)=0
 
        else
C...All QCD processes:
          msub(11)=1
          msub(12)=1
          msub(13)=1
          msub(28)=1
          msub(53)=1
          msub(68)=1
          msub(91)=1
          msub(92)=1
          msub(93)=1
          msub(94)=1
          msub(95)=1
        endif
 
      elseif(msel.ge.4.and.msel.le.8) then
C...Heavy quark production.
        msub(81)=1
        msub(82)=1
        msub(84)=1
        do 130 j=1,min(8,mdcy(21,3))
          mdme(mdcy(21,2)+j-1,1)=0
  130   continue
        mdme(mdcy(21,2)+msel-1,1)=1
        msub(85)=1
        do 140 j=1,min(12,mdcy(22,3))
          mdme(mdcy(22,2)+j-1,1)=0
  140   continue
        mdme(mdcy(22,2)+msel-1,1)=1
 
      elseif(msel.eq.10) then
C...Prompt photon production:
        msub(14)=1
        msub(18)=1
        msub(29)=1
 
      elseif(msel.eq.11) then
C...Z0/gamma* production:
        msub(1)=1
 
      elseif(msel.eq.12) then
C...W+/- production:
        msub(2)=1
 
      elseif(msel.eq.13) then
C...Z0 + jet:
        msub(15)=1
        msub(30)=1
 
      elseif(msel.eq.14) then
C...W+/- + jet:
        msub(16)=1
        msub(31)=1
 
      elseif(msel.eq.15) then
C...Z0 & W+/- pair production:
        msub(19)=1
        msub(20)=1
        msub(22)=1
        msub(23)=1
        msub(25)=1
 
      elseif(msel.eq.16) then
C...h0 production:
        msub(3)=1
        msub(102)=1
        msub(103)=1
        msub(123)=1
        msub(124)=1
 
      elseif(msel.eq.17) then
C...h0 & Z0 or W+/- pair production:
        msub(24)=1
        msub(26)=1
 
      elseif(msel.eq.18) then
C...h0 production; interesting processes in e+e-.
        msub(24)=1
        msub(103)=1
        msub(123)=1
        msub(124)=1
 
      elseif(msel.eq.19) then
C...h0, H0 and A0 production; interesting processes in e+e-.
        msub(24)=1
        msub(103)=1
        msub(123)=1
        msub(124)=1
        msub(153)=1
        msub(171)=1
        msub(173)=1
        msub(174)=1
        msub(158)=1
        msub(176)=1
        msub(178)=1
        msub(179)=1
 
      elseif(msel.eq.21) then
C...Z'0 production:
        msub(141)=1
 
      elseif(msel.eq.22) then
C...W'+/- production:
        msub(142)=1
 
      elseif(msel.eq.23) then
C...H+/- production:
        msub(143)=1
 
      elseif(msel.eq.24) then
C...R production:
        msub(144)=1
 
      elseif(msel.eq.25) then
C...LQ (leptoquark) production.
        msub(145)=1
        msub(162)=1
        msub(163)=1
        msub(164)=1
 
      elseif(msel.ge.35.and.msel.le.38) then
C...Production of one heavy quark (W exchange):
        msub(83)=1
        do 150 j=1,min(8,mdcy(21,3))
          mdme(mdcy(21,2)+j-1,1)=0
  150   continue
        mdme(mdcy(21,2)+msel-31,1)=1
 
CMRENNA++Define SUSY alternatives.
      elseif(msel.eq.39) then
C...Turn on all SUSY processes.
        if(mint(43).eq.4) then
C...Hadron-hadron processes.
          do 160 i=201,280
            if(iset(i).ge.0) msub(i)=1
  160     continue
        elseif(mint(43).eq.1) then
C...Lepton-lepton processes: QED production of squarks.
          do 170 i=201,214
            msub(i)=1
  170     continue
          msub(210)=0
          msub(211)=0
          msub(212)=0
          do 180 i=216,228
            msub(i)=1
  180     continue
          do 190 i=261,263
            msub(i)=1
  190     continue
          msub(277)=1
          msub(278)=1
        endif
 
      elseif(msel.eq.40) then
C...Gluinos and squarks.
        if(mint(43).eq.4) then
          msub(243)=1
          msub(244)=1
          msub(258)=1
          msub(259)=1
          msub(261)=1
          msub(262)=1
          msub(264)=1
          msub(265)=1
          do 200 i=271,280
            msub(i)=1
  200     continue
        elseif(mint(43).eq.1) then
          msub(277)=1
          msub(278)=1
        endif
 
      elseif(msel.eq.41) then
C...Stop production.
        msub(261)=1
        msub(262)=1
        msub(263)=1
        if(mint(43).eq.4) then
          msub(264)=1
          msub(265)=1
        endif
 
      elseif(msel.eq.42) then
C...Slepton production.
        do 210 i=201,214
          msub(i)=1
  210   continue
        if(mint(43).ne.4) then
          msub(210)=0
          msub(211)=0
          msub(212)=0
        endif
 
      elseif(msel.eq.43) then
C...Neutralino/Chargino + Gluino/Squark.
        if(mint(43).eq.4) then
          do 220 i=237,242
            msub(i)=1
  220     continue
          do 230 i=246,257
            msub(i)=1
  230     continue
        endif
 
      elseif(msel.eq.44) then
C...Neutralino/Chargino pair production.
        if(mint(43).eq.4) then
          do 240 i=216,236
            msub(i)=1
  240     continue
        elseif(mint(43).eq.1) then
          do 250 i=216,228
            msub(i)=1
  250     continue
        endif
      endif
 
C...Find heaviest new quark flavour allowed in processes 81-84.
      kflqm=1
      do 260 i=1,min(8,mdcy(21,3))
        idc=i+mdcy(21,2)-1
        if(mdme(idc,1).le.0) goto 260
        kflqm=i
  260 continue
      if(mstp(7).ge.1.and.mstp(7).le.8.and.(msel.le.3.or.msel.ge.9))
     &kflqm=mstp(7)
      mint(55)=kflqm
      kfpr(81,1)=kflqm
      kfpr(81,2)=kflqm
      kfpr(82,1)=kflqm
      kfpr(82,2)=kflqm
      kfpr(83,1)=kflqm
      kfpr(84,1)=kflqm
      kfpr(84,2)=kflqm
 
C...Find heaviest new fermion flavour allowed in process 85.
      kflfm=1
      do 270 i=1,min(12,mdcy(22,3))
        idc=i+mdcy(22,2)-1
        if(mdme(idc,1).le.0) goto 270
        kflfm=kfdp(idc,1)
  270 continue
      if(((mstp(7).ge.1.and.mstp(7).le.8).or.(mstp(7).ge.11.and.
     &mstp(7).le.18)).and.(msel.le.3.or.msel.ge.9)) kflfm=mstp(7)
      mint(56)=kflfm
      kfpr(85,1)=kflfm
      kfpr(85,2)=kflfm
 
      return
      end
 
C*********************************************************************
 
C...PYXTOT
C...Parametrizes total, elastic and diffractive cross-sections
C...for different energies and beams. Donnachie-Landshoff for
C...total and Schuler-Sjostrand for elastic and diffractive.
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
 
      subroutine pjxtot
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint5/ngenpd,ngen(0:500,3),xsec(0:500,3)
      common/pjint7/sigt(0:6,0:6,0:5)
      save /jydat1/,/pjpars/,/pjint1/,/pjint5/,/pjint7/
C...Local arrays.
      dimension nproc(30),xpar(30),ypar(30),ihada(20),ihadb(20),
     &pmhad(4),bhad(4),betp(4),ifitsd(20),ifitdd(20),ceffs(10,8),
     &ceffd(10,9),sigtmp(6,0:5)
 
C...Common constants.
      data eps/0.0808d0/, eta/-0.4525d0/, alp/0.25d0/, cres/2d0/,
     &pmrc/1.062d0/, smp/0.880d0/, facel/0.0511d0/, facsd/0.0336d0/,
     &facdd/0.0084d0/
 
C...Number of multiple processes to be evaluated (= 0 : undefined).
      data nproc/7*1,3*0,6*1,4*0,4*3,2*6,4*0/
C...X and Y parameters of sigmatot = X * s**epsilon + Y * s**(-eta).
      data xpar/2*21.70d0,3*13.63d0,10.01d0,0.970d0,3*0d0,
     &8.56d0,6.29d0,0.609d0,4.62d0,0.447d0,0.0434d0,4*0d0,
     &0.0677d0,0.0534d0,0.0425d0,0.0335d0,2.11d-4,1.31d-4,4*0d0/
      data ypar/
     &56.08d0,98.39d0,27.56d0,36.02d0,31.79d0,-1.51d0,-0.146d0,3*0d0,
     &13.08d0,-0.62d0,-0.060d0,0.030d0,-0.0028d0,0.00028d0,4*0d0,
     &0.129d0,0.115d0,0.081d0,0.072d0,2.15d-4,1.70d-4,4*0d0/
 
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
      data ((ceffs(j1,j2),j2=1,8),j1=1,10)/
     &0.213d0, 0.0d0, -0.47d0, 150d0, 0.213d0, 0.0d0, -0.47d0, 150d0,
     &0.213d0, 0.0d0, -0.47d0, 150d0, 0.267d0, 0.0d0, -0.47d0, 100d0,
     &0.213d0, 0.0d0, -0.47d0, 150d0, 0.232d0, 0.0d0, -0.47d0, 110d0,
     &0.213d0, 7.0d0, -0.55d0, 800d0, 0.115d0, 0.0d0, -0.47d0, 110d0,
     &0.267d0, 0.0d0, -0.46d0,  75d0, 0.267d0, 0.0d0, -0.46d0,  75d0,
     &0.232d0, 0.0d0, -0.46d0,  85d0, 0.267d0, 0.0d0, -0.48d0, 100d0,
     &0.115d0, 0.0d0, -0.50d0,  90d0, 0.267d0, 6.0d0, -0.56d0, 420d0,
     &0.232d0, 0.0d0, -0.48d0, 110d0, 0.232d0, 0.0d0, -0.48d0, 110d0,
     &0.115d0, 0.0d0, -0.52d0, 120d0, 0.232d0, 6.0d0, -0.56d0, 470d0,
     &0.115d0, 5.5d0, -0.58d0, 570d0, 0.115d0, 5.5d0, -0.58d0, 570d0/
      data ((ceffd(j1,j2),j2=1,9),j1=1,10)/
     &3.11d0, -7.34d0,  9.71d0, 0.068d0, -0.42d0,  1.31d0,
     &-1.37d0,  35.0d0,  118d0,  3.11d0, -7.10d0,  10.6d0,
     &0.073d0, -0.41d0, 1.17d0, -1.41d0,  31.6d0,   95d0,
     &3.12d0, -7.43d0,  9.21d0, 0.067d0, -0.44d0,  1.41d0,
     &-1.35d0,  36.5d0,  132d0,  3.13d0, -8.18d0, -4.20d0,
     &0.056d0, -0.71d0, 3.12d0, -1.12d0,  55.2d0, 1298d0,
     &3.11d0, -6.90d0,  11.4d0, 0.078d0, -0.40d0,  1.05d0,
     &-1.40d0,  28.4d0,   78d0,  3.11d0, -7.13d0,  10.0d0,
     &0.071d0, -0.41d0, 1.23d0, -1.34d0,  33.1d0,  105d0,
     &3.12d0, -7.90d0, -1.49d0, 0.054d0, -0.64d0,  2.72d0,
     &-1.13d0,  53.1d0,  995d0,  3.11d0, -7.39d0,  8.22d0,
     &0.065d0, -0.44d0, 1.45d0, -1.36d0,  38.1d0,  148d0,
     &3.18d0, -8.95d0, -3.37d0, 0.057d0, -0.76d0,  3.32d0,
     &-1.12d0,  55.6d0, 1472d0,  4.18d0, -29.2d0,  56.2d0,
     &0.074d0, -1.36d0, 6.67d0, -1.14d0, 116.2d0, 6532d0/
 
C...Parameters. Combinations of the energy.
      aem=paru(101)
      pmth=parp(102)
      s=vint(2)
      srt=vint(1)
      seps=s**eps
      seta=s**eta
      slog=log(s)
 
C...Ratio of gamma/pi (for rescaling in parton distributions).
      vint(281)=(xpar(22)*seps+ypar(22)*seta)/
     &(xpar(5)*seps+ypar(5)*seta)
      if(mint(50).ne.1) return
 
C...Order flavours of incoming particles: KF1 < KF2.
      if(iabs(mint(11)).le.iabs(mint(12))) then
        kf1=iabs(mint(11))
        kf2=iabs(mint(12))
        iord=1
      else
        kf1=iabs(mint(12))
        kf2=iabs(mint(11))
        iord=2
      endif
      isgn12=isign(1,mint(11)*mint(12))
 
C...Find process number (for lookup tables).
      if(kf1.gt.1000) then
        iproc=1
        if(isgn12.lt.0) iproc=2
      elseif(kf1.gt.100.and.kf2.gt.1000) then
        iproc=3
        if(isgn12.lt.0) iproc=4
        if(kf1.eq.111) iproc=5
      elseif(kf1.gt.100) then
        iproc=11
      elseif(kf2.gt.1000) then
        iproc=21
        if(mint(123).eq.2) iproc=22
      elseif(kf2.gt.100) then
        iproc=23
        if(mint(123).eq.2) iproc=24
      else
        iproc=25
        if(mint(123).eq.2) iproc=26
      endif
 
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
      if(mstp(31).le.0) return
      if(npr.eq.0) call pjerrm(26,
     &'(PYXTOT:) cross section for this process not yet implemented')
 
C...Parameters. Combinations of the energy.
      aem=paru(101)
      pmth=parp(102)
      s=vint(2)
      srt=vint(1)
      seps=s**eps
      seta=s**eta
      slog=log(s)
 
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
          sigtmp(i,j)=0d0
  100   continue
        if(srt.lt.pma+pmb+parp(104)) goto 110
 
C...Total cross-section. Elastic slope parameter and cross-section.
        sigtmp(i,0)=xpar(ipr)*seps+ypar(ipr)*seta
        bel=2d0*bha+2d0*bhb+4d0*seps-4.2d0
        sigtmp(i,1)=facel*sigtmp(i,0)**2/bel
 
C...Diffractive scattering A + B -> X + B.
        bsd=2d0*bhb
        sqml=(pma+pmth)**2
        sqmu=s*ceffs(isd,1)+ceffs(isd,2)
        sum1=log((bsd+2d0*alp*log(s/sqml))/
     &  (bsd+2d0*alp*log(s/sqmu)))/(2d0*alp)
        bxb=ceffs(isd,3)+ceffs(isd,4)/s
        sum2=cres*log(1d0+((pma+pmrc)/(pma+pmth))**2)/
     &  (bsd+2d0*alp*log(s/((pma+pmth)*(pma+pmrc)))+bxb)
        sigtmp(i,2)=facsd*xpar(ipr)*betp(ihb)*max(0d0,sum1+sum2)
 
C...Diffractive scattering A + B -> A + X.
        bsd=2d0*bha
        sqml=(pmb+pmth)**2
        sqmu=s*ceffs(isd,5)+ceffs(isd,6)
        sum1=log((bsd+2d0*alp*log(s/sqml))/
     &  (bsd+2d0*alp*log(s/sqmu)))/(2d0*alp)
        bax=ceffs(isd,7)+ceffs(isd,8)/s
        sum2=cres*log(1d0+((pmb+pmrc)/(pmb+pmth))**2)/
     &  (bsd+2d0*alp*log(s/((pmb+pmth)*(pmb+pmrc)))+bax)
        sigtmp(i,3)=facsd*xpar(ipr)*betp(iha)*max(0d0,sum1+sum2)
 
C...Order single diffractive correctly.
        if(iord.eq.2) then
          sigsav=sigtmp(i,2)
          sigtmp(i,2)=sigtmp(i,3)
          sigtmp(i,3)=sigsav
        endif
 
C...Double diffractive scattering A + B -> X1 + X2.
        yeff=log(s*smp/((pma+pmth)*(pmb+pmth))**2)
        deff=ceffd(idd,1)+ceffd(idd,2)/slog+ceffd(idd,3)/slog**2
        sum1=deff+yeff*(log(max(1d-10,yeff/deff))-1d0)/(2d0*alp)
        if(yeff.le.0) sum1=0d0
        sqmu=s*(ceffd(idd,4)+ceffd(idd,5)/slog+ceffd(idd,6)/slog**2)
        slup=log(max(1.1d0,s/(alp*(pma+pmth)**2*(pmb+pmth)*(pmb+pmrc))))
        sldn=log(max(1.1d0,s/(alp*sqmu*(pmb+pmth)*(pmb+pmrc))))
        sum2=cres*log(1d0+((pmb+pmrc)/(pmb+pmth))**2)*log(slup/sldn)/
     &  (2d0*alp)
        slup=log(max(1.1d0,s/(alp*(pmb+pmth)**2*(pma+pmth)*(pma+pmrc))))
        sldn=log(max(1.1d0,s/(alp*sqmu*(pma+pmth)*(pma+pmrc))))
        sum3=cres*log(1d0+((pma+pmrc)/(pma+pmth))**2)*log(slup/sldn)/
     &  (2d0*alp)
        bxx=ceffd(idd,7)+ceffd(idd,8)/srt+ceffd(idd,9)/s
        slrr=log(s/(alp*(pma+pmth)*(pma+pmrc)*(pmb+pmth)*(pmb*pmrc)))
        sum4=cres**2*log(1d0+((pma+pmrc)/(pma+pmth))**2)*
     &  log(1d0+((pmb+pmrc)/(pmb+pmth))**2)/max(0.1d0,2d0*alp*slrr+bxx)
        sigtmp(i,4)=facdd*xpar(ipr)*max(0d0,sum1+sum2+sum3+sum4)
 
C...Non-diffractive by unitarity.
        sigtmp(i,5)=sigtmp(i,0)-sigtmp(i,1)-sigtmp(i,2)-sigtmp(i,3)-
     &  sigtmp(i,4)
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
  130     continue
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
  160     continue
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
  190       continue
  200     continue
  210   continue
        do 230 j=0,5
          do 220 i=1,4
            sigt(i,0,j)=sigt(i,1,j)+sigt(i,2,j)+sigt(i,3,j)+sigt(i,4,j)
            sigt(0,i,j)=sigt(1,i,j)+sigt(2,i,j)+sigt(3,i,j)+sigt(4,i,j)
  220     continue
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
  240       continue
  250     continue
  260   continue
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYMAXI
C...Finds optimal set of coefficients for kinematical variable selection
C...and the maximum of the part of the differential cross-section used
C...in the event weighting.
 
      subroutine pjmaxi
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Parameter statement to help give large particle numbers.
      parameter (ksusy1=1000000,ksusy2=2000000,kexcit=4000000)
C...Commonblocks.
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
      common/pjint6/proc(0:500)
      character proc*28
      common/pjint7/sigt(0:6,0:6,0:5)
      save /jydat1/,/jydat2/,/jydat3/,/pjsubs/,/pjpars/,/pjint1/,
     &/pjint2/,/pjint3/,/pjint4/,/pjint5/,/pjint6/,/pjint7/
C...Local arrays, character variables and data.
      character cvar(4)*4
      dimension npts(4),mvarpt(500,4),vintpt(500,30),sigspt(500),
     &narel(7),wtrel(7),wtmat(7,7),wtreln(7),coefu(7),coefo(7),
     &iaccmx(4),sigsmx(4),sigssm(3),pmmn(2)
      data cvar/'tau ','tau''','y*  ','cth '/
      data sigssm/3*0d0/
 
C...Select subprocess to study: skip cases not applicable.
      nposi=0
      vint(143)=1d0
      vint(144)=1d0
      xsec(0,1)=0d0
      do 460 isub=1,500
        mint(51)=0
        if(iset(isub).eq.11) then
          xsec(isub,1)=1.00001d0*coef(isub,1)
          nposi=nposi+1
          goto 450
        elseif(isub.ge.91.and.isub.le.95) then
          xsec(isub,1)=sigt(0,0,isub-90)
          if(msub(isub).ne.1) goto 460
          nposi=nposi+1
          goto 450
        elseif(isub.eq.96) then
          if(mint(50).eq.0) goto 460
          if(msub(95).ne.1.and.mstp(81).le.0.and.mstp(131).le.0)
     &    goto 460
          if(mint(49).eq.0.and.mstp(131).eq.0) goto 460
        elseif(isub.eq.11.or.isub.eq.12.or.isub.eq.13.or.isub.eq.28.or.
     &    isub.eq.53.or.isub.eq.68) then
          if(msub(isub).ne.1.or.msub(95).eq.1) goto 460
        else
          if(msub(isub).ne.1) goto 460
        endif
        mint(1)=isub
        istsb=iset(isub)
        if(isub.eq.96) istsb=2
        if(mstp(122).ge.2) write(mstu(11),5000) isub
        mwtxs=0
        if(mstp(142).ge.1.and.isub.ne.96.and.msub(91)+msub(92)+msub(93)+
     &  msub(94)+msub(95).eq.0) mwtxs=1
 
C...Find resonances (explicit or implicit in cross-section).
        mint(72)=0
        kfr1=0
        if(istsb.eq.1.or.istsb.eq.3.or.istsb.eq.5) then
          kfr1=kfpr(isub,1)
        elseif(isub.eq.24.or.isub.eq.25.or.isub.eq.110.or.isub.eq.165
     &    .or.isub.eq.171.or.isub.eq.176) then
          kfr1=23
        elseif(isub.eq.23.or.isub.eq.26.or.isub.eq.166.or.isub.eq.172
     &    .or.isub.eq.177) then
          kfr1=24
        elseif(isub.ge.71.and.isub.le.77) then
          kfr1=25
          if(mstp(46).eq.5) then
            kfr1=30
            pmas(30,1)=parp(45)
            pmas(30,2)=parp(45)**3/(96d0*paru(1)*parp(47)**2)
          endif
        elseif(isub.eq.194) then
          kfr1=54
        endif
        ckmx=ckin(2)
        if(ckmx.le.0d0) ckmx=vint(1)
        kcr1=jamcomp(kfr1)
        if(kfr1.ne.0) then
          if(ckin(1).gt.pmas(kcr1,1)+20d0*pmas(kcr1,2).or.
     &    ckmx.lt.pmas(kcr1,1)-20d0*pmas(kcr1,2)) kfr1=0
        endif
        if(kfr1.ne.0) then
          taur1=pmas(kcr1,1)**2/vint(2)
          gamr1=pmas(kcr1,1)*pmas(kcr1,2)/vint(2)
          mint(72)=1
          mint(73)=kfr1
          vint(73)=taur1
          vint(74)=gamr1
        endif
        kfr2=0
        if(isub.eq.141.or.isub.eq.194) then
          kfr2=23
          if(isub.eq.194) kfr2=56
          kcr2=jamcomp(kfr2)
          taur2=pmas(kcr2,1)**2/vint(2)
          gamr2=pmas(kcr2,1)*pmas(kcr2,2)/vint(2)
          if(ckin(1).gt.pmas(kcr2,1)+20d0*pmas(kcr2,2).or.
     &    ckmx.lt.pmas(kcr2,1)-20d0*pmas(kcr2,2)) kfr2=0
          if(kfr2.ne.0.and.kfr1.ne.0) then
            mint(72)=2
            mint(74)=kfr2
            vint(75)=taur2
            vint(76)=gamr2
          elseif(kfr2.ne.0) then
            kfr1=kfr2
            taur1=taur2
            gamr1=gamr2
            mint(72)=1
            mint(73)=kfr1
            vint(73)=taur1
            vint(74)=gamr1
            kfr2=0
          endif
        endif
 
C...Find product masses and minimum pT of process.
        sqm3=0d0
        sqm4=0d0
        mint(71)=0
        vint(71)=ckin(3)
        vint(80)=1d0
        if(istsb.eq.2.or.istsb.eq.4) then
          nbw=0
          do 110 i=1,2
            pmmn(i)=0d0
            if(kfpr(isub,i).eq.0) then
            elseif(mstp(42).le.0.or.pmas(jamcomp(kfpr(isub,i)),2).lt.
     &        parp(41)) then
              if(i.eq.1) sqm3=pmas(jamcomp(kfpr(isub,i)),1)**2
              if(i.eq.2) sqm4=pmas(jamcomp(kfpr(isub,i)),1)**2
            else
              nbw=nbw+1
C...This prevents SUSY/t particles from becoming too light.
              kflw=kfpr(isub,i)
              if(kflw/ksusy1.eq.1.or.kflw/ksusy1.eq.2) then
                kcw=jamcomp(kflw)
                pmmn(i)=pmas(kcw,1)
                do 100 idc=mdcy(kcw,2),mdcy(kcw,2)+mdcy(kcw,3)-1
                  if(mdme(idc,1).gt.0.and.brat(idc).gt.1e-4) then
                    pmsum=pmas(jamcomp(kfdp(idc,1)),1)+
     &              pmas(jamcomp(kfdp(idc,2)),1)
                    if(kfdp(idc,3).ne.0) pmsum=pmsum+
     &              pmas(jamcomp(kfdp(idc,3)),1)
                    pmmn(i)=min(pmmn(i),pmsum)
                  endif
  100           continue
              elseif(kflw.eq.6) then
                pmmn(i)=pmas(24,1)+pmas(5,1)
              endif
            endif
  110     continue
          if(nbw.ge.1) then
            ckin41=ckin(41)
            ckin43=ckin(43)
            ckin(41)=max(pmmn(1),ckin(41))
            ckin(43)=max(pmmn(2),ckin(43))
            call pjofsh(3,0,kfpr(isub,1),kfpr(isub,2),0d0,pqm3,pqm4)
            ckin(41)=ckin41
            ckin(43)=ckin43
            if(mint(51).eq.1) then
              write(mstu(11),5100) isub
              msub(isub)=0
              goto 460
            endif
            sqm3=pqm3**2
            sqm4=pqm4**2
          endif
          if(min(sqm3,sqm4).lt.ckin(6)**2) mint(71)=1
          if(mint(71).eq.1) vint(71)=max(ckin(3),ckin(5))
          if(isub.eq.96.and.mstp(82).le.1) vint(71)=parp(81)
          if(isub.eq.96.and.mstp(82).ge.2) vint(71)=0.08d0*parp(82)
        endif
        vint(63)=sqm3
        vint(64)=sqm4
 
C...Prepare for additional variable choices in 2 -> 3.
        if(istsb.eq.5) then
          vint(201)=0d0
          if(kfpr(isub,2).gt.0) vint(201)=pmas(jamcomp(kfpr(isub,2)),1)
          vint(206)=vint(201)
          vint(204)=pmas(23,1)
          if(isub.eq.124) vint(204)=pmas(24,1)
          if(isub.eq.121.or.isub.eq.122.or.isub.eq.181.or.isub.eq.182
     &    .or.isub.eq.186.or.isub.eq.187) vint(204)=vint(201)
          vint(209)=vint(204)
        endif
 
C...Number of points for each variable: tau, tau', y*, cos(theta-hat).
        npts(1)=2+2*mint(72)
        if(mint(47).eq.1) then
          if(istsb.eq.1.or.istsb.eq.2) npts(1)=1
        elseif(mint(47).eq.5) then
          if(istsb.le.2.or.istsb.gt.5) npts(1)=npts(1)+1
        endif
        npts(2)=1
        if(istsb.ge.3.and.istsb.le.5) then
          if(mint(47).ge.2) npts(2)=2
          if(mint(47).eq.5) npts(2)=3
        endif
        npts(3)=1
        if(mint(47).ge.4) npts(3)=3
        if(mint(45).eq.3) npts(3)=npts(3)+1
        if(mint(46).eq.3) npts(3)=npts(3)+1
        npts(4)=1
        if(istsb.eq.2.or.istsb.eq.4) npts(4)=5
        ntry=npts(1)*npts(2)*npts(3)*npts(4)
 
C...Reset coefficients of cross-section weighting.
        do 120 j=1,20
          coef(isub,j)=0d0
  120   continue
        coef(isub,1)=1d0
        coef(isub,8)=0.5d0
        coef(isub,9)=0.5d0
        coef(isub,13)=1d0
        coef(isub,18)=1d0
        mcth=0
        mtaup=0
        metaup=0
        vint(23)=0d0
        vint(26)=0d0
        sigsam=0d0
 
C...Find limits and select tau, y*, cos(theta-hat) and tau' values,
C...in grid of phase space points.
        call pjklim(1)
        metau=mint(51)
        nacc=0
        do 150 itry=1,ntry
          mint(51)=0
          if(metau.eq.1) goto 150
          if(mod(itry-1,npts(2)*npts(3)*npts(4)).eq.0) then
            mtau=1+(itry-1)/(npts(2)*npts(3)*npts(4))
            if(mtau.gt.2+2*mint(72)) mtau=7
            rtau=0.5d0
C...Special case when both resonances have same mass,
C...as is often the case in process 194.
            if(mint(72).eq.2) then
              if(abs(pmas(kcr2,1)-pmas(kcr1,1)).lt.
     &        0.01d0*(pmas(kcr2,1)+pmas(kcr1,1))) then
                if(mtau.eq.3.or.mtau.eq.4) then
                  rtau=0.4d0
                elseif(mtau.eq.5.or.mtau.eq.6) then
                  rtau=0.6d0
                endif
              endif
            endif
            call pjkmap(1,mtau,rtau)
            if(istsb.ge.3.and.istsb.le.5) call pjklim(4)
            metaup=mint(51)
          endif
          if(metaup.eq.1) goto 150
          if(istsb.ge.3.and.istsb.le.5.and.mod(itry-1,npts(3)*npts(4))
     &    .eq.0) then
            mtaup=1+mod((itry-1)/(npts(3)*npts(4)),npts(2))
            call pjkmap(4,mtaup,0.5d0)
          endif
          if(mod(itry-1,npts(3)*npts(4)).eq.0) then
            call pjklim(2)
            meyst=mint(51)
          endif
          if(meyst.eq.1) goto 150
          if(mod(itry-1,npts(4)).eq.0) then
            myst=1+mod((itry-1)/npts(4),npts(3))
            if(myst.eq.4.and.mint(45).ne.3) myst=5
            call pjkmap(2,myst,0.5d0)
            call pjklim(3)
            mecth=mint(51)
          endif
          if(mecth.eq.1) goto 150
          if(istsb.eq.2.or.istsb.eq.4) then
            mcth=1+mod(itry-1,npts(4))
            call pjkmap(3,mcth,0.5d0)
          endif
          if(isub.eq.96) vint(25)=vint(21)*(1d0-vint(23)**2)
 
C...Store position and limits.
          mint(51)=0
          call pjklim(0)
          if(mint(51).eq.1) goto 150
          nacc=nacc+1
          mvarpt(nacc,1)=mtau
          mvarpt(nacc,2)=mtaup
          mvarpt(nacc,3)=myst
          mvarpt(nacc,4)=mcth
          do 130 j=1,30
            vintpt(nacc,j)=vint(10+j)
  130     continue
 
C...Normal case: calculate cross-section.
          if(istsb.ne.5) then
            call pjsigh(nchn,sigs)
            if(mwtxs.eq.1) then
              call pjevwt(wtxs)
              sigs=wtxs*sigs
            endif
 
C..2 -> 3: find highest value out of a number of tries.
          else
            sigs=0d0
            do 140 ikin3=1,mstp(129)
              call pjkmap(5,0,0d0)
              if(mint(51).eq.1) goto 140
              call pjsigh(nchn,sigtmp)
              if(mwtxs.eq.1) then
                call pjevwt(wtxs)
                sigtmp=wtxs*sigtmp
              endif
              if(sigtmp.gt.sigs) sigs=sigtmp
  140       continue
          endif
 
C...Store cross-section.
          sigspt(nacc)=sigs
          if(sigs.gt.sigsam) sigsam=sigs
          if(mstp(122).ge.2) write(mstu(11),5200) mtau,myst,mcth,mtaup,
     &    vint(21),vint(22),vint(23),vint(26),sigs
  150   continue
        if(nacc.eq.0) then
          write(mstu(11),5100) isub
          msub(isub)=0
          goto 460
        elseif(sigsam.eq.0d0) then
          write(mstu(11),5300) isub
          msub(isub)=0
          goto 460
        endif
        if(isub.ne.96) nposi=nposi+1
 
C...Calculate integrals in tau over maximal phase space limits.
        taumin=vint(11)
        taumax=vint(31)
        atau1=log(taumax/taumin)
        if(npts(1).ge.2) then
          atau2=(taumax-taumin)/(taumax*taumin)
        endif
        if(npts(1).ge.4) then
          atau3=log(taumax/taumin*(taumin+taur1)/(taumax+taur1))/taur1
          atau4=(atan((taumax-taur1)/gamr1)-atan((taumin-taur1)/gamr1))/
     &    gamr1
        endif
        if(npts(1).ge.6) then
          atau5=log(taumax/taumin*(taumin+taur2)/(taumax+taur2))/taur2
          atau6=(atan((taumax-taur2)/gamr2)-atan((taumin-taur2)/gamr2))/
     &    gamr2
        endif
        if(npts(1).gt.2+2*mint(72)) then
          atau7=log(max(2d-6,1d0-taumin)/max(2d-6,1d0-taumax))
        endif
 
C...Reset. Sum up cross-sections in points calculated.
        do 320 ivar=1,4
          if(npts(ivar).eq.1) goto 320
          if(isub.eq.96.and.ivar.eq.4) goto 320
          nbin=npts(ivar)
          do 170 j1=1,nbin
            narel(j1)=0
            wtrel(j1)=0d0
            coefu(j1)=0d0
            do 160 j2=1,nbin
              wtmat(j1,j2)=0d0
  160       continue
  170     continue
          do 180 iacc=1,nacc
            ibin=mvarpt(iacc,ivar)
            if(ivar.eq.1.and.ibin.eq.7) ibin=3+2*mint(72)
            if(ivar.eq.3.and.ibin.eq.5.and.mint(45).ne.3) ibin=4
            narel(ibin)=narel(ibin)+1
            wtrel(ibin)=wtrel(ibin)+sigspt(iacc)
 
C...Sum up tau cross-section pieces in points used.
            if(ivar.eq.1) then
              tau=vintpt(iacc,11)
              wtmat(ibin,1)=wtmat(ibin,1)+1d0
              wtmat(ibin,2)=wtmat(ibin,2)+(atau1/atau2)/tau
              if(nbin.ge.4) then
                wtmat(ibin,3)=wtmat(ibin,3)+(atau1/atau3)/(tau+taur1)
                wtmat(ibin,4)=wtmat(ibin,4)+(atau1/atau4)*tau/
     &          ((tau-taur1)**2+gamr1**2)
              endif
              if(nbin.ge.6) then
                wtmat(ibin,5)=wtmat(ibin,5)+(atau1/atau5)/(tau+taur2)
                wtmat(ibin,6)=wtmat(ibin,6)+(atau1/atau6)*tau/
     &          ((tau-taur2)**2+gamr2**2)
              endif
              if(nbin.gt.2+2*mint(72)) then
                wtmat(ibin,nbin)=wtmat(ibin,nbin)+(atau1/atau7)*
     &          tau/max(2d-6,1d0-tau)
              endif
 
C...Sum up tau' cross-section pieces in points used.
            elseif(ivar.eq.2) then
              tau=vintpt(iacc,11)
              taup=vintpt(iacc,16)
              taupmn=vintpt(iacc,6)
              taupmx=vintpt(iacc,26)
              ataup1=log(taupmx/taupmn)
              ataup2=((1d0-tau/taupmx)**4-(1d0-tau/taupmn)**4)/(4d0*tau)
              wtmat(ibin,1)=wtmat(ibin,1)+1d0
              wtmat(ibin,2)=wtmat(ibin,2)+(ataup1/ataup2)*
     &        (1d0-tau/taup)**3/taup
              if(nbin.ge.3) then
                ataup3=log(max(2d-6,1d0-taupmn)/max(2d-6,1d0-taupmx))
                wtmat(ibin,3)=wtmat(ibin,3)+(ataup1/ataup3)*
     &          taup/max(2d-6,1d0-taup)
              endif
 
C...Sum up y* cross-section pieces in points used.
            elseif(ivar.eq.3) then
              yst=vintpt(iacc,12)
              ystmin=vintpt(iacc,2)
              ystmax=vintpt(iacc,22)
              ayst0=ystmax-ystmin
              ayst1=0.5d0*(ystmax-ystmin)**2
              ayst2=ayst1
              ayst3=2d0*(atan(exp(ystmax))-atan(exp(ystmin)))
              wtmat(ibin,1)=wtmat(ibin,1)+(ayst0/ayst1)*(yst-ystmin)
              wtmat(ibin,2)=wtmat(ibin,2)+(ayst0/ayst2)*(ystmax-yst)
              wtmat(ibin,3)=wtmat(ibin,3)+(ayst0/ayst3)/cosh(yst)
              if(mint(45).eq.3) then
                taue=vintpt(iacc,11)
                if(istsb.ge.3.and.istsb.le.5) taue=vintpt(iacc,16)
                yst0=-0.5d0*log(taue)
                ayst4=log(max(1d-6,exp(yst0-ystmin)-1d0)/
     &          max(1d-6,exp(yst0-ystmax)-1d0))
                wtmat(ibin,4)=wtmat(ibin,4)+(ayst0/ayst4)/
     &          max(1d-6,1d0-exp(yst-yst0))
              endif
              if(mint(46).eq.3) then
                taue=vintpt(iacc,11)
                if(istsb.ge.3.and.istsb.le.5) taue=vintpt(iacc,16)
                yst0=-0.5d0*log(taue)
                ayst5=log(max(1d-6,exp(yst0+ystmax)-1d0)/
     &          max(1d-6,exp(yst0+ystmin)-1d0))
                wtmat(ibin,nbin)=wtmat(ibin,nbin)+(ayst0/ayst5)/
     &          max(1d-6,1d0-exp(-yst-yst0))
              endif
 
C...Sum up cos(theta-hat) cross-section pieces in points used.
            else
              rm34=max(1d-20,2d0*sqm3*sqm4/(vintpt(iacc,11)*vint(2))**2)
              rsqm=1d0+rm34
              cthmax=sqrt(1d0-4d0*vint(71)**2/(taumax*vint(2)))
              cthmin=-cthmax
              if(cthmax.gt.0.9999d0) rm34=max(rm34,2d0*vint(71)**2/
     &        (taumax*vint(2)))
              acth1=cthmax-cthmin
              acth2=log(max(rm34,rsqm-cthmin)/max(rm34,rsqm-cthmax))
              acth3=log(max(rm34,rsqm+cthmax)/max(rm34,rsqm+cthmin))
              acth4=1d0/max(rm34,rsqm-cthmax)-1d0/max(rm34,rsqm-cthmin)
              acth5=1d0/max(rm34,rsqm+cthmin)-1d0/max(rm34,rsqm+cthmax)
              cth=vintpt(iacc,13)
              wtmat(ibin,1)=wtmat(ibin,1)+1d0
              wtmat(ibin,2)=wtmat(ibin,2)+(acth1/acth2)/
     &        max(rm34,rsqm-cth)
              wtmat(ibin,3)=wtmat(ibin,3)+(acth1/acth3)/
     &        max(rm34,rsqm+cth)
              wtmat(ibin,4)=wtmat(ibin,4)+(acth1/acth4)/
     &        max(rm34,rsqm-cth)**2
              wtmat(ibin,5)=wtmat(ibin,5)+(acth1/acth5)/
     &        max(rm34,rsqm+cth)**2
            endif
  180     continue
 
C...Check that equation system solvable.
          if(mstp(122).ge.2) write(mstu(11),5400) cvar(ivar)
          msolv=1
          wtrels=0d0
          do 190 ibin=1,nbin
            if(mstp(122).ge.2) write(mstu(11),5500) (wtmat(ibin,ired),
     &      ired=1,nbin),wtrel(ibin)
            if(narel(ibin).eq.0) msolv=0
            wtrels=wtrels+wtrel(ibin)
  190     continue
          if(abs(wtrels).lt.1d-20) msolv=0
 
C...Solve to find relative importance of cross-section pieces.
          if(msolv.eq.1) then
            do 200 ibin=1,nbin
              wtreln(ibin)=max(0.1d0,wtrel(ibin)/wtrels)
  200       continue
            do 230 ired=1,nbin-1
              do 220 ibin=ired+1,nbin
                if(abs(wtmat(ired,ired)).lt.1d-20) then
                  msolv=0
                  goto 260
                endif
                rqt=wtmat(ibin,ired)/wtmat(ired,ired)
                wtrel(ibin)=wtrel(ibin)-rqt*wtrel(ired)
                do 210 icoe=ired,nbin
                  wtmat(ibin,icoe)=wtmat(ibin,icoe)-rqt*wtmat(ired,icoe)
  210           continue
  220         continue
  230       continue
            do 250 ired=nbin,1,-1
              do 240 icoe=ired+1,nbin
                wtrel(ired)=wtrel(ired)-wtmat(ired,icoe)*coefu(icoe)
  240         continue
              coefu(ired)=wtrel(ired)/wtmat(ired,ired)
  250       continue
          endif
 
C...Share evenly if failure.
  260     if(msolv.eq.0) then
            do 270 ibin=1,nbin
              coefu(ibin)=1d0
              wtreln(ibin)=0.1d0
              if(wtrels.gt.0d0) wtreln(ibin)=max(0.1d0,
     &        wtrel(ibin)/wtrels)
  270       continue
          endif
 
C...Normalize coefficients, with piece shared democratically.
          coefsu=0d0
          wtrels=0d0
          do 280 ibin=1,nbin
            coefu(ibin)=max(0d0,coefu(ibin))
            coefsu=coefsu+coefu(ibin)
            wtrels=wtrels+wtreln(ibin)
  280     continue
          if(coefsu.gt.0d0) then
            do 290 ibin=1,nbin
              coefo(ibin)=parp(122)/nbin+(1d0-parp(122))*0.5d0*
     &        (coefu(ibin)/coefsu+wtreln(ibin)/wtrels)
  290       continue
          else
            do 300 ibin=1,nbin
              coefo(ibin)=1d0/nbin
  300       continue
          endif
          if(ivar.eq.1) ioff=0
          if(ivar.eq.2) ioff=17
          if(ivar.eq.3) ioff=7
          if(ivar.eq.4) ioff=12
          do 310 ibin=1,nbin
            icof=ioff+ibin
            if(ivar.eq.1.and.ibin.gt.2+2*mint(72)) icof=7
            if(ivar.eq.3.and.ibin.eq.4.and.mint(45).ne.3) icof=icof+1
            coef(isub,icof)=coefo(ibin)
  310     continue
          if(mstp(122).ge.2) write(mstu(11),5600) cvar(ivar),
     &    (coefo(ibin),ibin=1,nbin)
  320   continue
 
C...Find two most promising maxima among points previously determined.
        do 330 j=1,4
          iaccmx(j)=0
          sigsmx(j)=0d0
  330   continue
        nmax=0
        do 390 iacc=1,nacc
          do 340 j=1,30
            vint(10+j)=vintpt(iacc,j)
  340     continue
          if(istsb.ne.5) then
            call pjsigh(nchn,sigs)
            if(mwtxs.eq.1) then
              call pjevwt(wtxs)
              sigs=wtxs*sigs
            endif
          else
            sigs=0d0
            do 350 ikin3=1,mstp(129)
              call pjkmap(5,0,0d0)
              if(mint(51).eq.1) goto 350
              call pjsigh(nchn,sigtmp)
              if(mwtxs.eq.1) then
                call pjevwt(wtxs)
                sigtmp=wtxs*sigtmp
              endif
              if(sigtmp.gt.sigs) sigs=sigtmp
  350       continue
          endif
          ieq=0
          do 360 imv=1,nmax
            if(abs(sigs-sigsmx(imv)).lt.1d-4*(sigs+sigsmx(imv))) ieq=imv
  360     continue
          if(ieq.eq.0) then
            do 370 imv=nmax,1,-1
              iin=imv+1
              if(sigs.le.sigsmx(imv)) goto 380
              iaccmx(imv+1)=iaccmx(imv)
              sigsmx(imv+1)=sigsmx(imv)
  370       continue
            iin=1
  380       iaccmx(iin)=iacc
            sigsmx(iin)=sigs
            if(nmax.le.1) nmax=nmax+1
          endif
  390   continue
 
C...Read out starting position for search.
        if(mstp(122).ge.2) write(mstu(11),5700)
        sigsam=sigsmx(1)
        do 440 imax=1,nmax
          iacc=iaccmx(imax)
          mtau=mvarpt(iacc,1)
          mtaup=mvarpt(iacc,2)
          myst=mvarpt(iacc,3)
          mcth=mvarpt(iacc,4)
          vtau=0.5d0
          vyst=0.5d0
          vcth=0.5d0
          vtaup=0.5d0
 
C...Starting point and step size in parameter space.
          do 430 irpt=1,2
            do 420 ivar=1,4
              if(npts(ivar).eq.1) goto 420
              if(ivar.eq.1) vvar=vtau
              if(ivar.eq.2) vvar=vtaup
              if(ivar.eq.3) vvar=vyst
              if(ivar.eq.4) vvar=vcth
              if(ivar.eq.1) mvar=mtau
              if(ivar.eq.2) mvar=mtaup
              if(ivar.eq.3) mvar=myst
              if(ivar.eq.4) mvar=mcth
              if(irpt.eq.1) vdel=0.1d0
              if(irpt.eq.2) vdel=max(0.01d0,min(0.05d0,vvar-0.02d0,
     &        0.98d0-vvar))
              if(irpt.eq.1) vmar=0.02d0
              if(irpt.eq.2) vmar=0.002d0
              imov0=1
              if(irpt.eq.1.and.ivar.eq.1) imov0=0
              do 410 imov=imov0,8
 
C...Define new point in parameter space.
                if(imov.eq.0) then
                  inew=2
                  vnew=vvar
                elseif(imov.eq.1) then
                  inew=3
                  vnew=vvar+vdel
                elseif(imov.eq.2) then
                  inew=1
                  vnew=vvar-vdel
                elseif(sigssm(3).ge.max(sigssm(1),sigssm(2)).and.
     &            vvar+2d0*vdel.lt.1d0-vmar) then
                  vvar=vvar+vdel
                  sigssm(1)=sigssm(2)
                  sigssm(2)=sigssm(3)
                  inew=3
                  vnew=vvar+vdel
                elseif(sigssm(1).ge.max(sigssm(2),sigssm(3)).and.
     &            vvar-2d0*vdel.gt.vmar) then
                  vvar=vvar-vdel
                  sigssm(3)=sigssm(2)
                  sigssm(2)=sigssm(1)
                  inew=1
                  vnew=vvar-vdel
                elseif(sigssm(3).ge.sigssm(1)) then
                  vdel=0.5d0*vdel
                  vvar=vvar+vdel
                  sigssm(1)=sigssm(2)
                  inew=2
                  vnew=vvar
                else
                  vdel=0.5d0*vdel
                  vvar=vvar-vdel
                  sigssm(3)=sigssm(2)
                  inew=2
                  vnew=vvar
                endif
 
C...Convert to relevant variables and find derived new limits.
                ilerr=0
                if(ivar.eq.1) then
                  vtau=vnew
                  call pjkmap(1,mtau,vtau)
                  if(istsb.ge.3.and.istsb.le.5) then
                    call pjklim(4)
                    if(mint(51).eq.1) ilerr=1
                  endif
                endif
                if(ivar.le.2.and.istsb.ge.3.and.istsb.le.5.and.
     &          ilerr.eq.0) then
                  if(ivar.eq.2) vtaup=vnew
                  call pjkmap(4,mtaup,vtaup)
                endif
                if(ivar.le.2.and.ilerr.eq.0) then
                  call pjklim(2)
                  if(mint(51).eq.1) ilerr=1
                endif
                if(ivar.le.3.and.ilerr.eq.0) then
                  if(ivar.eq.3) vyst=vnew
                  call pjkmap(2,myst,vyst)
                  call pjklim(3)
                  if(mint(51).eq.1) ilerr=1
                endif
                if((istsb.eq.2.or.istsb.eq.4.or.istsb.eq.6).and.
     &          ilerr.eq.0) then
                  if(ivar.eq.4) vcth=vnew
                  call pjkmap(3,mcth,vcth)
                endif
                if(isub.eq.96) vint(25)=vint(21)*(1.-vint(23)**2)
 
C...Evaluate cross-section. Save new maximum. Final maximum.
                if(ilerr.ne.0) then
                   sigs=0.
                elseif(istsb.ne.5) then
                  call pjsigh(nchn,sigs)
                  if(mwtxs.eq.1) then
                    call pjevwt(wtxs)
                    sigs=wtxs*sigs
                  endif
                else
                  sigs=0d0
                  do 400 ikin3=1,mstp(129)
                    call pjkmap(5,0,0d0)
                    if(mint(51).eq.1) goto 400
                    call pjsigh(nchn,sigtmp)
                    if(mwtxs.eq.1) then
                        call pjevwt(wtxs)
                        sigtmp=wtxs*sigtmp
                    endif
                    if(sigtmp.gt.sigs) sigs=sigtmp
  400             continue
                endif
                sigssm(inew)=sigs
                if(sigs.gt.sigsam) sigsam=sigs
                if(mstp(122).ge.2) write(mstu(11),5800) imax,ivar,mvar,
     &          imov,vnew,vint(21),vint(22),vint(23),vint(26),sigs
  410         continue
  420       continue
  430     continue
  440   continue
        if(mstp(121).eq.1) sigsam=parp(121)*sigsam
        xsec(isub,1)=1.05d0*sigsam
  450   continue
        if(mstp(173).eq.1.and.isub.ne.96) xsec(isub,1)=
     &  parp(174)*xsec(isub,1)
        if(isub.ne.96) xsec(0,1)=xsec(0,1)+xsec(isub,1)
  460 continue
      mint(51)=0
 
C...Print summary table.
      if(nposi.eq.0) then
        write(mstu(11),5900)
        stop
      endif
      if(mstp(122).ge.1) then
        write(mstu(11),6000)
        write(mstu(11),6100)
        do 470 isub=1,500
          if(msub(isub).ne.1.and.isub.ne.96) goto 470
          if(isub.eq.96.and.mint(50).eq.0) goto 470
          if(isub.eq.96.and.msub(95).ne.1.and.mstp(81).le.0) goto 470
          if(isub.eq.96.and.mint(49).eq.0.and.mstp(131).eq.0) goto 470
          if(msub(95).eq.1.and.(isub.eq.11.or.isub.eq.12.or.isub.eq.13
     &    .or.isub.eq.28.or.isub.eq.53.or.isub.eq.68)) goto 470
          write(mstu(11),6200) isub,proc(isub),xsec(isub,1)
  470   continue
        write(mstu(11),6300)
      endif
 
C...Format statements for maximization results.
 5000 format(/1x,'Coefficient optimization and maximum search for ',
     &'subprocess no',i4/1x,'Coefficient modes     tau',10x,'y*',9x,
     &'cth',9x,'tau''',7x,'sigma')
 5100 format(1x,'Warning: requested subprocess ',i3,' has no allowed ',
     &'phase space.'/1x,'Process switched off!')
 5200 format(1x,4i4,f12.8,f12.6,f12.7,f12.8,1p,d12.4)
 5300 format(1x,'Warning: requested subprocess ',i3,' has vanishing ',
     &'cross-section.'/1x,'Process switched off!')
 5400 format(1x,'Coefficients of equation system to be solved for ',a4)
 5500 format(1x,1p,8d11.3)
 5600 format(1x,'Result for ',a4,':',7f9.4)
 5700 format(1x,'Maximum search for given coefficients'/2x,'MAX VAR ',
     &'MOD MOV   VNEW',7x,'tau',7x,'y*',8x,'cth',7x,'tau''',7x,'sigma')
 5800 format(1x,4i4,f8.4,f11.7,f9.3,f11.6,f11.7,1p,d12.4)
 5900 format(1x,'Error: no requested process has non-vanishing ',
     &'cross-section.'/1x,'Execution stopped!')
 6000 format(/1x,8('*'),1x,'PYMAXI: summary of differential ',
     &'cross-section maximum search',1x,8('*'))
 6100 format(/11x,58('=')/11x,'I',38x,'I',17x,'I'/11x,'I  ISUB  ',
     &'Subprocess name',15x,'I  Maximum value  I'/11x,'I',38x,'I',
     &17x,'I'/11x,58('=')/11x,'I',38x,'I',17x,'I')
 6200 format(11x,'I',2x,i3,3x,a28,2x,'I',2x,1p,d12.4,3x,'I')
 6300 format(11x,'I',38x,'I',17x,'I'/11x,58('='))
 
      return
      end
 
C*********************************************************************
 
C...PYPILE
C...Initializes multiplicity distribution and selects mutliplicity
C...of pileup events, i.e. several events occuring at the same
C...beam crossing.
 
      subroutine pjpile(mpile)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint7/sigt(0:6,0:6,0:5)
      save /jydat1/,/pjpars/,/pjint1/,/pjint7/
C...Local arrays and saved variables.
      dimension wti(0:200)
      save imin,imax,wti,wts
 
C...Sum of allowed cross-sections for pileup events.
      if(mpile.eq.1) then
        vint(131)=sigt(0,0,5)
        if(mstp(132).ge.2) vint(131)=vint(131)+sigt(0,0,4)
        if(mstp(132).ge.3) vint(131)=vint(131)+sigt(0,0,2)+sigt(0,0,3)
        if(mstp(132).ge.4) vint(131)=vint(131)+sigt(0,0,1)
        if(mstp(133).le.0) return
 
C...Initialize multiplicity distribution at maximum.
        xnave=vint(131)*parp(131)
        if(xnave.gt.120d0) write(mstu(11),5000) xnave
        inave=max(1,min(200,nint(xnave)))
        wti(inave)=1d0
        wts=wti(inave)
        wtn=wti(inave)*inave
 
C...Find shape of multiplicity distribution below maximum.
        imin=inave
        do 100 i=inave-1,1,-1
          if(mstp(133).eq.1) wti(i)=wti(i+1)*(i+1)/xnave
          if(mstp(133).ge.2) wti(i)=wti(i+1)*i/xnave
          if(wti(i).lt.1d-6) goto 110
          wts=wts+wti(i)
          wtn=wtn+wti(i)*i
          imin=i
  100   continue
 
C...Find shape of multiplicity distribution above maximum.
  110   imax=inave
        do 120 i=inave+1,200
          if(mstp(133).eq.1) wti(i)=wti(i-1)*xnave/i
          if(mstp(133).ge.2) wti(i)=wti(i-1)*xnave/(i-1)
          if(wti(i).lt.1d-6) goto 130
          wts=wts+wti(i)
          wtn=wtn+wti(i)*i
          imax=i
  120   continue
  130   vint(132)=xnave
        vint(133)=wtn/wts
        if(mstp(133).eq.1.and.imin.eq.1) vint(134)=
     &  wts/(wts+wti(1)/xnave)
        if(mstp(133).eq.1.and.imin.gt.1) vint(134)=1d0
        if(mstp(133).ge.2) vint(134)=xnave
 
C...Pick multiplicity of pileup events.
      else
        if(mstp(133).le.0) then
          mint(81)=max(1,mstp(134))
        else
          wtr=wts*pjr(0)
          do 140 i=imin,imax
            mint(81)=i
            wtr=wtr-wti(i)
            if(wtr.le.0d0) goto 150
  140     continue
  150     continue
        endif
      endif
 
C...Format statement for error message.
 5000 format(1x,'Warning: requested average number of events per bunch',
     &'crossing too large, ',1p,d12.4)
 
      return
      end
 
C*********************************************************************
 
C...PYSAVE
C...Saves and restores parameter and cross section values for the
C...3 gamma-p and 6 gamma-gamma alnternatives. Also makes random
C...choice between alternatives.
 
      subroutine pjsave(isave,iga)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      common/pjint5/ngenpd,ngen(0:500,3),xsec(0:500,3)
      save /pjsubs/,/pjpars/,/pjint1/,/pjint2/,/pjint5/
C...Local arrays and saved variables.
      dimension ncp(10),nsubcp(10,20),msubcp(10,20),coefcp(10,20,20),
     &ngencp(10,0:20,3),xseccp(10,0:20,3),intcp(10,20),recp(10,20)
      save ncp,nsubcp,msubcp,coefcp,ngencp,xseccp,intcp,recp
 
C...Save list of subprocesses and cross-section information.
      if(isave.eq.1) then
        icp=0
        do 120 i=1,500
          if(msub(i).eq.0.and.i.ne.96.and.i.ne.97) goto 120
          icp=icp+1
          nsubcp(iga,icp)=i
          msubcp(iga,icp)=msub(i)
          do 100 j=1,20
            coefcp(iga,icp,j)=coef(i,j)
  100     continue
          do 110 j=1,3
            ngencp(iga,icp,j)=ngen(i,j)
            xseccp(iga,icp,j)=xsec(i,j)
  110     continue
  120   continue
        ncp(iga)=icp
        do 130 j=1,3
          ngencp(iga,0,j)=ngen(0,j)
          xseccp(iga,0,j)=xsec(0,j)
  130   continue
C...Save various common process variables.
        do 140 j=1,10
          intcp(iga,j)=mint(40+j)
  140   continue
        intcp(iga,11)=mint(101)
        intcp(iga,12)=mint(102)
        intcp(iga,13)=mint(107)
        intcp(iga,14)=mint(108)
        intcp(iga,15)=mint(123)
        recp(iga,1)=ckin(3)
 
C...Save cross-section information only.
      elseif(isave.eq.2) then
        do 160 icp=1,ncp(iga)
          i=nsubcp(iga,icp)
          do 150 j=1,3
            ngencp(iga,icp,j)=ngen(i,j)
            xseccp(iga,icp,j)=xsec(i,j)
  150     continue
  160   continue
        do 170 j=1,3
          ngencp(iga,0,j)=ngen(0,j)
          xseccp(iga,0,j)=xsec(0,j)
  170   continue
 
C...Choose between allowed alternatives.
      elseif(isave.eq.3.or.isave.eq.4) then
        if(isave.eq.4) then
          xsumcp=0d0
          do 180 ig=1,mint(121)
            xsumcp=xsumcp+xseccp(ig,0,1)
  180     continue
          xsumcp=xsumcp*pjr(0)
          do 190 ig=1,mint(121)
            iga=ig
            xsumcp=xsumcp-xseccp(ig,0,1)
            if(xsumcp.le.0d0) goto 200
  190     continue
  200     continue
        endif
 
C...Restore cross-section information.
        do 210 i=1,500
          msub(i)=0
  210   continue
        do 240 icp=1,ncp(iga)
          i=nsubcp(iga,icp)
          msub(i)=msubcp(iga,icp)
          do 220 j=1,20
            coef(i,j)=coefcp(iga,icp,j)
  220     continue
          do 230 j=1,3
            ngen(i,j)=ngencp(iga,icp,j)
            xsec(i,j)=xseccp(iga,icp,j)
  230     continue
  240   continue
        do 250 j=1,3
          ngen(0,j)=ngencp(iga,0,j)
          xsec(0,j)=xseccp(iga,0,j)
  250   continue
 
C...Restore various common process variables.
        do 260 j=1,10
          mint(40+j)=intcp(iga,j)
  260   continue
        mint(101)=intcp(iga,11)
        mint(102)=intcp(iga,12)
        mint(107)=intcp(iga,13)
        mint(108)=intcp(iga,14)
        mint(123)=intcp(iga,15)
        ckin(3)=recp(iga,1)
        ckin(1)=2d0*ckin(3)
 
C...Sum up cross-section info (for PYSTAT).
      elseif(isave.eq.5) then
        do 270 i=1,500
          msub(i)=0
          ngen(i,1)=0
          ngen(i,3)=0
          xsec(i,3)=0d0
  270   continue
        ngen(0,1)=0
        ngen(0,2)=0
        ngen(0,3)=0
        xsec(0,3)=0
        do 290 ig=1,mint(121)
          do 280 icp=1,ncp(ig)
            i=nsubcp(ig,icp)
            if(msubcp(ig,icp).eq.1) msub(i)=1
            ngen(i,1)=ngen(i,1)+ngencp(ig,icp,1)
            ngen(i,3)=ngen(i,3)+ngencp(ig,icp,3)
            xsec(i,3)=xsec(i,3)+xseccp(ig,icp,3)
  280     continue
          ngen(0,1)=ngen(0,1)+ngencp(ig,0,1)
          ngen(0,2)=ngen(0,2)+ngencp(ig,0,2)
          ngen(0,3)=ngen(0,3)+ngencp(ig,0,3)
          xsec(0,3)=xsec(0,3)+xseccp(ig,0,3)
  290   continue
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYRAND
C...Generates quantities characterizing the high-pT scattering at the
C...parton level according to the matrix elements. Chooses incoming,
C...reacting partons, their momentum fractions and one of the possible
C...subprocesses.
 
      subroutine pjrand
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Parameter statement to help give large particle numbers.
      parameter (ksusy1=1000000,ksusy2=2000000,kexcit=4000000)
C...Commonblocks.
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
      common/pjuppr/nup,kup(20,7),nfup,ifup(10,2),pup(20,5),q2up(0:10)
      common/pjmssm/imss(0:99),rmss(0:99)
      save /jydat1/,/jydat2/,/jydat3/,/pjsubs/,/pjpars/,/pjint1/,
     &/pjint2/,/pjint3/,/pjint4/,/pjint5/,/pjint7/,/pjuppr/,/pjmssm/
C...Local arrays.
      dimension xpq(-25:25),pmm(2),pdif(4),bhad(4),pmmn(2)
 
C...Parameters and data used in elastic/diffractive treatment.
      data eps/0.0808d0/, alp/0.25d0/, cres/2d0/, pmrc/1.062d0/,
     &smp/0.880d0/, bhad/2.3d0,1.4d0,1.4d0,0.23d0/
 
C...Initial values, specifically for (first) semihard interaction.
      mint(10)=0
      mint(17)=0
      mint(18)=0
      vint(143)=1d0
      vint(144)=1d0
      mfail=0
      if(mstp(171).eq.1.and.mstp(172).eq.2) mfail=1
      isub=0
      loop=0
  100 loop=loop+1
      mint(51)=0
 
C...Choice of process type - first event of pileup.
      if(mint(82).eq.1.and.(isub.le.90.or.isub.gt.96)) then
 
C...For gamma-p or gamma-gamma first pick between alternatives.
        if(mint(121).gt.1) call pjsave(4,iga)
        mint(122)=iga
 
C...For gamma + gamma with different nature, flip at random.
        if(mint(11).eq.22.and.mint(12).eq.22.and.mint(123).ge.4.and.
     &  pjr(0).gt.0.5d0) then
          mintsv=mint(41)
          mint(41)=mint(42)
          mint(42)=mintsv
          mintsv=mint(45)
          mint(45)=mint(46)
          mint(46)=mintsv
          mintsv=mint(107)
          mint(107)=mint(108)
          mint(108)=mintsv
          if(mint(47).eq.2.or.mint(47).eq.3) mint(47)=5-mint(47)
        endif
 
C...Pick process type.
        rsub=xsec(0,1)*pjr(0)
        do 110 i=1,500
          if(msub(i).ne.1) goto 110
          isub=i
          rsub=rsub-xsec(i,1)
          if(rsub.le.0d0) goto 120
  110   continue
  120   if(isub.eq.95) isub=96
        if(isub.eq.96) call pjmult(2)
 
C...Choice of inclusive process type - pileup events.
      elseif(mint(82).ge.2.and.isub.eq.0) then
        rsub=vint(131)*pjr(0)
        isub=96
        if(rsub.gt.sigt(0,0,5)) isub=94
        if(rsub.gt.sigt(0,0,5)+sigt(0,0,4)) isub=93
        if(rsub.gt.sigt(0,0,5)+sigt(0,0,4)+sigt(0,0,3)) isub=92
        if(rsub.gt.sigt(0,0,5)+sigt(0,0,4)+sigt(0,0,3)+sigt(0,0,2))
     &  isub=91
        if(isub.eq.96) call pjmult(2)
      endif
      if(mint(82).eq.1) ngen(0,1)=ngen(0,1)+1
      if(mint(82).eq.1) ngen(isub,1)=ngen(isub,1)+1
      if(isub.eq.96.and.loop.eq.1.and.mint(82).eq.1)
     &ngen(97,1)=ngen(97,1)+1
      mint(1)=isub
      istsb=iset(isub)
 
C...Random choice of flavour for some SUSY processes.
      if(isub.ge.201.and.isub.le.280) then
C...~e_L ~nu_e or ~mu_L ~nu_mu.
        if(isub.eq.210) then
          kfpr(isub,1)=ksusy1+11+2*int(0.5d0+pjr(0))
          kfpr(isub,2)=kfpr(isub,1)+1
C...~nu_e ~nu_e(bar) or ~nu_mu ~nu_mu(bar).
        elseif(isub.eq.213) then
          kfpr(isub,1)=ksusy1+12+2*int(0.5d0+pjr(0))
          kfpr(isub,2)=kfpr(isub,1)
C...~q ~chi/~g; ~q = ~d, ~u, ~s, ~c or ~b.
        elseif(isub.ge.246.and.isub.le.259) then
          if(mod(isub,2).eq.0) then
            kfpr(isub,1)=ksusy1+1+int(5d0*pjr(0))
          else
            kfpr(isub,1)=ksusy2+1+int(5d0*pjr(0))
          endif
C...~q1 ~q2; ~q = ~d, ~u, ~s, ~c or ~b.
        elseif(isub.ge.271.and.isub.le.276) then
          if(isub.eq.271.or.isub.eq.274) then
            ksu1=ksusy1
            ksu2=ksusy1
          elseif(isub.eq.272.or.isub.eq.275) then
            ksu1=ksusy2
            ksu2=ksusy2
          elseif(pjr(0).lt.0.5d0) then
            ksu1=ksusy1
            ksu2=ksusy2
          else
            ksu1=ksusy2
            ksu2=ksusy1
          endif
          kfpr(isub,1)=ksu1+1+int(5d0*pjr(0))
          kfpr(isub,2)=ksu2+1+int(5d0*pjr(0))
C...~q ~q(bar);  ~q = ~d, ~u, ~s, ~c or ~b.
        elseif(isub.eq.277.or.isub.eq.279) then
          kfpr(isub,1)=ksusy1+1+int(5d0*pjr(0))
          kfpr(isub,2)=kfpr(isub,1)
        elseif(isub.eq.278.or.isub.eq.280) then
          kfpr(isub,1)=ksusy2+1+int(5d0*pjr(0))
          kfpr(isub,2)=kfpr(isub,1)
        endif
      endif
 
C...Find resonances (explicit or implicit in cross-section).
      mint(72)=0
      kfr1=0
      if(istsb.eq.1.or.istsb.eq.3.or.istsb.eq.5) then
        kfr1=kfpr(isub,1)
      elseif(isub.eq.24.or.isub.eq.25.or.isub.eq.110.or.isub.eq.165.or.
     &  isub.eq.171.or.isub.eq.176) then
        kfr1=23
      elseif(isub.eq.23.or.isub.eq.26.or.isub.eq.166.or.isub.eq.172.or.
     &  isub.eq.177) then
        kfr1=24
      elseif(isub.ge.71.and.isub.le.77) then
        kfr1=25
        if(mstp(46).eq.5) then
          kfr1=30
          pmas(30,1)=parp(45)
          pmas(30,2)=parp(45)**3/(96d0*paru(1)*parp(47)**2)
        endif
      elseif(isub.eq.194) then
        kfr1=54
      endif
      ckmx=ckin(2)
      if(ckmx.le.0d0) ckmx=vint(1)
      kcr1=jamcomp(kfr1)
      if(kfr1.ne.0) then
        if(ckin(1).gt.pmas(kcr1,1)+20d0*pmas(kcr1,2).or.
     &  ckmx.lt.pmas(kcr1,1)-20d0*pmas(kcr1,2)) kfr1=0
      endif
      if(kfr1.ne.0) then
        taur1=pmas(kcr1,1)**2/vint(2)
        gamr1=pmas(kcr1,1)*pmas(kcr1,2)/vint(2)
        mint(72)=1
        mint(73)=kfr1
        vint(73)=taur1
        vint(74)=gamr1
      endif
      if(isub.eq.141.or.isub.eq.194) then
        kfr2=23
        if(isub.eq.194) kfr2=56
        kcr2=jamcomp(kfr2)
        taur2=pmas(kcr2,1)**2/vint(2)
        gamr2=pmas(kcr2,1)*pmas(kcr2,2)/vint(2)
        if(ckin(1).gt.pmas(kcr2,1)+20d0*pmas(kcr2,2).or.
     &  ckmx.lt.pmas(kcr2,1)-20d0*pmas(kcr2,2)) kfr2=0
        if(kfr2.ne.0.and.kfr1.ne.0) then
          mint(72)=2
          mint(74)=kfr2
          vint(75)=taur2
          vint(76)=gamr2
        elseif(kfr2.ne.0) then
          kfr1=kfr2
          taur1=taur2
          gamr1=gamr2
          mint(72)=1
          mint(73)=kfr1
          vint(73)=taur1
          vint(74)=gamr1
        endif
      endif
 
C...Find product masses and minimum pT of process,
C...optionally with broadening according to a truncated Breit-Wigner.
      vint(63)=0d0
      vint(64)=0d0
      mint(71)=0
      vint(71)=ckin(3)
      if(mint(82).ge.2) vint(71)=0d0
      vint(80)=1d0
      if(istsb.eq.2.or.istsb.eq.4) then
        nbw=0
        do 140 i=1,2
          pmmn(i)=0d0
          if(kfpr(isub,i).eq.0) then
          elseif(mstp(42).le.0.or.pmas(jamcomp(kfpr(isub,i)),2).lt.
     &      parp(41)) then
            vint(62+i)=pmas(jamcomp(kfpr(isub,i)),1)**2
          else
            nbw=nbw+1
C...This prevents SUSY/t particles from becoming too light.
            kflw=kfpr(isub,i)
            if(kflw/ksusy1.eq.1.or.kflw/ksusy1.eq.2) then
              kcw=jamcomp(kflw)
              pmmn(i)=pmas(kcw,1)
              do 130 idc=mdcy(kcw,2),mdcy(kcw,2)+mdcy(kcw,3)-1
                if(mdme(idc,1).gt.0.and.brat(idc).gt.1e-4) then
                  pmsum=pmas(jamcomp(kfdp(idc,1)),1)+
     &            pmas(jamcomp(kfdp(idc,2)),1)
                  if(kfdp(idc,3).ne.0) pmsum=pmsum+
     &            pmas(jamcomp(kfdp(idc,3)),1)
                  pmmn(i)=min(pmmn(i),pmsum)
                endif
  130         continue
            elseif(kflw.eq.6) then
              pmmn(i)=pmas(24,1)+pmas(5,1)
            endif
          endif
  140   continue
        if(nbw.ge.1) then
          ckin41=ckin(41)
          ckin43=ckin(43)
          ckin(41)=max(pmmn(1),ckin(41))
          ckin(43)=max(pmmn(2),ckin(43))
          call pjofsh(4,0,kfpr(isub,1),kfpr(isub,2),0d0,pqm3,pqm4)
          ckin(41)=ckin41
          ckin(43)=ckin43
          if(mint(51).eq.1) then
            if(mint(121).gt.1) call pjsave(2,iga)
            if(mfail.eq.1) then
              msti(61)=1
              return
            endif
            goto 100
          endif
          vint(63)=pqm3**2
          vint(64)=pqm4**2
        endif
        if(min(vint(63),vint(64)).lt.ckin(6)**2) mint(71)=1
        if(mint(71).eq.1) vint(71)=max(ckin(3),ckin(5))
      endif
 
C...Prepare for additional variable choices in 2 -> 3.
      if(istsb.eq.5) then
        vint(201)=0d0
        if(kfpr(isub,2).gt.0) vint(201)=pmas(jamcomp(kfpr(isub,2)),1)
        vint(206)=vint(201)
        vint(204)=pmas(23,1)
        if(isub.eq.124) vint(204)=pmas(24,1)
        if(isub.eq.121.or.isub.eq.122.or.isub.eq.181.or.isub.eq.182.or.
     &  isub.eq.186.or.isub.eq.187) vint(204)=vint(201)
        vint(209)=vint(204)
      endif
 
C...Select incoming VDM particle (rho/omega/phi/J/psi).
      if(istsb.ne.0.and.(mint(101).ge.2.or.mint(102).ge.2).and.
     &(mint(123).eq.2.or.mint(123).eq.5.or.mint(123).eq.7)) then
        vrn=pjr(0)*sigt(0,0,5)
        if(mint(101).le.1) then
          i1mn=0
          i1mx=0
        else
          i1mn=1
          i1mx=mint(101)
        endif
        if(mint(102).le.1) then
          i2mn=0
          i2mx=0
        else
          i2mn=1
          i2mx=mint(102)
        endif
        do 160 i1=i1mn,i1mx
          kfv1=110*i1+3
          do 150 i2=i2mn,i2mx
            kfv2=110*i2+3
            vrn=vrn-sigt(i1,i2,5)
            if(vrn.le.0d0) goto 170
  150     continue
  160   continue
  170   if(mint(101).ge.2) mint(103)=kfv1
        if(mint(102).ge.2) mint(104)=kfv2
      endif
 
      if(istsb.eq.0) then
C...Elastic scattering or single or double diffractive scattering.
 
C...Select incoming particle (rho/omega/phi/J/psi for VDM) and mass.
        mint(103)=mint(11)
        mint(104)=mint(12)
        pmm(1)=vint(3)
        pmm(2)=vint(4)
        if(mint(101).ge.2.or.mint(102).ge.2) then
          jj=isub-90
          vrn=pjr(0)*sigt(0,0,jj)
          if(mint(101).le.1) then
            i1mn=0
            i1mx=0
          else
            i1mn=1
            i1mx=mint(101)
          endif
          if(mint(102).le.1) then
            i2mn=0
            i2mx=0
          else
            i2mn=1
            i2mx=mint(102)
          endif
          do 190 i1=i1mn,i1mx
            kfv1=110*i1+3
            do 180 i2=i2mn,i2mx
              kfv2=110*i2+3
              vrn=vrn-sigt(i1,i2,jj)
              if(vrn.le.0d0) goto 200
  180       continue
  190     continue
  200     if(mint(101).ge.2) then
            mint(103)=kfv1
            pmm(1)=pjmass(kfv1)
          endif
          if(mint(102).ge.2) then
            mint(104)=kfv2
            pmm(2)=pjmass(kfv2)
          endif
        endif
 
C...Side/sides of diffractive system.
        mint(17)=0
        mint(18)=0
        if(isub.eq.92.or.isub.eq.94) mint(17)=1
        if(isub.eq.93.or.isub.eq.94) mint(18)=1
 
C...Find masses of particles and minimal masses of diffractive states.
        do 210 jt=1,2
          pdif(jt)=pmm(jt)
          vint(66+jt)=pdif(jt)
          if(mint(16+jt).eq.1) pdif(jt)=pdif(jt)+parp(102)
  210   continue
        sh=vint(2)
        sqm1=pmm(1)**2
        sqm2=pmm(2)**2
        sqm3=pdif(1)**2
        sqm4=pdif(2)**2
        smres1=(pmm(1)+pmrc)**2
        smres2=(pmm(2)+pmrc)**2
 
C...Find elastic slope and lower limit diffractive slope.
        iha=max(2,iabs(mint(103))/110)
        if(iha.ge.5) iha=1
        ihb=max(2,iabs(mint(104))/110)
        if(ihb.ge.5) ihb=1
        if(isub.eq.91) then
          bmn=2d0*bhad(iha)+2d0*bhad(ihb)+4d0*sh**eps-4.2d0
        elseif(isub.eq.92) then
          bmn=max(2d0,2d0*bhad(ihb))
        elseif(isub.eq.93) then
          bmn=max(2d0,2d0*bhad(iha))
        elseif(isub.eq.94) then
          bmn=2d0*alp*4d0
        endif
 
C...Determine maximum possible t range and coefficient of generation.
        sqla12=(sh-sqm1-sqm2)**2-4d0*sqm1*sqm2
        sqla34=(sh-sqm3-sqm4)**2-4d0*sqm3*sqm4
        tha=sh-(sqm1+sqm2+sqm3+sqm4)+(sqm1-sqm2)*(sqm3-sqm4)/sh
        thb=sqrt(max(0d0,sqla12))*sqrt(max(0d0,sqla34))/sh
        thc=(sqm3-sqm1)*(sqm4-sqm2)+(sqm1+sqm4-sqm2-sqm3)*
     &  (sqm1*sqm4-sqm2*sqm3)/sh
        thl=-0.5d0*(tha+thb)
        thu=thc/thl
        thrnd=exp(max(-50d0,bmn*(thl-thu)))-1d0
 
C...Select diffractive mass/masses according to dm^2/m^2.
  220   do 230 jt=1,2
          if(mint(16+jt).eq.0) then
            pdif(2+jt)=pdif(jt)
          else
            pmmin=pdif(jt)
            pmmax=max(vint(2+jt),vint(1)-pdif(3-jt))
            pdif(2+jt)=pmmin*(pmmax/pmmin)**pjr(0)
          endif
  230   continue
        sqm3=pdif(3)**2
        sqm4=pdif(4)**2
 
C..Additional mass factors, including resonance enhancement.
        if(pdif(3)+pdif(4).ge.vint(1)) goto 220
        if(isub.eq.92) then
          fsd=(1d0-sqm3/sh)*(1d0+cres*smres1/(smres1+sqm3))
          if(fsd.lt.pjr(0)*(1d0+cres)) goto 220
        elseif(isub.eq.93) then
          fsd=(1d0-sqm4/sh)*(1d0+cres*smres2/(smres2+sqm4))
          if(fsd.lt.pjr(0)*(1d0+cres)) goto 220
        elseif(isub.eq.94) then
          fdd=(1d0-(pdif(3)+pdif(4))**2/sh)*(sh*smp/
     &    (sh*smp+sqm3*sqm4))*(1d0+cres*smres1/(smres1+sqm3))*
     &    (1d0+cres*smres2/(smres2+sqm4))
          if(fdd.lt.pjr(0)*(1d0+cres)**2) goto 220
        endif
 
C...Select t according to exp(Bmn*t) and correct to right slope.
        th=thu+log(1d0+thrnd*pjr(0))/bmn
        if(isub.ge.92) then
          if(isub.eq.92) then
            badd=2d0*alp*log(sh/sqm3)
            if(bhad(ihb).lt.1d0) badd=max(0d0,badd+2d0*bhad(ihb)-2d0)
          elseif(isub.eq.93) then
            badd=2d0*alp*log(sh/sqm4)
            if(bhad(iha).lt.1d0) badd=max(0d0,badd+2d0*bhad(iha)-2d0)
          elseif(isub.eq.94) then
            badd=2d0*alp*(log(exp(4d0)+sh/(alp*sqm3*sqm4))-4d0)
          endif
          if(exp(max(-50d0,badd*(th-thu))).lt.pjr(0)) goto 220
        endif
 
C...Check whether m^2 and t choices are consistent.
        sqla34=(sh-sqm3-sqm4)**2-4d0*sqm3*sqm4
        tha=sh-(sqm1+sqm2+sqm3+sqm4)+(sqm1-sqm2)*(sqm3-sqm4)/sh
        thb=sqrt(max(0d0,sqla12))*sqrt(max(0d0,sqla34))/sh
        if(thb.le.1d-8) goto 220
        thc=(sqm3-sqm1)*(sqm4-sqm2)+(sqm1+sqm4-sqm2-sqm3)*
     &  (sqm1*sqm4-sqm2*sqm3)/sh
        thlm=-0.5d0*(tha+thb)
        thum=thc/thlm
        if(th.lt.thlm.or.th.gt.thum) goto 220
 
C...Information to output.
        vint(21)=1d0
        vint(22)=0d0
        vint(23)=min(1d0,max(-1d0,(tha+2d0*th)/thb))
        vint(45)=th
        vint(59)=2d0*sqrt(max(0d0,-(thc+tha*th+th**2)))/thb
        vint(63)=pdif(3)**2
        vint(64)=pdif(4)**2
 
C...Note: in the following, by In is meant the integral over the
C...quantity multiplying coefficient cn.
C...Choose tau according to h1(tau)/tau, where
C...h1(tau) = c1 + I1/I2*c2*1/tau + I1/I3*c3*1/(tau+tau_R) +
C...I1/I4*c4*tau/((s*tau-m^2)^2+(m*Gamma)^2) +
C...I1/I5*c5*1/(tau+tau_R') +
C...I1/I6*c6*tau/((s*tau-m'^2)^2+(m'*Gamma')^2) +
C...I1/I7*c7*tau/(1.-tau), and
C...c1 + c2 + c3 + c4 + c5 + c6 + c7 = 1.
      elseif(istsb.ge.1.and.istsb.le.5) then
        call pjklim(1)
        if(mint(51).ne.0) then
          if(mint(121).gt.1) call pjsave(2,iga)
          if(mfail.eq.1) then
            msti(61)=1
            return
          endif
          goto 100
        endif
        rtau=pjr(0)
        mtau=1
        if(rtau.gt.coef(isub,1)) mtau=2
        if(rtau.gt.coef(isub,1)+coef(isub,2)) mtau=3
        if(rtau.gt.coef(isub,1)+coef(isub,2)+coef(isub,3)) mtau=4
        if(rtau.gt.coef(isub,1)+coef(isub,2)+coef(isub,3)+coef(isub,4))
     &  mtau=5
        if(rtau.gt.coef(isub,1)+coef(isub,2)+coef(isub,3)+coef(isub,4)+
     &  coef(isub,5)) mtau=6
        if(rtau.gt.coef(isub,1)+coef(isub,2)+coef(isub,3)+coef(isub,4)+
     &  coef(isub,5)+coef(isub,6)) mtau=7
        call pjkmap(1,mtau,pjr(0))
 
C...2 -> 3, 4 processes:
C...Choose tau' according to h4(tau,tau')/tau', where
C...h4(tau,tau') = c1 + I1/I2*c2*(1 - tau/tau')^3/tau' +
C...I1/I3*c3*1/(1 - tau'), and c1 + c2 + c3 = 1.
        if(istsb.ge.3.and.istsb.le.5) then
          call pjklim(4)
          if(mint(51).ne.0) then
            if(mint(121).gt.1) call pjsave(2,iga)
            if(mfail.eq.1) then
              msti(61)=1
              return
            endif
            goto 100
          endif
          rtaup=pjr(0)
          mtaup=1
          if(rtaup.gt.coef(isub,18)) mtaup=2
          if(rtaup.gt.coef(isub,18)+coef(isub,19)) mtaup=3
          call pjkmap(4,mtaup,pjr(0))
        endif
 
C...Choose y* according to h2(y*), where
C...h2(y*) = I0/I1*c1*(y*-y*min) + I0/I2*c2*(y*max-y*) +
C...I0/I3*c3*1/cosh(y*) + I0/I4*c4*1/(1-exp(y*-y*max)) +
C...I0/I5*c5*1/(1-exp(-y*-y*min)), I0 = y*max-y*min,
C...and c1 + c2 + c3 + c4 + c5 = 1.
        call pjklim(2)
        if(mint(51).ne.0) then
          if(mint(121).gt.1) call pjsave(2,iga)
          if(mfail.eq.1) then
            msti(61)=1
            return
          endif
          goto 100
        endif
        ryst=pjr(0)
        myst=1
        if(ryst.gt.coef(isub,8)) myst=2
        if(ryst.gt.coef(isub,8)+coef(isub,9)) myst=3
        if(ryst.gt.coef(isub,8)+coef(isub,9)+coef(isub,10)) myst=4
        if(ryst.gt.coef(isub,8)+coef(isub,9)+coef(isub,10)+
     &  coef(isub,11)) myst=5
        call pjkmap(2,myst,pjr(0))
 
C...2 -> 2 processes:
C...Choose cos(theta-hat) (cth) according to h3(cth), where
C...h3(cth) = c0 + I0/I1*c1*1/(A - cth) + I0/I2*c2*1/(A + cth) +
C...I0/I3*c3*1/(A - cth)^2 + I0/I4*c4*1/(A + cth)^2,
C...A = 1 + 2*(m3*m4/sh)^2 (= 1 for massless products),
C...and c0 + c1 + c2 + c3 + c4 = 1.
        call pjklim(3)
        if(mint(51).ne.0) then
          if(mint(121).gt.1) call pjsave(2,iga)
          if(mfail.eq.1) then
            msti(61)=1
            return
          endif
          goto 100
        endif
        if(istsb.eq.2.or.istsb.eq.4) then
          rcth=pjr(0)
          mcth=1
          if(rcth.gt.coef(isub,13)) mcth=2
          if(rcth.gt.coef(isub,13)+coef(isub,14)) mcth=3
          if(rcth.gt.coef(isub,13)+coef(isub,14)+coef(isub,15)) mcth=4
          if(rcth.gt.coef(isub,13)+coef(isub,14)+coef(isub,15)+
     &    coef(isub,16)) mcth=5
          call pjkmap(3,mcth,pjr(0))
        endif
 
C...2 -> 3 : select pT1, phi1, pT2, phi2, y3 for 3 outgoing.
        if(istsb.eq.5) then
          call pjkmap(5,0,0d0)
          if(mint(51).ne.0) then
            if(mint(121).gt.1) call pjsave(2,iga)
            if(mfail.eq.1) then
              msti(61)=1
              return
            endif
            goto 100
          endif
        endif
 
C...Low-pT or multiple interactions (first semihard interaction).
      elseif(istsb.eq.9) then
        call pjmult(3)
        isub=mint(1)
 
C...Generate user-defined process: kinematics plus weight.
      elseif(istsb.eq.11) then
        msti(51)=0
        call pjupev(isub,sigs)
        if(nup.le.0) then
          mint(51)=2
          msti(51)=1
          if(mint(82).eq.1) then
            ngen(0,1)=ngen(0,1)-1
            ngen(0,2)=ngen(0,2)-1
            ngen(isub,1)=ngen(isub,1)-1
          endif
          if(mint(121).gt.1) call pjsave(2,iga)
          return
        endif
 
C...Construct 'trivial' kinematical variables needed.
        kfl1=kup(1,2)
        kfl2=kup(2,2)
        vint(41)=2d0*pup(1,4)/vint(1)
        vint(42)=2d0*pup(2,4)/vint(1)
        vint(21)=vint(41)*vint(42)
        vint(22)=0.5d0*log(vint(41)/vint(42))
        vint(44)=vint(21)*vint(2)
        vint(43)=sqrt(max(0d0,vint(44)))
        vint(56)=q2up(0)
        vint(55)=sqrt(max(0d0,vint(56)))
 
C...Construct other kinematical variables needed (approximately).
        vint(23)=0d0
        vint(26)=vint(21)
        vint(45)=-0.5d0*vint(44)
        vint(46)=-0.5d0*vint(44)
        vint(49)=vint(43)
        vint(50)=vint(44)
        vint(51)=vint(55)
        vint(52)=vint(56)
        vint(53)=vint(55)
        vint(54)=vint(56)
        vint(25)=0d0
        vint(48)=0d0
        do 240 iup=3,nup
          if(kup(iup,1).eq.1) vint(25)=vint(25)+2d0*(pup(iup,5)**2+
     &    pup(iup,1)**2+pup(iup,2)**2)/vint(1)
          if(kup(iup,1).eq.1) vint(48)=vint(48)+0.5d0*(pup(iup,1)**2+
     &    pup(iup,2)**2)
  240   continue
        vint(47)=sqrt(vint(48))
 
C...Calculate parton distribution weights.
        if(mint(47).ge.2) then
          do 260 i=3-min(2,mint(45)),min(2,mint(46))
            mint(105)=mint(102+i)
            mint(109)=mint(106+i)
            if(mstp(57).le.1) then
              call pjpdfu(mint(10+i),vint(40+i),q2up(0),xpq)
            else
              call pjpdfl(mint(10+i),vint(40+i),q2up(0),xpq)
            endif
            do 250 kfl=-25,25
              xsfx(i,kfl)=xpq(kfl)
  250       continue
  260     continue
        endif
      endif
 
C...Choose azimuthal angle.
      vint(24)=paru(2)*pjr(0)
 
C...Check against user cuts on kinematics at parton level.
      mint(51)=0
      if((isub.le.90.or.isub.gt.100).and.istsb.le.10) call pjklim(0)
      if(mint(51).ne.0) then
        if(mint(121).gt.1) call pjsave(2,iga)
        if(mfail.eq.1) then
          msti(61)=1
          return
        endif
        goto 100
      endif
      if(mint(82).eq.1.and.mstp(141).ge.1.and.istsb.le.10) then
        mcut=0
        if(msub(91)+msub(92)+msub(93)+msub(94)+msub(95).eq.0)
     &  call pjkcut(mcut)
        if(mcut.ne.0) then
          if(mint(121).gt.1) call pjsave(2,iga)
          if(mfail.eq.1) then
            msti(61)=1
            return
          endif
          goto 100
        endif
      endif
 
C...Calculate differential cross-section for different subprocesses.
      if(istsb.le.10) call pjsigh(nchn,sigs)
      sigsor=sigs
      siglpt=sigt(0,0,5)
 
C...Multiply cross-section by user-defined weights.
      if(mstp(173).eq.1) then
        sigs=parp(173)*sigs
        do 270 ichn=1,nchn
          sigh(ichn)=parp(173)*sigh(ichn)
  270   continue
        siglpt=parp(173)*siglpt
      endif
      wtxs=1d0
      sigswt=sigs
      vint(99)=1d0
      vint(100)=1d0
      if(mint(82).eq.1.and.mstp(142).ge.1) then
        if(isub.ne.96.and.msub(91)+msub(92)+msub(93)+msub(94)+
     &  msub(95).eq.0) call pjevwt(wtxs)
        sigswt=wtxs*sigs
        vint(99)=wtxs
        if(mstp(142).eq.1) vint(100)=1d0/wtxs
      endif
 
C...Calculations for Monte Carlo estimate of all cross-sections.
      if(mint(82).eq.1.and.isub.le.90.or.isub.ge.96) then
        if(mstp(142).le.1) then
          xsec(isub,2)=xsec(isub,2)+sigs
        else
          xsec(isub,2)=xsec(isub,2)+sigswt
        endif
      elseif(mint(82).eq.1) then
        xsec(isub,2)=xsec(isub,2)+sigs
      endif
      if((isub.eq.95.or.isub.eq.96).and.loop.eq.1.and.mint(82).eq.1)
     &xsec(97,2)=xsec(97,2)+siglpt
 
C...Multiple interactions: store results of cross-section calculation.
      if(mint(50).eq.1.and.mstp(82).ge.3) then
        vint(153)=sigsor
        call pjmult(4)
      endif
 
C...Check that weight not negative.
      viol=sigswt/xsec(isub,1)
      if(isub.eq.96.and.mstp(173).eq.1) viol=viol/parp(174)
      if(mstp(123).le.0) then
        if(viol.lt.-1d-3) then
          write(mstu(11),5000) viol,ngen(0,3)+1
          if(mstp(122).ge.1) write(mstu(11),5100) isub,vint(21),
     &    vint(22),vint(23),vint(26)
          stop
        endif
      else
        if(viol.lt.min(-1d-3,vint(109))) then
          vint(109)=viol
          write(mstu(11),5200) viol,ngen(0,3)+1
          if(mstp(122).ge.1) write(mstu(11),5100) isub,vint(21),
     &    vint(22),vint(23),vint(26)
        endif
      endif
 
C...Weighting using estimate of maximum of differential cross-section.
      if(mfail.eq.0) then
        if(viol.lt.pjr(0)) then
          if(mint(121).gt.1) call pjsave(2,iga)
          goto 100
        endif
      elseif(isub.ne.95.and.isub.ne.96) then
        if(viol.lt.pjr(0)) then
          msti(61)=1
          if(mint(121).gt.1) call pjsave(2,iga)
          return
        endif
      else
        ratnd=siglpt/xsec(95,1)
        if(loop.eq.1.and.ratnd.lt.pjr(0)) then
          msti(61)=1
          if(mint(121).gt.1) call pjsave(2,iga)
          return
        endif
        viol=viol/ratnd
        if(viol.lt.pjr(0)) then
          if(mint(121).gt.1) call pjsave(2,iga)
          goto 100
        endif
      endif
 
C...Check for possible violation of estimated maximum of differential
C...cross-section used in weighting.
      if(mstp(123).le.0) then
        if(viol.gt.1d0) then
          write(mstu(11),5300) viol,ngen(0,3)+1
          if(mstp(122).ge.2) write(mstu(11),5100) isub,vint(21),
     &    vint(22),vint(23),vint(26)
          stop
        endif
      elseif(mstp(123).eq.1) then
        if(viol.gt.vint(108)) then
          vint(108)=viol
          if(viol.gt.1d0) then
            mint(10)=1
            write(mstu(11),5400) viol,ngen(0,3)+1
            if(mstp(122).ge.2) write(mstu(11),5100) isub,vint(21),
     &      vint(22),vint(23),vint(26)
          endif
        endif
      elseif(viol.gt.vint(108)) then
        vint(108)=viol
        if(viol.gt.1d0) then
          mint(10)=1
          xdif=xsec(isub,1)*(viol-1d0)
          xsec(isub,1)=xsec(isub,1)+xdif
          if(msub(isub).eq.1.and.(isub.le.90.or.isub.gt.96))
     &    xsec(0,1)=xsec(0,1)+xdif
          write(mstu(11),5400) viol,ngen(0,3)+1
          if(mstp(122).ge.2) write(mstu(11),5100) isub,vint(21),
     &    vint(22),vint(23),vint(26)
          if(isub.le.9) then
            write(mstu(11),5500) isub,xsec(isub,1)
          elseif(isub.le.99) then
            write(mstu(11),5600) isub,xsec(isub,1)
          else
            write(mstu(11),5700) isub,xsec(isub,1)
          endif
          vint(108)=1d0
        endif
      endif
 
C...Multiple interactions: choose impact parameter.
      vint(148)=1d0
      if(mint(50).eq.1.and.(isub.le.90.or.isub.ge.96).and.
     &mstp(82).ge.3) then
        call pjmult(5)
        if(vint(150).lt.pjr(0)) then
          if(mint(121).gt.1) call pjsave(2,iga)
          if(mfail.eq.1) then
            msti(61)=1
            return
          endif
          goto 100
        endif
      endif
      if(mint(82).eq.1) ngen(0,2)=ngen(0,2)+1
      if(mint(82).eq.1.and.msub(95).eq.1) then
        if(isub.le.90.or.isub.ge.95) ngen(95,1)=ngen(95,1)+1
        if(isub.le.90.or.isub.ge.96) ngen(96,2)=ngen(96,2)+1
      endif
      if(isub.le.90.or.isub.ge.96) mint(31)=mint(31)+1
 
C...Choose flavour of reacting partons (and subprocess).
      if(istsb.ge.11) goto 290
      rsigs=sigs*pjr(0)
      qt2=vint(48)
      rqqbar=parp(87)*(1d0-(qt2/(qt2+(parp(88)*parp(82))**2))**2)
      if(isub.ne.95.and.(isub.ne.96.or.mstp(82).le.1.or.
     &pjr(0).gt.rqqbar)) then
        do 280 ichn=1,nchn
          kfl1=isig(ichn,1)
          kfl2=isig(ichn,2)
          mint(2)=isig(ichn,3)
          rsigs=rsigs-sigh(ichn)
          if(rsigs.le.0d0) goto 290
  280   continue
 
C...Multiple interactions: choose qqbar preferentially at small pT.
      elseif(isub.eq.96) then
        mint(105)=mint(103)
        mint(109)=mint(107)
        call pjspli(mint(11),21,kfl1,kfldum)
        mint(105)=mint(104)
        mint(109)=mint(108)
        call pjspli(mint(12),21,kfl2,kfldum)
        mint(1)=11
        mint(2)=1
        if(kfl1.eq.kfl2.and.pjr(0).lt.0.5d0) mint(2)=2
 
C...Low-pT: choose string drawing configuration.
      else
        kfl1=21
        kfl2=21
        rsigs=6d0*pjr(0)
        mint(2)=1
        if(rsigs.gt.1d0) mint(2)=2
        if(rsigs.gt.2d0) mint(2)=3
      endif
 
C...Reassign QCD process. Partons before initial state radiation.
  290 if(mint(2).gt.10) then
        mint(1)=mint(2)/10
        mint(2)=mod(mint(2),10)
      endif
      if(mint(82).eq.1.and.mstp(111).ge.0) ngen(mint(1),2)=
     &ngen(mint(1),2)+1
      mint(15)=kfl1
      mint(16)=kfl2
      mint(13)=mint(15)
      mint(14)=mint(16)
      vint(141)=vint(41)
      vint(142)=vint(42)
      vint(151)=0d0
      vint(152)=0d0
 
C...Calculate x value of photon for parton inside photon inside e.
      do 320 jt=1,2
        mint(18+jt)=0
        vint(154+jt)=0d0
        mspli=0
        if(jt.eq.1.and.mint(43).le.2) mspli=1
        if(jt.eq.2.and.mod(mint(43),2).eq.1) mspli=1
        if(iabs(mint(14+jt)).le.8.or.mint(14+jt).eq.21) mspli=mspli+1
        if(mspli.eq.2) then
          kflh=mint(14+jt)
          xhrd=vint(140+jt)
          q2hrd=vint(54)
          mint(105)=mint(102+jt)
          mint(109)=mint(106+jt)
          if(mstp(57).le.1) then
            call pjpdfu(22,xhrd,q2hrd,xpq)
          else
            call pjpdfl(22,xhrd,q2hrd,xpq)
          endif
          wtmx=4d0*xpq(kflh)
          if(mstp(13).eq.2) then
            q2pms=q2hrd/pmas(11,1)**2
            wtmx=wtmx*log(max(2d0,q2pms*(1d0-xhrd)/xhrd**2))
          endif
  300     xe=xhrd**pjr(0)
          xg=min(0.999999d0,xhrd/xe)
          if(mstp(57).le.1) then
            call pjpdfu(22,xg,q2hrd,xpq)
          else
            call pjpdfl(22,xg,q2hrd,xpq)
          endif
          wt=(1d0+(1d0-xe)**2)*xpq(kflh)
          if(mstp(13).eq.2) wt=wt*log(max(2d0,q2pms*(1d0-xe)/xe**2))
          if(wt.lt.pjr(0)*wtmx) goto 300
          mint(18+jt)=1
          vint(154+jt)=xe
          do 310 kfls=-25,25
            xsfx(jt,kfls)=xpq(kfls)
  310     continue
        endif
  320 continue
 
C...Pick scale where photon is resolved.
      if(mint(107).eq.3) vint(283)=parp(15)**2*
     &(vint(54)/parp(15)**2)**pjr(0)
      if(mint(108).eq.3) vint(284)=parp(15)**2*
     &(vint(54)/parp(15)**2)**pjr(0)
      if(mint(121).gt.1) call pjsave(2,iga)
 
C...Format statements for differential cross-section maximum violations.
 5000 format(/1x,'Error: negative cross-section fraction',1p,d11.3,1x,
     &'in event',1x,i7,'D0'/1x,'Execution stopped!')
 5100 format(1x,'ISUB = ',i3,'; Point of violation:'/1x,'tau =',1p,
     &d11.3,', y* =',d11.3,', cthe = ',0p,f11.7,', tau'' =',1p,d11.3)
 5200 format(/1x,'Warning: negative cross-section fraction',1p,d11.3,1x,
     &'in event',1x,i7)
 5300 format(/1x,'Error: maximum violated by',1p,d11.3,1x,
     &'in event',1x,i7,'D0'/1x,'Execution stopped!')
 5400 format(/1x,'Advisory warning: maximum violated by',1p,d11.3,1x,
     &'in event',1x,i7)
 5500 format(1x,'XSEC(',i1,',1) increased to',1p,d11.3)
 5600 format(1x,'XSEC(',i2,',1) increased to',1p,d11.3)
 5700 format(1x,'XSEC(',i3,',1) increased to',1p,d11.3)
 
      return
      end
 
C*********************************************************************
 
C...PYSCAT
C...Finds outgoing flavours and event type; sets up the kinematics
C...and colour flow of the hard scattering
 
      subroutine pjscat
 
C...Double precision and integer declarations
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
      common/pjuppr/nup,kup(20,7),nfup,ifup(10,2),pup(20,5),q2up(0:10)
      common/pjssmt/zmix(4,4),umix(2,2),vmix(2,2),smz(4),smw(2),
     &sfmix(16,4)
      save /jyjets/,/jydat1/,/jydat2/,/jydat3/,/pjsubs/,/pjpars/,
     &/pjint1/,/pjint2/,/pjint3/,/pjint4/,/pjint5/,/pjuppr/,/pjssmt/
C...Local arrays and saved variables
      dimension wdtp(0:200),wdte(0:200,0:5),pmq(2),z(2),cthe(2),
     &phi(2),kuppo(20),vintsv(41:66)
      save vintsv
 
C...Read out process
      isub=mint(1)
      isubsv=isub
 
C...Restore information for low-pT processes
      if(isub.eq.95.and.mint(57).ge.1) then
        do 100 j=41,66
        vint(j)=vintsv(j)
  100   continue
      endif
 
C...Convert H' or A process into equivalent H one
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
 
C...Choice of subprocess, number of documentation lines
      idoc=6+iset(isub)
      if(isub.eq.95) idoc=8
      if(iset(isub).eq.5) idoc=9
      if(iset(isub).eq.11) idoc=4+nup
      mint(3)=idoc-6
      if(idoc.ge.9.and.iset(isub).le.4) idoc=idoc+2
      mint(4)=idoc
      ipu1=mint(84)+1
      ipu2=mint(84)+2
      ipu3=mint(84)+3
      ipu4=mint(84)+4
      ipu5=mint(84)+5
      ipu6=mint(84)+6
 
C...Reset K, P and V vectors. Store incoming particles
      do 120 jt=1,mstp(126)+20
        i=mint(83)+jt
        do 110 j=1,5
          k(i,j)=0
          p(i,j)=0d0
          v(i,j)=0d0
  110   continue
  120 continue
      do 140 jt=1,2
        i=mint(83)+jt
        k(i,1)=21
        k(i,2)=mint(10+jt)
        do 130 j=1,5
          p(i,j)=vint(285+5*jt+j)
  130   continue
  140 continue
      mint(6)=2
      kfres=0
 
C...Store incoming partons in their CM-frame
      sh=vint(44)
      shr=sqrt(sh)
      shp=vint(26)*vint(2)
      shpr=sqrt(shp)
      shuser=shr
      if(iset(isub).ge.3.and.iset(isub).le.5) shuser=shpr
      do 150 jt=1,2
        i=mint(84)+jt
        k(i,1)=14
        k(i,2)=mint(14+jt)
        k(i,3)=mint(83)+2+jt
        p(i,3)=0.5d0*shuser*(-1d0)**(jt-1)
        p(i,4)=0.5d0*shuser
  150 continue
 
C...Copy incoming partons to documentation lines
      do 170 jt=1,2
        i1=mint(83)+4+jt
        i2=mint(84)+jt
        k(i1,1)=21
        k(i1,2)=k(i2,2)
        k(i1,3)=i1-2
        do 160 j=1,5
          p(i1,j)=p(i2,j)
  160   continue
  170 continue
 
C...Choose new quark/lepton flavour for relevant annihilation graphs
      if(isub.eq.12.or.isub.eq.53.or.isub.eq.54.or.isub.eq.58) then
        iglga=21
        if(isub.eq.58) iglga=22
        call pjwidt(iglga,sh,wdtp,wdte)
  180   rkfl=(wdte(0,1)+wdte(0,2)+wdte(0,4))*pjr(0)
        do 190 i=1,mdcy(iglga,3)
          kflf=kfdp(i+mdcy(iglga,2)-1,1)
          rkfl=rkfl-(wdte(i,1)+wdte(i,2)+wdte(i,4))
          if(rkfl.le.0d0) goto 200
  190   continue
  200   continue
        if(isub.eq.12.and.mstp(5).eq.1.and.iabs(mint(15)).le.2.and.
     &  iabs(kflf).ge.3) then
          facqqb=vint(58)**2*4d0/9d0*(vint(45)**2+vint(46)**2)/
     &    vint(44)**2
          faccib=vint(46)**2/paru(155)**4
          if(facqqb/(facqqb+faccib).lt.pjr(0)) goto 180
        elseif(isub.eq.54) then
          if((kchg(jamcomp(kflf),1)/2d0)**2.lt.pjr(0)) goto 180
        elseif(isub.eq.58) then
          if((kchg(jamcomp(kflf),1)/3d0)**2.lt.pjr(0)) goto 180
        endif
      endif
 
C...Final state flavours and colour flow: default values
      js=1
      mint(21)=mint(15)
      mint(22)=mint(16)
      mint(23)=0
      mint(24)=0
      kcc=20
      kcs=isign(1,mint(15))
 
      if(iset(isub).eq.11) then
C...User-defined processes: find products
        irup=0
        do 210 iup=3,nup
          if(kup(iup,1).ne.1) then
          elseif(irup.le.5) then
            irup=irup+1
            mint(20+irup)=kup(iup,2)
          endif
  210   continue
 
      elseif(isub.le.10) then
        if(isub.eq.1) then
C...f + fbar -> gamma*/Z0
          kfres=23
 
        elseif(isub.eq.2) then
C...f + fbar' -> W+/-
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          kfres=isign(24,kch1+kch2)
 
        elseif(isub.eq.3) then
C...f + fbar -> h0 (or H0, or A0)
          kfres=kfhigg
 
        elseif(isub.eq.4) then
C...gamma + W+/- -> W+/-
 
        elseif(isub.eq.5) then
C...Z0 + Z0 -> h0
          xh=sh/shp
          mint(21)=mint(15)
          mint(22)=mint(16)
          pmq(1)=pjmass(mint(21))
          pmq(2)=pjmass(mint(22))
  220     jt=int(1.5d0+pjr(0))
          zmin=2d0*pmq(jt)/shpr
          zmax=1d0-pmq(3-jt)/shpr-(sh-pmq(jt)**2)/
     &    (shpr*(shpr-pmq(3-jt)))
          zmax=min(1d0-xh,zmax)
          z(jt)=zmin+(zmax-zmin)*pjr(0)
          if(-1d0+(1d0+xh)/(1d0-z(jt))-xh/(1d0-z(jt))**2.lt.
     &    (1d0-xh)**2/(4d0*xh)*pjr(0)) goto 220
          sqc1=1d0-4d0*pmq(jt)**2/(z(jt)**2*shp)
          if(sqc1.lt.1.d-8) goto 220
          c1=sqrt(sqc1)
          c2=1d0+2d0*(pmas(23,1)**2-pmq(jt)**2)/(z(jt)*shp)
          cthe(jt)=(c2-(c2**2-c1**2)/(c2+(2d0*pjr(0)-1d0)*c1))/c1
          cthe(jt)=min(1d0,max(-1d0,cthe(jt)))
          z(3-jt)=1d0-xh/(1d0-z(jt))
          sqc1=1d0-4d0*pmq(3-jt)**2/(z(3-jt)**2*shp)
          if(sqc1.lt.1.d-8) goto 220
          c1=sqrt(sqc1)
          c2=1d0+2d0*(pmas(23,1)**2-pmq(3-jt)**2)/(z(3-jt)*shp)
          cthe(3-jt)=(c2-(c2**2-c1**2)/(c2+(2d0*pjr(0)-1d0)*c1))/c1
          cthe(3-jt)=min(1d0,max(-1d0,cthe(3-jt)))
          phir=paru(2)*pjr(0)
          cphi=cos(phir)
          ang=cthe(1)*cthe(2)-sqrt(1d0-cthe(1)**2)*
     &    sqrt(1d0-cthe(2)**2)*cphi
          z1=2d0-z(jt)
          z2=ang*sqrt(z(jt)**2-4d0*pmq(jt)**2/shp)
          z3=1d0-z(jt)-xh+(pmq(1)**2+pmq(2)**2)/shp
          z(3-jt)=2d0/(z1**2-z2**2)*(z1*z3+z2*sqrt(z3**2-(z1**2-z2**2)*
     &    pmq(3-jt)**2/shp))
          zmin=2d0*pmq(3-jt)/shpr
          zmax=1d0-pmq(jt)/shpr-(sh-pmq(3-jt)**2)/(shpr*(shpr-pmq(jt)))
          zmax=min(1d0-xh,zmax)
          if(z(3-jt).lt.zmin.or.z(3-jt).gt.zmax) goto 220
          kcc=22
          kfres=25
 
        elseif(isub.eq.6) then
C...Z0 + W+/- -> W+/-
 
        elseif(isub.eq.7) then
C...W+ + W- -> Z0
 
        elseif(isub.eq.8) then
C...W+ + W- -> h0
          xh=sh/shp
  230     do 260 jt=1,2
            i=mint(14+jt)
            ia=iabs(i)
            if(ia.le.10) then
              rvckm=vint(180+i)*pjr(0)
              do 240 j=1,mstp(1)
                ib=2*j-1+mod(ia,2)
                ipm=(5-isign(1,i))/2
                idc=j+mdcy(ia,2)+2
                if(mdme(idc,1).ne.1.and.mdme(idc,1).ne.ipm) goto 240
                mint(20+jt)=isign(ib,i)
                rvckm=rvckm-vckm((ia+1)/2,(ib+1)/2)
                if(rvckm.le.0d0) goto 250
  240         continue
            else
              ib=2*((ia+1)/2)-1+mod(ia,2)
              mint(20+jt)=isign(ib,i)
            endif
  250       pmq(jt)=pjmass(mint(20+jt))
  260     continue
          jt=int(1.5d0+pjr(0))
          zmin=2d0*pmq(jt)/shpr
          zmax=1d0-pmq(3-jt)/shpr-(sh-pmq(jt)**2)/
     &    (shpr*(shpr-pmq(3-jt)))
          zmax=min(1d0-xh,zmax)
          if(zmin.ge.zmax) goto 230
          z(jt)=zmin+(zmax-zmin)*pjr(0)
          if(-1d0+(1d0+xh)/(1d0-z(jt))-xh/(1d0-z(jt))**2.lt.
     &    (1d0-xh)**2/(4d0*xh)*pjr(0)) goto 230
          sqc1=1d0-4d0*pmq(jt)**2/(z(jt)**2*shp)
          if(sqc1.lt.1.d-8) goto 230
          c1=sqrt(sqc1)
          c2=1d0+2d0*(pmas(24,1)**2-pmq(jt)**2)/(z(jt)*shp)
          cthe(jt)=(c2-(c2**2-c1**2)/(c2+(2d0*pjr(0)-1d0)*c1))/c1
          cthe(jt)=min(1d0,max(-1d0,cthe(jt)))
          z(3-jt)=1d0-xh/(1d0-z(jt))
          sqc1=1d0-4d0*pmq(3-jt)**2/(z(3-jt)**2*shp)
          if(sqc1.lt.1.d-8) goto 230
          c1=sqrt(sqc1)
          c2=1d0+2d0*(pmas(24,1)**2-pmq(3-jt)**2)/(z(3-jt)*shp)
          cthe(3-jt)=(c2-(c2**2-c1**2)/(c2+(2d0*pjr(0)-1d0)*c1))/c1
          cthe(3-jt)=min(1d0,max(-1d0,cthe(3-jt)))
          phir=paru(2)*pjr(0)
          cphi=cos(phir)
          ang=cthe(1)*cthe(2)-sqrt(1d0-cthe(1)**2)*
     &    sqrt(1d0-cthe(2)**2)*cphi
          z1=2d0-z(jt)
          z2=ang*sqrt(z(jt)**2-4d0*pmq(jt)**2/shp)
          z3=1d0-z(jt)-xh+(pmq(1)**2+pmq(2)**2)/shp
          z(3-jt)=2d0/(z1**2-z2**2)*(z1*z3+z2*sqrt(z3**2-(z1**2-z2**2)*
     &    pmq(3-jt)**2/shp))
          zmin=2d0*pmq(3-jt)/shpr
          zmax=1d0-pmq(jt)/shpr-(sh-pmq(3-jt)**2)/(shpr*(shpr-pmq(jt)))
          zmax=min(1d0-xh,zmax)
          if(z(3-jt).lt.zmin.or.z(3-jt).gt.zmax) goto 230
          kcc=22
          kfres=25
 
        elseif(isub.eq.10) then
C...f + f' -> f + f' (gamma/Z/W exchange); th = (p(f)-p(f))**2
          if(mint(2).eq.1) then
            kcc=22
          else
C...W exchange: need to mix flavours according to CKM matrix
            do 280 jt=1,2
              i=mint(14+jt)
              ia=iabs(i)
              if(ia.le.10) then
                rvckm=vint(180+i)*pjr(0)
                do 270 j=1,mstp(1)
                  ib=2*j-1+mod(ia,2)
                  ipm=(5-isign(1,i))/2
                  idc=j+mdcy(ia,2)+2
                  if(mdme(idc,1).ne.1.and.mdme(idc,1).ne.ipm) goto 270
                  mint(20+jt)=isign(ib,i)
                  rvckm=rvckm-vckm((ia+1)/2,(ib+1)/2)
                  if(rvckm.le.0d0) goto 280
  270           continue
              else
                ib=2*((ia+1)/2)-1+mod(ia,2)
                mint(20+jt)=isign(ib,i)
              endif
  280       continue
            kcc=22
          endif
        endif
 
      elseif(isub.le.20) then
        if(isub.eq.11) then
C...f + f' -> f + f' (g exchange); th = (p(f)-p(f))**2
          kcc=mint(2)
          if(mint(15)*mint(16).lt.0) kcc=kcc+2
 
        elseif(isub.eq.12) then
C...f + fbar -> f' + fbar'; th = (p(f)-p(f'))**2
          mint(21)=isign(kflf,mint(15))
          mint(22)=-mint(21)
          kcc=4
 
        elseif(isub.eq.13) then
C...f + fbar -> g + g; th arbitrary
          mint(21)=21
          mint(22)=21
          kcc=mint(2)+4
 
        elseif(isub.eq.14) then
C...f + fbar -> g + gamma; th arbitrary
          if(pjr(0).gt.0.5d0) js=2
          mint(20+js)=21
          mint(23-js)=22
          kcc=17+js
 
        elseif(isub.eq.15) then
C...f + fbar -> g + Z0; th arbitrary
          if(pjr(0).gt.0.5d0) js=2
          mint(20+js)=21
          mint(23-js)=23
          kcc=17+js
 
        elseif(isub.eq.16) then
C...f + fbar' -> g + W+/-; th = (p(f)-p(W-))**2 or (p(fbar')-p(W+))**2
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          if(mint(15)*(kch1+kch2).lt.0) js=2
          mint(20+js)=21
          mint(23-js)=isign(24,kch1+kch2)
          kcc=17+js
 
        elseif(isub.eq.17) then
C...f + fbar -> g + h0; th arbitrary
          if(pjr(0).gt.0.5d0) js=2
          mint(20+js)=21
          mint(23-js)=25
          kcc=17+js
 
        elseif(isub.eq.18) then
C...f + fbar -> gamma + gamma; th arbitrary
          mint(21)=22
          mint(22)=22
 
        elseif(isub.eq.19) then
C...f + fbar -> gamma + Z0; th arbitrary
          if(pjr(0).gt.0.5d0) js=2
          mint(20+js)=22
          mint(23-js)=23
 
        elseif(isub.eq.20) then
C...f + fbar' -> gamma + W+/-; th = (p(f)-p(W-))**2 or
C...(p(fbar')-p(W+))**2
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          if(mint(15)*(kch1+kch2).lt.0) js=2
          mint(20+js)=22
          mint(23-js)=isign(24,kch1+kch2)
        endif
 
      elseif(isub.le.30) then
        if(isub.eq.21) then
C...f + fbar -> gamma + h0; th arbitrary
          if(pjr(0).gt.0.5d0) js=2
          mint(20+js)=22
          mint(23-js)=25
 
        elseif(isub.eq.22) then
C...f + fbar -> Z0 + Z0; th arbitrary
          mint(21)=23
          mint(22)=23
 
        elseif(isub.eq.23) then
C...f + fbar' -> Z0 + W+/-; th = (p(f)-p(W-))**2 or (p(fbar')-p(W+))**2
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          if(mint(15)*(kch1+kch2).lt.0) js=2
          mint(20+js)=23
          mint(23-js)=isign(24,kch1+kch2)
 
        elseif(isub.eq.24) then
C...f + fbar -> Z0 + h0 (or H0, or A0); th arbitrary
          if(pjr(0).gt.0.5d0) js=2
          mint(20+js)=23
          mint(23-js)=kfhigg
 
        elseif(isub.eq.25) then
C...f + fbar -> W+ + W-; th = (p(f)-p(W-))**2
          mint(21)=-isign(24,mint(15))
          mint(22)=-mint(21)
 
        elseif(isub.eq.26) then
C...f + fbar' -> W+/- + h0 (or H0, or A0);
C...th = (p(f)-p(W-))**2 or (p(fbar')-p(W+))**2
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          if(mint(15)*(kch1+kch2).gt.0) js=2
          mint(20+js)=isign(24,kch1+kch2)
          mint(23-js)=kfhigg
 
        elseif(isub.eq.27) then
C...f + fbar -> h0 + h0
 
        elseif(isub.eq.28) then
C...f + g -> f + g; th = (p(f)-p(f))**2
          kcc=mint(2)+6
          if(mint(15).eq.21) kcc=kcc+2
          if(mint(15).ne.21) kcs=isign(1,mint(15))
          if(mint(16).ne.21) kcs=isign(1,mint(16))
 
        elseif(isub.eq.29) then
C...f + g -> f + gamma; th = (p(f)-p(f))**2
          if(mint(15).eq.21) js=2
          mint(23-js)=22
          kcc=15+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.30) then
C...f + g -> f + Z0; th = (p(f)-p(f))**2
          if(mint(15).eq.21) js=2
          mint(23-js)=23
          kcc=15+js
          kcs=isign(1,mint(14+js))
        endif
 
      elseif(isub.le.40) then
        if(isub.eq.31) then
C...f + g -> f' + W+/-; th = (p(f)-p(f'))**2; choose flavour f'
          if(mint(15).eq.21) js=2
          i=mint(14+js)
          ia=iabs(i)
          mint(23-js)=isign(24,kchg(ia,1)*i)
          rvckm=vint(180+i)*pjr(0)
          do 290 j=1,mstp(1)
            ib=2*j-1+mod(ia,2)
            ipm=(5-isign(1,i))/2
            idc=j+mdcy(ia,2)+2
            if(mdme(idc,1).ne.1.and.mdme(idc,1).ne.ipm) goto 290
            mint(20+js)=isign(ib,i)
            rvckm=rvckm-vckm((ia+1)/2,(ib+1)/2)
            if(rvckm.le.0d0) goto 300
  290     continue
  300     kcc=15+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.32) then
C...f + g -> f + h0; th = (p(f)-p(f))**2
          if(mint(15).eq.21) js=2
          mint(23-js)=25
          kcc=15+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.33) then
C...f + gamma -> f + g; th=(p(f)-p(f))**2
          if(mint(15).eq.22) js=2
          mint(23-js)=21
          kcc=24+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.34) then
C...f + gamma -> f + gamma; th=(p(f)-p(f))**2
          if(mint(15).eq.22) js=2
          kcc=22
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.35) then
C...f + gamma -> f + Z0; th=(p(f)-p(f))**2
          if(mint(15).eq.22) js=2
          mint(23-js)=23
          kcc=22
 
        elseif(isub.eq.36) then
C...f + gamma -> f' + W+/-; th=(p(f)-p(f'))**2
          if(mint(15).eq.22) js=2
          i=mint(14+js)
          ia=iabs(i)
          mint(23-js)=isign(24,kchg(ia,1)*i)
          if(ia.le.10) then
            rvckm=vint(180+i)*pjr(0)
            do 310 j=1,mstp(1)
              ib=2*j-1+mod(ia,2)
              ipm=(5-isign(1,i))/2
              idc=j+mdcy(ia,2)+2
              if(mdme(idc,1).ne.1.and.mdme(idc,1).ne.ipm) goto 310
              mint(20+js)=isign(ib,i)
              rvckm=rvckm-vckm((ia+1)/2,(ib+1)/2)
              if(rvckm.le.0d0) goto 320
  310       continue
          else
            ib=2*((ia+1)/2)-1+mod(ia,2)
            mint(20+js)=isign(ib,i)
          endif
  320     kcc=22
 
        elseif(isub.eq.37) then
C...f + gamma -> f + h0
 
        elseif(isub.eq.38) then
C...f + Z0 -> f + g
 
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
C...f + W+/- -> f' + g
 
        elseif(isub.eq.44) then
C...f + W+/- -> f' + gamma
 
        elseif(isub.eq.45) then
C...f + W+/- -> f' + Z0
 
        elseif(isub.eq.46) then
C...f + W+/- -> f' + W+/-
 
        elseif(isub.eq.47) then
C...f + W+/- -> f' + h0
 
        elseif(isub.eq.48) then
C...f + h0 -> f + g
 
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
C...g + g -> f + fbar; th arbitrary
          kcs=(-1)**int(1.5d0+pjr(0))
          mint(21)=isign(kflf,kcs)
          mint(22)=-mint(21)
          kcc=mint(2)+10
 
        elseif(isub.eq.54) then
C...g + gamma -> f + fbar; th arbitrary
          kcs=(-1)**int(1.5d0+pjr(0))
          mint(21)=isign(kflf,kcs)
          mint(22)=-mint(21)
          kcc=27
          if(mint(16).eq.21) kcc=28
 
        elseif(isub.eq.55) then
C...g + Z0 -> f + fbar
 
        elseif(isub.eq.56) then
C...g + W+/- -> f + fbar'
 
        elseif(isub.eq.57) then
C...g + h0 -> f + fbar
 
        elseif(isub.eq.58) then
C...gamma + gamma -> f + fbar; th arbitrary
          kcs=(-1)**int(1.5d0+pjr(0))
          mint(21)=isign(kflf,kcs)
          mint(22)=-mint(21)
          kcc=21
 
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
C...g + g -> g + g; th arbitrary
          kcc=mint(2)+12
          kcs=(-1)**int(1.5d0+pjr(0))
 
        elseif(isub.eq.69) then
C...gamma + gamma -> W+ + W-; th arbitrary
          mint(21)=24
          mint(22)=-24
          kcc=21
 
        elseif(isub.eq.70) then
C...gamma + W+/- -> Z0 + W+/-; th=(p(W)-p(W))**2
          if(mint(15).eq.22) mint(21)=23
          if(mint(16).eq.22) mint(22)=23
          kcc=21
        endif
 
      elseif(isub.le.80) then
        if(isub.eq.71.or.isub.eq.72) then
C...Z0 + Z0 -> Z0 + Z0; Z0 + Z0 -> W+ + W-
          xh=sh/shp
          mint(21)=mint(15)
          mint(22)=mint(16)
          pmq(1)=pjmass(mint(21))
          pmq(2)=pjmass(mint(22))
  330     jt=int(1.5d0+pjr(0))
          zmin=2d0*pmq(jt)/shpr
          zmax=1d0-pmq(3-jt)/shpr-(sh-pmq(jt)**2)/
     &    (shpr*(shpr-pmq(3-jt)))
          zmax=min(1d0-xh,zmax)
          z(jt)=zmin+(zmax-zmin)*pjr(0)
          if(-1d0+(1d0+xh)/(1d0-z(jt))-xh/(1d0-z(jt))**2.lt.
     &    (1d0-xh)**2/(4d0*xh)*pjr(0)) goto 330
          sqc1=1d0-4d0*pmq(jt)**2/(z(jt)**2*shp)
          if(sqc1.lt.1.d-8) goto 330
          c1=sqrt(sqc1)
          c2=1d0+2d0*(pmas(23,1)**2-pmq(jt)**2)/(z(jt)*shp)
          cthe(jt)=(c2-(c2**2-c1**2)/(c2+(2d0*pjr(0)-1d0)*c1))/c1
          cthe(jt)=min(1d0,max(-1d0,cthe(jt)))
          z(3-jt)=1d0-xh/(1d0-z(jt))
          sqc1=1d0-4d0*pmq(3-jt)**2/(z(3-jt)**2*shp)
          if(sqc1.lt.1.d-8) goto 330
          c1=sqrt(sqc1)
          c2=1d0+2d0*(pmas(23,1)**2-pmq(3-jt)**2)/(z(3-jt)*shp)
          cthe(3-jt)=(c2-(c2**2-c1**2)/(c2+(2d0*pjr(0)-1d0)*c1))/c1
          cthe(3-jt)=min(1d0,max(-1d0,cthe(3-jt)))
          phir=paru(2)*pjr(0)
          cphi=cos(phir)
          ang=cthe(1)*cthe(2)-sqrt(1d0-cthe(1)**2)*
     &    sqrt(1d0-cthe(2)**2)*cphi
          z1=2d0-z(jt)
          z2=ang*sqrt(z(jt)**2-4d0*pmq(jt)**2/shp)
          z3=1d0-z(jt)-xh+(pmq(1)**2+pmq(2)**2)/shp
          z(3-jt)=2d0/(z1**2-z2**2)*(z1*z3+z2*sqrt(z3**2-(z1**2-z2**2)*
     &    pmq(3-jt)**2/shp))
          zmin=2d0*pmq(3-jt)/shpr
          zmax=1d0-pmq(jt)/shpr-(sh-pmq(3-jt)**2)/(shpr*(shpr-pmq(jt)))
          zmax=min(1d0-xh,zmax)
          if(z(3-jt).lt.zmin.or.z(3-jt).gt.zmax) goto 330
          kcc=22
 
        elseif(isub.eq.73) then
C...Z0 + W+/- -> Z0 + W+/-
          js=mint(2)
          xh=sh/shp
  340     jt=3-mint(2)
          i=mint(14+jt)
          ia=iabs(i)
          if(ia.le.10) then
            rvckm=vint(180+i)*pjr(0)
            do 350 j=1,mstp(1)
              ib=2*j-1+mod(ia,2)
              ipm=(5-isign(1,i))/2
              idc=j+mdcy(ia,2)+2
              if(mdme(idc,1).ne.1.and.mdme(idc,1).ne.ipm) goto 350
              mint(20+jt)=isign(ib,i)
              rvckm=rvckm-vckm((ia+1)/2,(ib+1)/2)
              if(rvckm.le.0d0) goto 360
  350       continue
          else
            ib=2*((ia+1)/2)-1+mod(ia,2)
            mint(20+jt)=isign(ib,i)
          endif
  360     pmq(jt)=pjmass(mint(20+jt))
          mint(23-jt)=mint(17-jt)
          pmq(3-jt)=pjmass(mint(23-jt))
          jt=int(1.5d0+pjr(0))
          zmin=2d0*pmq(jt)/shpr
          zmax=1d0-pmq(3-jt)/shpr-(sh-pmq(jt)**2)/
     &    (shpr*(shpr-pmq(3-jt)))
          zmax=min(1d0-xh,zmax)
          if(zmin.ge.zmax) goto 340
          z(jt)=zmin+(zmax-zmin)*pjr(0)
          if(-1d0+(1d0+xh)/(1d0-z(jt))-xh/(1d0-z(jt))**2.lt.
     &    (1d0-xh)**2/(4d0*xh)*pjr(0)) goto 340
          sqc1=1d0-4d0*pmq(jt)**2/(z(jt)**2*shp)
          if(sqc1.lt.1.d-8) goto 340
          c1=sqrt(sqc1)
          c2=1d0+2d0*(pmas(23,1)**2-pmq(jt)**2)/(z(jt)*shp)
          cthe(jt)=(c2-(c2**2-c1**2)/(c2+(2d0*pjr(0)-1d0)*c1))/c1
          cthe(jt)=min(1d0,max(-1d0,cthe(jt)))
          z(3-jt)=1d0-xh/(1d0-z(jt))
          sqc1=1d0-4d0*pmq(3-jt)**2/(z(3-jt)**2*shp)
          if(sqc1.lt.1.d-8) goto 340
          c1=sqrt(sqc1)
          c2=1d0+2d0*(pmas(23,1)**2-pmq(3-jt)**2)/(z(3-jt)*shp)
          cthe(3-jt)=(c2-(c2**2-c1**2)/(c2+(2d0*pjr(0)-1d0)*c1))/c1
          cthe(3-jt)=min(1d0,max(-1d0,cthe(3-jt)))
          phir=paru(2)*pjr(0)
          cphi=cos(phir)
          ang=cthe(1)*cthe(2)-sqrt(1d0-cthe(1)**2)*
     &    sqrt(1d0-cthe(2)**2)*cphi
          z1=2d0-z(jt)
          z2=ang*sqrt(z(jt)**2-4d0*pmq(jt)**2/shp)
          z3=1d0-z(jt)-xh+(pmq(1)**2+pmq(2)**2)/shp
          z(3-jt)=2d0/(z1**2-z2**2)*(z1*z3+z2*sqrt(z3**2-(z1**2-z2**2)*
     &    pmq(3-jt)**2/shp))
          zmin=2d0*pmq(3-jt)/shpr
          zmax=1d0-pmq(jt)/shpr-(sh-pmq(3-jt)**2)/(shpr*(shpr-pmq(jt)))
          zmax=min(1d0-xh,zmax)
          if(z(3-jt).lt.zmin.or.z(3-jt).gt.zmax) goto 340
          kcc=22
 
        elseif(isub.eq.74) then
C...Z0 + h0 -> Z0 + h0
 
        elseif(isub.eq.75) then
C...W+ + W- -> gamma + gamma
 
        elseif(isub.eq.76.or.isub.eq.77) then
C...W+ + W- -> Z0 + Z0; W+ + W- -> W+ + W-
          xh=sh/shp
  370     do 400 jt=1,2
            i=mint(14+jt)
            ia=iabs(i)
            if(ia.le.10) then
              rvckm=vint(180+i)*pjr(0)
              do 380 j=1,mstp(1)
                ib=2*j-1+mod(ia,2)
                ipm=(5-isign(1,i))/2
                idc=j+mdcy(ia,2)+2
                if(mdme(idc,1).ne.1.and.mdme(idc,1).ne.ipm) goto 380
                mint(20+jt)=isign(ib,i)
                rvckm=rvckm-vckm((ia+1)/2,(ib+1)/2)
                if(rvckm.le.0d0) goto 390
  380         continue
            else
              ib=2*((ia+1)/2)-1+mod(ia,2)
              mint(20+jt)=isign(ib,i)
            endif
  390       pmq(jt)=pjmass(mint(20+jt))
  400     continue
          jt=int(1.5d0+pjr(0))
          zmin=2d0*pmq(jt)/shpr
          zmax=1d0-pmq(3-jt)/shpr-(sh-pmq(jt)**2)/
     &    (shpr*(shpr-pmq(3-jt)))
          zmax=min(1d0-xh,zmax)
          if(zmin.ge.zmax) goto 370
          z(jt)=zmin+(zmax-zmin)*pjr(0)
          if(-1d0+(1d0+xh)/(1d0-z(jt))-xh/(1d0-z(jt))**2.lt.
     &    (1d0-xh)**2/(4d0*xh)*pjr(0)) goto 370
          sqc1=1d0-4d0*pmq(jt)**2/(z(jt)**2*shp)
          if(sqc1.lt.1.d-8) goto 370
          c1=sqrt(sqc1)
          c2=1d0+2d0*(pmas(24,1)**2-pmq(jt)**2)/(z(jt)*shp)
          cthe(jt)=(c2-(c2**2-c1**2)/(c2+(2d0*pjr(0)-1d0)*c1))/c1
          cthe(jt)=min(1d0,max(-1d0,cthe(jt)))
          z(3-jt)=1d0-xh/(1d0-z(jt))
          sqc1=1d0-4d0*pmq(3-jt)**2/(z(3-jt)**2*shp)
          if(sqc1.lt.1.d-8) goto 370
          c1=sqrt(sqc1)
          c2=1d0+2d0*(pmas(24,1)**2-pmq(3-jt)**2)/(z(3-jt)*shp)
          cthe(3-jt)=(c2-(c2**2-c1**2)/(c2+(2d0*pjr(0)-1d0)*c1))/c1
          cthe(3-jt)=min(1d0,max(-1d0,cthe(3-jt)))
          phir=paru(2)*pjr(0)
          cphi=cos(phir)
          ang=cthe(1)*cthe(2)-sqrt(1d0-cthe(1)**2)*
     &    sqrt(1d0-cthe(2)**2)*cphi
          z1=2d0-z(jt)
          z2=ang*sqrt(z(jt)**2-4d0*pmq(jt)**2/shp)
          z3=1d0-z(jt)-xh+(pmq(1)**2+pmq(2)**2)/shp
          z(3-jt)=2d0/(z1**2-z2**2)*(z1*z3+z2*sqrt(z3**2-(z1**2-z2**2)*
     &    pmq(3-jt)**2/shp))
          zmin=2d0*pmq(3-jt)/shpr
          zmax=1d0-pmq(jt)/shpr-(sh-pmq(3-jt)**2)/(shpr*(shpr-pmq(jt)))
          zmax=min(1d0-xh,zmax)
          if(z(3-jt).lt.zmin.or.z(3-jt).gt.zmax) goto 370
          kcc=22
 
        elseif(isub.eq.78) then
C...W+/- + h0 -> W+/- + h0
 
        elseif(isub.eq.79) then
C...h0 + h0 -> h0 + h0
 
        elseif(isub.eq.80) then
C...q + gamma -> q' + pi+/-; th=(p(q)-p(q'))**2
          if(mint(15).eq.22) js=2
          i=mint(14+js)
          ia=iabs(i)
          mint(23-js)=isign(211,kchg(ia,1)*i)
          ib=3-ia
          mint(20+js)=isign(ib,i)
          kcc=22
        endif
 
      elseif(isub.le.90) then
        if(isub.eq.81) then
C...q + qbar -> Q + Qbar; th = (p(q)-p(Q))**2
          mint(21)=isign(mint(55),mint(15))
          mint(22)=-mint(21)
          kcc=4
 
        elseif(isub.eq.82) then
C...g + g -> Q + Qbar; th arbitrary
          kcs=(-1)**int(1.5d0+pjr(0))
          mint(21)=isign(mint(55),kcs)
          mint(22)=-mint(21)
          kcc=mint(2)+10
 
        elseif(isub.eq.83) then
C...f + q -> f' + Q; th = (p(f) - p(f'))**2
          kfold=mint(16)
          if(mint(2).eq.2) kfold=mint(15)
          kfaold=iabs(kfold)
          if(kfaold.gt.10) then
            kfanew=kfaold+2*mod(kfaold,2)-1
          else
            rckm=vint(180+kfold)*pjr(0)
            ipm=(5-isign(1,kfold))/2
            kfanew=-mod(kfaold+1,2)
  410       kfanew=kfanew+2
            idc=mdcy(kfaold,2)+(kfanew+1)/2+2
            if(mdme(idc,1).eq.1.or.mdme(idc,1).eq.ipm) then
              if(mod(kfaold,2).eq.0) rckm=rckm-
     &        vckm(kfaold/2,(kfanew+1)/2)
              if(mod(kfaold,2).eq.1) rckm=rckm-
     &        vckm(kfanew/2,(kfaold+1)/2)
            endif
            if(kfanew.le.6.and.rckm.gt.0d0) goto 410
          endif
          if(mint(2).eq.1) then
            mint(21)=isign(mint(55),mint(15))
            mint(22)=isign(kfanew,mint(16))
          else
            mint(21)=isign(kfanew,mint(15))
            mint(22)=isign(mint(55),mint(16))
            js=2
          endif
          kcc=22
 
        elseif(isub.eq.84) then
C...g + gamma -> Q + Qbar; th arbitary
          kcs=(-1)**int(1.5d0+pjr(0))
          mint(21)=isign(mint(55),kcs)
          mint(22)=-mint(21)
          kcc=27
          if(mint(16).eq.21) kcc=28
 
        elseif(isub.eq.85) then
C...gamma + gamma -> F + Fbar; th arbitary
          kcs=(-1)**int(1.5d0+pjr(0))
          mint(21)=isign(mint(56),kcs)
          mint(22)=-mint(21)
          kcc=21
 
        elseif(isub.ge.86.and.isub.le.89) then
C...g + g -> (J/Psi, chi_0c, chi_1c or chi_2c) + g
          mint(21)=kfpr(isub,1)
          mint(22)=kfpr(isub,2)
          kcc=24
          kcs=(-1)**int(1.5d0+pjr(0))
        endif
 
      elseif(isub.le.100) then
        if(isub.eq.95) then
C...Low-pT ( = energyless g + g -> g + g)
          kcc=mint(2)+12
          kcs=(-1)**int(1.5d0+pjr(0))
 
        elseif(isub.eq.96) then
C...Multiple interactions (should be reassigned to QCD process)
        endif
 
      elseif(isub.le.110) then
        if(isub.eq.101) then
C...g + g -> gamma*/Z0
          kcc=21
          kfres=22
 
        elseif(isub.eq.102) then
C...g + g -> h0 (or H0, or A0)
          kcc=21
          kfres=kfhigg
 
        elseif(isub.eq.103) then
C...gamma + gamma -> h0 (or H0, or A0)
          kcc=21
          kfres=kfhigg
 
        elseif(isub.eq.106) then
C...g + g -> J/Psi + gamma
          mint(21)=kfpr(isub,1)
          mint(22)=kfpr(isub,2)
          kcc=21
 
        elseif(isub.eq.107) then
C...g + gamma -> J/Psi + g
          mint(21)=kfpr(isub,1)
          mint(22)=kfpr(isub,2)
          kcc=22
          if(mint(16).eq.22) kcc=33
 
        elseif(isub.eq.108) then
C...gamma + gamma -> J/Psi + gamma
          mint(21)=kfpr(isub,1)
          mint(22)=kfpr(isub,2)
 
        elseif(isub.eq.110) then
C...f + fbar -> gamma + h0; th arbitrary
          if(pjr(0).gt.0.5d0) js=2
          mint(20+js)=22
          mint(23-js)=kfhigg
        endif
 
      elseif(isub.le.120) then
        if(isub.eq.111) then
C...f + fbar -> g + h0; th arbitrary
          if(pjr(0).gt.0.5d0) js=2
          mint(20+js)=21
          mint(23-js)=25
          kcc=17+js
 
        elseif(isub.eq.112) then
C...f + g -> f + h0; th = (p(f) - p(f))**2
          if(mint(15).eq.21) js=2
          mint(23-js)=25
          kcc=15+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.113) then
C...g + g -> g + h0; th arbitrary
          if(pjr(0).gt.0.5d0) js=2
          mint(23-js)=25
          kcc=22+js
          kcs=(-1)**int(1.5d0+pjr(0))
 
        elseif(isub.eq.114) then
C...g + g -> gamma + gamma; th arbitrary
          if(pjr(0).gt.0.5d0) js=2
          mint(21)=22
          mint(22)=22
          kcc=21
 
        elseif(isub.eq.115) then
C...g + g -> g + gamma; th arbitrary
          if(pjr(0).gt.0.5d0) js=2
          mint(23-js)=22
          kcc=22+js
          kcs=(-1)**int(1.5d0+pjr(0))
 
        elseif(isub.eq.116) then
C...g + g -> gamma + Z0
 
        elseif(isub.eq.117) then
C...g + g -> Z0 + Z0
 
        elseif(isub.eq.118) then
C...g + g -> W+ + W-
        endif
 
      elseif(isub.le.140) then
        if(isub.eq.121) then
C...g + g -> Q + Qbar + h0
          kcs=(-1)**int(1.5d0+pjr(0))
          mint(21)=isign(kfpr(isubsv,2),kcs)
          mint(22)=-mint(21)
          kcc=11+int(0.5d0+pjr(0))
          kfres=kfhigg
 
        elseif(isub.eq.122) then
C...q + qbar -> Q + Qbar + h0
          mint(21)=isign(kfpr(isubsv,2),mint(15))
          mint(22)=-mint(21)
          kcc=4
          kfres=kfhigg
 
        elseif(isub.eq.123) then
C...f + f' -> f + f' + h0 (or H0, or A0) (Z0 + Z0 -> h0 as
C...inner process)
          kcc=22
          kfres=kfhigg
 
        elseif(isub.eq.124) then
C...f + f' -> f" + f"' + h0 (or H0, or A) (W+ + W- -> h0 as
C...inner process)
          do 430 jt=1,2
            i=mint(14+jt)
            ia=iabs(i)
            if(ia.le.10) then
              rvckm=vint(180+i)*pjr(0)
              do 420 j=1,mstp(1)
                ib=2*j-1+mod(ia,2)
                ipm=(5-isign(1,i))/2
                idc=j+mdcy(ia,2)+2
                if(mdme(idc,1).ne.1.and.mdme(idc,1).ne.ipm) goto 420
                mint(20+jt)=isign(ib,i)
                rvckm=rvckm-vckm((ia+1)/2,(ib+1)/2)
                if(rvckm.le.0d0) goto 430
  420         continue
            else
              ib=2*((ia+1)/2)-1+mod(ia,2)
              mint(20+jt)=isign(ib,i)
            endif
  430     continue
          kcc=22
          kfres=kfhigg
 
        elseif(isub.eq.131) then
C...g + g -> Z0 + q + qbar
        endif
 
      elseif(isub.le.160) then
        if(isub.eq.141) then
C...f + fbar -> gamma*/Z0/Z'0
          kfres=32
 
        elseif(isub.eq.142) then
C...f + fbar' -> W'+/-
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          kfres=isign(34,kch1+kch2)
 
        elseif(isub.eq.143) then
C...f + fbar' -> H+/-
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          kfres=isign(37,kch1+kch2)
 
        elseif(isub.eq.144) then
C...f + fbar' -> R
          kfres=isign(40,mint(15)+mint(16))
 
        elseif(isub.eq.145) then
C...q + l -> LQ (leptoquark)
          if(iabs(mint(16)).le.8) js=2
          kfres=isign(39,mint(14+js))
          kcc=28+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.147.or.isub.eq.148) then
C...q + g -> q* (excited quark)
          if(mint(15).eq.21) js=2
          kfres=isign(kfpr(isub,1),mint(14+js))
          kcc=30+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.149) then
C...g + g -> eta_techni
          kfres=38
          kcc=23
          kcs=(-1)**int(1.5d0+pjr(0))
        endif
 
      elseif(isub.le.200) then
        if(isub.eq.161) then
C...f + g -> f' + H+/-; th = (p(f)-p(f'))**2
          if(mint(15).eq.21) js=2
          i=mint(14+js)
          ia=iabs(i)
          mint(23-js)=isign(37,kchg(ia,1)*i)
          ib=ia+mod(ia,2)-mod(ia+1,2)
          mint(20+js)=isign(ib,i)
          kcc=15+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.162) then
C...q + g -> LQ + lbar; LQ=leptoquark; th=(p(q)-p(LQ))^2
          if(mint(15).eq.21) js=2
          mint(20+js)=isign(39,mint(14+js))
          kflql=kfdp(mdcy(39,2),2)
          mint(23-js)=-isign(kflql,mint(14+js))
          kcc=15+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.163) then
C...g + g -> LQ + LQbar; LQ=leptoquark; th arbitrary
          kcs=(-1)**int(1.5d0+pjr(0))
          mint(21)=isign(39,kcs)
          mint(22)=-mint(21)
          kcc=mint(2)+10
 
        elseif(isub.eq.164) then
C...q + qbar -> LQ + LQbar; LQ=leptoquark; th=(p(q)-p(LQ))**2
          mint(21)=isign(39,mint(15))
          mint(22)=-mint(21)
          kcc=4
 
        elseif(isub.eq.165) then
C...q + qbar -> l- + l+; th=(p(q)-p(l-))**2
          mint(21)=isign(kfpr(isub,1),mint(15))
          mint(22)=-mint(21)
 
        elseif(isub.eq.166) then
C...q + qbar' -> l + nu; th=(p(u)-p(nu))**2 or (p(ubar)-p(nubar))**2
          if(mod(mint(15),2).eq.0) then
            mint(21)=isign(kfpr(isub,1)+1,mint(15))
            mint(22)=isign(kfpr(isub,1),mint(16))
          else
            mint(21)=isign(kfpr(isub,1),mint(15))
            mint(22)=isign(kfpr(isub,1)+1,mint(16))
          endif
 
        elseif(isub.eq.167.or.isub.eq.168) then
C...q + q' -> q" + q* (excited quark)
          kfqstr=kfpr(isub,2)
          kfqexc=mod(kfqstr,kexcit)
          js=mint(2)
          mint(20+js)=isign(kfqstr,mint(14+js))
          if(iabs(mint(15)).ne.kfqexc.and.iabs(mint(16)).ne.kfqexc)
     &    mint(23-js)=isign(kfqexc,mint(17-js))
          kcc=22
 
        elseif(isub.eq.191) then
C...f + fbar -> rho_tech0.
          kfres=54
 
        elseif(isub.eq.192) then
C...f + fbar' -> rho_tech+/-
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          kfres=isign(55,kch1+kch2)
 
        elseif(isub.eq.193) then
C...f + fbar -> omega_tech0.
          kfres=56
 
        elseif(isub.eq.194) then
C...f + fbar -> f' + fbar' via mixture of s-channel
C...rho_tech and omega_tech; th=(p(f)-p(f'))**2
          mint(21)=isign(kfpr(isub,1),mint(15))
          mint(22)=-mint(21)
         endif
 
CMRENNA++
      elseif(isub.le.215) then
        if(isub.eq.201) then
C...f + fbar -> ~e_L + ~e_Lbar
          mint(21)=isign(ksusy1+11,kcs)
          mint(22)=-mint(21)
 
        elseif(isub.eq.202) then
C...f + fbar -> ~e_R + ~e_Rbar
          mint(21)=isign(ksusy2+11,kcs)
          mint(22)=-mint(21)
 
        elseif(isub.eq.203) then
C...f + fbar -> ~e_R + ~e_Lbar
          kcs=1
          if(mint(2).eq.2) kcs=-1
          mint(21)=isign(ksusy1+11,kcs)
          mint(22)=-isign(ksusy2+11,kcs)
 
        elseif(isub.eq.204) then
C...f + fbar -> ~mu_L + ~mu_Lbar
          mint(21)=isign(ksusy1+13,kcs)
          mint(22)=-mint(21)
 
        elseif(isub.eq.205) then
C...f + fbar -> ~mu_R + ~mu_Rbar
          mint(21)=isign(ksusy2+13,kcs)
          mint(22)=-mint(21)
 
        elseif(isub.eq.206) then
C...f + fbar -> ~mu_L + ~mu_Rbar
          kcs=1
          if(mint(2).eq.2) kcs=-1
          mint(21)=isign(ksusy1+13,kcs)
          mint(22)=-isign(ksusy2+13,kcs)
 
        elseif(isub.eq.207) then
C...f + fbar -> ~tau_1 + ~tau_1bar
          mint(21)=isign(ksusy1+15,kcs)
          mint(22)=-mint(21)
 
        elseif(isub.eq.208) then
C...f + fbar -> ~tau_2 + ~tau_2bar
          mint(21)=isign(ksusy2+15,kcs)
          mint(22)=-mint(21)
 
        elseif(isub.eq.209) then
C...f + fbar -> ~tau_1 + ~tau_2bar
          kcs=1
          if(mint(2).eq.2) kcs=-1
          mint(21)=isign(ksusy1+15,kcs)
          mint(22)=-isign(ksusy2+15,kcs)
 
        elseif(isub.eq.210) then
C...q + qbar' -> ~l_L + ~nulbar; th arbitrary
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          mint(21)=-isign(kfpr(isub,1),kch1+kch2)
          mint(22)=isign(kfpr(isub,2),kch1+kch2)
 
        elseif(isub.eq.211) then
C...q + qbar'-> ~tau_1 + ~nutaubar; th arbitrary
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          mint(21)=-isign(ksusy1+15,kch1+kch2)
          mint(22)=isign(ksusy1+16,kch1+kch2)
 
        elseif(isub.eq.212) then
C...q + qbar'-> ~tau_2 + ~nutaubar; th arbitrary
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          mint(21)=-isign(ksusy2+15,kch1+kch2)
          mint(22)=isign(ksusy1+16,kch1+kch2)
 
        elseif(isub.eq.213) then
C...f + fbar -> ~nul + ~nulbar
          mint(21)=isign(kfpr(isub,1),kcs)
          mint(22)=-mint(21)
 
        elseif(isub.eq.214) then
C...f + fbar -> ~nutau + ~nutaubar
          mint(21)=isign(ksusy1+16,kcs)
          mint(22)=-mint(21)
        endif
 
      elseif(isub.le.225) then
        if(isub.eq.216) then
C...f + fbar -> ~chi01 + ~chi01
          mint(21)=ksusy1+22
          mint(22)=ksusy1+22
 
        elseif(isub.eq.217) then
C...f + fbar -> ~chi02 + ~chi02
          mint(21)=ksusy1+23
          mint(22)=ksusy1+23
 
        elseif(isub.eq.218 ) then
C...f + fbar -> ~chi03 + ~chi03
          mint(21)=ksusy1+25
          mint(22)=ksusy1+25
 
        elseif(isub.eq.219 ) then
C...f + fbar -> ~chi04 + ~chi04
          mint(21)=ksusy1+35
          mint(22)=ksusy1+35
 
        elseif(isub.eq.220 ) then
C...f + fbar -> ~chi01 + ~chi02
          if(pjr(0).gt.0.5d0) js=2
          mint(20+js)=ksusy1+22
          mint(23-js)=ksusy1+23
 
        elseif(isub.eq.221 ) then
C...f + fbar -> ~chi01 + ~chi03
          if(pjr(0).gt.0.5d0) js=2
          mint(20+js)=ksusy1+22
          mint(23-js)=ksusy1+25
 
        elseif(isub.eq.222) then
C...f + fbar -> ~chi01 + ~chi04
          if(pjr(0).gt.0.5d0) js=2
          mint(20+js)=ksusy1+22
          mint(23-js)=ksusy1+35
 
        elseif(isub.eq.223) then
C...f + fbar -> ~chi02 + ~chi03
          if(pjr(0).gt.0.5d0) js=2
          mint(20+js)=ksusy1+23
          mint(23-js)=ksusy1+25
 
        elseif(isub.eq.224) then
C...f + fbar -> ~chi02 + ~chi04
          if(pjr(0).gt.0.5d0) js=2
          mint(20+js)=ksusy1+23
          mint(23-js)=ksusy1+35
 
        elseif(isub.eq.225) then
C...f + fbar -> ~chi03 + ~chi04
          if(pjr(0).gt.0.5d0) js=2
          mint(20+js)=ksusy1+25
          mint(23-js)=ksusy1+35
        endif
 
      elseif(isub.le.236) then
        if(isub.eq.226) then
C...f + fbar -> ~chi+-1 + ~chi-+1
C...th=(p(q)-p(chi+))**2 or (p(qbar)-p(chi-))**2
          mint(21)=isign(ksusy1+24,mint(15))
          mint(22)=-mint(21)
 
        elseif(isub.eq.227) then
C...f + fbar -> ~chi+-2 + ~chi-+2
          mint(21)=isign(ksusy1+37,mint(15))
          mint(22)=-mint(21)
 
        elseif(isub.eq.228) then
C...f + fbar -> ~chi+-1 + ~chi-+2
C...th=(p(q)-p(chi1+))**2 or th=(p(qbar)-p(chi1-))**2
C...js=1 if pjr<.5, js=2 if pjr>.5
C...if 15=q, 16=qbar and js=1, chi1+ + chi2-, th=(q-chi1+)**2
C...if 15=qbar, 16=q and js=1, chi2- + chi1+, th=(q-chi1+)**2
C...if 15=q, 16=qbar and js=2, chi1- + chi2+, th=(qbar-chi1-)**2
C...if 15=qbar, 16=q and js=2, chi2+ + chi1-, th=(q-chi1-)**2
          kch1=isign(1,mint(15))
          kch2=int(1-kch1)/2
          if(mint(2).eq.1) then
            mint(22-kch2)= -(ksusy1+24)
            mint(21+kch2)= ksusy1+37
            if(kch2.eq.0) js=2
          else
            mint(21+kch2)= ksusy1+24
            mint(22-kch2)= -(ksusy1+37)
            if(kch2.eq.1) js=2
          endif
 
        elseif(isub.eq.229) then
C...q + qbar' -> ~chi01 + ~chi+-1
C...th=(p(u)-p(chi+))**2 or (p(ubar)-p(chi-))**2
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
C...CHECK THIS
          if(mod(mint(15),2).ne.0) js=2
          mint(20+js)=ksusy1+22
          mint(23-js)=isign(ksusy1+24,kch1+kch2)
 
        elseif(isub.eq.230) then
C...q + qbar' -> ~chi02 + ~chi+-1
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          if(mod(mint(15),2).ne.0) js=2
          mint(20+js)=ksusy1+23
          mint(23-js)=isign(ksusy1+24,kch1+kch2)
 
        elseif(isub.eq.231) then
C...q + qbar' -> ~chi03 + ~chi+-1
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          if(mod(mint(15),2).ne.0) js=2
          mint(20+js)=ksusy1+25
          mint(23-js)=isign(ksusy1+24,kch1+kch2)
 
        elseif(isub.eq.232) then
C...q + qbar' -> ~chi04 + ~chi+-1
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          if(mod(mint(15),2).ne.0) js=2
          mint(20+js)=ksusy1+35
          mint(23-js)=isign(ksusy1+24,kch1+kch2)
 
        elseif(isub.eq.233) then
C...q + qbar' -> ~chi01 + ~chi+-2
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          if(mod(mint(15),2).ne.0) js=2
          mint(20+js)=ksusy1+22
          mint(23-js)=isign(ksusy1+37,kch1+kch2)
 
        elseif(isub.eq.234) then
C...q + qbar' -> ~chi02 + ~chi+-2
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          if(mod(mint(15),2).ne.0) js=2
          mint(20+js)=ksusy1+23
          mint(23-js)=isign(ksusy1+37,kch1+kch2)
 
        elseif(isub.eq.235) then
C...q + qbar' -> ~chi03 + ~chi+-2
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          if(mod(mint(15),2).ne.0) js=2
          mint(20+js)=ksusy1+25
          mint(23-js)=isign(ksusy1+37,kch1+kch2)
 
        elseif(isub.eq.236) then
C...q + qbar' -> ~chi04 + ~chi+-2
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          if(mod(mint(15),2).ne.0) js=2
          mint(20+js)=ksusy1+35
          mint(23-js)=isign(ksusy1+37,kch1+kch2)
        endif
 
      elseif(isub.le.245) then
        if(isub.eq.237) then
C...q + qbar -> ~chi01 + ~g
C...th arbitrary
          if(pjr(0).gt.0.5d0) js=2
          mint(20+js)=ksusy1+21
          mint(23-js)=ksusy1+22
          kcc=17+js
 
        elseif(isub.eq.238) then
C...q + qbar -> ~chi02 + ~g
C...th arbitrary
          if(pjr(0).gt.0.5d0) js=2
          mint(20+js)=ksusy1+21
          mint(23-js)=ksusy1+23
          kcc=17+js
 
        elseif(isub.eq.239) then
C...q + qbar -> ~chi03 + ~g
C...th arbitrary
          if(pjr(0).gt.0.5d0) js=2
          mint(20+js)=ksusy1+21
          mint(23-js)=ksusy1+25
          kcc=17+js
 
        elseif(isub.eq.240) then
C...q + qbar -> ~chi04 + ~g
C...th arbitrary
          if(pjr(0).gt.0.5d0) js=2
          mint(20+js)=ksusy1+21
          mint(23-js)=ksusy1+35
          kcc=17+js
 
        elseif(isub.eq.241) then
C...q + qbar' -> ~chi+-1 + ~g
C...if 15=u, 16=dbar, then (kch1+kch2)>0, js=1, chi+
C...if 15=d, 16=ubar, then (kch1+kch2)<0, js=2, chi-
C...if 15=ubar, 16=d, then (kch1+kch2)<0, js=1, chi-
C...if 15=dbar, 16=u, then (kch1+kch2)>0, js=2, chi+
C...th=(p(q)-p(chi+))**2 or (p(qbar')-p(chi-))**2
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          js=1
          if(mint(15)*(kch1+kch2).gt.0) js=2
          mint(20+js)=ksusy1+21
          mint(23-js)=isign(ksusy1+24,kch1+kch2)
          kcc=17+js
 
        elseif(isub.eq.242) then
C...q + qbar' -> ~chi+-2 + ~g
C...if 15=u, 16=dbar, then (kch1+kch2)>0, js=1, chi+
C...if 15=d, 16=ubar, then (kch1+kch2)<0, js=2, chi-
C...if 15=ubar, 16=d, then (kch1+kch2)<0, js=1, chi-
C...if 15=dbar, 16=u, then (kch1+kch2)>0, js=2, chi+
C...th=(p(q)-p(chi+))**2 or (p(qbar')-p(chi-))**2
          kch1=kchg(iabs(mint(15)),1)*isign(1,mint(15))
          kch2=kchg(iabs(mint(16)),1)*isign(1,mint(16))
          js=1
          if(mint(15)*(kch1+kch2).gt.0) js=2
          mint(20+js)=ksusy1+21
          mint(23-js)=isign(ksusy1+37,kch1+kch2)
          kcc=17+js
 
        elseif(isub.eq.243) then
C...q + qbar -> ~g + ~g ; th arbitrary
          mint(21)=ksusy1+21
          mint(22)=ksusy1+21
          kcc=mint(2)+4
 
        elseif(isub.eq.244) then
C...g + g -> ~g + ~g ; th arbitrary
          kcc=mint(2)+12
          kcs=(-1)**int(1.5d0+pjr(0))
          mint(21)=ksusy1+21
          mint(22)=ksusy1+21
        endif
 
      elseif(isub.le.260) then
        if(isub.eq.246) then
C...qj + g -> ~qj_L + ~chi01
          if(mint(15).eq.21) js=2
          i=mint(14+js)
          ia=iabs(i)
          mint(20+js)=isign(ksusy1+ia,i)
          mint(23-js)=ksusy1+22
          kcc=15+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.247) then
C...qj + g -> ~qj_R + ~chi01
          if(mint(15).eq.21) js=2
          i=mint(14+js)
          ia=iabs(i)
          mint(20+js)=isign(ksusy2+ia,i)
          mint(23-js)=ksusy1+22
          kcc=15+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.248) then
C...qj + g -> ~qj_L + ~chi02
          if(mint(15).eq.21) js=2
          i=mint(14+js)
          ia=iabs(i)
          mint(20+js)=isign(ksusy1+ia,i)
          mint(23-js)=ksusy1+23
          kcc=15+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.249) then
C...qj + g -> ~qj_R + ~chi02
          if(mint(15).eq.21) js=2
          i=mint(14+js)
          ia=iabs(i)
          mint(20+js)=isign(ksusy2+ia,i)
          mint(23-js)=ksusy1+23
          kcc=15+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.250) then
C...qj + g -> ~qj_L + ~chi03
          if(mint(15).eq.21) js=2
          i=mint(14+js)
          ia=iabs(i)
          mint(20+js)=isign(ksusy1+ia,i)
          mint(23-js)=ksusy1+25
          kcc=15+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.251) then
C...qj + g -> ~qj_R + ~chi03
          if(mint(15).eq.21) js=2
          i=mint(14+js)
          ia=iabs(i)
          mint(20+js)=isign(ksusy2+ia,i)
          mint(23-js)=ksusy1+25
          kcc=15+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.252) then
C...qj + g -> ~qj_L + ~chi04
          if(mint(15).eq.21) js=2
          i=mint(14+js)
          ia=iabs(i)
          mint(20+js)=isign(ksusy1+ia,i)
          mint(23-js)=ksusy1+35
          kcc=15+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.253) then
C...qj + g -> ~qj_R + ~chi04
          if(mint(15).eq.21) js=2
          i=mint(14+js)
          ia=iabs(i)
          mint(20+js)=isign(ksusy2+ia,i)
          mint(23-js)=ksusy1+35
          kcc=15+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.254) then
C...qj + g -> ~qk_L + ~chi+-1
          if(mint(15).eq.21) js=2
          i=mint(14+js)
          ia=iabs(i)
          mint(23-js)=isign(ksusy1+24,kchg(ia,1)*i)
          ib=-ia+int((ia+1)/2)*4-1
          mint(20+js)=isign(ksusy1+ib,i)
          kcc=15+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.255) then
C...qj + g -> ~qk_L + ~chi+-1
          if(mint(15).eq.21) js=2
          i=mint(14+js)
          ia=iabs(i)
          mint(23-js)=isign(ksusy1+24,kchg(ia,1)*i)
          ib=-ia+int((ia+1)/2)*4-1
          mint(20+js)=isign(ksusy2+ib,i)
          kcc=15+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.256) then
C...qj + g -> ~qk_L + ~chi+-2
          if(mint(15).eq.21) js=2
          i=mint(14+js)
          ia=iabs(i)
          ib=-ia+int((ia+1)/2)*4-1
          mint(20+js)=isign(ksusy1+ib,i)
          mint(23-js)=isign(ksusy1+37,kchg(ia,1)*i)
          kcc=15+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.257) then
C...qj + g -> ~qk_R + ~chi+-2
          if(mint(15).eq.21) js=2
          i=mint(14+js)
          ia=iabs(i)
          ib=-ia+int((ia+1)/2)*4-1
          mint(20+js)=isign(ksusy2+ib,i)
          mint(23-js)=isign(ksusy1+37,kchg(ia,1)*i)
          kcc=15+js
          kcs=isign(1,mint(14+js))
 
        elseif(isub.eq.258) then
C...qj + g -> ~qj_L + ~g
          if(mint(15).eq.21) js=2
          i=mint(14+js)
          ia=iabs(i)
          mint(20+js)=isign(ksusy1+ia,i)
          mint(23-js)=ksusy1+21
          kcc=mint(2)+6
          if(js.eq.2) kcc=kcc+2
          kcs=isign(1,i)
 
        elseif(isub.eq.259) then
C...qj + g -> ~qj_R + ~g
          if(mint(15).eq.21) js=2
          i=mint(14+js)
          ia=iabs(i)
          mint(20+js)=isign(ksusy2+ia,i)
          mint(23-js)=ksusy1+21
          kcc=mint(2)+6
          if(js.eq.2) kcc=kcc+2
          kcs=isign(1,i)
        endif
 
      elseif(isub.le.270) then
        if(isub.eq.261) then
C...f + fbar -> ~t_1 + ~t_1bar; th = (p(q)-p(sq))**2
          mint(21)=isign(kfpr(isub,1),kcs)
          mint(22)=-mint(21)
C...Correct color combination
          if(mint(43).eq.4) kcc=4
 
        elseif(isub.eq.262) then
C...f + fbar -> ~t_2 + ~t_2bar; th = (p(q)-p(sq))**2
          mint(21)=isign(kfpr(isub,1),kcs)
          mint(22)=-mint(21)
C...Correct color combination
          if(mint(43).eq.4) kcc=4
 
        elseif(isub.eq.263) then
C...f + fbar -> ~t_1 + ~t_2bar; th = (p(q)-p(sq))**2
          if((kcs.gt.0.and.mint(2).eq.1).or.
     &    (kcs.lt.0.and.mint(2).eq.2)) then
            mint(21)=isign(kfpr(isub,1),kcs)
            mint(22)=-isign(kfpr(isub,2),kcs)
          else
            js=2
            mint(21)=isign(kfpr(isub,2),kcs)
            mint(22)=-isign(kfpr(isub,1),kcs)
          endif
C...Correct color combination
          if(mint(43).eq.4) kcc=4
 
        elseif(isub.eq.264) then
C...g + g -> ~t_1 + ~t_1bar; th arbitrary
          kcs=(-1)**int(1.5d0+pjr(0))
          mint(21)=isign(kfpr(isub,1),kcs)
          mint(22)=-mint(21)
          kcc=mint(2)+10
 
        elseif(isub.eq.265) then
C...g + g -> ~t_2 + ~t_2bar; th arbitrary
          kcs=(-1)**int(1.5d0+pjr(0))
          mint(21)=isign(kfpr(isub,1),kcs)
          mint(22)=-mint(21)
          kcc=mint(2)+10
        endif
 
      elseif(isub.le.280) then
        if(isub.eq.271) then
C...qi + qj -> ~qi_L + ~qj_L
          kcc=mint(2)
          if(mint(15)*mint(16).lt.0) kcc=kcc+2
          mint(21)=isign(ksusy1+iabs(mint(15)),mint(15))
          mint(22)=isign(ksusy1+iabs(mint(16)),mint(16))
 
        elseif(isub.eq.272) then
C...qi + qj -> ~qi_R + ~qj_R
          kcc=mint(2)
          if(mint(15)*mint(16).lt.0) kcc=kcc+2
          mint(21)=isign(ksusy2+iabs(mint(15)),mint(15))
          mint(22)=isign(ksusy2+iabs(mint(16)),mint(16))
 
        elseif(isub.eq.273) then
C...qi + qj -> ~qi_L + ~qj_R
          mint(21)=isign(kfpr(isub,1),mint(15))
          mint(22)=isign(kfpr(isub,2),mint(16))
          kcc=mint(2)
          if(mint(15)*mint(16).lt.0) kcc=kcc+2
 
        elseif(isub.eq.274) then
C...qi + qjbar -> ~qi_L + ~qj_Lbar; th = (p(f)-p(sf'))**2
          mint(21)=isign(ksusy1+iabs(mint(15)),mint(15))
          mint(22)=isign(ksusy1+iabs(mint(16)),mint(16))
          kcc=mint(2)
          if(mint(15)*mint(16).lt.0) kcc=kcc+2
 
        elseif(isub.eq.275) then
C...qi + qjbar -> ~qi_R + ~qj_Rbar ; th = (p(f)-p(sf'))**2
          mint(21)=isign(ksusy2+iabs(mint(15)),mint(15))
          mint(22)=isign(ksusy2+iabs(mint(16)),mint(16))
          kcc=mint(2)
          if(mint(15)*mint(16).lt.0) kcc=kcc+2
 
        elseif(isub.eq.276) then
C...qi + qjbar -> ~qi_L + ~qj_Rbar ; th = (p(f)-p(sf'))**2
          mint(21)=isign(kfpr(isub,1),mint(15))
          mint(22)=isign(kfpr(isub,2),mint(16))
          kcc=mint(2)
          if(mint(15)*mint(16).lt.0) kcc=kcc+2
 
        elseif(isub.eq.277) then
C...f + fbar -> ~qi_L + ~qi_Lbar ; th = (p(q)-p(sq))**2
          isgn=1
          if(mint(43).eq.1.and.pjr(0).gt.0.5d0) isgn=-1
          mint(21)=isgn*isign(kfpr(isub,1),kcs)
          mint(22)=-mint(21)
          if(mint(43).eq.4) kcc=4
 
        elseif(isub.eq.278) then
C...f + fbar -> ~qi_R + ~qi_Rbar; th = (p(q)-p(sq))**2
          isgn=1
          if(mint(43).eq.1.and.pjr(0).gt.0.5d0) isgn=-1
          mint(21)=isgn*isign(kfpr(isub,1),kcs)
          mint(22)=-mint(21)
          if(mint(43).eq.4) kcc=4
 
        elseif(isub.eq.279) then
C...g + g -> ~qi_L + ~qi_Lbar ; th arbitrary
C...pure LL + RR
          kcs=(-1)**int(1.5d0+pjr(0))
          mint(21)=isign(kfpr(isub,1),kcs)
          mint(22)=-mint(21)
          kcc=mint(2)+10
 
        elseif(isub.eq.280) then
C...g + g -> ~qi_R + ~qi_Rbar ; th arbitrary
          kcs=(-1)**int(1.5d0+pjr(0))
          mint(21)=isign(kfpr(isub,1),kcs)
          mint(22)=-mint(21)
          kcc=mint(2)+10
        endif
 
CMRENNA--
      endif
 
      if(iset(isub).eq.11) then
C...Store documentation for user-defined processes
        bezup=(pup(1,4)-pup(2,4))/(pup(1,4)+pup(2,4))
        kuppo(1)=mint(83)+5
        kuppo(2)=mint(83)+6
        i=mint(83)+6
        do 450 iup=3,nup
          kuppo(iup)=0
          if(mstp(128).ge.2.and.kup(iup,3).ne.0) then
            idoc=idoc-1
            mint(4)=mint(4)-1
            goto 450
          endif
          i=i+1
          kuppo(iup)=i
          k(i,1)=21
          k(i,2)=kup(iup,2)
          k(i,3)=0
          if(kup(iup,3).ne.0) k(i,3)=kuppo(kup(iup,3))
          k(i,4)=0
          k(i,5)=0
          do 440 j=1,5
            p(i,j)=pup(iup,j)
  440     continue
  450   continue
        call pjrobo(mint(83)+7,mint(83)+4+nup,0d0,vint(24),0d0,0d0,
     &  -bezup)
 
C...Store final state partons for user-defined processes
        n=ipu2
        do 470 iup=3,nup
          n=n+1
          k(n,1)=1
          if(kup(iup,1).ne.1) k(n,1)=11
          k(n,2)=kup(iup,2)
          if(mstp(128).le.0.or.kup(iup,3).eq.0) then
            k(n,3)=kuppo(iup)
          else
            k(n,3)=mint(84)+kup(iup,3)
          endif
          k(n,4)=0
          k(n,5)=0
          do 460 j=1,5
            p(n,j)=pup(iup,j)
  460     continue
  470   continue
        call pjrobo(ipu3,n,0d0,vint(24),0d0,0d0,-bezup)
 
C...Arrange colour flow for user-defined processes
        n=mint(84)
        do 480 iup=1,nup
          n=n+1
          if(kchg(jamcomp(k(n,2)),2).eq.0) goto 480
          if(k(n,1).eq.1) k(n,1)=3
          if(k(n,1).eq.11) k(n,1)=14
          if(kup(iup,4).ne.0) k(n,4)=k(n,4)+mstu(5)*(kup(iup,4)+
     &    mint(84))
          if(kup(iup,5).ne.0) k(n,5)=k(n,5)+mstu(5)*(kup(iup,5)+
     &    mint(84))
          if(kup(iup,6).ne.0) k(n,4)=k(n,4)+kup(iup,6)+mint(84)
          if(kup(iup,7).ne.0) k(n,5)=k(n,5)+kup(iup,7)+mint(84)
  480   continue
 
      elseif(idoc.eq.7) then
C...Resonance not decaying; store kinematics
        i=mint(83)+7
        k(ipu3,1)=1
        k(ipu3,2)=kfres
        k(ipu3,3)=i
        p(ipu3,4)=shuser
        p(ipu3,5)=shuser
        k(i,1)=21
        k(i,2)=kfres
        p(i,4)=shuser
        p(i,5)=shuser
        n=ipu3
        mint(21)=kfres
        mint(22)=0
 
C...Special cases: colour flow in coloured resonances
        kcres=jamcomp(kfres)
        if(kchg(kcres,2).ne.0) then
          k(ipu3,1)=3
          do 490 j=1,2
            jc=j
            if(kcs.eq.-1) jc=3-j
            if(icol(kcc,1,jc).ne.0.and.k(ipu1,1).eq.14) k(ipu1,j+3)=
     &      mint(84)+icol(kcc,1,jc)
            if(icol(kcc,2,jc).ne.0.and.k(ipu2,1).eq.14) k(ipu2,j+3)=
     &      mint(84)+icol(kcc,2,jc)
            if(icol(kcc,3,jc).ne.0.and.k(ipu3,1).eq.3) k(ipu3,j+3)=
     &      mstu(5)*(mint(84)+icol(kcc,3,jc))
  490     continue
        else
          k(ipu1,4)=ipu2
          k(ipu1,5)=ipu2
          k(ipu2,4)=ipu1
          k(ipu2,5)=ipu1
        endif
 
      elseif(idoc.eq.8) then
C...2 -> 2 processes: store outgoing partons in their CM-frame
        do 500 jt=1,2
          i=mint(84)+2+jt
          kca=jamcomp(mint(20+jt))
          k(i,1)=1
          if(kchg(kca,2).ne.0) k(i,1)=3
          k(i,2)=mint(20+jt)
          k(i,3)=mint(83)+idoc+jt-2
          kfaa=iabs(k(i,2))
          if(mwid(kca).ne.0.and.kfpr(isubsv,1).ne.0) then
            p(i,5)=sqrt(vint(63+mod(js+jt,2)))
          elseif(mwid(kca).ne.0.and.kfpr(isubsv,2).ne.0) then
            p(i,5)=sqrt(vint(64))
          else
            p(i,5)=pjmass(k(i,2))
          endif
          if((kfaa.eq.6.or.kfaa.eq.7.or.kfaa.eq.8).and.
     &    p(i,5).lt.parp(42)) p(i,5)=pjmass(k(i,2))
  500   continue
        if(p(ipu3,5)+p(ipu4,5).ge.shr) then
          kfa1=iabs(mint(21))
          kfa2=iabs(mint(22))
          if((kfa1.gt.3.and.kfa1.ne.21).or.(kfa2.gt.3.and.kfa2.ne.21))
     &    then
            mint(51)=1
            return
          endif
          p(ipu3,5)=0d0
          p(ipu4,5)=0d0
        endif
        p(ipu3,4)=0.5d0*(shr+(p(ipu3,5)**2-p(ipu4,5)**2)/shr)
        p(ipu3,3)=sqrt(max(0d0,p(ipu3,4)**2-p(ipu3,5)**2))
        p(ipu4,4)=shr-p(ipu3,4)
        p(ipu4,3)=-p(ipu3,3)
        n=ipu4
        mint(7)=mint(83)+7
        mint(8)=mint(83)+8
 
C...Rotate outgoing partons using cos(theta)=(th-uh)/lam(sh,sqm3,sqm4)
        call pjrobo(ipu3,ipu4,acos(vint(23)),vint(24),0d0,0d0,0d0)
 
      elseif(idoc.eq.9) then
C...2 -> 3 processes: store outgoing partons in their CM frame
        do 510 jt=1,2
          i=mint(84)+2+jt
          kca=jamcomp(mint(20+jt))
          k(i,1)=1
          if(kchg(kca,2).ne.0) k(i,1)=3
          k(i,2)=mint(20+jt)
          k(i,3)=mint(83)+idoc+jt-3
          if(iabs(k(i,2)).le.22) then
            p(i,5)=pjmass(k(i,2))
          else
            p(i,5)=sqrt(vint(63+mod(js+jt,2)))
          endif
          pt=sqrt(max(0d0,vint(197+5*jt)-p(i,5)**2+vint(196+5*jt)**2))
          p(i,1)=pt*cos(vint(198+5*jt))
          p(i,2)=pt*sin(vint(198+5*jt))
  510   continue
        k(ipu5,1)=1
        k(ipu5,2)=kfres
        k(ipu5,3)=mint(83)+idoc
        p(ipu5,5)=shr
        p(ipu5,1)=-p(ipu3,1)-p(ipu4,1)
        p(ipu5,2)=-p(ipu3,2)-p(ipu4,2)
        pms1=p(ipu3,5)**2+p(ipu3,1)**2+p(ipu3,2)**2
        pms2=p(ipu4,5)**2+p(ipu4,1)**2+p(ipu4,2)**2
        pms3=p(ipu5,5)**2+p(ipu5,1)**2+p(ipu5,2)**2
        pmt3=sqrt(pms3)
        p(ipu5,3)=pmt3*sinh(vint(211))
        p(ipu5,4)=pmt3*cosh(vint(211))
        pms12=(shpr-p(ipu5,4))**2-p(ipu5,3)**2
        sql12=(pms12-pms1-pms2)**2-4d0*pms1*pms2
        if(sql12.le.0d0) then
          mint(51)=1
          return
        endif
        p(ipu3,3)=(-p(ipu5,3)*(pms12+pms1-pms2)+
     &  vint(213)*(shpr-p(ipu5,4))*sqrt(sql12))/(2d0*pms12)
        p(ipu4,3)=-p(ipu3,3)-p(ipu5,3)
        p(ipu3,4)=sqrt(pms1+p(ipu3,3)**2)
        p(ipu4,4)=sqrt(pms2+p(ipu4,3)**2)
        mint(23)=kfres
        n=ipu5
        mint(7)=mint(83)+7
        mint(8)=mint(83)+8
 
      elseif(idoc.eq.11) then
C...Z0 + Z0 -> h0, W+ + W- -> h0: store Higgs and outgoing partons
        phi(1)=paru(2)*pjr(0)
        phi(2)=phi(1)-phir
        do 520 jt=1,2
          i=mint(84)+2+jt
          k(i,1)=1
          if(kchg(jamcomp(mint(20+jt)),2).ne.0) k(i,1)=3
          k(i,2)=mint(20+jt)
          k(i,3)=mint(83)+idoc+jt-2
          p(i,5)=pjmass(k(i,2))
          if(0.5d0*shpr*z(jt).le.p(i,5)) then
            mint(51)=1
            return
          endif
          pabs=sqrt(max(0d0,(0.5d0*shpr*z(jt))**2-p(i,5)**2))
          ptabs=pabs*sqrt(max(0d0,1d0-cthe(jt)**2))
          p(i,1)=ptabs*cos(phi(jt))
          p(i,2)=ptabs*sin(phi(jt))
          p(i,3)=pabs*cthe(jt)*(-1)**(jt+1)
          p(i,4)=0.5d0*shpr*z(jt)
          izw=mint(83)+6+jt
          k(izw,1)=21
          k(izw,2)=23
          if(isub.eq.8) k(izw,2)=isign(24,jamchge(mint(14+jt)))
          k(izw,3)=izw-2
          p(izw,1)=-p(i,1)
          p(izw,2)=-p(i,2)
          p(izw,3)=(0.5d0*shpr-pabs*cthe(jt))*(-1)**(jt+1)
          p(izw,4)=0.5d0*shpr*(1d0-z(jt))
          p(izw,5)=-sqrt(max(0d0,p(izw,3)**2+ptabs**2-p(izw,4)**2))
  520   continue
        i=mint(83)+9
        k(ipu5,1)=1
        k(ipu5,2)=kfres
        k(ipu5,3)=i
        p(ipu5,5)=shr
        p(ipu5,1)=-p(ipu3,1)-p(ipu4,1)
        p(ipu5,2)=-p(ipu3,2)-p(ipu4,2)
        p(ipu5,3)=-p(ipu3,3)-p(ipu4,3)
        p(ipu5,4)=shpr-p(ipu3,4)-p(ipu4,4)
        k(i,1)=21
        k(i,2)=kfres
        do 530 j=1,5
          p(i,j)=p(ipu5,j)
  530   continue
        n=ipu5
        mint(23)=kfres
 
      elseif(idoc.eq.12) then
C...Z0 and W+/- scattering: store bosons and outgoing partons
        phi(1)=paru(2)*pjr(0)
        phi(2)=phi(1)-phir
        jtran=int(1.5d0+pjr(0))
        do 540 jt=1,2
          i=mint(84)+2+jt
          k(i,1)=1
          if(kchg(jamcomp(mint(20+jt)),2).ne.0) k(i,1)=3
          k(i,2)=mint(20+jt)
          k(i,3)=mint(83)+idoc+jt-2
          p(i,5)=pjmass(k(i,2))
          if(0.5d0*shpr*z(jt).le.p(i,5)) p(i,5)=0d0
          pabs=sqrt(max(0d0,(0.5d0*shpr*z(jt))**2-p(i,5)**2))
          ptabs=pabs*sqrt(max(0d0,1d0-cthe(jt)**2))
          p(i,1)=ptabs*cos(phi(jt))
          p(i,2)=ptabs*sin(phi(jt))
          p(i,3)=pabs*cthe(jt)*(-1)**(jt+1)
          p(i,4)=0.5d0*shpr*z(jt)
          izw=mint(83)+6+jt
          k(izw,1)=21
          if(mint(14+jt).eq.mint(20+jt)) then
            k(izw,2)=23
          else
            k(izw,2)=isign(24,jamchge(mint(14+jt))-jamchge(mint(20+jt)))
          endif
          k(izw,3)=izw-2
          p(izw,1)=-p(i,1)
          p(izw,2)=-p(i,2)
          p(izw,3)=(0.5d0*shpr-pabs*cthe(jt))*(-1)**(jt+1)
          p(izw,4)=0.5d0*shpr*(1d0-z(jt))
          p(izw,5)=-sqrt(max(0d0,p(izw,3)**2+ptabs**2-p(izw,4)**2))
          ipu=mint(84)+4+jt
          k(ipu,1)=3
          k(ipu,2)=kfpr(isub,jt)
          if(isub.eq.72.and.jt.eq.jtran) k(ipu,2)=-k(ipu,2)
          if(isub.eq.73.or.isub.eq.77) k(ipu,2)=k(izw,2)
          k(ipu,3)=mint(83)+8+jt
          if(iabs(k(ipu,2)).le.10.or.k(ipu,2).eq.21) then
            p(ipu,5)=pjmass(k(ipu,2))
          else
            p(ipu,5)=sqrt(vint(63+mod(js+jt,2)))
          endif
          mint(22+jt)=k(ipu,2)
  540   continue
C...Find rotation and boost for hard scattering subsystem
        i1=mint(83)+7
        i2=mint(83)+8
        bexcm=(p(i1,1)+p(i2,1))/(p(i1,4)+p(i2,4))
        beycm=(p(i1,2)+p(i2,2))/(p(i1,4)+p(i2,4))
        bezcm=(p(i1,3)+p(i2,3))/(p(i1,4)+p(i2,4))
        gamcm=(p(i1,4)+p(i2,4))/shr
        bepcm=bexcm*p(i1,1)+beycm*p(i1,2)+bezcm*p(i1,3)
        px=p(i1,1)+gamcm*(gamcm/(1d0+gamcm)*bepcm-p(i1,4))*bexcm
        py=p(i1,2)+gamcm*(gamcm/(1d0+gamcm)*bepcm-p(i1,4))*beycm
        pz=p(i1,3)+gamcm*(gamcm/(1d0+gamcm)*bepcm-p(i1,4))*bezcm
        thecm=pjangl(pz,sqrt(px**2+py**2))
        phicm=pjangl(px,py)
C...Store hard scattering subsystem. Rotate and boost it
        sqlam=(sh-p(ipu5,5)**2-p(ipu6,5)**2)**2-4d0*p(ipu5,5)**2*
     &  p(ipu6,5)**2
        pabs=sqrt(max(0d0,sqlam/(4d0*sh)))
        cthwz=vint(23)
        sthwz=sqrt(max(0d0,1d0-cthwz**2))
        phiwz=vint(24)-phicm
        p(ipu5,1)=pabs*sthwz*cos(phiwz)
        p(ipu5,2)=pabs*sthwz*sin(phiwz)
        p(ipu5,3)=pabs*cthwz
        p(ipu5,4)=sqrt(pabs**2+p(ipu5,5)**2)
        p(ipu6,1)=-p(ipu5,1)
        p(ipu6,2)=-p(ipu5,2)
        p(ipu6,3)=-p(ipu5,3)
        p(ipu6,4)=sqrt(pabs**2+p(ipu6,5)**2)
        call pjrobo(ipu5,ipu6,thecm,phicm,bexcm,beycm,bezcm)
        do 560 jt=1,2
          i1=mint(83)+8+jt
          i2=mint(84)+4+jt
          k(i1,1)=21
          k(i1,2)=k(i2,2)
          do 550 j=1,5
            p(i1,j)=p(i2,j)
  550     continue
  560   continue
        n=ipu6
        mint(7)=mint(83)+9
        mint(8)=mint(83)+10
      endif
 
      if(iset(isub).eq.11) then
      elseif(idoc.ge.8) then
C...Store colour connection indices
        do 570 j=1,2
          jc=j
          if(kcs.eq.-1) jc=3-j
          if(icol(kcc,1,jc).ne.0.and.k(ipu1,1).eq.14) k(ipu1,j+3)=
     &    k(ipu1,j+3)+mint(84)+icol(kcc,1,jc)
          if(icol(kcc,2,jc).ne.0.and.k(ipu2,1).eq.14) k(ipu2,j+3)=
     &    k(ipu2,j+3)+mint(84)+icol(kcc,2,jc)
          if(icol(kcc,3,jc).ne.0.and.k(ipu3,1).eq.3) k(ipu3,j+3)=
     &    mstu(5)*(mint(84)+icol(kcc,3,jc))
          if(icol(kcc,4,jc).ne.0.and.k(ipu4,1).eq.3) k(ipu4,j+3)=
     &    mstu(5)*(mint(84)+icol(kcc,4,jc))
  570   continue
 
C...Copy outgoing partons to documentation lines
        imax=2
        if(idoc.eq.9) imax=3
        do 590 i=1,imax
          i1=mint(83)+idoc-imax+i
          i2=mint(84)+2+i
          k(i1,1)=21
          k(i1,2)=k(i2,2)
          if(idoc.le.9) k(i1,3)=0
          if(idoc.ge.11) k(i1,3)=mint(83)+2+i
          do 580 j=1,5
            p(i1,j)=p(i2,j)
  580     continue
  590   continue
 
      elseif(idoc.eq.9) then
C...Store colour connection indices
        do 600 j=1,2
          jc=j
          if(kcs.eq.-1) jc=3-j
          if(icol(kcc,1,jc).ne.0.and.k(ipu1,1).eq.14) k(ipu1,j+3)=
     &    k(ipu1,j+3)+mint(84)+icol(kcc,1,jc)+
     &    max(0,min(1,icol(kcc,1,jc)-2))
          if(icol(kcc,2,jc).ne.0.and.k(ipu2,1).eq.14) k(ipu2,j+3)=
     &    k(ipu2,j+3)+mint(84)+icol(kcc,2,jc)+
     &    max(0,min(1,icol(kcc,2,jc)-2))
          if(icol(kcc,3,jc).ne.0.and.k(ipu4,1).eq.3) k(ipu4,j+3)=
     &    mstu(5)*(mint(84)+icol(kcc,3,jc))
          if(icol(kcc,4,jc).ne.0.and.k(ipu5,1).eq.3) k(ipu5,j+3)=
     &    mstu(5)*(mint(84)+icol(kcc,4,jc))
  600   continue
 
C...Copy outgoing partons to documentation lines
        do 620 i=1,3
          i1=mint(83)+idoc-3+i
          i2=mint(84)+2+i
          k(i1,1)=21
          k(i1,2)=k(i2,2)
          k(i1,3)=0
          do 610 j=1,5
            p(i1,j)=p(i2,j)
  610     continue
  620   continue
      endif
 
C...Low-pT events: remove gluons used for string drawing purposes
      if(isub.eq.95) then
        k(ipu3,1)=k(ipu3,1)+10
        k(ipu4,1)=k(ipu4,1)+10
        do 630 j=41,66
          vintsv(j)=vint(j)
          vint(j)=0d0
  630   continue
        do 650 i=mint(83)+5,mint(83)+8
          do 640 j=1,5
            p(i,j)=0d0
  640     continue
  650   continue
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYSSPA
C...Generates spacelike parton showers.
 
      subroutine pjsspa(ipu1,ipu2)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      common/pjint3/xsfx(2,-40:40),isig(1000,3),sigh(1000)
      save /jyjets/,/jydat1/,/jydat2/,/pjsubs/,/pjpars/,/pjint1/,
     &/pjint2/,/pjint3/
C...Local arrays and data.
      dimension kfls(4),is(2),xs(2),zs(2),q2s(2),tevcsv(2),tevesv(2),
     &xfs(2,-25:25),xfa(-25:25),xfb(-25:25),xfn(-25:25),wtapc(-25:25),
     &wtape(-25:25),wtsf(-25:25),the2(2),alam(2),dq2(3),dpc(3),dpd(4),
     &dpb(4),robo(5),more(2),kfbeam(2),q2mncs(2),kcfi(2),nfis(2),
     &thefis(2,2),isfi(2)
      data is/2*0/
 
C...Read out basic information; set global Q^2 scale.
      ipus1=ipu1
      ipus2=ipu2
      isub=mint(1)
      q2mx=vint(56)
      if(iset(isub).eq.2) q2mx=parp(67)*vint(56)
 
C...Initialize QCD evolution and check phase space.
      q2mnc=parp(62)**2
      q2mncs(1)=q2mnc
      if(mstp(66).eq.1.and.mint(107).eq.3)
     &q2mncs(1)=max(q2mnc,vint(283))
      q2mncs(2)=q2mnc
      if(mstp(66).eq.1.and.mint(108).eq.3)
     &q2mncs(2)=max(q2mnc,vint(284))
      mcev=0
      xec0=2d0*parp(65)/vint(1)
      alams=paru(112)
      paru(112)=parp(61)
      fq2c=1d0
      tcmx=0d0
      if(mint(47).ge.2.and.(mint(47).ne.5.or.mstp(12).ge.1)) then
        mcev=1
        if(mstp(64).eq.1) fq2c=parp(63)
        if(mstp(64).eq.2) fq2c=parp(64)
        tcmx=log(fq2c*q2mx/parp(61)**2)
        if(q2mx.lt.max(q2mnc,2d0*parp(61)**2).or.tcmx.lt.0.2d0)
     &  mcev=0
      endif
 
C...Initialize QED evolution and check phase space.
      q2mne=parp(68)**2
      meev=0
      xee=1d-6
      spme=pmas(11,1)**2
      temx=0d0
      fwte=10d0
      if(mint(45).eq.3.or.mint(46).eq.3) then
        meev=1
        temx=log(q2mx/spme)
        if(q2mx.le.q2mne.or.temx.lt.0.2d0) meev=0
      endif
      if(mcev.eq.0.and.meev.eq.0) return
 
C...Initial values: flavours, momenta, virtualities.
      ns=n
  100 n=ns
      do 120 jt=1,2
        more(jt)=1
        kfbeam(jt)=mint(10+jt)
        if(mint(18+jt).eq.1)kfbeam(jt)=22
        kfls(jt)=mint(14+jt)
        kfls(jt+2)=kfls(jt)
        xs(jt)=vint(40+jt)
        if(mint(18+jt).eq.1) xs(jt)=vint(40+jt)/vint(154+jt)
        zs(jt)=1d0
        q2s(jt)=q2mx
        tevcsv(jt)=tcmx
        alam(jt)=parp(61)
        the2(jt)=100d0
        tevesv(jt)=temx
        do 110 kfl=-25,25
          xfs(jt,kfl)=xsfx(jt,kfl)
  110   continue
  120 continue
      dsh=vint(44)
      if(iset(isub).ge.3.and.iset(isub).le.5) dsh=vint(26)*vint(2)
 
C...Find if interference with final state partons.
      mfis=0
      if(mstp(67).ge.1.and.mstp(67).le.3) mfis=mstp(67)
      if(mfis.ne.0) then
        do 140 i=1,2
          kcfi(i)=0
          kca=jamcomp(iabs(kfls(i)))
          if(kca.ne.0) kcfi(i)=kchg(kca,2)*isign(1,kfls(i))
          nfis(i)=0
          if(kcfi(i).ne.0) then
            if(i.eq.1) ipfs=ipus1
            if(i.eq.2) ipfs=ipus2
            do 130 j=1,2
              icsi=mod(k(ipfs,3+j),mstu(5))
              if(icsi.gt.0.and.icsi.ne.ipus1.and.icsi.ne.ipus2.and.
     &        (kcfi(i).eq.(-1)**(j+1).or.kcfi(i).eq.2)) then
                nfis(i)=nfis(i)+1
                thefis(i,nfis(i))=pjangl(p(icsi,3),sqrt(p(icsi,1)**2+
     &          p(icsi,2)**2))
                if(i.eq.2) thefis(i,nfis(i))=paru(1)-thefis(i,nfis(i))
              endif
  130       continue
          endif
  140   continue
        if(nfis(1)+nfis(2).eq.0) mfis=0
      endif
 
C...Pick up leg with highest virtuality.
  150 n=n+1
      jt=1
      if(n.gt.ns+1.and.q2s(2).gt.q2s(1)) jt=2
      if(more(jt).eq.0) jt=3-jt
      kflb=kfls(jt)
      xb=xs(jt)
      do 160 kfl=-25,25
        xfb(kfl)=xfs(jt,kfl)
  160 continue
      dshr=2d0*sqrt(dsh)
      dshz=dsh/zs(jt)
 
C...Check if allowed to branch.
      mcev=0
      if(iabs(kflb).le.10.or.kflb.eq.21) then
        mcev=1
        xec=max(xec0,xb*(1d0/(1d0-parp(66))-1d0))
        if(xb.ge.1d0-2d0*xec) mcev=0
      endif
      meev=0
      if(mint(44+jt).eq.3) then
        meev=1
        if(xb.ge.1d0-2d0*xee) meev=0
        if((iabs(kflb).le.10.or.kflb.eq.21).and.xb.ge.1d0-2d0*xec)
     &  meev=0
C***Currently kill QED shower for resolved photoproduction.
        if(mint(18+jt).eq.1) meev=0
C***Currently kill shower for W inside electron.
        if(iabs(kflb).eq.24) then
          mcev=0
          meev=0
        endif
      endif
      if(mcev.eq.0.and.meev.eq.0) then
        q2b=0d0
        goto 250
      endif
 
C...Maximum Q2 with or without Q2 ordering. Effective Lambda and n_f.
      q2b=q2s(jt)
      tevcb=tevcsv(jt)
      teveb=tevesv(jt)
      if(mstp(62).le.1) then
        if(zs(jt).gt.0.99999d0) then
          q2b=q2s(jt)
        else
          q2b=0.5d0*(1d0/zs(jt)+1d0)*q2s(jt)+0.5d0*(1d0/zs(jt)-1d0)*
     &    (q2s(3-jt)-dsh+sqrt((dsh+q2s(1)+q2s(2))**2+
     &    8d0*q2s(1)*q2s(2)*zs(jt)/(1d0-zs(jt))))
        endif
        if(mcev.eq.1) tevcb=log(fq2c*q2b/alam(jt)**2)
        if(meev.eq.1) teveb=log(q2b/spme)
      endif
      if(mcev.eq.1) then
        alsdum=pjalps(fq2c*q2b)
        tevcb=tevcb+2d0*log(alam(jt)/paru(117))
        alam(jt)=paru(117)
        b0=(33d0-2d0*mstu(118))/6d0
      endif
      tevcbs=tevcb
      tevebs=teveb
 
C...Select side for interference with final state partons.
      if(mfis.ge.1.and.n.le.ns+2) then
        ifi=n-ns
        isfi(ifi)=0
        if(iabs(kcfi(ifi)).eq.1.and.nfis(ifi).eq.1) then
          isfi(ifi)=1
        elseif(kcfi(ifi).eq.2.and.nfis(ifi).eq.1) then
          if(pjr(0).gt.0.5d0) isfi(ifi)=1
        elseif(kcfi(ifi).eq.2.and.nfis(ifi).eq.2) then
          isfi(ifi)=1
          if(pjr(0).gt.0.5d0) isfi(ifi)=2
        endif
      endif
 
C...Calculate Altarelli-Parisi weights.
      do 170 kfl=-25,25
        wtapc(kfl)=0d0
        wtape(kfl)=0d0
        wtsf(kfl)=0d0
  170 continue
C...q -> q, g -> q.
      if(iabs(kflb).le.10) then
        wtapc(kflb)=(8d0/3d0)*log((1d0-xec-xb)*(xb+xec)/(xec*(1d0-xec)))
        wtapc(21)=0.5d0*(xb/(xb+xec)-xb/(1d0-xec))
C...f -> f, gamma -> f.
      elseif(iabs(kflb).le.20) then
        wtapf1=log((1d0-xee-xb)*(xb+xee)/(xee*(1d0-xee)))
        wtapf2=log((1d0-xee-xb)*(1d0-xee)/(xee*(xb+xee)))
        wtape(kflb)=2d0*(wtapf1+wtapf2)
        if(mstp(12).ge.1) wtape(22)=xb/(xb+xee)-xb/(1d0-xee)
C...f -> g, g -> g.
      elseif(kflb.eq.21) then
        wtapq=(16d0/3d0)*(sqrt((1d0-xec)/xb)-sqrt((xb+xec)/xb))
        do 180 kfl=1,mstp(58)
          wtapc(kfl)=wtapq
          wtapc(-kfl)=wtapq
  180   continue
        wtapc(21)=6d0*log((1d0-xec-xb)/xec)
C...f -> gamma, W+, W-.
      elseif(kflb.eq.22) then
        wtapf=log((1d0-xee-xb)*(1d0-xee)/(xee*(xb+xee)))/xb
        wtape(11)=wtapf
        wtape(-11)=wtapf
      elseif(kflb.eq.24) then
        wtape(-11)=1d0/(4d0*paru(102))*log((1d0-xee-xb)*(1d0-xee)/
     &  (xee*(xb+xee)))/xb
      elseif(kflb.eq.-24) then
        wtape(11)=1d0/(4d0*paru(102))*log((1d0-xee-xb)*(1d0-xee)/
     &  (xee*(xb+xee)))/xb
      endif
 
C...Calculate parton distribution weights and sum.
      ntry=0
  190 ntry=ntry+1
      if(ntry.gt.500) then
        mint(51)=1
        return
      endif
      wtsumc=0d0
      wtsume=0d0
      xfbo=max(1d-10,xfb(kflb))
      do 200 kfl=-25,25
        wtsf(kfl)=xfb(kfl)/xfbo
        wtsumc=wtsumc+wtapc(kfl)*wtsf(kfl)
        wtsume=wtsume+wtape(kfl)*wtsf(kfl)
  200 continue
      wtsumc=max(0.0001d0,wtsumc)
      wtsume=max(0.0001d0/fwte,wtsume)
 
C...Choose new t: fix alpha_s, alpha_s(Q^2), alpha_s(k_T^2).
      ntry2=0
  210 ntry2=ntry2+1
      if(ntry2.gt.500) then
        mint(51)=1
        return
      endif
      if(mcev.eq.1) then
        if(mstp(64).le.0) then
          tevcb=tevcb+log(pjr(0))*paru(2)/(paru(111)*wtsumc)
        elseif(mstp(64).eq.1) then
          tevcb=tevcb*exp(max(-50d0,log(pjr(0))*b0/wtsumc))
        else
          tevcb=tevcb*exp(max(-50d0,log(pjr(0))*b0/(5d0*wtsumc)))
        endif
      endif
      if(meev.eq.1) then
        teveb=teveb*exp(max(-50d0,log(pjr(0))*paru(2)/
     &  (paru(101)*fwte*wtsume*temx)))
      endif
 
C...Translate t into Q2 scale; choose between QCD and QED evolution.
  220 if(mcev.eq.1) q2cb=alam(jt)**2*exp(max(-50d0,tevcb))/fq2c
      if(meev.eq.1) q2eb=spme*exp(max(-50d0,teveb))
      mce=0
      if(mcev.eq.0.and.meev.eq.0) then
      elseif(mcev.eq.1.and.meev.eq.0) then
        if(q2cb.gt.q2mncs(jt)) mce=1
      elseif(mcev.eq.0.and.meev.eq.1) then
        if(q2eb.gt.q2mne) mce=2
      elseif(q2mncs(jt).gt.q2mne) then
        mce=1
        if(q2eb.gt.q2cb.or.q2cb.le.q2mncs(jt)) mce=2
        if(mce.eq.2.and.q2eb.le.q2mne) mce=0
      else
        mce=2
        if(q2cb.gt.q2eb.or.q2eb.le.q2mne) mce=1
        if(mce.eq.1.and.q2cb.le.q2mncs(jt)) mce=0
      endif
 
C...Evolution possibly ended. Update t values.
      if(mce.eq.0) then
        q2b=0d0
        goto 250
      elseif(mce.eq.1) then
        q2b=q2cb
        q2ref=fq2c*q2b
        if(meev.eq.1) teveb=log(q2b/spme)
      else
        q2b=q2eb
        q2ref=q2b
        if(mcev.eq.1) tevcb=log(fq2c*q2b/alam(jt)**2)
      endif
 
C...Select flavour for branching parton.
      if(mce.eq.1) wtran=pjr(0)*wtsumc
      if(mce.eq.2) wtran=pjr(0)*wtsume
      kfla=-25
  230 kfla=kfla+1
      if(mce.eq.1) wtran=wtran-wtapc(kfla)*wtsf(kfla)
      if(mce.eq.2) wtran=wtran-wtape(kfla)*wtsf(kfla)
      if(kfla.le.24.and.wtran.gt.0d0) goto 230
      if(kfla.eq.25) then
        q2b=0d0
        goto 250
      endif
 
C...Choose z value and corrective weight.
      wtz=0d0
C...q -> q + g.
      if(iabs(kfla).le.10.and.iabs(kflb).le.10) then
        z=1d0-((1d0-xb-xec)/(1d0-xec))*
     &  (xec*(1d0-xec)/((xb+xec)*(1d0-xb-xec)))**pjr(0)
        wtz=0.5d0*(1d0+z**2)
C...q -> g + q.
      elseif(iabs(kfla).le.10.and.kflb.eq.21) then
        z=xb/(sqrt(xb+xec)+pjr(0)*(sqrt(1d0-xec)-sqrt(xb+xec)))**2
        wtz=0.5d0*(1d0+(1d0-z)**2)*sqrt(z)
C...f -> f + gamma.
      elseif(iabs(kfla).le.20.and.iabs(kflb).le.20) then
        if(wtapf1.gt.pjr(0)*(wtapf1+wtapf2)) then
          z=1d0-((1d0-xb-xee)/(1d0-xee))*
     &    (xee*(1d0-xee)/((xb+xee)*(1d0-xb-xee)))**pjr(0)
        else
          z=xb+xb*(xee/(1d0-xee))*
     &    ((1d0-xb-xee)*(1d0-xee)/(xee*(xb+xee)))**pjr(0)
        endif
        wtz=0.5d0*(1d0+z**2)*(z-xb)/(1d0-xb)
C...f -> gamma + f.
      elseif(iabs(kfla).le.20.and.kflb.eq.22) then
        z=xb+xb*(xee/(1d0-xee))*
     &  ((1d0-xb-xee)*(1d0-xee)/(xee*(xb+xee)))**pjr(0)
        wtz=0.5d0*(1d0+(1d0-z)**2)*xb*(z-xb)/z
C...f -> W+- + f'.
      elseif(iabs(kfla).le.20.and.iabs(kflb).eq.24) then
        z=xb+xb*(xee/(1d0-xee))*
     &  ((1d0-xb-xee)*(1d0-xee)/(xee*(xb+xee)))**pjr(0)
        wtz=0.5d0*(1d0+(1d0-z)**2)*(xb*(z-xb)/z)*
     &  (q2b/(q2b+pmas(24,1)**2))
C...g -> q + qbar.
      elseif(kfla.eq.21.and.iabs(kflb).le.10) then
        z=xb/(1d0-xec)+pjr(0)*(xb/(xb+xec)-xb/(1d0-xec))
        wtz=1d0-2d0*z*(1d0-z)
C...g -> g + g.
      elseif(kfla.eq.21.and.kflb.eq.21) then
        z=1d0/(1d0+((1d0-xec-xb)/xb)*(xec/(1d0-xec-xb))**pjr(0))
        wtz=(1d0-z*(1d0-z))**2
C...gamma -> f + fbar.
      elseif(kfla.eq.22.and.iabs(kflb).le.20) then
        z=xb/(1d0-xee)+pjr(0)*(xb/(xb+xee)-xb/(1d0-xee))
        wtz=1d0-2d0*z*(1d0-z)
      endif
      if(mce.eq.2) wtz=(wtz/fwte)*(teveb/temx)
 
C...Option with resummation of soft gluon emission as effective z shift.
      if(mce.eq.1) then
        if(mstp(65).ge.1) then
          rsoft=6d0
          if(kflb.ne.21) rsoft=8d0/3d0
          z=z*(tevcb/tevcsv(jt))**(rsoft*xec/((xb+xec)*b0))
          if(z.le.xb) goto 210
        endif
 
C...Option with alpha_s(k_T^2): demand k_T^2 > cutoff, reweight.
        if(mstp(64).ge.2) then
          if((1d0-z)*q2b.lt.q2mncs(jt)) goto 210
          alprat=tevcb/(tevcb+log(1d0-z))
          if(alprat.lt.5d0*pjr(0)) goto 210
          if(alprat.gt.5d0) wtz=wtz*alprat/5d0
        endif
 
C...Impose angular constraint in first branching from interference
C...with final state partons.
        if(mfis.ge.1.and.n.le.ns+2.and.ntry2.lt.200) then
          the2d=(4d0*q2b)/(dsh*(1d0-z))
          if(n.eq.ns+1.and.isfi(1).ge.1) then
            if(the2d.gt.thefis(1,isfi(1))**2) goto 210
          elseif(n.eq.ns+2.and.isfi(2).ge.1) then
            if(the2d.gt.thefis(2,isfi(2))**2) goto 210
          endif
        endif
 
C...Option with angular ordering requirement.
        if(mstp(62).ge.3.and.ntry2.lt.200) then
          the2t=(4d0*z**2*q2b)/(vint(2)*(1d0-z)*xb**2)
          if(the2t.gt.the2(jt)) goto 210
        endif
      endif
 
C...Weighting with new parton distributions.
      mint(105)=mint(102+jt)
      mint(109)=mint(106+jt)
      if(mstp(57).le.1) then
        call pjpdfu(kfbeam(jt),xb,q2ref,xfn)
      else
        call pjpdfl(kfbeam(jt),xb,q2ref,xfn)
      endif
      xfbn=xfn(kflb)
      if(xfbn.lt.1d-20) then
        if(kfla.eq.kflb) then
          tevcb=tevcbs
          teveb=tevebs
          wtapc(kflb)=0d0
          wtape(kflb)=0d0
          goto 190
        elseif(mce.eq.1.and.tevcbs-tevcb.gt.0.2d0) then
          tevcb=0.5d0*(tevcbs+tevcb)
          goto 220
        elseif(mce.eq.2.and.tevebs-teveb.gt.0.2d0) then
          teveb=0.5d0*(tevebs+teveb)
          goto 220
        else
          xfbn=1d-10
          xfn(kflb)=xfbn
        endif
      endif
      do 240 kfl=-25,25
        xfb(kfl)=xfn(kfl)
  240 continue
      xa=xb/z
      if(mstp(57).le.1) then
        call pjpdfu(kfbeam(jt),xa,q2ref,xfa)
      else
        call pjpdfl(kfbeam(jt),xa,q2ref,xfa)
      endif
      xfan=xfa(kfla)
      if(xfan.lt.1d-20) goto 190
      wtsfa=wtsf(kfla)
      if(wtz*xfan/xfbn.lt.pjr(0)*wtsfa) goto 190
 
C...Define two hard scatterers in their CM-frame.
  250 if(n.eq.ns+2) then
        dq2(jt)=q2b
        dplcm=sqrt((dsh+dq2(1)+dq2(2))**2-4d0*dq2(1)*dq2(2))/dshr
        do 270 jr=1,2
          i=ns+jr
          if(jr.eq.1) ipo=ipus1
          if(jr.eq.2) ipo=ipus2
          do 260 j=1,5
            k(i,j)=0
            p(i,j)=0d0
            v(i,j)=0d0
  260     continue
          k(i,1)=14
          k(i,2)=kfls(jr+2)
          k(i,4)=ipo
          k(i,5)=ipo
          p(i,3)=dplcm*(-1)**(jr+1)
          p(i,4)=(dsh+dq2(3-jr)-dq2(jr))/dshr
          p(i,5)=-sqrt(dq2(jr))
          k(ipo,1)=14
          k(ipo,3)=i
          k(ipo,4)=mod(k(ipo,4),mstu(5))+mstu(5)*i
          k(ipo,5)=mod(k(ipo,5),mstu(5))+mstu(5)*i
  270   continue
 
C...Find maximum allowed mass of timelike parton.
      elseif(n.gt.ns+2) then
        jr=3-jt
        dq2(3)=q2b
        dpc(1)=p(is(1),4)
        dpc(2)=p(is(2),4)
        dpc(3)=0.5d0*(abs(p(is(1),3))+abs(p(is(2),3)))
        dpd(1)=dsh+dq2(jr)+dq2(jt)
        dpd(2)=dshz+dq2(jr)+dq2(3)
        dpd(3)=sqrt(dpd(1)**2-4d0*dq2(jr)*dq2(jt))
        dpd(4)=sqrt(dpd(2)**2-4d0*dq2(jr)*dq2(3))
        ikin=0
        if(q2s(jr).ge.0.25d0*q2mnc.and.dpd(1)-dpd(3).ge.
     &  1d-10*dpd(1)) ikin=1
        if(ikin.eq.0) dmsma=(dq2(jt)/zs(jt)-dq2(3))*
     &  (dsh/(dsh+dq2(jt))-dsh/(dshz+dq2(3)))
        if(ikin.eq.1) dmsma=(dpd(1)*dpd(2)-dpd(3)*dpd(4))/
     &  (2d0*dq2(jr))-dq2(jt)-dq2(3)
 
C...Generate timelike parton shower (if required).
        it=n
        do 280 j=1,5
          k(it,j)=0
          p(it,j)=0d0
          v(it,j)=0d0
  280   continue
        k(it,1)=3
C...f -> f + g (gamma).
        if(iabs(kflb).le.20.and.iabs(kfls(jt+2)).le.20) then
          k(it,2)=21
          if(iabs(kflb).ge.11) k(it,2)=22
C...f -> g (gamma, W+-) + f.
        elseif(iabs(kflb).le.20.and.iabs(kfls(jt+2)).gt.20) then
          k(it,2)=kflb
          if(kfls(jt+2).eq.24) then
            k(it,2)=-12
          elseif(kfls(jt+2).eq.-24) then
            k(it,2)=12
          endif
C...g (gamma) -> f + fbar, g + g.
        else
          k(it,2)=-kfls(jt+2)
          if(kfls(jt+2).gt.20) k(it,2)=kfls(jt+2)
        endif
        p(it,5)=pjmass(k(it,2))
        if(dmsma.le.p(it,5)**2) goto 100
        if(mstp(63).ge.1.and.mce.eq.1) then
          mstj48=mstj(48)
          parj85=parj(85)
          p(it,4)=(dshz-dsh-p(it,5)**2)/dshr
          p(it,3)=sqrt(p(it,4)**2-p(it,5)**2)
          if(mstp(63).eq.1) then
            q2tim=dmsma
          elseif(mstp(63).eq.2) then
            q2tim=min(dmsma,parp(71)*q2s(jt))
          else
            q2tim=dmsma
            mstj(48)=1
            if(ikin.eq.0) dpt2=dmsma*(dshz+dq2(3))/(dsh+dq2(jt))
            if(ikin.eq.1) dpt2=dmsma*(0.5d0*dpd(1)*dpd(2)+0.5d0*dpd(3)*
     &      dpd(4)-dq2(jr)*(dq2(jt)+dq2(3)))/(4d0*dsh*dpc(3)**2)
            parj(85)=sqrt(max(0d0,dpt2))*
     &      (1d0/p(it,4)+1d0/p(is(jt),4))
          endif
          call pjshow(it,0,sqrt(q2tim))
          mstj(48)=mstj48
          parj(85)=parj85
          if(n.ge.it+1) p(it,5)=p(it+1,5)
        endif
 
C...Reconstruct kinematics of branching: timelike parton shower.
        dms=p(it,5)**2
        if(ikin.eq.0) dpt2=(dmsma-dms)*(dshz+dq2(3))/(dsh+dq2(jt))
        if(ikin.eq.1) dpt2=(dmsma-dms)*(0.5d0*dpd(1)*dpd(2)+
     &  0.5d0*dpd(3)*dpd(4)-dq2(jr)*(dq2(jt)+dq2(3)+dms))/
     &  (4d0*dsh*dpc(3)**2)
        if(dpt2.lt.0d0) goto 100
        dpb(1)=(0.5d0*dpd(2)-dpc(jr)*(dshz+dq2(jr)-dq2(jt)-dms)/
     &  dshr)/dpc(3)-dpc(3)
        p(it,1)=sqrt(dpt2)
        p(it,3)=dpb(1)*(-1)**(jt+1)
        p(it,4)=sqrt(dpt2+dpb(1)**2+dms)
        if(n.ge.it+1) then
          dpb(1)=sqrt(dpb(1)**2+dpt2)
          dpb(2)=sqrt(dpb(1)**2+dms)
          dpb(3)=p(it+1,3)
          dpb(4)=sqrt(dpb(3)**2+dms)
          dbez=(dpb(4)*dpb(1)-dpb(3)*dpb(2))/(dpb(4)*dpb(2)-dpb(3)*
     &    dpb(1))
          call pjrobo(it+1,n,0d0,0d0,0d0,0d0,dbez)
          the=pjangl(p(it,3),p(it,1))
          call pjrobo(it+1,n,the,0d0,0d0,0d0,0d0)
        endif
 
C...Reconstruct kinematics of branching: spacelike parton.
        do 290 j=1,5
          k(n+1,j)=0
          p(n+1,j)=0d0
          v(n+1,j)=0d0
  290   continue
        k(n+1,1)=14
        k(n+1,2)=kflb
        p(n+1,1)=p(it,1)
        p(n+1,3)=p(it,3)+p(is(jt),3)
        p(n+1,4)=p(it,4)+p(is(jt),4)
        p(n+1,5)=-sqrt(dq2(3))
 
C...Define colour flow of branching.
        k(is(jt),3)=n+1
        k(it,3)=n+1
        im1=n+1
        im2=n+1
C...f -> f + gamma (Z, W).
        if(iabs(k(it,2)).ge.22) then
          k(it,1)=1
          id1=is(jt)
          id2=is(jt)
C...f -> gamma (Z, W) + f.
        elseif(iabs(k(is(jt),2)).ge.22) then
          id1=it
          id2=it
C...gamma -> q + qbar, g + g.
        elseif(k(n+1,2).eq.22) then
          id1=is(jt)
          id2=it
          im1=id2
          im2=id1
C...q -> q + g.
        elseif(k(n+1,2).gt.0.and.k(n+1,2).ne.21.and.k(it,2).eq.21) then
          id1=it
          id2=is(jt)
C...q -> g + q.
        elseif(k(n+1,2).gt.0.and.k(n+1,2).ne.21) then
          id1=is(jt)
          id2=it
C...qbar -> qbar + g.
        elseif(k(n+1,2).lt.0.and.k(it,2).eq.21) then
          id1=is(jt)
          id2=it
C...qbar -> g + qbar.
        elseif(k(n+1,2).lt.0) then
          id1=it
          id2=is(jt)
C...g -> g + g; g -> q + qbar.
        elseif((k(it,2).eq.21.and.pjr(0).gt.0.5d0).or.k(it,2).lt.0) then
          id1=is(jt)
          id2=it
        else
          id1=it
          id2=is(jt)
        endif
        if(im1.eq.n+1) k(im1,4)=k(im1,4)+id1
        if(im2.eq.n+1) k(im2,5)=k(im2,5)+id2
        k(id1,4)=k(id1,4)+mstu(5)*im1
        k(id2,5)=k(id2,5)+mstu(5)*im2
        if(id1.ne.id2) then
          k(id1,5)=k(id1,5)+mstu(5)*id2
          k(id2,4)=k(id2,4)+mstu(5)*id1
        endif
        n=n+1
 
C...Boost to new CM-frame.
        dbsvx=(p(n,1)+p(is(jr),1))/(p(n,4)+p(is(jr),4))
        dbsvz=(p(n,3)+p(is(jr),3))/(p(n,4)+p(is(jr),4))
        if(dbsvx**2+dbsvz**2.ge.1d0) goto 100
        call pjrobo(ns+1,n,0d0,0d0,-dbsvx,0d0,-dbsvz)
        ir=n+(jt-1)*(is(1)-n)
        call pjrobo(ns+1,n,-pjangl(p(ir,3),p(ir,1)),paru(2)*pjr(0),
     &  0d0,0d0,0d0)
      endif
 
C...Update kinematics variables.
      is(jt)=n
      dq2(jt)=q2b
      if(mstp(62).ge.3) the2(jt)=the2t
      dsh=dshz
 
C...Save quantities; loop back.
      q2s(jt)=q2b
      if((mcev.eq.1.and.q2b.ge.0.25d0*q2mnc).or.
     &(meev.eq.1.and.q2b.ge.q2mne)) then
        kfls(jt+2)=kfls(jt)
        kfls(jt)=kfla
        xs(jt)=xa
        zs(jt)=z
        do 300 kfl=-25,25
          xfs(jt,kfl)=xfa(kfl)
  300   continue
        tevcsv(jt)=tevcb
        tevesv(jt)=teveb
      else
        more(jt)=0
        if(jt.eq.1) ipu1=n
        if(jt.eq.2) ipu2=n
      endif
      if(n.gt.mstu(4)-mstu(32)-10) then
        call pjerrm(11,'(PYSSPA:) no more memory left in PYJETS')
        if(mstu(21).ge.1) n=ns
        if(mstu(21).ge.1) return
      endif
      if(more(1).eq.1.or.more(2).eq.1) goto 150
 
C...Boost hard scattering partons to frame of shower initiators.
      do 310 j=1,3
        robo(j+2)=(p(ns+1,j)+p(ns+2,j))/(p(ns+1,4)+p(ns+2,4))
  310 continue
      k(n+2,1)=1
      do 320 j=1,5
        p(n+2,j)=p(ns+1,j)
  320 continue
      robot=robo(3)**2+robo(4)**2+robo(5)**2
      if(robot.ge.0.999999d0) then
        robot=1.00001d0*sqrt(robot)
        robo(3)=robo(3)/robot
        robo(4)=robo(4)/robot
        robo(5)=robo(5)/robot
      endif
      call pjrobo(n+2,n+2,0d0,0d0,-robo(3),-robo(4),-robo(5))
      robo(2)=pjangl(p(n+2,1),p(n+2,2))
      robo(1)=pjangl(p(n+2,3),sqrt(p(n+2,1)**2+p(n+2,2)**2))
      call pjrobo(mint(83)+5,ns,robo(1),robo(2),robo(3),robo(4),
     &robo(5))
 
C...Store user information. Reset Lambda value.
      k(ipu1,3)=mint(83)+3
      k(ipu2,3)=mint(83)+4
      do 330 jt=1,2
        mint(12+jt)=kfls(jt)
        vint(140+jt)=xs(jt)
        if(mint(18+jt).eq.1) vint(140+jt)=vint(154+jt)*xs(jt)
  330 continue
      paru(112)=alams
 
      return
      end
 
C*********************************************************************
 
C...PYRESD
C...Allows resonances to decay (including parton showers for hadronic
C...channels).
 
      subroutine pjresd(ires)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Parameter statement to help give large particle numbers.
      parameter (ksusy1=1000000,ksusy2=2000000,kexcit=4000000)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/jydat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      common/pjint4/mwid(500),wids(500,5)
      save /jyjets/,/jydat1/,/jydat2/,/jydat3/,/pjsubs/,/pjpars/,
     &/pjint1/,/pjint2/,/pjint4/
C...Local arrays and complex and character variables.
      dimension iref(50,8),kdcy(3),kfl1(3),kfl2(3),kfl3(3),keql(3),
     &kcqm(3),kcq1(3),kcq2(3),kcq3(3),nsd(3),pmmn(3),ilin(6),
     &hgz(3,3),coup(6,4),corl(2,2,2),pk(6,4),pkk(6,6),cthe(3),
     &phi(3),wdtp(0:200),wdte(0:200,0:5),dbezqq(3),dpmo(5),xm(5)
      complex fgk,ha(6,6),hc(6,6)
      real tir,uir
      character code*9,mass*9
 
C...The F, Xi and Xj functions of Gunion and Kunszt
C...(Phys. Rev. D33, 665, plus errata from the authors).
      fgk(i1,i2,i3,i4,i5,i6)=4.*ha(i1,i3)*hc(i2,i6)*(ha(i1,i5)*
     &hc(i1,i4)+ha(i3,i5)*hc(i3,i4))
      digk(dt,du)=-4d0*d34*d56+dt*(3d0*dt+4d0*du)+dt**2*(dt*du/
     &(d34*d56)-2d0*(1d0/d34+1d0/d56)*(dt+du)+2d0*(d34/d56+d56/d34))
      djgk(dt,du)=8d0*(d34+d56)**2-8d0*(d34+d56)*(dt+du)-6d0*dt*du-
     &2d0*dt*du*(dt*du/(d34*d56)-2d0*(1d0/d34+1d0/d56)*(dt+du)+
     &2d0*(d34/d56+d56/d34))
 
C...Some general constants.
      xw=paru(102)
      xwv=xw
      if(mstp(8).ge.2) xw=1d0-(pmas(24,1)/pmas(23,1))**2
      xw1=1d0-xw
      sqmz=pmas(23,1)**2
      gmmz=pmas(23,1)*pmas(23,2)
      sqmw=pmas(24,1)**2
      gmmw=pmas(24,1)*pmas(24,2)
      sh=vint(44)
 
C...Reset original resonance configuration.
      do 100 jt=1,8
        iref(1,jt)=0
  100 continue
 
C...Define initial one, two or three objects for subprocess.
      if(ires.eq.0) then
        isub=mint(1)
        if(iset(isub).eq.1.or.iset(isub).eq.3) then
          iref(1,1)=mint(84)+2+iset(isub)
          iref(1,4)=mint(83)+6+iset(isub)
        elseif(iset(isub).eq.2.or.iset(isub).eq.4) then
          iref(1,1)=mint(84)+1+iset(isub)
          iref(1,2)=mint(84)+2+iset(isub)
          iref(1,4)=mint(83)+5+iset(isub)
          iref(1,5)=mint(83)+6+iset(isub)
        elseif(iset(isub).eq.5) then
          iref(1,1)=mint(84)+3
          iref(1,2)=mint(84)+4
          iref(1,3)=mint(84)+5
          iref(1,4)=mint(83)+7
          iref(1,5)=mint(83)+8
          iref(1,6)=mint(83)+9
        endif
 
C...Define original resonance for odd cases.
      else
        isub=0
        iref(1,1)=ires
      endif
 
C...Check if initial resonance has been moved (in resonance + jet).
      do 120 jt=1,3
        if(iref(1,jt).gt.0) then
          if(k(iref(1,jt),1).gt.10) then
            kfa=iabs(k(iref(1,jt),2))
            if(kfa.ge.6.and.kchg(jamcomp(kfa),2).ne.0) then
              do 110 i=iref(1,jt)+1,n
                if(k(i,1).le.10.and.k(i,2).eq.k(iref(1,jt),2))
     &          iref(1,jt)=i
  110         continue
            else
              kda=mod(k(iref(1,jt),4),mstu(4))
              if(mwid(jamcomp(kfa)).ne.0.and.kda.gt.1) iref(1,jt)=kda
            endif
          endif
        endif
  120 continue
 
C...Loop over decay history.
      np=1
      ip=0
  130 ip=ip+1
      ninh=0
      jtmax=2
      if(iref(ip,2).eq.0) jtmax=1
      if(iref(ip,3).ne.0) jtmax=3
      it4=0
      nsav=n
 
C...Start treatment of one, two or three resonances in parallel.
  140 n=nsav
      do 220 jt=1,jtmax
        id=iref(ip,jt)
        kdcy(jt)=0
        kfl1(jt)=0
        kfl2(jt)=0
        kfl3(jt)=0
        keql(jt)=0
        nsd(jt)=id
 
C...Check whether particle can/is allowed to decay.
        if(id.eq.0) goto 210
        kfa=iabs(k(id,2))
        kca=jamcomp(kfa)
        if(mwid(kca).eq.0) goto 210
        if(k(id,1).gt.10.or.mdcy(kca,1).eq.0) goto 210
        if(kfa.eq.6.or.kfa.eq.7.or.kfa.eq.8.or.kfa.eq.17.or.
     &  kfa.eq.18) it4=it4+1
        k(id,4)=mstu(5)*(k(id,4)/mstu(5))
        k(id,5)=mstu(5)*(k(id,5)/mstu(5))
 
C...Info for selection of decay channel: sign, pairings.
        if(kchg(kca,3).eq.0) then
          ipm=2
        else
          ipm=(5-isign(1,k(id,2)))/2
        endif
        kfb=0
        if(jtmax.eq.2) then
          kfb=iabs(k(iref(ip,3-jt),2))
        elseif(jtmax.eq.3) then
          jt2=jt+1-3*(jt/3)
          kfb=iabs(k(iref(ip,jt2),2))
          if(kfb.ne.kfa) then
            jt2=jt+2-3*((jt+1)/3)
            kfb=iabs(k(iref(ip,jt2),2))
          endif
        endif
 
C...Select decay channel.
        if(isub.eq.1.or.isub.eq.15.or.isub.eq.19.or.isub.eq.22.or.
     &  isub.eq.30.or.isub.eq.35.or.isub.eq.141) mint(61)=1
        call pjwidt(kfa,p(id,5)**2,wdtp,wdte)
        wdte0s=wdte(0,1)+wdte(0,ipm)+wdte(0,4)
        if(kfb.eq.kfa) wdte0s=wdte0s+wdte(0,5)
        if(wdte0s.le.0d0) goto 210
        rkfl=wdte0s*pjr(0)
        idl=0
  150   idl=idl+1
        idc=idl+mdcy(kca,2)-1
        rkfl=rkfl-(wdte(idl,1)+wdte(idl,ipm)+wdte(idl,4))
        if(kfb.eq.kfa) rkfl=rkfl-wdte(idl,5)
        if(idl.lt.mdcy(kca,3).and.rkfl.gt.0d0) goto 150
 
C...Read out flavours and colour charges of decay channel chosen.
        kcqm(jt)=kchg(kca,2)*isign(1,k(id,2))
        if(kcqm(jt).eq.-2) kcqm(jt)=2
        kfl1(jt)=kfdp(idc,1)*isign(1,k(id,2))
        kfc1a=jamcomp(iabs(kfl1(jt)))
        if(kchg(kfc1a,3).eq.0) kfl1(jt)=iabs(kfl1(jt))
        kcq1(jt)=kchg(kfc1a,2)*isign(1,kfl1(jt))
        if(kcq1(jt).eq.-2) kcq1(jt)=2
        kfl2(jt)=kfdp(idc,2)*isign(1,k(id,2))
        kfc2a=jamcomp(iabs(kfl2(jt)))
        if(kchg(kfc2a,3).eq.0) kfl2(jt)=iabs(kfl2(jt))
        kcq2(jt)=kchg(kfc2a,2)*isign(1,kfl2(jt))
        if(kcq2(jt).eq.-2) kcq2(jt)=2
        kfl3(jt)=kfdp(idc,3)*isign(1,k(id,2))
        if(kfl3(jt).ne.0) then
          kfc3a=jamcomp(iabs(kfl3(jt)))
          if(kchg(kfc3a,3).eq.0) kfl3(jt)=iabs(kfl3(jt))
          kcq3(jt)=kchg(kfc3a,2)*isign(1,kfl3(jt))
          if(kcq3(jt).eq.-2) kcq3(jt)=2
        endif
 
C...Set/save further info on channel.
        kdcy(jt)=1
        if(kfb.eq.kfa) keql(jt)=mdme(idc,1)
        nsd(jt)=n
        hgz(jt,1)=vint(111)
        hgz(jt,2)=vint(112)
        hgz(jt,3)=vint(114)
 
C...Select masses; to begin with assume resonances narrow.
        do 170 i=1,3
          p(n+i,5)=0d0
          pmmn(i)=0d0
          if(i.eq.1) then
            kflw=iabs(kfl1(jt))
            kcw=kfc1a
          elseif(i.eq.2) then
            kflw=iabs(kfl2(jt))
            kcw=kfc2a
          elseif(i.eq.3) then
            if(kfl3(jt).eq.0) goto 170
            kflw=iabs(kfl3(jt))
            kcw=kfc3a
          endif
          p(n+i,5)=pmas(kcw,1)
CMRENNA++
C...This prevents SUSY/t particles from becoming too light.
          if(kflw/ksusy1.eq.1.or.kflw/ksusy1.eq.2) then
            pmmn(i)=pmas(kcw,1)
            do 160 idc=mdcy(kcw,2),mdcy(kcw,2)+mdcy(kcw,3)-1
              if(mdme(idc,1).gt.0.and.brat(idc).gt.1e-4) then
                pmsum=pmas(jamcomp(kfdp(idc,1)),1)+
     &          pmas(jamcomp(kfdp(idc,2)),1)
                if(kfdp(idc,3).ne.0) pmsum=pmsum+
     &          pmas(jamcomp(kfdp(idc,3)),1)
                pmmn(i)=min(pmmn(i),pmsum)
              endif
  160       continue
CMRENNA--
          elseif(kflw.eq.6) then
            pmmn(i)=pmas(24,1)+pmas(5,1)
          endif
  170   continue
 
C...Check which two out of three are widest.
        iwid1=1
        iwid2=2
        pwid1=pmas(kfc1a,2)
        pwid2=pmas(kfc2a,2)
        kflw1=iabs(kfl1(jt))
        kflw2=iabs(kfl2(jt))
        if(kfl3(jt).ne.0) then
          pwid3=pmas(kfc3a,2)
          if(pwid3.gt.pwid1.and.pwid2.ge.pwid1) then
            iwid1=3
            pwid1=pwid3
            kflw1=iabs(kfl3(jt))
          elseif(pwid3.gt.pwid2) then
            iwid2=3
            pwid2=pwid3
            kflw2=iabs(kfl3(jt))
          endif
        endif
 
C...If all narrow then only check that masses consistent.
        if(mstp(42).le.0.or.(pwid1.lt.parp(41).and.
     &  pwid2.lt.parp(41))) then
CMRENNA++
C....Handle near degeneracy cases.
          if(kfa/ksusy1.eq.1.or.kfa/ksusy1.eq.2) then
            if(p(n+1,5)+p(n+2,5)+p(n+3,5).gt.p(id,5)) then
              p(n+1,5)=p(id,5)-p(n+2,5)-0.5d0
              if(p(n+1,5).lt.0d0) p(n+1,5)=0d0
            endif
          endif
CMRENNA--
          if(p(n+1,5)+p(n+2,5)+p(n+3,5)+parj(64).gt.p(id,5)) then
            call pjerrm(13,'(PYRESD:) daughter masses too large')
            mint(51)=1
            return
          endif
 
C...For three wide resonances select narrower of three
C...according to BW decoupled from rest.
        else
          pmtot=p(id,5)
          if(kfl3(jt).ne.0) then
            iwid3=6-iwid1-iwid2
            kflw3=iabs(kfl1(jt))+iabs(kfl2(jt))+iabs(kfl3(jt))-
     &      kflw1-kflw2
            loop=0
  180       loop=loop+1
            p(n+iwid3,5)=pjmass(kflw3)
            if(loop.le.10.and. p(n+iwid3,5).le.pmmn(iwid3)) goto 180
            pmtot=pmtot-p(n+iwid3,5)
          endif
C...Select other two correlated within remaining phase space.
          if(ip.eq.1) then
            ckin45=ckin(45)
            ckin47=ckin(47)
            ckin(45)=max(pmmn(iwid1),ckin(45))
            ckin(47)=max(pmmn(iwid2),ckin(47))
            call pjofsh(2,kfa,kflw1,kflw2,pmtot,p(n+iwid1,5),
     &      p(n+iwid2,5))
            ckin(45)=ckin45
            ckin(47)=ckin47
          else
            ckin(49)=pmmn(iwid1)
            ckin(50)=pmmn(iwid2)
            call pjofsh(5,kfa,kflw1,kflw2,pmtot,p(n+iwid1,5),
     &      p(n+iwid2,5))
            ckin(49)=0d0
            ckin(50)=0d0
          endif
          if(mint(51).eq.1) return
        endif
 
C...Begin fill decay products, with colour flow for coloured objects.
        mstu10=mstu(10)
        mstu(10)=1
        mstu(19)=1
 
CMRENNA++
C...1) Three-body decays of SUSY particles (plus special case top).
        if(kfl3(jt).ne.0) then
c...jam
          call pjerrm(30,'(pjresd:)susy remove')
          do 200 i=n+1,n+3
            do 190 j=1,5
              k(i,j)=0
              v(i,j)=0d0
  190       continue
  200     continue
          xm(1)=p(n+1,5)
          xm(2)=p(n+2,5)
          xm(3)=p(n+3,5)
          xm(5)=p(id,5)
cjam      call pjtbdy(xm)
          k(n+1,1)=1
          k(n+1,2)=kfl1(jt)
          k(n+2,1)=1
          k(n+2,2)=kfl2(jt)
          k(n+3,1)=1
          k(n+3,2)=kfl3(jt)

C...Set colour flow for t -> W + b + Z.
          if(kfa.eq.6) then
            k(n+2,1)=3
            isid=4
            if(kcqm(jt).eq.-1) isid=5
            idau=n+2
            k(id,isid)=k(id,isid)+idau
            k(idau,isid)=mstu(5)*id
            
C...Set colour flow in three-body decays - programmed as special cases.
          elseif(kfc2a.le.6) then
            k(n+2,1)=3
            k(n+3,1)=3
            isid=4
            if(kfl2(jt).lt.0) isid=5
            k(n+2,isid)=mstu(5)*(n+3)
            k(n+3,9-isid)=mstu(5)*(n+2)
          endif
          if(kfl1(jt).eq.ksusy1+21) then
            k(n+1,1)=3
            k(n+2,1)=3
            k(n+3,1)=3
            isid=4
            if(kfl2(jt).lt.0) isid=5
            k(n+1,isid)=mstu(5)*(n+2)
            k(n+1,9-isid)=mstu(5)*(n+3)
            k(n+2,isid)=mstu(5)*(n+1)
            k(n+3,9-isid)=mstu(5)*(n+1)
          endif
          if(kfa.eq.ksusy1+21) then
            k(n+2,1)=3
            k(n+3,1)=3
            isid=4
            if(kfl2(jt).lt.0) isid=5
            k(id,isid)=k(id,isid)+(n+2)
            k(id,9-isid)=k(id,9-isid)+(n+3)
            k(n+2,isid)=mstu(5)*id
            k(n+3,9-isid)=mstu(5)*id
          endif
          n=n+3
CMRENNA--
 
C...2) Everything else two-body decay.
        else
          call pj2ent(n+1,kfl1(jt),kfl2(jt),p(id,5))
C...First set colour flow as if mother colour singlet.
          if(kcq1(jt).ne.0) then
            k(n-1,1)=3
            if(kcq1(jt).ne.-1) k(n-1,4)=mstu(5)*n
            if(kcq1(jt).ne.1) k(n-1,5)=mstu(5)*n
          endif
          if(kcq2(jt).ne.0) then
            k(n,1)=3
            if(kcq2(jt).ne.-1) k(n,4)=mstu(5)*(n-1)
            if(kcq2(jt).ne.1) k(n,5)=mstu(5)*(n-1)
          endif
C...Then redirect colour flow if mother (anti)triplet.
          if(kcqm(jt).eq.0) then
          elseif(kcqm(jt).ne.2) then
            isid=4
            if(kcqm(jt).eq.-1) isid=5
            idau=n-1
            if(kcq1(jt).eq.0.or.kcq2(jt).eq.2) idau=n
            k(id,isid)=k(id,isid)+idau
            k(idau,isid)=mstu(5)*id
C...Then redirect colour flow if mother octet.
          elseif(kcq1(jt).eq.0.or.kcq2(jt).eq.0) then
            idau=n-1
            if(kcq1(jt).eq.0) idau=n
            k(id,4)=k(id,4)+idau
            k(id,5)=k(id,5)+idau
            k(idau,4)=mstu(5)*id
            k(idau,5)=mstu(5)*id
          else
            isid=4
            if(kcq1(jt).eq.-1) isid=5
            if(kcq1(jt).eq.2) isid=int(4.5d0+pjr(0))
            k(id,isid)=k(id,isid)+(n-1)
            k(id,9-isid)=k(id,9-isid)+n
            k(n-1,isid)=mstu(5)*id
            k(n,9-isid)=mstu(5)*id
          endif
        endif
 
C...End loop over resonances for daughter flavour and mass selection.
        mstu(10)=mstu10
  210   if(mwid(kca).ne.0.and.(kfl1(jt).eq.0.or.kfl3(jt).ne.0))
     &  ninh=ninh+1
        if(ires.gt.0.and.mwid(kca).ne.0.and.kfl1(jt).eq.0) then
          write(code,'(I9)') k(id,2)
          write(mass,'(F9.3)') p(id,5)
          call pjerrm(3,'(PYRESD:) Failed to decay particle'//
     &    code//' with mass'//mass)
          mint(51)=1
          return
        endif
  220 continue
 
C...Check for allowed combinations. Skip if no decays.
      if(jtmax.eq.1) then
        if(kdcy(1).eq.0) goto 560
      elseif(jtmax.eq.2) then
        if(kdcy(1).eq.0.and.kdcy(2).eq.0) goto 560
        if(keql(1).eq.4.and.keql(2).eq.4) goto 140
        if(keql(1).eq.5.and.keql(2).eq.5) goto 140
      elseif(jtmax.eq.3) then
        if(kdcy(1).eq.0.and.kdcy(2).eq.0.and.kdcy(3).eq.0) goto 560
        if(keql(1).eq.4.and.keql(2).eq.4) goto 140
        if(keql(1).eq.4.and.keql(3).eq.4) goto 140
        if(keql(2).eq.4.and.keql(3).eq.4) goto 140
        if(keql(1).eq.5.and.keql(2).eq.5) goto 140
        if(keql(1).eq.5.and.keql(3).eq.5) goto 140
        if(keql(2).eq.5.and.keql(3).eq.5) goto 140
      endif
 
C...Special case: matrix element option for Z0 decay to quarks.
      if(mstp(48).eq.1.and.isub.eq.1.and.jtmax.eq.1.and.
     &iabs(mint(11)).eq.11.and.iabs(kfl1(1)).le.5) then
 
C...Check consistency of MSTJ options set.
        if(mstj(109).eq.2.and.mstj(110).ne.1) then
          call pjerrm(6,
     &    '(PYRESD:) MSTJ(109) value requires MSTJ(110) = 1')
          mstj(110)=1
        endif
        if(mstj(109).eq.2.and.mstj(111).ne.0) then
          call pjerrm(6,
     &    '(PYRESD) MSTJ(109) value requires MSTJ(111) = 0')
          mstj(111)=0
        endif
 
C...Select alpha_strong behaviour.
        mst111=mstu(111)
        par112=paru(112)
        mstu(111)=mstj(108)
        if(mstj(108).eq.2.and.(mstj(101).eq.0.or.mstj(101).eq.1))
     &  mstu(111)=1
        paru(112)=parj(121)
        if(mstu(111).eq.2) paru(112)=parj(122)
 
C...Find axial fraction in total cross section for scalar gluon model.
        parj(171)=0d0
        if((iabs(mstj(101)).eq.1.and.mstj(109).eq.1).or.
     &  (mstj(101).eq.5.and.mstj(49).eq.1)) then
          poll=1d0-parj(131)*parj(132)
          sff=1d0/(16d0*xw*xw1)
          sfw=p(id,5)**4/((p(id,5)**2-parj(123)**2)**2+
     &    (parj(123)*parj(124))**2)
          sfi=sfw*(1d0-(parj(123)/p(id,5))**2)
          ve=4d0*xw-1d0
          hf1i=sfi*sff*(ve*poll+parj(132)-parj(131))
          hf1w=sfw*sff**2*((ve**2+1d0)*poll+2d0*ve*
     &    (parj(132)-parj(131)))
          kflc=iabs(kfl1(1))
          pmq=pjmass(kflc)
          qf=kchg(kflc,1)/3d0
          vq=1d0
          if(mod(mstj(103),2).eq.1) vq=sqrt(max(0d0,
     &    1d0-(2d0*pmq/p(id,5))**2))
          vf=sign(1d0,qf)-4d0*qf*xw
          rfv=0.5d0*vq*(3d0-vq**2)*(qf**2*poll-2d0*qf*vf*hf1i+
     &    vf**2*hf1w)+vq**3*hf1w
          if(rfv.gt.0d0) parj(171)=min(1d0,vq**3*hf1w/rfv)
        endif
 
C...Choice of jet configuration.
        call pjxjet(p(id,5),njet,cut)
        kflc=iabs(kfl1(1))
        kfln=21
        if(njet.eq.4) then
          call pjx4jt(njet,cut,kflc,p(id,5),kfln,x1,x2,x4,x12,x14)
        elseif(njet.eq.3) then
          call pjx3jt(njet,cut,kflc,p(id,5),x1,x3)
        else
          mstj(120)=1
        endif
 
C...Fill jet configuration; return if incorrect kinematics.
        nc=n-2
        if(njet.eq.2.and.mstj(101).ne.5) then
          call pj2ent(nc+1,kflc,-kflc,p(id,5))
        elseif(njet.eq.2) then
          call pj2ent(-(nc+1),kflc,-kflc,p(id,5))
        elseif(njet.eq.3) then
          call pj3ent(nc+1,kflc,21,-kflc,p(id,5),x1,x3)
        elseif(kfln.eq.21) then
          call pj4ent(nc+1,kflc,kfln,kfln,-kflc,p(id,5),x1,x2,x4,
     &    x12,x14)
        else
          call pj4ent(nc+1,kflc,-kfln,kfln,-kflc,p(id,5),x1,x2,x4,
     &    x12,x14)
        endif
        if(mstu(24).ne.0) then
          mint(51)=1
          mstu(111)=mst111
          paru(112)=par112
          return
        endif
 
C...Angular orientation according to matrix element.
        if(mstj(106).eq.1) then
          call pjxdif(nc,njet,kflc,p(id,5),chi,the,phi)
          if(mint(11).lt.0) the=paru(1)-the
          cthe(1)=cos(the)
          call pjrobo(nc+1,n,0d0,chi,0d0,0d0,0d0)
          call pjrobo(nc+1,n,the,phi,0d0,0d0,0d0)
        endif
 
C...Boost partons to Z0 rest frame.
        call pjrobo(nc+1,n,0d0,0d0,p(id,1)/p(id,4),
     &  p(id,2)/p(id,4),p(id,3)/p(id,4))
 
C...Mark decayed resonance and add documentation lines,
        k(id,1)=k(id,1)+10
        idoc=mint(83)+mint(4)
        do 240 i=nc+1,n
          i1=mint(83)+mint(4)+1
          k(i,3)=i1
          if(mstp(128).ge.1) k(i,3)=id
          if(mstp(128).le.1.and.mint(4).lt.mstp(126)) then
            mint(4)=mint(4)+1
            k(i1,1)=21
            k(i1,2)=k(i,2)
            k(i1,3)=iref(ip,4)
            do 230 j=1,5
              p(i1,j)=p(i,j)
  230       continue
          endif
  240   continue
 
C...Generate parton shower.
        if(mstj(101).eq.5) call pjshow(n-1,n,p(id,5))
 
C... End special case for Z0: skip ahead.
        mstu(111)=mst111
        paru(112)=par112
        goto 550
      endif
 
C...Order incoming partons and outgoing resonances.
      if(jtmax.eq.2.and.mstp(47).ge.1.and.ninh.eq.0) then
        ilin(1)=mint(84)+1
        if(k(mint(84)+1,2).gt.0) ilin(1)=mint(84)+2
        if(k(ilin(1),2).eq.21) ilin(1)=2*mint(84)+3-ilin(1)
        ilin(2)=2*mint(84)+3-ilin(1)
        imin=1
        if(iref(ip,7).eq.25.or.iref(ip,7).eq.35.or.iref(ip,7)
     &  .eq.36) imin=3
        imax=2
        iord=1
        if(k(iref(ip,1),2).eq.23) iord=2
        if(k(iref(ip,1),2).eq.24.and.k(iref(ip,2),2).eq.-24) iord=2
        iakipd=iabs(k(iref(ip,iord),2))
        if(iakipd.eq.25.or.iakipd.eq.35.or.iakipd.eq.36) iord=3-iord
        if(kdcy(iord).eq.0) iord=3-iord
 
C...Order decay products of resonances.
        do 250 jt=iord,3-iord,3-2*iord
          if(kdcy(jt).eq.0) then
            ilin(imax+1)=nsd(jt)
            imax=imax+1
          elseif(k(nsd(jt)+1,2).gt.0) then
            ilin(imax+1)=n+2*jt-1
            ilin(imax+2)=n+2*jt
            imax=imax+2
            k(n+2*jt-1,2)=k(nsd(jt)+1,2)
            k(n+2*jt,2)=k(nsd(jt)+2,2)
          else
            ilin(imax+1)=n+2*jt
            ilin(imax+2)=n+2*jt-1
            imax=imax+2
            k(n+2*jt-1,2)=k(nsd(jt)+1,2)
            k(n+2*jt,2)=k(nsd(jt)+2,2)
          endif
  250   continue
 
C...Find charge, isospin, left- and righthanded couplings.
        do 270 i=imin,imax
          do 260 j=1,4
            coup(i,j)=0d0
  260     continue
          kfa=iabs(k(ilin(i),2))
          if(kfa.eq.0.or.kfa.gt.20) goto 270
          coup(i,1)=kchg(kfa,1)/3d0
          coup(i,2)=(-1)**mod(kfa,2)
          coup(i,4)=-2d0*coup(i,1)*xwv
          coup(i,3)=coup(i,2)+coup(i,4)
  270   continue
 
C...Full propagator dependence and flavour correlations for 2 gamma*/Z.
        if(isub.eq.22) then
          do 300 i=3,5,2
            i1=iord
            if(i.eq.5) i1=3-iord
            do 290 j1=1,2
              do 280 j2=1,2
                corl(i/2,j1,j2)=coup(1,1)**2*hgz(i1,1)*coup(i,1)**2/
     &          16d0+coup(1,1)*coup(1,j1+2)*hgz(i1,2)*coup(i,1)*
     &          coup(i,j2+2)/4d0+coup(1,j1+2)**2*hgz(i1,3)*
     &          coup(i,j2+2)**2
  280         continue
  290       continue
  300     continue
          cowt12=(corl(1,1,1)+corl(1,1,2))*(corl(2,1,1)+corl(2,1,2))+
     &    (corl(1,2,1)+corl(1,2,2))*(corl(2,2,1)+corl(2,2,2))
          comx12=(corl(1,1,1)+corl(1,1,2)+corl(1,2,1)+corl(1,2,2))*
     &    (corl(2,1,1)+corl(2,1,2)+corl(2,2,1)+corl(2,2,2))
          if(cowt12.lt.pjr(0)*comx12) goto 140
        endif
      endif
 
C...Select angular orientation type - Z'/W' only.
      mzpwp=0
      if(isub.eq.141) then
        if(pjr(0).lt.paru(130)) mzpwp=1
        if(ip.eq.2) then
          if(iabs(k(iref(2,1),2)).eq.37) mzpwp=2
          iakir=iabs(k(iref(2,2),2))
          if(iakir.eq.25.or.iakir.eq.35.or.iakir.eq.36) mzpwp=2
        endif
        if(ip.ge.3) mzpwp=2
      elseif(isub.eq.142) then
        if(pjr(0).lt.paru(136)) mzpwp=1
        if(ip.eq.2) then
          iakir=iabs(k(iref(2,2),2))
          if(iakir.eq.25.or.iakir.eq.35.or.iakir.eq.36) mzpwp=2
        endif
        if(ip.ge.3) mzpwp=2
      endif
 
C...Select random angles (begin of weighting procedure).
  310 do 320 jt=1,jtmax
        if(kdcy(jt).eq.0) goto 320
        if(jtmax.eq.1.and.isub.ne.0) then
          cthe(jt)=vint(13)+(vint(33)-vint(13)+vint(34)-vint(14))*pjr(0)
          if(cthe(jt).gt.vint(33)) cthe(jt)=cthe(jt)+vint(14)-vint(33)
          phi(jt)=vint(24)
        else
          cthe(jt)=2d0*pjr(0)-1d0
          phi(jt)=paru(2)*pjr(0)
        endif
  320 continue
 
      if(jtmax.eq.2.and.mstp(47).ge.1.and.ninh.eq.0) then
C...Construct massless four-vectors.
        do 340 i=n+1,n+4
          k(i,1)=1
          do 330 j=1,5
            p(i,j)=0d0
            v(i,j)=0d0
  330     continue
  340   continue
        do 350 jt=1,jtmax
          if(kdcy(jt).eq.0) goto 350
          id=iref(ip,jt)
          p(n+2*jt-1,3)=0.5d0*p(id,5)
          p(n+2*jt-1,4)=0.5d0*p(id,5)
          p(n+2*jt,3)=-0.5d0*p(id,5)
          p(n+2*jt,4)=0.5d0*p(id,5)
          call pjrobo(n+2*jt-1,n+2*jt,acos(cthe(jt)),phi(jt),
     &    p(id,1)/p(id,4),p(id,2)/p(id,4),p(id,3)/p(id,4))
  350   continue
 
C...Store incoming and outgoing momenta, with random rotation to
C...avoid accidental zeroes in HA expressions.
        do 370 i=1,imax
          k(n+4+i,1)=1
          p(n+4+i,4)=sqrt(p(ilin(i),1)**2+p(ilin(i),2)**2+
     &    p(ilin(i),3)**2+p(ilin(i),5)**2)
          p(n+4+i,5)=p(ilin(i),5)
          do 360 j=1,3
            p(n+4+i,j)=p(ilin(i),j)
  360     continue
  370   continue
  380   therr=acos(2d0*pjr(0)-1d0)
        phirr=paru(2)*pjr(0)
        call pjrobo(n+5,n+4+imax,therr,phirr,0d0,0d0,0d0)
        do 400 i=1,imax
          if(p(n+4+i,1)**2+p(n+4+i,2)**2.lt.1d-4*p(n+4+i,4)**2) goto 380
          do 390 j=1,4
            pk(i,j)=p(n+4+i,j)
  390     continue
  400   continue
 
C...Calculate internal products.
        if(isub.eq.22.or.isub.eq.23.or.isub.eq.25.or.isub.eq.141.or.
     &  isub.eq.142) then
          do 420 i1=imin,imax-1
            do 410 i2=i1+1,imax
              ha(i1,i2)=sngl(sqrt((pk(i1,4)-pk(i1,3))*(pk(i2,4)+
     &        pk(i2,3))/(1d-20+pk(i1,1)**2+pk(i1,2)**2)))*
     &        cmplx(sngl(pk(i1,1)),sngl(pk(i1,2)))-
     &        sngl(sqrt((pk(i1,4)+pk(i1,3))*(pk(i2,4)-pk(i2,3))/
     &        (1d-20+pk(i2,1)**2+pk(i2,2)**2)))*
     &        cmplx(sngl(pk(i2,1)),sngl(pk(i2,2)))
              hc(i1,i2)=conjg(ha(i1,i2))
              if(i1.le.2) ha(i1,i2)=cmplx(0.,1.)*ha(i1,i2)
              if(i1.le.2) hc(i1,i2)=cmplx(0.,1.)*hc(i1,i2)
              ha(i2,i1)=-ha(i1,i2)
              hc(i2,i1)=-hc(i1,i2)
  410       continue
  420     continue
        endif
        do 440 i=1,2
          do 430 j=1,4
            pk(i,j)=-pk(i,j)
  430     continue
  440   continue
        do 460 i1=imin,imax-1
          do 450 i2=i1+1,imax
            pkk(i1,i2)=2d0*(pk(i1,4)*pk(i2,4)-pk(i1,1)*pk(i2,1)-
     &      pk(i1,2)*pk(i2,2)-pk(i1,3)*pk(i2,3))
            pkk(i2,i1)=pkk(i1,i2)
  450     continue
  460   continue
      endif
 
      kfagm=iabs(iref(ip,7))
      if(mstp(47).le.0.or.ninh.ne.0) then
C...Isotropic decay selected by user.
        wt=1d0
        wtmax=1d0

      elseif(jtmax.eq.3) then
C...Isotropic decay when three mother particles.
        wt=1d0
        wtmax=1d0      
 
      elseif(it4.ge.1) then
C... Isotropic decay t -> b + W etc for 4th generation q and l.
        wt=1d0
        wtmax=1d0
 
      elseif(iref(ip,7).eq.25.or.iref(ip,7).eq.35.or.
     &  iref(ip,7).eq.36) then
C...Angular weight for h0 -> Z0 + Z0 or W+ + W- -> 4 quarks/leptons.
        if(ip.eq.1) wtmax=sh**2
        if(ip.ge.2) wtmax=p(iref(ip,8),5)**4
        kfa=iabs(k(iref(ip,1),2))
        if(kfa.eq.23) then
          kflf1a=iabs(kfl1(1))
          ef1=kchg(kflf1a,1)/3d0
          af1=sign(1d0,ef1+0.1d0)
          vf1=af1-4d0*ef1*xwv
          kflf2a=iabs(kfl1(2))
          ef2=kchg(kflf2a,1)/3d0
          af2=sign(1d0,ef2+0.1d0)
          vf2=af2-4d0*ef2*xwv
          va12as=4d0*vf1*af1*vf2*af2/((vf1**2+af1**2)*(vf2**2+af2**2))
          wt=8d0*(1d0+va12as)*pkk(3,5)*pkk(4,6)+
     &    8d0*(1d0-va12as)*pkk(3,6)*pkk(4,5)
        elseif(kfa.eq.24) then
          wt=16d0*pkk(3,5)*pkk(4,6)
        else
          wt=wtmax
        endif
 
      elseif((kfagm.eq.6.or.kfagm.eq.7.or.kfagm.eq.8.or.
     &  kfagm.eq.17.or.kfagm.eq.18).and.iabs(k(iref(ip,1),2)).eq.24)
     &  then
C...Angular correlation in f -> f' + W -> f' + 2 quarks/leptons.
        i1=iref(ip,8)
        if(mod(kfagm,2).eq.0) then
          i2=n+1
          i3=n+2
        else
          i2=n+2
          i3=n+1
        endif
        i4=iref(ip,2)
        wt=(p(i1,4)*p(i2,4)-p(i1,1)*p(i2,1)-p(i1,2)*p(i2,2)-
     &  p(i1,3)*p(i2,3))*(p(i3,4)*p(i4,4)-p(i3,1)*p(i4,1)-
     &  p(i3,2)*p(i4,2)-p(i3,3)*p(i4,3))
        wtmax=(p(i1,5)**4-p(iref(ip,1),5)**4)/8d0
 
      elseif(isub.eq.1) then
C...Angular weight for gamma*/Z0 -> 2 quarks/leptons.
        ei=kchg(iabs(mint(15)),1)/3d0
        ai=sign(1d0,ei+0.1d0)
        vi=ai-4d0*ei*xwv
        ef=kchg(iabs(kfl1(1)),1)/3d0
        af=sign(1d0,ef+0.1d0)
        vf=af-4d0*ef*xwv
        rmf=min(1d0,4d0*pmas(iabs(kfl1(1)),1)**2/sh)
        wt1=ei**2*vint(111)*ef**2+ei*vi*vint(112)*ef*vf+
     &  (vi**2+ai**2)*vint(114)*(vf**2+(1d0-rmf)*af**2)
        wt2=rmf*(ei**2*vint(111)*ef**2+ei*vi*vint(112)*ef*vf+
     &  (vi**2+ai**2)*vint(114)*vf**2)
        wt3=sqrt(1d0-rmf)*(ei*ai*vint(112)*ef*af+
     &  4d0*vi*ai*vint(114)*vf*af)
        wt=wt1*(1d0+cthe(1)**2)+wt2*(1d0-cthe(1)**2)+
     &  2d0*wt3*cthe(1)*isign(1,mint(15)*kfl1(1))
        wtmax=2d0*(wt1+abs(wt3))
 
      elseif(isub.eq.2) then
C...Angular weight for W+/- -> 2 quarks/leptons.
        wt=(1d0+cthe(1)*isign(1,mint(15)*kfl1(1)))**2
        wtmax=4d0
 
      elseif(isub.eq.15.or.isub.eq.19) then
C...Angular weight for f + fbar -> gluon/gamma + (gamma*/Z0) ->
C...-> gluon/gamma + 2 quarks/leptons.
        clilf=coup(1,1)**2*hgz(2,1)*coup(3,1)**2/16d0+
     &  coup(1,1)*coup(1,3)*hgz(2,2)*coup(3,1)*coup(3,3)/4d0+
     &  coup(1,3)**2*hgz(2,3)*coup(3,3)**2
        clirf=coup(1,1)**2*hgz(2,1)*coup(3,1)**2/16d0+
     &  coup(1,1)*coup(1,3)*hgz(2,2)*coup(3,1)*coup(3,4)/4d0+
     &  coup(1,3)**2*hgz(2,3)*coup(3,4)**2
        crilf=coup(1,1)**2*hgz(2,1)*coup(3,1)**2/16d0+
     &  coup(1,1)*coup(1,4)*hgz(2,2)*coup(3,1)*coup(3,3)/4d0+
     &  coup(1,4)**2*hgz(2,3)*coup(3,3)**2
        crirf=coup(1,1)**2*hgz(2,1)*coup(3,1)**2/16d0+
     &  coup(1,1)*coup(1,4)*hgz(2,2)*coup(3,1)*coup(3,4)/4d0+
     &  coup(1,4)**2*hgz(2,3)*coup(3,4)**2
        wt=(clilf+crirf)*(pkk(1,3)**2+pkk(2,4)**2)+
     &  (clirf+crilf)*(pkk(1,4)**2+pkk(2,3)**2)
        wtmax=(clilf+clirf+crilf+crirf)*
     &  ((pkk(1,3)+pkk(1,4))**2+(pkk(2,3)+pkk(2,4))**2)
 
      elseif(isub.eq.16.or.isub.eq.20) then
C...Angular weight for f + fbar' -> gluon/gamma + W+/- ->
C...-> gluon/gamma + 2 quarks/leptons.
        wt=pkk(1,3)**2+pkk(2,4)**2
        wtmax=(pkk(1,3)+pkk(1,4))**2+(pkk(2,3)+pkk(2,4))**2
 
      elseif(isub.eq.22) then
C...Angular weight for f + fbar -> Z0 + Z0 -> 4 quarks/leptons.
        s34=p(iref(ip,iord),5)**2
        s56=p(iref(ip,3-iord),5)**2
        ti=pkk(1,3)+pkk(1,4)+s34
        ui=pkk(1,5)+pkk(1,6)+s56
        tir=real(ti)
        uir=real(ui)
        fgk135=abs(fgk(1,2,3,4,5,6)/tir+fgk(1,2,5,6,3,4)/uir)**2
        fgk145=abs(fgk(1,2,4,3,5,6)/tir+fgk(1,2,5,6,4,3)/uir)**2
        fgk136=abs(fgk(1,2,3,4,6,5)/tir+fgk(1,2,6,5,3,4)/uir)**2
        fgk146=abs(fgk(1,2,4,3,6,5)/tir+fgk(1,2,6,5,4,3)/uir)**2
        fgk253=abs(fgk(2,1,5,6,3,4)/tir+fgk(2,1,3,4,5,6)/uir)**2
        fgk263=abs(fgk(2,1,6,5,3,4)/tir+fgk(2,1,3,4,6,5)/uir)**2
        fgk254=abs(fgk(2,1,5,6,4,3)/tir+fgk(2,1,4,3,5,6)/uir)**2
        fgk264=abs(fgk(2,1,6,5,4,3)/tir+fgk(2,1,4,3,6,5)/uir)**2
        wt=
     &  corl(1,1,1)*corl(2,1,1)*fgk135+corl(1,1,2)*corl(2,1,1)*fgk145+
     &  corl(1,1,1)*corl(2,1,2)*fgk136+corl(1,1,2)*corl(2,1,2)*fgk146+
     &  corl(1,2,1)*corl(2,2,1)*fgk253+corl(1,2,2)*corl(2,2,1)*fgk263+
     &  corl(1,2,1)*corl(2,2,2)*fgk254+corl(1,2,2)*corl(2,2,2)*fgk264
        wtmax=16d0*((corl(1,1,1)+corl(1,1,2))*(corl(2,1,1)+corl(2,1,2))+
     &  (corl(1,2,1)+corl(1,2,2))*(corl(2,2,1)+corl(2,2,2)))*s34*s56*
     &  ((ti**2+ui**2+2d0*sh*(s34+s56))/(ti*ui)-s34*s56*(1d0/ti**2+
     &  1d0/ui**2))
 
      elseif(isub.eq.23) then
C...Angular weight for f + fbar' -> Z0 + W+/- -> 4 quarks/leptons.
        d34=p(iref(ip,iord),5)**2
        d56=p(iref(ip,3-iord),5)**2
        dt=pkk(1,3)+pkk(1,4)+d34
        du=pkk(1,5)+pkk(1,6)+d56
        facbw=1d0/((sh-sqmw)**2+sqmw*pmas(24,2)**2)
        cawz=coup(2,3)/dt-2d0*xw1*coup(1,2)*(sh-sqmw)*facbw
        cbwz=coup(1,3)/du+2d0*xw1*coup(1,2)*(sh-sqmw)*facbw
        fgk135=abs(real(cawz)*fgk(1,2,3,4,5,6)+
     &  real(cbwz)*fgk(1,2,5,6,3,4))
        fgk136=abs(real(cawz)*fgk(1,2,3,4,6,5)+
     &  real(cbwz)*fgk(1,2,6,5,3,4))
        wt=(coup(5,3)*fgk135)**2+(coup(5,4)*fgk136)**2
        wtmax=4d0*d34*d56*(coup(5,3)**2+coup(5,4)**2)*(cawz**2*
     &  digk(dt,du)+cbwz**2*digk(du,dt)+cawz*cbwz*djgk(dt,du))
 
      elseif(isub.eq.24.or.isub.eq.171.or.isub.eq.176) then
C...Angular weight for f + fbar -> Z0 + h0 -> 2 quarks/leptons + h0
C...(or H0, or A0).
        wt=((coup(1,3)*coup(3,3))**2+(coup(1,4)*coup(3,4))**2)*
     &  pkk(1,3)*pkk(2,4)+((coup(1,3)*coup(3,4))**2+(coup(1,4)*
     &  coup(3,3))**2)*pkk(1,4)*pkk(2,3)
        wtmax=(coup(1,3)**2+coup(1,4)**2)*(coup(3,3)**2+coup(3,4)**2)*
     &  (pkk(1,3)+pkk(1,4))*(pkk(2,3)+pkk(2,4))
 
      elseif(isub.eq.25) then
C...Angular weight for f + fbar -> W+ + W- -> 4 quarks/leptons.
        d34=p(iref(ip,iord),5)**2
        d56=p(iref(ip,3-iord),5)**2
        dt=pkk(1,3)+pkk(1,4)+d34
        du=pkk(1,5)+pkk(1,6)+d56
        facbw=1d0/((sh-sqmz)**2+sqmz*pmas(23,2)**2)
        cdww=(coup(1,3)*sqmz*(sh-sqmz)*facbw+coup(1,2))/sh
        caww=cdww+0.5d0*(coup(1,2)+1d0)/dt
        cbww=cdww+0.5d0*(coup(1,2)-1d0)/du
        ccww=coup(1,4)*sqmz*(sh-sqmz)*facbw/sh
        fgk135=abs(real(caww)*fgk(1,2,3,4,5,6)-
     &  real(cbww)*fgk(1,2,5,6,3,4))
        fgk253=abs(fgk(2,1,5,6,3,4)-fgk(2,1,3,4,5,6))
        wt=fgk135**2+(ccww*fgk253)**2
        wtmax=4d0*d34*d56*(caww**2*digk(dt,du)+cbww**2*digk(du,dt)-caww*
     &  cbww*djgk(dt,du)+ccww**2*(digk(dt,du)+digk(du,dt)-djgk(dt,du)))
 
      elseif(isub.eq.26.or.isub.eq.172.or.isub.eq.177) then
C...Angular weight for f + fbar' -> W+/- + h0 -> 2 quarks/leptons + h0
C...(or H0, or A0).
        wt=pkk(1,3)*pkk(2,4)
        wtmax=(pkk(1,3)+pkk(1,4))*(pkk(2,3)+pkk(2,4))
 
      elseif(isub.eq.30.or.isub.eq.35) then
C...Angular weight for f + g/gamma -> f + (gamma*/Z0)
C...-> f + 2 quarks/leptons.
        clilf=coup(1,1)**2*hgz(2,1)*coup(3,1)**2/16d0+
     &  coup(1,1)*coup(1,3)*hgz(2,2)*coup(3,1)*coup(3,3)/4d0+
     &  coup(1,3)**2*hgz(2,3)*coup(3,3)**2
        clirf=coup(1,1)**2*hgz(2,1)*coup(3,1)**2/16d0+
     &  coup(1,1)*coup(1,3)*hgz(2,2)*coup(3,1)*coup(3,4)/4d0+
     &  coup(1,3)**2*hgz(2,3)*coup(3,4)**2
        crilf=coup(1,1)**2*hgz(2,1)*coup(3,1)**2/16d0+
     &  coup(1,1)*coup(1,4)*hgz(2,2)*coup(3,1)*coup(3,3)/4d0+
     &  coup(1,4)**2*hgz(2,3)*coup(3,3)**2
        crirf=coup(1,1)**2*hgz(2,1)*coup(3,1)**2/16d0+
     &  coup(1,1)*coup(1,4)*hgz(2,2)*coup(3,1)*coup(3,4)/4d0+
     &  coup(1,4)**2*hgz(2,3)*coup(3,4)**2
        if(k(ilin(1),2).gt.0) wt=(clilf+crirf)*(pkk(1,4)**2+
     &  pkk(3,5)**2)+(clirf+crilf)*(pkk(1,3)**2+pkk(4,5)**2)
        if(k(ilin(1),2).lt.0) wt=(clilf+crirf)*(pkk(1,3)**2+
     &  pkk(4,5)**2)+(clirf+crilf)*(pkk(1,4)**2+pkk(3,5)**2)
        wtmax=(clilf+clirf+crilf+crirf)*
     &  ((pkk(1,3)+pkk(1,4))**2+(pkk(3,5)+pkk(4,5))**2)
 
      elseif(isub.eq.31) then
C...Angular weight for f + g -> f' + W+/- -> f' + 2 quarks/leptons.
        if(k(ilin(1),2).gt.0) wt=pkk(1,4)**2+pkk(3,5)**2
        if(k(ilin(1),2).lt.0) wt=pkk(1,3)**2+pkk(4,5)**2
        wtmax=(pkk(1,3)+pkk(1,4))**2+(pkk(3,5)+pkk(4,5))**2
 
      elseif(isub.eq.71.or.isub.eq.72.or.isub.eq.73.or.isub.eq.76.or.
     &  isub.eq.77) then
C...Angular weight for V_L1 + V_L2 -> V_L3 + V_L4 (V = Z/W).
        wt=16d0*pkk(3,5)*pkk(4,6)
        wtmax=sh**2
 
      elseif(isub.eq.110) then
C...Angular weight for f + fbar -> gamma + h0 -> gamma + X is isotropic.
        wt=1d0
        wtmax=1d0
 
      elseif(isub.eq.141) then
        if(ip.eq.1.and.iabs(kfl1(1)).lt.20) then
C...Angular weight for f + fbar -> gamma*/Z0/Z'0 -> 2 quarks/leptons.
C...Couplings of incoming flavour.
          kfai=iabs(mint(15))
          ei=kchg(kfai,1)/3d0
          ai=sign(1d0,ei+0.1d0)
          vi=ai-4d0*ei*xwv
          kfaic=1
          if(kfai.le.10.and.mod(kfai,2).eq.0) kfaic=2
          if(kfai.gt.10.and.mod(kfai,2).ne.0) kfaic=3
          if(kfai.gt.10.and.mod(kfai,2).eq.0) kfaic=4
          vpi=paru(119+2*kfaic)
          api=paru(120+2*kfaic)
C...Couplings of final flavour.
          kfaf=iabs(kfl1(1))
          ef=kchg(kfaf,1)/3d0
          af=sign(1d0,ef+0.1d0)
          vf=af-4d0*ef*xwv
          kfafc=1
          if(kfaf.le.10.and.mod(kfaf,2).eq.0) kfafc=2
          if(kfaf.gt.10.and.mod(kfaf,2).ne.0) kfafc=3
          if(kfaf.gt.10.and.mod(kfaf,2).eq.0) kfafc=4
          vpf=paru(119+2*kfafc)
          apf=paru(120+2*kfafc)
C...Asymmetry and weight.
          asym=2d0*(ei*ai*vint(112)*ef*af+ei*api*vint(113)*ef*apf+
     &    4d0*vi*ai*vint(114)*vf*af+(vi*api+vpi*ai)*vint(115)*
     &    (vf*apf+vpf*af)+4d0*vpi*api*vint(116)*vpf*apf)/
     &    (ei**2*vint(111)*ef**2+ei*vi*vint(112)*ef*vf+
     &    ei*vpi*vint(113)*ef*vpf+(vi**2+ai**2)*vint(114)*
     &    (vf**2+af**2)+(vi*vpi+ai*api)*vint(115)*(vf*vpf+af*apf)+
     &    (vpi**2+api**2)*vint(116)*(vpf**2+apf**2))
          wt=1d0+asym*cthe(1)*isign(1,mint(15)*kfl1(1))+cthe(1)**2
          wtmax=2d0+abs(asym)
        elseif(ip.eq.1.and.iabs(kfl1(1)).eq.24) then
C...Angular weight for f + fbar -> Z' -> W+ + W-.
          rm1=p(nsd(1)+1,5)**2/sh
          rm2=p(nsd(1)+2,5)**2/sh
          ccos2=-(1d0/16d0)*((1d0-rm1-rm2)**2-4d0*rm1*rm2)*
     &    (1d0-2d0*rm1-2d0*rm2+rm1**2+rm2**2+10d0*rm1*rm2)
          cflat=-ccos2+0.5d0*(rm1+rm2)*(1d0-2d0*rm1-2d0*rm2+
     &    (rm2-rm1)**2)
          wt=cflat+ccos2*cthe(1)**2
          wtmax=cflat+max(0d0,ccos2)
        elseif(ip.eq.1.and.(kfl1(1).eq.25.or.kfl1(1).eq.35.or.
     &    iabs(kfl1(1)).eq.37)) then
C...Angular weight for f + fbar -> Z' -> h0 + A0, H0 + A0, H+ + H-.
          wt=1d0-cthe(1)**2
          wtmax=1d0
        elseif(ip.eq.1.and.kfl2(1).eq.25) then
C...Angular weight for f + fbar -> Z' -> Z0 + h0.
          rm1=p(nsd(1)+1,5)**2/sh
          rm2=p(nsd(1)+2,5)**2/sh
          flam2=max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2)
          wt=1d0+flam2*(1d0-cthe(1)**2)/(8d0*rm1)
          wtmax=1d0+flam2/(8d0*rm1)
        elseif(mzpwp.eq.0) then
C...Angular weight for f + fbar -> Z' -> W+ + W- -> 4 quarks/leptons
C...(W:s like if intermediate Z).
          d34=p(iref(ip,iord),5)**2
          d56=p(iref(ip,3-iord),5)**2
          dt=pkk(1,3)+pkk(1,4)+d34
          du=pkk(1,5)+pkk(1,6)+d56
          fgk135=abs(fgk(1,2,3,4,5,6)-fgk(1,2,5,6,3,4))
          fgk253=abs(fgk(2,1,5,6,3,4)-fgk(2,1,3,4,5,6))
          wt=(coup(1,3)*fgk135)**2+(coup(1,4)*fgk253)**2
          wtmax=4d0*d34*d56*(coup(1,3)**2+coup(1,4)**2)*
     &    (digk(dt,du)+digk(du,dt)-djgk(dt,du))
        elseif(mzpwp.eq.1) then
C...Angular weight for f + fbar -> Z' -> W+ + W- -> 4 quarks/leptons
C...(W:s approximately longitudinal, like if intermediate H).
          wt=16d0*pkk(3,5)*pkk(4,6)
          wtmax=sh**2
        else
C...Angular weight for f + fbar -> Z' -> H+ + H-, Z0 + h0, h0 + A0,
C...H0 + A0 -> 4 quarks/leptons.
          wt=1d0
          wtmax=1d0
        endif
 
      elseif(isub.eq.142) then
        if(ip.eq.1.and.iabs(kfl1(1)).lt.20) then
C...Angular weight for f + fbar' -> W'+/- -> 2 quarks/leptons.
          kfai=iabs(mint(15))
          kfaic=1
          if(kfai.gt.10) kfaic=2
          vi=paru(129+2*kfaic)
          ai=paru(130+2*kfaic)
          kfaf=iabs(kfl1(1))
          kfafc=1
          if(kfaf.gt.10) kfafc=2
          vf=paru(129+2*kfafc)
          af=paru(130+2*kfafc)
          asym=8d0*vi*ai*vf*af/((vi**2+ai**2)*(vf**2+af**2))
          wt=1d0+asym*cthe(1)*isign(1,mint(15)*kfl1(1))+cthe(1)**2
          wtmax=2d0+abs(asym)
        elseif(ip.eq.1.and.iabs(kfl2(1)).eq.23) then
C...Angular weight for f + fbar' -> W'+/- -> W+/- + Z0.
          rm1=p(nsd(1)+1,5)**2/sh
          rm2=p(nsd(1)+2,5)**2/sh
          ccos2=-(1d0/16d0)*((1d0-rm1-rm2)**2-4d0*rm1*rm2)*
     &    (1d0-2d0*rm1-2d0*rm2+rm1**2+rm2**2+10d0*rm1*rm2)
          cflat=-ccos2+0.5d0*(rm1+rm2)*(1d0-2d0*rm1-2d0*rm2+
     &    (rm2-rm1)**2)
          wt=cflat+ccos2*cthe(1)**2
          wtmax=cflat+max(0d0,ccos2)
        elseif(ip.eq.1.and.kfl2(1).eq.25) then
C...Angular weight for f + fbar -> W'+/- -> W+/- + h0.
          rm1=p(nsd(1)+1,5)**2/sh
          rm2=p(nsd(1)+2,5)**2/sh
          flam2=max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2)
          wt=1d0+flam2*(1d0-cthe(1)**2)/(8d0*rm1)
          wtmax=1d0+flam2/(8d0*rm1)
        elseif(mzpwp.eq.0) then
C...Angular weight for f + fbar' -> W' -> W + Z0 -> 4 quarks/leptons
C...(W/Z like if intermediate W).
          d34=p(iref(ip,iord),5)**2
          d56=p(iref(ip,3-iord),5)**2
          dt=pkk(1,3)+pkk(1,4)+d34
          du=pkk(1,5)+pkk(1,6)+d56
          fgk135=abs(fgk(1,2,3,4,5,6)-fgk(1,2,5,6,3,4))
          fgk136=abs(fgk(1,2,3,4,6,5)-fgk(1,2,6,5,3,4))
          wt=(coup(5,3)*fgk135)**2+(coup(5,4)*fgk136)**2
          wtmax=4d0*d34*d56*(coup(5,3)**2+coup(5,4)**2)*
     &    (digk(dt,du)+digk(du,dt)-djgk(dt,du))
        elseif(mzpwp.eq.1) then
C...Angular weight for f + fbar' -> W' -> W + Z0 -> 4 quarks/leptons
C...(W/Z approximately longitudinal, like if intermediate H).
          wt=16d0*pkk(3,5)*pkk(4,6)
          wtmax=sh**2
        else
C...Angular weight for f + fbar -> W' -> W + h0 -> whatever.
          wt=1d0
          wtmax=1d0
        endif
 
      elseif(isub.eq.145.or.isub.eq.162.or.isub.eq.163.or.isub.eq.164)
     &  then
C...Isotropic decay of leptoquarks (assumed spin 0).
        wt=1d0
        wtmax=1d0
 
      elseif(isub.eq.147.or.isub.eq.148) then
C...Decays of (spin 1/2) q* -> q + (g,gamma) or (Z0,W+-).
        side=1d0
        if(mint(16).eq.21) side=-1d0
        if(ip.eq.1.and.(kfl1(1).eq.21.or.kfl1(1).eq.22)) then
          wt=1d0+side*cthe(1)
          wtmax=2d0
        elseif(ip.eq.1) then
          rm1=p(nsd(1)+1,5)**2/sh
          wt=1d0+side*cthe(1)*(1d0-0.5d0*rm1)/(1d0+0.5d0*rm1)
          wtmax=1d0+(1d0-0.5d0*rm1)/(1d0+0.5d0*rm1)
        else
C...W/Z decay assumed isotropic, since not known.
          wt=1d0
          wtmax=1d0
        endif
 
      elseif(isub.eq.149) then
C...Isotropic decay of techni-eta.
        wt=1d0
        wtmax=1d0
 
      elseif(isub.eq.191) then
        if(ip.eq.1.and.iabs(kfl1(1)).gt.21) then
C...Angular weight for f + fbar -> rho_tech0 -> W+ W-,
C...W+ pi_tech-, pi_tech+ W- or pi_tech+ pi_tech-.
          wt=1d0-cthe(1)**2
          wtmax=1d0
        elseif(ip.eq.1) then
C...Angular weight for f + fbar -> rho_tech0 -> f fbar.
          cthesg=cthe(1)*isign(1,mint(15))
          xwrht=(1d0-2d0*xw)/(4d0*xw*(1d0-xw))
          bwzr=xwrht*sh*(sh-sqmz)/((sh-sqmz)**2+gmmz**2)
          bwzi=xwrht*sh*gmmz/((sh-sqmz)**2+gmmz**2)
          kfai=iabs(mint(15))
          ei=kchg(kfai,1)/3d0
          ai=sign(1d0,ei+0.1d0)
          vi=ai-4d0*ei*xwv
          vali=0.5d0*(vi+ai)
          vari=0.5d0*(vi-ai)
          alefti=(ei+vali*bwzr)**2+(vali*bwzi)**2
          arighi=(ei+vari*bwzr)**2+(vari*bwzi)**2
          kfaf=iabs(kfl1(1))
          ef=kchg(kfaf,1)/3d0
          af=sign(1d0,ef+0.1d0)
          vf=af-4d0*ef*xwv
          valf=0.5d0*(vf+af)
          varf=0.5d0*(vf-af)
          aleftf=(ef+valf*bwzr)**2+(valf*bwzi)**2
          arighf=(ef+varf*bwzr)**2+(varf*bwzi)**2
          asame=alefti*aleftf+arighi*arighf
          aflip=alefti*arighf+arighi*aleftf
          wt=asame*(1d0+cthesg)**2+aflip*(1d0-cthesg)**2
          wtmax=4d0*max(asame,aflip)
        else
C...Isotropic decay of W/pi_tech produced in rho_tech decay.
          wt=1d0
          wtmax=1d0
        endif
 
      elseif(isub.eq.192) then
        if(ip.eq.1.and.iabs(kfl1(1)).gt.21) then
C...Angular weight for f + fbar' -> rho_tech+ -> W+ Z0,
C...W+ pi_tech0, pi_tech+ Z0 or pi_tech+ pi_tech0.
          wt=1d0-cthe(1)**2
          wtmax=1d0
        elseif(ip.eq.1) then
C...Angular weight for f + fbar' -> rho_tech+ -> f fbar'.
          cthesg=cthe(1)*isign(1,mint(15))
          wt=(1d0+cthesg)**2
          wtmax=4d0
        else
C...Isotropic decay of W/Z/pi_tech produced in rho_tech+ decay.
          wt=1d0
          wtmax=1d0
        endif
 
      elseif(isub.eq.193) then
        if(ip.eq.1.and.iabs(kfl1(1)).gt.21) then
C...Angular weight for f + fbar -> omega_tech0 ->
C...gamma pi_tech0 or Z0 pi_tech0.
          wt=1d0+cthe(1)**2
          wtmax=2d0
        elseif(ip.eq.1) then
C...Angular weight for f + fbar -> omega_tech0 -> f fbar.
          cthesg=cthe(1)*isign(1,mint(15))
          bwzr=(0.5d0/(1d0-xw))*sh*(sh-sqmz)/((sh-sqmz)**2+gmmz**2)
          bwzi=(0.5d0/(1d0-xw))*sh*gmmz/((sh-sqmz)**2+gmmz**2)
          kfai=iabs(mint(15))
          ei=kchg(kfai,1)/3d0
          ai=sign(1d0,ei+0.1d0)
          vi=ai-4d0*ei*xwv
          vali=0.5d0*(vi+ai)
          vari=0.5d0*(vi-ai)
          blefti=(ei-vali*bwzr)**2+(vali*bwzi)**2
          brighi=(ei-vari*bwzr)**2+(vari*bwzi)**2
          kfaf=iabs(kfl1(1))
          ef=kchg(kfaf,1)/3d0
          af=sign(1d0,ef+0.1d0)
          vf=af-4d0*ef*xwv
          valf=0.5d0*(vf+af)
          varf=0.5d0*(vf-af)
          bleftf=(ef-valf*bwzr)**2+(valf*bwzi)**2
          brighf=(ef-varf*bwzr)**2+(varf*bwzi)**2
          bsame=blefti*bleftf+brighi*brighf
          bflip=blefti*brighf+brighi*bleftf
          wt=bsame*(1d0+cthesg)**2+bflip*(1d0-cthesg)**2
          wtmax=4d0*max(bsame,bflip)
        else
C...Isotropic decay of Z/pi_tech produced in omega_tech decay.
          wt=1d0
          wtmax=1d0
        endif
 
C...Obtain correct angular distribution by rejection techniques.
      else
        wt=1d0
        wtmax=1d0
      endif
      if(wt.lt.pjr(0)*wtmax) goto 310
 
C...Construct massive four-vectors using angles chosen.
  470 do 540 jt=1,jtmax
        if(kdcy(jt).eq.0) goto 540
        id=iref(ip,jt)
        do 480 j=1,5
          dpmo(j)=p(id,j)
  480   continue
        dpmo(4)=sqrt(dpmo(1)**2+dpmo(2)**2+dpmo(3)**2+dpmo(5)**2)
CMRENNA++
        if(kfl3(jt).eq.0) then
          call pjrobo(nsd(jt)+1,nsd(jt)+2,acos(cthe(jt)),phi(jt),
     &    dpmo(1)/dpmo(4),dpmo(2)/dpmo(4),dpmo(3)/dpmo(4))
        else
          call pjrobo(nsd(jt)+1,nsd(jt)+3,acos(cthe(jt)),phi(jt),
     &    dpmo(1)/dpmo(4),dpmo(2)/dpmo(4),dpmo(3)/dpmo(4))
        endif
CMRENNA--
 
C...Mark decayed resonances; trace history.
        k(id,1)=k(id,1)+10
        kfa=iabs(k(id,2))
        kca=jamcomp(kfa)
        if(kcqm(jt).ne.0) then
C...Do not kill colour flow through coloured resonance!
        else
          k(id,4)=nsd(jt)+1
          k(id,5)=nsd(jt)+2
          if(kfl3(jt).ne.0) k(id,5)=nsd(jt)+3
        endif
 
C...Add documentation lines.
        if(isub.ne.0) then
          idoc=mint(83)+mint(4)
CMRENNA+++
          ihi=nsd(jt)+2
          if(kfl3(jt).ne.0) ihi=ihi+1
          do 500 i=nsd(jt)+1,ihi
CMRENNA---
            i1=mint(83)+mint(4)+1
            k(i,3)=i1
            if(mstp(128).ge.1) k(i,3)=id
            if(mstp(128).le.1.and.mint(4).lt.mstp(126)) then
              mint(4)=mint(4)+1
              k(i1,1)=21
              k(i1,2)=k(i,2)
              k(i1,3)=iref(ip,jt+3)
              do 490 j=1,5
                p(i1,j)=p(i,j)
  490         continue
            endif
  500     continue
        else
          k(nsd(jt)+1,3)=id
          k(nsd(jt)+2,3)=id
          if(kfl3(jt).ne.0) k(nsd(jt)+3,3)=id 
        endif
 
C...Do showering if any of the two/three products can shower.
        nshbef=n
        if(mstp(71).ge.1) then
          ishow1=0
          kfl1a=iabs(kfl1(jt))
          if(kfl1a.le.22) ishow1=1
          ishow2=0
          kfl2a=iabs(kfl2(jt))
          if(kfl2a.le.22) ishow2=1
          ishow3=0
          if(kfl3(jt).ne.0) then
            kfl3a=iabs(kfl3(jt))
            if(kfl3a.le.22) ishow3=1
          endif
          if(ishow1.eq.0.and.ishow2.eq.0.and.ishow3.eq.0) then
          elseif(kfl3(jt).eq.0) then
            call pjshow(nsd(jt)+1,nsd(jt)+2,p(id,5))
          else
            nsd1=nsd(jt)+1
            nsd2=nsd(jt)+2
            if(ishow1.eq.0.and.ishow3.ne.0) then
              nsd1=nsd(jt)+3
            elseif(ishow2.eq.0.and.ishow3.ne.0) then
              nsd2=nsd(jt)+3
            endif
            pmshow=sqrt(max(0d0,(p(nsd1,4)+p(nsd2,4))**2-
     &      (p(nsd1,1)+p(nsd2,1))**2-(p(nsd1,2)+p(nsd2,2))**2-
     &      (p(nsd1,3)+p(nsd2,3))**2))
            call pjshow(nsd1,nsd2,pmshow)
          endif
        endif
        nshaft=n
        if(jt.eq.1) naft1=n
 
C...Check if decay products moved by shower.
        nsd1=nsd(jt)+1
        nsd2=nsd(jt)+2
        nsd3=nsd(jt)+3
        if(nshaft.gt.nshbef) then
          if(k(nsd1,1).gt.10) then
            do 510 i=nshbef+1,nshaft
              if(k(i,1).lt.10.and.k(i,2).eq.k(nsd1,2)) nsd1=i
  510       continue
          endif
          if(k(nsd2,1).gt.10) then
            do 520 i=nshbef+1,nshaft
              if(k(i,1).lt.10.and.k(i,2).eq.k(nsd2,2).and.
     &        i.ne.nsd1) nsd2=i
  520       continue
          endif
          if(kfl3(jt).ne.0.and.k(nsd3,1).gt.10) then
            do 530 i=nshbef+1,nshaft
              if(k(i,1).lt.10.and.k(i,2).eq.k(nsd3,2).and.
     &        i.ne.nsd1.and.i.ne.nsd2) nsd3=i
  530       continue
          endif
        endif
 
C...Store decay products for further treatment.
        np=np+1
        iref(np,1)=nsd1
        iref(np,2)=nsd2
        iref(np,3)=0
        if(kfl3(jt).ne.0) iref(np,3)=nsd3
        iref(np,4)=idoc+1
        iref(np,5)=idoc+2
        iref(np,6)=0
        if(kfl3(jt).ne.0) iref(np,6)=idoc+3
        iref(np,7)=k(iref(ip,jt),2)
        iref(np,8)=iref(ip,jt)
  540 continue
 
C...Fill information for 2 -> 1 -> 2.
  550 if(jtmax.eq.1.and.kdcy(1).ne.0.and.isub.ne.0) then
        mint(7)=mint(83)+6+2*iset(isub)
        mint(8)=mint(83)+7+2*iset(isub)
        mint(25)=kfl1(1)
        mint(26)=kfl2(1)
        vint(23)=cthe(1)
        rm3=p(n-1,5)**2/sh
        rm4=p(n,5)**2/sh
        be34=sqrt(max(0d0,(1d0-rm3-rm4)**2-4d0*rm3*rm4))
        vint(45)=-0.5d0*sh*(1d0-rm3-rm4-be34*cthe(1))
        vint(46)=-0.5d0*sh*(1d0-rm3-rm4+be34*cthe(1))
        vint(48)=0.25d0*sh*be34**2*max(0d0,1d0-cthe(1)**2)
        vint(47)=sqrt(vint(48))
      endif
 
C...Possibility of colour rearrangement in W+W- events.
      if(isub.eq.25.and.mstp(115).ge.1) then
        iakf1=iabs(kfl1(1))
        iakf2=iabs(kfl1(2))
        iakf3=iabs(kfl2(1))
        iakf4=iabs(kfl2(2))
        if(min(iakf1,iakf2,iakf3,iakf4).ge.1.and.
     &  max(iakf1,iakf2,iakf3,iakf4).le.5) call
     &  pjreco(iref(1,1),iref(1,2),nsd(1),naft1)
      endif
 
C...Loop back if needed.
  560 if(ip.lt.np) goto 130
 
      return
      end
 
C*********************************************************************
 
C...PYMULT
C...Initializes treatment of multiple interactions, selects kinematics
C...of hardest interaction if low-pT physics included in run, and
C...generates all non-hardest interactions.
 
      subroutine pjmult(mmul)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      common/pjint3/xsfx(2,-40:40),isig(1000,3),sigh(1000)
      common/pjint5/ngenpd,ngen(0:500,3),xsec(0:500,3)
      common/pjint7/sigt(0:6,0:6,0:5)
      save /jyjets/,/jydat1/,/jydat2/,/pjsubs/,/pjpars/,/pjint1/,
     &/pjint2/,/pjint3/,/pjint5/,/pjint7/
C...Local arrays and saved variables.
      dimension nmul(20),sigm(20),kstr(500,2),vintsv(80)
      save xt2,xt2fac,xc2,xts,irbin,rbin,nmul,sigm
c...JAM:
      common/jampyda1/qmult(10),imult(2,10)
 
C...Initialization of multiple interaction treatment.
      if(mmul.eq.1) then
        if(mstp(122).ge.1) write(mstu(11),5000) mstp(82)
        isub=96
        mint(1)=96
        vint(63)=0d0
        vint(64)=0d0
        vint(143)=1d0
        vint(144)=1d0
 
C...Loop over phase space points: xT2 choice in 20 bins.
  100   sigsum=0d0
        do 120 ixt2=1,20
          nmul(ixt2)=mstp(83)
          sigm(ixt2)=0d0
          do 110 itry=1,mstp(83)
            rsca=0.05d0*((21-ixt2)-pjr(0))
            xt2=vint(149)*(1d0+vint(149))/(vint(149)+rsca)-vint(149)
            xt2=max(0.01d0*vint(149),xt2)
            vint(25)=xt2
 
C...Choose tau and y*. Calculate cos(theta-hat).
            if(pjr(0).le.coef(isub,1)) then
              taut=(2d0*(1d0+sqrt(1d0-xt2))/xt2-1d0)**pjr(0)
              tau=xt2*(1d0+taut)**2/(4d0*taut)
            else
              tau=xt2*(1d0+tan(pjr(0)*atan(sqrt(1d0/xt2-1d0)))**2)
            endif
            vint(21)=tau
            call pjklim(2)
            ryst=pjr(0)
            myst=1
            if(ryst.gt.coef(isub,8)) myst=2
            if(ryst.gt.coef(isub,8)+coef(isub,9)) myst=3
            call pjkmap(2,myst,pjr(0))
            vint(23)=sqrt(max(0d0,1d0-xt2/tau))*(-1)**int(1.5d0+pjr(0))
 
C...Calculate differential cross-section.
            vint(71)=0.5d0*vint(1)*sqrt(xt2)
            call pjsigh(nchn,sigs)
            sigm(ixt2)=sigm(ixt2)+sigs
  110     continue
          sigsum=sigsum+sigm(ixt2)
  120   continue
        sigsum=sigsum/(20d0*mstp(83))
 
C...Reject result if sigma(parton-parton) is smaller than hadronic one.
        if(sigsum.lt.1.1d0*sigt(0,0,5)) then
          if(mstp(122).ge.1) write(mstu(11),5100) parp(82),sigsum
          parp(82)=0.9d0*parp(82)
          vint(149)=4d0*parp(82)**2/vint(2)
          goto 100
        endif
        if(mstp(122).ge.1) write(mstu(11),5200) parp(82), sigsum
 
C...Start iteration to find k factor.
        yke=sigsum/sigt(0,0,5)
        so=0.5d0
        xi=0d0
        yi=0d0
        xf=0d0
        yf=0d0
        xk=0.5d0
        iit=0
  130   if(iit.eq.0) then
          xk=2d0*xk
        elseif(iit.eq.1) then
          xk=0.5d0*xk
        else
          xk=xi+(yke-yi)*(xf-xi)/(yf-yi)
        endif
 
C...Evaluate overlap integrals.
        if(mstp(82).eq.2) then
          sp=0.5d0*paru(1)*(1d0-exp(-xk))
          sop=sp/paru(1)
        else
          if(mstp(82).eq.3) deltab=0.02d0
          if(mstp(82).eq.4) deltab=min(0.01d0,0.05d0*parp(84))
          sp=0d0
          sop=0d0
          b=-0.5d0*deltab
  140     b=b+deltab
          if(mstp(82).eq.3) then
            ov=exp(-b**2)/paru(2)
          else
            cq2=parp(84)**2
            ov=((1d0-parp(83))**2*exp(-min(50d0,b**2))+
     &      2d0*parp(83)*(1d0-parp(83))*2d0/(1d0+cq2)*
     &      exp(-min(50d0,b**2*2d0/(1d0+cq2)))+
     &      parp(83)**2/cq2*exp(-min(50d0,b**2/cq2)))/paru(2)
          endif
          pacc=1d0-exp(-min(50d0,paru(1)*xk*ov))
          sp=sp+paru(2)*b*deltab*pacc
          sop=sop+paru(2)*b*deltab*ov*pacc
          if(b.lt.1d0.or.b*pacc.gt.1d-6) goto 140
        endif
        yk=paru(1)*xk*so/sp
 
C...Continue iteration until convergence.
        if(yk.lt.yke) then
          xi=xk
          yi=yk
          if(iit.eq.1) iit=2
        else
          xf=xk
          yf=yk
          if(iit.eq.0) iit=1
        endif
        if(abs(yk-yke).ge.1d-5*yke) goto 130
 
C...Store some results for subsequent use.
        vint(145)=sigsum
        vint(146)=sop/so
        vint(147)=sop/sp
 
C...Initialize iteration in xT2 for hardest interaction.
      elseif(mmul.eq.2) then
        if(mstp(82).le.0) then
        elseif(mstp(82).eq.1) then
          xt2=1d0
          xt2fac=xsec(96,1)/sigt(0,0,5)*vint(149)/(1d0-vint(149))
        elseif(mstp(82).eq.2) then
          xt2=1d0
          xt2fac=vint(146)*xsec(96,1)/sigt(0,0,5)*vint(149)*
     &    (1d0+vint(149))
        else
          xc2=4d0*ckin(3)**2/vint(2)
          if(ckin(3).le.ckin(5).or.mint(82).ge.2) xc2=0d0
        endif
 
      elseif(mmul.eq.3) then
C...Low-pT or multiple interactions (first semihard interaction):
C...choose xT2 according to dpT2/pT2**2*exp(-(sigma above pT2)/norm)
C...or (MSTP(82)>=2) dpT2/(pT2+pT0**2)**2*exp(-....).
        isub=mint(1)
        if(mstp(82).le.0) then
          xt2=0d0
        elseif(mstp(82).eq.1) then
          xt2=xt2fac*xt2/(xt2fac-xt2*log(pjr(0)))
        elseif(mstp(82).eq.2) then
          if(xt2.lt.1d0.and.exp(-xt2fac*xt2/(vint(149)*(xt2+
     &    vint(149)))).gt.pjr(0)) xt2=1d0
          if(xt2.ge.1d0) then
            xt2=(1d0+vint(149))*xt2fac/(xt2fac-(1d0+vint(149))*log(1d0-
     &      pjr(0)*(1d0-exp(-xt2fac/(vint(149)*(1d0+vint(149)))))))-
     &      vint(149)
          else
            xt2=-xt2fac/log(exp(-xt2fac/(xt2+vint(149)))+pjr(0)*
     &      (exp(-xt2fac/vint(149))-exp(-xt2fac/(xt2+vint(149)))))-
     &      vint(149)
          endif
          xt2=max(0.01d0*vint(149),xt2)
        else
          xt2=(xc2+vint(149))*(1d0+vint(149))/(1d0+vint(149)-
     &    pjr(0)*(1d0-xc2))-vint(149)
          xt2=max(0.01d0*vint(149),xt2)
        endif
        vint(25)=xt2
 
C...Low-pT: choose xT2, tau, y* and cos(theta-hat) fixed.
        if(mstp(82).le.1.and.xt2.lt.vint(149)) then
          if(mint(82).eq.1) ngen(0,1)=ngen(0,1)-1
          if(mint(82).eq.1) ngen(isub,1)=ngen(isub,1)-1
          isub=95
          mint(1)=isub
          vint(21)=0.01d0*vint(149)
          vint(22)=0d0
          vint(23)=0d0
          vint(25)=0.01d0*vint(149)
 
        else
C...Multiple interactions (first semihard interaction).
C...Choose tau and y*. Calculate cos(theta-hat).
          if(pjr(0).le.coef(isub,1)) then
            taut=(2d0*(1d0+sqrt(1d0-xt2))/xt2-1d0)**pjr(0)
            tau=xt2*(1d0+taut)**2/(4d0*taut)
          else
            tau=xt2*(1d0+tan(pjr(0)*atan(sqrt(1d0/xt2-1d0)))**2)
          endif
          vint(21)=tau
          call pjklim(2)
          ryst=pjr(0)
          myst=1
          if(ryst.gt.coef(isub,8)) myst=2
          if(ryst.gt.coef(isub,8)+coef(isub,9)) myst=3
          call pjkmap(2,myst,pjr(0))
          vint(23)=sqrt(max(0d0,1d0-xt2/tau))*(-1)**int(1.5d0+pjr(0))
        endif
        vint(71)=0.5d0*vint(1)*sqrt(vint(25))
 
C...Store results of cross-section calculation.
      elseif(mmul.eq.4) then
        isub=mint(1)
        xts=vint(25)
        if(iset(isub).eq.1) xts=vint(21)
        if(iset(isub).eq.2)
     &  xts=(4d0*vint(48)+2d0*vint(63)+2d0*vint(64))/vint(2)
        if(iset(isub).ge.3.and.iset(isub).le.5) xts=vint(26)
        rbin=max(0.000001d0,min(0.999999d0,xts*(1d0+vint(149))/
     &  (xts+vint(149))))
        irbin=int(1d0+20d0*rbin)
        if(isub.eq.96.and.mstp(171).eq.0) then
          nmul(irbin)=nmul(irbin)+1
          sigm(irbin)=sigm(irbin)+vint(153)
        endif
 
C...Choose impact parameter.
      elseif(mmul.eq.5) then
        if(mstp(82).eq.3) then
          vint(148)=pjr(0)/(paru(2)*vint(147))
        else
          rtype=pjr(0)
          cq2=parp(84)**2
          if(rtype.lt.(1d0-parp(83))**2) then
            b2=-log(pjr(0))
          elseif(rtype.lt.1d0-parp(83)**2) then
            b2=-0.5d0*(1d0+cq2)*log(pjr(0))
          else
            b2=-cq2*log(pjr(0))
          endif
          vint(148)=((1d0-parp(83))**2*exp(-min(50d0,b2))+2d0*parp(83)*
     &    (1d0-parp(83))*2d0/(1d0+cq2)*exp(-min(50d0,b2*2d0/(1d0+cq2)))+
     &    parp(83)**2/cq2*exp(-min(50d0,b2/cq2)))/(paru(2)*vint(147))
        endif
 
C...Multiple interactions (variable impact parameter) : reject with
C...probability exp(-overlap*cross-section above pT/normalization).
        rncor=(irbin-20d0*rbin)*nmul(irbin)
        sigcor=(irbin-20d0*rbin)*sigm(irbin)
        do 150 ibin=irbin+1,20
          rncor=rncor+nmul(ibin)
          sigcor=sigcor+sigm(ibin)
  150   continue
        sigabv=(sigcor/rncor)*vint(149)*(1d0-xts)/(xts+vint(149))
        if(mstp(171).eq.1) sigabv=sigabv*vint(2)/vint(289)
        vint(150)=exp(-min(50d0,vint(146)*vint(148)*
     &  sigabv/sigt(0,0,5)))
 
C...Generate additional multiple semihard interactions.
      elseif(mmul.eq.6) then
        isubsv=mint(1)
        do 160 j=11,80
          vintsv(j)=vint(j)
  160   continue
        isub=96
        mint(1)=96
 
C...Reconstruct strings in hard scattering.
        nmax=mint(84)+4
        if(iset(isubsv).eq.1) nmax=mint(84)+2
        if(iset(isubsv).eq.11) nmax=mint(84)+2+mint(3)
        nstr=0
        do 180 i=mint(84)+1,nmax
          kcs=kchg(jamcomp(k(i,2)),2)*isign(1,k(i,2))
          if(kcs.eq.0) goto 180
 
          do 170 j=1,4
            if(kcs.eq.1.and.(j.eq.2.or.j.eq.4)) goto 170
            if(kcs.eq.-1.and.(j.eq.1.or.j.eq.3)) goto 170
            if(j.le.2) then
              ist=mod(k(i,j+3)/mstu(5),mstu(5))
            else
              ist=mod(k(i,j+1),mstu(5))
            endif
            if(ist.lt.mint(84).or.ist.gt.i) goto 170
            if(kchg(jamcomp(k(ist,2)),2).eq.0) goto 170
            nstr=nstr+1
            if(j.eq.1.or.j.eq.4) then
              kstr(nstr,1)=i
              kstr(nstr,2)=ist
            else
              kstr(nstr,1)=ist
              kstr(nstr,2)=i
            endif
  170     continue
  180   continue
 
C...Set up starting values for iteration in xT2.
        xt2=vint(25)
        if(iset(isubsv).eq.1) xt2=vint(21)
        if(iset(isubsv).eq.2)
     &  xt2=(4d0*vint(48)+2d0*vint(63)+2d0*vint(64))/vint(2)
        if(iset(isubsv).ge.3.and.iset(isubsv).le.5) xt2=vint(26)
        if(mstp(82).le.1) then
          xt2fac=xsec(isub,1)*vint(149)/((1d0-vint(149))*sigt(0,0,5))
        else
          xt2fac=vint(146)*vint(148)*xsec(isub,1)/sigt(0,0,5)*
     &    vint(149)*(1d0+vint(149))
        endif
        vint(63)=0d0
        vint(64)=0d0
        vint(143)=1d0-vint(141)
        vint(144)=1d0-vint(142)
 
C...Iterate downwards in xT2.
  190   if(mstp(82).le.1) then
          xt2=xt2fac*xt2/(xt2fac-xt2*log(pjr(0)))
          if(xt2.lt.vint(149)) goto 240
        else
          if(xt2.le.0.01001d0*vint(149)) goto 240
          xt2=xt2fac*(xt2+vint(149))/(xt2fac-(xt2+vint(149))*
     &    log(pjr(0)))-vint(149)
          if(xt2.le.0d0) goto 240
          xt2=max(0.01d0*vint(149),xt2)
        endif
        vint(25)=xt2
 
C...Choose tau and y*. Calculate cos(theta-hat).
        if(pjr(0).le.coef(isub,1)) then
          taut=(2d0*(1d0+sqrt(1d0-xt2))/xt2-1d0)**pjr(0)
          tau=xt2*(1d0+taut)**2/(4d0*taut)
        else
          tau=xt2*(1d0+tan(pjr(0)*atan(sqrt(1d0/xt2-1d0)))**2)
        endif
        vint(21)=tau
        call pjklim(2)
        ryst=pjr(0)
        myst=1
        if(ryst.gt.coef(isub,8)) myst=2
        if(ryst.gt.coef(isub,8)+coef(isub,9)) myst=3
        call pjkmap(2,myst,pjr(0))
        vint(23)=sqrt(max(0d0,1d0-xt2/tau))*(-1)**int(1.5d0+pjr(0))
 
C...Check that x not used up. Accept or reject kinematical variables.
        x1m=sqrt(tau)*exp(vint(22))
        x2m=sqrt(tau)*exp(-vint(22))
        if(vint(143)-x1m.lt.0.01d0.or.vint(144)-x2m.lt.0.01d0) goto 190
        vint(71)=0.5d0*vint(1)*sqrt(xt2)
        call pjsigh(nchn,sigs)
        if(sigs.lt.xsec(isub,1)*pjr(0)) goto 190
 
C...Reset K, P and V vectors. Select some variables.
        do 210 i=n+1,n+2
          do 200 j=1,5
            k(i,j)=0
            p(i,j)=0d0
            v(i,j)=0d0
  200     continue
  210   continue
        rflav=pjr(0)
        pt=0.5d0*vint(1)*sqrt(xt2)
        phi=paru(2)*pjr(0)
        cth=vint(23)
 
C...Add first parton to event record.
        k(n+1,1)=3
        k(n+1,2)=21
        if(rflav.ge.max(parp(85),parp(86))) k(n+1,2)=
     &  1+int((2d0+parj(2))*pjr(0))
        p(n+1,1)=pt*cos(phi)
        p(n+1,2)=pt*sin(phi)
        p(n+1,3)=0.25d0*vint(1)*(vint(41)*(1d0+cth)-vint(42)*(1d0-cth))
        p(n+1,4)=0.25d0*vint(1)*(vint(41)*(1d0+cth)+vint(42)*(1d0-cth))
        p(n+1,5)=0d0
 
C...Add second parton to event record.
        k(n+2,1)=3
        k(n+2,2)=21
        if(k(n+1,2).ne.21) k(n+2,2)=-k(n+1,2)
        p(n+2,1)=-p(n+1,1)
        p(n+2,2)=-p(n+1,2)
        p(n+2,3)=0.25d0*vint(1)*(vint(41)*(1d0-cth)-vint(42)*(1d0+cth))
        p(n+2,4)=0.25d0*vint(1)*(vint(41)*(1d0-cth)+vint(42)*(1d0+cth))
        p(n+2,5)=0d0
 
        if(rflav.lt.parp(85).and.nstr.ge.1) then
C....Choose relevant string pieces to place gluons on.
          do 230 i=n+1,n+2
            dmin=1d8
            do 220 istr=1,nstr
              i1=kstr(istr,1)
              i2=kstr(istr,2)
              dist=(p(i,4)*p(i1,4)-p(i,1)*p(i1,1)-p(i,2)*p(i1,2)-
     &        p(i,3)*p(i1,3))*(p(i,4)*p(i2,4)-p(i,1)*p(i2,1)-
     &        p(i,2)*p(i2,2)-p(i,3)*p(i2,3))/max(1d0,p(i1,4)*p(i2,4)-
     &        p(i1,1)*p(i2,1)-p(i1,2)*p(i2,2)-p(i1,3)*p(i2,3))
              if(istr.eq.1.or.dist.lt.dmin) then
                dmin=dist
                ist1=i1
                ist2=i2
                istm=istr
              endif
  220       continue
 
C....Colour flow adjustments, new string pieces.
            if(k(ist1,4)/mstu(5).eq.ist2) k(ist1,4)=mstu(5)*i+
     &      mod(k(ist1,4),mstu(5))
            if(mod(k(ist1,5),mstu(5)).eq.ist2) k(ist1,5)=
     &      mstu(5)*(k(ist1,5)/mstu(5))+i
            k(i,5)=mstu(5)*ist1
            k(i,4)=mstu(5)*ist2
            if(k(ist2,5)/mstu(5).eq.ist1) k(ist2,5)=mstu(5)*i+
     &      mod(k(ist2,5),mstu(5))
            if(mod(k(ist2,4),mstu(5)).eq.ist1) k(ist2,4)=
     &      mstu(5)*(k(ist2,4)/mstu(5))+i
            kstr(istm,2)=i
            kstr(nstr+1,1)=i
            kstr(nstr+1,2)=ist2
            nstr=nstr+1
  230     continue
 
C...String drawing and colour flow for gluon loop.
        elseif(k(n+1,2).eq.21) then
          k(n+1,4)=mstu(5)*(n+2)
          k(n+1,5)=mstu(5)*(n+2)
          k(n+2,4)=mstu(5)*(n+1)
          k(n+2,5)=mstu(5)*(n+1)
          kstr(nstr+1,1)=n+1
          kstr(nstr+1,2)=n+2
          kstr(nstr+2,1)=n+2
          kstr(nstr+2,2)=n+1
          nstr=nstr+2
 
C...String drawing and colour flow for qqbar pair.
        else
          k(n+1,4)=mstu(5)*(n+2)
          k(n+2,5)=mstu(5)*(n+1)
          kstr(nstr+1,1)=n+1
          kstr(nstr+1,2)=n+2
          nstr=nstr+1
        endif

c...jam:save multiple scattered partons.
        imult(1,mint(31))=n+1
        imult(2,mint(31))=n+2
        qmult(mint(31))=pt
 
C...Update remaining energy; iterate.
        n=n+2
        if(n.gt.mstu(4)-mstu(32)-10) then
          call pjerrm(11,'(PYMULT:) no more memory left in PYJETS')
          if(mstu(21).ge.1) return
        endif
        mint(31)=mint(31)+1
        vint(151)=vint(151)+vint(41)
        vint(152)=vint(152)+vint(42)
        vint(143)=vint(143)-vint(41)
        vint(144)=vint(144)-vint(42)
        if(mint(31).lt.240) goto 190
  240   continue
        mint(1)=isubsv
        do 250 j=11,80
          vint(j)=vintsv(j)
  250   continue
      endif
 
C...Format statements for printout.
 5000 format(/1x,'****** PYMULT: initialization of multiple inter',
     &'actions for MSTP(82) =',i2,' ******')
 5100 format(8x,'pT0 =',f5.2,' GeV gives sigma(parton-parton) =',1p,
     &d9.2,' mb: rejected')
 5200 format(8x,'pT0 =',f5.2,' GeV gives sigma(parton-parton) =',1p,
     &d9.2,' mb: accepted')
 
      return
      end
 
C*********************************************************************
 
C...PYREMN
C...Adds on target remnants (one or two from each side) and
C...includes primordial kT for hadron beams.
 
      subroutine pjremn(ipu1,ipu2)
 
C...Double precision declarations.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
cp    common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
cp    common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      save /jyjets/,/pjpars/,/pjint1/
C...Local arrays.
      dimension kflch(2),kflsp(2),chi(2),pms(0:6),is(2),isn(2),robo(5),
     &psys(0:2,5),pmin(0:2),qold(4),qnew(4),dbe(3),psum(4)

      real*8 jamrnd2
 
C...Find event type and remaining energy.
      isub=mint(1)
      ns=n
      if(mint(50).eq.0.or.mstp(81).le.0) then
        vint(143)=1d0-vint(141)
        vint(144)=1d0-vint(142)
      endif
 
C...Define initial partons.
      ntry=0
  100 ntry=ntry+1
      do 130 jt=1,2
        i=mint(83)+jt+2
        if(jt.eq.1) ipu=ipu1
        if(jt.eq.2) ipu=ipu2
        k(i,1)=21
        k(i,2)=k(ipu,2)
        k(i,3)=i-2
        pms(jt)=0d0
        vint(156+jt)=0d0
        vint(158+jt)=0d0
        if(mint(47).eq.1) then
          do 110 j=1,5
            p(i,j)=p(i-2,j)
  110     continue
        elseif(isub.eq.95) then
          k(i,2)=21
        else
          p(i,5)=p(ipu,5)
 
C...No primordial kT, or chosen according to truncated Gaussian or
C...exponential, or (for photon) predetermined or power law.
  120     if(mint(40+jt).eq.2.and.mint(10+jt).ne.22) then
            if(mstp(91).le.0) then
              pt=0d0
            elseif(mstp(91).eq.1) then
              pt=parp(91)*sqrt(-log(pjr(0)))
            elseif(mstp(91).eq.2) then
              rpt1=pjr(0)
              rpt2=pjr(0)
              pt=-parp(92)*log(rpt1*rpt2)
            else
            endif
            if(pt.gt.parp(93)) goto 120
          elseif(mint(106+jt).eq.3) then
            pt=sqrt(vint(282+jt))
            pt=pt*0.8d0**mint(57)
            if(ntry.gt.10) pt=pt*0.8d0**(ntry-10)
          elseif(iabs(mint(14+jt)).le.8.or.mint(14+jt).eq.21) then
            if(mstp(93).le.0) then
              pt=0d0
            elseif(mstp(93).eq.1) then
              pt=parp(99)*sqrt(-log(pjr(0)))
            elseif(mstp(93).eq.2) then
              rpt1=pjr(0)
              rpt2=pjr(0)
              pt=-parp(99)*log(rpt1*rpt2)
            elseif(mstp(93).eq.3) then
              ha=parp(99)**2
              hb=parp(100)**2
              pt=sqrt(max(0d0,ha*(ha+hb)/(ha+hb-pjr(0)*hb)-ha))
            else
              ha=parp(99)**2
              hb=parp(100)**2
              if(mstp(93).eq.5) hb=min(vint(48),parp(100)**2)
              pt=sqrt(max(0d0,ha*((ha+hb)/ha)**pjr(0)-ha))
            endif
            if(pt.gt.parp(100)) goto 120
          else
            pt=0d0
          endif
          vint(156+jt)=pt
          phi=paru(2)*pjr(0)
          p(i,1)=pt*cos(phi)
          p(i,2)=pt*sin(phi)
          pms(jt)=p(i,5)**2+p(i,1)**2+p(i,2)**2
        endif
  130 continue
      if(mint(47).eq.1) return
 
C...Kinematics construction for initial partons.
      i1=mint(83)+3
      i2=mint(83)+4
      if(isub.eq.95) then
        shs=0d0
        shr=0d0
      else
        shs=vint(141)*vint(142)*vint(2)+(p(i1,1)+p(i2,1))**2+
     &  (p(i1,2)+p(i2,2))**2
        shr=sqrt(max(0d0,shs))
        if((shs-pms(1)-pms(2))**2-4d0*pms(1)*pms(2).le.0d0) goto 100
        p(i1,4)=0.5d0*(shr+(pms(1)-pms(2))/shr)
        p(i1,3)=sqrt(max(0d0,p(i1,4)**2-pms(1)))
        p(i2,4)=shr-p(i1,4)
        p(i2,3)=-p(i1,3)
 
C...Transform partons to overall CM-frame.
        robo(3)=(p(i1,1)+p(i2,1))/shr
        robo(4)=(p(i1,2)+p(i2,2))/shr
        call pjrobo(i1,i2,0d0,0d0,-robo(3),-robo(4),0d0)
        robo(2)=pjangl(p(i1,1),p(i1,2))
        call pjrobo(i1,i2,0d0,-robo(2),0d0,0d0,0d0)
        robo(1)=pjangl(p(i1,3),p(i1,1))
        call pjrobo(i1,i2,-robo(1),0d0,0d0,0d0,0d0)
        call pjrobo(i1,mint(52),robo(1),robo(2),robo(3),robo(4),0d0)
        robo(5)=max(-0.999999d0,min(0.999999d0,(vint(141)-vint(142))/
     &  (vint(141)+vint(142))))
        call pjrobo(i1,mint(52),0d0,0d0,0d0,0d0,robo(5))
      endif
 
C...Optionally fix up x and Q2 definitions for leptoproduction.
      idisxq=0
      if((mint(43).eq.2.or.mint(43).eq.3).and.((isub.eq.10.and.
     &mstp(23).ge.1).or.(isub.eq.83.and.mstp(23).ge.2))) idisxq=1
      if(idisxq.eq.1) then
 
C...Find where incoming and outgoing leptons/partons are sitting.
        lesd=1
        if(mint(42).eq.1) lesd=2
        lpin=mint(83)+3-lesd
        lein=mint(84)+lesd
        lqin=mint(84)+3-lesd
        leout=mint(84)+2+lesd
        lqout=mint(84)+5-lesd
        if(k(lein,3).gt.lein) lein=k(lein,3)
        if(k(lqin,3).gt.lqin) lqin=k(lqin,3)
        lscms=0
        do 140 i=mint(84)+5,n
          if(k(i,2).eq.94) then
            lscms=i
            leout=i+lesd
            lqout=i+3-lesd
          endif
  140   continue
        lqbg=ipu1
        if(lesd.eq.1) lqbg=ipu2
 
C...Calculate actual and wanted momentum transfer.
        xnom=vint(43-lesd)
        q2nom=-vint(45)
        hpk=2d0*(p(lpin,4)*p(lein,4)-p(lpin,1)*p(lein,1)-
     &  p(lpin,2)*p(lein,2)-p(lpin,3)*p(lein,3))*
     &  (p(mint(83)+lesd,4)*vint(40+lesd)/p(lein,4))
        hpt2=max(0d0,q2nom*(1d0-q2nom/(xnom*hpk)))
        fac=sqrt(hpt2/(p(leout,1)**2+p(leout,2)**2))
        p(n+1,1)=fac*p(leout,1)
        p(n+1,2)=fac*p(leout,2)
        p(n+1,3)=0.25d0*((hpk-q2nom/xnom)/p(lpin,4)-
     &  q2nom/(p(mint(83)+lesd,4)*vint(40+lesd)))*(-1)**(lesd+1)
        p(n+1,4)=sqrt(p(leout,5)**2+p(n+1,1)**2+p(n+1,2)**2+
     &  p(n+1,3)**2)
        do 150 j=1,4
          qold(j)=p(lein,j)-p(leout,j)
          qnew(j)=p(lein,j)-p(n+1,j)
  150   continue
 
C...Boost outgoing electron and daughters.
        if(lscms.eq.0) then
          do 160 j=1,4
            p(leout,j)=p(n+1,j)
  160     continue
        else
          do 170 j=1,3
            p(n+2,j)=(p(n+1,j)-p(leout,j))/(p(n+1,4)+p(leout,4))
  170     continue
          pinv=2d0/(1d0+p(n+2,1)**2+p(n+2,2)**2+p(n+2,3)**2)
          do 180 j=1,3
            dbe(j)=pinv*p(n+2,j)
  180     continue
          do 200 i=lscms+1,n
            iorig=i
  190       iorig=k(iorig,3)
            if(iorig.gt.leout) goto 190
            if(i.eq.leout.or.iorig.eq.leout)
     &      call pjrobo(i,i,0d0,0d0,dbe(1),dbe(2),dbe(3))
  200     continue
        endif
 
C...Copy shower initiator and all outgoing partons.
        ncop=n+1
        k(ncop,3)=lqbg
        do 210 j=1,5
          p(ncop,j)=p(lqbg,j)
  210   continue
        do 240 i=mint(84)+1,n
          icop=0
          if(k(i,1).gt.10) goto 240
          if(i.eq.lqbg.or.i.eq.lqout) then
            icop=i
          else
            iorig=i
  220       iorig=k(iorig,3)
            if(iorig.eq.lqbg.or.iorig.eq.lqout) then
              icop=iorig
            elseif(iorig.gt.mint(84).and.iorig.le.n) then
              goto 220
            endif
          endif
          if(icop.ne.0) then
            ncop=ncop+1
            k(ncop,3)=i
            do 230 j=1,5
              p(ncop,j)=p(i,j)
  230       continue
          endif
  240   continue
 
C...Calculate relative rescaling factors.
        slc=3-2*lesd
        plcsum=0d0
        do 250 i=n+2,ncop
          plcsum=plcsum+(p(i,4)+slc*p(i,3))
  250   continue
        do 260 i=n+2,ncop
          v(i,1)=(p(i,4)+slc*p(i,3))/plcsum
  260   continue
 
C...Transfer extra three-momentum of current.
        do 280 i=n+2,ncop
          do 270 j=1,3
            p(i,j)=p(i,j)+v(i,1)*(qnew(j)-qold(j))
  270     continue
          p(i,4)=sqrt(p(i,5)**2+p(i,1)**2+p(i,2)**2+p(i,3)**2)
  280   continue
 
C...Iterate change of initiator momentum to get energy right.
        iter=0
  290   iter=iter+1
        peex=-p(n+1,4)-qnew(4)
        pemv=-p(n+1,3)/p(n+1,4)
        do 300 i=n+2,ncop
          peex=peex+p(i,4)
          pemv=pemv+v(i,1)*p(i,3)/p(i,4)
  300   continue
        if(abs(pemv).lt.1d-10) then
          mint(51)=1
          mint(57)=mint(57)+1
          return
        endif
        pzch=-peex/pemv
        p(n+1,3)=p(n+1,3)+pzch
        p(n+1,4)=sqrt(p(n+1,5)**2+p(n+1,1)**2+p(n+1,2)**2+p(n+1,3)**2)
        do 310 i=n+2,ncop
          p(i,3)=p(i,3)+v(i,1)*pzch
          p(i,4)=sqrt(p(i,5)**2+p(i,1)**2+p(i,2)**2+p(i,3)**2)
  310   continue
        if(iter.lt.10.and.abs(peex).gt.1d-6*p(n+1,4)) goto 290
 
C...Modify momenta in event record.
        hbe=2d0*(p(n+1,4)+p(lqbg,4))*(p(n+1,3)-p(lqbg,3))/
     &  ((p(n+1,4)+p(lqbg,4))**2+(p(n+1,3)-p(lqbg,3))**2)
        if(abs(hbe).gt.0.999999d0) then
          mint(51)=1
          mint(57)=mint(57)+1
          return
        endif
        i=mint(83)+5-lesd
        call pjrobo(i,i,0d0,0d0,0d0,0d0,hbe)
        do 330 i=n+1,ncop
          icop=k(i,3)
          do 320 j=1,4
            p(icop,j)=p(i,j)
  320     continue
  330   continue
      endif
 
C...Check minimum invariant mass of remnant system(s).
      psys(0,4)=p(i1,4)+p(i2,4)+0.5d0*vint(1)*(vint(151)+vint(152))
      psys(0,3)=p(i1,3)+p(i2,3)+0.5d0*vint(1)*(vint(151)-vint(152))
      pms(0)=max(0d0,psys(0,4)**2-psys(0,3)**2)
      pmin(0)=sqrt(pms(0))
      do 340 jt=1,2
        psys(jt,4)=0.5d0*vint(1)*vint(142+jt)
        psys(jt,3)=psys(jt,4)*(-1)**(jt-1)
        pmin(jt)=0d0
        if(mint(44+jt).eq.1) goto 340
        mint(105)=mint(102+jt)
        mint(109)=mint(106+jt)
        call pjspli(mint(10+jt),mint(12+jt),kflch(jt),kflsp(jt))
        if(kflch(jt).ne.0) pmin(jt)=pmin(jt)+pjmass(kflch(jt))
        if(kflsp(jt).ne.0) pmin(jt)=pmin(jt)+pjmass(kflsp(jt))
        if(kflch(jt)*kflsp(jt).ne.0) pmin(jt)=pmin(jt)+0.5d0*parp(111)
        pmin(jt)=sqrt(pmin(jt)**2+p(mint(83)+jt+2,1)**2+
     &  p(mint(83)+jt+2,2)**2)
  340 continue
      if(pmin(0)+pmin(1)+pmin(2).gt.vint(1).or.(mint(45).ge.2.and.
     &pmin(1).gt.psys(1,4)).or.(mint(46).ge.2.and.pmin(2).gt.
     &psys(2,4))) then
        mint(51)=1
        mint(57)=mint(57)+1
        return
      endif
 
C...Loop over two remnants; skip if none there.
      i=ns
      do 410 jt=1,2
        isn(jt)=0
        if(mint(44+jt).eq.1) goto 410
        if(jt.eq.1) ipu=ipu1
        if(jt.eq.2) ipu=ipu2
 
C...Store first remnant parton.
        i=i+1
        is(jt)=i
        isn(jt)=1
        do 350 j=1,5
          k(i,j)=0
          p(i,j)=0d0
          v(i,j)=0d0
  350   continue
        k(i,1)=1
        k(i,2)=kflsp(jt)
        k(i,3)=mint(83)+jt
        p(i,5)=pjmass(k(i,2))
 
C...First parton colour connections and kinematics.
        kcol=kchg(jamcomp(kflsp(jt)),2)
        if(kcol.eq.2) then
          k(i,1)=3
          k(i,4)=mstu(5)*ipu+ipu
          k(i,5)=mstu(5)*ipu+ipu
          k(ipu,4)=mod(k(ipu,4),mstu(5))+mstu(5)*i
          k(ipu,5)=mod(k(ipu,5),mstu(5))+mstu(5)*i
        elseif(kcol.ne.0) then
          k(i,1)=3
          kfls=(3-kcol*isign(1,kflsp(jt)))/2
          k(i,kfls+3)=ipu
          k(ipu,6-kfls)=mod(k(ipu,6-kfls),mstu(5))+mstu(5)*i
        endif
        if(kflch(jt).eq.0) then
          p(i,1)=-p(mint(83)+jt+2,1)
          p(i,2)=-p(mint(83)+jt+2,2)
          pms(jt)=p(i,5)**2+p(i,1)**2+p(i,2)**2
          psys(jt,3)=sqrt(max(0d0,psys(jt,4)**2-pms(jt)))*(-1)**(jt-1)
          p(i,3)=psys(jt,3)
          p(i,4)=psys(jt,4)
 
C...When extra remnant parton or hadron: store extra remnant.
        else
          i=i+1
          isn(jt)=2
          do 360 j=1,5
            k(i,j)=0
            p(i,j)=0d0
            v(i,j)=0d0
  360     continue
          k(i,1)=1
          k(i,2)=kflch(jt)
          k(i,3)=mint(83)+jt
          p(i,5)=pjmass(k(i,2))
 
C...Find parton colour connections of extra remnant.
          kcol=kchg(jamcomp(kflch(jt)),2)
          if(kcol.eq.2) then
            k(i,1)=3
            k(i,4)=mstu(5)*ipu+ipu
            k(i,5)=mstu(5)*ipu+ipu
            k(ipu,4)=mod(k(ipu,4),mstu(5))+mstu(5)*i
            k(ipu,5)=mod(k(ipu,5),mstu(5))+mstu(5)*i
          elseif(kcol.ne.0) then
            k(i,1)=3
            kfls=(3-kcol*isign(1,kflch(jt)))/2
            k(i,kfls+3)=ipu
            k(ipu,6-kfls)=mod(k(ipu,6-kfls),mstu(5))+mstu(5)*i
          endif
 
C...Relative transverse momentum when two remnants.
          loop=0
  370     loop=loop+1
c+JAM
          iii=0
          if(iii.eq.0) then
            call pjptdi(1,p(i-1,1),p(i-1,2))
          else if(iii.eq.1) then
            pt=0.5d0*sqrt(-log(max(1d-10,pjr(0))))
            pt=paru(2)*pjr(0)
            p(i-1,1)=pt*cos(phix)
            p(i-1,2)=pt*sin(phix)

          else

c.....HIJING parametrization.
            pkcmx=sqrt(vint(2)/4d0)
c           pkc=jamrnd2(3,0.0d0,pkcmx)
            pkc=sqrt(hirnd2(3,0.0d0,pkcmx**2))
            if(pkc.gt.parc(68))
     &        pkc=parc(66)*sqrt(-log(exp(-parc(68)**2/parc(66)**2)
     &         -rn(0)*(exp(-parc(68)**2/parc(66)**2)-
     &         exp(-pkcmx**2/parc(66)**2))))
            phix=paru(2)*pjr(0)
            p(i-1,1)=pkc*cos(phix)
            p(i-1,2)=pkc*sin(phix)
            print *,'pkc2',pkc
 
          endif
c-JAM
          if(iabs(mint(10+jt)).lt.20) then
            p(i-1,1)=0d0
            p(i-1,2)=0d0
          endif
          pms(jt+2)=p(i-1,5)**2+p(i-1,1)**2+p(i-1,2)**2
          p(i,1)=-p(mint(83)+jt+2,1)-p(i-1,1)
          p(i,2)=-p(mint(83)+jt+2,2)-p(i-1,2)
          pms(jt+4)=p(i,5)**2+p(i,1)**2+p(i,2)**2
 
C...Meson or baryon; photon as meson. For splitup below.
          imb=1
          if(mod(mint(10+jt)/1000,10).ne.0) imb=2
 
C***Relative distribution for electron into two electrons. Temporary!
          if(iabs(mint(10+jt)).lt.20.and.mint(14+jt).eq.-mint(10+jt))
     &    then
            chi(jt)=pjr(0)
 
C...Relative distribution of electron energy into electron plus parton.
          elseif(iabs(mint(10+jt)).lt.20) then
            xhrd=vint(140+jt)
            xe=vint(154+jt)
            chi(jt)=(xe-xhrd)/(1d0-xhrd)
 
C...Relative distribution of energy for particle into two jets.
          elseif(iabs(kflch(jt)).le.10.or.kflch(jt).eq.21) then
            chik=parp(92+2*imb)
            if(mstp(92).le.1) then
              if(imb.eq.1) chi(jt)=pjr(0)
              if(imb.eq.2) chi(jt)=1d0-sqrt(pjr(0))
            elseif(mstp(92).eq.2) then
              chi(jt)=1d0-pjr(0)**(1d0/(1d0+chik))
            elseif(mstp(92).eq.3) then
              cut=2d0*0.3d0/vint(1)
  380         chi(jt)=pjr(0)**2
              if((chi(jt)**2/(chi(jt)**2+cut**2))**0.25d0*
     &        (1d0-chi(jt))**chik.lt.pjr(0)) goto 380
            elseif(mstp(92).eq.4) then
              cut=2d0*0.3d0/vint(1)
              cutr=(1d0+sqrt(1d0+cut**2))/cut
  390         chir=cut*cutr**pjr(0)
              chi(jt)=(chir**2-cut**2)/(2d0*chir)
              if((1d0-chi(jt))**chik.lt.pjr(0)) goto 390
            else
              cut=2d0*0.3d0/vint(1)
              cuta=cut**(1d0-parp(98))
              cutb=(1d0+cut)**(1d0-parp(98))
cp400         chi(jt)=(cuta+pjr(0)*(cutb-cuta))**(1d0/(1d0-parp(98)))
  400         chi(jt)=(cuta+pjr(0)*(cutb-cuta))**(1d0/(1d0-parp(98)))
     $        -cut
              if(((chi(jt)+cut)**2/(2d0*(chi(jt)**2+cut**2)))**
     &        (0.5d0*parp(98))*(1d0-chi(jt))**chik.lt.pjr(0)) goto 400
            endif
 
C...Relative distribution of energy for particle into jet plus particle.
          else
            if(mstp(94).le.1) then
              if(imb.eq.1) chi(jt)=pjr(0)
              if(imb.eq.2) chi(jt)=1d0-sqrt(pjr(0))
              if(mod(kflch(jt)/1000,10).ne.0) chi(jt)=1d0-chi(jt)
            elseif(mstp(94).eq.2) then
              chi(jt)=1d0-pjr(0)**(1d0/(1d0+parp(93+2*imb)))
              if(mod(kflch(jt)/1000,10).ne.0) chi(jt)=1d0-chi(jt)
            elseif(mstp(94).eq.3) then
              call pjzdis(1,0,pms(jt+4),zz)
              chi(jt)=zz
            else
              call pjzdis(1000,0,pms(jt+4),zz)
              chi(jt)=zz
            endif
          endif
 
C...Construct total transverse mass; reject if too large.
          pms(jt)=pms(jt+4)/chi(jt)+pms(jt+2)/(1d0-chi(jt))
          if(pms(jt).gt.psys(jt,4)**2) then
            if(loop.lt.10) then
              goto 370
            else
              mint(51)=1
              mint(57)=mint(57)+1
              return
            endif
          endif
          psys(jt,3)=sqrt(max(0d0,psys(jt,4)**2-pms(jt)))*(-1)**(jt-1)
          vint(158+jt)=chi(jt)
 
C...Subdivide longitudinal momentum according to value selected above.
          pw1=chi(jt)*(psys(jt,4)+abs(psys(jt,3)))
          p(is(jt)+1,4)=0.5d0*(pw1+pms(jt+4)/pw1)
          p(is(jt)+1,3)=0.5d0*(pw1-pms(jt+4)/pw1)*(-1)**(jt-1)
          p(is(jt),4)=psys(jt,4)-p(is(jt)+1,4)
          p(is(jt),3)=psys(jt,3)-p(is(jt)+1,3)
        endif
  410 continue
      n=i
 
C...Check if longitudinal boosts needed - if so pick two systems.
      pdev=abs(psys(0,4)+psys(1,4)+psys(2,4)-vint(1))+
     &abs(psys(0,3)+psys(1,3)+psys(2,3))
      if(pdev.le.1d-6*vint(1)) return
      if(isn(1).eq.0) then
        ir=0
        il=2
      elseif(isn(2).eq.0) then
        ir=1
        il=0
      elseif(vint(143).gt.0.2d0.and.vint(144).gt.0.2d0) then
        ir=1
        il=2
      elseif(vint(143).gt.0.2d0) then
        ir=1
        il=0
      elseif(vint(144).gt.0.2d0) then
        ir=0
        il=2
      elseif(pms(1)/psys(1,4)**2.gt.pms(2)/psys(2,4)**2) then
        ir=1
        il=0
      else
        ir=0
        il=2
      endif
      ig=3-ir-il
 
C...E+-pL wanted for system to be modified.
      if((ig.eq.1.and.isn(1).eq.0).or.(ig.eq.2.and.isn(2).eq.0)) then
        ppb=vint(1)
        pnb=vint(1)
      else
        ppb=vint(1)-(psys(ig,4)+psys(ig,3))
        pnb=vint(1)-(psys(ig,4)-psys(ig,3))
      endif
 
C...To keep x and Q2 in leptoproduction: do not count scattered lepton.
      if(idisxq.eq.1.and.ig.ne.0) then
        pmtb=ppb*pnb
        pmtr=pms(ir)
        pmtl=pms(il)
        sqlam=sqrt(max(0d0,(pmtb-pmtr-pmtl)**2-4d0*pmtr*pmtl))
        sqsgn=sign(1d0,psys(ir,3)*psys(il,4)-psys(il,3)*psys(ir,4))
        rkr=(pmtb+pmtr-pmtl+sqlam*sqsgn)/(2d0*(psys(ir,4)+psys(ir,3))
     &  *pnb)
        rkl=(pmtb+pmtl-pmtr+sqlam*sqsgn)/(2d0*(psys(il,4)-psys(il,3))
     &  *ppb)
        ber=(rkr**2-1d0)/(rkr**2+1d0)
        bel=-(rkl**2-1d0)/(rkl**2+1d0)
        ppb=ppb-(psys(0,4)+psys(0,3))
        pnb=pnb-(psys(0,4)-psys(0,3))
        do 420 j=1,4
          psys(0,j)=0d0
  420   continue
        do 450 i=mint(84)+1,ns
          if(k(i,1).gt.10) goto 450
          incl=0
          iorig=i
  430     if(iorig.eq.lqout.or.iorig.eq.lpin+2) incl=1
          iorig=k(iorig,3)
          if(iorig.gt.lpin) goto 430
          if(incl.eq.0) goto 450
          do 440 j=1,4
            psys(0,j)=psys(0,j)+p(i,j)
  440     continue
  450   continue
        pms(0)=max(0d0,psys(0,4)**2-psys(0,3)**2)
        ppb=ppb+(psys(0,4)+psys(0,3))
        pnb=pnb+(psys(0,4)-psys(0,3))
      endif
 
C...Construct longitudinal boosts.
      dpmtb=ppb*pnb
      dpmtr=pms(ir)
      dpmtl=pms(il)
      dsqlam=sqrt(max(0d0,(dpmtb-dpmtr-dpmtl)**2-4d0*dpmtr*dpmtl))
      if(dsqlam.le.1d-6*dpmtb) then
        mint(51)=1
        mint(57)=mint(57)+1
        return
      endif
      dsqsgn=sign(1d0,psys(ir,3)*psys(il,4)-psys(il,3)*psys(ir,4))
      drkr=(dpmtb+dpmtr-dpmtl+dsqlam*dsqsgn)/
     &(2d0*(psys(ir,4)+psys(ir,3))*pnb)
      drkl=(dpmtb+dpmtl-dpmtr+dsqlam*dsqsgn)/
     &(2d0*(psys(il,4)-psys(il,3))*ppb)
      dber=(drkr**2-1d0)/(drkr**2+1d0)
      dbel=-(drkl**2-1d0)/(drkl**2+1d0)
 
C...Perform longitudinal boosts.
      if(ir.eq.1.and.isn(1).eq.1.and.dber.le.-0.99999999d0) then
        p(is(1),3)=0d0
        p(is(1),4)=sqrt(p(is(1),5)**2+p(is(1),1)**2+p(is(1),2)**2)
      elseif(ir.eq.1) then
        call pjrobo(is(1),is(1)+isn(1)-1,0d0,0d0,0d0,0d0,dber)
      elseif(idisxq.eq.1) then
        do 470 i=i1,ns
          incl=0
          iorig=i
  460     if(iorig.eq.lqout.or.iorig.eq.lpin+2) incl=1
          iorig=k(iorig,3)
          if(iorig.gt.lpin) goto 460
          if(incl.eq.1) call pjrobo(i,i,0d0,0d0,0d0,0d0,dber)
  470   continue
      else
        call pjrobo(i1,ns,0d0,0d0,0d0,0d0,dber)
      endif
      if(il.eq.2.and.isn(2).eq.1.and.dbel.ge.0.99999999d0) then
        p(is(2),3)=0d0
        p(is(2),4)=sqrt(p(is(2),5)**2+p(is(2),1)**2+p(is(2),2)**2)
      elseif(il.eq.2) then
        call pjrobo(is(2),is(2)+isn(2)-1,0d0,0d0,0d0,0d0,dbel)
      elseif(idisxq.eq.1) then
        do 490 i=i1,ns
          incl=0
          iorig=i
  480     if(iorig.eq.lqout.or.iorig.eq.lpin+2) incl=1
          iorig=k(iorig,3)
          if(iorig.gt.lpin) goto 480
          if(incl.eq.1) call pjrobo(i,i,0d0,0d0,0d0,0d0,dbel)
  490   continue
      else
        call pjrobo(i1,ns,0d0,0d0,0d0,0d0,dbel)
      endif
 
C...Final check that energy-momentum conservation worked.
      pesum=0d0
      pzsum=0d0
      do 500 i=mint(84)+1,n
        if(k(i,1).gt.10) goto 500
        pesum=pesum+p(i,4)
        pzsum=pzsum+p(i,3)
  500 continue
      pdev=abs(pesum-vint(1))+abs(pzsum)
      if(pdev.gt.1d-4*vint(1)) then
        mint(51)=1
        mint(57)=mint(57)+1
        return
      endif
 
C...Calculate rotation and boost from overall CM frame to
C...hadronic CM frame in leptoproduction.
      mint(91)=0
      if(mint(82).eq.1.and.(mint(43).eq.2.or.mint(43).eq.3)) then
        mint(91)=1
        lesd=1
        if(mint(42).eq.1) lesd=2
        lpin=mint(83)+3-lesd
 
C...Sum upp momenta of everything not lepton or photon to define boost.
        do 510 j=1,4
          psum(j)=0d0
  510   continue
        do 530 i=1,n
          if(k(i,1).le.0.or.k(i,1).gt.10) goto 530
          if(iabs(k(i,2)).ge.11.and.iabs(k(i,2)).le.20) goto 530
          if(k(i,2).eq.22) goto 530
          do 520 j=1,4
            psum(j)=psum(j)+p(i,j)
  520     continue
  530   continue
        vint(223)=-psum(1)/psum(4)
        vint(224)=-psum(2)/psum(4)
        vint(225)=-psum(3)/psum(4)
 
C...Boost incoming hadron to hadronic CM frame to determine rotations.
        k(n+1,1)=1
        do 540 j=1,5
          p(n+1,j)=p(lpin,j)
          v(n+1,j)=v(lpin,j)
  540   continue
        call pjrobo(n+1,n+1,0d0,0d0,vint(223),vint(224),vint(225))
        vint(222)=-pjangl(p(n+1,1),p(n+1,2))
        call pjrobo(n+1,n+1,0d0,vint(222),0d0,0d0,0d0)
        if(lesd.eq.2) then
          vint(221)=-pjangl(p(n+1,3),p(n+1,1))
        else
          vint(221)=pjangl(-p(n+1,3),p(n+1,1))
        endif
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYDIFF
C...Handles diffractive and elastic scattering.
 
      subroutine pjdiff
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      save /jyjets/,/jydat1/,/pjpars/,/pjint1/
 
C...Reset K, P and V vectors. Store incoming particles.
      do 110 jt=1,mstp(126)+10
        i=mint(83)+jt
        do 100 j=1,5
          k(i,j)=0
          p(i,j)=0d0
          v(i,j)=0d0
  100   continue
  110 continue
      n=mint(84)
      mint(3)=0
      mint(21)=0
      mint(22)=0
      mint(23)=0
      mint(24)=0
      mint(4)=4
      do 130 jt=1,2
        i=mint(83)+jt
        k(i,1)=21
        k(i,2)=mint(10+jt)
        do 120 j=1,5
          p(i,j)=vint(285+5*jt+j)
  120   continue
  130 continue
      mint(6)=2
 
C...Subprocess; kinematics.
      sqlam=(vint(2)-vint(63)-vint(64))**2-4d0*vint(63)*vint(64)
      pz=sqrt(sqlam)/(2d0*vint(1))
      do 200 jt=1,2
        i=mint(83)+jt
        pe=(vint(2)+vint(62+jt)-vint(65-jt))/(2d0*vint(1))
        kfh=mint(102+jt)
 
C...Elastically scattered particle.
        if(mint(16+jt).le.0) then
          n=n+1
          k(n,1)=1
          k(n,2)=kfh
          k(n,3)=i+2
          p(n,3)=pz*(-1)**(jt+1)
          p(n,4)=pe
          p(n,5)=sqrt(vint(62+jt))
 
C...Decay rho from elastic scattering of gamma with sin**2(theta)
C...distribution of decay products (in rho rest frame).
          if(kfh.eq.113.and.mint(10+jt).eq.22.and.mstp(102).eq.1) then
            nsav=n
            dbetaz=p(n,3)/sqrt(p(n,3)**2+p(n,5)**2)
            p(n,3)=0d0
            p(n,4)=p(n,5)
            call pjdecy(nsav,icon)
            if(n.eq.nsav+2.and.iabs(k(nsav+1,2)).eq.211) then
              phi=pjangl(p(nsav+1,1),p(nsav+1,2))
              call pjrobo(nsav+1,nsav+2,0d0,-phi,0d0,0d0,0d0)
              the=pjangl(p(nsav+1,3),p(nsav+1,1))
              call pjrobo(nsav+1,nsav+2,-the,0d0,0d0,0d0,0d0)
  140         cthe=2d0*pjr(0)-1d0
              if(1d0-cthe**2.lt.pjr(0)) goto 140
              call pjrobo(nsav+1,nsav+2,acos(cthe),phi,0d0,0d0,0d0)
            endif
            call pjrobo(nsav,nsav+2,0d0,0d0,0d0,0d0,dbetaz)
          endif
 
C...Diffracted particle: low-mass system to two particles.
        elseif(vint(62+jt).lt.(vint(66+jt)+parp(103))**2) then
          n=n+2
          k(n-1,1)=1
          k(n,1)=1
          k(n-1,3)=i+2
          k(n,3)=i+2
          pmmas=sqrt(vint(62+jt))
          ntry=0
  150     ntry=ntry+1
          if(ntry.lt.20) then
            mint(105)=mint(102+jt)
            mint(109)=mint(106+jt)
            call pjspli(kfh,21,kfl1,kfl2)
            call pjkfdi(kfl1,0,kfl3,kf1)
            if(kf1.eq.0) goto 150
            call pjkfdi(kfl2,-kfl3,kfldum,kf2)
            if(kf2.eq.0) goto 150
          else
            kf1=kfh
            kf2=111
          endif
          pm1=pjmass(kf1)
          pm2=pjmass(kf2)
          if(pm1+pm2+parj(64).gt.pmmas) goto 150
          k(n-1,2)=kf1
          k(n,2)=kf2
          p(n-1,5)=pm1
          p(n,5)=pm2
          pzp=sqrt(max(0d0,(pmmas**2-pm1**2-pm2**2)**2-
     &    4d0*pm1**2*pm2**2))/(2d0*pmmas)
          p(n-1,3)=pzp
          p(n,3)=-pzp
          p(n-1,4)=sqrt(pm1**2+pzp**2)
          p(n,4)=sqrt(pm2**2+pzp**2)
          call pjrobo(n-1,n,acos(2d0*pjr(0)-1d0),paru(2)*pjr(0),
     &    0d0,0d0,0d0)
          dbetaz=pz*(-1)**(jt+1)/sqrt(pz**2+pmmas**2)
          call pjrobo(n-1,n,0d0,0d0,0d0,0d0,dbetaz)
 
C...Diffracted particle: valence quark kicked out.
        elseif(mstp(101).eq.1.or.(mstp(101).eq.3.and.pjr(0).lt.
     &    parp(101))) then
          n=n+2
          k(n-1,1)=2
          k(n,1)=1
          k(n-1,3)=i+2
          k(n,3)=i+2
          mint(105)=mint(102+jt)
          mint(109)=mint(106+jt)
          call pjspli(kfh,21,k(n,2),k(n-1,2))
          p(n-1,5)=pjmass(k(n-1,2))
          p(n,5)=pjmass(k(n,2))
          sqlam=(vint(62+jt)-p(n-1,5)**2-p(n,5)**2)**2-
     &    4d0*p(n-1,5)**2*p(n,5)**2
          p(n-1,3)=(pe*sqrt(sqlam)+pz*(vint(62+jt)+p(n-1,5)**2-
     &    p(n,5)**2))/(2d0*vint(62+jt))*(-1)**(jt+1)
          p(n-1,4)=sqrt(p(n-1,3)**2+p(n-1,5)**2)
          p(n,3)=pz*(-1)**(jt+1)-p(n-1,3)
          p(n,4)=sqrt(p(n,3)**2+p(n,5)**2)
 
C...Diffracted particle: gluon kicked out.
        else
          n=n+3
          k(n-2,1)=2
          k(n-1,1)=2
          k(n,1)=1
          k(n-2,3)=i+2
          k(n-1,3)=i+2
          k(n,3)=i+2
          mint(105)=mint(102+jt)
          mint(109)=mint(106+jt)
          call pjspli(kfh,21,k(n,2),k(n-2,2))
          k(n-1,2)=21
          p(n-2,5)=pjmass(k(n-2,2))
          p(n-1,5)=0d0
          p(n,5)=pjmass(k(n,2))
C...Energy distribution for particle into two jets.
  160     imb=1
          if(mod(kfh/1000,10).ne.0) imb=2
          chik=parp(92+2*imb)
          if(mstp(92).le.1) then
            if(imb.eq.1) chi=pjr(0)
            if(imb.eq.2) chi=1d0-sqrt(pjr(0))
          elseif(mstp(92).eq.2) then
            chi=1d0-pjr(0)**(1d0/(1d0+chik))
          elseif(mstp(92).eq.3) then
            cut=2d0*0.3d0/vint(1)
  170       chi=pjr(0)**2
            if((chi**2/(chi**2+cut**2))**0.25d0*(1d0-chi)**chik.lt.
     &      pjr(0)) goto 170
          elseif(mstp(92).eq.4) then
            cut=2d0*0.3d0/vint(1)
            cutr=(1d0+sqrt(1d0+cut**2))/cut
  180       chir=cut*cutr**pjr(0)
            chi=(chir**2-cut**2)/(2d0*chir)
            if((1d0-chi)**chik.lt.pjr(0)) goto 180
          else
            cut=2d0*0.3d0/vint(1)
            cuta=cut**(1d0-parp(98))
            cutb=(1d0+cut)**(1d0-parp(98))
  190       chi=(cuta+pjr(0)*(cutb-cuta))**(1d0/(1d0-parp(98)))
            if(((chi+cut)**2/(2d0*(chi**2+cut**2)))**
     &      (0.5d0*parp(98))*(1d0-chi)**chik.lt.pjr(0)) goto 190
          endif
          if(chi.lt.p(n,5)**2/vint(62+jt).or.chi.gt.1d0-p(n-2,5)**2/
     &    vint(62+jt)) goto 160
          sqm=p(n-2,5)**2/(1d0-chi)+p(n,5)**2/chi
          if((sqrt(sqm)+parj(32))**2.ge.vint(62+jt)) goto 160
          pzi=(pe*(vint(62+jt)-sqm)+pz*(vint(62+jt)+sqm))/
     &    (2d0*vint(62+jt))
          pei=sqrt(pzi**2+sqm)
          pqqp=(1d0-chi)*(pei+pzi)
          p(n-2,3)=0.5d0*(pqqp-p(n-2,5)**2/pqqp)*(-1)**(jt+1)
          p(n-2,4)=sqrt(p(n-2,3)**2+p(n-2,5)**2)
          p(n-1,4)=0.5d0*(vint(62+jt)-sqm)/(pei+pzi)
          p(n-1,3)=p(n-1,4)*(-1)**jt
          p(n,3)=pzi*(-1)**(jt+1)-p(n-2,3)
          p(n,4)=sqrt(p(n,3)**2+p(n,5)**2)
        endif
 
C...Documentation lines.
        k(i+2,1)=21
        if(mint(16+jt).eq.0) k(i+2,2)=kfh
        if(mint(16+jt).ne.0) k(i+2,2)=10*(kfh/10)
        k(i+2,3)=i
        p(i+2,3)=pz*(-1)**(jt+1)
        p(i+2,4)=pe
        p(i+2,5)=sqrt(vint(62+jt))
  200 continue
 
C...Rotate outgoing partons/particles using cos(theta).
      if(vint(23).lt.0.9d0) then
        call pjrobo(mint(83)+3,n,acos(vint(23)),vint(24),0d0,0d0,0d0)
      else
        call pjrobo(mint(83)+3,n,asin(vint(59)),vint(24),0d0,0d0,0d0)
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYDOCU
C...Handles the documentation of the process in MSTI and PARI,
C...and also computes cross-sections based on accumulated statistics.
 
      subroutine pjdocu
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      common/pjint5/ngenpd,ngen(0:500,3),xsec(0:500,3)
      save /jyjets/,/jydat1/,/pjsubs/,/pjpars/,/pjint1/,/pjint2/,
     &/pjint5/
 
C...Calculate Monte Carlo estimates of cross-sections.
      isub=mint(1)
      if(mstp(111).ne.-1) ngen(isub,3)=ngen(isub,3)+1
      ngen(0,3)=ngen(0,3)+1
      xsec(0,3)=0d0
      do 100 i=1,500
        if(i.eq.96.or.i.eq.97) then
          xsec(i,3)=0d0
        elseif(msub(95).eq.1.and.(i.eq.11.or.i.eq.12.or.i.eq.13.or.
     &    i.eq.28.or.i.eq.53.or.i.eq.68)) then
          xsec(i,3)=xsec(96,2)*ngen(i,3)/max(1d0,dble(ngen(96,1))*
     &    dble(ngen(96,2)))
        elseif(msub(i).eq.0.or.ngen(i,1).eq.0) then
          xsec(i,3)=0d0
        elseif(ngen(i,2).eq.0) then
          xsec(i,3)=xsec(i,2)*ngen(0,3)/(dble(ngen(i,1))*
     &    dble(ngen(0,2)))
        else
          xsec(i,3)=xsec(i,2)*ngen(i,3)/(dble(ngen(i,1))*
     &    dble(ngen(i,2)))
        endif
        xsec(0,3)=xsec(0,3)+xsec(i,3)
  100 continue
 
C...Rescale to known low-pT cross-section for standard QCD processes.
      if(msub(95).eq.1) then
        xsech=xsec(11,3)+xsec(12,3)+xsec(13,3)+xsec(28,3)+xsec(53,3)+
     &  xsec(68,3)+xsec(95,3)
        xsecw=xsec(97,2)/max(1d0,dble(ngen(97,1)))
        if(xsech.gt.1d-10.and.xsecw.gt.1d-10) then
          fac=xsecw/xsech
          xsec(11,3)=fac*xsec(11,3)
          xsec(12,3)=fac*xsec(12,3)
          xsec(13,3)=fac*xsec(13,3)
          xsec(28,3)=fac*xsec(28,3)
          xsec(53,3)=fac*xsec(53,3)
          xsec(68,3)=fac*xsec(68,3)
          xsec(95,3)=fac*xsec(95,3)
          xsec(0,3)=xsec(0,3)-xsech+xsecw
        endif
      endif
 
C...Save information for gamma-p and gamma-gamma.
      if(mint(121).gt.1) then
        iga=mint(122)
        call pjsave(2,iga)
        call pjsave(5,0)
      endif
 
C...Reset information on hard interaction.
      do 110 j=1,200
        msti(j)=0
        pari(j)=0d0
  110 continue
 
C...Copy integer valued information from MINT into MSTI.
      do 120 j=1,32
        msti(j)=mint(j)
  120 continue
      if(mint(121).gt.1) msti(9)=mint(122)
 
C...Store cross-section variables in PARI.
      pari(1)=xsec(0,3)
      pari(2)=xsec(0,3)/mint(5)
      pari(9)=vint(99)
      pari(10)=vint(100)
      vint(98)=vint(98)+vint(100)
      if(mstp(142).eq.1) pari(2)=xsec(0,3)/vint(98)
 
C...Store kinematics variables in PARI.
      pari(11)=vint(1)
      pari(12)=vint(2)
      if(isub.ne.95) then
        do 130 j=13,26
          pari(j)=vint(30+j)
  130   continue
        pari(31)=vint(141)
        pari(32)=vint(142)
        pari(33)=vint(41)
        pari(34)=vint(42)
        pari(35)=pari(33)-pari(34)
        pari(36)=vint(21)
        pari(37)=vint(22)
        pari(38)=vint(26)
        pari(39)=vint(157)
        pari(40)=vint(158)
        pari(41)=vint(23)
        pari(42)=2d0*vint(47)/vint(1)
      endif
 
C...Store information on scattered partons in PARI.
      if(isub.ne.95.and.mint(7)*mint(8).ne.0) then
        do 140 is=7,8
          i=mint(is)
          pari(36+is)=p(i,3)/vint(1)
          pari(38+is)=p(i,4)/vint(1)
          pr=max(1d-20,p(i,5)**2+p(i,1)**2+p(i,2)**2)
          pari(40+is)=sign(log(min((sqrt(pr+p(i,3)**2)+abs(p(i,3)))/
     &    sqrt(pr),1d20)),p(i,3))
          pr=max(1d-20,p(i,1)**2+p(i,2)**2)
          pari(42+is)=sign(log(min((sqrt(pr+p(i,3)**2)+abs(p(i,3)))/
     &    sqrt(pr),1d20)),p(i,3))
          pari(44+is)=p(i,3)/sqrt(1d-20+p(i,1)**2+p(i,2)**2+p(i,3)**2)
          pari(46+is)=pjangl(p(i,3),sqrt(p(i,1)**2+p(i,2)**2))
          pari(48+is)=pjangl(p(i,1),p(i,2))
  140   continue
      endif
 
C...Store sum up transverse and longitudinal momenta.
      pari(65)=2d0*pari(17)
      if(isub.le.90.or.isub.ge.95) then
        do 150 i=mstp(126)+1,n
          if(k(i,1).le.0.or.k(i,1).gt.10) goto 150
          pt=sqrt(p(i,1)**2+p(i,2)**2)
          pari(69)=pari(69)+pt
          if(i.le.mint(52)) pari(66)=pari(66)+pt
          if(i.gt.mint(52).and.i.le.mint(53)) pari(68)=pari(68)+pt
  150   continue
        pari(67)=pari(68)
        pari(71)=vint(151)
        pari(72)=vint(152)
        pari(73)=vint(151)
        pari(74)=vint(152)
      else
        pari(66)=pari(65)
        pari(69)=pari(65)
      endif
 
C...Store various other pieces of information into PARI.
      pari(61)=vint(148)
      pari(75)=vint(155)
      pari(76)=vint(156)
      pari(77)=vint(159)
      pari(78)=vint(160)
      pari(81)=vint(138)
 
C...Set information for PYTABU.
      if(iset(isub).eq.1.or.iset(isub).eq.3) then
        mstu(161)=mint(21)
        mstu(162)=0
      elseif(iset(isub).eq.5) then
        mstu(161)=mint(23)
        mstu(162)=0
      else
        mstu(161)=mint(21)
        mstu(162)=mint(22)
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYFRAM
C...Performs transformations between different coordinate frames.
 
      subroutine pjfram(iframe)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      save /jydat1/,/pjpars/,/pjint1/
 
C...Check that transformation can and should be done.
      if(iframe.eq.1.or.iframe.eq.2.or.(iframe.eq.3.and.
     &mint(91).eq.1)) then
        if(iframe.eq.mint(6)) return
      else
        write(mstu(11),5000) iframe,mint(6)
        return
      endif
 
      if(mint(6).eq.1) then
C...Transform from fixed target or user specified frame to
C...overall CM frame.
        call pjrobo(0,0,0d0,0d0,-vint(8),-vint(9),-vint(10))
        call pjrobo(0,0,0d0,-vint(7),0d0,0d0,0d0)
        call pjrobo(0,0,-vint(6),0d0,0d0,0d0,0d0)
      elseif(mint(6).eq.3) then
C...Transform from hadronic CM frame in DIS to overall CM frame.
        call pjrobo(0,0,-vint(221),-vint(222),-vint(223),-vint(224),
     &  -vint(225))
      endif
 
      if(iframe.eq.1) then
C...Transform from overall CM frame to fixed target or user specified
C...frame.
        call pjrobo(0,0,vint(6),vint(7),vint(8),vint(9),vint(10))
      elseif(iframe.eq.3) then
C...Transform from overall CM frame to hadronic CM frame in DIS.
        call pjrobo(0,0,0d0,0d0,vint(223),vint(224),vint(225))
        call pjrobo(0,0,0d0,vint(222),0d0,0d0,0d0)
        call pjrobo(0,0,vint(221),0d0,0d0,0d0,0d0)
      endif
 
C...Set information about new frame.
      mint(6)=iframe
      msti(6)=iframe
 
 5000 format(1x,'Error: illegal values in subroutine PYFRAM.',1x,
     &'No transformation performed.'/1x,'IFRAME =',1x,i5,'; MINT(6) =',
     &1x,i5)
 
      return
      end
 
C*********************************************************************
 
C...PYWIDT
C...Calculates full and partial widths of resonances.
 
      subroutine pjwidt(kflr,sh,wdtp,wdte)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Parameter statement to help give large particle numbers.
      parameter (ksusy1=1000000,ksusy2=2000000,kexcit=4000000)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/jydat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint4/mwid(500),wids(500,5)
      common/pjmssm/imss(0:99),rmss(0:99)
      common/pjssmt/zmix(4,4),umix(2,2),vmix(2,2),smz(4),smw(2),
     &sfmix(16,4)
      save /jydat1/,/jydat2/,/jydat3/,/pjsubs/,/pjpars/,/pjint1/,
     &/pjint4/,/pjmssm/,/pjssmt/
C...Local arrays and saved variables.
      dimension wdtp(0:200),wdte(0:200,0:5),mofsv(3,2),widwsv(3,2),
     &wid2sv(3,2)
      save mofsv,widwsv,wid2sv
      data mofsv/6*0/,widwsv/6*0d0/,wid2sv/6*0d0/
 
C...Compressed code and sign; mass.
      kfla=iabs(kflr)
      kfls=isign(1,kflr)
      kc=jamcomp(kfla)
      shr=sqrt(sh)
      pmr=pmas(kc,1)
 
C...Reset width information.
      do 110 i=0,200
        wdtp(i)=0d0
        do 100 j=0,5
          wdte(i,j)=0d0
  100   continue
  110 continue
 
C...Not to be treated as a resonance: return.
      if((mwid(kc).le.0.or.mwid(kc).ge.4).and.kfla.ne.21.and.
     &kfla.ne.22) then
        wdtp(0)=1d0
        wdte(0,0)=1d0
        mint(61)=0
        mint(62)=0
        mint(63)=0
        return
 
C...Treatment as a resonance based on tabulated branching ratios.
      elseif(mwid(kc).eq.2.or.(mwid(kc).eq.3.and.mint(63).eq.0)) then
C...Loop over possible decay channels; skip irrelevant ones.
        do 120 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 120
 
C...Read out decay products and nominal masses.
          kfd1=kfdp(idc,1)
          kfc1=jamcomp(kfd1)
          if(kchg(kfc1,3).eq.1) kfd1=kfls*kfd1
          pm1=pmas(kfc1,1)
          kfd2=kfdp(idc,2)
          kfc2=jamcomp(kfd2)
          if(kchg(kfc2,3).eq.1) kfd2=kfls*kfd2
          pm2=pmas(kfc2,1)
          kfd3=kfdp(idc,3)
          pm3=0d0
          if(kfd3.ne.0) then
            kfc3=jamcomp(kfd3)
            if(kchg(kfc3,3).eq.1) kfd3=kfls*kfd3
            pm3=pmas(kfc3,1)
          endif
 
C...Naive partial width and alternative threshold factors.
          wdtp(i)=pmas(kc,2)*brat(idc)*(shr/pmr)
          if(mdme(idc,2).ge.51.and.mdme(idc,2).le.53.and.
     &    pm1+pm2+pm3.ge.shr) then
             wdtp(i)=0d0
          elseif(mdme(idc,2).eq.52.and.kfd3.eq.0) then
            wdtp(i)=wdtp(i)*sqrt(max(0d0,(sh-pm1**2-pm2**2)**2-
     &      4d0*pm1**2*pm2**2))/sh
          elseif(mdme(idc,2).eq.52) then
            pma=max(pm1,pm2,pm3)
            pmc=min(pm1,pm2,pm3)
            pmb=pm1+pm2+pm3-pma-pmc
            pmbc=pmb+pmc+0.5d0*(shr-pma-pmc-pmc)
            pman=pma**2/sh
            pmbn=pmb**2/sh
            pmcn=pmc**2/sh
            pmbcn=pmbc**2/sh
            wdtp(i)=wdtp(i)*sqrt(max(0d0,
     &      ((1d0-pman-pmbcn)**2-4d0*pman*pmbcn)*
     &      ((pmbcn-pmbn-pmcn)**2-4d0*pmbn*pmcn)))*
     &      ((shr-pma)**2-(pmb+pmc)**2)*
     &      (1d0+0.25d0*(pma+pmb+pmc)/shr)/
     &      ((1d0-pmbcn)*pmbcn*sh)
          elseif(mdme(idc,2).eq.53.and.kfd3.eq.0) then
            wdtp(i)=wdtp(i)*sqrt(
     &      max(0d0,(sh-pm1**2-pm2**2)**2-4d0*pm1**2*pm2**2)/
     &      max(1d-4,(pmr**2-pm1**2-pm2**2)**2-4d0*pm1**2*pm2**2))
          elseif(mdme(idc,2).eq.53) then
            pma=max(pm1,pm2,pm3)
            pmc=min(pm1,pm2,pm3)
            pmb=pm1+pm2+pm3-pma-pmc
            pmbc=pmb+pmc+0.5d0*(shr-pma-pmb-pmc)
            pman=pma**2/sh
            pmbn=pmb**2/sh
            pmcn=pmc**2/sh
            pmbcn=pmbc**2/sh
            facact=sqrt(max(0d0,
     &      ((1d0-pman-pmbcn)**2-4d0*pman*pmbcn)*
     &      ((pmbcn-pmbn-pmcn)**2-4d0*pmbn*pmcn)))*
     &      ((shr-pma)**2-(pmb+pmc)**2)*
     &      (1d0+0.25d0*(pma+pmb+pmc)/shr)/
     &      ((1d0-pmbcn)*pmbcn*sh)
            pmbc=pmb+pmc+0.5d0*(pmr-pma-pmb-pmc)
            pman=pma**2/pmr**2
            pmbn=pmb**2/pmr**2
            pmcn=pmc**2/pmr**2
            pmbcn=pmbc**2/pmr**2
            facnom=sqrt(max(0d0,
     &      ((1d0-pman-pmbcn)**2-4d0*pman*pmbcn)*
     &      ((pmbcn-pmbn-pmcn)**2-4d0*pmbn*pmcn)))*
     &      ((pmr-pma)**2-(pmb+pmc)**2)*
     &      (1d0+0.25d0*(pma+pmb+pmc)/pmr)/
     &      ((1d0-pmbcn)*pmbcn*pmr**2)
            wdtp(i)=wdtp(i)*facact/max(1d-6,facnom)
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
 
C...Calculate secondary width (at most two identical/opposite).
          if(mdme(idc,1).gt.0) then
            if(kfd2.eq.kfd1) then
              if(kchg(kfc1,3).eq.0) then
                wid2=wids(kfc1,1)
              elseif(kfd1.gt.0) then
                wid2=wids(kfc1,4)
              else
                wid2=wids(kfc1,5)
              endif
              if(kfd3.gt.0) then
                wid2=wid2*wids(kfc3,2)
              elseif(kfd3.lt.0) then
                wid2=wid2*wids(kfc3,3)
              endif
            elseif(kfd2.eq.-kfd1) then
              wid2=wids(kfc1,1)
              if(kfd3.gt.0) then
                wid2=wid2*wids(kfc3,2)
              elseif(kfd3.lt.0) then
                wid2=wid2*wids(kfc3,3)
              endif
            elseif(kfd3.eq.kfd1) then
              if(kchg(kfc1,3).eq.0) then
                wid2=wids(kfc1,1)
              elseif(kfd1.gt.0) then
                wid2=wids(kfc1,4)
              else
                wid2=wids(kfc1,5)
              endif
              if(kfd2.gt.0) then
                wid2=wid2*wids(kfc2,2)
              elseif(kfd2.lt.0) then
                wid2=wid2*wids(kfc2,3)
              endif
            elseif(kfd3.eq.-kfd1) then
              wid2=wids(kfc1,1)
              if(kfd2.gt.0) then
                wid2=wid2*wids(kfc2,2)
              elseif(kfd2.lt.0) then
                wid2=wid2*wids(kfc2,3)
              endif
            elseif(kfd3.eq.kfd2) then
              if(kchg(kfc2,3).eq.0) then
                wid2=wids(kfc2,1)
              elseif(kfd2.gt.0) then
                wid2=wids(kfc2,4)
              else
                wid2=wids(kfc2,5)
              endif
              if(kfd1.gt.0) then
                wid2=wid2*wids(kfc1,2)
              elseif(kfd1.lt.0) then
                wid2=wid2*wids(kfc1,3)
              endif
            elseif(kfd3.eq.-kfd2) then
              wid2=wids(kfc2,1)
              if(kfd1.gt.0) then
                wid2=wid2*wids(kfc1,2)
              elseif(kfd1.lt.0) then
                wid2=wid2*wids(kfc1,3)
              endif
            else
              if(kfd1.gt.0) then
                wid2=wids(kfc1,2)
              else
                wid2=wids(kfc1,3)
              endif
              if(kfd2.gt.0) then
                wid2=wid2*wids(kfc2,2)
              else
                wid2=wid2*wids(kfc2,3)
              endif
              if(kfd3.gt.0) then
                wid2=wid2*wids(kfc3,2)
              elseif(kfd3.lt.0) then
                wid2=wid2*wids(kfc3,3)
              endif
            endif
 
C...Store effective widths according to case.
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  120   continue
C...Return.
        mint(61)=0
        mint(62)=0
        mint(63)=0
        return
      endif
 
C...Here begins detailed dynamical calculation of resonance widths.
C...Shared treatment of Higgs states.
      kfhigg=25
      ihigg=1
      if(kfla.eq.35.or.kfla.eq.36) then
        kfhigg=kfla
        ihigg=kfla-33
      endif
 
C...Common electroweak and strong constants.
      xw=paru(102)
      xwv=xw
      if(mstp(8).ge.2) xw=1d0-(pmas(24,1)/pmas(23,1))**2
      xw1=1d0-xw
      aem=pjalem(sh)
      if(mstp(8).ge.1) aem=sqrt(2d0)*paru(105)*pmas(24,1)**2*xw/paru(1)
      as=pjalps(sh)
      radc=1d0+as/paru(1)
 
      if(kfla.eq.6) then
C...t quark.
        fac=(aem/(16d0*xw))*(sh/pmas(24,1)**2)*shr
        radct=1d0-2.5d0*as/paru(1)
        do 130 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 130
          rm1=pmas(jamcomp(kfdp(idc,1)),1)**2/sh
          rm2=pmas(jamcomp(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 130
          if(i.ge.4.and.i.le.7) then
C...t -> W + q; including approximate QCD correction factor.
            wdtp(i)=fac*vckm(3,i-3)*radct*
     &      sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))*
     &      ((1d0-rm2)**2+(1d0+rm2)*rm1-2d0*rm1**2)
            if(kflr.gt.0) then
              wid2=wids(24,2)
              if(i.eq.7) wid2=wid2*wids(7,2)
            else
              wid2=wids(24,3)
              if(i.eq.7) wid2=wid2*wids(7,3)
            endif
          elseif(i.eq.9) then
C...t -> H + b.
            wdtp(i)=fac*sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))*
     &      ((1d0+rm2-rm1)*(rm2*paru(141)**2+1d0/paru(141)**2)+4d0*rm2)
            wid2=wids(37,2)
            if(kflr.lt.0) wid2=wids(37,3)
CMRENNA++
          elseif(i.ge.10.and.i.le.13.and.imss(1).ne.0) then
C...t -> ~t + ~chi_i0, i = 1, 2, 3 or 4.
            beta=atan(rmss(5))
            sinb=sin(beta)
            tanw=sqrt(paru(102)/(1d0-paru(102)))
            et=kchg(6,1)/3d0
            t3l=sign(0.5d0,et)
            kfc1=jamcomp(kfdp(idc,1))
            kfc2=jamcomp(kfdp(idc,2))
            pmnchi=pmas(kfc1,1)
            pmstop=pmas(kfc2,1)
            if(shr.gt.pmnchi+pmstop) then
              iz=i-9
              al=shr*zmix(iz,4)/(2.0d0*pmas(24,1)*sinb)
              ar=-et*zmix(iz,1)*tanw
              bl=t3l*(zmix(iz,2)-zmix(iz,1)*tanw)-ar
              br=al
              fl=sfmix(6,1)*al+sfmix(6,2)*ar
              fr=sfmix(6,1)*bl+sfmix(6,2)*br
              pcm=sqrt((sh-(pmnchi+pmstop)**2)*
     &        (sh-(pmnchi-pmstop)**2))/(2d0*shr)
              wdtp(i)=(0.5d0*pjalem(sh)/paru(102))*pcm*((fl**2+fr**2)*
     &        (sh+pmnchi**2-pmstop**2)+smz(iz)*4d0*shr*fl*fr)/sh
              if(kflr.gt.0) then
                wid2=wids(kfc1,2)*wids(kfc2,2)
              else
                wid2=wids(kfc1,2)*wids(kfc2,3)
              endif
            endif
CMRENNA--
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  130   continue
 
      elseif(kfla.eq.7) then
C...b' quark.
        fac=(aem/(16d0*xw))*(sh/pmas(24,1)**2)*shr
        do 140 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 140
          rm1=pmas(jamcomp(kfdp(idc,1)),1)**2/sh
          rm2=pmas(jamcomp(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 140
          if(i.ge.4.and.i.le.7) then
C...b' -> W + q.
            wdtp(i)=fac*vckm(i-3,4)*
     &      sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))*
     &      ((1d0-rm2)**2+(1d0+rm2)*rm1-2d0*rm1**2)
            if(kflr.gt.0) then
              wid2=wids(24,3)
              if(i.eq.6) wid2=wid2*wids(6,2)
              if(i.eq.7) wid2=wid2*wids(8,2)
            else
              wid2=wids(24,2)
              if(i.eq.6) wid2=wid2*wids(6,3)
              if(i.eq.7) wid2=wid2*wids(8,3)
            endif
            wid2=wids(24,3)
            if(kflr.lt.0) wid2=wids(24,2)
          elseif(i.eq.9.or.i.eq.10) then
C...b' -> H + q.
            wdtp(i)=fac*sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))*
     &      ((1d0+rm2-rm1)*(paru(141)**2+rm2/paru(141)**2)+4d0*rm2)
            if(kflr.gt.0) then
              wid2=wids(37,3)
              if(i.eq.10) wid2=wid2*wids(6,2)
            else
              wid2=wids(37,2)
              if(i.eq.10) wid2=wid2*wids(6,3)
            endif
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  140   continue
 
      elseif(kfla.eq.8) then
C...t' quark.
        fac=(aem/(16d0*xw))*(sh/pmas(24,1)**2)*shr
        do 150 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 150
          rm1=pmas(jamcomp(kfdp(idc,1)),1)**2/sh
          rm2=pmas(jamcomp(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 150
          if(i.ge.4.and.i.le.7) then
C...t' -> W + q.
            wdtp(i)=fac*vckm(4,i-3)*
     &      sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))*
     &      ((1d0-rm2)**2+(1d0+rm2)*rm1-2d0*rm1**2)
            if(kflr.gt.0) then
              wid2=wids(24,2)
              if(i.eq.7) wid2=wid2*wids(7,2)
            else
              wid2=wids(24,3)
              if(i.eq.7) wid2=wid2*wids(7,3)
            endif
          elseif(i.eq.9.or.i.eq.10) then
C...t' -> H + q.
            wdtp(i)=fac*sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))*
     &      ((1d0+rm2-rm1)*(rm2*paru(141)**2+1d0/paru(141)**2)+4d0*rm2)
            if(kflr.gt.0) then
              wid2=wids(37,2)
              if(i.eq.10) wid2=wid2*wids(7,2)
            else
              wid2=wids(37,3)
              if(i.eq.10) wid2=wid2*wids(7,3)
            endif
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  150   continue
 
      elseif(kfla.eq.17) then
C...tau' lepton.
        fac=(aem/(16d0*xw))*(sh/pmas(24,1)**2)*shr
        do 160 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 160
          rm1=pmas(jamcomp(kfdp(idc,1)),1)**2/sh
          rm2=pmas(jamcomp(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 160
          if(i.eq.3) then
C...tau' -> W + nu'_tau.
            wdtp(i)=fac*sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))*
     &      ((1d0-rm2)**2+(1d0+rm2)*rm1-2d0*rm1**2)
            if(kflr.gt.0) then
              wid2=wids(24,3)
              wid2=wid2*wids(18,2)
            else
              wid2=wids(24,2)
              wid2=wid2*wids(18,3)
            endif
          elseif(i.eq.5) then
C...tau' -> H + nu'_tau.
            wdtp(i)=fac*sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))*
     &      ((1d0+rm2-rm1)*(paru(141)**2+rm2/paru(141)**2)+4d0*rm2)
            if(kflr.gt.0) then
              wid2=wids(37,3)
              wid2=wid2*wids(18,2)
            else
              wid2=wids(37,2)
              wid2=wid2*wids(18,3)
            endif
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  160   continue
 
      elseif(kfla.eq.18) then
C...nu'_tau neutrino.
        fac=(aem/(16d0*xw))*(sh/pmas(24,1)**2)*shr
        do 170 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 170
          rm1=pmas(jamcomp(kfdp(idc,1)),1)**2/sh
          rm2=pmas(jamcomp(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 170
          if(i.eq.2) then
C...nu'_tau -> W + tau'.
            wdtp(i)=fac*sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))*
     &      ((1d0-rm2)**2+(1d0+rm2)*rm1-2d0*rm1**2)
            if(kflr.gt.0) then
              wid2=wids(24,2)
              wid2=wid2*wids(17,2)
            else
              wid2=wids(24,3)
              wid2=wid2*wids(17,3)
            endif
          elseif(i.eq.3) then
C...nu'_tau -> H + tau'.
            wdtp(i)=fac*sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))*
     &      ((1d0+rm2-rm1)*(rm2*paru(141)**2+1d0/paru(141)**2)+4d0*rm2)
            if(kflr.gt.0) then
              wid2=wids(37,2)
              wid2=wid2*wids(17,2)
            else
              wid2=wids(37,3)
              wid2=wid2*wids(17,3)
            endif
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  170   continue
 
      elseif(kfla.eq.21) then
C...QCD:
C***Note that widths are not given in dimensional quantities here.
        do 180 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 180
          rm1=pmas(iabs(kfdp(idc,1)),1)**2/sh
          rm2=pmas(iabs(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 180
          wid2=1d0
          if(i.le.8) then
C...QCD -> q + qbar
            wdtp(i)=(1d0+2d0*rm1)*sqrt(max(0d0,1d0-4d0*rm1))
            if(i.eq.6) wid2=wids(6,1)
            if((i.eq.7.or.i.eq.8)) wid2=wids(i,1)
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  180   continue
 
      elseif(kfla.eq.22) then
C...QED photon.
C***Note that widths are not given in dimensional quantities here.
        do 190 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 190
          rm1=pmas(iabs(kfdp(idc,1)),1)**2/sh
          rm2=pmas(iabs(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 190
          wid2=1d0
          if(i.le.8) then
C...QED -> q + qbar.
            ef=kchg(i,1)/3d0
            fcof=3d0*radc
            if(i.ge.6.and.mstp(35).ge.1) fcof=fcof*pjhfth(sh,sh*rm1,1d0)
            wdtp(i)=fcof*ef**2*(1d0+2d0*rm1)*sqrt(max(0d0,1d0-4d0*rm1))
            if(i.eq.6) wid2=wids(6,1)
            if((i.eq.7.or.i.eq.8)) wid2=wids(i,1)
          elseif(i.le.12) then
C...QED -> l+ + l-.
            ef=kchg(9+2*(i-8),1)/3d0
            wdtp(i)=ef**2*(1d0+2d0*rm1)*sqrt(max(0d0,1d0-4d0*rm1))
            if(i.eq.12) wid2=wids(17,1)
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  190   continue
 
      elseif(kfla.eq.23) then
C...Z0:
        icase=1
        xwc=1d0/(16d0*xw*xw1)
        fac=(aem*xwc/3d0)*shr
  200   continue
        if(mint(61).ge.1.and.icase.eq.2) then
          vint(111)=0d0
          vint(112)=0d0
          vint(114)=0d0
        endif
        if(mint(61).eq.1.and.icase.eq.2) then
          kfi=iabs(mint(15))
          if(kfi.gt.20) kfi=iabs(mint(16))
          ei=kchg(kfi,1)/3d0
          ai=sign(1d0,ei)
          vi=ai-4d0*ei*xwv
          sqmz=pmas(23,1)**2
          hz=shr*wdtp(0)
          if(mstp(43).eq.1.or.mstp(43).eq.3) vint(111)=1d0
          if(mstp(43).eq.3) vint(112)=
     &    2d0*xwc*sh*(sh-sqmz)/((sh-sqmz)**2+hz**2)
          if(mstp(43).eq.2.or.mstp(43).eq.3) vint(114)=
     &    xwc**2*sh**2/((sh-sqmz)**2+hz**2)
        endif
        do 210 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 210
          rm1=pmas(iabs(kfdp(idc,1)),1)**2/sh
          rm2=pmas(iabs(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 210
          wid2=1d0
          if(i.le.8) then
C...Z0 -> q + qbar
            ef=kchg(i,1)/3d0
            af=sign(1d0,ef+0.1d0)
            vf=af-4d0*ef*xwv
            fcof=3d0*radc
            if(i.ge.6.and.mstp(35).ge.1) fcof=fcof*pjhfth(sh,sh*rm1,1d0)
            if(i.eq.6) wid2=wids(6,1)
            if((i.eq.7.or.i.eq.8)) wid2=wids(i,1)
          elseif(i.le.16) then
C...Z0 -> l+ + l-, nu + nubar
            ef=kchg(i+2,1)/3d0
            af=sign(1d0,ef+0.1d0)
            vf=af-4d0*ef*xwv
            fcof=1d0
            if((i.eq.15.or.i.eq.16)) wid2=wids(2+i,1)
          endif
          be34=sqrt(max(0d0,1d0-4d0*rm1))
          if(icase.eq.1) then
            wdtp(i)=fac*fcof*(vf**2*(1d0+2d0*rm1)+af**2*(1d0-4d0*rm1))*
     &      be34
          elseif(mint(61).eq.1.and.icase.eq.2) then
            wdtp(i)=fac*fcof*((ei**2*vint(111)*ef**2+ei*vi*vint(112)*
     &      ef*vf+(vi**2+ai**2)*vint(114)*vf**2)*(1d0+2d0*rm1)+
     &      (vi**2+ai**2)*vint(114)*af**2*(1d0-4d0*rm1))*be34
          elseif(mint(61).eq.2.and.icase.eq.2) then
            fggf=fcof*ef**2*(1d0+2d0*rm1)*be34
            fgzf=fcof*ef*vf*(1d0+2d0*rm1)*be34
            fzzf=fcof*(vf**2*(1d0+2d0*rm1)+af**2*(1d0-4d0*rm1))*be34
          endif
          if(icase.eq.1) wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            if((icase.eq.1.and.mint(61).ne.1).or.
     &      (icase.eq.2.and.mint(61).eq.1)) then
              wdte(i,mdme(idc,1))=wdtp(i)*wid2
              wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+
     &        wdte(i,mdme(idc,1))
              wdte(i,0)=wdte(i,mdme(idc,1))
              wdte(0,0)=wdte(0,0)+wdte(i,0)
            endif
            if(mint(61).eq.2.and.icase.eq.2) then
              if(mstp(43).eq.1.or.mstp(43).eq.3) vint(111)=
     &        vint(111)+fggf*wid2
              if(mstp(43).eq.3) vint(112)=vint(112)+fgzf*wid2
              if(mstp(43).eq.2.or.mstp(43).eq.3) vint(114)=
     &        vint(114)+fzzf*wid2
            endif
          endif
  210   continue
        if(mint(61).ge.1) icase=3-icase
        if(icase.eq.2) goto 200
 
      elseif(kfla.eq.24) then
C...W+/-:
        fac=(aem/(24d0*xw))*shr
        do 220 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 220
          rm1=pmas(iabs(kfdp(idc,1)),1)**2/sh
          rm2=pmas(iabs(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 220
          wid2=1d0
          if(i.le.16) then
C...W+/- -> q + qbar'
            fcof=3d0*radc*vckm((i-1)/4+1,mod(i-1,4)+1)
            if(kflr.gt.0) then
              if(mod(i,4).eq.3) wid2=wids(6,2)
              if(mod(i,4).eq.0) wid2=wids(8,2)
              if(i.ge.13) wid2=wid2*wids(7,3)
            else
              if(mod(i,4).eq.3) wid2=wids(6,3)
              if(mod(i,4).eq.0) wid2=wids(8,3)
              if(i.ge.13) wid2=wid2*wids(7,2)
            endif
          elseif(i.le.20) then
C...W+/- -> l+/- + nu
            fcof=1d0
            if(kflr.gt.0) then
              if(i.eq.20) wid2=wids(17,3)*wids(18,2)
            else
              if(i.eq.20) wid2=wids(17,2)*wids(18,3)
            endif
          endif
          wdtp(i)=fac*fcof*(2d0-rm1-rm2-(rm1-rm2)**2)*
     &    sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  220   continue
 
      elseif(kfla.eq.25.or.kfla.eq.35.or.kfla.eq.36) then
C...h0 (or H0, or A0):
        if(mstp(49).eq.0) then
          fac=(aem/(8d0*xw))*(sh/pmas(24,1)**2)*shr
        else
          fac=(aem/(8d0*xw))*(pmas(kfhigg,1)/pmas(24,1))**2*shr
        endif
        do 260 i=1,mdcy(kfhigg,3)
          idc=i+mdcy(kfhigg,2)-1
          if(mdme(idc,1).lt.0) goto 260
          kfc1=jamcomp(kfdp(idc,1))
          kfc2=jamcomp(kfdp(idc,2))
          rm1=pmas(kfc1,1)**2/sh
          rm2=pmas(kfc2,1)**2/sh
          if(i.ne.16.and.i.ne.17.and.sqrt(rm1)+sqrt(rm2).gt.1d0)
     &    goto 260
          wid2=1d0
 
          if(i.le.8) then
C...h0 -> q + qbar
            wdtp(i)=fac*3d0*rm1*(1d0-4d0*rm1)*sqrt(max(0d0,
     &      1d0-4d0*rm1))*radc
            if(mstp(37).eq.1.and.mstp(2).ge.1) wdtp(i)=wdtp(i)*
     &      (log(max(4d0,parp(37)**2*rm1*sh/paru(117)**2))/
     &      log(max(4d0,sh/paru(117)**2)))**(24d0/(33d0-2d0*mstu(118)))
            if(mstp(4).ge.1.or.ihigg.ge.2) then
              if(mod(i,2).eq.1) wdtp(i)=wdtp(i)*paru(151+10*ihigg)**2
              if(mod(i,2).eq.0) wdtp(i)=wdtp(i)*paru(152+10*ihigg)**2
            endif
            if(i.eq.6) wid2=wids(6,1)
            if((i.eq.7.or.i.eq.8)) wid2=wids(i,1)
 
          elseif(i.le.12) then
C...h0 -> l+ + l-
            wdtp(i)=fac*rm1*(1d0-4d0*rm1)*sqrt(max(0d0,1d0-4d0*rm1))
            if(mstp(4).ge.1.or.ihigg.ge.2) wdtp(i)=wdtp(i)*
     &      paru(153+10*ihigg)**2
            if(i.eq.12) wid2=wids(17,1)
 
          elseif(i.eq.13) then
C...h0 -> g + g; quark loop contribution only
            etare=0d0
            etaim=0d0
            do 230 j=1,2*mstp(1)
              eps=(2d0*pmas(j,1))**2/sh
C...Loop integral; function of eps=4m^2/shat; different for A0.
              if(eps.le.1d0) then
                if(eps.gt.1.d-4) then
                  root=sqrt(1d0-eps)
                  rln=log((1d0+root)/(1d0-root))
                else
                  rln=log(4d0/eps-2d0)
                endif
                phire=-0.25d0*(rln**2-paru(1)**2)
                phiim=0.5d0*paru(1)*rln
              else
                phire=(asin(1d0/sqrt(eps)))**2
                phiim=0d0
              endif
              if(ihigg.le.2) then
                etarej=-0.5d0*eps*(1d0+(1d0-eps)*phire)
                etaimj=-0.5d0*eps*(1d0-eps)*phiim
              else
                etarej=-0.5d0*eps*phire
                etaimj=-0.5d0*eps*phiim
              endif
C...Couplings (=1 for standard model Higgs).
              if(mstp(4).ge.1.or.ihigg.ge.2) then
                if(mod(j,2).eq.1) then
                  etarej=etarej*paru(151+10*ihigg)
                  etaimj=etaimj*paru(151+10*ihigg)
                else
                  etarej=etarej*paru(152+10*ihigg)
                  etaimj=etaimj*paru(152+10*ihigg)
                endif
              endif
              etare=etare+etarej
              etaim=etaim+etaimj
  230       continue
            eta2=etare**2+etaim**2
            wdtp(i)=fac*(as/paru(1))**2*eta2
 
          elseif(i.eq.14) then
C...h0 -> gamma + gamma; quark, lepton, W+- and H+- loop contributions
            etare=0d0
            etaim=0d0
            jmax=3*mstp(1)+1
            if(mstp(4).ge.1.or.ihigg.ge.2) jmax=jmax+1
            do 240 j=1,jmax
              if(j.le.2*mstp(1)) then
                ej=kchg(j,1)/3d0
                eps=(2d0*pmas(j,1))**2/sh
              elseif(j.le.3*mstp(1)) then
                jl=2*(j-2*mstp(1))-1
                ej=kchg(10+jl,1)/3d0
                eps=(2d0*pmas(10+jl,1))**2/sh
              elseif(j.eq.3*mstp(1)+1) then
                eps=(2d0*pmas(24,1))**2/sh
              else
                eps=(2d0*pmas(37,1))**2/sh
              endif
C...Loop integral; function of eps=4m^2/shat.
              if(eps.le.1d0) then
                if(eps.gt.1.d-4) then
                  root=sqrt(1d0-eps)
                  rln=log((1d0+root)/(1d0-root))
                else
                  rln=log(4d0/eps-2d0)
                endif
                phire=-0.25d0*(rln**2-paru(1)**2)
                phiim=0.5d0*paru(1)*rln
              else
                phire=(asin(1d0/sqrt(eps)))**2
                phiim=0d0
              endif
              if(j.le.3*mstp(1)) then
C...Fermion loops: loop integral different for A0; charges.
                if(ihigg.le.2) then
                  phipre=-0.5d0*eps*(1d0+(1d0-eps)*phire)
                  phipim=-0.5d0*eps*(1d0-eps)*phiim
                else
                  phipre=-0.5d0*eps*phire
                  phipim=-0.5d0*eps*phiim
                endif
                if(j.le.2*mstp(1).and.mod(j,2).eq.1) then
                  ejc=3d0*ej**2
                  ejh=paru(151+10*ihigg)
                elseif(j.le.2*mstp(1)) then
                  ejc=3d0*ej**2
                  ejh=paru(152+10*ihigg)
                else
                  ejc=ej**2
                  ejh=paru(153+10*ihigg)
                endif
                if(mstp(4).eq.0.and.ihigg.eq.1) ejh=1d0
                etarej=ejc*ejh*phipre
                etaimj=ejc*ejh*phipim
              elseif(j.eq.3*mstp(1)+1) then
C...W loops: loop integral and charges.
                etarej=0.5d0+0.75d0*eps*(1d0+(2d0-eps)*phire)
                etaimj=0.75d0*eps*(2d0-eps)*phiim
                if(mstp(4).ge.1.or.ihigg.ge.2) then
                  etarej=etarej*paru(155+10*ihigg)
                  etaimj=etaimj*paru(155+10*ihigg)
                endif
              else
C...Charged H loops: loop integral and charges.
                fachhh=(pmas(24,1)/pmas(37,1))**2*
     &          paru(158+10*ihigg+2*(ihigg/3))
                etarej=eps*(1d0-eps*phire)*fachhh
                etaimj=-eps**2*phiim*fachhh
              endif
              etare=etare+etarej
              etaim=etaim+etaimj
  240       continue
            eta2=etare**2+etaim**2
            wdtp(i)=fac*(aem/paru(1))**2*0.5d0*eta2
 
          elseif(i.eq.15) then
C...h0 -> gamma + Z0; quark, lepton, W and H+- loop contributions
            etare=0d0
            etaim=0d0
            jmax=3*mstp(1)+1
            if(mstp(4).ge.1.or.ihigg.ge.2) jmax=jmax+1
            do 250 j=1,jmax
              if(j.le.2*mstp(1)) then
                ej=kchg(j,1)/3d0
                aj=sign(1d0,ej+0.1d0)
                vj=aj-4d0*ej*xwv
                eps=(2d0*pmas(j,1))**2/sh
                epsp=(2d0*pmas(j,1)/pmas(23,1))**2
              elseif(j.le.3*mstp(1)) then
                jl=2*(j-2*mstp(1))-1
                ej=kchg(10+jl,1)/3d0
                aj=sign(1d0,ej+0.1d0)
                vj=aj-4d0*ej*xwv
                eps=(2d0*pmas(10+jl,1))**2/sh
                epsp=(2d0*pmas(10+jl,1)/pmas(23,1))**2
              else
                eps=(2d0*pmas(24,1))**2/sh
                epsp=(2d0*pmas(24,1)/pmas(23,1))**2
              endif
C...Loop integrals; functions of eps=4m^2/shat and eps'=4m^2/m_Z^2.
              if(eps.le.1d0) then
                root=sqrt(1d0-eps)
                if(eps.gt.1.d-4) then
                  rln=log((1d0+root)/(1d0-root))
                else
                  rln=log(4d0/eps-2d0)
                endif
                phire=-0.25d0*(rln**2-paru(1)**2)
                phiim=0.5d0*paru(1)*rln
                psire=0.5d0*root*rln
                psiim=-0.5d0*root*paru(1)
              else
                phire=(asin(1d0/sqrt(eps)))**2
                phiim=0d0
                psire=sqrt(eps-1d0)*asin(1d0/sqrt(eps))
                psiim=0d0
              endif
              if(epsp.le.1d0) then
                root=sqrt(1d0-epsp)
                if(epsp.gt.1.d-4) then
                  rln=log((1d0+root)/(1d0-root))
                else
                  rln=log(4d0/epsp-2d0)
                endif
                phirep=-0.25d0*(rln**2-paru(1)**2)
                phiimp=0.5d0*paru(1)*rln
                psirep=0.5d0*root*rln
                psiimp=-0.5d0*root*paru(1)
              else
                phirep=(asin(1d0/sqrt(epsp)))**2
                phiimp=0d0
                psirep=sqrt(epsp-1d0)*asin(1d0/sqrt(epsp))
                psiimp=0d0
              endif
              fxyre=eps*epsp/(8d0*(eps-epsp))*(1d0+eps*epsp/(eps-epsp)*
     &        (phire-phirep)+2d0*eps/(eps-epsp)*(psire-psirep))
              fxyim=eps**2*epsp/(8d0*(eps-epsp)**2)*
     &        (epsp*(phiim-phiimp)+2d0*(psiim-psiimp))
              f1re=-eps*epsp/(2d0*(eps-epsp))*(phire-phirep)
              f1im=-eps*epsp/(2d0*(eps-epsp))*(phiim-phiimp)
              if(j.le.3*mstp(1)) then
C...Fermion loops: loop integral different for A0; charges.
                if(ihigg.eq.3) fxyre=0d0
                if(ihigg.eq.3) fxyim=0d0
                if(j.le.2*mstp(1).and.mod(j,2).eq.1) then
                  ejc=-3d0*ej*vj
                  ejh=paru(151+10*ihigg)
                elseif(j.le.2*mstp(1)) then
                  ejc=-3d0*ej*vj
                  ejh=paru(152+10*ihigg)
                else
                  ejc=-ej*vj
                  ejh=paru(153+10*ihigg)
                endif
                if(mstp(4).eq.0.and.ihigg.eq.1) ejh=1d0
                etarej=ejc*ejh*(fxyre-0.25d0*f1re)
                etaimj=ejc*ejh*(fxyim-0.25d0*f1im)
              elseif(j.eq.3*mstp(1)+1) then
C...W loops: loop integral and charges.
                heps=(1d0+2d0/eps)*xw/xw1-(5d0+2d0/eps)
                etarej=-xw1*((3d0-xw/xw1)*f1re+heps*fxyre)
                etaimj=-xw1*((3d0-xw/xw1)*f1im+heps*fxyim)
                if(mstp(4).ge.1.or.ihigg.ge.2) then
                  etarej=etarej*paru(155+10*ihigg)
                  etaimj=etaimj*paru(155+10*ihigg)
                endif
              else
C...Charged H loops: loop integral and charges.
                fachhh=(pmas(24,1)/pmas(37,1))**2*(1d0-2d0*xw)*
     &          paru(158+10*ihigg+2*(ihigg/3))
                etarej=fachhh*fxyre
                etaimj=fachhh*fxyim
              endif
              etare=etare+etarej
              etaim=etaim+etaimj
  250       continue
            eta2=(etare**2+etaim**2)/(xw*xw1)
            wdtp(i)=fac*(aem/paru(1))**2*(1d0-pmas(23,1)**2/sh)**3*eta2
            wid2=wids(23,2)
 
          elseif(i.le.17) then
C...h0 -> Z0 + Z0, W+ + W-
            pm1=pmas(iabs(kfdp(idc,1)),1)
            pg1=pmas(iabs(kfdp(idc,1)),2)
            if(mint(62).ge.1) then
              if(mstp(42).eq.0.or.(4d0*(pm1+10d0*pg1)**2.lt.sh.and.
     &        ckin(46).lt.ckin(45).and.ckin(48).lt.ckin(47).and.
     &        max(ckin(45),ckin(47)).lt.pm1-10d0*pg1)) then
                mofsv(ihigg,i-15)=0
                widw=(1d0-4d0*rm1+12d0*rm1**2)*sqrt(max(0d0,
     &          1d0-4d0*rm1))
                wid2=1d0
              else
                mofsv(ihigg,i-15)=1
                rmas=sqrt(max(0d0,sh))
                call pjofsh(1,kfla,kfdp(idc,1),kfdp(idc,2),rmas,widw,
     &          wid2)
                widwsv(ihigg,i-15)=widw
                wid2sv(ihigg,i-15)=wid2
              endif
            else
              if(mofsv(ihigg,i-15).eq.0) then
                widw=(1d0-4d0*rm1+12d0*rm1**2)*sqrt(max(0d0,
     &          1d0-4d0*rm1))
                wid2=1d0
              else
                widw=widwsv(ihigg,i-15)
                wid2=wid2sv(ihigg,i-15)
              endif
            endif
            wdtp(i)=fac*widw/(2d0*(18-i))
            if(mstp(4).ge.1.or.ihigg.ge.2) wdtp(i)=wdtp(i)*
     &      paru(138+i+10*ihigg)**2
            wid2=wid2*wids(7+i,1)
 
          elseif(i.eq.18.and.kfla.eq.35) then
C***H0 -> Z0 + h0 (not yet implemented).
 
          elseif(i.eq.19.and.kfla.eq.35) then
C...H0 -> h0 + h0.
            wdtp(i)=fac*paru(176)**2*0.25d0*pmas(23,1)**4/sh**2*
     &      sqrt(max(0d0,1d0-4d0*rm1))
            wid2=wids(25,2)**2
 
          elseif(i.eq.20.and.kfla.eq.35) then
C...H0 -> A0 + A0.
            wdtp(i)=fac*paru(177)**2*0.25d0*pmas(23,1)**4/sh**2*
     &      sqrt(max(0d0,1d0-4d0*rm1))
            wid2=wids(36,2)**2
 
          elseif(i.eq.18.and.kfla.eq.36) then
C...A0 -> Z0 + h0.
            wdtp(i)=fac*paru(186)**2*0.5d0*sqrt(max(0d0,
     &      (1d0-rm1-rm2)**2-4d0*rm1*rm2))**3
            wid2=wids(23,2)*wids(25,2)
 
CMRENNA++
          else
C...Add in SUSY decays (two-body) by rescaling by phase space factor.
            rm10=rm1*sh/pmr**2
            rm20=rm2*sh/pmr**2
            wfac0=1d0+rm10**2+rm20**2-2d0*(rm10+rm20+rm10*rm20)
            wfac=1d0+rm1**2+rm2**2-2d0*(rm1+rm2+rm1*rm2)
            if(wfac.le.0d0 .or. wfac0.le.0d0) then
              wfac=0d0
            else
              wfac=wfac/wfac0
            endif
            wdtp(i)=pmas(kfla,2)*brat(idc)*(shr/pmr)*sqrt(wfac)
CMRENNA--
            if(kfc2.eq.kfc1) then
              wid2=wids(kfc1,1)
            else
              ksgn1=2
              if(kfdp(idc,1).lt.0) ksgn1=3
              ksgn2=2
              if(kfdp(idc,2).lt.0) ksgn2=3
              wid2=wids(kfc1,ksgn1)*wids(kfc2,ksgn2)
            endif
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  260   continue
 
      elseif(kfla.eq.32) then
C...Z'0:
        icase=1
        xwc=1d0/(16d0*xw*xw1)
        fac=(aem*xwc/3d0)*shr
        vint(117)=0d0
  270   continue
        if(mint(61).ge.1.and.icase.eq.2) then
          vint(111)=0d0
          vint(112)=0d0
          vint(113)=0d0
          vint(114)=0d0
          vint(115)=0d0
          vint(116)=0d0
        endif
        if(mint(61).eq.1.and.icase.eq.2) then
          kfai=iabs(mint(15))
          ei=kchg(kfai,1)/3d0
          ai=sign(1d0,ei+0.1d0)
          vi=ai-4d0*ei*xwv
          kfaic=1
          if(kfai.le.10.and.mod(kfai,2).eq.0) kfaic=2
          if(kfai.gt.10.and.mod(kfai,2).ne.0) kfaic=3
          if(kfai.gt.10.and.mod(kfai,2).eq.0) kfaic=4
          vpi=paru(119+2*kfaic)
          api=paru(120+2*kfaic)
          sqmz=pmas(23,1)**2
          hz=shr*fac*vint(117)
          sqmzp=pmas(32,1)**2
          hzp=shr*fac*wdtp(0)
          if(mstp(44).eq.1.or.mstp(44).eq.4.or.mstp(44).eq.5.or.
     &    mstp(44).eq.7) vint(111)=1d0
          if(mstp(44).eq.4.or.mstp(44).eq.7) vint(112)=
     &    2d0*xwc*sh*(sh-sqmz)/((sh-sqmz)**2+hz**2)
          if(mstp(44).eq.5.or.mstp(44).eq.7) vint(113)=
     &    2d0*xwc*sh*(sh-sqmzp)/((sh-sqmzp)**2+hzp**2)
          if(mstp(44).eq.2.or.mstp(44).eq.4.or.mstp(44).eq.6.or.
     &    mstp(44).eq.7) vint(114)=xwc**2*sh**2/((sh-sqmz)**2+hz**2)
          if(mstp(44).eq.6.or.mstp(44).eq.7) vint(115)=
     &    2d0*xwc**2*sh**2*((sh-sqmz)*(sh-sqmzp)+hz*hzp)/
     &    (((sh-sqmz)**2+hz**2)*((sh-sqmzp)**2+hzp**2))
          if(mstp(44).eq.3.or.mstp(44).eq.5.or.mstp(44).eq.6.or.
     &    mstp(44).eq.7) vint(116)=xwc**2*sh**2/((sh-sqmzp)**2+hzp**2)
        endif
        do 280 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 280
          rm1=pmas(jamcomp(kfdp(idc,1)),1)**2/sh
          rm2=pmas(jamcomp(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0.or.mdme(idc,1).lt.0) goto 280
          wid2=1d0
          if(i.le.16) then
            if(i.le.8) then
C...Z'0 -> q + qbar
              ef=kchg(i,1)/3d0
              af=sign(1d0,ef+0.1d0)
              vf=af-4d0*ef*xwv
              vpf=paru(123-2*mod(i,2))
              apf=paru(124-2*mod(i,2))
              fcof=3d0*radc
              if(i.ge.6.and.mstp(35).ge.1) fcof=fcof*
     &        pjhfth(sh,sh*rm1,1d0)
              if(i.eq.6) wid2=wids(6,1)
              if((i.eq.7.or.i.eq.8)) wid2=wids(i,1)
            elseif(i.le.16) then
C...Z'0 -> l+ + l-, nu + nubar
              ef=kchg(i+2,1)/3d0
              af=sign(1d0,ef+0.1d0)
              vf=af-4d0*ef*xwv
              vpf=paru(127-2*mod(i,2))
              apf=paru(128-2*mod(i,2))
              fcof=1d0
              if((i.eq.15.or.i.eq.16)) wid2=wids(2+i,1)
            endif
            be34=sqrt(max(0d0,1d0-4d0*rm1))
            if(icase.eq.1) then
              wdtpz=fcof*(vf**2*(1d0+2d0*rm1)+af**2*(1d0-4d0*rm1))*be34
              wdtp(i)=fac*fcof*(vpf**2*(1d0+2d0*rm1)+
     &        apf**2*(1d0-4d0*rm1))*be34
            elseif(mint(61).eq.1.and.icase.eq.2) then
              wdtp(i)=fac*fcof*((ei**2*vint(111)*ef**2+ei*vi*vint(112)*
     &        ef*vf+ei*vpi*vint(113)*ef*vpf+(vi**2+ai**2)*vint(114)*
     &        vf**2+(vi*vpi+ai*api)*vint(115)*vf*vpf+(vpi**2+api**2)*
     &        vint(116)*vpf**2)*(1d0+2d0*rm1)+((vi**2+ai**2)*vint(114)*
     &        af**2+(vi*vpi+ai*api)*vint(115)*af*apf+(vpi**2+api**2)*
     &        vint(116)*apf**2)*(1d0-4d0*rm1))*be34
            elseif(mint(61).eq.2) then
              fggf=fcof*ef**2*(1d0+2d0*rm1)*be34
              fgzf=fcof*ef*vf*(1d0+2d0*rm1)*be34
              fgzpf=fcof*ef*vpf*(1d0+2d0*rm1)*be34
              fzzf=fcof*(vf**2*(1d0+2d0*rm1)+af**2*(1d0-4d0*rm1))*be34
              fzzpf=fcof*(vf*vpf*(1d0+2d0*rm1)+af*apf*(1d0-4d0*rm1))*
     &        be34
              fzpzpf=fcof*(vpf**2*(1d0+2d0*rm1)+apf**2*(1d0-4d0*rm1))*
     &        be34
            endif
          elseif(i.eq.17) then
C...Z'0 -> W+ + W-
            wdtpzp=paru(129)**2*xw1**2*
     &      sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))**3*
     &      (1d0+10d0*rm1+10d0*rm2+rm1**2+rm2**2+10d0*rm1*rm2)
            if(icase.eq.1) then
              wdtpz=0d0
              wdtp(i)=fac*wdtpzp
            elseif(mint(61).eq.1.and.icase.eq.2) then
              wdtp(i)=fac*(vpi**2+api**2)*vint(116)*wdtpzp
            elseif(mint(61).eq.2) then
              fggf=0d0
              fgzf=0d0
              fgzpf=0d0
              fzzf=0d0
              fzzpf=0d0
              fzpzpf=wdtpzp
            endif
            wid2=wids(24,1)
          elseif(i.eq.18) then
C...Z'0 -> H+ + H-
            czc=2d0*(1d0-2d0*xw)
            be34c=(1d0-4d0*rm1)*sqrt(max(0d0,1d0-4d0*rm1))
            if(icase.eq.1) then
              wdtpz=0.25d0*paru(142)**2*czc**2*be34c
              wdtp(i)=fac*0.25d0*paru(143)**2*czc**2*be34c
            elseif(mint(61).eq.1.and.icase.eq.2) then
              wdtp(i)=fac*0.25d0*(ei**2*vint(111)+paru(142)*ei*vi*
     &        vint(112)*czc+paru(143)*ei*vpi*vint(113)*czc+paru(142)**2*
     &        (vi**2+ai**2)*vint(114)*czc**2+paru(142)*paru(143)*
     &        (vi*vpi+ai*api)*vint(115)*czc**2+paru(143)**2*
     &        (vpi**2+api**2)*vint(116)*czc**2)*be34c
            elseif(mint(61).eq.2) then
              fggf=0.25d0*be34c
              fgzf=0.25d0*paru(142)*czc*be34c
              fgzpf=0.25d0*paru(143)*czc*be34c
              fzzf=0.25d0*paru(142)**2*czc**2*be34c
              fzzpf=0.25d0*paru(142)*paru(143)*czc**2*be34c
              fzpzpf=0.25d0*paru(143)**2*czc**2*be34c
            endif
            wid2=wids(37,1)
          elseif(i.eq.19) then
C...Z'0 -> Z0 + gamma.
          elseif(i.eq.20) then
C...Z'0 -> Z0 + h0
            flam=sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))
            wdtpzp=paru(145)**2*4d0*abs(1d0-2d0*xw)*
     &      (3d0*rm1+0.25d0*flam**2)*flam
            if(icase.eq.1) then
              wdtpz=0d0
              wdtp(i)=fac*wdtpzp
            elseif(mint(61).eq.1.and.icase.eq.2) then
              wdtp(i)=fac*(vpi**2+api**2)*vint(116)*wdtpzp
            elseif(mint(61).eq.2) then
              fggf=0d0
              fgzf=0d0
              fgzpf=0d0
              fzzf=0d0
              fzzpf=0d0
              fzpzpf=wdtpzp
            endif
            wid2=wids(23,2)*wids(25,2)
          elseif(i.eq.21.or.i.eq.22) then
C...Z' -> h0 + A0 or H0 + A0.
            be34c=sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))**3
            if(i.eq.21) then
              czah=paru(186)
              czpah=paru(188)
            else
              czah=paru(187)
              czpah=paru(189)
            endif
            if(icase.eq.1) then
              wdtpz=czah**2*be34c
              wdtp(i)=fac*czpah**2*be34c
            elseif(mint(61).eq.1.and.icase.eq.2) then
              wdtp(i)=fac*(czah**2*(vi**2+ai**2)*vint(114)+czah*czpah*
     &        (vi*vpi+ai*api)*vint(115)+czpah**2*(vpi**2+api**2)*
     &        vint(116))*be34c
            elseif(mint(61).eq.2) then
              fggf=0d0
              fgzf=0d0
              fgzpf=0d0
              fzzf=czah**2*be34c
              fzzpf=czah*czpah*be34c
              fzpzpf=czpah**2*be34c
            endif
            if(i.eq.21) wid2=wids(25,2)*wids(36,2)
            if(i.eq.22) wid2=wids(35,2)*wids(36,2)
          endif
          if(icase.eq.1) then
            vint(117)=vint(117)+wdtpz
            wdtp(0)=wdtp(0)+wdtp(i)
          endif
          if(mdme(idc,1).gt.0) then
            if((icase.eq.1.and.mint(61).ne.1).or.
     &      (icase.eq.2.and.mint(61).eq.1)) then
              wdte(i,mdme(idc,1))=wdtp(i)*wid2
              wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+
     &        wdte(i,mdme(idc,1))
              wdte(i,0)=wdte(i,mdme(idc,1))
              wdte(0,0)=wdte(0,0)+wdte(i,0)
            endif
            if(mint(61).eq.2.and.icase.eq.2) then
              if(mstp(44).eq.1.or.mstp(44).eq.4.or.mstp(44).eq.5.or.
     &        mstp(44).eq.7) vint(111)=vint(111)+fggf*wid2
              if(mstp(44).eq.4.or.mstp(44).eq.7) vint(112)=vint(112)+
     &        fgzf*wid2
              if(mstp(44).eq.5.or.mstp(44).eq.7) vint(113)=vint(113)+
     &        fgzpf*wid2
              if(mstp(44).eq.2.or.mstp(44).eq.4.or.mstp(44).eq.6.or.
     &        mstp(44).eq.7) vint(114)=vint(114)+fzzf*wid2
              if(mstp(44).eq.6.or.mstp(44).eq.7) vint(115)=vint(115)+
     &        fzzpf*wid2
              if(mstp(44).eq.3.or.mstp(44).eq.5.or.mstp(44).eq.6.or.
     &        mstp(44).eq.7) vint(116)=vint(116)+fzpzpf*wid2
            endif
          endif
  280   continue
        if(mint(61).ge.1) icase=3-icase
        if(icase.eq.2) goto 270
 
      elseif(kfla.eq.34) then
C...W'+/-:
        fac=(aem/(24d0*xw))*shr
        do 290 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 290
          rm1=pmas(jamcomp(kfdp(idc,1)),1)**2/sh
          rm2=pmas(jamcomp(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 290
          wid2=1d0
          if(i.le.20) then
            if(i.le.16) then
C...W'+/- -> q + qbar'
              fcof=3d0*radc*(paru(131)**2+paru(132)**2)*
     &        vckm((i-1)/4+1,mod(i-1,4)+1)
              if(kflr.gt.0) then
                if(mod(i,4).eq.3) wid2=wids(6,2)
                if(mod(i,4).eq.0) wid2=wids(8,2)
                if(i.ge.13) wid2=wid2*wids(7,3)
              else
                if(mod(i,4).eq.3) wid2=wids(6,3)
                if(mod(i,4).eq.0) wid2=wids(8,3)
                if(i.ge.13) wid2=wid2*wids(7,2)
              endif
            elseif(i.le.20) then
C...W'+/- -> l+/- + nu
              fcof=paru(133)**2+paru(134)**2
              if(kflr.gt.0) then
                if(i.eq.20) wid2=wids(17,3)*wids(18,2)
              else
                if(i.eq.20) wid2=wids(17,2)*wids(18,3)
              endif
            endif
            wdtp(i)=fac*fcof*0.5d0*(2d0-rm1-rm2-(rm1-rm2)**2)*
     &      sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))
          elseif(i.eq.21) then
C...W'+/- -> W+/- + Z0
            wdtp(i)=fac*paru(135)**2*0.5d0*xw1*(rm1/rm2)*
     &      sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))**3*
     &      (1d0+10d0*rm1+10d0*rm2+rm1**2+rm2**2+10d0*rm1*rm2)
            if(kflr.gt.0) wid2=wids(24,2)*wids(23,2)
            if(kflr.lt.0) wid2=wids(24,3)*wids(23,2)
          elseif(i.eq.23) then
C...W'+/- -> W+/- + h0
            flam=sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))
            wdtp(i)=fac*paru(146)**2*2d0*(3d0*rm1+0.25d0*flam**2)*flam
            if(kflr.gt.0) wid2=wids(24,2)*wids(25,2)
            if(kflr.lt.0) wid2=wids(24,3)*wids(25,2)
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  290   continue
 
      elseif(kfla.eq.37) then
C...H+/-:
        fac=(aem/(8d0*xw))*(sh/pmas(24,1)**2)*shr
        do 300 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 300
          kfc1=jamcomp(kfdp(idc,1))
          kfc2=jamcomp(kfdp(idc,2))
          rm1=pmas(kfc1,1)**2/sh
          rm2=pmas(kfc2,1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 300
          wid2=1d0
          if(i.le.4) then
C...H+/- -> q + qbar'
            rm1r=rm1
            if(mstp(37).eq.1.and.mstp(2).ge.1) rm1r=rm1*
     &      (log(max(4d0,parp(37)**2*rm1*sh/paru(117)**2))/
     &      log(max(4d0,sh/paru(117)**2)))**(24d0/(33d0-2d0*mstu(118)))
            wdtp(i)=fac*3d0*radc*((rm1r*paru(141)**2+rm2/paru(141)**2)*
     &      (1d0-rm1r-rm2)-4d0*rm1r*rm2)*
     &      sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))
            if(kflr.gt.0) then
              if(i.eq.3) wid2=wids(6,2)
              if(i.eq.4) wid2=wids(7,3)*wids(8,2)
            else
              if(i.eq.3) wid2=wids(6,3)
              if(i.eq.4) wid2=wids(7,2)*wids(8,3)
            endif
          elseif(i.le.8) then
C...H+/- -> l+/- + nu
            wdtp(i)=fac*((rm1*paru(141)**2+rm2/paru(141)**2)*
     &      (1d0-rm1-rm2)-4d0*rm1*rm2)*sqrt(max(0d0,(1d0-rm1-rm2)**2-
     &      4d0*rm1*rm2))
            if(kflr.gt.0) then
              if(i.eq.8) wid2=wids(17,3)*wids(18,2)
            else
              if(i.eq.8) wid2=wids(17,2)*wids(18,3)
            endif
          elseif(i.eq.9) then
C...H+/- -> W+/- + h0.
            wdtp(i)=fac*paru(195)**2*0.5d0*sqrt(max(0d0,
     &      (1d0-rm1-rm2)**2-4d0*rm1*rm2))**3
            if(kflr.gt.0) wid2=wids(24,2)*wids(25,2)
            if(kflr.lt.0) wid2=wids(24,3)*wids(25,2)
 
CMRENNA++
          else
C...Add in SUSY decays (two-body) by rescaling by phase space factor.
            rm10=rm1*sh/pmr**2
            rm20=rm2*sh/pmr**2
            wfac0=1d0+rm10**2+rm20**2-2d0*(rm10+rm20+rm10*rm20)
            wfac=1d0+rm1**2+rm2**2-2d0*(rm1+rm2+rm1*rm2)
            if(wfac.le.0d0 .or. wfac0.le.0d0) then
              wfac=0d0
            else
              wfac=wfac/wfac0
            endif
            wdtp(i)=pmas(kc,2)*brat(idc)*(shr/pmr)*sqrt(wfac)
CMRENNA--
            ksgn1=2
            if(kfls*kfdp(idc,1).lt.0.and.kchg(kfc1,3).eq.1) ksgn1=3
            ksgn2=2
            if(kfls*kfdp(idc,2).lt.0.and.kchg(kfc2,3).eq.1) ksgn2=3
            wid2=wids(kfc1,ksgn1)*wids(kfc2,ksgn2)
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  300   continue
 
      elseif(kfla.eq.38) then
C...Techni-eta.
        fac=(sh/parp(46)**2)*shr
        do 310 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 310
          rm1=pmas(jamcomp(kfdp(idc,1)),1)**2/sh
          rm2=pmas(jamcomp(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 310
          wid2=1d0
          if(i.le.2) then
            wdtp(i)=fac*rm1*sqrt(max(0d0,1d0-4d0*rm1))/(4d0*paru(1))
            if(i.eq.2) wid2=wids(6,1)
          else
            wdtp(i)=fac*5d0*as**2/(96d0*paru(1)**3)
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  310   continue
 
      elseif(kfla.eq.39) then
C...LQ (leptoquark).
        fac=(aem/4d0)*paru(151)*shr
        do 320 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 320
          rm1=pmas(jamcomp(kfdp(idc,1)),1)**2/sh
          rm2=pmas(jamcomp(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 320
          wdtp(i)=fac*sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))**3
          wid2=1d0
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  320   continue
 
      elseif(kfla.eq.40) then
C...R:
        fac=(aem/(12d0*xw))*shr
        do 330 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 330
          rm1=pmas(jamcomp(kfdp(idc,1)),1)**2/sh
          rm2=pmas(jamcomp(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 330
          wid2=1d0
          if(i.le.6) then
C...R -> q + qbar'
            fcof=3d0*radc
          elseif(i.le.9) then
C...R -> l+ + l'-
            fcof=1d0
          endif
          wdtp(i)=fac*fcof*(2d0-rm1-rm2-(rm1-rm2)**2)*
     &    sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))
          if(kflr.gt.0) then
            if(i.eq.4) wid2=wids(6,3)
            if(i.eq.5) wid2=wids(7,3)
            if(i.eq.6) wid2=wids(6,2)*wids(8,3)
            if(i.eq.9) wid2=wids(17,3)
          else
            if(i.eq.4) wid2=wids(6,2)
            if(i.eq.5) wid2=wids(7,2)
            if(i.eq.6) wid2=wids(6,3)*wids(8,2)
            if(i.eq.9) wid2=wids(17,2)
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  330   continue
 
      elseif(kfla.eq.51.or.kfla.eq.52) then
C...Techni-pi0 and techni-pi+-:
        fac=(3d0/(32d0*paru(1)*parp(142)**2))*shr
        do 340 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 340
          pm1=pmas(jamcomp(kfdp(idc,1)),1)
          pm2=pmas(jamcomp(kfdp(idc,2)),1)
          rm1=pm1**2/sh
          rm2=pm2**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 340
          wid2=1d0
C...pi_tech -> f + f'.
          fcof=1d0
          if(iabs(kfdp(idc,1)).lt.10) fcof=3d0*radc
          wdtp(i)=fac*fcof*(pm1+pm2)**2*
     &    sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  340   continue
 
      elseif(kfla.eq.53) then
C...Techni-pi'0 not yet implemented.
 
      elseif(kfla.eq.54) then
C...Techni-rho0:
        alprht=2.91d0*(3d0/parp(144))
        fac=(alprht/12d0)*shr
        facf=(1d0/6d0)*(aem**2/alprht)*(pmas(kfla,1)**4/shr**3)
        sqmz=pmas(23,1)**2
        gmmz=pmas(23,1)*pmas(23,2)
        xwrht=(1d0-2d0*xw)/(4d0*xw*(1d0-xw))
        bwzr=xwrht*sh*(sh-sqmz)/((sh-sqmz)**2+gmmz**2)
        bwzi=xwrht*sh*gmmz/((sh-sqmz)**2+gmmz**2)
        do 350 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 350
          rm1=pmas(jamcomp(kfdp(idc,1)),1)**2/sh
          rm2=pmas(jamcomp(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 350
          if(i.eq.1) then
C...rho_tech0 -> W+ + W-.
            wdtp(i)=fac*parp(141)**4*
     &      sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))**3
            wid2=wids(24,1)
          elseif(i.eq.2) then
C...rho_tech0 -> W+ + pi_tech-.
            wdtp(i)=fac*parp(141)**2*(1d0-parp(141)**2)*
     &      sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))**3
            wid2=wids(24,2)*wids(52,3)
          elseif(i.eq.3) then
C...rho_tech0 -> pi_tech+ + W-.
            wdtp(i)=fac*parp(141)**2*(1d0-parp(141)**2)*
     &      sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))**3
            wid2=wids(52,2)*wids(24,3)
          elseif(i.eq.4) then
C...rho_tech0 -> pi_tech+ + pi_tech-.
            wdtp(i)=fac*(1d0-parp(141)**2)**2*
     &      sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))**3
            wid2=wids(52,1)
          else
C...rho_tech0 -> f + fbar.
            wid2=1d0
            if(i.le.12) then
              ia=i-4
              fcof=3d0*radc
              if(ia.ge.6.and.ia.le.8) wid2=wids(ia,1)
            else
              ia=i-2
              fcof=1d0
              if(ia.ge.17) wid2=wids(ia,1)
            endif
            ei=kchg(ia,1)/3d0
            ai=sign(1d0,ei+0.1d0)
            vi=ai-4d0*ei*xwv
            vali=0.5d0*(vi+ai)
            vari=0.5d0*(vi-ai)
            wdtp(i)=facf*fcof*(1d0-rm1)*sqrt(max(0d0,1d0-4d0*rm1))*
     &      ((ei+vali*bwzr)**2+(vali*bwzi)**2+
     &      (ei+vari*bwzr)**2+(vari*bwzi)**2)
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  350   continue
 
      elseif(kfla.eq.55) then
C...Techni-rho+/-:
        alprht=2.91d0*(3d0/parp(144))
        fac=(alprht/12d0)*shr
        sqmw=pmas(24,1)**2
        gmmw=pmas(24,1)*pmas(24,2)
        facf=(1d0/6d0)*(aem**2/alprht)*(pmas(kfla,1)**4/shr**3)*
     &  (0.25d0/xw**2)*sh**2/((sh-sqmw)**2+gmmw**2)
        do 360 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 360
          rm1=pmas(jamcomp(kfdp(idc,1)),1)**2/sh
          rm2=pmas(jamcomp(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 360
          if(i.eq.1) then
C...rho_tech+ -> W+ + Z0.
            wdtp(i)=fac*parp(141)**4*
     &      sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))**3
            if(kflr.gt.0) then
              wid2=wids(24,2)*wids(23,2)
            else
              wid2=wids(24,3)*wids(23,2)
            endif
          elseif(i.eq.2) then
C...rho_tech+ -> W+ + pi_tech0.
            wdtp(i)=fac*parp(141)**2*(1d0-parp(141)**2)*
     &      sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))**3
            if(kflr.gt.0) then
              wid2=wids(24,2)*wids(51,2)
            else
              wid2=wids(24,3)*wids(51,2)
            endif
          elseif(i.eq.3) then
C...rho_tech+ -> pi_tech+ + Z0.
            wdtp(i)=fac*parp(141)**2*(1d0-parp(141)**2)*
     &      sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))**3
            if(kflr.gt.0) then
              wid2=wids(52,2)*wids(23,2)
            else
              wid2=wids(52,3)*wids(23,2)
            endif
          elseif(i.eq.4) then
C...rho_tech+ -> pi_tech+ + pi_tech0.
            wdtp(i)=fac*(1d0-parp(141)**2)**2*
     &      sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))**3
            if(kflr.gt.0) then
              wid2=wids(52,2)*wids(51,2)
            else
              wid2=wids(52,3)*wids(51,2)
            endif
          else
C...rho_tech+ -> f + fbar'.
            ia=i-4
            wid2=1d0
            if(ia.le.16) then
              fcof=3d0*radc*vckm((ia-1)/4+1,mod(ia-1,4)+1)
              if(kflr.gt.0) then
                if(mod(ia,4).eq.3) wid2=wids(6,2)
                if(mod(ia,4).eq.0) wid2=wids(8,2)
                if(ia.ge.13) wid2=wid2*wids(7,3)
              else
                if(mod(ia,4).eq.3) wid2=wids(6,3)
                if(mod(ia,4).eq.0) wid2=wids(8,3)
                if(ia.ge.13) wid2=wid2*wids(7,2)
              endif
            else
              fcof=1d0
              if(kflr.gt.0) then
                if(ia.eq.20) wid2=wids(17,3)*wids(18,2)
              else
                if(ia.eq.20) wid2=wids(17,2)*wids(18,3)
              endif
            endif
            wdtp(i)=facf*fcof*0.5d0*(2d0-rm1-rm2-(rm1-rm2)**2)*
     &      sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  360   continue
 
      elseif(kfla.eq.56) then
C...Techni-omega:
        alprht=2.91d0*(3d0/parp(144))
        fac=(aem/24d0)*(shr**3/parp(145)**2)
        facf=(1d0/6d0)*(aem**2/alprht)*(pmas(kfla,1)**4/shr**3)*
     &  (2d0*parp(143)-1d0)**2
        sqmz=pmas(23,1)**2
        gmmz=pmas(23,1)*pmas(23,2)
        bwzr=(0.5d0/(1d0-xw))*sh*(sh-sqmz)/((sh-sqmz)**2+gmmz**2)
        bwzi=(0.5d0/(1d0-xw))*sh*gmmz/((sh-sqmz)**2+gmmz**2)
        do 370 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 370
          rm1=pmas(jamcomp(kfdp(idc,1)),1)**2/sh
          rm2=pmas(jamcomp(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 370
          if(i.eq.1) then
C...omega_tech0 -> gamma + pi_tech0.
            wdtp(i)=fac*
     &      sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))**3
            wid2=wids(51,2)
          elseif(i.eq.2) then
C...omega_tech0 -> Z0 + pi_tech0 not known.
            wdtp(i)=0d0
            wid2=wids(23,2)*wids(51,2)
          else
C...omega_tech0 -> f + fbar.
            wid2=1d0
            if(i.le.10) then
              ia=i-2
              fcof=3d0*radc
              if(ia.ge.6.and.ia.le.8) wid2=wids(ia,1)
            else
              ia=i
              fcof=1d0
              if(ia.ge.17) wid2=wids(ia,1)
            endif
            ei=kchg(ia,1)/3d0
            ai=sign(1d0,ei+0.1d0)
            vi=ai-4d0*ei*xwv
            vali=0.5d0*(vi+ai)
            vari=0.5d0*(vi-ai)
            wdtp(i)=facf*fcof*(1d0-rm1)*sqrt(max(0d0,1d0-4d0*rm1))*
     &      ((ei-vali*bwzr)**2+(vali*bwzi)**2+
     &      (ei-vari*bwzr)**2+(vari*bwzi)**2)
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  370   continue
 
      elseif(kfla.eq.kexcit+1) then
C...d* excited quark.
        fac=(sh/paru(155)**2)*shr
        do 380 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 380
          rm1=pmas(jamcomp(kfdp(idc,1)),1)**2/sh
          rm2=pmas(jamcomp(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 380
          if(i.eq.1) then
C...d* -> g + d.
            wdtp(i)=fac*as*paru(159)**2/3d0
            wid2=1d0
          elseif(i.eq.2) then
C...d* -> gamma + d.
            qf=-paru(157)/2d0+paru(158)/6d0
            wdtp(i)=fac*aem*qf**2/4d0
            wid2=1d0
          elseif(i.eq.3) then
C...d* -> Z0 + d.
            qf=-paru(157)*xw1/2d0-paru(158)*xw/6d0
            wdtp(i)=fac*aem*qf**2/(8d0*xw*xw1)*
     &      (1d0-rm1)**2*(2d0+rm1)
            wid2=wids(23,2)
          elseif(i.eq.4) then
C...d* -> W- + u.
            wdtp(i)=fac*aem*paru(157)**2/(16d0*xw)*
     &      (1d0-rm1)**2*(2d0+rm1)
            if(kflr.gt.0) wid2=wids(24,3)
            if(kflr.lt.0) wid2=wids(24,2)
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  380   continue
 
      elseif(kfla.eq.kexcit+2) then
C...u* excited quark.
        fac=(sh/paru(155)**2)*shr
        do 390 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 390
          rm1=pmas(jamcomp(kfdp(idc,1)),1)**2/sh
          rm2=pmas(jamcomp(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 390
          if(i.eq.1) then
C...u* -> g + u.
            wdtp(i)=fac*as*paru(159)**2/3d0
            wid2=1d0
          elseif(i.eq.2) then
C...u* -> gamma + u.
            qf=paru(157)/2d0+paru(158)/6d0
            wdtp(i)=fac*aem*qf**2/4d0
            wid2=1d0
          elseif(i.eq.3) then
C...u* -> Z0 + u.
            qf=paru(157)*xw1/2d0-paru(158)*xw/6d0
            wdtp(i)=fac*aem*qf**2/(8d0*xw*xw1)*
     &      (1d0-rm1)**2*(2d0+rm1)
            wid2=wids(23,2)
          elseif(i.eq.4) then
C...u* -> W+ + d.
            wdtp(i)=fac*aem*paru(157)**2/(16d0*xw)*
     &      (1d0-rm1)**2*(2d0+rm1)
            if(kflr.gt.0) wid2=wids(24,2)
            if(kflr.lt.0) wid2=wids(24,3)
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  390   continue
 
      elseif(kfla.eq.kexcit+11) then
C...e* excited lepton.
        fac=(sh/paru(155)**2)*shr
        do 400 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 400
          rm1=pmas(jamcomp(kfdp(idc,1)),1)**2/sh
          rm2=pmas(jamcomp(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 400
          if(i.eq.1) then
C...e* -> gamma + e.
            qf=-paru(157)/2d0-paru(158)/2d0
            wdtp(i)=fac*aem*qf**2/4d0
            wid2=1d0
          elseif(i.eq.2) then
C...e* -> Z0 + e.
            qf=-paru(157)*xw1/2d0+paru(158)*xw/2d0
            wdtp(i)=fac*aem*qf**2/(8d0*xw*xw1)*
     &      (1d0-rm1)**2*(2d0+rm1)
            wid2=wids(23,2)
          elseif(i.eq.3) then
C...e* -> W- + nu.
            wdtp(i)=fac*aem*paru(157)**2/(16d0*xw)*
     &      (1d0-rm1)**2*(2d0+rm1)
            if(kflr.gt.0) wid2=wids(24,3)
            if(kflr.lt.0) wid2=wids(24,2)
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  400   continue
 
      elseif(kfla.eq.kexcit+12) then
C...nu*_e excited neutrino.
        fac=(sh/paru(155)**2)*shr
        do 410 i=1,mdcy(kc,3)
          idc=i+mdcy(kc,2)-1
          if(mdme(idc,1).lt.0) goto 410
          rm1=pmas(jamcomp(kfdp(idc,1)),1)**2/sh
          rm2=pmas(jamcomp(kfdp(idc,2)),1)**2/sh
          if(sqrt(rm1)+sqrt(rm2).gt.1d0) goto 410
          if(i.eq.1) then
C...nu*_e -> Z0 + nu*_e.
            qf=paru(157)*xw1/2d0+paru(158)*xw/2d0
            wdtp(i)=fac*aem*qf**2/(8d0*xw*xw1)*
     &      (1d0-rm1)**2*(2d0+rm1)
            wid2=wids(23,2)
          elseif(i.eq.2) then
C...nu*_e -> W+ + e.
            wdtp(i)=fac*aem*paru(157)**2/(16d0*xw)*
     &      (1d0-rm1)**2*(2d0+rm1)
            if(kflr.gt.0) wid2=wids(24,2)
            if(kflr.lt.0) wid2=wids(24,3)
          endif
          wdtp(0)=wdtp(0)+wdtp(i)
          if(mdme(idc,1).gt.0) then
            wdte(i,mdme(idc,1))=wdtp(i)*wid2
            wdte(0,mdme(idc,1))=wdte(0,mdme(idc,1))+wdte(i,mdme(idc,1))
            wdte(i,0)=wdte(i,mdme(idc,1))
            wdte(0,0)=wdte(0,0)+wdte(i,0)
          endif
  410   continue
 
      endif
      mint(61)=0
      mint(62)=0
      mint(63)=0
 
      return
      end
 
C***********************************************************************
 
C...PYOFSH
C...Calculates partial width and differential cross-section maxima
C...of channels/processes not allowed on mass-shell, and selects
C...masses in such channels/processes.
 
      subroutine pjofsh(mofsh,kfmo,kfd1,kfd2,pmmo,ret1,ret2)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/jydat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      common/pjint5/ngenpd,ngen(0:500,3),xsec(0:500,3)
      save /jydat1/,/jydat2/,/jydat3/,/pjsubs/,/pjpars/,/pjint1/,
     &/pjint2/,/pjint5/
C...Local arrays.
      dimension kfd(2),mbw(2),pmd(2),pgd(2),pmg(2),pml(2),pmu(2),
     &pmh(2),atl(2),atu(2),ath(2),rmg(2),inx1(100),xpt1(100),
     &fpt1(100),inx2(100),xpt2(100),fpt2(100)
 
C...Find if particles equal, maximum mass, matrix elements, etc.
      mint(51)=0
      isub=mint(1)
      kfd(1)=iabs(kfd1)
      kfd(2)=iabs(kfd2)
      meql=0
      if(kfd(1).eq.kfd(2)) meql=1
      mlm=0
      if(mofsh.ge.2.and.meql.eq.1) mlm=int(1.5d0+pjr(0))
      if(mofsh.le.2.or.mofsh.eq.5) then
        noff=44
        pmmx=pmmo
      else
        noff=40
        pmmx=vint(1)
        if(ckin(2).gt.ckin(1)) pmmx=min(ckin(2),vint(1))
      endif
      mmed=0
      if((kfmo.eq.25.or.kfmo.eq.35.or.kfmo.eq.36).and.meql.eq.1.and.
     &(kfd(1).eq.23.or.kfd(1).eq.24)) mmed=1
      if((kfmo.eq.32.or.iabs(kfmo).eq.34).and.(kfd(1).eq.23.or.
     &kfd(1).eq.24).and.(kfd(2).eq.23.or.kfd(2).eq.24)) mmed=2
      if((kfmo.eq.32.or.iabs(kfmo).eq.34).and.(kfd(2).eq.25.or.
     &kfd(2).eq.35.or.kfd(2).eq.36)) mmed=3
      loop=1
 
C...Find where Breit-Wigners are required, else select discrete masses.
  100 do 110 i=1,2
        kfca=jamcomp(kfd(i))
        if(kfca.gt.0) then
          pmd(i)=pmas(kfca,1)
          pgd(i)=pmas(kfca,2)
        else
          pmd(i)=0d0
          pgd(i)=0d0
        endif
        if(mstp(42).le.0.or.pgd(i).lt.parp(41)) then
          mbw(i)=0
          pmg(i)=pmd(i)
          rmg(i)=(pmg(i)/pmmx)**2
        else
          mbw(i)=1
        endif
  110 continue
 
C...Find allowed mass range and Breit-Wigner parameters.
      do 120 i=1,2
        if(mofsh.eq.1.and.loop.eq.1.and.mbw(i).eq.1) then
          pml(i)=parp(42)
          pmu(i)=pmmx-parp(42)
          if(mbw(3-i).eq.0) pmu(i)=min(pmu(i),pmmx-pmd(3-i))
          if(pmu(i).lt.pml(i)+parj(64)) mbw(i)=-1
        elseif(mbw(i).eq.1.and.mofsh.ne.5) then
          ilm=i
          if(mlm.eq.2) ilm=3-i
          pml(i)=max(ckin(noff+2*ilm-1),parp(42))
          pmu(i)=pmmx-max(ckin(noff+5-2*ilm),parp(42))
          if(ckin(noff+2*ilm).gt.ckin(noff+2*ilm-1)) pmu(i)=min(pmu(i),
     &    ckin(noff+2*ilm))
          if(mbw(3-i).eq.0) pmu(i)=min(pmu(i),pmmx-pmd(3-i))
          if(i.eq.mlm) pmu(i)=min(pmu(i),0.5d0*pmmx)
          if(meql.eq.0) pmh(i)=min(pmu(i),0.5d0*pmmx)
          if(pmu(i).lt.pml(i)+parj(64)) mbw(i)=-1
          if(mbw(i).eq.1) then
            atl(i)=atan((pml(i)**2-pmd(i)**2)/(pmd(i)*pgd(i)))
            atu(i)=atan((pmu(i)**2-pmd(i)**2)/(pmd(i)*pgd(i)))
            if(meql.eq.0) ath(i)=atan((pmh(i)**2-pmd(i)**2)/(pmd(i)*
     &      pgd(i)))
          endif
        elseif(mbw(i).eq.1.and.mofsh.eq.5) then
          ilm=i
          if(mlm.eq.2) ilm=3-i
          pml(i)=max(ckin(48+i),parp(42))
          pmu(i)=pmmx-max(ckin(51-i),parp(42))
          if(mbw(3-i).eq.0) pmu(i)=min(pmu(i),pmmx-pmd(3-i))
          if(i.eq.mlm) pmu(i)=min(pmu(i),0.5d0*pmmx)
          if(meql.eq.0) pmh(i)=min(pmu(i),0.5d0*pmmx)
          if(pmu(i).lt.pml(i)+parj(64)) mbw(i)=-1
          if(mbw(i).eq.1) then
            atl(i)=atan((pml(i)**2-pmd(i)**2)/(pmd(i)*pgd(i)))
            atu(i)=atan((pmu(i)**2-pmd(i)**2)/(pmd(i)*pgd(i)))
            if(meql.eq.0) ath(i)=atan((pmh(i)**2-pmd(i)**2)/(pmd(i)*
     &      pgd(i)))
          endif
        endif
  120 continue
      if(mbw(1).lt.0.or.mbw(2).lt.0.or.(mbw(1).eq.0.and.mbw(2).eq.0))
     &then
        call pjerrm(3,'(PYOFSH:) no allowed decay product masses')
        mint(51)=1
        return
      endif
 
C...Calculation of partial width of resonance.
      if(mofsh.eq.1) then
 
C..If only one integration, pick that to be the inner.
        if(mbw(1).eq.0) then
          pm2=pmd(1)
          pmd(1)=pmd(2)
          pgd(1)=pgd(2)
          pml(1)=pml(2)
          pmu(1)=pmu(2)
        elseif(mbw(2).eq.0) then
          pm2=pmd(2)
        endif
 
C...Start outer loop of integration.
        if(mbw(1).eq.1.and.mbw(2).eq.1) then
          atl2=atan((pml(2)**2-pmd(2)**2)/(pmd(2)*pgd(2)))
          atu2=atan((pmu(2)**2-pmd(2)**2)/(pmd(2)*pgd(2)))
          npt2=1
          xpt2(1)=1d0
          inx2(1)=0
          fmax2=0d0
        endif
  130   if(mbw(1).eq.1.and.mbw(2).eq.1) then
          pm2s=pmd(2)**2+pmd(2)*pgd(2)*tan(atl2+xpt2(npt2)*(atu2-atl2))
          pm2=min(pmu(2),max(pml(2),sqrt(max(0d0,pm2s))))
        endif
        rm2=(pm2/pmmx)**2
 
C...Start inner loop of integration.
        pml1=pml(1)
        pmu1=min(pmu(1),pmmx-pm2)
        if(meql.eq.1) pmu1=min(pmu1,pm2)
        atl1=atan((pml1**2-pmd(1)**2)/(pmd(1)*pgd(1)))
        atu1=atan((pmu1**2-pmd(1)**2)/(pmd(1)*pgd(1)))
        if(pml1+parj(64).ge.pmu1.or.atl1+1d-7.ge.atu1) then
          func2=0d0
          goto 180
        endif
        npt1=1
        xpt1(1)=1d0
        inx1(1)=0
        fmax1=0d0
  140   pm1s=pmd(1)**2+pmd(1)*pgd(1)*tan(atl1+xpt1(npt1)*(atu1-atl1))
        pm1=min(pmu1,max(pml1,sqrt(max(0d0,pm1s))))
        rm1=(pm1/pmmx)**2
 
C...Evaluate function value - inner loop.
        func1=sqrt(max(0d0,(1d0-rm1-rm2)**2-4d0*rm1*rm2))
        if(mmed.eq.1) func1=func1*((1d0-rm1-rm2)**2+8d0*rm1*rm2)
        if(mmed.eq.2) func1=func1**3*(1d0+10d0*rm1+10d0*rm2+rm1**2+
     &  rm2**2+10d0*rm1*rm2)
        if(func1.gt.fmax1) fmax1=func1
        fpt1(npt1)=func1
 
C...Go to next position in inner loop.
        if(npt1.eq.1) then
          npt1=npt1+1
          xpt1(npt1)=0d0
          inx1(npt1)=1
          goto 140
        elseif(npt1.le.8) then
          npt1=npt1+1
          if(npt1.le.4.or.npt1.eq.6) ish1=1
          ish1=ish1+1
          xpt1(npt1)=0.5d0*(xpt1(ish1)+xpt1(inx1(ish1)))
          inx1(npt1)=inx1(ish1)
          inx1(ish1)=npt1
          goto 140
        elseif(npt1.lt.100) then
          isn1=ish1
  150     ish1=ish1+1
          if(ish1.gt.npt1) ish1=2
          if(ish1.eq.isn1) goto 160
          dfpt1=abs(fpt1(ish1)-fpt1(inx1(ish1)))
          if(dfpt1.lt.parp(43)*fmax1) goto 150
          npt1=npt1+1
          xpt1(npt1)=0.5d0*(xpt1(ish1)+xpt1(inx1(ish1)))
          inx1(npt1)=inx1(ish1)
          inx1(ish1)=npt1
          goto 140
        endif
 
C...Calculate integral over inner loop.
  160   fsum1=0d0
        do 170 ipt1=2,npt1
          fsum1=fsum1+0.5d0*(fpt1(ipt1)+fpt1(inx1(ipt1)))*
     &    (xpt1(inx1(ipt1))-xpt1(ipt1))
  170   continue
        func2=fsum1*(atu1-atl1)/paru(1)
  180   if(mbw(1).eq.1.and.mbw(2).eq.1) then
          if(func2.gt.fmax2) fmax2=func2
          fpt2(npt2)=func2
 
C...Go to next position in outer loop.
          if(npt2.eq.1) then
            npt2=npt2+1
            xpt2(npt2)=0d0
            inx2(npt2)=1
            goto 130
          elseif(npt2.le.8) then
            npt2=npt2+1
            if(npt2.le.4.or.npt2.eq.6) ish2=1
            ish2=ish2+1
            xpt2(npt2)=0.5d0*(xpt2(ish2)+xpt2(inx2(ish2)))
            inx2(npt2)=inx2(ish2)
            inx2(ish2)=npt2
            goto 130
          elseif(npt2.lt.100) then
            isn2=ish2
  190       ish2=ish2+1
            if(ish2.gt.npt2) ish2=2
            if(ish2.eq.isn2) goto 200
            dfpt2=abs(fpt2(ish2)-fpt2(inx2(ish2)))
            if(dfpt2.lt.parp(43)*fmax2) goto 190
            npt2=npt2+1
            xpt2(npt2)=0.5d0*(xpt2(ish2)+xpt2(inx2(ish2)))
            inx2(npt2)=inx2(ish2)
            inx2(ish2)=npt2
            goto 130
          endif
 
C...Calculate integral over outer loop.
  200     fsum2=0d0
          do 210 ipt2=2,npt2
            fsum2=fsum2+0.5d0*(fpt2(ipt2)+fpt2(inx2(ipt2)))*
     &      (xpt2(inx2(ipt2))-xpt2(ipt2))
  210     continue
          fsum2=fsum2*(atu2-atl2)/paru(1)
          if(meql.eq.1) fsum2=2d0*fsum2
        else
          fsum2=func2
        endif
 
C...Save result; second integration for user-selected mass range.
        if(loop.eq.1) widw=fsum2
        wid2=fsum2
        if(loop.eq.1.and.(ckin(46).ge.ckin(45).or.ckin(48).ge.ckin(47)
     &  .or.max(ckin(45),ckin(47)).ge.1.01d0*parp(42))) then
          loop=2
          goto 100
        endif
        ret1=widw
        ret2=wid2/widw
 
C...Select two decay product masses of a resonance.
      elseif(mofsh.eq.2.or.mofsh.eq.5) then
  220   do 230 i=1,2
          if(mbw(i).eq.0) goto 230
          pmbw=pmd(i)**2+pmd(i)*pgd(i)*tan(atl(i)+pjr(0)*
     &    (atu(i)-atl(i)))
          pmg(i)=min(pmu(i),max(pml(i),sqrt(max(0d0,pmbw))))
          rmg(i)=(pmg(i)/pmmx)**2
  230   continue
        if((meql.eq.1.and.pmg(max(1,mlm)).gt.pmg(min(2,3-mlm))).or.
     &  pmg(1)+pmg(2)+parj(64).gt.pmmx) goto 220
 
C...Weight with matrix element (if none known, use beta factor).
        flam=sqrt(max(0d0,(1d0-rmg(1)-rmg(2))**2-4d0*rmg(1)*rmg(2)))
        if(mmed.eq.1) then
          wtbe=flam*((1d0-rmg(1)-rmg(2))**2+8d0*rmg(1)*rmg(2))
        elseif(mmed.eq.2) then
          wtbe=flam**3*(1d0+10d0*rmg(1)+10d0*rmg(2)+rmg(1)**2+
     &    rmg(2)**2+10d0*rmg(1)*rmg(2))
        elseif(mmed.eq.3) then
          wtbe=flam*(rmg(1)+flam**2/12d0)
        else
          wtbe=flam
        endif
        if(wtbe.lt.pjr(0)) goto 220
        ret1=pmg(1)
        ret2=pmg(2)
 
C...Find suitable set of masses for initialization of 2 -> 2 processes.
      elseif(mofsh.eq.3) then
        if(mbw(1).ne.0.and.mbw(2).eq.0) then
          pmg(1)=min(pmd(1),0.5d0*(pml(1)+pmu(1)))
          pmg(2)=pmd(2)
        elseif(mbw(2).ne.0.and.mbw(1).eq.0) then
          pmg(1)=pmd(1)
          pmg(2)=min(pmd(2),0.5d0*(pml(2)+pmu(2)))
        else
          idiv=-1
  240     idiv=idiv+1
          pmg(1)=min(pmd(1),0.1d0*(idiv*pml(1)+(10-idiv)*pmu(1)))
          pmg(2)=min(pmd(2),0.1d0*(idiv*pml(2)+(10-idiv)*pmu(2)))
          if(idiv.le.9.and.pmg(1)+pmg(2).gt.0.9d0*pmmx) goto 240
        endif
        ret1=pmg(1)
        ret2=pmg(2)
 
C...Evaluate importance of excluded tails of Breit-Wigners.
        if(meql.eq.0.and.mbw(1).eq.1.and.mbw(2).eq.1.and.pmd(1)+pmd(2)
     &  .gt.pmmx.and.pmh(1).gt.pml(1).and.pmh(2).gt.pml(2)) meql=2
        if(meql.le.1) then
          vint(80)=1d0
          do 250 i=1,2
            if(mbw(i).ne.0) vint(80)=vint(80)*1.25d0*(atu(i)-atl(i))/
     &      paru(1)
  250     continue
        else
          vint(80)=(1.25d0/paru(1))**2*max((atu(1)-atl(1))*
     &    (ath(2)-atl(2)),(ath(1)-atl(1))*(atu(2)-atl(2)))
        endif
        if((isub.eq.15.or.isub.eq.19.or.isub.eq.30.or.isub.eq.35).and.
     &  mstp(43).ne.2) vint(80)=2d0*vint(80)
        if(isub.eq.22.and.mstp(43).ne.2) vint(80)=4d0*vint(80)
        if(meql.ge.1) vint(80)=2d0*vint(80)
 
C...Pick one particle to be the lighter (if improves efficiency).
      elseif(mofsh.eq.4) then
        if(meql.eq.0.and.mbw(1).eq.1.and.mbw(2).eq.1.and.pmd(1)+pmd(2)
     &  .gt.pmmx.and.pmh(1).gt.pml(1).and.pmh(2).gt.pml(2)) meql=2
  260   if(meql.eq.2) mlm=int(1.5d0+pjr(0))
 
C...Select two masses according to Breit-Wigner + flat in s + 1/s.
        do 270 i=1,2
          if(mbw(i).eq.0) goto 270
          pmv=pmu(i)
          if(meql.eq.2.and.i.eq.mlm) pmv=pmh(i)
          atv=atu(i)
          if(meql.eq.2.and.i.eq.mlm) atv=ath(i)
          rbr=pjr(0)
          if((isub.eq.15.or.isub.eq.19.or.isub.eq.22.or.isub.eq.30.or.
     &    isub.eq.35).and.mstp(43).ne.2) rbr=2d0*rbr
          if(rbr.lt.0.8d0) then
            pmsr=pmd(i)**2+pmd(i)*pgd(i)*tan(atl(i)+pjr(0)*(atv-atl(i)))
            pmg(i)=min(pmv,max(pml(i),sqrt(max(0d0,pmsr))))
          elseif(rbr.lt.0.9d0) then
            pmg(i)=sqrt(max(0d0,pml(i)**2+pjr(0)*(pmv**2-pml(i)**2)))
          elseif(rbr.lt.1.5d0) then
            pmg(i)=pml(i)*(pmv/pml(i))**pjr(0)
          else
            pmg(i)=sqrt(max(0d0,pml(i)**2*pmv**2/(pml(i)**2+pjr(0)*
     &      (pmv**2-pml(i)**2))))
          endif
  270   continue
        if((meql.ge.1.and.pmg(max(1,mlm)).gt.pmg(min(2,3-mlm))).or.
     &  pmg(1)+pmg(2)+parj(64).gt.pmmx) then
          if(mint(48).eq.1) then
            ngen(0,1)=ngen(0,1)+1
            ngen(mint(1),1)=ngen(mint(1),1)+1
            goto 260
          else
            mint(51)=1
            return
          endif
        endif
        ret1=pmg(1)
        ret2=pmg(2)
 
C...Give weight for selected mass distribution.
        vint(80)=1d0
        do 280 i=1,2
          if(mbw(i).eq.0) goto 280
          pmv=pmu(i)
          if(meql.eq.2.and.i.eq.mlm) pmv=pmh(i)
          atv=atu(i)
          if(meql.eq.2.and.i.eq.mlm) atv=ath(i)
          f0=pmd(i)*pgd(i)/((pmg(i)**2-pmd(i)**2)**2+
     &    (pmd(i)*pgd(i))**2)/paru(1)
          f1=1d0
          f2=1d0/pmg(i)**2
          f3=1d0/pmg(i)**4
          fi0=(atv-atl(i))/paru(1)
          fi1=pmv**2-pml(i)**2
          fi2=2d0*log(pmv/pml(i))
          fi3=1d0/pml(i)**2-1d0/pmv**2
          if((isub.eq.15.or.isub.eq.19.or.isub.eq.22.or.isub.eq.30.or.
     &    isub.eq.35).and.mstp(43).ne.2) then
            vint(80)=vint(80)*20d0/(8d0+(fi0/f0)*(f1/fi1+6d0*f2/fi2+
     &      5d0*f3/fi3))
          else
            vint(80)=vint(80)*10d0/(8d0+(fi0/f0)*(f1/fi1+f2/fi2))
          endif
          vint(80)=vint(80)*fi0
  280   continue
        if(meql.ge.1) vint(80)=2d0*vint(80)
      endif
 
      return
      end
 
C***********************************************************************
 
C...PYRECO
C...Handles the possibility of colour reconnection in W+W- events,
C...Based on the main scenarios of the Sjostrand and Khoze study:
C...I, II, II', intermediate and instantaneous; plus one model
C...along the lines of the Gustafson and Hakkinen: GH.
 
      subroutine pjreco(iw1,iw2,nsd1,naft1)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Parameter value; number of points in MC integration.
      parameter (npt=100)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      save /jyjets/,/jydat1/,/jydat2/,/pjpars/,/pjint1/
C...Local arrays.
      dimension nbeg(2),nend(2),inp(50),inm(50),beww(3),xp(3),xm(3),
     &v1(3),v2(3),betp(50,4),dirp(50,3),betm(50,4),dirm(50,3),
     &xd(4),xb(4),iap(npt),iam(npt),wta(npt),v1p(3),v2p(3),v1m(3),
     &v2m(3),q(4,3),xpp(3),xmm(3),ipc(20),imc(20),tc(0:20),tpc(20),
     &tmc(20),ijoin(100)
 
C...Functions to give four-product and to do determinants.
      four(i,j)=p(i,4)*p(j,4)-p(i,1)*p(j,1)-p(i,2)*p(j,2)-p(i,3)*p(j,3)
      deter(i,j,l)=q(i,1)*q(j,2)*q(l,3)-q(i,1)*q(l,2)*q(j,3)+
     &q(j,1)*q(l,2)*q(i,3)-q(j,1)*q(i,2)*q(l,3)+
     &q(l,1)*q(i,2)*q(j,3)-q(l,1)*q(j,2)*q(i,3)
 
C...Only allow fraction of recoupling for GH, intermediate and
C...instantaneous.
      if(mstp(115).eq.5.or.mstp(115).eq.11.or.mstp(115).eq.12) then
        if(pjr(0).gt.parp(120)) return
      endif
 
C...Common part for scenarios I, II, II', and GH.
      if(mstp(115).eq.1.or.mstp(115).eq.2.or.mstp(115).eq.3.or.
     &mstp(115).eq.5) then
 
C...Read out frequently-used parameters.
        pi=paru(1)
        hbar=paru(3)
        pmw=pmas(24,1)
        pgw=pmas(24,2)
        tfrag=parp(115)
        rhad=parp(116)
        fact=parp(117)
        blowr=parp(118)
        blowt=parp(119)
 
C...Find range of decay products of the W's.
C...Background: the W's are stored in IW1 and IW2.
C...Their direct decay products in NSD1+1 through NSD1+4.
C...Products after shower (if any) in NSD1+5 through NAFT1
C...for first W and in NAFT1+1 through N for the second.
        if(k(iw1,2).gt.0) then
          jt=1
        else
          jt=2
        endif
        jr=3-jt
        if(naft1.gt.nsd1+4) then
          nbeg(jt)=nsd1+5
          nend(jt)=naft1
        else
          nbeg(jt)=nsd1+1
          nend(jt)=nsd1+2
        endif
        if(n.gt.naft1) then
          nbeg(jr)=naft1+1
          nend(jr)=n
        else
          nbeg(jr)=nsd1+3
          nend(jr)=nsd1+4
        endif
 
C...Rearrange parton shower products along strings.
        nold=n
        call pjprep(nsd1+1)
 
C...Find partons pointing back to W+ and W-; store them with quark
C...end of string first.
        nnp=0
        nnm=0
        isgp=0
        isgm=0
        do 120 i=nold+1,n
          if(k(i,1).ne.1.and.k(i,1).ne.2) goto 120
          if(iabs(k(i,2)).ge.22) goto 120
          if(k(i,3).ge.nbeg(1).and.k(i,3).le.nend(1)) then
            if(isgp.eq.0) isgp=isign(1,k(i,2))
            nnp=nnp+1
            if(isgp.eq.1) then
              inp(nnp)=i
            else
              do 100 i1=nnp,2,-1
                inp(i1)=inp(i1-1)
  100         continue
              inp(1)=i
            endif
            if(k(i,1).eq.1) isgp=0
          elseif(k(i,3).ge.nbeg(2).and.k(i,3).le.nend(2)) then
            if(isgm.eq.0) isgm=isign(1,k(i,2))
            nnm=nnm+1
            if(isgm.eq.1) then
              inm(nnm)=i
            else
              do 110 i1=nnm,2,-1
                inm(i1)=inm(i1-1)
  110         continue
              inm(1)=i
            endif
            if(k(i,1).eq.1) isgm=0
          endif
  120   continue
 
C...Boost to W+W- rest frame (not strictly needed).
        do 130 j=1,3
          beww(j)=(p(iw1,j)+p(iw2,j))/(p(iw1,4)+p(iw2,4))
  130   continue
        call pjrobo(iw1,iw1,0d0,0d0,-beww(1),-beww(2),-beww(3))
        call pjrobo(iw2,iw2,0d0,0d0,-beww(1),-beww(2),-beww(3))
        call pjrobo(nold+1,n,0d0,0d0,-beww(1),-beww(2),-beww(3))
 
C...Select decay vertices of W+ and W-.
        tp=hbar*(-log(pjr(0)))*p(iw1,4)/
     &  sqrt((p(iw1,5)**2-pmw**2)**2+(p(iw1,5)**2*pgw/pmw)**2)
        tm=hbar*(-log(pjr(0)))*p(iw2,4)/
     &  sqrt((p(iw2,5)**2-pmw**2)**2+(p(iw2,5)**2*pgw/pmw)**2)
        gtmax=max(tp,tm)
        do 140 j=1,3
          xp(j)=tp*p(iw1,j)/p(iw1,4)
          xm(j)=tm*p(iw2,j)/p(iw2,4)
  140   continue
 
C...Begin scenario I specifics.
        if(mstp(115).eq.1) then
 
C...Reconstruct velocity and direction of W+ string pieces.
          do 170 iip=1,nnp-1
            if(k(inp(iip),2).lt.0) goto 170
            i1=inp(iip)
            i2=inp(iip+1)
            p1a=sqrt(p(i1,1)**2+p(i1,2)**2+p(i1,3)**2)
            p2a=sqrt(p(i2,1)**2+p(i2,2)**2+p(i2,3)**2)
            do 150 j=1,3
              v1(j)=p(i1,j)/p1a
              v2(j)=p(i2,j)/p2a
              betp(iip,j)=0.5d0*(v1(j)+v2(j))
              dirp(iip,j)=v1(j)-v2(j)
  150       continue
            betp(iip,4)=1d0/sqrt(1d0-betp(iip,1)**2-betp(iip,2)**2-
     &      betp(iip,3)**2)
            dirl=sqrt(dirp(iip,1)**2+dirp(iip,2)**2+dirp(iip,3)**2)
            do 160 j=1,3
              dirp(iip,j)=dirp(iip,j)/dirl
  160       continue
  170     continue
 
C...Reconstruct velocity and direction of W- string pieces.
          do 200 iim=1,nnm-1
            if(k(inm(iim),2).lt.0) goto 200
            i1=inm(iim)
            i2=inm(iim+1)
            p1a=sqrt(p(i1,1)**2+p(i1,2)**2+p(i1,3)**2)
            p2a=sqrt(p(i2,1)**2+p(i2,2)**2+p(i2,3)**2)
            do 180 j=1,3
              v1(j)=p(i1,j)/p1a
              v2(j)=p(i2,j)/p2a
              betm(iim,j)=0.5d0*(v1(j)+v2(j))
              dirm(iim,j)=v1(j)-v2(j)
  180       continue
            betm(iim,4)=1d0/sqrt(1d0-betm(iim,1)**2-betm(iim,2)**2-
     &      betm(iim,3)**2)
            dirl=sqrt(dirm(iim,1)**2+dirm(iim,2)**2+dirm(iim,3)**2)
            do 190 j=1,3
              dirm(iim,j)=dirm(iim,j)/dirl
  190       continue
  200     continue
 
C...Loop over number of space-time points.
          nacc=0
          sum=0d0
          do 250 ipt=1,npt
 
C...Pick x,y,z,t Gaussian (width RHAD and TFRAG, respectively).
            r=sqrt(-log(pjr(0)))
            phi=2d0*pi*pjr(0)
            x=blowr*rhad*r*cos(phi)
            y=blowr*rhad*r*sin(phi)
            r=sqrt(-log(pjr(0)))
            phi=2d0*pi*pjr(0)
            z=blowr*rhad*r*cos(phi)
            t=gtmax+blowt*sqrt(0.5d0)*tfrag*r*abs(sin(phi))
 
C...Weight for sample distribution.
            wtsmp=exp(-(x**2+y**2+z**2)/(blowr*rhad)**2)*
     &      exp(-2d0*(t-gtmax)**2/(blowt*tfrag)**2)
 
C...Loop over W+ string pieces and find one with largest weight.
            imaxp=0
            wtmaxp=1d-10
            xd(1)=x-xp(1)
            xd(2)=y-xp(2)
            xd(3)=z-xp(3)
            xd(4)=t-tp
            do 220 iip=1,nnp-1
              if(k(inp(iip),2).lt.0) goto 220
              bed=betp(iip,1)*xd(1)+betp(iip,2)*xd(2)+betp(iip,3)*xd(3)
              bedg=betp(iip,4)*(betp(iip,4)*bed/(1d0+betp(iip,4))-xd(4))
              do 210 j=1,3
                xb(j)=xd(j)+bedg*betp(iip,j)
  210         continue
              xb(4)=betp(iip,4)*(xd(4)-bed)
              sr2=xb(1)**2+xb(2)**2+xb(3)**2
              sz2=(dirp(iip,1)*xb(1)+dirp(iip,2)*xb(2)+
     &        dirp(iip,3)*xb(3))**2
              wtp=exp(-(sr2-sz2)/(2d0*rhad**2))*exp(-(xb(4)**2-sz2)/
     &        tfrag**2)
              if(xb(4)-sqrt(sr2).lt.0d0) wtp=0d0
              if(wtp.gt.wtmaxp) then
                imaxp=iip
                wtmaxp=wtp
              endif
  220       continue
 
C...Loop over W- string pieces and find one with largest weight.
            imaxm=0
            wtmaxm=1d-10
            xd(1)=x-xm(1)
            xd(2)=y-xm(2)
            xd(3)=z-xm(3)
            xd(4)=t-tm
            do 240 iim=1,nnm-1
              if(k(inm(iim),2).lt.0) goto 240
              bed=betm(iim,1)*xd(1)+betm(iim,2)*xd(2)+betm(iim,3)*xd(3)
              bedg=betm(iim,4)*(betm(iim,4)*bed/(1d0+betm(iim,4))-xd(4))
              do 230 j=1,3
                xb(j)=xd(j)+bedg*betm(iim,j)
  230         continue
              xb(4)=betm(iim,4)*(xd(4)-bed)
              sr2=xb(1)**2+xb(2)**2+xb(3)**2
              sz2=(dirm(iim,1)*xb(1)+dirm(iim,2)*xb(2)+
     &        dirm(iim,3)*xb(3))**2
              wtm=exp(-(sr2-sz2)/(2d0*rhad**2))*exp(-(xb(4)**2-sz2)/
     &        tfrag**2)
              if(xb(4)-sqrt(sr2).lt.0d0) wtm=0d0
              if(wtm.gt.wtmaxm) then
                imaxm=iim
                wtmaxm=wtm
              endif
  240       continue
 
C...Result of integration.
            wt=0d0
            if(imaxp.ne.0.and.imaxm.ne.0) then
              wt=wtmaxp*wtmaxm/wtsmp
              sum=sum+wt
              nacc=nacc+1
              iap(nacc)=imaxp
              iam(nacc)=imaxm
              wta(nacc)=wt
            endif
  250     continue
          res=blowr**3*blowt*sum/npt
 
C...Decide whether to reconnect and, if so, where.
          iacc=0
          prec=1d0-exp(-fact*res)
          if(prec.gt.pjr(0)) then
            rsum=pjr(0)*sum
            do 260 ia=1,nacc
              iacc=ia
              rsum=rsum-wta(ia)
              if(rsum.le.0d0) goto 270
  260       continue
  270       iip=iap(iacc)
            iim=iam(iacc)
          endif
 
C...Begin scenario II and II' specifics.
        elseif(mstp(115).eq.2.or.mstp(115).eq.3) then
 
C...Loop through all string pieces, one from W+ and one from W-.
          ncross=0
          tc(0)=0d0
          do 340 iip=1,nnp-1
            if(k(inp(iip),2).lt.0) goto 340
            i1p=inp(iip)
            i2p=inp(iip+1)
            do 330 iim=1,nnm-1
              if(k(inm(iim),2).lt.0) goto 330
              i1m=inm(iim)
              i2m=inm(iim+1)
 
C...Find endpoint velocity vectors.
              do 280 j=1,3
                v1p(j)=p(i1p,j)/p(i1p,4)
                v2p(j)=p(i2p,j)/p(i2p,4)
                v1m(j)=p(i1m,j)/p(i1m,4)
                v2m(j)=p(i2m,j)/p(i2m,4)
  280         continue
 
C...Define q matrix and find t.
              do 290 j=1,3
                q(1,j)=v2p(j)-v1p(j)
                q(2,j)=-(v2m(j)-v1m(j))
                q(3,j)=xp(j)-xm(j)-tp*v1p(j)+tm*v1m(j)
                q(4,j)=v1p(j)-v1m(j)
  290         continue
              t=-deter(1,2,3)/deter(1,2,4)
 
C...Find alpha and beta; i.e. coordinates of crossing point.
              s11=q(1,1)*(t-tp)
              s12=q(2,1)*(t-tm)
              s13=q(3,1)+q(4,1)*t
              s21=q(1,2)*(t-tp)
              s22=q(2,2)*(t-tm)
              s23=q(3,2)+q(4,2)*t
              den=s11*s22-s12*s21
              alp=(s12*s23-s22*s13)/den
              bet=(s21*s13-s11*s23)/den
 
C...Check if solution acceptable.
              iansw=1
              if(t.lt.gtmax) iansw=0
              if(alp.lt.0d0.or.alp.gt.1d0) iansw=0
              if(bet.lt.0d0.or.bet.gt.1d0) iansw=0
 
C...Find point of crossing and check that not inconsistent.
              do 300 j=1,3
                xpp(j)=xp(j)+(v1p(j)+alp*(v2p(j)-v1p(j)))*(t-tp)
                xmm(j)=xm(j)+(v1m(j)+bet*(v2m(j)-v1m(j)))*(t-tm)
  300         continue
              d2pm=(xpp(1)-xmm(1))**2+(xpp(2)-xmm(2))**2+
     &        (xpp(3)-xmm(3))**2
              d2p=xpp(1)**2+xpp(2)**2+xpp(3)**2
              d2m=xmm(1)**2+xmm(2)**2+xmm(3)**2
              if(d2pm.gt.1d-4*(d2p+d2m)) iansw=-1
 
C...Find string eigentimes at crossing.
              if(iansw.eq.1) then
                taup=sqrt(max(0d0,(t-tp)**2-(xpp(1)-xp(1))**2-
     &          (xpp(2)-xp(2))**2-(xpp(3)-xp(3))**2))
                taum=sqrt(max(0d0,(t-tm)**2-(xmm(1)-xm(1))**2-
     &          (xmm(2)-xm(2))**2-(xmm(3)-xm(3))**2))
              else
                taup=0d0
                taum=0d0
              endif
 
C...Order crossings by time. End loop over crossings.
              if(iansw.eq.1.and.ncross.lt.20) then
                ncross=ncross+1
                do 310 i1=ncross,1,-1
                  if(t.gt.tc(i1-1).or.i1.eq.1) then
                    ipc(i1)=iip
                    imc(i1)=iim
                    tc(i1)=t
                    tpc(i1)=taup
                    tmc(i1)=taum
                    goto 320
                  else
                    ipc(i1)=ipc(i1-1)
                    imc(i1)=imc(i1-1)
                    tc(i1)=tc(i1-1)
                    tpc(i1)=tpc(i1-1)
                    tmc(i1)=tmc(i1-1)
                  endif
  310           continue
  320           continue
              endif
  330       continue
  340     continue
 
C...Loop over crossings; find first (if any) acceptable one.
          iacc=0
          if(ncross.ge.1) then
            do 350 ic=1,ncross
              pnfrag=exp(-(tpc(ic)**2+tmc(ic)**2)/tfrag**2)
              if(pnfrag.gt.pjr(0)) then
C...Scenario II: only compare with fragmentation time.
                if(mstp(115).eq.2) then
                  iacc=ic
                  iip=ipc(iacc)
                  iim=imc(iacc)
                  goto 360
C...Scenario II': also require that string length decreases.
                else
                  iip=ipc(ic)
                  iim=imc(ic)
                  i1p=inp(iip)
                  i2p=inp(iip+1)
                  i1m=inm(iim)
                  i2m=inm(iim+1)
                  elold=four(i1p,i2p)*four(i1m,i2m)
                  elnew=four(i1p,i2m)*four(i1m,i2p)
                  if(elnew.lt.elold) then
                    iacc=ic
                    iip=ipc(iacc)
                    iim=imc(iacc)
                    goto 360
                  endif
                endif
              endif
  350       continue
  360       continue
          endif
 
C...Begin scenario GH specifics.
        elseif(mstp(115).eq.5) then
 
C...Loop through all string pieces, one from W+ and one from W-.
          iacc=0
          elmin=1d0
          do 380 iip=1,nnp-1
            if(k(inp(iip),2).lt.0) goto 380
            i1p=inp(iip)
            i2p=inp(iip+1)
            do 370 iim=1,nnm-1
              if(k(inm(iim),2).lt.0) goto 370
              i1m=inm(iim)
              i2m=inm(iim+1)
 
C...Look for largest decrease of (exponent of) Lambda measure.
              elold=four(i1p,i2p)*four(i1m,i2m)
              elnew=four(i1p,i2m)*four(i1m,i2p)
              eldif=elnew/max(1d-10,elold)
              if(eldif.lt.elmin) then
                iacc=iip+iim
                elmin=eldif
                ipc(1)=iip
                imc(1)=iim
              endif
  370       continue
  380     continue
          iip=ipc(1)
          iim=imc(1)
        endif
 
C...Common for scenarios I, II, II' and GH: reconnect strings.
        if(iacc.ne.0) then
          mint(32)=1
          njoin=0
          do 390 is=1,nnp+nnm
            njoin=njoin+1
            if(is.le.iip) then
              i=inp(is)
            elseif(is.le.iip+nnm-iim) then
              i=inm(is-iip+iim)
            elseif(is.le.iip+nnm) then
              i=inm(is-iip-nnm+iim)
            else
              i=inp(is-nnm)
            endif
            ijoin(njoin)=i
            if(k(i,2).lt.0) then
              call pjjoin(njoin,ijoin)
              njoin=0
            endif
  390     continue
 
C...Restore original event record if no reconnection.
        else
          do 400 i=nsd1+1,nold
            if(k(i,1).eq.13.or.k(i,1).eq.14) then
              k(i,4)=mod(k(i,4),mstu(5)**2)
              k(i,5)=mod(k(i,5),mstu(5)**2)
            endif
  400     continue
          do 410 i=nold+1,n
            k(k(i,3),1)=3
  410     continue
          n=nold
        endif
 
C...Boost back system.
        call pjrobo(iw1,iw1,0d0,0d0,beww(1),beww(2),beww(3))
        call pjrobo(iw2,iw2,0d0,0d0,beww(1),beww(2),beww(3))
        if(n.gt.nold) call pjrobo(nold+1,n,0d0,0d0,
     &  beww(1),beww(2),beww(3))
 
C...Common part for intermediate and instantaneous scenarios.
      elseif(mstp(115).eq.11.or.mstp(115).eq.12) then
        mint(32)=1
 
C...Remove old shower products and reset showering ones.
        n=nsd1+4
        do 420 i=nsd1+1,nsd1+4
          k(i,1)=3
          k(i,4)=mod(k(i,4),mstu(5)**2)
          k(i,5)=mod(k(i,5),mstu(5)**2)
  420   continue
 
C...Identify quark-antiquark pairs.
        iq1=nsd1+1
        iq2=nsd1+2
        iq3=nsd1+3
        if(k(iq1,2)*k(iq3,2).lt.0) iq3=nsd1+4
        iq4=2*nsd1+7-iq3
 
C...Reconnect strings.
        ijoin(1)=iq1
        ijoin(2)=iq4
        call pjjoin(2,ijoin)
        ijoin(1)=iq3
        ijoin(2)=iq2
        call pjjoin(2,ijoin)
 
C...Do new parton showers in intermediate scenario.
        if(mstp(71).ge.1.and.mstp(115).eq.11) then
          mstj50=mstj(50)
          mstj(50)=0
          call pjshow(iq1,iq2,p(iw1,5))
          call pjshow(iq3,iq4,p(iw2,5))
          mstj(50)=mstj50
 
C...Do new parton showers in instantaneous scenario.
        elseif(mstp(71).ge.1.and.mstp(115).eq.12) then
          ppm2=(p(iq1,4)+p(iq4,4))**2-(p(iq1,1)+p(iq4,1))**2-
     &    (p(iq1,2)+p(iq4,2))**2-(p(iq1,3)+p(iq4,3))**2
          ppm=sqrt(max(0d0,ppm2))
          call pjshow(iq1,iq4,ppm)
          ppm2=(p(iq3,4)+p(iq2,4))**2-(p(iq3,1)+p(iq2,1))**2-
     &    (p(iq3,2)+p(iq2,2))**2-(p(iq3,3)+p(iq2,3))**2
          ppm=sqrt(max(0d0,ppm2))
          call pjshow(iq3,iq2,ppm)
        endif
      endif
 
      return
      end
 
C***********************************************************************
 
C...PYKLIM
C...Checks generated variables against pre-set kinematical limits;
C...also calculates limits on variables used in generation.
 
      subroutine pjklim(ilim)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/jydat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      save /jyjets/,/jydat1/,/jydat2/,/jydat3/,/pjsubs/,/pjpars/,
     &/pjint1/,/pjint2/
 
C...Common kinematical expressions.
      mint(51)=0
      isub=mint(1)
      istsb=iset(isub)
      if(isub.eq.96) goto 100
      sqm3=vint(63)
      sqm4=vint(64)
      if(ilim.ne.0) then
        if(abs(sqm3).lt.1d-4.and.abs(sqm4).lt.1d-4) then
          ckin09=max(ckin(9),ckin(13))
          ckin10=min(ckin(10),ckin(14))
          ckin11=max(ckin(11),ckin(15))
          ckin12=min(ckin(12),ckin(16))
        else
          ckin09=max(ckin(9),min(0d0,ckin(13)))
          ckin10=min(ckin(10),max(0d0,ckin(14)))
          ckin11=max(ckin(11),min(0d0,ckin(15)))
          ckin12=min(ckin(12),max(0d0,ckin(16)))
        endif
      endif
      if(ilim.ne.1) then
        tau=vint(21)
        rm3=sqm3/(tau*vint(2))
        rm4=sqm4/(tau*vint(2))
        be34=sqrt(max(1d-20,(1d0-rm3-rm4)**2-4d0*rm3*rm4))
      endif
      pthmin=ckin(3)
      if(min(sqm3,sqm4).lt.ckin(6)**2.and.istsb.ne.1.and.istsb.ne.3)
     &pthmin=max(ckin(3),ckin(5))
 
      if(ilim.eq.0) then
C...Check generated values of tau, y*, cos(theta-hat), and tau' against
C...pre-set kinematical limits.
        yst=vint(22)
        cth=vint(23)
        taup=vint(26)
        taue=tau
        if(istsb.ge.3.and.istsb.le.5) taue=taup
        x1=sqrt(taue)*exp(yst)
        x2=sqrt(taue)*exp(-yst)
        xf=x1-x2
        if(mint(47).ne.1) then
          if(tau*vint(2).lt.ckin(1)**2) mint(51)=1
          if(ckin(2).ge.0d0.and.tau*vint(2).gt.ckin(2)**2) mint(51)=1
          if(yst.lt.ckin(7).or.yst.gt.ckin(8)) mint(51)=1
          if(xf.lt.ckin(25).or.xf.gt.ckin(26)) mint(51)=1
        endif
        if(mint(45).ne.1) then
          if(x1.lt.ckin(21).or.x1.gt.ckin(22)) mint(51)=1
        endif
        if(mint(46).ne.1) then
          if(x2.lt.ckin(23).or.x2.gt.ckin(24)) mint(51)=1
        endif
        if(istsb.eq.2.or.istsb.eq.4) then
          pth=0.5d0*be34*sqrt(tau*vint(2)*max(0d0,1d0-cth**2))
          expy3=max(1.d-10,(1d0+rm3-rm4+be34*cth)/
     &    max(1.d-10,(1d0+rm3-rm4-be34*cth)))
          expy4=max(1.d-10,(1d0-rm3+rm4-be34*cth)/
     &    max(1.d-10,(1d0-rm3+rm4+be34*cth)))
          y3=yst+0.5d0*log(expy3)
          y4=yst+0.5d0*log(expy4)
          ylarge=max(y3,y4)
          ysmall=min(y3,y4)
          etalar=10d0
          etasma=-10d0
          sth=sqrt(max(0d0,1d0-cth**2))
          exsq3=sqrt(max(1d-20,((1d0+rm3-rm4)*cosh(yst)+be34*sinh(yst)*
     &    cth)**2-4d0*rm3))
          exsq4=sqrt(max(1d-20,((1d0-rm3+rm4)*cosh(yst)-be34*sinh(yst)*
     &    cth)**2-4d0*rm4))
          if(sth.ge.1.d-6) then
            expet3=((1d0+rm3-rm4)*sinh(yst)+be34*cosh(yst)*cth+exsq3)/
     &      (be34*sth)
            expet4=((1d0-rm3+rm4)*sinh(yst)-be34*cosh(yst)*cth+exsq4)/
     &      (be34*sth)
            eta3=log(min(1.d10,max(1.d-10,expet3)))
            eta4=log(min(1.d10,max(1.d-10,expet4)))
            etalar=max(eta3,eta4)
            etasma=min(eta3,eta4)
          endif
          cts3=((1d0+rm3-rm4)*sinh(yst)+be34*cosh(yst)*cth)/exsq3
          cts4=((1d0-rm3+rm4)*sinh(yst)-be34*cosh(yst)*cth)/exsq4
          ctslar=min(1d0,max(cts3,cts4))
          ctssma=max(-1d0,min(cts3,cts4))
          sh=tau*vint(2)
          rpts=4d0*vint(71)**2/sh
          be34l=sqrt(max(0d0,(1d0-rm3-rm4)**2-4d0*rm3*rm4-rpts))
          rm34=max(1d-20,2d0*rm3*rm4)
          if(2d0*vint(71)**2/(vint(21)*vint(2)).lt.0.0001d0)
     &    rm34=max(rm34,2d0*vint(71)**2/(vint(21)*vint(2)))
          rthm=(4d0*rm3*rm4+rpts)/(1d0-rm3-rm4+be34l)
          tha=0.5d0*sh*max(rthm,1d0-rm3-rm4-be34*cth)
          uha=0.5d0*sh*max(rthm,1d0-rm3-rm4+be34*cth)
          if(pth.lt.pthmin) mint(51)=1
          if(ckin(4).ge.0d0.and.pth.gt.ckin(4)) mint(51)=1
          if(ylarge.lt.ckin(9).or.ylarge.gt.ckin(10)) mint(51)=1
          if(ysmall.lt.ckin(11).or.ysmall.gt.ckin(12)) mint(51)=1
          if(etalar.lt.ckin(13).or.etalar.gt.ckin(14)) mint(51)=1
          if(etasma.lt.ckin(15).or.etasma.gt.ckin(16)) mint(51)=1
          if(ctslar.lt.ckin(17).or.ctslar.gt.ckin(18)) mint(51)=1
          if(ctssma.lt.ckin(19).or.ctssma.gt.ckin(20)) mint(51)=1
          if(cth.lt.ckin(27).or.cth.gt.ckin(28)) mint(51)=1
          if(tha.lt.ckin(35)) mint(51)=1
          if(ckin(36).ge.0d0.and.tha.gt.ckin(36)) mint(51)=1
          if(uha.lt.ckin(37)) mint(51)=1
          if(ckin(38).ge.0d0.and.uha.gt.ckin(38)) mint(51)=1
        endif
        if(istsb.ge.3.and.istsb.le.5) then
          if(taup*vint(2).lt.ckin(31)**2) mint(51)=1
          if(ckin(32).ge.0d0.and.taup*vint(2).gt.ckin(32)**2) mint(51)=1
        endif
 
C...Additional cuts on W2 (approximately) in DIS.
        if(isub.eq.10) then
          xbj=x2
          if(iabs(mint(12)).lt.20) xbj=x1
          q2bj=tha
          w2bj=q2bj*(1d0-xbj)/xbj
          if(w2bj.lt.ckin(39)) mint(51)=1
          if(ckin(40).gt.0d0.and.w2bj.gt.ckin(40)) mint(51)=1
        endif
 
      elseif(ilim.eq.1) then
C...Calculate limits on tau
C...0) due to definition
        taumn0=0d0
        taumx0=1d0
C...1) due to limits on subsystem mass
        taumn1=ckin(1)**2/vint(2)
        taumx1=1d0
        if(ckin(2).ge.0d0) taumx1=ckin(2)**2/vint(2)
C...2) due to limits on pT-hat (and non-overlapping rapidity intervals)
        tm3=sqrt(sqm3+pthmin**2)
        tm4=sqrt(sqm4+pthmin**2)
        ydcosh=1d0
        if(ckin09.gt.ckin12) ydcosh=cosh(ckin09-ckin12)
        taumn2=(tm3**2+2d0*tm3*tm4*ydcosh+tm4**2)/vint(2)
        taumx2=1d0
C...3) due to limits on pT-hat and cos(theta-hat)
        cth2mn=min(ckin(27)**2,ckin(28)**2)
        cth2mx=max(ckin(27)**2,ckin(28)**2)
        taumn3=0d0
        if(ckin(27)*ckin(28).gt.0d0) taumn3=
     &  (sqrt(sqm3+pthmin**2/(1d0-cth2mn))+
     &  sqrt(sqm4+pthmin**2/(1d0-cth2mn)))**2/vint(2)
        taumx3=1d0
        if(ckin(4).ge.0d0.and.cth2mx.lt.1d0) taumx3=
     &  (sqrt(sqm3+ckin(4)**2/(1d0-cth2mx))+
     &  sqrt(sqm4+ckin(4)**2/(1d0-cth2mx)))**2/vint(2)
C...4) due to limits on x1 and x2
        taumn4=ckin(21)*ckin(23)
        taumx4=ckin(22)*ckin(24)
C...5) due to limits on xF
        taumn5=0d0
        taumx5=max(1d0-ckin(25),1d0+ckin(26))
C...6) due to limits on that and uhat
        taumn6=(sqm3+sqm4+ckin(35)+ckin(37))/vint(2)
        taumx6=1d0
        if(ckin(36).gt.0d0.and.ckin(38).gt.0d0) taumx6=
     &  (sqm3+sqm4+ckin(36)+ckin(38))/vint(2)
 
C...Net effect of all separate limits.
        vint(11)=max(taumn0,taumn1,taumn2,taumn3,taumn4,taumn5,taumn6)
        vint(31)=min(taumx0,taumx1,taumx2,taumx3,taumx4,taumx5,taumx6)
        if(mint(47).eq.1.and.(istsb.eq.1.or.istsb.eq.2)) then
          vint(11)=0.99999d0
          vint(31)=1.00001d0
        elseif(mint(47).eq.5) then
          vint(31)=min(vint(31),0.999998d0)
        endif
        if(vint(31).le.vint(11)) mint(51)=1
 
      elseif(ilim.eq.2) then
C...Calculate limits on y*
        taue=tau
        if(istsb.ge.3.and.istsb.le.5) taue=vint(26)
        taurt=sqrt(taue)
C...0) due to kinematics
        ystmn0=log(taurt)
        ystmx0=-ystmn0
C...1) due to explicit limits
        ystmn1=ckin(7)
        ystmx1=ckin(8)
C...2) due to limits on x1
        ystmn2=log(max(taue,ckin(21))/taurt)
        ystmx2=log(max(taue,ckin(22))/taurt)
C...3) due to limits on x2
        ystmn3=-log(max(taue,ckin(24))/taurt)
        ystmx3=-log(max(taue,ckin(23))/taurt)
C...4) due to limits on xF
        yepmn4=0.5d0*abs(ckin(25))/taurt
        ystmn4=sign(log(max(1d-20,sqrt(1d0+yepmn4**2)+yepmn4)),ckin(25))
        yepmx4=0.5d0*abs(ckin(26))/taurt
        ystmx4=sign(log(max(1d-20,sqrt(1d0+yepmx4**2)+yepmx4)),ckin(26))
C...5) due to simultaneous limits on y-large and y-small
        yepsmn=(rm3-rm4)*sinh(ckin09-ckin11)
        yepsmx=(rm3-rm4)*sinh(ckin10-ckin12)
        ydifmn=abs(log(max(1d-20,sqrt(1d0+yepsmn**2)-yepsmn)))
        ydifmx=abs(log(max(1d-20,sqrt(1d0+yepsmx**2)-yepsmx)))
        ystmn5=0.5d0*(ckin09+ckin11-ydifmn)
        ystmx5=0.5d0*(ckin10+ckin12+ydifmx)
C...6) due to simultaneous limits on cos(theta-hat) and y-large or
C...   y-small
        cthlim=sqrt(max(0d0,1d0-4d0*pthmin**2/(be34**2*taue*vint(2))))
        rzmn=be34*max(ckin(27),-cthlim)
        rzmx=be34*min(ckin(28),cthlim)
        yex3mx=(1d0+rm3-rm4+rzmx)/max(1d-10,1d0+rm3-rm4-rzmx)
        yex4mx=(1d0+rm4-rm3-rzmn)/max(1d-10,1d0+rm4-rm3+rzmn)
        yex3mn=max(1d-10,1d0+rm3-rm4+rzmn)/(1d0+rm3-rm4-rzmn)
        yex4mn=max(1d-10,1d0+rm4-rm3-rzmx)/(1d0+rm4-rm3+rzmx)
        ystmn6=ckin09-0.5d0*log(max(yex3mx,yex4mx))
        ystmx6=ckin12-0.5d0*log(min(yex3mn,yex4mn))
 
C...Net effect of all separate limits.
        vint(12)=max(ystmn0,ystmn1,ystmn2,ystmn3,ystmn4,ystmn5,ystmn6)
        vint(32)=min(ystmx0,ystmx1,ystmx2,ystmx3,ystmx4,ystmx5,ystmx6)
        if(mint(47).eq.1) then
          vint(12)=-0.00001d0
          vint(32)=0.00001d0
        elseif(mint(47).eq.2) then
          vint(12)=0.99999d0*ystmx0
          vint(32)=1.00001d0*ystmx0
        elseif(mint(47).eq.3) then
          vint(12)=-1.00001d0*ystmx0
          vint(32)=-0.99999d0*ystmx0
        elseif(mint(47).eq.5) then
          ystee=log(0.999999d0/taurt)
          vint(12)=max(vint(12),-ystee)
          vint(32)=min(vint(32),ystee)
        endif
        if(vint(32).le.vint(12)) mint(51)=1
 
      elseif(ilim.eq.3) then
C...Calculate limits on cos(theta-hat)
        yst=vint(22)
C...0) due to definition
        ctnmn0=-1d0
        ctnmx0=0d0
        ctpmn0=0d0
        ctpmx0=1d0
C...1) due to explicit limits
        ctnmn1=min(0d0,ckin(27))
        ctnmx1=min(0d0,ckin(28))
        ctpmn1=max(0d0,ckin(27))
        ctpmx1=max(0d0,ckin(28))
C...2) due to limits on pT-hat
        ctnmn2=-sqrt(max(0d0,1d0-4d0*pthmin**2/(be34**2*tau*vint(2))))
        ctpmx2=-ctnmn2
        ctnmx2=0d0
        ctpmn2=0d0
        if(ckin(4).ge.0d0) then
          ctnmx2=-sqrt(max(0d0,1d0-4d0*ckin(4)**2/
     &    (be34**2*tau*vint(2))))
          ctpmn2=-ctnmx2
        endif
C...3) due to limits on y-large and y-small
        ctnmn3=min(0d0,max((1d0+rm3-rm4)/be34*tanh(ckin11-yst),
     &  -(1d0-rm3+rm4)/be34*tanh(ckin10-yst)))
        ctnmx3=min(0d0,(1d0+rm3-rm4)/be34*tanh(ckin12-yst),
     &  -(1d0-rm3+rm4)/be34*tanh(ckin09-yst))
        ctpmn3=max(0d0,(1d0+rm3-rm4)/be34*tanh(ckin09-yst),
     &  -(1d0-rm3+rm4)/be34*tanh(ckin12-yst))
        ctpmx3=max(0d0,min((1d0+rm3-rm4)/be34*tanh(ckin10-yst),
     &  -(1d0-rm3+rm4)/be34*tanh(ckin11-yst)))
C...4) due to limits on that
        ctnmn4=-1d0
        ctnmx4=0d0
        ctpmn4=0d0
        ctpmx4=1d0
        sh=tau*vint(2)
        if(ckin(35).gt.0d0) then
          ctlim=(1d0-rm3-rm4-2d0*ckin(35)/sh)/be34
          if(ctlim.gt.0d0) then
            ctpmx4=ctlim
          else
            ctpmx4=0d0
            ctnmx4=ctlim
          endif
        endif
        if(ckin(36).gt.0d0) then
          ctlim=(1d0-rm3-rm4-2d0*ckin(36)/sh)/be34
          if(ctlim.lt.0d0) then
            ctnmn4=ctlim
          else
            ctnmn4=0d0
            ctpmn4=ctlim
          endif
        endif
C...5) due to limits on uhat
        ctnmn5=-1d0
        ctnmx5=0d0
        ctpmn5=0d0
        ctpmx5=1d0
        if(ckin(37).gt.0d0) then
          ctlim=(2d0*ckin(37)/sh-(1d0-rm3-rm4))/be34
          if(ctlim.lt.0d0) then
            ctnmn5=ctlim
          else
            ctnmn5=0d0
            ctpmn5=ctlim
          endif
        endif
        if(ckin(38).gt.0d0) then
          ctlim=(2d0*ckin(38)/sh-(1d0-rm3-rm4))/be34
          if(ctlim.gt.0d0) then
            ctpmx5=ctlim
          else
            ctpmx5=0d0
            ctnmx5=ctlim
          endif
        endif
 
C...Net effect of all separate limits.
        vint(13)=max(ctnmn0,ctnmn1,ctnmn2,ctnmn3,ctnmn4,ctnmn5)
        vint(33)=min(ctnmx0,ctnmx1,ctnmx2,ctnmx3,ctnmx4,ctnmx5)
        vint(14)=max(ctpmn0,ctpmn1,ctpmn2,ctpmn3,ctpmn4,ctpmn5)
        vint(34)=min(ctpmx0,ctpmx1,ctpmx2,ctpmx3,ctpmx4,ctpmx5)
        if(vint(33).le.vint(13).and.vint(34).le.vint(14)) mint(51)=1
 
      elseif(ilim.eq.4) then
C...Calculate limits on tau'
C...0) due to kinematics
        tapmn0=tau
        if(istsb.eq.5.and.kfpr(isub,2).gt.0) then
          pqrat=2d0*pmas(jamcomp(kfpr(isub,2)),1)/vint(1)
          tapmn0=(sqrt(tau)+pqrat)**2
        endif
        tapmx0=1d0
C...1) due to explicit limits
        tapmn1=ckin(31)**2/vint(2)
        tapmx1=1d0
        if(ckin(32).ge.0d0) tapmx1=ckin(32)**2/vint(2)
 
C...Net effect of all separate limits.
        vint(16)=max(tapmn0,tapmn1)
        vint(36)=min(tapmx0,tapmx1)
        if(mint(47).eq.1) then
          vint(16)=0.99999d0
          vint(36)=1.00001d0
        endif
        if(vint(36).le.vint(16)) mint(51)=1
 
      endif
      return
 
C...Special case for low-pT and multiple interactions:
C...effective kinematical limits for tau, y*, cos(theta-hat).
  100 if(ilim.eq.0) then
      elseif(ilim.eq.1) then
        if(mstp(82).le.1) vint(11)=4d0*parp(81)**2/vint(2)
        if(mstp(82).ge.2) vint(11)=parp(82)**2/vint(2)
        vint(31)=1d0
      elseif(ilim.eq.2) then
        vint(12)=0.5d0*log(vint(21))
        vint(32)=-vint(12)
      elseif(ilim.eq.3) then
        if(mstp(82).le.1) st2eff=4d0*parp(81)**2/(vint(21)*vint(2))
        if(mstp(82).ge.2) st2eff=0.01d0*parp(82)**2/(vint(21)*vint(2))
        vint(13)=-sqrt(max(0d0,1d0-st2eff))
        vint(33)=0d0
        vint(14)=0d0
        vint(34)=-vint(13)
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYKMAP
C...Maps a uniform distribution into a distribution of a kinematical
C...variable according to one of the possibilities allowed. It is
C...assumed that kinematical limits have been set by a PYKLIM call.
 
      subroutine pjkmap(ivar,mvar,vvar)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/pjsubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      save /jydat1/,/jydat2/,/pjsubs/,/pjpars/,/pjint1/,/pjint2/
 
C...Convert VVAR to tau variable.
      isub=mint(1)
      istsb=iset(isub)
      if(ivar.eq.1) then
        taumin=vint(11)
        taumax=vint(31)
        if(mvar.eq.3.or.mvar.eq.4) then
          taure=vint(73)
          gamre=vint(74)
        elseif(mvar.eq.5.or.mvar.eq.6) then
          taure=vint(75)
          gamre=vint(76)
        endif
        if(mint(47).eq.1.and.(istsb.eq.1.or.istsb.eq.2)) then
          tau=1d0
        elseif(mvar.eq.1) then
          tau=taumin*(taumax/taumin)**vvar
        elseif(mvar.eq.2) then
          tau=taumax*taumin/(taumin+(taumax-taumin)*vvar)
        elseif(mvar.eq.3.or.mvar.eq.5) then
          ratgen=(taure+taumax)/(taure+taumin)*taumin/taumax
          tau=taure*taumin/((taure+taumin)*ratgen**vvar-taumin)
        elseif(mvar.eq.4.or.mvar.eq.6) then
          aupp=atan((taumax-taure)/gamre)
          alow=atan((taumin-taure)/gamre)
          tau=taure+gamre*tan(alow+(aupp-alow)*vvar)
        else
          aupp=log(max(2d-6,1d0-taumax))
          alow=log(max(2d-6,1d0-taumin))
          tau=1d0-exp(aupp+vvar*(alow-aupp))
        endif
        vint(21)=min(taumax,max(taumin,tau))
 
C...Convert VVAR to y* variable.
      elseif(ivar.eq.2) then
        ystmin=vint(12)
        ystmax=vint(32)
        taue=vint(21)
        if(istsb.ge.3.and.istsb.le.5) taue=vint(26)
        if(mint(47).eq.1) then
          yst=0d0
        elseif(mint(47).eq.2) then
          yst=-0.5d0*log(taue)
        elseif(mint(47).eq.3) then
          yst=0.5d0*log(taue)
        elseif(mvar.eq.1) then
          yst=ystmin+(ystmax-ystmin)*sqrt(vvar)
        elseif(mvar.eq.2) then
          yst=ystmax-(ystmax-ystmin)*sqrt(1d0-vvar)
        elseif(mvar.eq.3) then
          aupp=atan(exp(ystmax))
          alow=atan(exp(ystmin))
          yst=log(tan(alow+(aupp-alow)*vvar))
        elseif(mvar.eq.4) then
          yst0=-0.5d0*log(taue)
          aupp=log(max(1d-6,exp(yst0-ystmin)-1d0))
          alow=log(max(1d-6,exp(yst0-ystmax)-1d0))
          yst=yst0-log(1d0+exp(alow+vvar*(aupp-alow)))
        else
          yst0=-0.5d0*log(taue)
          aupp=log(max(1d-6,exp(yst0+ystmin)-1d0))
          alow=log(max(1d-6,exp(yst0+ystmax)-1d0))
          yst=log(1d0+exp(aupp+vvar*(alow-aupp)))-yst0
        endif
        vint(22)=min(ystmax,max(ystmin,yst))
 
C...Convert VVAR to cos(theta-hat) variable.
      elseif(ivar.eq.3) then
        rm34=max(1d-20,2d0*vint(63)*vint(64)/(vint(21)*vint(2))**2)
        rsqm=1d0+rm34
        if(2d0*vint(71)**2/(vint(21)*vint(2)).lt.0.0001d0)
     &  rm34=max(rm34,2d0*vint(71)**2/(vint(21)*vint(2)))
        ctnmin=vint(13)
        ctnmax=vint(33)
        ctpmin=vint(14)
        ctpmax=vint(34)
        if(mvar.eq.1) then
          aneg=ctnmax-ctnmin
          apos=ctpmax-ctpmin
          if(aneg.gt.0d0.and.vvar*(aneg+apos).le.aneg) then
            vctn=vvar*(aneg+apos)/aneg
            cth=ctnmin+(ctnmax-ctnmin)*vctn
          else
            vctp=(vvar*(aneg+apos)-aneg)/apos
            cth=ctpmin+(ctpmax-ctpmin)*vctp
          endif
        elseif(mvar.eq.2) then
          rmnmin=max(rm34,rsqm-ctnmin)
          rmnmax=max(rm34,rsqm-ctnmax)
          rmpmin=max(rm34,rsqm-ctpmin)
          rmpmax=max(rm34,rsqm-ctpmax)
          aneg=log(rmnmin/rmnmax)
          apos=log(rmpmin/rmpmax)
          if(aneg.gt.0d0.and.vvar*(aneg+apos).le.aneg) then
            vctn=vvar*(aneg+apos)/aneg
            cth=rsqm-rmnmin*(rmnmax/rmnmin)**vctn
          else
            vctp=(vvar*(aneg+apos)-aneg)/apos
            cth=rsqm-rmpmin*(rmpmax/rmpmin)**vctp
          endif
        elseif(mvar.eq.3) then
          rmnmin=max(rm34,rsqm+ctnmin)
          rmnmax=max(rm34,rsqm+ctnmax)
          rmpmin=max(rm34,rsqm+ctpmin)
          rmpmax=max(rm34,rsqm+ctpmax)
          aneg=log(rmnmax/rmnmin)
          apos=log(rmpmax/rmpmin)
          if(aneg.gt.0d0.and.vvar*(aneg+apos).le.aneg) then
            vctn=vvar*(aneg+apos)/aneg
            cth=rmnmin*(rmnmax/rmnmin)**vctn-rsqm
          else
            vctp=(vvar*(aneg+apos)-aneg)/apos
            cth=rmpmin*(rmpmax/rmpmin)**vctp-rsqm
          endif
        elseif(mvar.eq.4) then
          rmnmin=max(rm34,rsqm-ctnmin)
          rmnmax=max(rm34,rsqm-ctnmax)
          rmpmin=max(rm34,rsqm-ctpmin)
          rmpmax=max(rm34,rsqm-ctpmax)
          aneg=1d0/rmnmax-1d0/rmnmin
          apos=1d0/rmpmax-1d0/rmpmin
          if(aneg.gt.0d0.and.vvar*(aneg+apos).le.aneg) then
            vctn=vvar*(aneg+apos)/aneg
            cth=rsqm-1d0/(1d0/rmnmin+aneg*vctn)
          else
            vctp=(vvar*(aneg+apos)-aneg)/apos
            cth=rsqm-1d0/(1d0/rmpmin+apos*vctp)
          endif
        elseif(mvar.eq.5) then
          rmnmin=max(rm34,rsqm+ctnmin)
          rmnmax=max(rm34,rsqm+ctnmax)
          rmpmin=max(rm34,rsqm+ctpmin)
          rmpmax=max(rm34,rsqm+ctpmax)
          aneg=1d0/rmnmin-1d0/rmnmax
          apos=1d0/rmpmin-1d0/rmpmax
          if(aneg.gt.0d0.and.vvar*(aneg+apos).le.aneg) then
            vctn=vvar*(aneg+apos)/aneg
            cth=1d0/(1d0/rmnmin-aneg*vctn)-rsqm
          else
            vctp=(vvar*(aneg+apos)-aneg)/apos
            cth=1d0/(1d0/rmpmin-apos*vctp)-rsqm
          endif
        endif
        if(cth.lt.0d0) cth=min(ctnmax,max(ctnmin,cth))
        if(cth.gt.0d0) cth=min(ctpmax,max(ctpmin,cth))
        vint(23)=cth
 
C...Convert VVAR to tau' variable.
      elseif(ivar.eq.4) then
        tau=vint(21)
        taupmn=vint(16)
        taupmx=vint(36)
        if(mint(47).eq.1) then
          taup=1d0
        elseif(mvar.eq.1) then
          taup=taupmn*(taupmx/taupmn)**vvar
        elseif(mvar.eq.2) then
          aupp=(1d0-tau/taupmx)**4
          alow=(1d0-tau/taupmn)**4
          taup=tau/max(1d-7,1d0-(alow+(aupp-alow)*vvar)**0.25d0)
        else
          aupp=log(max(2d-6,1d0-taupmx))
          alow=log(max(2d-6,1d0-taupmn))
          taup=1d0-exp(aupp+vvar*(alow-aupp))
        endif
        vint(26)=min(taupmx,max(taupmn,taup))
 
C...Selection of extra variables needed in 2 -> 3 process:
C...pT1, pT2, phi1, phi2, y3 for three outgoing particles.
C...Since no options are available, the functions of PYKLIM
C...and PYKMAP are joint for these choices.
      elseif(ivar.eq.5) then
 
C...Read out total energy and particle masses.
        mint(51)=0
        mptpk=1
        if(isub.eq.123.or.isub.eq.124.or.isub.eq.173.or.isub.eq.174
     &  .or.isub.eq.178.or.isub.eq.179) mptpk=2
        shp=vint(26)*vint(2)
        shpr=sqrt(shp)
        pm1=vint(201)
        pm2=vint(206)
        pm3=sqrt(vint(21))*vint(1)
        if(pm1+pm2+pm3.gt.0.9999d0*shpr) then
          mint(51)=1
          return
        endif
        pmrs1=vint(204)**2
        pmrs2=vint(209)**2
 
C...Specify coefficients of pT choice; upper and lower limits.
        if(mptpk.eq.1) then
          hwt1=0.4d0
          hwt2=0.4d0
        else
          hwt1=0.05d0
          hwt2=0.05d0
        endif
        hwt3=1d0-hwt1-hwt2
        ptsmx1=((shp-pm1**2-(pm2+pm3)**2)**2-(2d0*pm1*(pm2+pm3))**2)/
     &  (4d0*shp)
        if(ckin(52).gt.0d0) ptsmx1=min(ptsmx1,ckin(52)**2)
        ptsmn1=ckin(51)**2
        ptsmx2=((shp-pm2**2-(pm1+pm3)**2)**2-(2d0*pm2*(pm1+pm3))**2)/
     &  (4d0*shp)
        if(ckin(54).gt.0d0) ptsmx2=min(ptsmx2,ckin(54)**2)
        ptsmn2=ckin(53)**2
 
C...Select transverse momenta according to
C...dp_T^2 * (a + b/(M^2 + p_T^2) + c/(M^2 + p_T^2)^2).
        hmx=pmrs1+ptsmx1
        hmn=pmrs1+ptsmn1
        if(hmx.lt.1.0001d0*hmn) then
          mint(51)=1
          return
        endif
        hde=ptsmx1-ptsmn1
        rpt=pjr(0)
        if(rpt.lt.hwt1) then
          pts1=ptsmn1+pjr(0)*hde
        elseif(rpt.lt.hwt1+hwt2) then
          pts1=max(ptsmn1,hmn*(hmx/hmn)**pjr(0)-pmrs1)
        else
          pts1=max(ptsmn1,hmn*hmx/(hmn+pjr(0)*hde)-pmrs1)
        endif
        wtpts1=hde/(hwt1+hwt2*hde/(log(hmx/hmn)*(pmrs1+pts1))+
     &  hwt3*hmn*hmx/(pmrs1+pts1)**2)
        hmx=pmrs2+ptsmx2
        hmn=pmrs2+ptsmn2
        if(hmx.lt.1.0001d0*hmn) then
          mint(51)=1
          return
        endif
        hde=ptsmx2-ptsmn2
        rpt=pjr(0)
        if(rpt.lt.hwt1) then
          pts2=ptsmn2+pjr(0)*hde
        elseif(rpt.lt.hwt1+hwt2) then
          pts2=max(ptsmn2,hmn*(hmx/hmn)**pjr(0)-pmrs2)
        else
          pts2=max(ptsmn2,hmn*hmx/(hmn+pjr(0)*hde)-pmrs2)
        endif
        wtpts2=hde/(hwt1+hwt2*hde/(log(hmx/hmn)*(pmrs2+pts2))+
     &  hwt3*hmn*hmx/(pmrs2+pts2)**2)
 
C...Select azimuthal angles and check pT choice.
        phi1=paru(2)*pjr(0)
        phi2=paru(2)*pjr(0)
        phir=phi2-phi1
        pts3=max(0d0,pts1+pts2+2d0*sqrt(pts1*pts2)*cos(phir))
        if(pts3.lt.ckin(55)**2.or.(ckin(56).gt.0d0.and.pts3.gt.
     &  ckin(56)**2)) then
          mint(51)=1
          return
        endif
 
C...Calculate transverse masses and check phase space not closed.
        pms1=pm1**2+pts1
        pms2=pm2**2+pts2
        pms3=pm3**2+pts3
        pmt1=sqrt(pms1)
        pmt2=sqrt(pms2)
        pmt3=sqrt(pms3)
        pm12=(pmt1+pmt2)**2
        if(pmt1+pmt2+pmt3.gt.0.9999d0*shpr) then
          mint(51)=1
          return
        endif
 
C...Select rapidity for particle 3 and check phase space not closed.
        y3max=log((shp+pms3-pm12+sqrt(max(0d0,(shp-pms3-pm12)**2-
     &  4d0*pms3*pm12)))/(2d0*shpr*pmt3))
        if(y3max.lt.1d-6) then
          mint(51)=1
          return
        endif
        y3=(2d0*pjr(0)-1d0)*0.999999d0*y3max
        pz3=pmt3*sinh(y3)
        pe3=pmt3*cosh(y3)
 
C...Find momentum transfers in two mirror solutions (in 1-2 frame).
        pz12=-pz3
        pe12=shpr-pe3
        pms12=pe12**2-pz12**2
        sql12=sqrt(max(0d0,(pms12-pms1-pms2)**2-4d0*pms1*pms2))
        if(sql12.lt.1d-6*shp) then
          mint(51)=1
          return
        endif
        pmm1=pms12+pms1-pms2
        pmm2=pms12+pms2-pms1
        tfac=-shpr/(2d0*pms12)
        t1p=tfac*(pe12-pz12)*(pmm1-sql12)
        t1n=tfac*(pe12-pz12)*(pmm1+sql12)
        t2p=tfac*(pe12+pz12)*(pmm2-sql12)
        t2n=tfac*(pe12+pz12)*(pmm2+sql12)
 
C...Construct relative mirror weights and make choice.
        if(mptpk.eq.1) then
          wtpu=1d0
          wtnu=1d0
        else
          wtpu=1d0/((t1p-pmrs1)*(t2p-pmrs2))**2
          wtnu=1d0/((t1n-pmrs1)*(t2n-pmrs2))**2
        endif
        wtp=wtpu/(wtpu+wtnu)
        wtn=wtnu/(wtpu+wtnu)
        eps=1d0
        if(wtn.gt.pjr(0)) eps=-1d0
 
C...Store result of variable choice and associated weights.
        vint(202)=pts1
        vint(207)=pts2
        vint(203)=phi1
        vint(208)=phi2
        vint(205)=wtpts1
        vint(210)=wtpts2
        vint(211)=y3
        vint(212)=y3max
        vint(213)=eps
        if(eps.gt.0d0) then
          vint(214)=1d0/wtp
          vint(215)=t1p
          vint(216)=t2p
        else
          vint(214)=1d0/wtn
          vint(215)=t1n
          vint(216)=t2n
        endif
        vint(217)=-0.5d0*tfac*(pe12-pz12)*(pmm2+eps*sql12)
        vint(218)=-0.5d0*tfac*(pe12+pz12)*(pmm1+eps*sql12)
        vint(219)=0.5d0*(pms12-pts3)
        vint(220)=sql12
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYHFTH
C...Gives threshold attractive/repulsive factor for heavy flavour
C...production.
 
      function pjhfth(sh,sqm,fratt)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      save /jydat1/,/pjpars/,/pjint1/
 
C...Value for alpha_strong.
      if(mstp(35).le.1) then
        alssg=parp(35)
      else
        mst115=mstu(115)
        mstu(115)=mstp(36)
        q2bn=sqrt(max(1d0,sqm*((sqrt(sh)-2d0*sqrt(sqm))**2+
     &  parp(36)**2)))
        alssg=pjalps(q2bn)
        mstu(115)=mst115
      endif
 
C...Evaluate attractive and repulsive factors.
      xattr=4d0*paru(1)*alssg/(3d0*sqrt(max(1d-20,1d0-4d0*sqm/sh)))
      fattr=xattr/(1d0-exp(-min(50d0,xattr)))
      xrepu=paru(1)*alssg/(6d0*sqrt(max(1d-20,1d0-4d0*sqm/sh)))
      frepu=xrepu/(exp(min(50d0,xrepu))-1d0)
      pjhfth=fratt*fattr+(1d0-fratt)*frepu
      vint(138)=pjhfth
 
      return
      end
 
C*********************************************************************
 
C...PYSPLI:modifyed by JAM.
C...Splits a hadron remnant into two (partons or hadron + parton)
C...in case it is more complicated than just a quark or a diquark.
 
      subroutine pjspli(kf0,kflin,kflch,kflsp)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...JAM Commonblocks.
      include 'jam2.inc'
C...Commonblocks.
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      save /pjpars/,/pjint1/
C...Local array.
      dimension kfl(3)
 
cjam++
      kf=kf0
      kfs=isign(1,kf)
      kc=jamcomp(kf)
      id0=kchg(kc,5)
      iz0=kchg(kc,1)/3
      if(id0.eq.id_nucls) then
        if(iz0.eq. 0) kf=2112*kfs
        if(iz0.eq. 1) kf=2212*kfs
      else if(id0.eq.id_delts) then
        if(iz0.eq.-1) kf=1114*kfs
        if(iz0.eq. 0) kf=2114*kfs
        if(iz0.eq. 1) kf=2214*kfs
        if(iz0.eq. 2) kf=2224*kfs
      else if(id0.eq.id_lambs) then
        kf=3122*kfs
      else if(id0.eq.id_sigms) then
        if(iz0.eq.-1) kf=3112*kfs
        if(iz0.eq. 0) kf=3212*kfs
        if(iz0.eq. 1) kf=3222*kfs
      else if(id0.eq.id_xis) then
        if(iz0.eq.-1) kf=3312*kfs
        if(iz0.eq. 0) kf=3322*kfs
      endif
cjam--

C...Preliminaries. Parton composition.
      kfa=iabs(kf)
      kfl(1)=mod(kfa/1000,10)
      kfl(2)=mod(kfa/100,10)
      kfl(3)=mod(kfa/10,10)

      if(kfa.eq.22.and.mint(109).eq.2) then
        kfl(2)=int(1.5d0+pjr(0))
        if(mint(105).eq.333) kfl(2)=3
        if(mint(105).eq.443) kfl(2)=4
        kfl(3)=kfl(2)
cjam++
c     elseif((kfa.eq.111.or.kfa.eq.113).and.pjr(0).gt.0.5d0) then
      else if(kfl(1).eq.0) then !mesons.
        kfma=kfl(2)*10+kfl(3)
        if(kfma.eq.11.and.pjr(0).gt.0.5d0)then
          kfl(2)=2
          kfl(3)=2
        else if(kfma.eq.22.and.pjr(0).gt.0.5)then
          kfl(2)=1
          kfl(3)=1
        endif
cjam--
      endif
      if(kflin.ne.21.and.kflin.ne.22.and.kflin.ne.23) then
        kflr=kflin*kfs
      else
        kflr=kflin
      endif
      kflch=0
 
C...Subdivide lepton.
      if(kfa.ge.11.and.kfa.le.18) then
        if(kflr.eq.kfa) then
          kflsp=kfs*22
        elseif(kflr.eq.22) then
          kflsp=kfa
        elseif(kflr.eq.-24.and.mod(kfa,2).eq.1) then
          kflsp=kfa+1
        elseif(kflr.eq.24.and.mod(kfa,2).eq.0) then
          kflsp=kfa-1
        elseif(kflr.eq.21) then
          kflsp=kfa
          kflch=kfs*21
        else
          kflsp=kfa
          kflch=-kflr
        endif
 
C...Subdivide photon.
      elseif(kfa.eq.22.and.mint(109).ne.2) then
        if(kflr.ne.21) then
          kflsp=-kflr
        else
          ragr=0.75d0*pjr(0)
          kflsp=1
          if(ragr.gt.0.125d0) kflsp=2
          if(ragr.gt.0.625d0) kflsp=3
          if(pjr(0).gt.0.5d0) kflsp=-kflsp
          kflch=-kflsp
        endif
 
C...Subdivide Reggeon or Pomeron.
      elseif(kfa.eq.28.or.kfa.eq.29) then
        if(kflin.eq.21) then
          kflsp=kfs*21
        else
          kflsp=-kflin
        endif
 
C...Subdivide meson.
      elseif(kfl(1).eq.0) then
        kfl(2)=kfl(2)*(-1)**kfl(2)
        kfl(3)=-kfl(3)*(-1)**iabs(kfl(2))
        if(kflr.eq.kfl(2)) then
          kflsp=kfl(3)
        elseif(kflr.eq.kfl(3)) then
          kflsp=kfl(2)
        elseif(kflr.eq.21.and.pjr(0).gt.0.5d0) then
          kflsp=kfl(2)
          kflch=kfl(3)
        elseif(kflr.eq.21) then
          kflsp=kfl(3)
          kflch=kfl(2)
        elseif(kflr*kfl(2).gt.0) then
          call pjkfdi(-kflr,kfl(2),kfdump,kflch)
          kflsp=kfl(3)
        else
          call pjkfdi(-kflr,kfl(3),kfdump,kflch)
          kflsp=kfl(2)
        endif
 
C...Subdivide baryon.
      else
        nagr=0
        do 100 j=1,3
          if(kflr.eq.kfl(j)) nagr=nagr+1
  100   continue
        if(nagr.ge.1) then
          ragr=0.00001d0+(nagr-0.00002d0)*pjr(0)
          iagr=0
          do 110 j=1,3
            if(kflr.eq.kfl(j)) ragr=ragr-1d0
            if(iagr.eq.0.and.ragr.le.0d0) iagr=j
  110     continue
        else
          iagr=1.00001d0+2.99998d0*pjr(0)
        endif
        id1=1
        if(iagr.eq.1) id1=2
        if(iagr.eq.1.and.kfl(3).gt.kfl(2)) id1=3
        id2=6-iagr-id1
        ksp=3
        if(mod(kfa,10).eq.2.and.kfl(1).eq.kfl(2)) then
          if(iagr.ne.3.and.pjr(0).gt.0.25d0) ksp=1
        elseif(mod(kfa,10).eq.2.and.kfl(2).ge.kfl(3)) then
          if(iagr.ne.1.and.pjr(0).gt.0.25d0) ksp=1
        elseif(mod(kfa,10).eq.2) then
          if(iagr.eq.1) ksp=1
          if(iagr.ne.1.and.pjr(0).gt.0.75d0) ksp=1
        endif
        kflsp=1000*kfl(id1)+100*kfl(id2)+ksp
        if(kflr.eq.21) then
          kflch=kfl(iagr)
        elseif(nagr.eq.0.and.kflr.gt.0) then
          call pjkfdi(-kflr,kfl(iagr),kfdump,kflch)
        elseif(nagr.eq.0) then
          call pjkfdi(10000+kflsp,-kflr,kfdump,kflch)
          kflsp=kfl(iagr)
        endif
      endif
 
cjam++ check
      if(jamcomp(kflsp).eq.0)then
        write(check(1),8000)kf,kflin,kflch,kflsp
        write(check(2),8100)kflr,nagr
        call jamerrm(30,2,'(pjspli:)kflsp=0')
      endif
cjam--

C...Add on correct sign for result.
      kflch=kflch*kfs
      kflsp=kflsp*kfs

      return
 8000 format('kf=',i9,' kflin=',i9,' kflch=',i9,' kfslp=',i9)
 8100 format('kflr=',i9,' nagr=',i9)
 
      end
 
C*********************************************************************
 
C...PYGAMM
C...Gives ordinary Gamma function Gamma(x) for positive, real arguments;
C...see M. Abramowitz, I. A. Stegun: Handbook of Mathematical Functions
C...(Dover, 1965) 6.1.36.
 
      function pjgamm(x)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Local array and data.
      dimension b(8)
      data b/-0.577191652d0,0.988205891d0,-0.897056937d0,0.918206857d0,
     &-0.756704078d0,0.482199394d0,-0.193527818d0,0.035868343d0/
 
      nx=int(x)
      dx=x-nx
 
      pjgamm=1d0
      dxp=1d0
      do 100 i=1,8
        dxp=dxp*dx
        pjgamm=pjgamm+b(i)*dxp
  100 continue
      if(x.lt.1d0) then
        pjgamm=pjgamm/x
      else
        do 110 ix=1,nx-1
          pjgamm=(x-ix)*pjgamm
  110   continue
      endif
 
      return
      end
 
C***********************************************************************
 
C...PYWAUX
C...Calculates real and imaginary parts of the auxiliary functions W1
C...and W2; see R. K. Ellis, I. Hinchliffe, M. Soldate and J. J. van
C...der Bij, Nucl. Phys. B297 (1988) 221.
 
      subroutine pjwaux(iaux,eps,wre,wim)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      save /jydat1/
 
      asinh(x)=log(x+sqrt(x**2+1d0))
      acosh(x)=log(x+sqrt(x**2-1d0))
 
      if(eps.lt.0d0) then
        if(iaux.eq.1) wre=2d0*sqrt(1d0-eps)*asinh(sqrt(-1d0/eps))
        if(iaux.eq.2) wre=4d0*(asinh(sqrt(-1d0/eps)))**2
        wim=0d0
      elseif(eps.lt.1d0) then
        if(iaux.eq.1) wre=2d0*sqrt(1d0-eps)*acosh(sqrt(1d0/eps))
        if(iaux.eq.2) wre=4d0*(acosh(sqrt(1d0/eps)))**2-paru(1)**2
        if(iaux.eq.1) wim=-paru(1)*sqrt(1d0-eps)
        if(iaux.eq.2) wim=-4d0*paru(1)*acosh(sqrt(1d0/eps))
      else
        if(iaux.eq.1) wre=2d0*sqrt(eps-1d0)*asin(sqrt(1d0/eps))
        if(iaux.eq.2) wre=-4d0*(asin(sqrt(1d0/eps)))**2
        wim=0d0
      endif
 
      return
      end
 
C***********************************************************************
 
C...PYI3AU
C...Calculates real and imaginary parts of the auxiliary function I3;
C...see R. K. Ellis, I. Hinchliffe, M. Soldate and J. J. van der Bij,
C...Nucl. Phys. B297 (1988) 221.
 
      subroutine pji3au(eps,rat,y3re,y3im)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      save /jydat1/
 
      be=0.5d0*(1d0+sqrt(1d0+rat*eps))
      if(eps.lt.1d0) ga=0.5d0*(1d0+sqrt(1d0-eps))
 
      if(eps.lt.0d0) then
        if(abs(eps).lt.1.d-4.and.abs(rat*eps).lt.1.d-4) then
          f3re=pjspen(-0.25d0*eps/(1d0+0.25d0*(rat-1d0)*eps),0d0,1)-
     &    pjspen((1d0-0.25d0*eps)/(1d0+0.25d0*(rat-1d0)*eps),0d0,1)+
     &    pjspen(0.25d0*(rat+1d0)*eps/(1d0+0.25d0*rat*eps),0d0,1)-
     &    pjspen((rat+1d0)/rat,0d0,1)+0.5d0*(log(1d0+0.25d0*rat*eps)**2-
     &    log(0.25d0*rat*eps)**2)+log(1d0-0.25d0*eps)*
     &    log((1d0+0.25d0*(rat-1d0)*eps)/(1d0+0.25d0*rat*eps))+
     &    log(-0.25d0*eps)*log(0.25d0*rat*eps/(1d0+0.25d0*(rat-1d0)*
     &    eps))
        elseif(abs(eps).lt.1.d-4.and.abs(rat*eps).ge.1.d-4) then
          f3re=pjspen(-0.25d0*eps/(be-0.25d0*eps),0d0,1)-
     &    pjspen((1d0-0.25d0*eps)/(be-0.25d0*eps),0d0,1)+
     &    pjspen((be-1d0+0.25d0*eps)/be,0d0,1)-
     &    pjspen((be-1d0+0.25d0*eps)/(be-1d0),0d0,1)+
     &    0.5d0*(log(be)**2-log(be-1d0)**2)+
     &    log(1d0-0.25d0*eps)*log((be-0.25d0*eps)/be)+
     &    log(-0.25d0*eps)*log((be-1d0)/(be-0.25d0*eps))
        elseif(abs(eps).ge.1.d-4.and.abs(rat*eps).lt.1.d-4) then
          f3re=pjspen((ga-1d0)/(ga+0.25d0*rat*eps),0d0,1)-
     &    pjspen(ga/(ga+0.25d0*rat*eps),0d0,1)+
     &    pjspen((1d0+0.25d0*rat*eps-ga)/(1d0+0.25d0*rat*eps),0d0,1)-
     &    pjspen((1d0+0.25d0*rat*eps-ga)/(0.25d0*rat*eps),0d0,1)+
     &    0.5d0*(log(1d0+0.25d0*rat*eps)**2-log(0.25d0*rat*eps)**2)+
     &    log(ga)*log((ga+0.25d0*rat*eps)/(1d0+0.25d0*rat*eps))+
     &    log(ga-1d0)*log(0.25d0*rat*eps/(ga+0.25d0*rat*eps))
        else
          f3re=pjspen((ga-1d0)/(ga+be-1d0),0d0,1)-
     &    pjspen(ga/(ga+be-1d0),0d0,1)+pjspen((be-ga)/be,0d0,1)-
     &    pjspen((be-ga)/(be-1d0),0d0,1)+0.5d0*(log(be)**2-
     &    log(be-1d0)**2)+log(ga)*log((ga+be-1d0)/be)+
     &    log(ga-1d0)*log((be-1d0)/(ga+be-1d0))
        endif
        f3im=0d0
      elseif(eps.lt.1d0) then
        if(abs(eps).lt.1.d-4.and.abs(rat*eps).lt.1.d-4) then
          f3re=pjspen(-0.25d0*eps/(1d0+0.25d0*(rat-1d0)*eps),0d0,1)-
     &    pjspen((1d0-0.25d0*eps)/(1d0+0.25d0*(rat-1d0)*eps),0d0,1)+
     &    pjspen((1d0-0.25d0*eps)/(-0.25d0*(rat+1d0)*eps),0d0,1)-
     &    pjspen(1d0/(rat+1d0),0d0,1)+log((1d0-0.25d0*eps)/
     &    (0.25d0*eps))*log((1d0+0.25d0*(rat-1d0)*eps)/
     &    (0.25d0*(rat+1d0)*eps))
          f3im=-paru(1)*log((1d0+0.25d0*(rat-1d0)*eps)/
     &    (0.25d0*(rat+1d0)*eps))
        elseif(abs(eps).lt.1.d-4.and.abs(rat*eps).ge.1.d-4) then
          f3re=pjspen(-0.25d0*eps/(be-0.25d0*eps),0d0,1)-
     &    pjspen((1d0-0.25d0*eps)/(be-0.25d0*eps),0d0,1)+
     &    pjspen((1d0-0.25d0*eps)/(1d0-0.25d0*eps-be),0d0,1)-
     &    pjspen(-0.25d0*eps/(1d0-0.25d0*eps-be),0d0,1)+
     &    log((1d0-0.25d0*eps)/(0.25d0*eps))*
     &    log((be-0.25d0*eps)/(be-1d0+0.25d0*eps))
          f3im=-paru(1)*log((be-0.25d0*eps)/(be-1d0+0.25d0*eps))
        elseif(abs(eps).ge.1.d-4.and.abs(rat*eps).lt.1.d-4) then
          f3re=pjspen((ga-1d0)/(ga+0.25d0*rat*eps),0d0,1)-
     &    pjspen(ga/(ga+0.25d0*rat*eps),0d0,1)+
     &    pjspen(ga/(ga-1d0-0.25d0*rat*eps),0d0,1)-
     &    pjspen((ga-1d0)/(ga-1d0-0.25d0*rat*eps),0d0,1)+
     &    log(ga/(1d0-ga))*log((ga+0.25d0*rat*eps)/
     &    (1d0+0.25d0*rat*eps-ga))
          f3im=-paru(1)*log((ga+0.25d0*rat*eps)/
     &    (1d0+0.25d0*rat*eps-ga))
        else
          f3re=pjspen((ga-1d0)/(ga+be-1d0),0d0,1)-
     &    pjspen(ga/(ga+be-1d0),0d0,1)+pjspen(ga/(ga-be),0d0,1)-
     &    pjspen((ga-1d0)/(ga-be),0d0,1)+log(ga/(1d0-ga))*
     &    log((ga+be-1d0)/(be-ga))
          f3im=-paru(1)*log((ga+be-1d0)/(be-ga))
        endif
      else
        rsq=eps/(eps-1d0+(2d0*be-1d0)**2)
        rcthe=rsq*(1d0-2d0*be/eps)
        rsthe=sqrt(max(0d0,rsq-rcthe**2))
        rcphi=rsq*(1d0+2d0*(be-1d0)/eps)
        rsphi=sqrt(max(0d0,rsq-rcphi**2))
        r=sqrt(rsq)
        the=acos(max(-0.999999d0,min(0.999999d0,rcthe/r)))
        phi=acos(max(-0.999999d0,min(0.999999d0,rcphi/r)))
        f3re=pjspen(rcthe,rsthe,1)+pjspen(rcthe,-rsthe,1)-
     &  pjspen(rcphi,rsphi,1)-pjspen(rcphi,-rsphi,1)+
     &  (phi-the)*(phi+the-paru(1))
        f3im=pjspen(rcthe,rsthe,2)+pjspen(rcthe,-rsthe,2)-
     &  pjspen(rcphi,rsphi,2)-pjspen(rcphi,-rsphi,2)
      endif
 
      y3re=2d0/(2d0*be-1d0)*f3re
      y3im=2d0/(2d0*be-1d0)*f3im
 
      return
      end
 
C***********************************************************************
 
C...PYSPEN
C...Calculates real and imaginary part of Spence function; see
C...G. 't Hooft and M. Veltman, Nucl. Phys. B153 (1979) 365.
 
      function pjspen(xrein,ximin,ireim)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      save /jydat1/
C...Local array and data.
      dimension b(0:14)
      data b/
     &1.000000d+00,        -5.000000d-01,         1.666667d-01,
     &0.000000d+00,        -3.333333d-02,         0.000000d+00,
     &2.380952d-02,         0.000000d+00,        -3.333333d-02,
     &0.000000d+00,         7.575757d-02,         0.000000d+00,
     &-2.531135d-01,         0.000000d+00,         1.166667d+00/
 
      xre=xrein
      xim=ximin
      if(abs(1d0-xre).lt.1.d-6.and.abs(xim).lt.1.d-6) then
        if(ireim.eq.1) pjspen=paru(1)**2/6d0
        if(ireim.eq.2) pjspen=0d0
        return
      endif
 
      xmod=sqrt(xre**2+xim**2)
      if(xmod.lt.1.d-6) then
        if(ireim.eq.1) pjspen=0d0
        if(ireim.eq.2) pjspen=0d0
        return
      endif
 
      xarg=sign(acos(xre/xmod),xim)
      sp0re=0d0
      sp0im=0d0
      sgn=1d0
      if(xmod.gt.1d0) then
        algxre=log(xmod)
        algxim=xarg-sign(paru(1),xarg)
        sp0re=-paru(1)**2/6d0-(algxre**2-algxim**2)/2d0
        sp0im=-algxre*algxim
        sgn=-1d0
        xmod=1d0/xmod
        xarg=-xarg
        xre=xmod*cos(xarg)
        xim=xmod*sin(xarg)
      endif
      if(xre.gt.0.5d0) then
        algxre=log(xmod)
        algxim=xarg
        xre=1d0-xre
        xim=-xim
        xmod=sqrt(xre**2+xim**2)
        xarg=sign(acos(xre/xmod),xim)
        algyre=log(xmod)
        algyim=xarg
        sp0re=sp0re+sgn*(paru(1)**2/6d0-(algxre*algyre-algxim*algyim))
        sp0im=sp0im-sgn*(algxre*algyim+algxim*algyre)
        sgn=-sgn
      endif
 
      xre=1d0-xre
      xim=-xim
      xmod=sqrt(xre**2+xim**2)
      xarg=sign(acos(xre/xmod),xim)
      zre=-log(xmod)
      zim=-xarg
 
      spre=0d0
      spim=0d0
      savere=1d0
      saveim=0d0
      do 100 i=0,14
        if(max(abs(savere),abs(saveim)).lt.1d-30) goto 110
        termre=(savere*zre-saveim*zim)/dble(i+1)
        termim=(savere*zim+saveim*zre)/dble(i+1)
        savere=termre
        saveim=termim
        spre=spre+b(i)*termre
        spim=spim+b(i)*termim
  100 continue
 
  110 if(ireim.eq.1) pjspen=sp0re+sgn*spre
      if(ireim.eq.2) pjspen=sp0im+sgn*spim
 
      return
      end
 
C***********************************************************************
 
C...PYQQBH
C...Calculates the matrix element for the processes
C...g + g or q + qbar -> Q + Qbar + H (normally with Q = t).
C...REDUCE output and part of the rest courtesy Z. Kunszt, see
C...Z. Kunszt, Nucl. Phys. B247 (1984) 339.
 
      subroutine pjqqbh(wtqqbh)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      save /jydat1/,/jydat2/,/pjpars/,/pjint1/,/pjint2/
C...Local arrays and function.
      dimension pp(15,4),clr(8,8),fm(10,10),rm(8,8),dx(8)
      dot(i,j)=pp(i,4)*pp(j,4)-pp(i,1)*pp(j,1)-pp(i,2)*pp(j,2)-
     &pp(i,3)*pp(j,3)
 
C...Mass parameters.
      wtqqbh=0d0
      isub=mint(1)
      shpr=sqrt(vint(26))*vint(1)
      pq=pmas(jamcomp(kfpr(isub,2)),1)
      ph=sqrt(vint(21))*vint(1)
      spq=pq**2
      sph=ph**2
 
C...Set up outgoing kinematics: 1=t, 2=tbar, 3=H.
      do 100 i=1,2
        pt=sqrt(max(0d0,vint(197+5*i)))
        pp(i,1)=pt*cos(vint(198+5*i))
        pp(i,2)=pt*sin(vint(198+5*i))
  100 continue
      pp(3,1)=-pp(1,1)-pp(2,1)
      pp(3,2)=-pp(1,2)-pp(2,2)
      pms1=spq+pp(1,1)**2+pp(1,2)**2
      pms2=spq+pp(2,1)**2+pp(2,2)**2
      pms3=sph+pp(3,1)**2+pp(3,2)**2
      pmt3=sqrt(pms3)
      pp(3,3)=pmt3*sinh(vint(211))
      pp(3,4)=pmt3*cosh(vint(211))
      pms12=(shpr-pp(3,4))**2-pp(3,3)**2
      pp(1,3)=(-pp(3,3)*(pms12+pms1-pms2)+
     &vint(213)*(shpr-pp(3,4))*vint(220))/(2d0*pms12)
      pp(2,3)=-pp(1,3)-pp(3,3)
      pp(1,4)=sqrt(pms1+pp(1,3)**2)
      pp(2,4)=sqrt(pms2+pp(2,3)**2)
 
C...Set up incoming kinematics and derived momentum combinations.
      do 110 i=4,5
        pp(i,1)=0d0
        pp(i,2)=0d0
        pp(i,3)=-0.5d0*shpr*(-1)**i
        pp(i,4)=-0.5d0*shpr
  110 continue
      do 120 j=1,4
        pp(6,j)=pp(1,j)+pp(2,j)
        pp(7,j)=pp(1,j)+pp(3,j)
        pp(8,j)=pp(1,j)+pp(4,j)
        pp(9,j)=pp(1,j)+pp(5,j)
        pp(10,j)=-pp(2,j)-pp(3,j)
        pp(11,j)=-pp(2,j)-pp(4,j)
        pp(12,j)=-pp(2,j)-pp(5,j)
        pp(13,j)=-pp(4,j)-pp(5,j)
  120 continue
 
C...Derived kinematics invariants.
      x1=dot(1,2)
      x2=dot(1,3)
      x3=dot(1,4)
      x4=dot(1,5)
      x5=dot(2,3)
      x6=dot(2,4)
      x7=dot(2,5)
      x8=dot(3,4)
      x9=dot(3,5)
      x10=dot(4,5)
 
C...Propagators.
      ss1=dot(7,7)-spq
      ss2=dot(8,8)-spq
      ss3=dot(9,9)-spq
      ss4=dot(10,10)-spq
      ss5=dot(11,11)-spq
      ss6=dot(12,12)-spq
      ss7=dot(13,13)
      dx(1)=ss1*ss6
      dx(2)=ss2*ss6
      dx(3)=ss2*ss4
      dx(4)=ss1*ss5
      dx(5)=ss3*ss5
      dx(6)=ss3*ss4
      dx(7)=ss7*ss1
      dx(8)=ss7*ss4
 
C...Define colour coefficients for g + g -> Q + Qbar + H.
      if(isub.eq.121.or.isub.eq.181.or.isub.eq.186) then
        do 140 i=1,3
          do 130 j=1,3
            clr(i,j)=16d0/3d0
            clr(i+3,j+3)=16d0/3d0
            clr(i,j+3)=-2d0/3d0
            clr(i+3,j)=-2d0/3d0
  130     continue
  140   continue
        do 160 l=1,2
          do 150 i=1,3
            clr(i,6+l)=-6d0
            clr(i+3,6+l)=6d0
            clr(6+l,i)=-6d0
            clr(6+l,i+3)=6d0
  150     continue
  160   continue
        do 180 k1=1,2
          do 170 k2=1,2
            clr(6+k1,6+k2)=12d0
  170     continue
  180   continue
 
C...Evaluate matrix elements for g + g -> Q + Qbar + H.
        fm(1,1)=64*pq**6+16*pq**4*ph**2+32*pq**4*(x1+2*x2+x4+x9+2*
     &  x7+x5)+8*pq**2*ph**2*(-x1-x4+2*x7)+16*pq**2*(x2*x9+4*x2*
     &  x7+x2*x5-2*x4*x7-2*x9*x7)+8*ph**2*x4*x7-16*x2*x9*x7
        fm(1,2)=16*pq**6+8*pq**4*(-2*x1+x2-2*x3-2*x4-4*x10+x9-x8+2
     &  *x7-4*x6+x5)+8*pq**2*(-2*x1*x2-2*x2*x4-2*x2*x10+x2*x7-2*
     &  x2*x6-2*x3*x7+2*x4*x7+4*x10*x7-x9*x7-x8*x7)+16*x2*x7*(x4+
     &  x10)
        fm(1,3)=16*pq**6-4*pq**4*ph**2+8*pq**4*(-2*x1+2*x2-2*x3-4*
     &  x4-8*x10+x9+x8-2*x7-4*x6+2*x5)-(4*pq**2*ph**2)*(x1+x4+x10
     &  +x6)+8*pq**2*(-2*x1*x2-2*x1*x10+x1*x9+x1*x8-2*x1*x5+x2**2
     &  -4*x2*x4-5*x2*x10+x2*x8-x2*x7-3*x2*x6+x2*x5+x3*x9+2*x3*x7
     &  -x3*x5+x4*x8+2*x4*x6-3*x4*x5-5*x10*x5+x9*x8+x9*x6+x9*x5+
     &  x8*x7-4*x6*x5+x5**2)-(16*x2*x5)*(x1+x4+x10+x6)
        fm(1,4)=16*pq**6+4*pq**4*ph**2+16*pq**4*(-x1+x2-x3-x4+x10-
     &  x9-x8+2*x7+2*x6-x5)+4*pq**2*ph**2*(x1+x3+x4+x10+2*x7+2*x6
     &  )+8*pq**2*(4*x1*x10+4*x1*x7+4*x1*x6+2*x2*x10-x2*x9-x2*x8+
     &  4*x2*x7+4*x2*x6-x2*x5+4*x10*x5+4*x7*x5+4*x6*x5)-(8*ph**2*
     &  x1)*(x10+x7+x6)+16*x2*x5*(x10+x7+x6)
        fm(1,5)=8*pq**4*(-2*x1-2*x4+x10-x9)+4*pq**2*(4*x1**2-2*x1*
     &  x2+8*x1*x3+6*x1*x10-2*x1*x9+4*x1*x8+4*x1*x7+4*x1*x6+2*x1*
     &  x5+x2*x10+4*x3*x4-x3*x9+2*x3*x7+3*x4*x8-2*x4*x6+2*x4*x5-4
     &  *x10*x7+3*x10*x5-3*x9*x6+3*x8*x7-4*x7**2+4*x7*x5)+8*(x1**
     &  2*x9-x1**2*x8-x1*x2*x7+x1*x2*x6+x1*x3*x9+x1*x3*x5-x1*x4*
     &  x8-x1*x4*x5+x1*x10*x9+x1*x9*x7+x1*x9*x6-x1*x8*x7-x2*x3*x7
     &  +x2*x4*x6-x2*x10*x7-x2*x7**2+x3*x7*x5-x4*x10*x5-x4*x7*x5-
     &  x4*x6*x5)
        fm(1,6)=16*pq**4*(-4*x1-x4+x9-x7)+4*pq**2*ph**2*(-2*x1-x4-
     &  x7)+16*pq**2*(-2*x1**2-3*x1*x2-2*x1*x4-3*x1*x9-2*x1*x7-3*
     &  x1*x5-2*x2*x4-2*x7*x5)-8*ph**2*x4*x7+8*(-x1*x2*x9-2*x1*x2
     &  *x5-x1*x9**2-x1*x9*x5+x2**2*x7-x2*x4*x5+x2*x9*x7-x2*x7*x5
     &  +x4*x9*x5+x4*x5**2)
        fm(1,7)=8*pq**4*(2*x3+x4+3*x10+x9+2*x8+3*x7+6*x6)+2*pq**2*
     &  ph**2*(-2*x3-x4+3*x10+3*x7+6*x6)+4*pq**2*(4*x1*x10+4*x1*
     &  x7+8*x1*x6+6*x2*x10+x2*x9+2*x2*x8+6*x2*x7+12*x2*x6-8*x3*
     &  x7+4*x4*x7+4*x4*x6+4*x10*x5+4*x9*x7+4*x9*x6-8*x8*x7+4*x7*
     &  x5+8*x6*x5)+4*ph**2*(-x1*x10-x1*x7-2*x1*x6+2*x3*x7-x4*x7-
     &  x4*x6)+8*x2*(x10*x5+x9*x7+x9*x6-2*x8*x7+x7*x5+2*x6*x5)
        fm(1,8)=8*pq**4*(2*x3+x4+3*x10+2*x9+x8+3*x7+6*x6)+2*pq**2*
     &  ph**2*(-2*x3-x4+2*x10+x7+2*x6)+4*pq**2*(4*x1*x10-2*x1*x9+
     &  2*x1*x8+4*x1*x7+8*x1*x6+5*x2*x10+2*x2*x9+x2*x8+4*x2*x7+8*
     &  x2*x6-x3*x9-8*x3*x7+2*x3*x5+2*x4*x9-x4*x8+4*x4*x7+4*x4*x6
     &  +4*x4*x5+5*x10*x5+x9**2-x9*x8+2*x9*x7+5*x9*x6+x9*x5-7*x8*
     &  x7+2*x8*x5+2*x7*x5+10*x6*x5)+2*ph**2*(-x1*x10+x3*x7-2*x4*
     &  x7+x4*x6)+4*(-x1*x9**2+x1*x9*x8-2*x1*x9*x5-x1*x8*x5+2*x2*
     &  x10*x5+x2*x9*x7+x2*x9*x6-2*x2*x8*x7+3*x2*x6*x5+x3*x9*x5+
     &  x3*x5**2+x4*x9*x5-2*x4*x8*x5+2*x4*x5**2)
        fm(2,2)=16*pq**6+16*pq**4*(-x1+x3-x4-x10+x7-x6)+16*pq**2*(
     &  x3*x10+x3*x7+x3*x6+x4*x7+x10*x7)-16*x3*x10*x7
        fm(2,3)=16*pq**6+8*pq**4*(-2*x1+x2+2*x3-4*x4-4*x10-x9+x8-2
     &  *x7-2*x6+x5)+8*pq**2*(-2*x1*x5+4*x3*x10-x3*x9-x3*x8-2*x3*
     &  x7+2*x3*x6+x3*x5-2*x4*x5-2*x10*x5-2*x6*x5)+16*x3*x5*(x10+
     &  x6)
        fm(2,4)=8*pq**4*(-2*x1-2*x3+x10-x8)+4*pq**2*(4*x1**2-2*x1*
     &  x2+8*x1*x4+6*x1*x10+4*x1*x9-2*x1*x8+4*x1*x7+4*x1*x6+2*x1*
     &  x5+x2*x10+4*x3*x4+3*x3*x9-2*x3*x7+2*x3*x5-x4*x8+2*x4*x6-4
     &  *x10*x6+3*x10*x5+3*x9*x6-3*x8*x7-4*x6**2+4*x6*x5)+8*(-x1
     &  **2*x9+x1**2*x8+x1*x2*x7-x1*x2*x6-x1*x3*x9-x1*x3*x5+x1*x4
     &  *x8+x1*x4*x5+x1*x10*x8-x1*x9*x6+x1*x8*x7+x1*x8*x6+x2*x3*
     &  x7-x2*x4*x6-x2*x10*x6-x2*x6**2-x3*x10*x5-x3*x7*x5-x3*x6*
     &  x5+x4*x6*x5)
        fm(2,5)=16*pq**4*x10+8*pq**2*(2*x1**2+2*x1*x3+2*x1*x4+2*x1
     &  *x10+2*x1*x7+2*x1*x6+x3*x7+x4*x6)+8*(-2*x1**3-2*x1**2*x3-
     &  2*x1**2*x4-2*x1**2*x10-2*x1**2*x7-2*x1**2*x6-2*x1*x3*x4-
     &  x1*x3*x10-2*x1*x3*x6-x1*x4*x10-2*x1*x4*x7-x1*x10**2-x1*
     &  x10*x7-x1*x10*x6-2*x1*x7*x6+x3**2*x7-x3*x4*x7-x3*x4*x6+x3
     &  *x10*x7+x3*x7**2-x3*x7*x6+x4**2*x6+x4*x10*x6-x4*x7*x6+x4*
     &  x6**2)
        fm(2,6)=8*pq**4*(-2*x1+x10-x9-2*x7)+4*pq**2*(4*x1**2+2*x1*
     &  x2+4*x1*x3+4*x1*x4+6*x1*x10-2*x1*x9+4*x1*x8+8*x1*x6-2*x1*
     &  x5+4*x2*x4+3*x2*x10+2*x2*x7-3*x3*x9-2*x3*x7-4*x4**2-4*x4*
     &  x10+3*x4*x8+2*x4*x6+x10*x5-x9*x6+3*x8*x7+4*x7*x6)+8*(x1**
     &  2*x9-x1**2*x8-x1*x2*x7+x1*x2*x6+x1*x3*x9+x1*x3*x5+x1*x4*
     &  x9-x1*x4*x8-x1*x4*x5+x1*x10*x9+x1*x9*x6-x1*x8*x7-x2*x3*x7
     &  -x2*x4*x7+x2*x4*x6-x2*x10*x7+x3*x7*x5-x4**2*x5-x4*x10*x5-
     &  x4*x6*x5)
        fm(2,7)=8*pq**4*(x3+2*x4+3*x10+x7+2*x6)+4*pq**2*(-4*x1*x3-
     &  2*x1*x4-2*x1*x10+x1*x9-x1*x8-4*x1*x7-2*x1*x6+x2*x3+2*x2*
     &  x4+3*x2*x10+x2*x7+2*x2*x6-6*x3*x4-6*x3*x10-2*x3*x9-2*x3*
     &  x7-4*x3*x6-x3*x5-6*x4**2-6*x4*x10-3*x4*x9-x4*x8-4*x4*x7-2
     &  *x4*x6-2*x4*x5-3*x10*x9-3*x10*x8-6*x10*x7-6*x10*x6+x10*x5
     &  +x9*x7-2*x8*x7-2*x8*x6-6*x7*x6+x7*x5-6*x6**2+2*x6*x5)+4*(
     &  -x1**2*x9+x1**2*x8-2*x1*x2*x10-3*x1*x2*x7-3*x1*x2*x6+x1*
     &  x3*x9-x1*x3*x5+x1*x4*x9+x1*x4*x8+x1*x4*x5+x1*x10*x9+x1*
     &  x10*x8-x1*x9*x6+x1*x8*x6+x2*x3*x7-3*x2*x4*x7-x2*x4*x6-3*
     &  x2*x10*x7-3*x2*x10*x6-3*x2*x7*x6-3*x2*x6**2-2*x3*x4*x5-x3
     &  *x10*x5-x3*x6*x5-x4**2*x5-x4*x10*x5+x4*x6*x5)
        fm(2,8)=8*pq**4*(x3+2*x4+3*x10+x7+2*x6)+4*pq**2*(-4*x1*x3-
     &  2*x1*x4-2*x1*x10-x1*x9+x1*x8-4*x1*x7-2*x1*x6+x2*x3+2*x2*
     &  x4+x2*x10-x2*x7-2*x2*x6-6*x3*x4-6*x3*x10-2*x3*x9+x3*x8-2*
     &  x3*x7-4*x3*x6+x3*x5-6*x4**2-6*x4*x10-2*x4*x9-4*x4*x7-2*x4
     &  *x6+2*x4*x5-3*x10*x9-3*x10*x8-6*x10*x7-6*x10*x6+3*x10*x5-
     &  x9*x6-2*x8*x7-3*x8*x6-6*x7*x6+x7*x5-6*x6**2+2*x6*x5)+4*(
     &  x1**2*x9-x1**2*x8-x1*x2*x7+x1*x2*x6-3*x1*x3*x5+x1*x4*x9-
     &  x1*x4*x8-3*x1*x4*x5+x1*x10*x9+x1*x10*x8-2*x1*x10*x5+x1*x9
     &  *x6+x1*x8*x7+x1*x8*x6-x2*x4*x7+x2*x4*x6-x2*x10*x7-x2*x10*
     &  x6-2*x2*x7*x6-x2*x6**2-3*x3*x4*x5-3*x3*x10*x5+x3*x7*x5-3*
     &  x3*x6*x5-3*x4**2*x5-3*x4*x10*x5-x4*x6*x5)
        fm(3,3)=64*pq**6+16*pq**4*ph**2+32*pq**4*(x1+x2+2*x3+x8+x6
     &  +2*x5)+8*pq**2*ph**2*(-x1+2*x3-x6)+16*pq**2*(x2*x5-2*x3*
     &  x8-2*x3*x6+4*x3*x5+x8*x5)+8*ph**2*x3*x6-16*x3*x8*x5
        fm(3,4)=16*pq**4*(-4*x1-x3+x8-x6)+4*pq**2*ph**2*(-2*x1-x3-
     &  x6)+16*pq**2*(-2*x1**2-3*x1*x2-2*x1*x3-3*x1*x8-2*x1*x6-3*
     &  x1*x5-2*x2*x3-2*x6*x5)-8*ph**2*x3*x6+8*(-x1*x2*x8-2*x1*x2
     &  *x5-x1*x8**2-x1*x8*x5+x2**2*x6-x2*x3*x5+x2*x8*x6-x2*x6*x5
     &  +x3*x8*x5+x3*x5**2)
        fm(3,5)=8*pq**4*(-2*x1+x10-x8-2*x6)+4*pq**2*(4*x1**2+2*x1*
     &  x2+4*x1*x3+4*x1*x4+6*x1*x10+4*x1*x9-2*x1*x8+8*x1*x7-2*x1*
     &  x5+4*x2*x3+3*x2*x10+2*x2*x6-4*x3**2-4*x3*x10+3*x3*x9+2*x3
     &  *x7-3*x4*x8-2*x4*x6+x10*x5+3*x9*x6-x8*x7+4*x7*x6)+8*(-x1
     &  **2*x9+x1**2*x8+x1*x2*x7-x1*x2*x6-x1*x3*x9+x1*x3*x8-x1*x3
     &  *x5+x1*x4*x8+x1*x4*x5+x1*x10*x8-x1*x9*x6+x1*x8*x7+x2*x3*
     &  x7-x2*x3*x6-x2*x4*x6-x2*x10*x6-x3**2*x5-x3*x10*x5-x3*x7*
     &  x5+x4*x6*x5)
        fm(3,6)=16*pq**6+4*pq**4*ph**2+16*pq**4*(-x1-x2+2*x3+2*x4+
     &  x10-x9-x8-x7-x6+x5)+4*pq**2*ph**2*(x1+2*x3+2*x4+x10+x7+x6
     &  )+8*pq**2*(4*x1*x3+4*x1*x4+4*x1*x10+4*x2*x3+4*x2*x4+4*x2*
     &  x10-x2*x5+4*x3*x5+4*x4*x5+2*x10*x5-x9*x5-x8*x5)-(8*ph**2*
     &  x1)*(x3+x4+x10)+16*x2*x5*(x3+x4+x10)
        fm(3,7)=8*pq**4*(3*x3+6*x4+3*x10+x9+2*x8+2*x7+x6)+2*pq**2*
     &  ph**2*(x3+2*x4+2*x10-2*x7-x6)+4*pq**2*(4*x1*x3+8*x1*x4+4*
     &  x1*x10+2*x1*x9-2*x1*x8+2*x2*x3+10*x2*x4+5*x2*x10+2*x2*x9+
     &  x2*x8+2*x2*x7+4*x2*x6-7*x3*x9+2*x3*x8-8*x3*x7+4*x3*x6+4*
     &  x3*x5+5*x4*x8+4*x4*x6+8*x4*x5+5*x10*x5-x9*x8-x9*x6+x9*x5+
     &  x8**2-x8*x7+2*x8*x6+2*x8*x5)+2*ph**2*(-x1*x10+x3*x7-2*x3*
     &  x6+x4*x6)+4*(-x1*x2*x9-2*x1*x2*x8+x1*x9*x8-x1*x8**2+x2**2
     &  *x7+2*x2**2*x6+3*x2*x4*x5+2*x2*x10*x5-2*x2*x9*x6+x2*x8*x7
     &  +x2*x8*x6-2*x3*x9*x5+x3*x8*x5+x4*x8*x5)
        fm(3,8)=8*pq**4*(3*x3+6*x4+3*x10+2*x9+x8+2*x7+x6)+2*pq**2*
     &  ph**2*(3*x3+6*x4+3*x10-2*x7-x6)+4*pq**2*(4*x1*x3+8*x1*x4+
     &  4*x1*x10+4*x2*x3+8*x2*x4+4*x2*x10-8*x3*x9+4*x3*x8-8*x3*x7
     &  +4*x3*x6+6*x3*x5+4*x4*x8+4*x4*x6+12*x4*x5+6*x10*x5+2*x9*
     &  x5+x8*x5)+4*ph**2*(-x1*x3-2*x1*x4-x1*x10+2*x3*x7-x3*x6-x4
     &  *x6)+8*x5*(x2*x3+2*x2*x4+x2*x10-2*x3*x9+x3*x8+x4*x8)
        fm(4,4)=64*pq**6+16*pq**4*ph**2+32*pq**4*(x1+2*x2+x3+x8+2*
     &  x6+x5)+8*pq**2*ph**2*(-x1-x3+2*x6)+16*pq**2*(x2*x8+4*x2*
     &  x6+x2*x5-2*x3*x6-2*x8*x6)+8*ph**2*x3*x6-16*x2*x8*x6
        fm(4,5)=16*pq**6+8*pq**4*(-2*x1+x2-2*x3-2*x4-4*x10-x9+x8-4
     &  *x7+2*x6+x5)+8*pq**2*(-2*x1*x2-2*x2*x3-2*x2*x10-2*x2*x7+
     &  x2*x6+2*x3*x6-2*x4*x6+4*x10*x6-x9*x6-x8*x6)+16*x2*x6*(x3+
     &  x10)
        fm(4,6)=16*pq**6-4*pq**4*ph**2+8*pq**4*(-2*x1+2*x2-4*x3-2*
     &  x4-8*x10+x9+x8-4*x7-2*x6+2*x5)-(4*pq**2*ph**2)*(x1+x3+x10
     &  +x7)+8*pq**2*(-2*x1*x2-2*x1*x10+x1*x9+x1*x8-2*x1*x5+x2**2
     &  -4*x2*x3-5*x2*x10+x2*x9-3*x2*x7-x2*x6+x2*x5+x3*x9+2*x3*x7
     &  -3*x3*x5+x4*x8+2*x4*x6-x4*x5-5*x10*x5+x9*x8+x9*x6+x8*x7+
     &  x8*x5-4*x7*x5+x5**2)-(16*x2*x5)*(x1+x3+x10+x7)
        fm(4,7)=8*pq**4*(-x3-2*x4-3*x10-2*x9-x8-6*x7-3*x6)+2*pq**2
     &  *ph**2*(x3+2*x4-3*x10-6*x7-3*x6)+4*pq**2*(-4*x1*x10-8*x1*
     &  x7-4*x1*x6-6*x2*x10-2*x2*x9-x2*x8-12*x2*x7-6*x2*x6-4*x3*
     &  x7-4*x3*x6+8*x4*x6-4*x10*x5+8*x9*x6-4*x8*x7-4*x8*x6-8*x7*
     &  x5-4*x6*x5)+4*ph**2*(x1*x10+2*x1*x7+x1*x6+x3*x7+x3*x6-2*
     &  x4*x6)+8*x2*(-x10*x5+2*x9*x6-x8*x7-x8*x6-2*x7*x5-x6*x5)
        fm(4,8)=8*pq**4*(-x3-2*x4-3*x10-x9-2*x8-6*x7-3*x6)+2*pq**2
     &  *ph**2*(x3+2*x4-2*x10-2*x7-x6)+4*pq**2*(-4*x1*x10-2*x1*x9
     &  +2*x1*x8-8*x1*x7-4*x1*x6-5*x2*x10-x2*x9-2*x2*x8-8*x2*x7-4
     &  *x2*x6+x3*x9-2*x3*x8-4*x3*x7-4*x3*x6-4*x3*x5+x4*x8+8*x4*
     &  x6-2*x4*x5-5*x10*x5+x9*x8+7*x9*x6-2*x9*x5-x8**2-5*x8*x7-2
     &  *x8*x6-x8*x5-10*x7*x5-2*x6*x5)+2*ph**2*(x1*x10-x3*x7+2*x3
     &  *x6-x4*x6)+4*(-x1*x9*x8+x1*x9*x5+x1*x8**2+2*x1*x8*x5-2*x2
     &  *x10*x5+2*x2*x9*x6-x2*x8*x7-x2*x8*x6-3*x2*x7*x5+2*x3*x9*
     &  x5-x3*x8*x5-2*x3*x5**2-x4*x8*x5-x4*x5**2)
        fm(5,5)=16*pq**6+16*pq**4*(-x1-x3+x4-x10-x7+x6)+16*pq**2*(
     &  x3*x6+x4*x10+x4*x7+x4*x6+x10*x6)-16*x4*x10*x6
        fm(5,6)=16*pq**6+8*pq**4*(-2*x1+x2-4*x3+2*x4-4*x10+x9-x8-2
     &  *x7-2*x6+x5)+8*pq**2*(-2*x1*x5-2*x3*x5+4*x4*x10-x4*x9-x4*
     &  x8+2*x4*x7-2*x4*x6+x4*x5-2*x10*x5-2*x7*x5)+16*x4*x5*(x10+
     &  x7)
        fm(5,7)=8*pq**4*(-2*x3-x4-3*x10-2*x7-x6)+4*pq**2*(2*x1*x3+
     &  4*x1*x4+2*x1*x10+x1*x9-x1*x8+2*x1*x7+4*x1*x6-2*x2*x3-x2*
     &  x4-3*x2*x10-2*x2*x7-x2*x6+6*x3**2+6*x3*x4+6*x3*x10+x3*x9+
     &  3*x3*x8+2*x3*x7+4*x3*x6+2*x3*x5+6*x4*x10+2*x4*x8+4*x4*x7+
     &  2*x4*x6+x4*x5+3*x10*x9+3*x10*x8+6*x10*x7+6*x10*x6-x10*x5+
     &  2*x9*x7+2*x9*x6-x8*x6+6*x7**2+6*x7*x6-2*x7*x5-x6*x5)+4*(-
     &  x1**2*x9+x1**2*x8+2*x1*x2*x10+3*x1*x2*x7+3*x1*x2*x6-x1*x3
     &  *x9-x1*x3*x8-x1*x3*x5-x1*x4*x8+x1*x4*x5-x1*x10*x9-x1*x10*
     &  x8-x1*x9*x7+x1*x8*x7+x2*x3*x7+3*x2*x3*x6-x2*x4*x6+3*x2*
     &  x10*x7+3*x2*x10*x6+3*x2*x7**2+3*x2*x7*x6+x3**2*x5+2*x3*x4
     &  *x5+x3*x10*x5-x3*x7*x5+x4*x10*x5+x4*x7*x5)
        fm(5,8)=8*pq**4*(-2*x3-x4-3*x10-2*x7-x6)+4*pq**2*(2*x1*x3+
     &  4*x1*x4+2*x1*x10-x1*x9+x1*x8+2*x1*x7+4*x1*x6-2*x2*x3-x2*
     &  x4-x2*x10+2*x2*x7+x2*x6+6*x3**2+6*x3*x4+6*x3*x10+2*x3*x8+
     &  2*x3*x7+4*x3*x6-2*x3*x5+6*x4*x10-x4*x9+2*x4*x8+4*x4*x7+2*
     &  x4*x6-x4*x5+3*x10*x9+3*x10*x8+6*x10*x7+6*x10*x6-3*x10*x5+
     &  3*x9*x7+2*x9*x6+x8*x7+6*x7**2+6*x7*x6-2*x7*x5-x6*x5)+4*(
     &  x1**2*x9-x1**2*x8-x1*x2*x7+x1*x2*x6+x1*x3*x9-x1*x3*x8+3*
     &  x1*x3*x5+3*x1*x4*x5-x1*x10*x9-x1*x10*x8+2*x1*x10*x5-x1*x9
     &  *x7-x1*x9*x6-x1*x8*x7-x2*x3*x7+x2*x3*x6+x2*x10*x7+x2*x10*
     &  x6+x2*x7**2+2*x2*x7*x6+3*x3**2*x5+3*x3*x4*x5+3*x3*x10*x5+
     &  x3*x7*x5+3*x4*x10*x5+3*x4*x7*x5-x4*x6*x5)
        fm(6,6)=64*pq**6+16*pq**4*ph**2+32*pq**4*(x1+x2+2*x4+x9+x7
     &  +2*x5)+8*pq**2*ph**2*(-x1+2*x4-x7)+16*pq**2*(x2*x5-2*x4*
     &  x9-2*x4*x7+4*x4*x5+x9*x5)+8*ph**2*x4*x7-16*x4*x9*x5
        fm(6,7)=8*pq**4*(-6*x3-3*x4-3*x10-2*x9-x8-x7-2*x6)+2*pq**2
     &  *ph**2*(-2*x3-x4-2*x10+x7+2*x6)+4*pq**2*(-8*x1*x3-4*x1*x4
     &  -4*x1*x10+2*x1*x9-2*x1*x8-10*x2*x3-2*x2*x4-5*x2*x10-x2*x9
     &  -2*x2*x8-4*x2*x7-2*x2*x6-5*x3*x9-4*x3*x7-8*x3*x5-2*x4*x9+
     &  7*x4*x8-4*x4*x7+8*x4*x6-4*x4*x5-5*x10*x5-x9**2+x9*x8-2*x9
     &  *x7+x9*x6-2*x9*x5+x8*x7-x8*x5)+2*ph**2*(x1*x10-x3*x7+2*x4
     &  *x7-x4*x6)+4*(2*x1*x2*x9+x1*x2*x8+x1*x9**2-x1*x9*x8-2*x2
     &  **2*x7-x2**2*x6-3*x2*x3*x5-2*x2*x10*x5-x2*x9*x7-x2*x9*x6+
     &  2*x2*x8*x7-x3*x9*x5-x4*x9*x5+2*x4*x8*x5)
        fm(6,8)=8*pq**4*(-6*x3-3*x4-3*x10-x9-2*x8-x7-2*x6)+2*pq**2
     &  *ph**2*(-6*x3-3*x4-3*x10+x7+2*x6)+4*pq**2*(-8*x1*x3-4*x1*
     &  x4-4*x1*x10-8*x2*x3-4*x2*x4-4*x2*x10-4*x3*x9-4*x3*x7-12*
     &  x3*x5-4*x4*x9+8*x4*x8-4*x4*x7+8*x4*x6-6*x4*x5-6*x10*x5-x9
     &  *x5-2*x8*x5)+4*ph**2*(2*x1*x3+x1*x4+x1*x10+x3*x7+x4*x7-2*
     &  x4*x6)+8*x5*(-2*x2*x3-x2*x4-x2*x10-x3*x9-x4*x9+2*x4*x8)
        fm(7,7)=72*pq**4*x10+18*pq**2*ph**2*x10+8*pq**2*(x1*x10+9*
     &  x2*x10+7*x3*x7+2*x3*x6+2*x4*x7+7*x4*x6+x10*x5+2*x9*x7+7*
     &  x9*x6+7*x8*x7+2*x8*x6)+2*ph**2*(-x1*x10-7*x3*x7-2*x3*x6-2
     &  *x4*x7-7*x4*x6)+4*x2*(x10*x5+2*x9*x7+7*x9*x6+7*x8*x7+2*x8
     &  *x6)
        fm(7,8)=72*pq**4*x10+2*pq**2*ph**2*x10+4*pq**2*(2*x1*x10+
     &  10*x2*x10+7*x3*x9+2*x3*x8+14*x3*x7+4*x3*x6+2*x4*x9+7*x4*
     &  x8+4*x4*x7+14*x4*x6+10*x10*x5+x9**2+7*x9*x8+2*x9*x7+7*x9*
     &  x6+x8**2+7*x8*x7+2*x8*x6)+2*ph**2*(7*x1*x10-7*x3*x7-2*x3*
     &  x6-2*x4*x7-7*x4*x6)+2*(-2*x1*x9**2-14*x1*x9*x8-2*x1*x8**2
     &  +2*x2*x10*x5+2*x2*x9*x7+7*x2*x9*x6+7*x2*x8*x7+2*x2*x8*x6+
     &  7*x3*x9*x5+2*x3*x8*x5+2*x4*x9*x5+7*x4*x8*x5)
        fm(8,8)=72*pq**4*x10+18*pq**2*ph**2*x10+8*pq**2*(x1*x10+x2
     &  *x10+7*x3*x9+2*x3*x8+7*x3*x7+2*x3*x6+2*x4*x9+7*x4*x8+2*x4
     &  *x7+7*x4*x6+9*x10*x5)+2*ph**2*(-x1*x10-7*x3*x7-2*x3*x6-2*
     &  x4*x7-7*x4*x6)+4*x5*(x2*x10+7*x3*x9+2*x3*x8+2*x4*x9+7*x4*
     &  x8)
        fm(9,9)=-4*pq**4*x10-pq**2*ph**2*x10+4*pq**2*(-x1*x10-x2*x10+
     &  x3*x7+x4*x6-x10*x5+x9*x6+x8*x7)+ph**2*(x1*x10-x3*x7-x4*x6
     &  )+2*x2*(-x10*x5+x9*x6+x8*x7)
        fm(9,10)=-4*pq**4*x10-pq**2*ph**2*x10+2*pq**2*(-2*x1*x10-2*x2*
     &  x10+2*x3*x9+2*x3*x7+2*x4*x6-2*x10*x5+x9*x8+2*x8*x7)+ph**2
     &  *(x1*x10-x3*x7-x4*x6)+2*(-x1*x9*x8-x2*x10*x5+x2*x8*x7+x3*
     &  x9*x5)
        fmxx=-4*pq**4*x10-pq**2*ph**2*x10+2*pq**2*(-2*x1*x10-2*x2*
     &  x10+2*x4*x8+2*x4*x6+2*x3*x7-2*x10*x5+x9*x8+2*x9*x6)+ph**2
     &  *(x1*x10-x3*x7-x4*x6)+2*(-x1*x9*x8-x2*x10*x5+x2*x9*x6+x4*
     &  x8*x5)
        fm(9,10)=0.5d0*(fmxx+fm(9,10))
        fm(10,10)=-4*pq**4*x10-pq**2*ph**2*x10+4*pq**2*(-x1*x10-x2*x10+
     &  x3*x7+x4*x6-x10*x5+x9*x3+x8*x4)+ph**2*(x1*x10-x3*x7-x4*x6
     &  )+2*x5*(-x10*x2+x9*x3+x8*x4)
 
C...Repackage matrix elements.
        do 200 i=1,8
          do 190 j=1,8
            rm(i,j)=fm(i,j)
  190     continue
  200   continue
        rm(7,7)=fm(7,7)-2d0*fm(9,9)
        rm(7,8)=fm(7,8)-2d0*fm(9,10)
        rm(8,8)=fm(8,8)-2d0*fm(10,10)
 
C...Produce final result: matrix elements * colours * propagators.
        do 220 i=1,8
          do 210 j=i,8
            fac=8d0
            if(i.eq.j)fac=4d0
            wtqqbh=wtqqbh+rm(i,j)*fac*clr(i,j)/(dx(i)*dx(j))
  210     continue
  220   continue
        wtqqbh=-wtqqbh/256d0
 
      else
C...Evaluate matrix elements for q + qbar -> Q + Qbar + H.
        a11=-8d0*pq**4*x10-2d0*pq**2*ph**2*x10-(8d0*pq**2)*(x2*x10+x3
     &  *x7+x4*x6+x9*x6+x8*x7)+2d0*ph**2*(x3*x7+x4*x6)-(4d0*x2)*(x9
     &  *x6+x8*x7)
        a12=-8d0*pq**4*x10+4d0*pq**2*(-x2*x10-x3*x9-2d0*x3*x7-x4*x8-
     &  2d0*x4*x6-x10*x5-x9*x8-x9*x6-x8*x7)+2d0*ph**2*(-x1*x10+x3*x7
     &  +x4*x6)+2d0*(2d0*x1*x9*x8-x2*x9*x6-x2*x8*x7-x3*x9*x5-x4*x8*
     &  x5)
        a22=-8d0*pq**4*x10-2d0*pq**2*ph**2*x10-(8d0*pq**2)*(x3*x9+x3*
     &  x7+x4*x8+x4*x6+x10*x5)+2d0*ph**2*(x3*x7+x4*x6)-(4d0*x5)*(x3
     &  *x9+x4*x8)
 
C...Produce final result: matrix elements * propagators.
        a11=a11/dx(7)**2
        a12=a12/(dx(7)*dx(8))
        a22=a22/dx(8)**2
        wtqqbh=-(a11+a22+2d0*a12)/8d0
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYKCUT
C...Dummy routine, which the user can replace in order to make cuts on
C...the kinematics on the parton level before the matrix elements are
C...evaluated and the event is generated. The cross-section estimates
C...will automatically take these cuts into account, so the given
C...values are for the allowed phase space region only. MCUT=0 means
C...that the event has passed the cuts, MCUT=1 that it has failed.
 
      subroutine pjkcut(mcut)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      save /jydat1/,/pjint1/,/pjint2/
 
C...Set default value (accepting event) for MCUT.
      mcut=0
 
C...Read out subprocess number.
      isub=mint(1)
      istsb=iset(isub)
 
C...Read out tau, y*, cos(theta), tau' (where defined, else =0).
      tau=vint(21)
      yst=vint(22)
      cth=0d0
      if(istsb.eq.2.or.istsb.eq.4) cth=vint(23)
      taup=0d0
      if(istsb.ge.3.and.istsb.le.5) taup=vint(26)
 
C...Calculate x_1, x_2, x_F.
      if(istsb.le.2.or.istsb.ge.5) then
        x1=sqrt(tau)*exp(yst)
        x2=sqrt(tau)*exp(-yst)
      else
        x1=sqrt(taup)*exp(yst)
        x2=sqrt(taup)*exp(-yst)
      endif
      xf=x1-x2
 
C...Calculate shat, that, uhat, p_T^2.
      shat=tau*vint(2)
      sqm3=vint(63)
      sqm4=vint(64)
      rm3=sqm3/shat
      rm4=sqm4/shat
      be34=sqrt(max(0d0,(1d0-rm3-rm4)**2-4d0*rm3*rm4))
      rpts=4d0*vint(71)**2/shat
      be34l=sqrt(max(0d0,(1d0-rm3-rm4)**2-4d0*rm3*rm4-rpts))
      rm34=2d0*rm3*rm4
      rsqm=1d0+rm34
      rthm=(4d0*rm3*rm4+rpts)/(1d0-rm3-rm4+be34l)
      that=-0.5d0*shat*max(rthm,1d0-rm3-rm4-be34*cth)
      uhat=-0.5d0*shat*max(rthm,1d0-rm3-rm4+be34*cth)
      pt2=max(vint(71)**2,0.25d0*shat*be34**2*(1d0-cth**2))
 
C...Decisions by user to be put here.
 
C...Stop program if this routine is ever called.
C...You should not copy these lines to your own routine.
      write(mstu(11),5000)
      if(pjr(0).lt.10d0) stop
 
C...Format for error printout.
 5000 format(1x,'Error: you did not link your PYKCUT routine ',
     &'correctly.'/1x,'Dummy routine in PYTHIA file called instead.'/
     &1x,'Execution stopped!')
 
      return
      end
 
C*********************************************************************
 
C...PYEVWT
C...Dummy routine, which the user can replace in order to multiply the
C...standard PYTHIA differential cross-section by a process- and
C...kinematics-dependent factor WTXS. For MSTP(142)=1 this corresponds
C...to generation of weighted events, with weight 1/WTXS, while for
C...MSTP(142)=2 it corresponds to a modification of the underlying
C...physics.
 
      subroutine pjevwt(wtxs)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/pjint1/mint(400),vint(400)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      save /jydat1/,/pjint1/,/pjint2/
 
C...Set default weight for WTXS.
      wtxs=1d0
 
C...Read out subprocess number.
      isub=mint(1)
      istsb=iset(isub)
 
C...Read out tau, y*, cos(theta), tau' (where defined, else =0).
      tau=vint(21)
      yst=vint(22)
      cth=0d0
      if(istsb.eq.2.or.istsb.eq.4) cth=vint(23)
      taup=0d0
      if(istsb.ge.3.and.istsb.le.5) taup=vint(26)
 
C...Read out x_1, x_2, x_F, shat, that, uhat, p_T^2.
      x1=vint(41)
      x2=vint(42)
      xf=x1-x2
      shat=vint(44)
      that=vint(45)
      uhat=vint(46)
      pt2=vint(48)
 
C...Modifications by user to be put here.
 
C...Stop program if this routine is ever called.
C...You should not copy these lines to your own routine.
      write(mstu(11),5000)
      if(pjr(0).lt.10d0) stop
 
C...Format for error printout.
 5000 format(1x,'Error: you did not link your PYEVWT routine ',
     &'correctly.'/1x,'Dummy routine in PYTHIA file called instead.'/
     &1x,'Execution stopped!')
 
      return
      end
 
C*********************************************************************
 
C...PYUPIN
C...Dummy copy of routine to be called by user to set up a user-defined
C...process.
 
      subroutine pjupin(isub,title,sigmax)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/pjint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      common/pjint6/proc(0:500)
      character proc*28
      save /jydat1/,/pjint2/,/pjint6/
C...Local character variable.
      character*(*) title
 
C...Check that subprocess number free.
      if(isub.lt.1.or.isub.gt.500.or.iset(isub).ge.0) then
        write(mstu(11),5000) isub
        stop
      endif
 
C...Fill information on new process.
      iset(isub)=11
      coef(isub,1)=sigmax
      proc(isub)=title//' '
 
C...Format for error output.
 5000 format(1x,'Error: user-defined subprocess code ',i4,
     &' not allowed.'//1x,'Execution stopped!')
 
      return
      end
 
C*********************************************************************
 
C...PYUPEV
C...Dummy routine, to be replaced by user. When called from PYTHIA
C...the subprocess number ISUB will be given, and PYUPEV is supposed
C...to generate an event of this type, to be stored in the PYUPPR
C...commonblock. SIGEV gives the differential cross-section associated
C...with the event, i.e. the acceptance probability of the event is
C...taken to be SIGEV/SIGMAX, where SIGMAX was given in the PYUPIN
C...call.
 
      subroutine pjupev(isub,sigev)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/pjuppr/nup,kup(20,7),nfup,ifup(10,2),pup(20,5),q2up(0:10)
      save /jydat1/,/pjuppr/
 
C...Stop program if this routine is ever called.
C...You should not copy these lines to your own routine.
      write(mstu(11),5000)
      if(pjr(0).lt.10d0) stop
      sigev=isub
 
C...Format for error printout.
 5000 format(1x,'Error: you did not link your PYUPEV routine ',
     &'correctly.'/1x,'Dummy routine in PYTHIA file called instead.'/
     &1x,'Execution stopped!')
 
      return
      end
 
 
C*********************************************************************
 
C...PYTAUD
C...Dummy routine, to be replaced by user, to handle the decay of a
C...polarized tau lepton.
C...Input:
C...ITAU is the position where the decaying tau is stored in /PYJETS/.
C...IORIG is the position where the mother of the tau is stored;
C...     is 0 when the mother is not stored.
C...KFORIG is the flavour of the mother of the tau;
C...     is 0 when the mother is not known.
C...Note that IORIG=0 does not necessarily imply KFORIG=0;
C...     e.g. in B hadron semileptonic decays the W  propagator
C...     is not explicitly stored but the W code is still unambiguous.
C...Output:
C...NDECAY is the number of decay products in the current tau decay.
C...These decay products should be added to the /PYJETS/ common block,
C...in positions N+1 through N+NDECAY. For each product I you must
C...give the flavour codes K(I,2) and the five-momenta P(I,1), P(I,2),
C...P(I,3), P(I,4) and P(I,5). The rest will be stored automatically.
 
      subroutine pjtaud(itau,iorig,kforig,ndecay)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      save /jyjets/,/jydat1/
 
C...Stop program if this routine is ever called.
C...You should not copy these lines to your own routine.
      ndecay=itau+iorig+kforig
      write(mstu(11),5000)
      if(pjr(0).lt.10d0) stop
 
C...Format for error printout.
 5000 format(1x,'Error: you did not link your PYTAUD routine ',
     &'correctly.'/1x,'Dummy routine in PYTHIA file called instead.'/
     &1x,'Execution stopped!')
 
      return
      end
 
