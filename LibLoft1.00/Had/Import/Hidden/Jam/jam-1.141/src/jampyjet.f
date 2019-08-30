C*********************************************************************
C*********************************************************************
C                                                                    *
C           PYTHIA jetset part: modified for JAM                     *
C                                                                    *
C  S   PY1ENT   to fill one entry (= parton or particle)             *
C  S   PY2ENT   to fill two entries                                  *
C  S   PY3ENT   to fill three entries                                *
C  S   PY4ENT   to fill four entries                                 *
C  S   PYJOIN   to connect entries with colour flow information      *
C  S   PYGIVE   to fill (or query) commonblock variables             *
C  S   PYEXEC   to administrate fragmentation and decay chain        *
C  S   PYPREP   to rearrange showered partons along strings          *
C  S   PYSTRF   to do string fragmentation of jet system             *
C  S   PYINDF   to do independent fragmentation of one or many jets  *
C  S   PYDECY   to do the decay of a particle                        *
C  S   PYDCYK   to select parton and hadron flavours in decays       *
C  S   PYKFDI   to select parton and hadron flavours in fragm        *
C  S   PYNMES   to select number of popcorn mesons                   *
C  S   PYKFIN   to calculate falvour prod. ratios from input params. *
C  S   PYPTDI   to select transverse momenta in fragm                *
C  S   PYZDIS   to select longitudinal scaling variable in fragm     *
C  S   PYSHOW   to do timelike parton shower evolution               *
C  S   PYBOEI   to include Bose-Einstein effects (crudely)           *
C                                                                    *
C  F   PYMASS   to give the mass of a particle or parton             *
C  S   PYNAME   to give the name of a particle or parton             *
C  S   PYERRM   to write error messages and abort faulty run         *
C  F   PYALEM   to give the alpha_electromagnetic value              *
C  F   PYALPS   to give the alpha_strong value                       *
C  F   PYANGL   to give the angle from known x and y components      *
C  F   PYRND    to provide a random number generator                 *
C  S   PYRGET   to save the state of the random number generator     *
C  S   PYRSET   to set the state of the random number generator      *
C  S   PYROBO   to rotate and/or boost an event                      *
C  S   PYEDIT   to remove unwanted entries from record               *
C  S   PYLIST   to list event record or particle data                *
C  S   PYLOGO   to write a logo                                      *
C  S   PYUPDA   to update particle data                              *
C  F   PYK      to provide integer-valued event information          *
C  F   PYP      to provide real-valued event information             *
C                                                                    *
C  S   PYEEVT   to administrate the generation of an e+e- event      *
C  S   PYXTEE   to give the total cross-section at given CM energy   *
C  S   PYRADK   to generate initial state photon radiation           *
C  S   PYXKFL   to select flavour of primary qqbar pair              *
C  S   PYXJET   to select (matrix element) jet multiplicity          *
C  S   PYX3JT   to select kinematics of three-jet event              *
C  S   PYX4JT   to select kinematics of four-jet event               *
C  S   PYXDIF   to select angular orientation of event               *
C  S   PYONIA   to perform generation of onium decay to gluons       *
c  f   pyr      to interface for random number generator for pythia  *
C                                                                    *
C                                                                    *
C*********************************************************************
 
C...PY1ENT
C...Stores one parton/particle in commonblock PYJETS.
 
      subroutine pj1ent(ip,kf,pe,the,phi)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jyjets/,/jydat1/,/jydat2/
 
C...Standard checks.
      mstu(28)=0
      if(mstu(12).ge.1) call pjlist(0)
      ipa=max(1,iabs(ip))
      if(ipa.gt.mstu(4)) call pjerrm(21,
     &'(PY1ENT:) writing outside PYJETS memory')
      kc=jamcomp(kf)
      if(kc.eq.0) call pjerrm(12,'(PY1ENT:) unknown flavour code')
 
C...Find mass. Reset K, P and V vectors.
      pm=0d0
      if(mstu(10).eq.1) pm=p(ipa,5)
      if(mstu(10).ge.2) pm=pjmass(kf)
      do 100 j=1,5
        k(ipa,j)=0
        p(ipa,j)=0d0
        v(ipa,j)=0d0
  100 continue
 
C...Store parton/particle in K and P vectors.
      k(ipa,1)=1
      if(ip.lt.0) k(ipa,1)=2
      k(ipa,2)=kf
      p(ipa,5)=pm
      p(ipa,4)=max(pe,pm)
      pa=sqrt(p(ipa,4)**2-p(ipa,5)**2)
      p(ipa,1)=pa*sin(the)*cos(phi)
      p(ipa,2)=pa*sin(the)*sin(phi)
      p(ipa,3)=pa*cos(the)
 
C...Set N. Optionally fragment/decay.
      n=ipa
      if(ip.eq.0) call pjexec
 
      return
      end
 
C*********************************************************************
 
C...PY2ENT
C...Stores two partons/particles in their CM frame,
C...with the first along the +z axis.
 
      subroutine pj2ent(ip,kf1,kf2,pecm)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jyjets/,/jydat1/,/jydat2/
 
C...Standard checks.
      mstu(28)=0
      if(mstu(12).ge.1) call pjlist(0)
      ipa=max(1,iabs(ip))
      if(ipa.gt.mstu(4)-1) call pjerrm(21,
     &'(PY2ENT:) writing outside PYJETS memory')
      kc1=jamcomp(kf1)
      kc2=jamcomp(kf2)
      if(kc1.eq.0.or.kc2.eq.0) call pjerrm(12,
     &'(PY2ENT:) unknown flavour code')
 
C...Find masses. Reset K, P and V vectors.
      pm1=0d0
      if(mstu(10).eq.1) pm1=p(ipa,5)
      if(mstu(10).ge.2) pm1=pjmass(kf1)
      pm2=0d0
      if(mstu(10).eq.1) pm2=p(ipa+1,5)
      if(mstu(10).ge.2) pm2=pjmass(kf2)
      do 110 i=ipa,ipa+1
        do 100 j=1,5
          k(i,j)=0
          p(i,j)=0d0
          v(i,j)=0d0
  100   continue
  110 continue
 
C...Check flavours.
      kq1=kchg(kc1,2)*isign(1,kf1)
      kq2=kchg(kc2,2)*isign(1,kf2)
      if(mstu(19).eq.1) then
        mstu(19)=0
      else
        if(kq1+kq2.ne.0.and.kq1+kq2.ne.4) call pjerrm(2,
     &  '(PY2ENT:) unphysical flavour combination')
      endif
      k(ipa,2)=kf1
      k(ipa+1,2)=kf2
 
C...Store partons/particles in K vectors for normal case.
      if(ip.ge.0) then
        k(ipa,1)=1
        if(kq1.ne.0.and.kq2.ne.0) k(ipa,1)=2
        k(ipa+1,1)=1
 
C...Store partons in K vectors for parton shower evolution.
      else
        k(ipa,1)=3
        k(ipa+1,1)=3
        k(ipa,4)=mstu(5)*(ipa+1)
        k(ipa,5)=k(ipa,4)
        k(ipa+1,4)=mstu(5)*ipa
        k(ipa+1,5)=k(ipa+1,4)
      endif
 
C...Check kinematics and store partons/particles in P vectors.
      if(pecm.le.pm1+pm2) call pjerrm(13,
     &'(PY2ENT:) energy smaller than sum of masses')
      pa=sqrt(max(0d0,(pecm**2-pm1**2-pm2**2)**2-(2d0*pm1*pm2)**2))/
     &(2d0*pecm)
      p(ipa,3)=pa
      p(ipa,4)=sqrt(pm1**2+pa**2)
      p(ipa,5)=pm1
      p(ipa+1,3)=-pa
      p(ipa+1,4)=sqrt(pm2**2+pa**2)
      p(ipa+1,5)=pm2
 
C...Set N. Optionally fragment/decay.
      n=ipa+1
      if(ip.eq.0) call pjexec
 
      return
      end
 
C*********************************************************************
 
C...PY3ENT
C...Stores three partons or particles in their CM frame,
C...with the first along the +z axis and the third in the (x,z)
C...plane with x > 0.
 
      subroutine pj3ent(ip,kf1,kf2,kf3,pecm,x1,x3)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jyjets/,/jydat1/,/jydat2/
 
C...Standard checks.
      mstu(28)=0
      if(mstu(12).ge.1) call pjlist(0)
      ipa=max(1,iabs(ip))
      if(ipa.gt.mstu(4)-2) call pjerrm(21,
     &'(PY3ENT:) writing outside PYJETS memory')
      kc1=jamcomp(kf1)
      kc2=jamcomp(kf2)
      kc3=jamcomp(kf3)
      if(kc1.eq.0.or.kc2.eq.0.or.kc3.eq.0) call pjerrm(12,
     &'(PY3ENT:) unknown flavour code')
 
C...Find masses. Reset K, P and V vectors.
      pm1=0d0
      if(mstu(10).eq.1) pm1=p(ipa,5)
      if(mstu(10).ge.2) pm1=pjmass(kf1)
      pm2=0d0
      if(mstu(10).eq.1) pm2=p(ipa+1,5)
      if(mstu(10).ge.2) pm2=pjmass(kf2)
      pm3=0d0
      if(mstu(10).eq.1) pm3=p(ipa+2,5)
      if(mstu(10).ge.2) pm3=pjmass(kf3)
      do 110 i=ipa,ipa+2
        do 100 j=1,5
          k(i,j)=0
          p(i,j)=0d0
          v(i,j)=0d0
  100   continue
  110 continue
 
C...Check flavours.
      kq1=kchg(kc1,2)*isign(1,kf1)
      kq2=kchg(kc2,2)*isign(1,kf2)
      kq3=kchg(kc3,2)*isign(1,kf3)
      if(mstu(19).eq.1) then
        mstu(19)=0
      elseif(kq1.eq.0.and.kq2.eq.0.and.kq3.eq.0) then
      elseif(kq1.ne.0.and.kq2.eq.2.and.(kq1+kq3.eq.0.or.
     &  kq1+kq3.eq.4)) then
      else
        call pjerrm(2,'(PY3ENT:) unphysical flavour combination')
      endif
      k(ipa,2)=kf1
      k(ipa+1,2)=kf2
      k(ipa+2,2)=kf3
 
C...Store partons/particles in K vectors for normal case.
      if(ip.ge.0) then
        k(ipa,1)=1
        if(kq1.ne.0.and.(kq2.ne.0.or.kq3.ne.0)) k(ipa,1)=2
        k(ipa+1,1)=1
        if(kq2.ne.0.and.kq3.ne.0) k(ipa+1,1)=2
        k(ipa+2,1)=1
 
C...Store partons in K vectors for parton shower evolution.
      else
        k(ipa,1)=3
        k(ipa+1,1)=3
        k(ipa+2,1)=3
        kcs=4
        if(kq1.eq.-1) kcs=5
        k(ipa,kcs)=mstu(5)*(ipa+1)
        k(ipa,9-kcs)=mstu(5)*(ipa+2)
        k(ipa+1,kcs)=mstu(5)*(ipa+2)
        k(ipa+1,9-kcs)=mstu(5)*ipa
        k(ipa+2,kcs)=mstu(5)*ipa
        k(ipa+2,9-kcs)=mstu(5)*(ipa+1)
      endif
 
C...Check kinematics.
      mkerr=0
      if(0.5d0*x1*pecm.le.pm1.or.0.5d0*(2d0-x1-x3)*pecm.le.pm2.or.
     &0.5d0*x3*pecm.le.pm3) mkerr=1
      pa1=sqrt(max(1d-10,(0.5d0*x1*pecm)**2-pm1**2))
      pa2=sqrt(max(1d-10,(0.5d0*(2d0-x1-x3)*pecm)**2-pm2**2))
      pa3=sqrt(max(1d-10,(0.5d0*x3*pecm)**2-pm3**2))
      cthe2=(pa3**2-pa1**2-pa2**2)/(2d0*pa1*pa2)
      cthe3=(pa2**2-pa1**2-pa3**2)/(2d0*pa1*pa3)
      if(abs(cthe2).ge.1.001d0.or.abs(cthe3).ge.1.001d0) mkerr=1
      cthe3=max(-1d0,min(1d0,cthe3))
      if(mkerr.ne.0) call pjerrm(13,
     &'(PY3ENT:) unphysical kinematical variable setup')
 
C...Store partons/particles in P vectors.
      p(ipa,3)=pa1
      p(ipa,4)=sqrt(pa1**2+pm1**2)
      p(ipa,5)=pm1
      p(ipa+2,1)=pa3*sqrt(1d0-cthe3**2)
      p(ipa+2,3)=pa3*cthe3
      p(ipa+2,4)=sqrt(pa3**2+pm3**2)
      p(ipa+2,5)=pm3
      p(ipa+1,1)=-p(ipa+2,1)
      p(ipa+1,3)=-p(ipa,3)-p(ipa+2,3)
      p(ipa+1,4)=sqrt(p(ipa+1,1)**2+p(ipa+1,3)**2+pm2**2)
      p(ipa+1,5)=pm2
 
C...Set N. Optionally fragment/decay.
      n=ipa+2
      if(ip.eq.0) call pjexec
 
      return
      end
 
C*********************************************************************
 
C...PY4ENT
C...Stores four partons or particles in their CM frame, with
C...the first along the +z axis, the last in the xz plane with x > 0
C...and the second having y < 0 and y > 0 with equal probability.
 
      subroutine pj4ent(ip,kf1,kf2,kf3,kf4,pecm,x1,x2,x4,x12,x14)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jyjets/,/jydat1/,/jydat2/
 
C...Standard checks.
      mstu(28)=0
      if(mstu(12).ge.1) call pjlist(0)
      ipa=max(1,iabs(ip))
      if(ipa.gt.mstu(4)-3) call pjerrm(21,
     &'(PY4ENT:) writing outside PYJETS momory')
      kc1=jamcomp(kf1)
      kc2=jamcomp(kf2)
      kc3=jamcomp(kf3)
      kc4=jamcomp(kf4)
      if(kc1.eq.0.or.kc2.eq.0.or.kc3.eq.0.or.kc4.eq.0) call pjerrm(12,
     &'(PY4ENT:) unknown flavour code')
 
C...Find masses. Reset K, P and V vectors.
      pm1=0d0
      if(mstu(10).eq.1) pm1=p(ipa,5)
      if(mstu(10).ge.2) pm1=pjmass(kf1)
      pm2=0d0
      if(mstu(10).eq.1) pm2=p(ipa+1,5)
      if(mstu(10).ge.2) pm2=pjmass(kf2)
      pm3=0d0
      if(mstu(10).eq.1) pm3=p(ipa+2,5)
      if(mstu(10).ge.2) pm3=pjmass(kf3)
      pm4=0d0
      if(mstu(10).eq.1) pm4=p(ipa+3,5)
      if(mstu(10).ge.2) pm4=pjmass(kf4)
      do 110 i=ipa,ipa+3
        do 100 j=1,5
          k(i,j)=0
          p(i,j)=0d0
          v(i,j)=0d0
  100   continue
  110 continue
 
C...Check flavours.
      kq1=kchg(kc1,2)*isign(1,kf1)
      kq2=kchg(kc2,2)*isign(1,kf2)
      kq3=kchg(kc3,2)*isign(1,kf3)
      kq4=kchg(kc4,2)*isign(1,kf4)
      if(mstu(19).eq.1) then
        mstu(19)=0
      elseif(kq1.eq.0.and.kq2.eq.0.and.kq3.eq.0.and.kq4.eq.0) then
      elseif(kq1.ne.0.and.kq2.eq.2.and.kq3.eq.2.and.(kq1+kq4.eq.0.or.
     &  kq1+kq4.eq.4)) then
      elseif(kq1.ne.0.and.kq1+kq2.eq.0.and.kq3.ne.0.and.kq3+kq4.eq.0d0)
     &  then
      else
        call pjerrm(2,'(PY4ENT:) unphysical flavour combination')
      endif
      k(ipa,2)=kf1
      k(ipa+1,2)=kf2
      k(ipa+2,2)=kf3
      k(ipa+3,2)=kf4
 
C...Store partons/particles in K vectors for normal case.
      if(ip.ge.0) then
        k(ipa,1)=1
        if(kq1.ne.0.and.(kq2.ne.0.or.kq3.ne.0.or.kq4.ne.0)) k(ipa,1)=2
        k(ipa+1,1)=1
        if(kq2.ne.0.and.kq1+kq2.ne.0.and.(kq3.ne.0.or.kq4.ne.0))
     &  k(ipa+1,1)=2
        k(ipa+2,1)=1
        if(kq3.ne.0.and.kq4.ne.0) k(ipa+2,1)=2
        k(ipa+3,1)=1
 
C...Store partons for parton shower evolution from q-g-g-qbar or
C...g-g-g-g event.
      elseif(kq1+kq2.ne.0) then
        k(ipa,1)=3
        k(ipa+1,1)=3
        k(ipa+2,1)=3
        k(ipa+3,1)=3
        kcs=4
        if(kq1.eq.-1) kcs=5
        k(ipa,kcs)=mstu(5)*(ipa+1)
        k(ipa,9-kcs)=mstu(5)*(ipa+3)
        k(ipa+1,kcs)=mstu(5)*(ipa+2)
        k(ipa+1,9-kcs)=mstu(5)*ipa
        k(ipa+2,kcs)=mstu(5)*(ipa+3)
        k(ipa+2,9-kcs)=mstu(5)*(ipa+1)
        k(ipa+3,kcs)=mstu(5)*ipa
        k(ipa+3,9-kcs)=mstu(5)*(ipa+2)
 
C...Store partons for parton shower evolution from q-qbar-q-qbar event.
      else
        k(ipa,1)=3
        k(ipa+1,1)=3
        k(ipa+2,1)=3
        k(ipa+3,1)=3
        k(ipa,4)=mstu(5)*(ipa+1)
        k(ipa,5)=k(ipa,4)
        k(ipa+1,4)=mstu(5)*ipa
        k(ipa+1,5)=k(ipa+1,4)
        k(ipa+2,4)=mstu(5)*(ipa+3)
        k(ipa+2,5)=k(ipa+2,4)
        k(ipa+3,4)=mstu(5)*(ipa+2)
        k(ipa+3,5)=k(ipa+3,4)
      endif
 
C...Check kinematics.
      mkerr=0
      if(0.5d0*x1*pecm.le.pm1.or.0.5d0*x2*pecm.le.pm2.or.
     &0.5d0*(2d0-x1-x2-x4)*pecm.le.pm3.or.0.5d0*x4*pecm.le.pm4)
     &mkerr=1
      pa1=sqrt(max(1d-10,(0.5d0*x1*pecm)**2-pm1**2))
      pa2=sqrt(max(1d-10,(0.5d0*x2*pecm)**2-pm2**2))
      pa4=sqrt(max(1d-10,(0.5d0*x4*pecm)**2-pm4**2))
      x24=x1+x2+x4-1d0-x12-x14+(pm3**2-pm1**2-pm2**2-pm4**2)/pecm**2
      cthe4=(x1*x4-2d0*x14)*pecm**2/(4d0*pa1*pa4)
      if(abs(cthe4).ge.1.002d0) mkerr=1
      cthe4=max(-1d0,min(1d0,cthe4))
      sthe4=sqrt(1d0-cthe4**2)
      cthe2=(x1*x2-2d0*x12)*pecm**2/(4d0*pa1*pa2)
      if(abs(cthe2).ge.1.002d0) mkerr=1
      cthe2=max(-1d0,min(1d0,cthe2))
      sthe2=sqrt(1d0-cthe2**2)
      cphi2=((x2*x4-2d0*x24)*pecm**2-4d0*pa2*cthe2*pa4*cthe4)/
     &max(1d-8*pecm**2,4d0*pa2*sthe2*pa4*sthe4)
      if(abs(cphi2).ge.1.05d0) mkerr=1
      cphi2=max(-1d0,min(1d0,cphi2))
      if(mkerr.eq.1) call pjerrm(13,
     &'(PY4ENT:) unphysical kinematical variable setup')
 
C...Store partons/particles in P vectors.
      p(ipa,3)=pa1
      p(ipa,4)=sqrt(pa1**2+pm1**2)
      p(ipa,5)=pm1
      p(ipa+3,1)=pa4*sthe4
      p(ipa+3,3)=pa4*cthe4
      p(ipa+3,4)=sqrt(pa4**2+pm4**2)
      p(ipa+3,5)=pm4
      p(ipa+1,1)=pa2*sthe2*cphi2
      p(ipa+1,2)=pa2*sthe2*sqrt(1d0-cphi2**2)*(-1d0)**int(pjr(0)+0.5d0)
      p(ipa+1,3)=pa2*cthe2
      p(ipa+1,4)=sqrt(pa2**2+pm2**2)
      p(ipa+1,5)=pm2
      p(ipa+2,1)=-p(ipa+1,1)-p(ipa+3,1)
      p(ipa+2,2)=-p(ipa+1,2)
      p(ipa+2,3)=-p(ipa,3)-p(ipa+1,3)-p(ipa+3,3)
      p(ipa+2,4)=sqrt(p(ipa+2,1)**2+p(ipa+2,2)**2+p(ipa+2,3)**2+pm3**2)
      p(ipa+2,5)=pm3
 
C...Set N. Optionally fragment/decay.
      n=ipa+3
      if(ip.eq.0) call pjexec
 
      return
      end
 
C*********************************************************************
 
C...PYJOIN
C...Connects a sequence of partons with colour flow indices,
C...as required for subsequent shower evolution (or other operations).
 
      subroutine pjjoin(njoin,ijoin)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jyjets/,/jydat1/,/jydat2/
C...Local array.
      dimension ijoin(*)
 
C...Check that partons are of right types to be connected.
      if(njoin.lt.2) goto 120
      kqsum=0
      do 100 ijn=1,njoin
        i=ijoin(ijn)
        if(i.le.0.or.i.gt.n) goto 120
        if(k(i,1).lt.1.or.k(i,1).gt.3) goto 120
        kc=jamcomp(k(i,2))
        if(kc.eq.0) goto 120
        kq=kchg(kc,2)*isign(1,k(i,2))
        if(kq.eq.0) goto 120
        if(ijn.ne.1.and.ijn.ne.njoin.and.kq.ne.2) goto 120
        if(kq.ne.2) kqsum=kqsum+kq
        if(ijn.eq.1) kqs=kq
  100 continue
      if(kqsum.ne.0) goto 120
 
C...Connect the partons sequentially (closing for gluon loop).
      kcs=(9-kqs)/2
      if(kqs.eq.2) kcs=int(4.5d0+pjr(0))
      do 110 ijn=1,njoin
        i=ijoin(ijn)
        k(i,1)=3
        if(ijn.ne.1) ip=ijoin(ijn-1)
        if(ijn.eq.1) ip=ijoin(njoin)
        if(ijn.ne.njoin) in=ijoin(ijn+1)
        if(ijn.eq.njoin) in=ijoin(1)
        k(i,kcs)=mstu(5)*in
        k(i,9-kcs)=mstu(5)*ip
        if(ijn.eq.1.and.kqs.ne.2) k(i,9-kcs)=0
        if(ijn.eq.njoin.and.kqs.ne.2) k(i,kcs)=0
  110 continue
 
C...Error exit: no action taken.
      return
  120 call pjerrm(12,
     &'(PYJOIN:) given entries can not be joined by one string')
 
      return
      end
 
C*********************************************************************
 
C...PYGIVE
C...Sets values of commonblock variables.
 
      subroutine pjgive(chin)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/jydat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
      common/jydat4/chaf(500,2)
      character chaf*16
      common/jydatr/mrpy(6),rrpy(100)
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
      common/pjint8/xpvmd(-6:6),xpanl(-6:6),xpanh(-6:6),xpbeh(-6:6),
     &xpdir(-6:6)
      common/pjmssm/imss(0:99),rmss(0:99)
      save /jyjets/,/jydat1/,/jydat2/,/jydat3/,/jydat4/,/jydatr/,
     &/pjsubs/,/pjpars/,/pjint1/,/pjint2/,/pjint3/,/pjint4/,
     &/pjint5/,/pjint6/,/pjint7/,/pjint8/,/pjmssm/
C...Local arrays and character variables.
      character chin*(*),chfix*104,chbit*104,chold*8,chnew*8,chold2*28,
     &chnew2*28,chnam*6,chvar(49)*6,chalp(2)*26,chind*8,chini*10,
     &chinr*16
      dimension msvar(49,8)
 
C...For each variable to be translated give: name,
C...integer/real/character, no. of indices, lower&upper index bounds.
      data chvar/'N','K','P','V','MSTU','PARU','MSTJ','PARJ','KCHG',
     &'PMAS','PARF','VCKM','MDCY','MDME','BRAT','KFDP','CHAF','MRPY',
     &'RRPY','MSEL','MSUB','KFIN','CKIN','MSTP','PARP','MSTI','PARI',
     &'MINT','VINT','ISET','KFPR','COEF','ICOL','XSFX','ISIG','SIGH',
     &'MWID','WIDS','NGEN','XSEC','PROC','SIGT','XPVMD','XPANL',
     &'XPANH','XPBEH','XPDIR','IMSS','RMSS'/
      data ((msvar(i,j),j=1,8),i=1,49)/ 1,7*0,  1,2,1,4000,1,5,2*0,
     &2,2,1,4000,1,5,2*0,  2,2,1,4000,1,5,2*0,  1,1,1,200,4*0,
     &2,1,1,200,4*0,  1,1,1,200,4*0,  2,1,1,200,4*0,
     &1,2,1,500,1,4,2*0,  2,2,1,500,1,4,2*0,  2,1,1,2000,4*0,
     &2,2,1,4,1,4,2*0,  1,2,1,500,1,3,2*0,  1,2,1,4000,1,2,2*0,
     &2,1,1,4000,4*0,  1,2,1,4000,1,5,2*0,  3,2,1,500,1,2,2*0,
     &1,1,1,6,4*0,  2,1,1,100,4*0,
     &1,7*0,  1,1,1,500,4*0,  1,2,1,2,-40,40,2*0,  2,1,1,200,4*0,
     &1,1,1,200,4*0,  2,1,1,200,4*0,  1,1,1,200,4*0,  2,1,1,200,4*0,
     &1,1,1,400,4*0,  2,1,1,400,4*0,  1,1,1,500,4*0,
     &1,2,1,500,1,2,2*0,  2,2,1,500,1,20,2*0,  1,3,1,40,1,4,1,2,
     &2,2,1,2,-40,40,2*0,  1,2,1,1000,1,3,2*0,  2,1,1,1000,4*0,
     &1,1,1,500,4*0,   2,2,1,500,1,5,2*0,   1,2,0,500,1,3,2*0,
     &2,2,0,500,1,3,2*0,   4,1,0,500,4*0,   2,3,0,6,0,6,0,5,
     &2,1,-6,6,4*0,     2,1,-6,6,4*0,    2,1,-6,6,4*0,
     &2,1,-6,6,4*0,  2,1,-6,6,4*0,  1,1,0,99,4*0,  2,1,0,99,4*0/
      data chalp/'abcdefghijklmnopqrstuvwxyz',
     &'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
 
C...Length of character variable. Subdivide it into instructions.
      if(mstu(12).ge.1) call pjlist(0)
      chbit=chin//' '
      lbit=101
  100 lbit=lbit-1
      if(chbit(lbit:lbit).eq.' ') goto 100
      ltot=0
      do 110 lcom=1,lbit
        if(chbit(lcom:lcom).eq.' ') goto 110
        ltot=ltot+1
        chfix(ltot:ltot)=chbit(lcom:lcom)
  110 continue
      llow=0
  120 lhig=llow+1
  130 lhig=lhig+1
      if(lhig.le.ltot.and.chfix(lhig:lhig).ne.';') goto 130
      lbit=lhig-llow-1
      chbit(1:lbit)=chfix(llow+1:lhig-1)
 
C...Identify commonblock variable.
      lnam=1
  140 lnam=lnam+1
      if(chbit(lnam:lnam).ne.'('.and.chbit(lnam:lnam).ne.'='.and.
     &lnam.le.6) goto 140
      chnam=chbit(1:lnam-1)//' '
      do 160 lcom=1,lnam-1
        do 150 lalp=1,26
          if(chnam(lcom:lcom).eq.chalp(1)(lalp:lalp)) chnam(lcom:lcom)=
     &    chalp(2)(lalp:lalp)
  150   continue
  160 continue
      ivar=0
      do 170 iv=1,49
        if(chnam.eq.chvar(iv)) ivar=iv
  170 continue
      if(ivar.eq.0) then
        call pjerrm(18,'(PYGIVE:) do not recognize variable '//chnam)
        llow=lhig
        if(llow.lt.ltot) goto 120
        return
      endif
 
C...Identify any indices.
      i1=0
      i2=0
      i3=0
      nindx=0
      if(chbit(lnam:lnam).eq.'(') then
        lind=lnam
  180   lind=lind+1
        if(chbit(lind:lind).ne.')'.and.chbit(lind:lind).ne.',') goto 180
        chind=' '
        if((chbit(lnam+1:lnam+1).eq.'C'.or.chbit(lnam+1:lnam+1).eq.'c')
     &  .and.(ivar.eq.9.or.ivar.eq.10.or.ivar.eq.13.or.ivar.eq.17))
     &  then
          chind(lnam-lind+11:8)=chbit(lnam+2:lind-1)
          read(chind,'(I8)') kf
          i1=jamcomp(kf)
        elseif(chbit(lnam+1:lnam+1).eq.'C'.or.chbit(lnam+1:lnam+1).eq.
     &    'c') then
          call pjerrm(18,'(PYGIVE:) not allowed to use C index for '//
     &    chnam)
          llow=lhig
          if(llow.lt.ltot) goto 120
          return
        else
          chind(lnam-lind+10:8)=chbit(lnam+1:lind-1)
          read(chind,'(I8)') i1
        endif
        lnam=lind
        if(chbit(lnam:lnam).eq.')') lnam=lnam+1
        nindx=1
      endif
      if(chbit(lnam:lnam).eq.',') then
        lind=lnam
  190   lind=lind+1
        if(chbit(lind:lind).ne.')'.and.chbit(lind:lind).ne.',') goto 190
        chind=' '
        chind(lnam-lind+10:8)=chbit(lnam+1:lind-1)
        read(chind,'(I8)') i2
        lnam=lind
        if(chbit(lnam:lnam).eq.')') lnam=lnam+1
        nindx=2
      endif
      if(chbit(lnam:lnam).eq.',') then
        lind=lnam
  200   lind=lind+1
        if(chbit(lind:lind).ne.')'.and.chbit(lind:lind).ne.',') goto 200
        chind=' '
        chind(lnam-lind+10:8)=chbit(lnam+1:lind-1)
        read(chind,'(I8)') i3
        lnam=lind+1
        nindx=3
      endif
 
C...Check that indices allowed.
      ierr=0
      if(nindx.ne.msvar(ivar,2)) ierr=1
      if(nindx.ge.1.and.(i1.lt.msvar(ivar,3).or.i1.gt.msvar(ivar,4)))
     &ierr=2
      if(nindx.ge.2.and.(i2.lt.msvar(ivar,5).or.i2.gt.msvar(ivar,6)))
     &ierr=3
      if(nindx.eq.3.and.(i3.lt.msvar(ivar,7).or.i3.gt.msvar(ivar,8)))
     &ierr=4
      if(chbit(lnam:lnam).ne.'=') ierr=5
      if(ierr.ge.1) then
        call pjerrm(18,'(PYGIVE:) unallowed indices for '//
     &  chbit(1:lnam-1))
        llow=lhig
        if(llow.lt.ltot) goto 120
        return
      endif
 
C...Save old value of variable.
      if(ivar.eq.1) then
        iold=n
      elseif(ivar.eq.2) then
        iold=k(i1,i2)
      elseif(ivar.eq.3) then
        rold=p(i1,i2)
      elseif(ivar.eq.4) then
        rold=v(i1,i2)
      elseif(ivar.eq.5) then
        iold=mstu(i1)
      elseif(ivar.eq.6) then
        rold=paru(i1)
      elseif(ivar.eq.7) then
        iold=mstj(i1)
      elseif(ivar.eq.8) then
        rold=parj(i1)
      elseif(ivar.eq.9) then
        iold=kchg(i1,i2)
      elseif(ivar.eq.10) then
        rold=pmas(i1,i2)
      elseif(ivar.eq.11) then
        rold=parf(i1)
      elseif(ivar.eq.12) then
        rold=vckm(i1,i2)
      elseif(ivar.eq.13) then
        iold=mdcy(i1,i2)
      elseif(ivar.eq.14) then
        iold=mdme(i1,i2)
      elseif(ivar.eq.15) then
        rold=brat(i1)
      elseif(ivar.eq.16) then
        iold=kfdp(i1,i2)
      elseif(ivar.eq.17) then
        chold=chaf(i1,i2)
      elseif(ivar.eq.18) then
        iold=mrpy(i1)
      elseif(ivar.eq.19) then
        rold=rrpy(i1)
      elseif(ivar.eq.20) then
        iold=msel
      elseif(ivar.eq.21) then
        iold=msub(i1)
      elseif(ivar.eq.22) then
        iold=kfin(i1,i2)
      elseif(ivar.eq.23) then
        rold=ckin(i1)
      elseif(ivar.eq.24) then
        iold=mstp(i1)
      elseif(ivar.eq.25) then
        rold=parp(i1)
      elseif(ivar.eq.26) then
        iold=msti(i1)
      elseif(ivar.eq.27) then
        rold=pari(i1)
      elseif(ivar.eq.28) then
        iold=mint(i1)
      elseif(ivar.eq.29) then
        rold=vint(i1)
      elseif(ivar.eq.30) then
        iold=iset(i1)
      elseif(ivar.eq.31) then
        iold=kfpr(i1,i2)
      elseif(ivar.eq.32) then
        rold=coef(i1,i2)
      elseif(ivar.eq.33) then
        iold=icol(i1,i2,i3)
      elseif(ivar.eq.34) then
        rold=xsfx(i1,i2)
      elseif(ivar.eq.35) then
        iold=isig(i1,i2)
      elseif(ivar.eq.36) then
        rold=sigh(i1)
      elseif(ivar.eq.37) then
        iold=mwid(i1)
      elseif(ivar.eq.38) then
        rold=wids(i1,i2)
      elseif(ivar.eq.39) then
        iold=ngen(i1,i2)
      elseif(ivar.eq.40) then
        rold=xsec(i1,i2)
      elseif(ivar.eq.41) then
        chold2=proc(i1)
      elseif(ivar.eq.42) then
        rold=sigt(i1,i2,i3)
      elseif(ivar.eq.43) then
        rold=xpvmd(i1)
      elseif(ivar.eq.44) then
        rold=xpanl(i1)
      elseif(ivar.eq.45) then
        rold=xpanh(i1)
      elseif(ivar.eq.46) then
        rold=xpbeh(i1)
      elseif(ivar.eq.47) then
        rold=xpdir(i1)
      elseif(ivar.eq.48) then
        iold=imss(i1)
      elseif(ivar.eq.49) then
        rold=rmss(i1)
      endif
 
C...Print current value of variable. Loop back.
      if(lnam.ge.lbit) then
        chbit(lnam:14)=' '
        chbit(15:60)=' has the value                                '
        if(msvar(ivar,1).eq.1) then
          write(chbit(51:60),'(I10)') iold
        elseif(msvar(ivar,1).eq.2) then
          write(chbit(47:60),'(F14.5)') rold
        elseif(msvar(ivar,1).eq.3) then
          chbit(53:60)=chold
        else
          chbit(33:60)=chold
        endif
        if(mstu(13).ge.1) write(mstu(11),5000) chbit(1:60)
        llow=lhig
        if(llow.lt.ltot) goto 120
        return
      endif
 
C...Read in new variable value.
      if(msvar(ivar,1).eq.1) then
        chini=' '
        chini(lnam-lbit+11:10)=chbit(lnam+1:lbit)
        read(chini,'(I10)') inew
      elseif(msvar(ivar,1).eq.2) then
        chinr=' '
        chinr(lnam-lbit+17:16)=chbit(lnam+1:lbit)
        read(chinr,*) rnew
      elseif(msvar(ivar,1).eq.3) then
        chnew=chbit(lnam+1:lbit)//' '
      else
        chnew2=chbit(lnam+1:lbit)//' '
      endif
 
C...Store new variable value.
      if(ivar.eq.1) then
        n=inew
      elseif(ivar.eq.2) then
        k(i1,i2)=inew
      elseif(ivar.eq.3) then
        p(i1,i2)=rnew
      elseif(ivar.eq.4) then
        v(i1,i2)=rnew
      elseif(ivar.eq.5) then
        mstu(i1)=inew
      elseif(ivar.eq.6) then
        paru(i1)=rnew
      elseif(ivar.eq.7) then
        mstj(i1)=inew
      elseif(ivar.eq.8) then
        parj(i1)=rnew
      elseif(ivar.eq.9) then
        kchg(i1,i2)=inew
      elseif(ivar.eq.10) then
        pmas(i1,i2)=rnew
      elseif(ivar.eq.11) then
        parf(i1)=rnew
      elseif(ivar.eq.12) then
        vckm(i1,i2)=rnew
      elseif(ivar.eq.13) then
        mdcy(i1,i2)=inew
      elseif(ivar.eq.14) then
        mdme(i1,i2)=inew
      elseif(ivar.eq.15) then
        brat(i1)=rnew
      elseif(ivar.eq.16) then
        kfdp(i1,i2)=inew
      elseif(ivar.eq.17) then
        chaf(i1,i2)=chnew
      elseif(ivar.eq.18) then
        mrpy(i1)=inew
      elseif(ivar.eq.19) then
        rrpy(i1)=rnew
      elseif(ivar.eq.20) then
        msel=inew
      elseif(ivar.eq.21) then
        msub(i1)=inew
      elseif(ivar.eq.22) then
        kfin(i1,i2)=inew
      elseif(ivar.eq.23) then
        ckin(i1)=rnew
      elseif(ivar.eq.24) then
        mstp(i1)=inew
      elseif(ivar.eq.25) then
        parp(i1)=rnew
      elseif(ivar.eq.26) then
        msti(i1)=inew
      elseif(ivar.eq.27) then
        pari(i1)=rnew
      elseif(ivar.eq.28) then
        mint(i1)=inew
      elseif(ivar.eq.29) then
        vint(i1)=rnew
      elseif(ivar.eq.30) then
        iset(i1)=inew
      elseif(ivar.eq.31) then
        kfpr(i1,i2)=inew
      elseif(ivar.eq.32) then
        coef(i1,i2)=rnew
      elseif(ivar.eq.33) then
        icol(i1,i2,i3)=inew
      elseif(ivar.eq.34) then
        xsfx(i1,i2)=rnew
      elseif(ivar.eq.35) then
        isig(i1,i2)=inew
      elseif(ivar.eq.36) then
        sigh(i1)=rnew
      elseif(ivar.eq.37) then
        mwid(i1)=inew
      elseif(ivar.eq.38) then
        wids(i1,i2)=rnew
      elseif(ivar.eq.39) then
        ngen(i1,i2)=inew
      elseif(ivar.eq.40) then
        xsec(i1,i2)=rnew
      elseif(ivar.eq.41) then
        proc(i1)=chnew2
      elseif(ivar.eq.42) then
        sigt(i1,i2,i3)=rnew
      elseif(ivar.eq.43) then
        xpvmd(i1)=rnew
      elseif(ivar.eq.44) then
        xpanl(i1)=rnew
      elseif(ivar.eq.45) then
        xpanh(i1)=rnew
      elseif(ivar.eq.46) then
        xpbeh(i1)=rnew
      elseif(ivar.eq.47) then
        xpdir(i1)=rnew
      elseif(ivar.eq.48) then
        imss(i1)=inew
      elseif(ivar.eq.49) then
        rmss(i1)=rnew
      endif
 
C...Write old and new value. Loop back.
      chbit(lnam:14)=' '
      chbit(15:60)=' changed from                to               '
      if(msvar(ivar,1).eq.1) then
        write(chbit(33:42),'(I10)') iold
        write(chbit(51:60),'(I10)') inew
        if(mstu(13).ge.1) write(mstu(11),5000) chbit(1:60)
      elseif(msvar(ivar,1).eq.2) then
        write(chbit(29:42),'(F14.5)') rold
        write(chbit(47:60),'(F14.5)') rnew
        if(mstu(13).ge.1) write(mstu(11),5000) chbit(1:60)
      elseif(msvar(ivar,1).eq.3) then
        chbit(35:42)=chold
        chbit(53:60)=chnew
        if(mstu(13).ge.1) write(mstu(11),5000) chbit(1:60)
      else
        chbit(15:88)=' changed from '//chold2//' to '//chnew2
        if(mstu(13).ge.1) write(mstu(11),5100) chbit(1:88)
      endif
      llow=lhig
      if(llow.lt.ltot) goto 120
 
C...Format statement for output on unit MSTU(11) (by default 6).
 5000 format(5x,a60)
 5100 format(5x,a88)
 
      return
      end
 
C*********************************************************************
 
C...PYEXEC
C...Administrates the fragmentation and decay chain.
 
      subroutine pjexec
 
C...Double precision.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/jydat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
      common/pjint4/mwid(500),wids(500,5)
      save /jyjets/,/jydat1/,/jydat2/,/jydat3/,/pjint4/
C...Local array.
      dimension ps(2,6),ijoin(100)
 
C...Initialize and reset.
      mstu(24)=0
      if(mstu(12).ge.1) call pjlist(0)
      mstu(31)=mstu(31)+1
      mstu(1)=0
      mstu(2)=0
      mstu(3)=0
      if(mstu(17).le.0) mstu(90)=0
      mcons=1
 
C...Sum up momentum, energy and charge for starting entries.
      nsav=n
      do 110 i=1,2
        do 100 j=1,6
          ps(i,j)=0d0
  100   continue
  110 continue
      do 130 i=1,n
        if(k(i,1).le.0.or.k(i,1).gt.10) goto 130
        do 120 j=1,4
          ps(1,j)=ps(1,j)+p(i,j)
  120   continue
        ps(1,6)=ps(1,6)+jamchge(k(i,2))
  130 continue
      paru(21)=ps(1,4)
 
C...Prepare system for subsequent fragmentation/decay.
      call pjprep(0)
 
C...Loop through jet fragmentation and particle decays.
      mbe=0
  140 mbe=mbe+1
      ip=0
  150 ip=ip+1
      kc=0
      if(k(ip,1).gt.0.and.k(ip,1).le.10) kc=jamcomp(k(ip,2))
      if(kc.eq.0) then
 
C...Deal with any remaining undecayed resonance
C...(normally the task of PYEVNT, so seldom used).
      elseif(mwid(kc).ne.0) then
        ibeg=ip
        if(kchg(kc,2).ne.0.and.k(i,1).ne.3) then
          ibeg=ip+1
  160     ibeg=ibeg-1
          if(ibeg.ge.2.and.k(ibeg,1).eq.2) goto 160
          if(k(ibeg,1).ne.2) ibeg=ibeg+1
          iend=ip-1
  170     iend=iend+1
          if(iend.lt.n.and.k(iend,1).eq.2) goto 170
          if(iend.lt.n.and.kchg(jamcomp(k(iend,2)),2).eq.0) goto 170
          njoin=0
          do 180 i=ibeg,iend
            if(kchg(jamcomp(k(iend,2)),2).ne.0) then
              njoin=njoin+1
              ijoin(njoin)=i
            endif
  180     continue
        endif
        call pjresd(ip)
        call pjprep(ibeg)
 
C...Particle decay if unstable and allowed. Save long-lived particle
C...decays until second pass after Bose-Einstein effects.
      elseif(kchg(kc,2).eq.0) then
        if(mstj(21).ge.1.and.mdcy(kc,1).ge.1.and.(mstj(51).le.0.or.mbe
     &  .eq.2.or.pmas(kc,2).ge.parj(91).or.iabs(k(ip,2)).eq.311))
     &  call pjdecy(ip,icon)
 
C...Decay products may develop a shower.
        if(mstj(92).gt.0) then
          ip1=mstj(92)
          qmax=sqrt(max(0d0,(p(ip1,4)+p(ip1+1,4))**2-(p(ip1,1)+p(ip1+1,
     &    1))**2-(p(ip1,2)+p(ip1+1,2))**2-(p(ip1,3)+p(ip1+1,3))**2))
          call pjshow(ip1,ip1+1,qmax)
          call pjprep(ip1)
          mstj(92)=0
        elseif(mstj(92).lt.0) then
          ip1=-mstj(92)
          call pjshow(ip1,-3,p(ip,5))
          call pjprep(ip1)
          mstj(92)=0
        endif
 
C...Jet fragmentation: string or independent fragmentation.
      elseif(k(ip,1).eq.1.or.k(ip,1).eq.2) then
        mfrag=mstj(1)
        if(mfrag.ge.1.and.k(ip,1).eq.1) mfrag=2
        if(mstj(21).ge.2.and.k(ip,1).eq.2.and.n.gt.ip) then
          if(k(ip+1,1).eq.1.and.k(ip+1,3).eq.k(ip,3).and.
     &    k(ip,3).gt.0.and.k(ip,3).lt.ip) then
            if(kchg(jamcomp(k(k(ip,3),2)),2).eq.0) mfrag=min(1,mfrag)
          endif
        endif
        if(mfrag.eq.1) call pjstrf(ip)
        if(mfrag.eq.2) call pjindf(ip)
        if(mfrag.eq.2.and.k(ip,1).eq.1) mcons=0
        if(mfrag.eq.2.and.(mstj(3).le.0.or.mod(mstj(3),5).eq.0)) mcons=0
      endif
 
C...Loop back if enough space left in PYJETS and no error abort.
      if(mstu(24).ne.0.and.mstu(21).ge.2) then
      elseif(ip.lt.n.and.n.lt.mstu(4)-20-mstu(32)) then
        goto 150
      elseif(ip.lt.n) then
        call pjerrm(11,'(PYEXEC:) no more memory left in PYJETS')
      endif
 
C...Include simple Bose-Einstein effect parametrization if desired.
      if(mbe.eq.1.and.mstj(51).ge.1) then
        call pjboei(nsav)
        goto 140
      endif
 
C...Check that momentum, energy and charge were conserved.
      do 200 i=1,n
        if(k(i,1).le.0.or.k(i,1).gt.10) goto 200
        do 190 j=1,4
          ps(2,j)=ps(2,j)+p(i,j)
  190   continue
        ps(2,6)=ps(2,6)+jamchge(k(i,2))
  200 continue
      pdev=(abs(ps(2,1)-ps(1,1))+abs(ps(2,2)-ps(1,2))+abs(ps(2,3)-
     &ps(1,3))+abs(ps(2,4)-ps(1,4)))/(1d0+abs(ps(2,4))+abs(ps(1,4)))
      if(mcons.eq.1.and.pdev.gt.paru(11)) call pjerrm(15,
     &'(PYEXEC:) four-momentum was not conserved')
      if(mcons.eq.1.and.abs(ps(2,6)-ps(1,6)).gt.0.1d0) call pjerrm(15,
     &'(PYEXEC:) charge was not conserved')
 
      return
      end
 
C*********************************************************************
 
C...PYPREP
C...Rearranges partons along strings. Allows small systems
C...to collapse into one or two particles and checks flavours.
 
      subroutine pjprep(ip)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
C...Commonblocks.
      common/jampos1/jqconst(2),kfcq(4),icq(4),icms
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      save /jyjets/,/jampos1/
c     common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
c     common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
c     common/jydat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
c     save /jyjets/,/jydat1/,/jydat2/,/jydat3/

C...Local arrays.
      dimension dps(5),dpc(5),ue(3)
 
C...Rearrange parton shower product listing along strings: begin loop.
      i1=n
      do 130 mqgst=1,2
        do 120 i=max(1,ip),n
          if(k(i,1).ne.3) goto 120
          kc=jamcomp(k(i,2))
          if(kc.eq.0) goto 120
          kq=kchg(kc,2)
          if(kq.eq.0.or.(mqgst.eq.1.and.kq.eq.2)) goto 120
 
C...Pick up loose string end.
          kcs=4
          if(kq*isign(1,k(i,2)).lt.0) kcs=5
          ia=i
          nstp=0
  100     nstp=nstp+1
          if(nstp.gt.4*n) then
            call pjerrm(14,'(PYPREP:) caught in infinite loop')
            return
          endif
 
C...Copy undecayed parton.
          if(k(ia,1).eq.3) then
            if(i1.ge.mstu(4)-mstu(32)-5) then
              call pjerrm(11,'(PYPREP:) no more memory left in PYJETS')
              return
            endif
            i1=i1+1
            k(i1,1)=2
            if(nstp.ge.2.and.kchg(jamcomp(k(ia,2)),2).ne.2) k(i1,1)=1
            k(i1,2)=k(ia,2)
            k(i1,3)=ia
            k(i1,4)=0
            k(i1,5)=0
            do 110 j=1,5
              p(i1,j)=p(ia,j)
              v(i1,j)=v(ia,j)
  110       continue
            k(ia,1)=k(ia,1)+10
            if(k(i1,1).eq.1) goto 120
          endif
 
C...Go to next parton in colour space.
          ib=ia
          if(mod(k(ib,kcs)/mstu(5)**2,2).eq.0.and.mod(k(ib,kcs),mstu(5))
     &    .ne.0) then
            ia=mod(k(ib,kcs),mstu(5))
            k(ib,kcs)=k(ib,kcs)+mstu(5)**2
            mrev=0
          else
            if(k(ib,kcs).ge.2*mstu(5)**2.or.mod(k(ib,kcs)/mstu(5),
     &      mstu(5)).eq.0) kcs=9-kcs
            ia=mod(k(ib,kcs)/mstu(5),mstu(5))
            k(ib,kcs)=k(ib,kcs)+2*mstu(5)**2
            mrev=1
          endif
          if(ia.le.0.or.ia.gt.n) then
            call pjerrm(12,'(PYPREP:) colour rearrangement failed')
            return
          endif
          if(mod(k(ia,4)/mstu(5),mstu(5)).eq.ib.or.mod(k(ia,5)/mstu(5),
     &    mstu(5)).eq.ib) then
            if(mrev.eq.1) kcs=9-kcs
            if(mod(k(ia,kcs)/mstu(5),mstu(5)).ne.ib) kcs=9-kcs
            k(ia,kcs)=k(ia,kcs)+2*mstu(5)**2
          else
            if(mrev.eq.0) kcs=9-kcs
            if(mod(k(ia,kcs),mstu(5)).ne.ib) kcs=9-kcs
            k(ia,kcs)=k(ia,kcs)+mstu(5)**2
          endif
          if(ia.ne.i) goto 100
          k(i1,1)=1
  120   continue
  130 continue
      n=i1
      if(mstj(14).lt.0) return
 
C...Find lowest-mass colour singlet jet system, OK if above threshold.
      if(mstj(14).eq.0) goto 320
      ns=n
  140 nsin=n-ns
      pdm=1d0+parj(32)
      ic=0
      do 190 i=max(1,ip),ns
        if(k(i,1).ne.1.and.k(i,1).ne.2) then
        elseif(k(i,1).eq.2.and.ic.eq.0) then
          nsin=nsin+1
          ic=i
          do 150 j=1,4
            dps(j)=p(i,j)
  150     continue
          mstj(93)=1
          dps(5)=pjmass(k(i,2))
        elseif(k(i,1).eq.2) then
          do 160 j=1,4
            dps(j)=dps(j)+p(i,j)
  160     continue
        elseif(ic.ne.0.and.kchg(jamcomp(k(i,2)),2).ne.0) then
          do 170 j=1,4
            dps(j)=dps(j)+p(i,j)
  170     continue
          mstj(93)=1
          dps(5)=dps(5)+pjmass(k(i,2))
          pd=sqrt(max(0d0,dps(4)**2-dps(1)**2-dps(2)**2-dps(3)**2))-
     &    dps(5)
          if(pd.lt.pdm) then
            pdm=pd
            do 180 j=1,5
              dpc(j)=dps(j)
  180       continue
            ic1=ic
            ic2=i
          endif
          ic=0
        else
          nsin=nsin+1
        endif
  190 continue
      if(pdm.ge.parj(32)) goto 320
 
C...Fill small-mass system as cluster.
      nsav=n
      pecm=sqrt(max(0d0,dpc(4)**2-dpc(1)**2-dpc(2)**2-dpc(3)**2))
      k(n+1,1)=11
      k(n+1,2)=91
      k(n+1,3)=ic1
      k(n+1,4)=n+2
      k(n+1,5)=n+3
      p(n+1,1)=dpc(1)
      p(n+1,2)=dpc(2)
      p(n+1,3)=dpc(3)
      p(n+1,4)=dpc(4)
      p(n+1,5)=pecm
 
C...Form two particles from flavours of lowest-mass system, if feasible.
      k(n+2,1)=1
      k(n+3,1)=1
      if(mstu(16).ne.2) then
        k(n+2,3)=n+1
        k(n+3,3)=n+1
      else
        k(n+2,3)=ic1
        k(n+3,3)=ic2
      endif
      k(n+2,4)=0
      k(n+3,4)=0
      k(n+2,5)=0
      k(n+3,5)=0
      if(iabs(k(ic1,2)).ne.21) then
        kc1=jamcomp(k(ic1,2))
        kc2=jamcomp(k(ic2,2))
        if(kc1.eq.0.or.kc2.eq.0) then
          call pjerrm(2,'(PYPREP:) kc1=0 or kc2=0')
          goto 320
        endif
        kq1=kchg(kc1,2)*isign(1,k(ic1,2))
        kq2=kchg(kc2,2)*isign(1,k(ic2,2))
        if(kq1+kq2.ne.0)then
          call pjerrm(2,'(PYPREP:) kq1+kq2 non zero')
          goto 320
        endif
C.. Start with qq, if there is one. Only allow for rank 1 popcorn meson
  200   k1=k(ic1,2)
        if(iabs(k(ic2,2)).gt.10) k1=k(ic2,2)
        mstu(125)=0
        call pjdcyk(k1,0,kfln,k(n+2,2))
        call pjdcyk(k(ic1,2)+k(ic2,2)-k1,-kfln,kfldmp,k(n+3,2))
        if(k(n+2,2).eq.0.or.k(n+3,2).eq.0) goto 200
      else
        if(iabs(k(ic2,2)).ne.21) then
          call pjerrm(2,'(PYPREP:) gg system??')
          goto 320
        endif
C.. No room for popcorn mesons in closed string -> 2 hadrons.
        mstu(125)=0
  210   call pjdcyk(1+int((2d0+parj(2))*pjr(0)),0,kfln,kfdmp)
        call pjdcyk(kfln,0,kflm,k(n+2,2))
        call pjdcyk(-kfln,-kflm,kfldmp,k(n+3,2))
        if(k(n+2,2).eq.0.or.k(n+3,2).eq.0) goto 210
      endif
      p(n+2,5)=pjmass(k(n+2,2))
      p(n+3,5)=pjmass(k(n+3,2))
      if(p(n+2,5)+p(n+3,5)+parj(64).ge.pecm.and.nsin.eq.1) then
        call pjerrm(2,
     $ '(PYPREP:) two particle decay impossible due to maas')
        goto 320
      endif
      if(p(n+2,5)+p(n+3,5)+parj(64).ge.pecm) goto 260
 
C...Perform two-particle decay of jet system, if possible.
      if(pecm.ge.0.02d0*dpc(4)) then
        pa=sqrt((pecm**2-(p(n+2,5)+p(n+3,5))**2)*(pecm**2-
     &  (p(n+2,5)-p(n+3,5))**2))/(2d0*pecm)
        ue(3)=2d0*pjr(0)-1d0
        phi=paru(2)*pjr(0)
        ue(1)=sqrt(1d0-ue(3)**2)*cos(phi)
        ue(2)=sqrt(1d0-ue(3)**2)*sin(phi)
        do 220 j=1,3
          p(n+2,j)=pa*ue(j)
          p(n+3,j)=-pa*ue(j)
  220   continue
        p(n+2,4)=sqrt(pa**2+p(n+2,5)**2)
        p(n+3,4)=sqrt(pa**2+p(n+3,5)**2)
        mstu(33)=1
        call pjrobo(n+2,n+3,0d0,0d0,dpc(1)/dpc(4),dpc(2)/dpc(4),
     &  dpc(3)/dpc(4))
      else
        np=0
        do 230 i=ic1,ic2
          if(k(i,1).eq.1.or.k(i,1).eq.2) np=np+1
  230   continue
        ha=p(ic1,4)*p(ic2,4)-p(ic1,1)*p(ic2,1)-p(ic1,2)*p(ic2,2)-
     &  p(ic1,3)*p(ic2,3)
        if(np.ge.3.or.ha.le.1.25d0*p(ic1,5)*p(ic2,5)) goto 260
        hd1=0.5d0*(p(n+2,5)**2-p(ic1,5)**2)
        hd2=0.5d0*(p(n+3,5)**2-p(ic2,5)**2)
        hr=sqrt(max(0d0,((ha-hd1-hd2)**2-(p(n+2,5)*p(n+3,5))**2)/
     &  (ha**2-(p(ic1,5)*p(ic2,5))**2)))-1d0
        hc=p(ic1,5)**2+2d0*ha+p(ic2,5)**2
        hk1=((p(ic2,5)**2+ha)*hr+hd1-hd2)/hc
        hk2=((p(ic1,5)**2+ha)*hr+hd2-hd1)/hc
        do 240 j=1,4
          p(n+2,j)=(1d0+hk1)*p(ic1,j)-hk2*p(ic2,j)
          p(n+3,j)=(1d0+hk2)*p(ic2,j)-hk1*p(ic1,j)
  240   continue
      endif
      do 250 j=1,4
        v(n+1,j)=v(ic1,j)
        v(n+2,j)=v(ic1,j)
        v(n+3,j)=v(ic2,j)
  250 continue
      v(n+1,5)=0d0
      v(n+2,5)=0d0
      v(n+3,5)=0d0
      n=n+3
c+JAM
      k(n+1,4)=jqconst(1)
      k(n+2,4)=jqconst(2)
c-JAM
      goto 300
 
C...Else form one particle from the flavours available, if possible.
  260 k(n+1,5)=n+2
      if(iabs(k(ic1,2)).gt.100.and.iabs(k(ic2,2)).gt.100) then
        goto 320
      elseif(iabs(k(ic1,2)).ne.21) then
        call pjkfdi(k(ic1,2),k(ic2,2),kfldmp,k(n+2,2))
      else
        kfln=1+int((2d0+parj(2))*pjr(0))
        call pjkfdi(kfln,-kfln,kfldmp,k(n+2,2))
      endif
      if(k(n+2,2).eq.0) goto 260
      p(n+2,5)=pjmass(k(n+2,2))
 
C...Find parton/particle which combines to largest extra mass.
      ir=0
      ha=0d0
      hsm=0d0
      do 280 mcomb=1,3
        if(ir.ne.0) goto 280
        do 270 i=max(1,ip),n
          if(k(i,1).le.0.or.k(i,1).gt.10.or.(i.ge.ic1.and.i.le.ic2
     &    .and.k(i,1).ge.1.and.k(i,1).le.2)) goto 270
          if(mcomb.eq.1) kci=jamcomp(k(i,2))
          if(mcomb.eq.1.and.kci.eq.0) goto 270
          if(mcomb.eq.1.and.kchg(kci,2).eq.0.and.i.le.ns) goto 270
          if(mcomb.eq.2.and.iabs(k(i,2)).gt.10.and.iabs(k(i,2)).le.100)
     &    goto 270
          hcr=dpc(4)*p(i,4)-dpc(1)*p(i,1)-dpc(2)*p(i,2)-dpc(3)*p(i,3)
          hsr=2d0*hcr+pecm**2-p(n+2,5)**2-2d0*p(n+2,5)*p(i,5)
          if(hsr.gt.hsm) then
            ir=i
            ha=hcr
            hsm=hsr
          endif
  270   continue
  280 continue
 
C...Shuffle energy and momentum to put new particle on mass shell.
      if(ir.ne.0) then
        hb=pecm**2+ha
        hc=p(n+2,5)**2+ha
        hd=p(ir,5)**2+ha
        hk2=0.5d0*(hb*sqrt(max(0d0,((hb+hc)**2-4d0*(hb+hd)*p(n+2,5)**2)/
     &  (ha**2-(pecm*p(ir,5))**2)))-(hb+hc))/(hb+hd)
        hk1=(0.5d0*(p(n+2,5)**2-pecm**2)+hd*hk2)/hb
        do 290 j=1,4
          p(n+2,j)=(1d0+hk1)*dpc(j)-hk2*p(ir,j)
          p(ir,j)=(1d0+hk2)*p(ir,j)-hk1*dpc(j)
          v(n+1,j)=v(ic1,j)
          v(n+2,j)=v(ic1,j)
  290   continue
        v(n+1,5)=0d0
        v(n+2,5)=0d0
        n=n+2
      else
        call pjerrm(3,'(PYPREP:) no match for collapsing cluster')
        return
      endif
 
C...Mark collapsed system and store daughter pointers. Iterate.
  300 do 310 i=ic1,ic2
        if((k(i,1).eq.1.or.k(i,1).eq.2)
     $               .and.kchg(jamcomp(k(i,2)),2).ne.0) then
          k(i,1)=k(i,1)+10
          if(mstu(16).ne.2) then
            k(i,4)=nsav+1
            k(i,5)=nsav+1
          else
            k(i,4)=nsav+2
            k(i,5)=n
          endif
        endif
  310 continue
      if(n.lt.mstu(4)-mstu(32)-5) goto 140
 
C...Check flavours and invariant masses in parton systems.
  320 np=0
      kfn=0
      kqs=0
      do 330 j=1,5
        dps(j)=0d0
  330 continue
      do 360 i=max(1,ip),n
        if(k(i,1).le.0.or.k(i,1).gt.10) goto 360
        kc=jamcomp(k(i,2))
        if(kc.eq.0) goto 360
        kq=kchg(kc,2)*isign(1,k(i,2))
        if(kq.eq.0) goto 360
        np=np+1
        if(kq.ne.2) then
          kfn=kfn+1
          kqs=kqs+kq
          mstj(93)=1
          dps(5)=dps(5)+pjmass(k(i,2))
        endif
        do 340 j=1,4
          dps(j)=dps(j)+p(i,j)
  340   continue

        if(k(i,1).eq.1) then
          if(np.ne.1.and.(kfn.eq.1.or.kfn.ge.3.or.kqs.ne.0)) call
     &    pjerrm(2,'(PYPREP:) unphysical flavour combination')

          if(np.ne.1.and.dps(4)**2-dps(1)**2-dps(2)**2-dps(3)**2.lt.
     &    (0.9d0*parj(32)+dps(5))**2) then
            write(check(1),'(''n kf1 kf2='',i4,1x,i9,1x,i9)')
     $                                             n,k(2,1),k(2,n)
            write(check(2),'(''mass='',g13.4)')
     $                sqrt(dps(4)**2-dps(1)**2-dps(2)**2-dps(3)**2)
            call jamerrm(3,2,'(PYPREP:) too small mass in jet system')
          endif

          np=0
          kfn=0
          kqs=0
          do 350 j=1,5
            dps(j)=0d0
  350     continue
        endif

  360 continue
 
      return
      end
 
C*********************************************************************
 
C...PYSTRF
C...Handles the fragmentation of an arbitrary colour singlet
C...jet system according to the Lund string fragmentation model.
 
      subroutine pjstrf(ip)
 
C...Double precision.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jyjets/,/jydat1/,/jydat2/
c...JAM:
      common/jampos1/jqconst(2),kfcq(4),icq(4),icms
      save /jampos1/
C...Local arrays. All MOPS variables ends with MO
      dimension dps(5),kfl(3),pmq(3),px(3),py(3),gam(3),ie(2),pr(2),
     &in(9),dhm(4),dhg(4),dp(5,5),irank(2),mju(4),iju(3),pju(5,5),
     &tju(5),kfjh(2),njs(2),kfjs(2),pjs(4,5),mstu9t(8),paru9t(8),
     &inmo(9),pm2qmo(2),xtmo(2)
c....JAM:Local arrays.
      dimension zposp(2),zposm(2),iqconst(2)
 
C...Function: four-product of two vectors.
      four(i,j)=p(i,4)*p(j,4)-p(i,1)*p(j,1)-p(i,2)*p(j,2)-p(i,3)*p(j,3)
      dfour(i,j)=dp(i,4)*dp(j,4)-dp(i,1)*dp(j,1)-dp(i,2)*dp(j,2)-
     &dp(i,3)*dp(j,3)
 
C...Reset counters. Identify parton system.
      mstj(91)=0
      nsav=n
      mstu90=mstu(90)
      np=0
      kqsum=0
      do 100 j=1,5
        dps(j)=0d0
  100 continue
      mju(1)=0
      mju(2)=0

      i=ip-1
  110 i=i+1
      if(i.gt.min(n,mstu(4)-mstu(32))) then
        call pjerrm(12,'(PYSTRF:) failed to reconstruct jet system')
        if(mstu(21).ge.1) return
      endif
      if(k(i,1).ne.1.and.k(i,1).ne.2.and.k(i,1).ne.41) goto 110
      kc=jamcomp(k(i,2))
      if(kc.eq.0) goto 110
      kq=kchg(kc,2)*isign(1,k(i,2))
      if(kq.eq.0) goto 110
      if(n+5*np+11.gt.mstu(4)-mstu(32)-5) then
        call pjerrm(11,'(PYSTRF:) no more memory left in PYJETS')
        if(mstu(21).ge.1) return
      endif
 
C...Take copy of partons to be considered. Check flavour sum.
      np=np+1
      do 120 j=1,5
        k(n+np,j)=k(i,j)
        p(n+np,j)=p(i,j)
        if(j.ne.4) dps(j)=dps(j)+p(i,j)
  120 continue
      dps(4)=dps(4)+sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
      k(n+np,3)=i
      if(kq.ne.2) kqsum=kqsum+kq
      if(k(i,1).eq.41) then
        kqsum=kqsum+2*kq
        if(kqsum.eq.kq) mju(1)=n+np
        if(kqsum.ne.kq) mju(2)=n+np
      endif
      if(k(i,1).eq.2.or.k(i,1).eq.41) goto 110

      if(kqsum.ne.0) then
        call pjerrm(12,'(PYSTRF:) unphysical flavour combination')
        if(mstu(21).ge.1) return
      endif

cJAM++
      wstr=sqrt(dps(4)**2-(dps(1)**2+dps(2)**2+dps(3)**2))
      mstj12=mstj(12)
c     if(wstr.le.2.5d0) mstj(12)=1
cJAM--
 
C...Boost copied system to CM frame (for better numerical precision).
cJAM++
      if(icms.eq.0) then
      if(abs(dps(3)).lt.0.99d0*dps(4)) then
        mbst=0
        mstu(33)=1
        call pjrobo(n+1,n+np,0d0,0d0,-dps(1)/dps(4),-dps(2)/dps(4),
     &  -dps(3)/dps(4))
      else
        mbst=1
        hhbz=sqrt(max(1d-6,dps(4)+dps(3))/max(1d-6,dps(4)-dps(3)))
        do 130 i=n+1,n+np
          hhpmt=p(i,1)**2+p(i,2)**2+p(i,5)**2
          if(p(i,3).gt.0d0) then
            hhpez=(p(i,4)+p(i,3))/hhbz
            p(i,3)=0.5d0*(hhpez-hhpmt/hhpez)
            p(i,4)=0.5d0*(hhpez+hhpmt/hhpez)
          else
            hhpez=(p(i,4)-p(i,3))*hhbz
            p(i,3)=-0.5d0*(hhpez-hhpmt/hhpez)
            p(i,4)=0.5d0*(hhpez+hhpmt/hhpez)
          endif
  130   continue
      endif
      endif
cJAM--
 
C...Search for very nearby partons that may be recombined.
      ntryr=0
      paru12=paru(12)  ! 0.09GeV^2
      paru13=paru(13)  ! 0.01 effective angular cut-off
      mju(3)=mju(1)
      mju(4)=mju(2)
      nr=np

  140 if(nr.ge.3) then
        pdrmin=2d0*paru12
        do 150 i=n+1,n+nr
          if(i.eq.n+nr.and.iabs(k(n+1,2)).ne.21) goto 150
          i1=i+1
          if(i.eq.n+nr) i1=n+1
          if(k(i,1).eq.41.or.k(i1,1).eq.41) goto 150
          if(mju(1).ne.0.and.i1.lt.mju(1).and.iabs(k(i1,2)).ne.21)
     &    goto 150
          if(mju(2).ne.0.and.i.gt.mju(2).and.iabs(k(i,2)).ne.21)
     &    goto 150
          pap=sqrt((p(i,1)**2+p(i,2)**2+p(i,3)**2)*(p(i1,1)**2+
     &    p(i1,2)**2+p(i1,3)**2))
          pvp=p(i,1)*p(i1,1)+p(i,2)*p(i1,2)+p(i,3)*p(i1,3)
          pdr=4d0*(pap-pvp)**2/max(1d-6,paru13**2*pap+2d0*(pap-pvp))
          if(pdr.lt.pdrmin) then
            ir=i
            pdrmin=pdr
          endif
  150   continue
 
C...Recombine very nearby partons to avoid machine precision problems.
        if(pdrmin.lt.paru12.and.ir.eq.n+nr) then
          do 160 j=1,4
            p(n+1,j)=p(n+1,j)+p(n+nr,j)
  160     continue
          p(n+1,5)=sqrt(max(0d0,p(n+1,4)**2-p(n+1,1)**2-p(n+1,2)**2-
     &    p(n+1,3)**2))
          nr=nr-1
          goto 140
        elseif(pdrmin.lt.paru12) then
          do 170 j=1,4
            p(ir,j)=p(ir,j)+p(ir+1,j)
  170     continue
          p(ir,5)=sqrt(max(0d0,p(ir,4)**2-p(ir,1)**2-p(ir,2)**2-
     &    p(ir,3)**2))
          do 190 i=ir+1,n+nr-1
            k(i,2)=k(i+1,2)
            do 180 j=1,5
              p(i,j)=p(i+1,j)
  180       continue
  190     continue
          if(ir.eq.n+nr-1) k(ir,2)=k(n+nr,2)
          nr=nr-1
          if(mju(1).gt.ir) mju(1)=mju(1)-1
          if(mju(2).gt.ir) mju(2)=mju(2)-1
          goto 140
        endif
      endif
      ntryr=ntryr+1
 
C...Reset particle counter. Skip ahead if no junctions are present;
C...this is usually the case!
      nrs=max(5*nr+11,np)
      ntry=0
  200 ntry=ntry+1
      if(ntry.gt.100.and.ntryr.le.4) then
        paru12=4d0*paru12
        paru13=2d0*paru13
        goto 140
      elseif(ntry.gt.100) then
        call pjerrm(14,'(PYSTRF:) 200 caught in infinite loop')
        if(mstu(21).ge.1) return
      endif

      i=n+nrs
      mstu(90)=mstu90

c...Junction string.
cJAM++
c     if(mju(1).eq.0.and.mju(2).eq.0) goto 580
      if(mju(1).ne.0.or.mju(2).ne.0) goto 5000

 
C...Open versus closed strings. Choose breakup region for latter.
  580 if(mju(1).ne.0.and.mju(2).ne.0) then
        ns=mju(2)-mju(1)
        nb=mju(1)-n
      elseif(mju(1).ne.0) then
        ns=n+nr-mju(1)
        nb=mju(1)-n
      elseif(mju(2).ne.0) then
        ns=mju(2)-n
        nb=1
      elseif(iabs(k(n+1,2)).ne.21) then
        ns=nr-1
        nb=1
      else
        ns=nr+1
        w2sum=0d0
        do 590 is=1,nr
          p(n+nr+is,1)=0.5d0*four(n+is,n+is+1-nr*(is/nr))
          w2sum=w2sum+p(n+nr+is,1)
  590   continue
        w2ran=pjr(0)*w2sum
        nb=0
  600   nb=nb+1
        w2sum=w2sum-p(n+nr+nb,1)
        if(w2sum.gt.w2ran.and.nb.lt.nr) goto 600
      endif
 
C...Find longitudinal string directions (i.e. lightlike four-vectors).
      do 630 is=1,ns
        is1=n+is+nb-1-nr*((is+nb-2)/nr)
        is2=n+is+nb-nr*((is+nb-1)/nr)
        do 610 j=1,5
          dp(1,j)=p(is1,j)
          if(iabs(k(is1,2)).eq.21) dp(1,j)=0.5d0*dp(1,j)
          if(is1.eq.mju(1)) dp(1,j)=pjs(1,j)-pjs(3,j)
          dp(2,j)=p(is2,j)
          if(iabs(k(is2,2)).eq.21) dp(2,j)=0.5d0*dp(2,j)
          if(is2.eq.mju(2)) dp(2,j)=pjs(2,j)-pjs(4,j)
  610   continue
        dp(3,5)=dfour(1,1)
        dp(4,5)=dfour(2,2)
        dhkc=dfour(1,2)
        if(dp(3,5)+2d0*dhkc+dp(4,5).le.0d0) then
          dp(3,5)=dp(1,5)**2
          dp(4,5)=dp(2,5)**2
          dp(1,4)=sqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2+dp(1,5)**2)
          dp(2,4)=sqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2+dp(2,5)**2)
          dhkc=dfour(1,2)
        endif
        dhks=sqrt(dhkc**2-dp(3,5)*dp(4,5))
        dhk1=0.5d0*((dp(4,5)+dhkc)/dhks-1d0)
        dhk2=0.5d0*((dp(3,5)+dhkc)/dhks-1d0)
        in1=n+nr+4*is-3
        p(in1,5)=sqrt(dp(3,5)+2d0*dhkc+dp(4,5))
        do 620 j=1,4
          p(in1,j)=(1d0+dhk1)*dp(1,j)-dhk2*dp(2,j)
          p(in1+1,j)=(1d0+dhk2)*dp(2,j)-dhk1*dp(1,j)
  620   continue
  630 continue
 
c================================================================
C...Begin initialization: sum up energy, set starting position.
c================================================================
      isav=i
      mstu91=mstu(90)
  640 ntry=ntry+1
      if(ntry.gt.100.and.ntryr.le.4) then
        paru12=4d0*paru12
        paru13=2d0*paru13
        goto 140
      elseif(ntry.gt.100) then
        call pjerrm(14,'(PYSTRF:) 640 caught in infinite loop')
        if(mstu(21).ge.1) return
      endif
      i=isav
      mstu(90)=mstu91
      do 660 j=1,4
        p(n+nrs,j)=0d0
        do 650 is=1,nr
          p(n+nrs,j)=p(n+nrs,j)+p(n+is,j)
  650   continue
  660 continue
      do 680 jt=1,2
        irank(jt)=0
        if(mju(jt).ne.0) irank(jt)=njs(jt)
        if(ns.gt.nr) irank(jt)=1
        ie(jt)=k(n+1+(jt/2)*(np-1),3)
        in(3*jt+1)=n+nr+1+4*(jt/2)*(ns-1)
        in(3*jt+2)=in(3*jt+1)+1
        in(3*jt+3)=n+nr+4*ns+2*jt-1
        do 670 in1=n+nr+2+jt,n+nr+4*ns-2+jt,4
          p(in1,1)=2-jt
          p(in1,2)=jt-1
          p(in1,3)=1d0
  670   continue
  680 continue

C.. MOPS variables and switches
      nrvmo=0
      xbmo=1d0
      mstu(121)=0
      mstu(122)=0
 
C...Initialize flavour and pT variables for open string.
      if(ns.lt.nr) then
        px(1)=0d0
        py(1)=0d0
        if(ns.eq.1.and.mju(1)+mju(2).eq.0) call pjptdi(0,px(1),py(1))
        px(2)=-px(1)
        py(2)=-py(1)
        do 690 jt=1,2
          kfl(jt)=k(ie(jt),2)
          if(mju(jt).ne.0) kfl(jt)=kfjs(jt)
          mstj(93)=1
          pmq(jt)=pjmass(kfl(jt))
          gam(jt)=0d0
  690   continue
 
C...Closed string: random initial breakup flavour, pT and vertex.
      else
        kfl(3)=int(1d0+(2d0+parj(2))*pjr(0))*(-1)**int(pjr(0)+0.5d0)
        ibmo=0
  700   call pjkfdi(kfl(3),0,kfl(1),kdump)
C.. Closed string: first vertex diq attempt => enforced second
C.. vertex diq
        if(iabs(kfl(1)).gt.10)then
           ibmo=1
           mstu(121)=0
           goto 700
        endif
        if(ibmo.eq.1) mstu(121)=-1
        kfl(2)=-kfl(1)
        call pjptdi(kfl(1),px(1),py(1))
        px(2)=-px(1)
        py(2)=-py(1)
        pr3=min(25d0,0.1d0*p(n+nr+1,5)**2)
  710   call pjzdis(kfl(1),kfl(2),pr3,z)
        zr=pr3/(z*p(n+nr+1,5)**2)
        if(zr.ge.1d0) goto 710
        do 720 jt=1,2
          mstj(93)=1
          pmq(jt)=pjmass(kfl(jt))
          gam(jt)=pr3*(1d0-z)/z
          in1=n+nr+3+4*(jt/2)*(ns-1)
          p(in1,jt)=1d0-z
          p(in1,3-jt)=jt-1
          p(in1,3)=(2-jt)*(1d0-z)+(jt-1)*z
          p(in1+1,jt)=zr
          p(in1+1,3-jt)=2-jt
          p(in1+1,3)=(2-jt)*(1d0-zr)+(jt-1)*zr
  720   continue
      endif

C.. MOPS variables
      do 730 jt=1,2
         xtmo(jt)=1d0
         pm2qmo(jt)=pmq(jt)**2
         if(iabs(kfl(jt)).gt.10) pm2qmo(jt)=0d0
  730 continue

cJAM:Position of vertex.
c     wstr=sqrt(four(n+nrs,n+nrs))
      zposp(1)=wstr
      zposm(1)=0.d0
      zposp(2)=0.d0
      zposm(2)=wstr

cJAM:Const. quark flag.
      iqconst(1)=jqconst(1)
      iqconst(2)=jqconst(2)
      kfcq(1)=0
      kfcq(2)=0
      kfcq(3)=0
      kfcq(4)=0
      icq(1)=0
      icq(2)=0
      icq(3)=0
      icq(4)=0
      ic=0

C...Find initial transverse directions (i.e. spacelike four-vectors).
      do 770 jt=1,2
        if(jt.eq.1.or.ns.eq.nr-1) then
          in1=in(3*jt+1)
          in3=in(3*jt+3)
          do 740 j=1,4
            dp(1,j)=p(in1,j)
            dp(2,j)=p(in1+1,j)
            dp(3,j)=0d0
            dp(4,j)=0d0
  740     continue
          dp(1,4)=sqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2)
          dp(2,4)=sqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2)
          dp(5,1)=dp(1,1)/dp(1,4)-dp(2,1)/dp(2,4)
          dp(5,2)=dp(1,2)/dp(1,4)-dp(2,2)/dp(2,4)
          dp(5,3)=dp(1,3)/dp(1,4)-dp(2,3)/dp(2,4)
          if(dp(5,1)**2.le.dp(5,2)**2+dp(5,3)**2) dp(3,1)=1d0
          if(dp(5,1)**2.gt.dp(5,2)**2+dp(5,3)**2) dp(3,3)=1d0
          if(dp(5,2)**2.le.dp(5,1)**2+dp(5,3)**2) dp(4,2)=1d0
          if(dp(5,2)**2.gt.dp(5,1)**2+dp(5,3)**2) dp(4,3)=1d0
          dhc12=dfour(1,2)
          dhcx1=dfour(3,1)/dhc12
          dhcx2=dfour(3,2)/dhc12
          dhcxx=1d0/sqrt(1d0+2d0*dhcx1*dhcx2*dhc12)
          dhcy1=dfour(4,1)/dhc12
          dhcy2=dfour(4,2)/dhc12
          dhcyx=dhcxx*(dhcx1*dhcy2+dhcx2*dhcy1)*dhc12
          dhcyy=1d0/sqrt(1d0+2d0*dhcy1*dhcy2*dhc12-dhcyx**2)
          do 750 j=1,4
            dp(3,j)=dhcxx*(dp(3,j)-dhcx2*dp(1,j)-dhcx1*dp(2,j))
            p(in3,j)=dp(3,j)
            p(in3+1,j)=dhcyy*(dp(4,j)-dhcy2*dp(1,j)-dhcy1*dp(2,j)-
     &      dhcyx*dp(3,j))
  750     continue
        else
          do 760 j=1,4
            p(in3+2,j)=p(in3,j)
            p(in3+3,j)=p(in3+1,j)
  760     continue
        endif
  770 continue
 
C...Remove energy used up in junction string fragmentation.
      if(mju(1)+mju(2).gt.0) then
        do 790 jt=1,2
          if(njs(jt).eq.0) goto 790
          do 780 j=1,4
            p(n+nrs,j)=p(n+nrs,j)-pjs(jt+2,j)
  780     continue
  790   continue
      endif
 
c=========================================
C...Produce new particle: side, origin.
c=========================================
      iprod=0 ! JAM
  800 i=i+1
      iprod=iprod+1   !JAM
      if(2*i-nsav.ge.mstu(4)-mstu(32)-5) then
        call pjerrm(11,'(PYSTRF:) no more memory left in PYJETS')
        if(mstu(21).ge.1) return
      endif

C.. New side priority for popcorn systems
      if(mstu(121).le.0)then
         jt=1.5d0+pjr(0)
         if(iabs(kfl(3-jt)).gt.10) jt=3-jt
         if(iabs(kfl(3-jt)).ge.4.and.iabs(kfl(3-jt)).le.8) jt=3-jt
      endif
      jr=3-jt
      js=3-2*jt
      irank(jt)=irank(jt)+1
      k(i,1)=1
      k(i,3)=ie(jt)
      k(i,4)=0
      k(i,5)=iprod ! JAM
 
C...Generate flavour, hadron and pT.
  810 continue
      call pjkfdi(kfl(jt),0,kfl(3),k(i,2))
      if(k(i,2).eq.0) goto 640
      mu90mo=mstu(90)
      if(mstu(121).eq.-1) goto 840

      if(irank(jt).eq.1) then
c...Extra baryon suppression.
        if(iabs(kfl(jt)).le.10.and.iabs(kfl(3)).gt.10) then
          if(pjr(0).gt.parj(19)) goto 810
        endif
cJAM++
c       if(jt.eq.2.and.iabs(kfl(jt)).ge.10
c    $     .and. iqconst(jt).ge.1.and.iqconst(jr).eq.0.
c    $      and. mod(iabs(k(i,2))/1000,10).ne.0) then
c          goto 810
c       endif
cJAM--
      endif

      p(i,5)=pjmass(k(i,2))
      call pjptdi(kfl(jt),px(3),py(3))
      pr(jt)=p(i,5)**2+(px(jt)+px(3))**2+(py(jt)+py(3))**2
 
C...Final hadrons for small invariant mass.
      mstj(93)=1
      pmq(3)=pjmass(kfl(3))
      parjst=parj(33)
      if(mstj(11).eq.2) parjst=parj(34)
      wmin=parjst+pmq(1)+pmq(2)+parj(36)*pmq(3)
      if(iabs(kfl(jt)).gt.10.and.iabs(kfl(3)).gt.10) wmin=
     & wmin-0.5d0*parj(36)*pmq(3)
      wrem2=four(n+nrs,n+nrs)
      if(wrem2.lt.0.10d0) goto 640  ! Try at the first.

cJAM++
      if(iqconst(jt).ge.1) then
c....This is a rank 0 diquark
        if(abs(kfl(jt)).gt.10) then
c.....diquark -> m + x
          if(mod(abs(k(i,2))/1000,10).eq.0) then
             iqconst(jt)=iqconst(jt)-1
             if(iqconst(jt).eq.0.or.iqconst(jt).eq.1) then
               ic=ic+1
               kfcq(ic)=kfl(jt)
               icq(ic)=i
               k(i,4)=1
             endif
c.....diquark -> B + x
          else
            iqconst(jt)=iqconst(jt)-2
            if(iqconst(jt).eq.0) then
              k(i,4)=2
              ic=ic+1
              kfcq(ic)=kfl(jt)
              icq(ic)=i
            else if(iqconst(jt).eq.-1) then
              k(i,4)=1
              ic=ic+1
              kfcq(ic)=kfl(jt)
              icq(ic)=i
            endif
          endif
c.....quark -> h + x
        else
          iqconst(jt)=iqconst(jt)-1
          if(iqconst(jt).eq.0) then
             k(i,4)=1
             ic=ic+1
             kfcq(ic)=kfl(jt)
             icq(ic)=i
          endif
        endif
      endif
cJAM--

c...Fragmentation finish.
      if(wrem2.lt.max(wmin*(1d0+(2d0*pjr(0)-1d0)*parj(37)),
     &parj(32)+pmq(1)+pmq(2))**2) goto 1010
 
C...Choose z, which gives Gamma. Shift z for heavy flavours.
      call pjzdis(kfl(jt),kfl(3),pr(jt),z)
      if(iabs(kfl(jt)).ge.4.and.iabs(kfl(jt)).le.8.and.
     &                                         mstu(90).lt.8) then
        mstu(90)=mstu(90)+1
        mstu(90+mstu(90))=i
        paru(90+mstu(90))=z
      endif

      kfl1a=iabs(kfl(1))
      kfl2a=iabs(kfl(2))

c....Heavy flavours.
      if(max(mod(kfl1a,10),mod(kfl1a/1000,10),mod(kfl2a,10),
     &mod(kfl2a/1000,10)).ge.4) then
        pr(jr)=(pmq(jr)+pmq(3))**2+(px(jr)-px(3))**2+(py(jr)-py(3))**2
        pw12=sqrt(max(0d0,(wrem2-pr(1)-pr(2))**2-4d0*pr(1)*pr(2)))
        z=(wrem2+pr(jt)-pr(jr)+pw12*(2d0*z-1d0))/(2d0*wrem2)
        pr(jr)=(pmq(jr)+parjst)**2+(px(jr)-px(3))**2+(py(jr)-py(3))**2
        if((1d0-z)*(wrem2-pr(jt)/z).lt.pr(jr)) goto 1010
      endif
      gam(3)=(1d0-z)*(gam(jt)+pr(jt)/z)
 
C.. MOPS baryon model modification
      xtmo3=(1d0-z)*xtmo(jt)
      if(iabs(kfl(3)).le.10) nrvmo=0
      if(iabs(kfl(3)).gt.10.and.mstj(12).ge.4) then
         gtstmo=1d0
         ptstmo=1d0
         rtstmo=pjr(0)
         if(iabs(kfl(jt)).le.10)then
            xbmo=min(xtmo3,1d0-(2d-10))
            gbmo=gam(3)
            pmmo=0d0
            pgmo=gbmo+log(1d0-xbmo)*pm2qmo(jt)
            gtstmo=1d0-parf(192)**pgmo
         else
            if(irank(jt).eq.1) then
               gbmo=gam(jt)
               pmmo=0d0
               xbmo=1d0
            endif
            if(xbmo.lt.1d0-(1d-10))then
               pgnmo=gbmo*xtmo3/xbmo+pm2qmo(jt)*log(1d0-xtmo3)
               gtstmo=(1d0-parf(192)**pgnmo)/(1d0-parf(192)**pgmo)
               pgmo=pgnmo
            endif
            if(mstj(12).ge.5)then
               pmnmo=sqrt((xbmo-xtmo3)*(gam(3)/xtmo3-gbmo/xbmo))
               pmmo=pmmo+pmas(jamcomp(k(i,2)),1)-pmas(jamcomp(k(i,2)),3)
               ptstmo=exp((pmmo-pmnmo)*parf(193))
               pmmo=pmnmo
            endif
         endif
 
C.. MOPS Accepting popcorn system hadron.
         if(ptstmo*gtstmo.gt.rtstmo) then
            if(irank(jt).eq.1.or.iabs(kfl(jt)).le.10) then
               nrvmo=i-n-nr
               if(i+nrvmo.gt.mstu(4)-mstu(32)-5) then
                  call pjerrm(11,
     &                 '(PYSTRF:) no more memory left in PYJETS')
                  if(mstu(21).ge.1) return
               endif
               imo=i
               kflmo=kfl(jt)
               pmqmo=pmq(jt)
               pxmo=px(jt)
               pymo=py(jt)
               gammo=gam(jt)
               irmo=irank(jt)
               xmo=xtmo(jt)
               iqmo=iqconst(jt) !JAM
               do 830 j=1,9
                  if(j.le.5) then
                     do 820 line=1,i-n-nr
                        p(mstu(4)-mstu(32)-line,j)=p(n+nr+line,j)
                        k(mstu(4)-mstu(32)-line,j)=k(n+nr+line,j)
  820                continue
                  endif
                  inmo(j)=in(j)
  830          continue
            endif
         else
C..Reject popcorn system, flag=-1 if enforcing new one
            mstu(121)=-1
            if(ptstmo.gt.rtstmo) mstu(121)=-2
         endif
      endif


C..Lift restoring string outside MOPS block
 840  if(mstu(121).lt.0) then
         if(mstu(121).eq.-2) mstu(121)=0
         mstu(90)=mu90mo
         nrvmo=0
         if(irank(jt).eq.1.or.iabs(kfl(jt)).le.10) goto 810
         i=imo
         kfl(jt)=kflmo
         pmq(jt)=pmqmo
         px(jt)=pxmo
         py(jt)=pymo
         gam(jt)=gammo
         irank(jt)=irmo
         xtmo(jt)=xmo
         iqconst(jt)=iqmo  !JAM
         do 860 j=1,9
            if(j.le.5) then
               do 850 line=1,i-n-nr
                  p(n+nr+line,j)=p(mstu(4)-mstu(32)-line,j)
                  k(n+nr+line,j)=k(mstu(4)-mstu(32)-line,j)
 850           continue
            endif
            in(j)=inmo(j)
 860     continue
         goto 810
      endif
      xtmo(jt)=xtmo3
C.. MOPS end of modification

 
      do 870 j=1,3
        in(j)=in(3*jt+j)
  870 continue
 
C...Stepping within or from 'low' string region easy.
      if(in(1)+1.eq.in(2).and.z*p(in(1)+2,3)*p(in(2)+2,3)*
     &                                    p(in(1),5)**2.ge.pr(jt)) then
        p(in(jt)+2,4)=z*p(in(jt)+2,3)
        p(in(jr)+2,4)=pr(jt)/(p(in(jt)+2,4)*p(in(1),5)**2)
        do 880 j=1,4
          p(i,j)=(px(jt)+px(3))*p(in(3),j)+(py(jt)+py(3))*p(in(3)+1,j)
  880   continue
        goto 970
      elseif(in(1)+1.eq.in(2)) then
        p(in(jr)+2,4)=p(in(jr)+2,3)
        p(in(jr)+2,jt)=1d0
        in(jr)=in(jr)+4*js
        if(js*in(jr).gt.js*in(4*jr)) goto 640
        if(four(in(1),in(2)).le.1d-2) then
          p(in(jt)+2,4)=p(in(jt)+2,3)
          p(in(jt)+2,jt)=0d0
          in(jt)=in(jt)+4*js
        endif
      endif
 
C...Find new transverse directions (i.e. spacelike string vectors).
  890 if(js*in(1).gt.js*in(3*jr+1).or.js*in(2).gt.js*in(3*jr+2).or.
     &                                       in(1).gt.in(2)) goto 640
      if(in(1).ne.in(3*jt+1).or.in(2).ne.in(3*jt+2)) then
        do 900 j=1,4
          dp(1,j)=p(in(1),j)
          dp(2,j)=p(in(2),j)
          dp(3,j)=0d0
          dp(4,j)=0d0
  900   continue
        dp(1,4)=sqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2)
        dp(2,4)=sqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2)
        dhc12=dfour(1,2)
        if(dhc12.le.1d-2) then
          p(in(jt)+2,4)=p(in(jt)+2,3)
          p(in(jt)+2,jt)=0d0
          in(jt)=in(jt)+4*js
          goto 890
        endif
        in(3)=n+nr+4*ns+5
        dp(5,1)=dp(1,1)/dp(1,4)-dp(2,1)/dp(2,4)
        dp(5,2)=dp(1,2)/dp(1,4)-dp(2,2)/dp(2,4)
        dp(5,3)=dp(1,3)/dp(1,4)-dp(2,3)/dp(2,4)
        if(dp(5,1)**2.le.dp(5,2)**2+dp(5,3)**2) dp(3,1)=1d0
        if(dp(5,1)**2.gt.dp(5,2)**2+dp(5,3)**2) dp(3,3)=1d0
        if(dp(5,2)**2.le.dp(5,1)**2+dp(5,3)**2) dp(4,2)=1d0
        if(dp(5,2)**2.gt.dp(5,1)**2+dp(5,3)**2) dp(4,3)=1d0
        dhcx1=dfour(3,1)/dhc12
        dhcx2=dfour(3,2)/dhc12
        dhcxx=1d0/sqrt(1d0+2d0*dhcx1*dhcx2*dhc12)
        dhcy1=dfour(4,1)/dhc12
        dhcy2=dfour(4,2)/dhc12
        dhcyx=dhcxx*(dhcx1*dhcy2+dhcx2*dhcy1)*dhc12
        dhcyy=1d0/sqrt(1d0+2d0*dhcy1*dhcy2*dhc12-dhcyx**2)
        do 910 j=1,4
          dp(3,j)=dhcxx*(dp(3,j)-dhcx2*dp(1,j)-dhcx1*dp(2,j))
          p(in(3),j)=dp(3,j)
          p(in(3)+1,j)=dhcyy*(dp(4,j)-dhcy2*dp(1,j)-dhcy1*dp(2,j)-
     &    dhcyx*dp(3,j))
  910   continue
C...Express pT with respect to new axes, if sensible.
        pxp=-(px(3)*four(in(3*jt+3),in(3))+py(3)*
     &  four(in(3*jt+3)+1,in(3)))
        pjp=-(px(3)*four(in(3*jt+3),in(3)+1)+py(3)*
     &  four(in(3*jt+3)+1,in(3)+1))
        if(abs(pxp**2+pjp**2-px(3)**2-py(3)**2).lt.0.01d0) then
          px(3)=pxp
          py(3)=pjp
        endif
      endif
 
C...Sum up known four-momentum. Gives coefficients for m2 expression.
      do 940 j=1,4
        dhg(j)=0d0
        p(i,j)=px(jt)*p(in(3*jt+3),j)+py(jt)*p(in(3*jt+3)+1,j)+
     &  px(3)*p(in(3),j)+py(3)*p(in(3)+1,j)
        do 920 in1=in(3*jt+1),in(1)-4*js,4*js
          p(i,j)=p(i,j)+p(in1+2,3)*p(in1,j)
  920   continue
        do 930 in2=in(3*jt+2),in(2)-4*js,4*js
          p(i,j)=p(i,j)+p(in2+2,3)*p(in2,j)
  930   continue
  940 continue
      dhm(1)=four(i,i)
      dhm(2)=2d0*four(i,in(1))
      dhm(3)=2d0*four(i,in(2))
      dhm(4)=2d0*four(in(1),in(2))
 
C...Find coefficients for Gamma expression.
      do 960 in2=in(1)+1,in(2),4
        do 950 in1=in(1),in2-1,4
          dhc=2d0*four(in1,in2)
          dhg(1)=dhg(1)+p(in1+2,jt)*p(in2+2,jt)*dhc
          if(in1.eq.in(1)) dhg(2)=dhg(2)-js*p(in2+2,jt)*dhc
          if(in2.eq.in(2)) dhg(3)=dhg(3)+js*p(in1+2,jt)*dhc
          if(in1.eq.in(1).and.in2.eq.in(2)) dhg(4)=dhg(4)-dhc
  950   continue
  960 continue
 
C...Solve (m2, Gamma) equation system for energies taken.
      dhs1=dhm(jr+1)*dhg(4)-dhm(4)*dhg(jr+1)
      if(abs(dhs1).lt.1d-4) goto 640
      dhs2=dhm(4)*(gam(3)-dhg(1))-dhm(jt+1)*dhg(jr+1)-dhg(4)*
     &(p(i,5)**2-dhm(1))+dhg(jt+1)*dhm(jr+1)
      dhs3=dhm(jt+1)*(gam(3)-dhg(1))-dhg(jt+1)*(p(i,5)**2-dhm(1))
      p(in(jr)+2,4)=0.5d0*(sqrt(max(0d0,dhs2**2-4d0*dhs1*dhs3))/
     &abs(dhs1)-dhs2/dhs1)
      if(dhm(jt+1)+dhm(4)*p(in(jr)+2,4).le.0d0) goto 640
      p(in(jt)+2,4)=(p(i,5)**2-dhm(1)-dhm(jr+1)*p(in(jr)+2,4))/
     &(dhm(jt+1)+dhm(4)*p(in(jr)+2,4))
 
C...Step to new region if necessary.
      if(p(in(jr)+2,4).gt.p(in(jr)+2,3)) then
        p(in(jr)+2,4)=p(in(jr)+2,3)
        p(in(jr)+2,jt)=1d0
        in(jr)=in(jr)+4*js
        if(js*in(jr).gt.js*in(4*jr)) goto 640
        if(four(in(1),in(2)).le.1d-2) then
          p(in(jt)+2,4)=p(in(jt)+2,3)
          p(in(jt)+2,jt)=0d0
          in(jt)=in(jt)+4*js
        endif
        goto 890
      elseif(p(in(jt)+2,4).gt.p(in(jt)+2,3)) then
        p(in(jt)+2,4)=p(in(jt)+2,3)
        p(in(jt)+2,jt)=0d0
        in(jt)=in(jt)+4*js
        goto 890
      endif
 
C...Four-momentum of particle. Remaining quantities. Loop back.
  970 do 980 j=1,4
        p(i,j)=p(i,j)+p(in(1)+2,4)*p(in(1),j)+p(in(2)+2,4)*p(in(2),j)
        p(n+nrs,j)=p(n+nrs,j)-p(i,j)
  980 continue
      if(p(i,4).lt.p(i,5)) goto 640

cJAM++...Calculate production point.
      zpsave=zposp(jt)
      zmsave=zposm(jt)
      emtr2=p(i,5)**2+p(i,1)**2+p(i,2)**2
c     emtr2=pr(jt)**2
      if(jt.eq.1) then
        zposm(jt)=zposm(jt)+emtr2/(z*zposp(jt))
        zposp(jt)=(1.d0-z)*zposp(jt)
      else
        zposp(jt)=zposp(jt)+emtr2/(z*zposm(jt))
        zposm(jt)=(1.d0-z)*zposm(jt)
      endif

c...Constituent formation point.
      if(mstj(10).eq.1.or.(irank(jt).eq.1.and.jt.eq.2)) then
         vp=zpsave
         vm=zmsave

c...yo-yo formation point.
      else if(mstj(10).eq.2) then
         if(jt.eq.1) then
           vp=zpsave
           vm=zposm(jt)
         else
           vp=zposp(jt)
           vm=zmsave
         endif

c....Mean of two-production points.
      else
         vp=0.5d0*(zpsave+zposp(jt))
         vm=0.5d0*(zmsave+zposm(jt))
      endif
      v(i,1)=0.0d0
      v(i,2)=0.0d0
      v(i,3)=0.5d0*(vp-vm)
      v(i,4)=0.5d0*(vp+vm)
      v(i,5)=sqrt(emtr2)
cJAM--

      kfl(jt)=-kfl(3)
      pmq(jt)=pmq(3)
      px(jt)=-px(3)
      py(jt)=-py(3)
      gam(jt)=gam(3)

      if(in(3).ne.in(3*jt+3)) then
        do 990 j=1,4
          p(in(3*jt+3),j)=p(in(3),j)
          p(in(3*jt+3)+1,j)=p(in(3)+1,j)
  990   continue
      endif
      do 1000 jq=1,2
        in(3*jt+jq)=in(jq)
        p(in(jq)+2,3)=p(in(jq)+2,3)-p(in(jq)+2,4)
        p(in(jq)+2,jt)=p(in(jq)+2,jt)-js*(3-2*jq)*p(in(jq)+2,4)
 1000 continue
c=================
      goto 800
c=================
 
C...Final hadron: side, flavour, hadron, mass.
 1010 i=i+1
      k(i,1)=1
      k(i,3)=ie(jr)
      k(i,4)=0
      k(i,5)=0
      call pjkfdi(kfl(jr),-kfl(3),kfldmp,k(i,2))
      if(k(i,2).eq.0) goto 640
      p(i,5)=pjmass(k(i,2))
      pr(jr)=p(i,5)**2+(px(jr)-px(3))**2+(py(jr)-py(3))**2
cJAM++
      irank(jr)=irank(jr)+1
c     if(irank(jr).eq.1) k(i,4)=kfl(jr)
      k(i,5)=iprod+1
cJAM--
 
C...Final two hadrons: find common setup of four-vectors.
      jq=1
      if(p(in(4)+2,3)*p(in(5)+2,3)*four(in(4),in(5)).lt.p(in(7),3)*
     &p(in(8),3)*four(in(7),in(8))) jq=2
      dhc12=four(in(3*jq+1),in(3*jq+2))
      dhr1=four(n+nrs,in(3*jq+2))/dhc12
      dhr2=four(n+nrs,in(3*jq+1))/dhc12
      if(in(4).ne.in(7).or.in(5).ne.in(8)) then
        px(3-jq)=-four(n+nrs,in(3*jq+3))-px(jq)
        py(3-jq)=-four(n+nrs,in(3*jq+3)+1)-py(jq)
        pr(3-jq)=p(i+(jt+jq-3)**2-1,5)**2+(px(3-jq)+(2*jq-3)*js*
     &  px(3))**2+(py(3-jq)+(2*jq-3)*js*py(3))**2
      endif
 
C...Solve kinematics for final two hadrons, if possible.
      wrem2=wrem2+(px(1)+px(2))**2+(py(1)+py(2))**2
      fd=(sqrt(pr(1))+sqrt(pr(2)))/sqrt(wrem2)

c......Is the transverse mass of the remainder-system smaller than
c......the sum of transverse masses of the final two hadrons?
      if(mju(1)+mju(2).ne.0.and.i.eq.isav+2.and.fd.ge.1d0) goto 200
      if(fd.ge.1d0) goto 640

      fa=wrem2+pr(jt)-pr(jr)
      if(mstj(11).ne.2) prev=0.5d0*exp(max(-50d0,log(fd)*parj(38)*
     &(pr(1)+pr(2))**2))
      if(mstj(11).eq.2) prev=0.5d0*fd**parj(39)

      fb=sign(sqrt(max(0d0,fa**2-4d0*wrem2*pr(jt))),js*(pjr(0)-prev))
      kfl1a=iabs(kfl(1))
      kfl2a=iabs(kfl(2))

      if(max(mod(kfl1a,10),mod(kfl1a/1000,10),mod(kfl2a,10),
     &mod(kfl2a/1000,10)).ge.6) fb=sign(sqrt(max(0d0,fa**2-
     &4d0*wrem2*pr(jt))),dble(js))
      do 1020 j=1,4
        p(i-1,j)=(px(jt)+px(3))*p(in(3*jq+3),j)+(py(jt)+py(3))*
     &  p(in(3*jq+3)+1,j)+0.5d0*(dhr1*(fa+fb)*p(in(3*jq+1),j)+
     &  dhr2*(fa-fb)*p(in(3*jq+2),j))/wrem2
        p(i,j)=p(n+nrs,j)-p(i-1,j)
 1020 continue

      if(p(i-1,4).lt.p(i-1,5).or.p(i,4).lt.p(i,5)) goto 640


cJAM++
      if(iqconst(jr).ge.1) then

c....This is a rank 0 diquark
        if(abs(kfl(jr)).gt.10) then
c.....diquark -> m + x
          if(mod(abs(k(i,2))/1000,10).eq.0) then

             iqconst(jr)=iqconst(jr)-1
             if(iqconst(jr).eq.0.or.iqconst(jr).eq.1) then
               k(i,4)=1
               ic=ic+1
               kfcq(ic)=kfl(jr)
               icq(ic)=i
               if(iqconst(jt).ge.1) then
                 k(i,4)=2
                 ic=ic+1
                 kfcq(ic)=kfl(jt)
                 icq(ic)=i
               endif
             endif

c.....diquark -> B + x
          else
            iqconst(jr)=iqconst(jr)-2
            if(iqconst(jr).eq.0) then
              k(i,4)=2
              ic=ic+1
              kfcq(ic)=kfl(jr)
              icq(ic)=i
            else if(iqconst(jr).eq.-1) then
              k(i,4)=1
              ic=ic+1
              kfcq(ic)=kfl(jr)
              icq(ic)=i
              if(iqconst(jt).ge.1) then
                 k(i,4)=3
                 ic=ic+1
                 kfcq(ic)=kfl(jt)
                 icq(ic)=i
              endif
            endif
          endif

c.....quark -> h + x
        else
          iqconst(jr)=iqconst(jr)-1
          if(iqconst(jr).eq.0) then
            k(i,4)=1
            ic=ic+1
            kfcq(ic)=kfl(jr)
            icq(ic)=i
            if(iqconst(jt).ge.1) then
              k(i,4)=2
              if(abs(kfl(3)).gt.10) k(i,4)=3
              ic=ic+1
              kfcq(ic)=kfl(jt)
              icq(ic)=i
            endif
          endif
        endif

c...Check jt
      else if(iqconst(jt).ge.1) then
c....This is a rank 0 diquark
        if(abs(kfl(jt)).gt.10) then
c.....diquark -> m + x
          if(mod(abs(k(i,2))/1000,10).eq.0) then

             iqconst(jt)=iqconst(jt)-1
             if(iqconst(jt).eq.0.or.iqconst(jt).eq.1) then
               k(i,4)=1
               ic=ic+1
               kfcq(ic)=kfl(jt)
               icq(ic)=i
             endif
c.....diquark -> B + x
          else
            iqconst(jt)=iqconst(jt)-2
            if(iqconst(jt).eq.0) then
              k(i,4)=2
              ic=ic+1
              kfcq(ic)=kfl(jt)
              icq(ic)=i
            else if(iqconst(jt).eq.-1) then
              k(i,4)=1
              ic=ic+1
              kfcq(ic)=kfl(jt)
              icq(ic)=i
            endif
          endif
c.....quark -> h + x
        else
          iqconst(jt)=iqconst(jt)-1
          if(iqconst(jt).eq.0) then
            k(i,4)=1
            ic=ic+1
            kfcq(ic)=kfl(jt)
            icq(ic)=i
          endif
        endif

      endif

c...Calculate positions of final hadrons.
c     if(jt.eq.2) then
        jj1=i-1
        jj2=i
c     else
c       jj1=i
c       jj2=i-1
c     endif

      pp1=p(jj1,4)+p(jj1,3)
      pm1=p(jj1,4)-p(jj1,3)
      pp2=p(jj2,4)+p(jj2,3)
      pm2=p(jj2,4)-p(jj2,3)

      zpnew=zposp(1)-pp1
      zmnew=pm1-zposm(1)

c...Constituent formation point.
      if(mstj(10).eq.1.or.(irank(jt).eq.1.and.jt.eq.2)) then
         vp1=zposp(1)
         vm1=zposm(1)
         vp2=zposp(2)
         vm2=zposm(2)

c...Yo-yo formation point.
      else if(mstj(10).eq.2) then
         vp1=zposp(1)
         vm1=zmnew
         vp2=zpnew
         vm2=zposm(2)

c....Mean of two-production points.
      else
         vp1=0.5d0*(zpnew+zposp(1))
         vm1=0.5d0*(zmnew+zposm(1))
         vp2=0.5d0*(zpnew+zposp(2))
         vm2=0.5d0*(zmnew+zposm(2))
      endif
      v(jj1,1)=0.0d0
      v(jj1,2)=0.0d0
      v(jj1,3)=0.5d0*(vp1-vm1)
      v(jj1,4)=0.5d0*(vp1+vm1)
      v(jj1,5)=sqrt(p(jj1,5)**2+p(jj1,1)**2+p(jj1,2)**2)
      v(jj2,1)=0.0d0
      v(jj2,2)=0.0d0
      v(jj2,3)=0.5d0*(vp2-vm2)
      v(jj2,4)=0.5d0*(vp2+vm2)
      v(jj2,5)=sqrt(p(jj2,5)**2+p(jj2,1)**2+p(jj2,2)**2)
cJAM--

c========================
c...End fragmentation.
c========================
 
C...Mark jets as fragmented and give daughter pointers.
      n=i-nrs+1
      do 1030 i=nsav+1,nsav+np
        im=k(i,3)
        k(im,1)=k(im,1)+10
        if(mstu(16).ne.2) then
          k(im,4)=nsav+1
          k(im,5)=nsav+1
        else
          k(im,4)=nsav+2
          k(im,5)=n
        endif
 1030 continue
 
C...Document string system. Move up particles.
      nsav=nsav+1
      k(nsav,1)=11
      k(nsav,2)=92
      k(nsav,3)=ip
      k(nsav,4)=nsav+1
      k(nsav,5)=n
      do 1040 j=1,4
        p(nsav,j)=dps(j)
        v(nsav,j)=v(ip,j)
 1040 continue
      p(nsav,5)=sqrt(max(0d0,dps(4)**2-dps(1)**2-dps(2)**2-dps(3)**2))
      v(nsav,5)=0d0

      do 1060 i=nsav+1,n
        do 1050 j=1,5
          k(i,j)=k(i+nrs-1,j)
          p(i,j)=p(i+nrs-1,j)
          v(i,j)=v(i+nrs-1,j) !JAM
 1050   continue
cJAM++
        do j=1,3
          if(icq(j).eq.i+nrs-1) then
             icq(j)=i
          endif
        enddo
cJAM--
 1060 continue

      mstu91=mstu(90)
      do 1070 iz=mstu90+1,mstu91
        mstu9t(iz)=mstu(90+iz)-nrs+1-nsav+n
        paru9t(iz)=paru(90+iz)
 1070 continue
      mstu(90)=mstu90
 
C...Order particles in rank along the chain. Update mother pointer.
      do 1090 i=nsav+1,n
        do 1080 j=1,5
          k(i-nsav+n,j)=k(i,j)
          p(i-nsav+n,j)=p(i,j)
          v(i-nsav+n,j)=v(i,j) !JAM
 1080   continue
 1090 continue

      i1=nsav
      do 1120 i=n+1,2*n-nsav
        if(k(i,3).ne.ie(1)) goto 1120
        i1=i1+1
        do 1100 j=1,5
          k(i1,j)=k(i,j)
          p(i1,j)=p(i,j)
          v(i1,j)=v(i,j) ! JAM
 1100   continue
        if(mstu(16).ne.2) k(i1,3)=nsav
        do 1110 iz=mstu90+1,mstu91
          if(mstu9t(iz).eq.i) then
            mstu(90)=mstu(90)+1
            mstu(90+mstu(90))=i1
            paru(90+mstu(90))=paru9t(iz)
          endif
 1110   continue
 1120 continue

      do 1150 i=2*n-nsav,n+1,-1
        if(k(i,3).eq.ie(1)) goto 1150
        i1=i1+1
        do 1130 j=1,5
          k(i1,j)=k(i,j)
          p(i1,j)=p(i,j)
          v(i1,j)=v(i,j) !JAM
 1130   continue
        if(mstu(16).ne.2) k(i1,3)=nsav
        do 1140 iz=mstu90+1,mstu91
          if(mstu9t(iz).eq.i) then
            mstu(90)=mstu(90)+1
            mstu(90+mstu(90))=i1
            paru(90+mstu(90))=paru9t(iz)
          endif
 1140   continue
 1150 continue
 
C...Boost back particle system. Set production vertices.
cJAM++
      if(icms.eq.0)then !JAM
      if(mbst.eq.0) then
cJAM    mstu(33)=1
        mstu(33)=0
        call pjrobo(nsav+1,n,0d0,0d0,dps(1)/dps(4),dps(2)/dps(4),
     &  dps(3)/dps(4))
      else
        dga=dps(4)/wstr !JAM
        do 1160 i=nsav+1,n
          hhpmt=p(i,1)**2+p(i,2)**2+p(i,5)**2
          if(p(i,3).gt.0d0) then
            hhpez=(p(i,4)+p(i,3))*hhbz
            p(i,3)=0.5d0*(hhpez-hhpmt/hhpez)
            p(i,4)=0.5d0*(hhpez+hhpmt/hhpez)
          else
            hhpez=(p(i,4)-p(i,3))/hhbz
            p(i,3)=-0.5d0*(hhpez-hhpmt/hhpez)
            p(i,4)=0.5d0*(hhpez+hhpmt/hhpez)
          endif
          dbv=dps(1)*v(i,1)+dps(2)*v(i,2)+dps(3)*v(i,3) 
          dgabv=dga*(dga*dbv/(1d0+dga)+v(i,4)) 
          v(i,1)=v(i,1)+dgabv*dps(1)/dps(4) !JAM
          v(i,2)=v(i,2)+dgabv*dps(1)/dps(4) !JAM
          v(i,3)=v(i,3)+dgabv*dps(1)/dps(4) !JAM
          v(i,4)=dga*(v(i,4)+dbv)           !JAM
 1160   continue
      endif
      endif
cJAM--

cJAM++
c     do 1180 i=nsav+1,n
c       do 1170 j=1,4
c         v(i,j)=v(ip,j)
c1170   continue
c1180 continue
cJAM--
      mstj(12)=mstj12

      return
c=======================================================================
c.....Junction string.
 5000 continue
      if(mstj(12).ge.4) call pjerrm(29,'(PYSTRF:) sorry,'//
     &     ' junction strings not handled by MSTJ(12)>3 options')
      do 570 jt=1,2
        njs(jt)=0
        if(mju(jt).eq.0) goto 570
        js=3-2*jt
 
C...Find and sum up momentum on three sides of junction. Check flavours.
        do 220 iu=1,3
          iju(iu)=0
          do 210 j=1,5
            pju(iu,j)=0d0
  210     continue
  220   continue
        iu=0
        do 240 i1=n+1+(jt-1)*(nr-1),n+nr+(jt-1)*(1-nr),js
          if(k(i1,2).ne.21.and.iu.le.2) then
            iu=iu+1
            iju(iu)=i1
          endif
          do 230 j=1,4
            pju(iu,j)=pju(iu,j)+p(i1,j)
  230     continue
  240   continue
        do 250 iu=1,3
          pju(iu,5)=sqrt(pju(iu,1)**2+pju(iu,2)**2+pju(iu,3)**2)
  250   continue
        if(k(iju(3),2)/100.ne.10*k(iju(1),2)+k(iju(2),2).and.
     &  k(iju(3),2)/100.ne.10*k(iju(2),2)+k(iju(1),2)) then
          call pjerrm(12,'(PYSTRF:) unphysical flavour combination')
          if(mstu(21).ge.1) return
        endif
 
C...Calculate (approximate) boost to rest frame of junction.
        t12=(pju(1,1)*pju(2,1)+pju(1,2)*pju(2,2)+pju(1,3)*pju(2,3))/
     &  (pju(1,5)*pju(2,5))
        t13=(pju(1,1)*pju(3,1)+pju(1,2)*pju(3,2)+pju(1,3)*pju(3,3))/
     &  (pju(1,5)*pju(3,5))
        t23=(pju(2,1)*pju(3,1)+pju(2,2)*pju(3,2)+pju(2,3)*pju(3,3))/
     &  (pju(2,5)*pju(3,5))
        t11=sqrt((2d0/3d0)*(1d0-t12)*(1d0-t13)/(1d0-t23))
        t22=sqrt((2d0/3d0)*(1d0-t12)*(1d0-t23)/(1d0-t13))
        tsq=sqrt((2d0*t11*t22+t12-1d0)*(1d0+t12))
        t1f=(tsq-t22*(1d0+t12))/(1d0-t12**2)
        t2f=(tsq-t11*(1d0+t12))/(1d0-t12**2)
        do 260 j=1,3
          tju(j)=-(t1f*pju(1,j)/pju(1,5)+t2f*pju(2,j)/pju(2,5))
  260   continue
        tju(4)=sqrt(1d0+tju(1)**2+tju(2)**2+tju(3)**2)
        do 270 iu=1,3
          pju(iu,5)=tju(4)*pju(iu,4)-tju(1)*pju(iu,1)-tju(2)*pju(iu,2)-
     &    tju(3)*pju(iu,3)
  270   continue
 
C...Put junction at rest if motion could give inconsistencies.
        if(pju(1,5)+pju(2,5).gt.pju(1,4)+pju(2,4)) then
          do 280 j=1,3
            tju(j)=0d0
  280     continue
          tju(4)=1d0
          pju(1,5)=pju(1,4)
          pju(2,5)=pju(2,4)
          pju(3,5)=pju(3,4)
        endif
 
C...Start preparing for fragmentation of two strings from junction.
        ista=i
        do 550 iu=1,2
          ns=iju(iu+1)-iju(iu)
 
C...Junction strings: find longitudinal string directions.
          do 310 is=1,ns
            is1=iju(iu)+is-1
            is2=iju(iu)+is
            do 290 j=1,5
              dp(1,j)=0.5d0*p(is1,j)
              if(is.eq.1) dp(1,j)=p(is1,j)
              dp(2,j)=0.5d0*p(is2,j)
              if(is.eq.ns) dp(2,j)=-pju(iu,j)
  290       continue
            if(is.eq.ns) dp(2,4)=sqrt(pju(iu,1)**2+pju(iu,2)**2+
     &      pju(iu,3)**2)
            if(is.eq.ns) dp(2,5)=0d0
            dp(3,5)=dfour(1,1)
            dp(4,5)=dfour(2,2)
            dhkc=dfour(1,2)
            if(dp(3,5)+2d0*dhkc+dp(4,5).le.0d0) then
              dp(1,4)=sqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2)
              dp(2,4)=sqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2)
              dp(3,5)=0d0
              dp(4,5)=0d0
              dhkc=dfour(1,2)
            endif
            dhks=sqrt(dhkc**2-dp(3,5)*dp(4,5))
            dhk1=0.5d0*((dp(4,5)+dhkc)/dhks-1d0)
            dhk2=0.5d0*((dp(3,5)+dhkc)/dhks-1d0)
            in1=n+nr+4*is-3
            p(in1,5)=sqrt(dp(3,5)+2d0*dhkc+dp(4,5))
            do 300 j=1,4
              p(in1,j)=(1d0+dhk1)*dp(1,j)-dhk2*dp(2,j)
              p(in1+1,j)=(1d0+dhk2)*dp(2,j)-dhk1*dp(1,j)
  300       continue
  310     continue
 
C...Junction strings: initialize flavour, momentum and starting pos.
          isav=i
          mstu91=mstu(90)
  320     ntry=ntry+1
          if(ntry.gt.100.and.ntryr.le.4) then
            paru12=4d0*paru12
            paru13=2d0*paru13
            goto 140
          elseif(ntry.gt.100) then
            call pjerrm(14,'(PYSTRF:) 320 caught in infinite loop')
            if(mstu(21).ge.1) return
          endif
          i=isav
          mstu(90)=mstu91
          irankj=0
          ie(1)=k(n+1+(jt/2)*(np-1),3)
          in(4)=n+nr+1
          in(5)=in(4)+1
          in(6)=n+nr+4*ns+1
          do 340 jq=1,2
            do 330 in1=n+nr+2+jq,n+nr+4*ns-2+jq,4
              p(in1,1)=2-jq
              p(in1,2)=jq-1
              p(in1,3)=1d0
  330       continue
  340     continue
          kfl(1)=k(iju(iu),2)
          px(1)=0d0
          py(1)=0d0
          gam(1)=0d0
          do 350 j=1,5
            pju(iu+3,j)=0d0
  350     continue
 
C...Junction strings: find initial transverse directions.
          do 360 j=1,4
            dp(1,j)=p(in(4),j)
            dp(2,j)=p(in(4)+1,j)
            dp(3,j)=0d0
            dp(4,j)=0d0
  360     continue
          dp(1,4)=sqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2)
          dp(2,4)=sqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2)
          dp(5,1)=dp(1,1)/dp(1,4)-dp(2,1)/dp(2,4)
          dp(5,2)=dp(1,2)/dp(1,4)-dp(2,2)/dp(2,4)
          dp(5,3)=dp(1,3)/dp(1,4)-dp(2,3)/dp(2,4)
          if(dp(5,1)**2.le.dp(5,2)**2+dp(5,3)**2) dp(3,1)=1d0
          if(dp(5,1)**2.gt.dp(5,2)**2+dp(5,3)**2) dp(3,3)=1d0
          if(dp(5,2)**2.le.dp(5,1)**2+dp(5,3)**2) dp(4,2)=1d0
          if(dp(5,2)**2.gt.dp(5,1)**2+dp(5,3)**2) dp(4,3)=1d0
          dhc12=dfour(1,2)
          dhcx1=dfour(3,1)/dhc12
          dhcx2=dfour(3,2)/dhc12
          dhcxx=1d0/sqrt(1d0+2d0*dhcx1*dhcx2*dhc12)
          dhcy1=dfour(4,1)/dhc12
          dhcy2=dfour(4,2)/dhc12
          dhcyx=dhcxx*(dhcx1*dhcy2+dhcx2*dhcy1)*dhc12
          dhcyy=1d0/sqrt(1d0+2d0*dhcy1*dhcy2*dhc12-dhcyx**2)
          do 370 j=1,4
            dp(3,j)=dhcxx*(dp(3,j)-dhcx2*dp(1,j)-dhcx1*dp(2,j))
            p(in(6),j)=dp(3,j)
            p(in(6)+1,j)=dhcyy*(dp(4,j)-dhcy2*dp(1,j)-dhcy1*dp(2,j)-
     &      dhcyx*dp(3,j))
  370     continue
 
C...Junction strings: produce new particle, origin.
  380     i=i+1
          if(2*i-nsav.ge.mstu(4)-mstu(32)-5) then
            call pjerrm(11,'(PYSTRF:) no more memory left in PYJETS')
            if(mstu(21).ge.1) return
          endif
          irankj=irankj+1
          k(i,1)=1
          k(i,3)=ie(1)
          k(i,4)=0
          k(i,5)=0
 
C...Junction strings: generate flavour, hadron, pT, z and Gamma.
  390     call pjkfdi(kfl(1),0,kfl(3),k(i,2))
          if(k(i,2).eq.0) goto 320
          if(irankj.eq.1.and.iabs(kfl(1)).le.10.and.
     &    iabs(kfl(3)).gt.10) then
            if(pjr(0).gt.parj(19)) goto 390
          endif
          p(i,5)=pjmass(k(i,2))
          call pjptdi(kfl(1),px(3),py(3))
          pr(1)=p(i,5)**2+(px(1)+px(3))**2+(py(1)+py(3))**2
          call pjzdis(kfl(1),kfl(3),pr(1),z)
          if(iabs(kfl(1)).ge.4.and.iabs(kfl(1)).le.8.and.
     &    mstu(90).lt.8) then
            mstu(90)=mstu(90)+1
            mstu(90+mstu(90))=i
            paru(90+mstu(90))=z
          endif
          gam(3)=(1d0-z)*(gam(1)+pr(1)/z)
          do 400 j=1,3
            in(j)=in(3+j)
  400     continue
 
C...Junction strings: stepping within or from 'low' string region easy.
          if(in(1)+1.eq.in(2).and.z*p(in(1)+2,3)*p(in(2)+2,3)*
     &    p(in(1),5)**2.ge.pr(1)) then
            p(in(1)+2,4)=z*p(in(1)+2,3)
            p(in(2)+2,4)=pr(1)/(p(in(1)+2,4)*p(in(1),5)**2)
            do 410 j=1,4
              p(i,j)=(px(1)+px(3))*p(in(3),j)+(py(1)+py(3))*p(in(3)+1,j)
  410       continue
            goto 500
          elseif(in(1)+1.eq.in(2)) then
            p(in(2)+2,4)=p(in(2)+2,3)
            p(in(2)+2,1)=1d0
            in(2)=in(2)+4
            if(in(2).gt.n+nr+4*ns) goto 320
            if(four(in(1),in(2)).le.1d-2) then
              p(in(1)+2,4)=p(in(1)+2,3)
              p(in(1)+2,1)=0d0
              in(1)=in(1)+4
            endif
          endif
 
C...Junction strings: find new transverse directions.
  420     if(in(1).gt.n+nr+4*ns.or.in(2).gt.n+nr+4*ns.or.
     &    in(1).gt.in(2)) goto 320
          if(in(1).ne.in(4).or.in(2).ne.in(5)) then
            do 430 j=1,4
              dp(1,j)=p(in(1),j)
              dp(2,j)=p(in(2),j)
              dp(3,j)=0d0
              dp(4,j)=0d0
  430       continue
            dp(1,4)=sqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2)
            dp(2,4)=sqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2)
            dhc12=dfour(1,2)
            if(dhc12.le.1d-2) then
              p(in(1)+2,4)=p(in(1)+2,3)
              p(in(1)+2,1)=0d0
              in(1)=in(1)+4
              goto 420
            endif
            in(3)=n+nr+4*ns+5
            dp(5,1)=dp(1,1)/dp(1,4)-dp(2,1)/dp(2,4)
            dp(5,2)=dp(1,2)/dp(1,4)-dp(2,2)/dp(2,4)
            dp(5,3)=dp(1,3)/dp(1,4)-dp(2,3)/dp(2,4)
            if(dp(5,1)**2.le.dp(5,2)**2+dp(5,3)**2) dp(3,1)=1d0
            if(dp(5,1)**2.gt.dp(5,2)**2+dp(5,3)**2) dp(3,3)=1d0
            if(dp(5,2)**2.le.dp(5,1)**2+dp(5,3)**2) dp(4,2)=1d0
            if(dp(5,2)**2.gt.dp(5,1)**2+dp(5,3)**2) dp(4,3)=1d0
            dhcx1=dfour(3,1)/dhc12
            dhcx2=dfour(3,2)/dhc12
            dhcxx=1d0/sqrt(1d0+2d0*dhcx1*dhcx2*dhc12)
            dhcy1=dfour(4,1)/dhc12
            dhcy2=dfour(4,2)/dhc12
            dhcyx=dhcxx*(dhcx1*dhcy2+dhcx2*dhcy1)*dhc12
            dhcyy=1d0/sqrt(1d0+2d0*dhcy1*dhcy2*dhc12-dhcyx**2)
            do 440 j=1,4
              dp(3,j)=dhcxx*(dp(3,j)-dhcx2*dp(1,j)-dhcx1*dp(2,j))
              p(in(3),j)=dp(3,j)
              p(in(3)+1,j)=dhcyy*(dp(4,j)-dhcy2*dp(1,j)-dhcy1*dp(2,j)-
     &        dhcyx*dp(3,j))
  440       continue
C...Express pT with respect to new axes, if sensible.
            pxp=-(px(3)*four(in(6),in(3))+py(3)*four(in(6)+1,in(3)))
            pjp=-(px(3)*four(in(6),in(3)+1)+py(3)*four(in(6)+1,in(3)+1))
            if(abs(pxp**2+pjp**2-px(3)**2-py(3)**2).lt.0.01d0) then
              px(3)=pxp
              py(3)=pjp
            endif
          endif
 
C...Junction strings: sum up known four-momentum, coefficients for m2.
          do 470 j=1,4
            dhg(j)=0d0
            p(i,j)=px(1)*p(in(6),j)+py(1)*p(in(6)+1,j)+px(3)*p(in(3),j)+
     &      py(3)*p(in(3)+1,j)
            do 450 in1=in(4),in(1)-4,4
              p(i,j)=p(i,j)+p(in1+2,3)*p(in1,j)
  450       continue
            do 460 in2=in(5),in(2)-4,4
              p(i,j)=p(i,j)+p(in2+2,3)*p(in2,j)
  460       continue
  470     continue
          dhm(1)=four(i,i)
          dhm(2)=2d0*four(i,in(1))
          dhm(3)=2d0*four(i,in(2))
          dhm(4)=2d0*four(in(1),in(2))
 
C...Junction strings: find coefficients for Gamma expression.
          do 490 in2=in(1)+1,in(2),4
            do 480 in1=in(1),in2-1,4
              dhc=2d0*four(in1,in2)
              dhg(1)=dhg(1)+p(in1+2,1)*p(in2+2,1)*dhc
              if(in1.eq.in(1)) dhg(2)=dhg(2)-p(in2+2,1)*dhc
              if(in2.eq.in(2)) dhg(3)=dhg(3)+p(in1+2,1)*dhc
              if(in1.eq.in(1).and.in2.eq.in(2)) dhg(4)=dhg(4)-dhc
  480       continue
  490     continue
 
C...Junction strings: solve (m2, Gamma) equation system for energies.
          dhs1=dhm(3)*dhg(4)-dhm(4)*dhg(3)
          if(abs(dhs1).lt.1d-4) goto 320
          dhs2=dhm(4)*(gam(3)-dhg(1))-dhm(2)*dhg(3)-dhg(4)*
     &    (p(i,5)**2-dhm(1))+dhg(2)*dhm(3)
          dhs3=dhm(2)*(gam(3)-dhg(1))-dhg(2)*(p(i,5)**2-dhm(1))
          p(in(2)+2,4)=0.5d0*(sqrt(max(0d0,dhs2**2-4d0*dhs1*dhs3))/
     &    abs(dhs1)-dhs2/dhs1)
          if(dhm(2)+dhm(4)*p(in(2)+2,4).le.0d0) goto 320
          p(in(1)+2,4)=(p(i,5)**2-dhm(1)-dhm(3)*p(in(2)+2,4))/
     &    (dhm(2)+dhm(4)*p(in(2)+2,4))
 
C...Junction strings: step to new region if necessary.
          if(p(in(2)+2,4).gt.p(in(2)+2,3)) then
            p(in(2)+2,4)=p(in(2)+2,3)
            p(in(2)+2,1)=1d0
            in(2)=in(2)+4
            if(in(2).gt.n+nr+4*ns) goto 320
            if(four(in(1),in(2)).le.1d-2) then
              p(in(1)+2,4)=p(in(1)+2,3)
              p(in(1)+2,1)=0d0
              in(1)=in(1)+4
            endif
            goto 420
          elseif(p(in(1)+2,4).gt.p(in(1)+2,3)) then
            p(in(1)+2,4)=p(in(1)+2,3)
            p(in(1)+2,1)=0d0
            in(1)=in(1)+js
            goto 890
          endif
 
C...Junction strings: particle four-momentum, remainder, loop back.
  500     do 510 j=1,4
            p(i,j)=p(i,j)+p(in(1)+2,4)*p(in(1),j)+
     &      p(in(2)+2,4)*p(in(2),j)
            pju(iu+3,j)=pju(iu+3,j)+p(i,j)
  510     continue
          if(p(i,4).lt.p(i,5)) goto 320
          pju(iu+3,5)=tju(4)*pju(iu+3,4)-tju(1)*pju(iu+3,1)-
     &    tju(2)*pju(iu+3,2)-tju(3)*pju(iu+3,3)
          if(pju(iu+3,5).lt.pju(iu,5)) then
            kfl(1)=-kfl(3)
            px(1)=-px(3)
            py(1)=-py(3)
            gam(1)=gam(3)
            if(in(3).ne.in(6)) then
              do 520 j=1,4
                p(in(6),j)=p(in(3),j)
                p(in(6)+1,j)=p(in(3)+1,j)
  520         continue
            endif
            do 530 jq=1,2
              in(3+jq)=in(jq)
              p(in(jq)+2,3)=p(in(jq)+2,3)-p(in(jq)+2,4)
              p(in(jq)+2,1)=p(in(jq)+2,1)-(3-2*jq)*p(in(jq)+2,4)
  530       continue
            goto 380
          endif
 
C...Junction strings: save quantities left after each string.
          if(iabs(kfl(1)).gt.10) goto 320
          i=i-1
          kfjh(iu)=kfl(1)
          do 540 j=1,4
            pju(iu+3,j)=pju(iu+3,j)-p(i+1,j)
  540     continue
  550   continue
 
C...Junction strings: put together to new effective string endpoint.
        njs(jt)=i-ista
        kfjs(jt)=k(k(mju(jt+2),3),2)
        kfls=2*int(pjr(0)+3d0*parj(4)/(1d0+3d0*parj(4)))+1
        if(kfjh(1).eq.kfjh(2)) kfls=3
        if(ista.ne.i) kfjs(jt)=isign(1000*max(iabs(kfjh(1)),
     &  iabs(kfjh(2)))+100*min(iabs(kfjh(1)),iabs(kfjh(2)))+
     &  kfls,kfjh(1))
        do 560 j=1,4
          pjs(jt,j)=pju(1,j)+pju(2,j)+p(mju(jt),j)
          pjs(jt+2,j)=pju(4,j)+pju(5,j)
  560   continue
        pjs(jt,5)=sqrt(max(0d0,pjs(jt,4)**2-pjs(jt,1)**2-pjs(jt,2)**2-
     &  pjs(jt,3)**2))
  570 continue
      goto 580

      end
 
C*********************************************************************
 
C...PYINDF
C...Handles the fragmentation of a jet system (or a single
C...jet) according to independent fragmentation models.
 
      subroutine pjindf(ip)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jyjets/,/jydat1/,/jydat2/
C...Local arrays.
      dimension dps(5),psi(4),nfi(3),nfl(3),ifet(3),kflf(3),
     &kflo(2),pxo(2),pyo(2),wo(2)
 
C.. MOPS error message
      if(mstj(12).gt.3) call pjerrm(9,'(PYINDF:) MSTJ(12)>3 options'//
     &' are not treated as expected in independent fragmentation')
 
C...Reset counters. Identify parton system and take copy. Check flavour.
      nsav=n
      mstu90=mstu(90)
      njet=0
      kqsum=0
      do 100 j=1,5
        dps(j)=0d0
  100 continue
      i=ip-1
  110 i=i+1
      if(i.gt.min(n,mstu(4)-mstu(32))) then
        call pjerrm(12,'(PYINDF:) failed to reconstruct jet system')
        if(mstu(21).ge.1) return
      endif
      if(k(i,1).ne.1.and.k(i,1).ne.2) goto 110
      kc=jamcomp(k(i,2))
      if(kc.eq.0) goto 110
      kq=kchg(kc,2)*isign(1,k(i,2))
      if(kq.eq.0) goto 110
      njet=njet+1
      if(kq.ne.2) kqsum=kqsum+kq
      do 120 j=1,5
        k(nsav+njet,j)=k(i,j)
        p(nsav+njet,j)=p(i,j)
        dps(j)=dps(j)+p(i,j)
  120 continue
      k(nsav+njet,3)=i
      if(k(i,1).eq.2.or.(mstj(3).le.5.and.n.gt.i.and.
     &k(i+1,1).eq.2)) goto 110
      if(njet.ne.1.and.kqsum.ne.0) then
        call pjerrm(12,'(PYINDF:) unphysical flavour combination')
        if(mstu(21).ge.1) return
      endif
 
C...Boost copied system to CM frame. Find CM energy and sum flavours.
      if(njet.ne.1) then
        mstu(33)=1
        call pjrobo(nsav+1,nsav+njet,0d0,0d0,-dps(1)/dps(4),
     &  -dps(2)/dps(4),-dps(3)/dps(4))
      endif
      pecm=0d0
      do 130 j=1,3
        nfi(j)=0
  130 continue
      do 140 i=nsav+1,nsav+njet
        pecm=pecm+p(i,4)
        kfa=iabs(k(i,2))
        if(kfa.le.3) then
          nfi(kfa)=nfi(kfa)+isign(1,k(i,2))
        elseif(kfa.gt.1000) then
          kfla=mod(kfa/1000,10)
          kflb=mod(kfa/100,10)
          if(kfla.le.3) nfi(kfla)=nfi(kfla)+isign(1,k(i,2))
          if(kflb.le.3) nfi(kflb)=nfi(kflb)+isign(1,k(i,2))
        endif
  140 continue
 
C...Loop over attempts made. Reset counters.
      ntry=0
  150 ntry=ntry+1
      if(ntry.gt.200) then
        call pjerrm(14,'(PYINDF:) caught in infinite loop')
        if(mstu(21).ge.1) return
      endif
      n=nsav+njet
      mstu(90)=mstu90
      do 160 j=1,3
        nfl(j)=nfi(j)
        ifet(j)=0
        kflf(j)=0
  160 continue
 
C...Loop over jets to be fragmented.
      do 230 ip1=nsav+1,nsav+njet
        mstj(91)=0
        nsav1=n
        mstu91=mstu(90)
 
C...Initial flavour and momentum values. Jet along +z axis.
        kflh=iabs(k(ip1,2))
        if(kflh.gt.10) kflh=mod(kflh/1000,10)
        kflo(2)=0
        wf=p(ip1,4)+sqrt(p(ip1,1)**2+p(ip1,2)**2+p(ip1,3)**2)
 
C...Initial values for quark or diquark jet.
  170   if(iabs(k(ip1,2)).ne.21) then
          nstr=1
          kflo(1)=k(ip1,2)
          call pjptdi(0,pxo(1),pyo(1))
          wo(1)=wf
 
C...Initial values for gluon treated like random quark jet.
        elseif(mstj(2).le.2) then
          nstr=1
          if(mstj(2).eq.2) mstj(91)=1
          kflo(1)=int(1d0+(2d0+parj(2))*pjr(0))*(-1)**int(pjr(0)+0.5d0)
          call pjptdi(0,pxo(1),pyo(1))
          wo(1)=wf
 
C...Initial values for gluon treated like quark-antiquark jet pair,
C...sharing energy according to Altarelli-Parisi splitting function.
        else
          nstr=2
          if(mstj(2).eq.4) mstj(91)=1
          kflo(1)=int(1d0+(2d0+parj(2))*pjr(0))*(-1)**int(pjr(0)+0.5d0)
          kflo(2)=-kflo(1)
          call pjptdi(0,pxo(1),pyo(1))
          pxo(2)=-pxo(1)
          pyo(2)=-pyo(1)
          wo(1)=wf*pjr(0)**(1d0/3d0)
          wo(2)=wf-wo(1)
        endif
 
C...Initial values for rank, flavour, pT and W+.
        do 220 istr=1,nstr
  180     i=n
          mstu(90)=mstu91
          irank=0
          kfl1=kflo(istr)
          px1=pxo(istr)
          py1=pyo(istr)
          w=wo(istr)
 
C...New hadron. Generate flavour and hadron species.
  190     i=i+1
          if(i.ge.mstu(4)-mstu(32)-njet-5) then
            call pjerrm(11,'(PYINDF:) no more memory left in PYJETS')
            if(mstu(21).ge.1) return
          endif
          irank=irank+1
          k(i,1)=1
          k(i,3)=ip1
          k(i,4)=0
          k(i,5)=0
  200     call pjkfdi(kfl1,0,kfl2,k(i,2))
          if(k(i,2).eq.0) goto 180
          if(irank.eq.1.and.iabs(kfl1).le.10.and.iabs(kfl2).gt.10) then
            if(pjr(0).gt.parj(19)) goto 200
          endif
 
C...Find hadron mass. Generate four-momentum.
          p(i,5)=pjmass(k(i,2))
          call pjptdi(kfl1,px2,py2)
          p(i,1)=px1+px2
          p(i,2)=py1+py2
          pr=p(i,5)**2+p(i,1)**2+p(i,2)**2
          call pjzdis(kfl1,kfl2,pr,z)
          mzsav=0
          if(iabs(kfl1).ge.4.and.iabs(kfl1).le.8.and.mstu(90).lt.8) then
            mzsav=1
            mstu(90)=mstu(90)+1
            mstu(90+mstu(90))=i
            paru(90+mstu(90))=z
          endif
          p(i,3)=0.5d0*(z*w-pr/max(1d-4,z*w))
          p(i,4)=0.5d0*(z*w+pr/max(1d-4,z*w))
          if(mstj(3).ge.1.and.irank.eq.1.and.kflh.ge.4.and.
     &    p(i,3).le.0.001d0) then
            if(w.ge.p(i,5)+0.5d0*parj(32)) goto 180
            p(i,3)=0.0001d0
            p(i,4)=sqrt(pr)
            z=p(i,4)/w
          endif
 
C...Remaining flavour and momentum.
          kfl1=-kfl2
          px1=-px2
          py1=-py2
          w=(1d0-z)*w
          do 210 j=1,5
            v(i,j)=0d0
  210     continue
 
C...Check if pL acceptable. Go back for new hadron if enough energy.
          if(mstj(3).ge.0.and.p(i,3).lt.0d0) then
            i=i-1
            if(mzsav.eq.1) mstu(90)=mstu(90)-1
          endif
          if(w.gt.parj(31)) goto 190
          n=i
  220   continue
        if(mod(mstj(3),5).eq.4.and.n.eq.nsav1) wf=wf+0.1d0*parj(32)
        if(mod(mstj(3),5).eq.4.and.n.eq.nsav1) goto 170
 
C...Rotate jet to new direction.
        the=pjangl(p(ip1,3),sqrt(p(ip1,1)**2+p(ip1,2)**2))
        phi=pjangl(p(ip1,1),p(ip1,2))
        mstu(33)=1
        call pjrobo(nsav1+1,n,the,phi,0d0,0d0,0d0)
        k(k(ip1,3),4)=nsav1+1
        k(k(ip1,3),5)=n
 
C...End of jet generation loop. Skip conservation in some cases.
  230 continue
      if(njet.eq.1.or.mstj(3).le.0) goto 490
      if(mod(mstj(3),5).ne.0.and.n-nsav-njet.lt.2) goto 150
 
C...Subtract off produced hadron flavours, finished if zero.
      do 240 i=nsav+njet+1,n
        kfa=iabs(k(i,2))
        kfla=mod(kfa/1000,10)
        kflb=mod(kfa/100,10)
        kflc=mod(kfa/10,10)
        if(kfla.eq.0) then
          if(kflb.le.3) nfl(kflb)=nfl(kflb)-isign(1,k(i,2))*(-1)**kflb
          if(kflc.le.3) nfl(kflc)=nfl(kflc)+isign(1,k(i,2))*(-1)**kflb
        else
          if(kfla.le.3) nfl(kfla)=nfl(kfla)-isign(1,k(i,2))
          if(kflb.le.3) nfl(kflb)=nfl(kflb)-isign(1,k(i,2))
          if(kflc.le.3) nfl(kflc)=nfl(kflc)-isign(1,k(i,2))
        endif
  240 continue
      nreq=(iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3))-iabs(nfl(1)+
     &nfl(2)+nfl(3)))/2+iabs(nfl(1)+nfl(2)+nfl(3))/3
      if(nreq.eq.0) goto 320
 
C...Take away flavour of low-momentum particles until enough freedom.
      nrem=0
  250 irem=0
      p2min=pecm**2
      do 260 i=nsav+njet+1,n
        p2=p(i,1)**2+p(i,2)**2+p(i,3)**2
        if(k(i,1).eq.1.and.p2.lt.p2min) irem=i
        if(k(i,1).eq.1.and.p2.lt.p2min) p2min=p2
  260 continue
      if(irem.eq.0) goto 150
      k(irem,1)=7
      kfa=iabs(k(irem,2))
      kfla=mod(kfa/1000,10)
      kflb=mod(kfa/100,10)
      kflc=mod(kfa/10,10)
      if(kfla.ge.4.or.kflb.ge.4) k(irem,1)=8
      if(k(irem,1).eq.8) goto 250
      if(kfla.eq.0) then
        isgn=isign(1,k(irem,2))*(-1)**kflb
        if(kflb.le.3) nfl(kflb)=nfl(kflb)+isgn
        if(kflc.le.3) nfl(kflc)=nfl(kflc)-isgn
      else
        if(kfla.le.3) nfl(kfla)=nfl(kfla)+isign(1,k(irem,2))
        if(kflb.le.3) nfl(kflb)=nfl(kflb)+isign(1,k(irem,2))
        if(kflc.le.3) nfl(kflc)=nfl(kflc)+isign(1,k(irem,2))
      endif
      nrem=nrem+1
      nreq=(iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3))-iabs(nfl(1)+
     &nfl(2)+nfl(3)))/2+iabs(nfl(1)+nfl(2)+nfl(3))/3
      if(nreq.gt.nrem) goto 250
      do 270 i=nsav+njet+1,n
        if(k(i,1).eq.8) k(i,1)=1
  270 continue
 
C...Find combination of existing and new flavours for hadron.
  280 nfet=2
      if(nfl(1)+nfl(2)+nfl(3).ne.0) nfet=3
      if(nreq.lt.nrem) nfet=1
      if(iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3)).eq.0) nfet=0
      do 290 j=1,nfet
        ifet(j)=1+(iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3)))*pjr(0)
        kflf(j)=isign(1,nfl(1))
        if(ifet(j).gt.iabs(nfl(1))) kflf(j)=isign(2,nfl(2))
        if(ifet(j).gt.iabs(nfl(1))+iabs(nfl(2))) kflf(j)=isign(3,nfl(3))
  290 continue
      if(nfet.eq.2.and.(ifet(1).eq.ifet(2).or.kflf(1)*kflf(2).gt.0))
     &goto 280
      if(nfet.eq.3.and.(ifet(1).eq.ifet(2).or.ifet(1).eq.ifet(3).or.
     &ifet(2).eq.ifet(3).or.kflf(1)*kflf(2).lt.0.or.kflf(1)*kflf(3)
     &.lt.0.or.kflf(1)*(nfl(1)+nfl(2)+nfl(3)).lt.0)) goto 280
      if(nfet.eq.0) kflf(1)=1+int((2d0+parj(2))*pjr(0))
      if(nfet.eq.0) kflf(2)=-kflf(1)
      if(nfet.eq.1) kflf(2)=isign(1+int((2d0+parj(2))*pjr(0)),-kflf(1))
      if(nfet.le.2) kflf(3)=0
      if(kflf(3).ne.0) then
        kflfc=isign(1000*max(iabs(kflf(1)),iabs(kflf(3)))+
     &  100*min(iabs(kflf(1)),iabs(kflf(3)))+1,kflf(1))
        if(kflf(1).eq.kflf(3).or.(1d0+3d0*parj(4))*pjr(0).gt.1d0)
     &  kflfc=kflfc+isign(2,kflfc)
      else
        kflfc=kflf(1)
      endif
      call pjkfdi(kflfc,kflf(2),kfldmp,kf)
      if(kf.eq.0) goto 280
      do 300 j=1,max(2,nfet)
        nfl(iabs(kflf(j)))=nfl(iabs(kflf(j)))-isign(1,kflf(j))
  300 continue
 
C...Store hadron at random among free positions.
      npos=min(1+int(pjr(0)*nrem),nrem)
      do 310 i=nsav+njet+1,n
        if(k(i,1).eq.7) npos=npos-1
        if(k(i,1).eq.1.or.npos.ne.0) goto 310
        k(i,1)=1
        k(i,2)=kf
        p(i,5)=pjmass(k(i,2))
        p(i,4)=sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
  310 continue
      nrem=nrem-1
      nreq=(iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3))-iabs(nfl(1)+
     &nfl(2)+nfl(3)))/2+iabs(nfl(1)+nfl(2)+nfl(3))/3
      if(nrem.gt.0) goto 280
 
C...Compensate for missing momentum in global scheme (3 options).
  320 if(mod(mstj(3),5).ne.0.and.mod(mstj(3),5).ne.4) then
        do 340 j=1,3
          psi(j)=0d0
          do 330 i=nsav+njet+1,n
            psi(j)=psi(j)+p(i,j)
  330     continue
  340   continue
        psi(4)=psi(1)**2+psi(2)**2+psi(3)**2
        pws=0d0
        do 350 i=nsav+njet+1,n
          if(mod(mstj(3),5).eq.1) pws=pws+p(i,4)
          if(mod(mstj(3),5).eq.2) pws=pws+sqrt(p(i,5)**2+(psi(1)*p(i,1)+
     &    psi(2)*p(i,2)+psi(3)*p(i,3))**2/psi(4))
          if(mod(mstj(3),5).eq.3) pws=pws+1d0
  350   continue
        do 370 i=nsav+njet+1,n
          if(mod(mstj(3),5).eq.1) pw=p(i,4)
          if(mod(mstj(3),5).eq.2) pw=sqrt(p(i,5)**2+(psi(1)*p(i,1)+
     &    psi(2)*p(i,2)+psi(3)*p(i,3))**2/psi(4))
          if(mod(mstj(3),5).eq.3) pw=1d0
          do 360 j=1,3
            p(i,j)=p(i,j)-psi(j)*pw/pws
  360     continue
          p(i,4)=sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
  370   continue
 
C...Compensate for missing momentum withing each jet separately.
      elseif(mod(mstj(3),5).eq.4) then
        do 390 i=n+1,n+njet
          k(i,1)=0
          do 380 j=1,5
            p(i,j)=0d0
  380     continue
  390   continue
        do 410 i=nsav+njet+1,n
          ir1=k(i,3)
          ir2=n+ir1-nsav
          k(ir2,1)=k(ir2,1)+1
          pls=(p(i,1)*p(ir1,1)+p(i,2)*p(ir1,2)+p(i,3)*p(ir1,3))/
     &    (p(ir1,1)**2+p(ir1,2)**2+p(ir1,3)**2)
          do 400 j=1,3
            p(ir2,j)=p(ir2,j)+p(i,j)-pls*p(ir1,j)
  400     continue
          p(ir2,4)=p(ir2,4)+p(i,4)
          p(ir2,5)=p(ir2,5)+pls
  410   continue
        pss=0d0
        do 420 i=n+1,n+njet
          if(k(i,1).ne.0) pss=pss+p(i,4)/(pecm*(0.8d0*p(i,5)+0.2d0))
  420   continue
        do 440 i=nsav+njet+1,n
          ir1=k(i,3)
          ir2=n+ir1-nsav
          pls=(p(i,1)*p(ir1,1)+p(i,2)*p(ir1,2)+p(i,3)*p(ir1,3))/
     &    (p(ir1,1)**2+p(ir1,2)**2+p(ir1,3)**2)
          do 430 j=1,3
            p(i,j)=p(i,j)-p(ir2,j)/k(ir2,1)+(1d0/(p(ir2,5)*pss)-1d0)*
     &      pls*p(ir1,j)
  430     continue
          p(i,4)=sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
  440   continue
      endif
 
C...Scale momenta for energy conservation.
      if(mod(mstj(3),5).ne.0) then
        pms=0d0
        pes=0d0
        pqs=0d0
        do 450 i=nsav+njet+1,n
          pms=pms+p(i,5)
          pes=pes+p(i,4)
          pqs=pqs+p(i,5)**2/p(i,4)
  450   continue
        if(pms.ge.pecm) goto 150
        neco=0
  460   neco=neco+1
        pfac=(pecm-pqs)/(pes-pqs)
        pes=0d0
        pqs=0d0
        do 480 i=nsav+njet+1,n
          do 470 j=1,3
            p(i,j)=pfac*p(i,j)
  470     continue
          p(i,4)=sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
          pes=pes+p(i,4)
          pqs=pqs+p(i,5)**2/p(i,4)
  480   continue
        if(neco.lt.10.and.abs(pecm-pes).gt.2d-6*pecm) goto 460
      endif
 
C...Origin of produced particles and parton daughter pointers.
  490 do 500 i=nsav+njet+1,n
        if(mstu(16).ne.2) k(i,3)=nsav+1
        if(mstu(16).eq.2) k(i,3)=k(k(i,3),3)
  500 continue
      do 510 i=nsav+1,nsav+njet
        i1=k(i,3)
        k(i1,1)=k(i1,1)+10
        if(mstu(16).ne.2) then
          k(i1,4)=nsav+1
          k(i1,5)=nsav+1
        else
          k(i1,4)=k(i1,4)-njet+1
          k(i1,5)=k(i1,5)-njet+1
          if(k(i1,5).lt.k(i1,4)) then
            k(i1,4)=0
            k(i1,5)=0
          endif
        endif
  510 continue
 
C...Document independent fragmentation system. Remove copy of jets.
      nsav=nsav+1
      k(nsav,1)=11
      k(nsav,2)=93
      k(nsav,3)=ip
      k(nsav,4)=nsav+1
      k(nsav,5)=n-njet+1
      do 520 j=1,4
        p(nsav,j)=dps(j)
        v(nsav,j)=v(ip,j)
  520 continue
      p(nsav,5)=sqrt(max(0d0,dps(4)**2-dps(1)**2-dps(2)**2-dps(3)**2))
      v(nsav,5)=0d0
      do 540 i=nsav+njet,n
        do 530 j=1,5
          k(i-njet+1,j)=k(i,j)
          p(i-njet+1,j)=p(i,j)
          v(i-njet+1,j)=v(i,j)
  530   continue
  540 continue
      n=n-njet+1
      do 550 iz=mstu90+1,mstu(90)
        mstu(90+iz)=mstu(90+iz)-njet+1
  550 continue
 
C...Boost back particle system. Set production vertices.
      if(njet.ne.1) call pjrobo(nsav+1,n,0d0,0d0,dps(1)/dps(4),
     &dps(2)/dps(4),dps(3)/dps(4))
      do 570 i=nsav+1,n
        do 560 j=1,4
          v(i,j)=v(ip,j)
  560   continue
  570 continue
 
      return
      end
 
C*********************************************************************
 
C...PYDECY
C...Handles the decay of unstable particles.
 
      subroutine pjdecy(ip,icon)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/jydat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
      save /jyjets/,/jydat1/,/jydat2/,/jydat3/
C...Local arrays.
      dimension vdcy(4),kflo(4),kfl1(4),pv(10,5),rord(10),ue(3),be(3),
     &wtcor(10),ptau(4),pcmtau(4),dbetau(3)
cp    character cidc*4
      character cidc*72
      data wtcor/2d0,5d0,15d0,60d0,250d0,1500d0,1.2d4,1.2d5,150d0,16d0/
 
C...Functions: momentum in two-particle decays and four-product.
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2d0*a)
      four(i,j)=p(i,4)*p(j,4)-p(i,1)*p(j,1)-p(i,2)*p(j,2)-p(i,3)*p(j,3)
 
C...Initial values.
      icon=0
      ntry=0
      nsav=n
      kfa=iabs(k(ip,2))
      kfs=isign(1,k(ip,2))
      kc=jamcomp(kfa)
      mstj(92)=0
C...Choose lifetime and determine decay vertex.
      if(k(ip,1).eq.5) then
        v(ip,5)=0d0
      elseif(k(ip,1).ne.4) then
        v(ip,5)=-pmas(kc,4)*log(pjr(0))
      endif
      do 100 j=1,4
        vdcy(j)=v(ip,j)+v(ip,5)*p(ip,j)/p(ip,5)
  100 continue
 
C...Determine whether decay allowed or not.
      mout=0
      if(mstj(22).eq.2) then
        if(pmas(kc,4).gt.parj(71)) mout=1
      elseif(mstj(22).eq.3) then
        if(vdcy(1)**2+vdcy(2)**2+vdcy(3)**2.gt.parj(72)**2) mout=1
      elseif(mstj(22).eq.4) then
        if(vdcy(1)**2+vdcy(2)**2.gt.parj(73)**2) mout=1
        if(abs(vdcy(3)).gt.parj(74)) mout=1
      endif
      if(mout.eq.1.and.k(ip,1).ne.5) then
        k(ip,1)=4
        return
      endif
 
C...Interface to external tau decay library (for tau polarization).
      if(kfa.eq.15.and.mstj(28).ge.1) then
 
C...Starting values for pointers and momenta.
        itau=ip
        do 110 j=1,4
          ptau(j)=p(itau,j)
          pcmtau(j)=p(itau,j)
  110   continue
 
C...Iterate to find position and code of mother of tau.
        imtau=itau
  120   imtau=k(imtau,3)
 
        if(imtau.eq.0) then
C...If no known origin then impossible to do anything further.
          kforig=0
          iorig=0
 
        elseif(k(imtau,2).eq.k(itau,2)) then
C...If tau -> tau + gamma then add gamma energy and loop.
          if(k(k(imtau,4),2).eq.22) then
            do 130 j=1,4
              pcmtau(j)=pcmtau(j)+p(k(imtau,4),j)
  130       continue
          elseif(k(k(imtau,5),2).eq.22) then
            do 140 j=1,4
              pcmtau(j)=pcmtau(j)+p(k(imtau,5),j)
  140       continue
          endif
          goto 120
 
        elseif(iabs(k(imtau,2)).gt.100) then
C...If coming from weak decay of hadron then W is not stored in record,
C...but can be reconstructed by adding neutrino momentum.
          kforig=-isign(24,k(itau,2))
          iorig=0
          do 160 ii=k(imtau,4),k(imtau,5)
            if(k(ii,2)*isign(1,k(itau,2)).eq.-16) then
              do 150 j=1,4
                pcmtau(j)=pcmtau(j)+p(ii,j)
  150         continue
            endif
  160     continue
 
        else
C...If coming from resonance decay then find latest copy of this
C...resonance (may not completely agree).
          kforig=k(imtau,2)
          iorig=imtau
          do 170 ii=imtau+1,ip-1
            if(k(ii,2).eq.kforig.and.k(ii,3).eq.iorig.and.
     &      abs(p(ii,5)-p(iorig,5)).lt.1d-5*p(iorig,5)) iorig=ii
  170     continue
          do 180 j=1,4
            pcmtau(j)=p(iorig,j)
  180     continue
        endif
 
C...Boost tau to rest frame of production process (where known)
C...and rotate it to sit along +z axis.
        do 190 j=1,3
          dbetau(j)=pcmtau(j)/pcmtau(4)
  190   continue
        if(kforig.ne.0) call pjrobo(itau,itau,0d0,0d0,-dbetau(1),
     &  -dbetau(2),-dbetau(3))
        phitau=pjangl(p(itau,1),p(itau,2))
        call pjrobo(itau,itau,0d0,-phitau,0d0,0d0,0d0)
        thetau=pjangl(p(itau,3),p(itau,1))
        call pjrobo(itau,itau,-thetau,0d0,0d0,0d0,0d0)
 
C...Call tau decay routine (if meaningful) and fill extra info.
        if(kforig.ne.0.or.mstj(28).eq.2) then
          call pjtaud(itau,iorig,kforig,ndecay)
          do 200 ii=nsav+1,nsav+ndecay
            k(ii,1)=1
            k(ii,3)=ip
            k(ii,4)=0
            k(ii,5)=0
  200     continue
          n=nsav+ndecay
        endif
 
C...Boost back decay tau and decay products.
        do 210 j=1,4
          p(itau,j)=ptau(j)
  210   continue
        if(kforig.ne.0.or.mstj(28).eq.2) then
          call pjrobo(nsav+1,n,thetau,phitau,0d0,0d0,0d0)
          if(kforig.ne.0) call pjrobo(nsav+1,n,0d0,0d0,dbetau(1),
     &    dbetau(2),dbetau(3))
 
C...Skip past ordinary tau decay treatment.
          mmat=0
          mbst=0
          nd=0
          goto 630
        endif
      endif
 
C...B-Bbar mixing: flip sign of meson appropriately.
      mmix=0
      if((kfa.eq.511.or.kfa.eq.531).and.mstj(26).ge.1) then
        xbbmix=parj(76)
        if(kfa.eq.531) xbbmix=parj(77)
        if(sin(0.5d0*xbbmix*v(ip,5)/pmas(kc,4))**2.gt.pjr(0)) mmix=1
        if(mmix.eq.1) kfs=-kfs
      endif
 
C...Check existence of decay channels. Particle/antiparticle rules.
      kca=kc
      if(mdcy(kc,2).gt.0) then
        mdmdcy=mdme(mdcy(kc,2),2)
        if(mdmdcy.gt.80.and.mdmdcy.le.90) kca=mdmdcy
      endif
      if(mdcy(kca,2).le.0.or.mdcy(kca,3).le.0) then
        call pjerrm(9,'(PYDECY:) no decay channel defined')
        return
      endif
      if(mod(kfa/1000,10).eq.0.and.kca.eq.85) kfs=-kfs
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
 
C...Sum branching ratios of allowed decay channels.
  220 nope=0
      brsu=0d0
      do 230 idl=mdcy(kca,2),mdcy(kca,2)+mdcy(kca,3)-1
        if(mdme(idl,1).ne.1.and.kfsp*mdme(idl,1).ne.2.and.
     &  kfsn*mdme(idl,1).ne.3) goto 230
        if(mdme(idl,2).gt.100) goto 230
        nope=nope+1
        brsu=brsu+brat(idl)
  230 continue
      if(nope.eq.0) then
        call pjerrm(2,'(PYDECY:) all decay channels closed by user')
        icon=1
        return
      endif
 
C...Select decay channel among allowed ones.
  240 rbr=brsu*pjr(0)
      idl=mdcy(kca,2)-1
  250 idl=idl+1
      if(mdme(idl,1).ne.1.and.kfsp*mdme(idl,1).ne.2.and.
     &kfsn*mdme(idl,1).ne.3) then
        if(idl.lt.mdcy(kca,2)+mdcy(kca,3)-1) goto 250
      elseif(mdme(idl,2).gt.100) then
        if(idl.lt.mdcy(kca,2)+mdcy(kca,3)-1) goto 250
      else
        idc=idl
        rbr=rbr-brat(idl)
        if(idl.lt.mdcy(kca,2)+mdcy(kca,3)-1.and.rbr.gt.0d0) goto 250
      endif
 
C...Start readout of decay channel: matrix element, reset counters.
      mmat=mdme(idc,2)
  260 ntry=ntry+1
      if(mod(ntry,200).eq.0) then
        write(cidc,'(i4,1x,i7,1x,g9.3)')ip,k(ip,2),p(ip,5)
cp      write(cidc,'(i4)') idc
        call pjerrm(4,'(PYDECY:) caught in loop for decay channel'//
     &  cidc)
        goto 240
      endif
      if(ntry.gt.2000) then
        write(mstu(11),*)ip,k(ip,2),p(ip,5)
cp      call pjerrm(14,'(PYDECY:) caught in infinite loop')
        call pjerrm(4,'(PYDECY:) caught in infinite loop')
        icon=2
        return
cp      if(mstu(21).ge.1) return
      endif
      i=n
      np=0
      nq=0
      mbst=0
      if(mmat.ge.11.and.p(ip,4).gt.20d0*p(ip,5)) mbst=1
      do 270 j=1,4
        pv(1,j)=0d0
        if(mbst.eq.0) pv(1,j)=p(ip,j)
  270 continue
      if(mbst.eq.1) pv(1,4)=p(ip,5)
      pv(1,5)=p(ip,5)
      ps=0d0
      psq=0d0
      mrem=0
      mhaddy=0
cynara2010/7/2  
c     if(kfa.gt.80) mhaddy=1
      if(kfa.eq.311) mhaddy=1
C.. Random flavour and popcorn system memory.
      irndmo=0
      jtmo=0
      mstu(121)=0
      mstu(125)=10
 
C...Read out decay products. Convert to standard flavour code.
      jtmax=5
      if(mdme(idc+1,2).eq.101) jtmax=10
      do 280 jt=1,jtmax
        if(jt.le.5) kp=kfdp(idc,jt)
        if(jt.ge.6) kp=kfdp(idc+1,jt-5)
        if(kp.eq.0) goto 280
        kpa=iabs(kp)
        kcp=jamcomp(kpa)
cynara2010/7/2  
c       if(kpa.gt.80) mhaddy=1
        if(kpa.eq.311) mhaddy=1
        if(kchg(kcp,3).eq.0.and.kpa.ne.81.and.kpa.ne.82) then
          kfp=kp
        elseif(kpa.ne.81.and.kpa.ne.82) then
          kfp=kfs*kp
        elseif(kpa.eq.81.and.mod(kfa/1000,10).eq.0) then
          kfp=-kfs*mod(kfa/10,10)
        elseif(kpa.eq.81.and.mod(kfa/100,10).ge.mod(kfa/10,10)) then
          kfp=kfs*(100*mod(kfa/10,100)+3)
        elseif(kpa.eq.81) then
          kfp=kfs*(1000*mod(kfa/10,10)+100*mod(kfa/100,10)+1)
        elseif(kp.eq.82) then
          call pjdcyk(-kfs*int(1d0+(2d0+parj(2))*pjr(0)),0,kfp,kdump)
          if(kfp.eq.0) goto 260
          kfp=-kfp
          irndmo=1
          mstj(93)=1
          if(pv(1,5).lt.parj(32)+2d0*pjmass(kfp)) goto 260
        elseif(kp.eq.-82) then
          kfp=mstu(124)
        endif
        if(kpa.eq.81.or.kpa.eq.82) kcp=jamcomp(kfp)
 
C...Add decay product to event record or to quark flavour list.
        kfpa=iabs(kfp)
        kqp=kchg(kcp,2)
        if(mmat.ge.11.and.mmat.le.30.and.kqp.ne.0) then
          nq=nq+1
          kflo(nq)=kfp
C...set rndmflav popcorn system pointer
          if(kp.eq.82.and.mstu(121).gt.0) jtmo=nq
          mstj(93)=2
          psq=psq+pjmass(kflo(nq))
        elseif((mmat.eq.42.or.mmat.eq.43.or.mmat.eq.48).and.np.eq.3.and.
     &    mod(nq,2).eq.1) then
          nq=nq-1
          ps=ps-p(i,5)
          k(i,1)=1
          kfi=k(i,2)
          call pjkfdi(kfp,kfi,kfldmp,k(i,2))
          if(k(i,2).eq.0) goto 260
          mstj(93)=1
          p(i,5)=pjmass(k(i,2))
          ps=ps+p(i,5)
        else
          i=i+1
          np=np+1
          if(mmat.ne.33.and.kqp.ne.0) nq=nq+1
          if(mmat.eq.33.and.kqp.ne.0.and.kqp.ne.2) nq=nq+1
          k(i,1)=1+mod(nq,2)
          if(mmat.eq.4.and.jt.le.2.and.kfp.eq.21) k(i,1)=2
          if(mmat.eq.4.and.jt.eq.3) k(i,1)=1
          k(i,2)=kfp
          k(i,3)=ip
          k(i,4)=0
          k(i,5)=0
          p(i,5)=pjmass(kfp)
          ps=ps+p(i,5)
        endif
  280 continue
 
C...Check masses for resonance decays.
      if(mhaddy.eq.0) then
        if(ps+parj(64).gt.pv(1,5)) goto 240
      endif
 
C...Choose decay multiplicity in phase space model.
  290 if(mmat.ge.11.and.mmat.le.30) then
        psp=ps
        cnde=parj(61)*log(max((pv(1,5)-ps-psq)/parj(62),1.1d0))
        if(mmat.eq.12) cnde=cnde+parj(63)
  300   ntry=ntry+1
C...Reset popcorn flags if new attempt. Re-select rndmflav if failed.
        if(irndmo.eq.0) then
           mstu(121)=0
           jtmo=0
        elseif(irndmo.eq.1) then
           irndmo=2
        else
           goto 260
        endif
        if(ntry.gt.1000) then
          call pjerrm(14,'(PYDECY:) caught in infinite loop')
          if(mstu(21).ge.1) return
        endif
        if(mmat.le.20) then
          gauss=sqrt(-2d0*cnde*log(max(1d-10,pjr(0))))*
     &    sin(paru(2)*pjr(0))
          nd=0.5d0+0.5d0*np+0.25d0*nq+cnde+gauss
          if(nd.lt.np+nq/2.or.nd.lt.2.or.nd.gt.10) goto 300
          if(mmat.eq.13.and.nd.eq.2) goto 300
          if(mmat.eq.14.and.nd.le.3) goto 300
          if(mmat.eq.15.and.nd.le.4) goto 300
        else
          nd=mmat-20
        endif
C.. Set maximum popcorn meson number. Test rndmflav popcorn size.
        mstu(125)=nd-nq/2
        if(mstu(121).gt.mstu(125)) goto 300
 
C...Form hadrons from flavour content.
        do 310 jt=1,4
          kfl1(jt)=kflo(jt)
  310   continue
        if(nd.eq.np+nq/2) goto 330
        do 320 i=n+np+1,n+nd-nq/2
C.. Stick to started popcorn system, else pick side at random
          jt=jtmo
          if(jt.eq.0) jt=1+int((nq-1)*pjr(0))
          call pjdcyk(kfl1(jt),0,kfl2,k(i,2))
          if(k(i,2).eq.0) goto 300
          mstu(125)=mstu(125)-1
          jtmo=0
          if(mstu(121).gt.0) jtmo=jt
          kfl1(jt)=-kfl2
  320   continue
  330   jt=2
        jt2=3
        jt3=4
        if(nq.eq.4.and.pjr(0).lt.parj(66)) jt=4
        if(jt.eq.4.and.isign(1,kfl1(1)*(10-iabs(kfl1(1))))*
     &  isign(1,kfl1(jt)*(10-iabs(kfl1(jt)))).gt.0) jt=3
        if(jt.eq.3) jt2=2
        if(jt.eq.4) jt3=2
        call pjdcyk(kfl1(1),kfl1(jt),kfldmp,k(n+nd-nq/2+1,2))
        if(k(n+nd-nq/2+1,2).eq.0) goto 300
        if(nq.eq.4) call pjdcyk(kfl1(jt2),kfl1(jt3),kfldmp,k(n+nd,2))
        if(nq.eq.4.and.k(n+nd,2).eq.0) goto 300
 
C...Check that sum of decay product masses not too large.
        ps=psp
        do 340 i=n+np+1,n+nd
          k(i,1)=1
          k(i,3)=ip
          k(i,4)=0
          k(i,5)=0
          p(i,5)=pjmass(k(i,2))
          ps=ps+p(i,5)
  340   continue
        if(ps+parj(64).gt.pv(1,5)) goto 300
 
C...Rescale energy to subtract off spectator quark mass.
      elseif((mmat.eq.31.or.mmat.eq.33.or.mmat.eq.44)
     &  .and.np.ge.3) then
        ps=ps-p(n+np,5)
        pqt=(p(n+np,5)+parj(65))/pv(1,5)
        do 350 j=1,5
          p(n+np,j)=pqt*pv(1,j)
          pv(1,j)=(1d0-pqt)*pv(1,j)
  350   continue
        if(ps+parj(64).gt.pv(1,5)) goto 260
        nd=np-1
        mrem=1
 
C...Fully specified final state: check mass broadening effects.
      else
        if(np.ge.2.and.ps+parj(64).gt.pv(1,5)) goto 260
        nd=np
      endif
 
C...Determine position of grandmother, number of sisters.
      nm=0
      kfas=0
      msgn=0
      if(mmat.eq.3) then
        im=k(ip,3)
        if(im.lt.0.or.im.ge.ip) im=0
        if(im.ne.0) kfam=iabs(k(im,2))
        if(im.ne.0) then
          do 360 il=max(ip-2,im+1),min(ip+2,n)
            if(k(il,3).eq.im) nm=nm+1
            if(k(il,3).eq.im.and.il.ne.ip) isis=il
  360     continue
          if(nm.ne.2.or.kfam.le.100.or.mod(kfam,10).ne.1.or.
     &    mod(kfam/1000,10).ne.0) nm=0
          if(nm.eq.2) then
            kfas=iabs(k(isis,2))
            if((kfas.le.100.or.mod(kfas,10).ne.1.or.
     &      mod(kfas/1000,10).ne.0).and.kfas.ne.22) nm=0
          endif
        endif
      endif
 
C...Kinematics of one-particle decays.
      if(nd.eq.1) then
        do 370 j=1,4
          p(n+1,j)=p(ip,j)
  370   continue
        goto 630
      endif
 
C...Calculate maximum weight ND-particle decay.
      pv(nd,5)=p(n+nd,5)
      if(nd.ge.3) then
        wtmax=1d0/wtcor(nd-2)
        pmax=pv(1,5)-ps+p(n+nd,5)
        pmin=0d0
        do 380 il=nd-1,1,-1
          pmax=pmax+p(n+il,5)
          pmin=pmin+p(n+il+1,5)
          wtmax=wtmax*pawt(pmax,pmin,p(n+il,5))
  380   continue
      endif
 
C...Find virtual gamma mass in Dalitz decay.
  390 if(nd.eq.2) then
      elseif(mmat.eq.2) then
        pmes=4d0*pmas(11,1)**2
        pmrho2=pmas(131,1)**2
        pgrho2=pmas(131,2)**2
  400   pmst=pmes*(p(ip,5)**2/pmes)**pjr(0)
        wt=(1+0.5d0*pmes/pmst)*sqrt(max(0d0,1d0-pmes/pmst))*
     &  (1d0-pmst/p(ip,5)**2)**3*(1d0+pgrho2/pmrho2)/
     &  ((1d0-pmst/pmrho2)**2+pgrho2/pmrho2)
        if(wt.lt.pjr(0)) goto 400
        pv(2,5)=max(2.00001d0*pmas(11,1),sqrt(pmst))
 
C...M-generator gives weight. If rejected, try again.
      else
  410   rord(1)=1d0
        do 440 il1=2,nd-1
          rsav=pjr(0)
          do 420 il2=il1-1,1,-1
            if(rsav.le.rord(il2)) goto 430
            rord(il2+1)=rord(il2)
  420     continue
  430     rord(il2+1)=rsav
  440   continue
        rord(nd)=0d0
        wt=1d0
        do 450 il=nd-1,1,-1
          pv(il,5)=pv(il+1,5)+p(n+il,5)+(rord(il)-rord(il+1))*
     &    (pv(1,5)-ps)
          wt=wt*pawt(pv(il,5),pv(il+1,5),p(n+il,5))
  450   continue
        if(wt.lt.pjr(0)*wtmax) goto 410
      endif
 
C...Perform two-particle decays in respective CM frame.
  460 do 480 il=1,nd-1
        pa=pawt(pv(il,5),pv(il+1,5),p(n+il,5))
        ue(3)=2d0*pjr(0)-1d0
        phi=paru(2)*pjr(0)
        ue(1)=sqrt(1d0-ue(3)**2)*cos(phi)
        ue(2)=sqrt(1d0-ue(3)**2)*sin(phi)
        do 470 j=1,3
          p(n+il,j)=pa*ue(j)
          pv(il+1,j)=-pa*ue(j)
  470   continue
        p(n+il,4)=sqrt(pa**2+p(n+il,5)**2)
        pv(il+1,4)=sqrt(pa**2+pv(il+1,5)**2)
  480 continue
 
C...Lorentz transform decay products to lab frame.
      do 490 j=1,4
        p(n+nd,j)=pv(nd,j)
  490 continue
      do 530 il=nd-1,1,-1
        do 500 j=1,3
          be(j)=pv(il,j)/pv(il,4)
  500   continue
        ga=pv(il,4)/pv(il,5)
        do 520 i=n+il,n+nd
          bep=be(1)*p(i,1)+be(2)*p(i,2)+be(3)*p(i,3)
          do 510 j=1,3
            p(i,j)=p(i,j)+ga*(ga*bep/(1d0+ga)+p(i,4))*be(j)
  510     continue
          p(i,4)=ga*(p(i,4)+bep)
  520   continue
  530 continue
 
C...Check that no infinite loop in matrix element weight.
      ntry=ntry+1
      if(ntry.gt.800) goto 560
 
C...Matrix elements for omega and phi decays.
      if(mmat.eq.1) then
        wt=(p(n+1,5)*p(n+2,5)*p(n+3,5))**2-(p(n+1,5)*four(n+2,n+3))**2
     &  -(p(n+2,5)*four(n+1,n+3))**2-(p(n+3,5)*four(n+1,n+2))**2
     &  +2d0*four(n+1,n+2)*four(n+1,n+3)*four(n+2,n+3)
        if(max(wt*wtcor(9)/p(ip,5)**6,0.001d0).lt.pjr(0)) goto 390
 
C...Matrix elements for pi0 or eta Dalitz decay to gamma e+ e-.
      elseif(mmat.eq.2) then
        four12=four(n+1,n+2)
        four13=four(n+1,n+3)
        wt=(pmst-0.5d0*pmes)*(four12**2+four13**2)+
     &  pmes*(four12*four13+four12**2+four13**2)
        if(wt.lt.pjr(0)*0.25d0*pmst*(p(ip,5)**2-pmst)**2) goto 460
 
C...Matrix element for S0 -> S1 + V1 -> S1 + S2 + S3 (S scalar,
C...V vector), of form cos**2(theta02) in V1 rest frame, and for
C...S0 -> gamma + V1 -> gamma + S2 + S3, of form sin**2(theta02).
      elseif(mmat.eq.3.and.nm.eq.2) then
        four10=four(ip,im)
        four12=four(ip,n+1)
        four02=four(im,n+1)
        pms1=p(ip,5)**2
        pms0=p(im,5)**2
        pms2=p(n+1,5)**2
        if(kfas.ne.22) hnum=(four10*four12-pms1*four02)**2
        if(kfas.eq.22) hnum=pms1*(2d0*four10*four12*four02-
     &  pms1*four02**2-pms0*four12**2-pms2*four10**2+pms1*pms0*pms2)
        hnum=max(1d-6*pms1**2*pms0*pms2,hnum)
        hden=(four10**2-pms1*pms0)*(four12**2-pms1*pms2)
        if(hnum.lt.pjr(0)*hden) goto 460
 
C...Matrix element for "onium" -> g + g + g or gamma + g + g.
      elseif(mmat.eq.4) then
        hx1=2d0*four(ip,n+1)/p(ip,5)**2
        hx2=2d0*four(ip,n+2)/p(ip,5)**2
        hx3=2d0*four(ip,n+3)/p(ip,5)**2
        wt=((1d0-hx1)/(hx2*hx3))**2+((1d0-hx2)/(hx1*hx3))**2+
     &  ((1d0-hx3)/(hx1*hx2))**2
        if(wt.lt.2d0*pjr(0)) goto 390
        if(k(ip+1,2).eq.22.and.(1d0-hx1)*p(ip,5)**2.lt.4d0*parj(32)**2)
     &  goto 390
 
C...Effective matrix element for nu spectrum in tau -> nu + hadrons.
      elseif(mmat.eq.41) then
        hx1=2d0*four(ip,n+1)/p(ip,5)**2
        hxm=min(0.75d0,2d0*(1d0-ps/p(ip,5)))
        if(hx1*(3d0-2d0*hx1).lt.pjr(0)*hxm*(3d0-2d0*hxm)) goto 390
 
C...Matrix elements for weak decays (only semileptonic for c and b)
      elseif((mmat.eq.42.or.mmat.eq.43.or.mmat.eq.44.or.mmat.eq.48)
     &  .and.nd.eq.3) then
        if(mbst.eq.0) wt=four(ip,n+1)*four(n+2,n+3)
        if(mbst.eq.1) wt=p(ip,5)*p(n+1,4)*four(n+2,n+3)
        if(wt.lt.pjr(0)*p(ip,5)*pv(1,5)**3/wtcor(10)) goto 390
      elseif(mmat.eq.42.or.mmat.eq.43.or.mmat.eq.44.or.mmat.eq.48) then
        do 550 j=1,4
          p(n+np+1,j)=0d0
          do 540 is=n+3,n+np
            p(n+np+1,j)=p(n+np+1,j)+p(is,j)
  540     continue
  550   continue
        if(mbst.eq.0) wt=four(ip,n+1)*four(n+2,n+np+1)
        if(mbst.eq.1) wt=p(ip,5)*p(n+1,4)*four(n+2,n+np+1)
        if(wt.lt.pjr(0)*p(ip,5)*pv(1,5)**3/wtcor(10)) goto 390
      endif
 
C...Scale back energy and reattach spectator.
  560 if(mrem.eq.1) then
        do 570 j=1,5
          pv(1,j)=pv(1,j)/(1d0-pqt)
  570   continue
        nd=nd+1
        mrem=0
      endif
 
C...Low invariant mass for system with spectator quark gives particle,
C...not two jets. Readjust momenta accordingly.
      if(mmat.eq.31.and.nd.eq.3) then
        mstj(93)=1
        pm2=pjmass(k(n+2,2))
        mstj(93)=1
        pm3=pjmass(k(n+3,2))
        if(p(n+2,5)**2+p(n+3,5)**2+2d0*four(n+2,n+3).ge.
     &  (parj(32)+pm2+pm3)**2) goto 630
        k(n+2,1)=1
        kftemp=k(n+2,2)
        call pjkfdi(kftemp,k(n+3,2),kfldmp,k(n+2,2))
        if(k(n+2,2).eq.0) goto 260
        p(n+2,5)=pjmass(k(n+2,2))
        ps=p(n+1,5)+p(n+2,5)
        pv(2,5)=p(n+2,5)
        mmat=0
        nd=2
        goto 460
      elseif(mmat.eq.44) then
        mstj(93)=1
        pm3=pjmass(k(n+3,2))
        mstj(93)=1
        pm4=pjmass(k(n+4,2))
        if(p(n+3,5)**2+p(n+4,5)**2+2d0*four(n+3,n+4).ge.
     &  (parj(32)+pm3+pm4)**2) goto 600
        k(n+3,1)=1
        kftemp=k(n+3,2)
        call pjkfdi(kftemp,k(n+4,2),kfldmp,k(n+3,2))
        if(k(n+3,2).eq.0) goto 260
        p(n+3,5)=pjmass(k(n+3,2))
        do 580 j=1,3
          p(n+3,j)=p(n+3,j)+p(n+4,j)
  580   continue
        p(n+3,4)=sqrt(p(n+3,1)**2+p(n+3,2)**2+p(n+3,3)**2+p(n+3,5)**2)
        ha=p(n+1,4)**2-p(n+2,4)**2
        hb=ha-(p(n+1,5)**2-p(n+2,5)**2)
        hc=(p(n+1,1)-p(n+2,1))**2+(p(n+1,2)-p(n+2,2))**2+
     &  (p(n+1,3)-p(n+2,3))**2
        hd=(pv(1,4)-p(n+3,4))**2
        he=ha**2-2d0*hd*(p(n+1,4)**2+p(n+2,4)**2)+hd**2
        hf=hd*hc-hb**2
        hg=hd*hc-ha*hb
        hh=(sqrt(hg**2+he*hf)-hg)/(2d0*hf)
        do 590 j=1,3
          pcor=hh*(p(n+1,j)-p(n+2,j))
          p(n+1,j)=p(n+1,j)+pcor
          p(n+2,j)=p(n+2,j)-pcor
  590   continue
        p(n+1,4)=sqrt(p(n+1,1)**2+p(n+1,2)**2+p(n+1,3)**2+p(n+1,5)**2)
        p(n+2,4)=sqrt(p(n+2,1)**2+p(n+2,2)**2+p(n+2,3)**2+p(n+2,5)**2)
        nd=nd-1
      endif
 
C...Check invariant mass of W jets. May give one particle or start over.
  600 if((mmat.eq.42.or.mmat.eq.43.or.mmat.eq.44.or.mmat.eq.48)
     &.and.iabs(k(n+1,2)).lt.10) then
        pmr=sqrt(max(0d0,p(n+1,5)**2+p(n+2,5)**2+2d0*four(n+1,n+2)))
        mstj(93)=1
        pm1=pjmass(k(n+1,2))
        mstj(93)=1
        pm2=pjmass(k(n+2,2))
        if(pmr.gt.parj(32)+pm1+pm2) goto 610
        kfldum=int(1.5d0+pjr(0))
        call pjkfdi(k(n+1,2),-isign(kfldum,k(n+1,2)),kfldmp,kf1)
        call pjkfdi(k(n+2,2),-isign(kfldum,k(n+2,2)),kfldmp,kf2)
        if(kf1.eq.0.or.kf2.eq.0) goto 260
        psm=pjmass(kf1)+pjmass(kf2)
        if((mmat.eq.42.or.mmat.eq.48).and.pmr.gt.parj(64)+psm) goto 610
        if(mmat.ge.43.and.pmr.gt.0.2d0*parj(32)+psm) goto 610
        if(mmat.eq.48) goto 390
        if(nd.eq.4.or.kfa.eq.15) goto 260
        k(n+1,1)=1
        kftemp=k(n+1,2)
        call pjkfdi(kftemp,k(n+2,2),kfldmp,k(n+1,2))
        if(k(n+1,2).eq.0) goto 260
        p(n+1,5)=pjmass(k(n+1,2))
        k(n+2,2)=k(n+3,2)
        p(n+2,5)=p(n+3,5)
        ps=p(n+1,5)+p(n+2,5)
        if(ps+parj(64).gt.pv(1,5)) goto 260
        pv(2,5)=p(n+3,5)
        mmat=0
        nd=2
        goto 460
      endif
 
C...Phase space decay of partons from W decay.
  610 if((mmat.eq.42.or.mmat.eq.48).and.iabs(k(n+1,2)).lt.10) then
        kflo(1)=k(n+1,2)
        kflo(2)=k(n+2,2)
        k(n+1,1)=k(n+3,1)
        k(n+1,2)=k(n+3,2)
        do 620 j=1,5
          pv(1,j)=p(n+1,j)+p(n+2,j)
          p(n+1,j)=p(n+3,j)
  620   continue
        pv(1,5)=pmr
        n=n+1
        np=0
        nq=2
        ps=0d0
        mstj(93)=2
        psq=pjmass(kflo(1))
        mstj(93)=2
        psq=psq+pjmass(kflo(2))
        mmat=11
        goto 290
      endif
 
C...Boost back for rapidly moving particle.
  630 n=n+nd
      if(mbst.eq.1) then
        do 640 j=1,3
          be(j)=p(ip,j)/p(ip,4)
  640   continue
        ga=p(ip,4)/p(ip,5)
        do 660 i=nsav+1,n
          bep=be(1)*p(i,1)+be(2)*p(i,2)+be(3)*p(i,3)
          do 650 j=1,3
            p(i,j)=p(i,j)+ga*(ga*bep/(1d0+ga)+p(i,4))*be(j)
  650     continue
          p(i,4)=ga*(p(i,4)+bep)
  660   continue
      endif
 
C...Fill in position of decay vertex.
      do 680 i=nsav+1,n
        do 670 j=1,4
          v(i,j)=vdcy(j)
  670   continue
        v(i,5)=0d0
  680 continue
 
C...Set up for parton shower evolution from jets.
      if(mstj(23).ge.1.and.mmat.eq.4.and.k(nsav+1,2).eq.21) then
        k(nsav+1,1)=3
        k(nsav+2,1)=3
        k(nsav+3,1)=3
        k(nsav+1,4)=mstu(5)*(nsav+2)
        k(nsav+1,5)=mstu(5)*(nsav+3)
        k(nsav+2,4)=mstu(5)*(nsav+3)
        k(nsav+2,5)=mstu(5)*(nsav+1)
        k(nsav+3,4)=mstu(5)*(nsav+1)
        k(nsav+3,5)=mstu(5)*(nsav+2)
        mstj(92)=-(nsav+1)
      elseif(mstj(23).ge.1.and.mmat.eq.4) then
        k(nsav+2,1)=3
        k(nsav+3,1)=3
        k(nsav+2,4)=mstu(5)*(nsav+3)
        k(nsav+2,5)=mstu(5)*(nsav+3)
        k(nsav+3,4)=mstu(5)*(nsav+2)
        k(nsav+3,5)=mstu(5)*(nsav+2)
        mstj(92)=nsav+2
      elseif(mstj(23).ge.1.and.(mmat.eq.32.or.mmat.eq.44).and.
     &  iabs(k(nsav+1,2)).le.10.and.iabs(k(nsav+2,2)).le.10) then
        k(nsav+1,1)=3
        k(nsav+2,1)=3
        k(nsav+1,4)=mstu(5)*(nsav+2)
        k(nsav+1,5)=mstu(5)*(nsav+2)
        k(nsav+2,4)=mstu(5)*(nsav+1)
        k(nsav+2,5)=mstu(5)*(nsav+1)
        mstj(92)=nsav+1
      elseif(mstj(23).ge.1.and.(mmat.eq.32.or.mmat.eq.44).and.
     &  iabs(k(nsav+1,2)).le.20.and.iabs(k(nsav+2,2)).le.20) then
        mstj(92)=nsav+1
      elseif(mstj(23).ge.1.and.mmat.eq.33.and.iabs(k(nsav+2,2)).eq.21)
     &  then
        k(nsav+1,1)=3
        k(nsav+2,1)=3
        k(nsav+3,1)=3
        kcp=jamcomp(k(nsav+1,2))
        kqp=kchg(kcp,2)*isign(1,k(nsav+1,2))
        jcon=4
        if(kqp.lt.0) jcon=5
        k(nsav+1,jcon)=mstu(5)*(nsav+2)
        k(nsav+2,9-jcon)=mstu(5)*(nsav+1)
        k(nsav+2,jcon)=mstu(5)*(nsav+3)
        k(nsav+3,9-jcon)=mstu(5)*(nsav+2)
        mstj(92)=nsav+1
      elseif(mstj(23).ge.1.and.mmat.eq.33) then
        k(nsav+1,1)=3
        k(nsav+3,1)=3
        k(nsav+1,4)=mstu(5)*(nsav+3)
        k(nsav+1,5)=mstu(5)*(nsav+3)
        k(nsav+3,4)=mstu(5)*(nsav+1)
        k(nsav+3,5)=mstu(5)*(nsav+1)
        mstj(92)=nsav+1
      endif
 
C...Mark decayed particle; special option for B-Bbar mixing.
      if(k(ip,1).eq.5) k(ip,1)=15
      if(k(ip,1).le.10) k(ip,1)=11
      if(mmix.eq.1.and.mstj(26).eq.2.and.k(ip,1).eq.11) k(ip,1)=12
      k(ip,4)=nsav+1
      k(ip,5)=n
 
      return
      end
 
C*********************************************************************
 
C...PYDCYK
C...Handles flavour production in the decay of unstable particles
C...and small string clusters.
 
      subroutine pjdcyk(kfl1,kfl2,kfl3,kf)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jydat1/,/jydat2/
 

C.. Call PYKFDI directly if no popcorn option is on
      if(mstj(12).lt.2) then
         call pjkfdi(kfl1,kfl2,kfl3,kf)
         mstu(124)=kfl3
         return
      endif
 
      kfl3=0
      kf=0
      if(kfl1.eq.0) return
      kf1a=iabs(kfl1)
      kf2a=iabs(kfl2)
 
      nsto=130
      nmax=min(mstu(125),10)
 
C.. Identify rank 0 cluster qq
      irank=1
      if(kf1a.gt.10.and.kf1a.lt.10000) irank=0
 
      if(kf2a.gt.0)then
C.. Join jets: Fails if store not empty
         if(mstu(121).gt.0) then
            mstu(121)=0
            return
         endif
         call pjkfdi(kfl1,kfl2,kfl3,kf)
      elseif(kf1a.gt.10.and.mstu(121).gt.0)then
C.. Pick popcorn meson from store, return same qq, decrease store
         kf=mstu(nsto+mstu(121))
         kfl3=-kfl1
         mstu(121)=mstu(121)-1
      else
C.. Generate new flavour. Then done if no diquark is generated
  100    call pjkfdi(kfl1,0,kfl3,kf)
         if(mstu(121).eq.-1) goto 100
         mstu(124)=kfl3
         if(kf.eq.0.or.iabs(kfl3).le.10) return
 
C.. Simple case if no dynamical popcorn suppressions are considered
         if(mstj(12).lt.4) then
            if(mstu(121).eq.0) return
            nmes=1
            kfprev=-kfl3
            call pjkfdi(kfprev,0,kfl3,kfm)
C.. Due to eta+eta' suppr., a qq->M+qq attempt might end as qq->B+q
            if(iabs(kfl3).le.10)then
               kfl3=-kfprev
               return
            endif
            goto 120
         endif
 
C test output qq against fake Gamma, then return if no popcorn.
         gb=2d0
         if(irank.ne.0)then
            call pjzdis(1,2103,5d0,z)
            gb=3d0*(1d0-z)/z
            if(1d0-parf(192)**gb.lt.pjr(0)) then
               mstu(121)=0
               goto 100
            endif
         endif      
         if(mstu(121).eq.0) return
 
C..Set store size memory. Pick fake dynamical variables of qq.
         nmes=mstu(121)
         call pjptdi(1,px3,py3)
         x=1d0
         popm=0d0
         g=gb
         popg=gb
 
C.. Pick next popcorn meson, test with fake dynamical variables
  110    kfprev=-kfl3
         px1=-px3
         py1=-py3
         call pjkfdi(kfprev,0,kfl3,kfm)
         if(mstu(121).eq.-1) goto 100
         call pjptdi(kfl3,px3,py3)
         pm=pjmass(kfm)**2+(px1+px3)**2+(py1+py3)**2
         call pjzdis(kfprev,kfl3,pm,z)
         g=(1d0-z)*(g+pm/z)
         x=(1d0-z)*x
 
         ptst=1d0
         gtst=1d0
         rtst=pjr(0)
         if(mstj(12).gt.4)then
            popmn=sqrt((1d0-x)*(g/x-gb))
            popm=popm+pmas(jamcomp(kfm),1)-pmas(jamcomp(kfm),3)
            ptst=exp((popm-popmn)*parf(193))
            popm=popmn
         endif
         if(irank.ne.0)then
            popgn=x*gb
            gtst=(1d0-parf(192)**popgn)/(1d0-parf(192)**popg)
            popg=popgn
         endif
         if(rtst.gt.ptst*gtst)then
            mstu(121)=0
            if(rtst.gt.ptst) mstu(121)=-1
            goto 100
         endif
 
C.. Store meson
  120    if(nmes.le.nmax) mstu(nsto+mstu(121)+1)=kfm
         if(mstu(121).gt.0) goto 110
 
C.. Test accepted system size. If OK set global popcorn size variable.
         if(nmes.gt.nmax)then
            kf=0
            kfl3=0
            return
         endif
         mstu(121)=nmes
      endif
 
      return
      end
 
C********************************************************************
 
C...PYKFDI
C...Generates a new flavour pair and combines off a hadron
 
      subroutine pjkfdi(kfl1,kfl2,kfl3,kf)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jydat1/,/jydat2/
      common/jydat4/chaf(500,2)
      character chaf*16
 
      if(mstu(123).eq.0.and.mstj(12).gt.0)  call pjkfin
 
C...Default flavour values. Input consistency checks.
      kf1a=iabs(kfl1)
      kf2a=iabs(kfl2)
      kfl3=0
      kf=0
      if(kf1a.eq.0) return
      if(kf2a.ne.0)then
        if(kf1a.le.10.and.kf2a.le.10.and.kfl1*kfl2.gt.0) return
        if(kf1a.gt.10.and.kf2a.gt.10) return
        if((kf1a.gt.10.or.kf2a.gt.10).and.kfl1*kfl2.lt.0) return
      endif
 
C...Check if tabulated flavour probabilities are to be used.
      if(mstj(15).eq.1) then
        if(mstj(12).ge.5)  call pjerrm(29,
     &        '(PYKFDI:) Sorry, option MSTJ(15)=1 not available' //
     &        ' together with MSTJ(12)>=5 modification')
        ktab1=-1
        if(kf1a.ge.1.and.kf1a.le.6) ktab1=kf1a
        kfl1a=mod(kf1a/1000,10)
        kfl1b=mod(kf1a/100,10)
        kfl1s=mod(kf1a,10)
        if(kfl1a.ge.1.and.kfl1a.le.4.and.kfl1b.ge.1.and.kfl1b.le.4)
     &  ktab1=6+kfl1a*(kfl1a-2)+2*kfl1b+(kfl1s-1)/2
        if(kfl1a.ge.1.and.kfl1a.le.4.and.kfl1a.eq.kfl1b) ktab1=ktab1-1
        if(kf1a.ge.1.and.kf1a.le.6) kfl1a=kf1a
        ktab2=0
        if(kf2a.ne.0) then
          ktab2=-1
          if(kf2a.ge.1.and.kf2a.le.6) ktab2=kf2a
          kfl2a=mod(kf2a/1000,10)
          kfl2b=mod(kf2a/100,10)
          kfl2s=mod(kf2a,10)
          if(kfl2a.ge.1.and.kfl2a.le.4.and.kfl2b.ge.1.and.kfl2b.le.4)
     &    ktab2=6+kfl2a*(kfl2a-2)+2*kfl2b+(kfl2s-1)/2
          if(kfl2a.ge.1.and.kfl2a.le.4.and.kfl2a.eq.kfl2b) ktab2=ktab2-1
        endif
        if(ktab1.ge.0.and.ktab2.ge.0) goto 140
      endif
 
C.. Recognize rank 0 diquark case
  100 irank=1
      kfdiq=max(kf1a,kf2a)
      if(kfdiq.gt.10.and.kfdiq.lt.10000) irank=0
 
C.. Join two flavours to meson or baryon. Test for popcorn.
      if(kf2a.gt.0)then
        mbary=0
        if(kfdiq.gt.10) then
          if(irank.eq.0.and.mstj(12).lt.5)
     &         call pjnmes(kfdiq)
          if(mstu(121).ne.0) return
          mbary=2
        endif
        kfqold=kf1a
        kfqver=kf2a
        goto 130
      endif
 
C.. Separate incoming flavours, curtain flavour consistency check
      kfin=kfl1
      kfqold=kf1a
      kfqpop=kf1a/10000
      if(kf1a.gt.10)then
         kfin=-kfl1
         kfl1a=mod(kf1a/1000,10)
         kfl1b=mod(kf1a/100,10)
         if(irank.eq.0)then
            qawt=1d0
            if(kfl1a.ge.3) qawt=parf(136+kfl1a/4)
            if(kfl1b.ge.3) qawt=qawt/parf(136+kfl1b/4)
            kfqpop=kfl1a+(kfl1b-kfl1a)*int(1d0/(qawt+1d0)+pjr(0))
         endif
         if(kfqpop.ne.kfl1b.and.kfqpop.ne.kfl1a) return
         kfqold=kfl1a+kfl1b-kfqpop
      endif
 
C...Meson/baryon choice. Set number of mesons if starting a popcorn
C...system.
  110 mbary=0
      if(kf1a.le.10.and.mstj(12).gt.0)then
         if(mstu(121).eq.-1.or.(1d0+parj(1))*pjr(0).gt.1d0)then
            mbary=1
            call pjnmes(0)
         endif
      elseif(kf1a.gt.10)then
         mbary=2
         if(irank.eq.0) call pjnmes(kf1a)
         if(mstu(121).gt.0) mbary=-1
      endif
 
C..x->H+q: Choose single vertex quark. Jump to form hadron.
      if(mbary.eq.0.or.mbary.eq.2)then
         kfqver=1+int((2d0+parj(2))*pjr(0))
         kfl3=isign(kfqver,-kfin)
         goto 130
      endif
 
C..x->H+qq: (IDW=proper PARF position for diquark weights)
      idw=160
C..   q->B+qq: Get curtain quark, different weights for q->B+B and
C..   q->B+M+...
      if(mbary.eq.1)then
         if(mstu(121).eq.0) idw=150
         sqwt=parf(idw+1)
         if(mstu(121).gt.0) sqwt=sqwt*parf(135)*parf(138)**mstu(121)
         kfqpop=1+int((2d0+sqwt)*pjr(0))
C..   Shift to s-curtain parameters if needed
         if(kfqpop.ge.3.and.mstj(12).ge.5)then
            parf(194)=parf(138)*parf(139)
            parf(193)=parj(8)+parj(9)
         endif
      endif
 
C.. x->H+qq: Get vertex quark
      if(mbary.eq.-1.and.mstj(12).ge.5)then
         idw=mstu(122)
         mstu(121)=mstu(121)-1
         if(idw.eq.170) then
            if(mstu(121).eq.0)then
               ipos=3*min(kfqpop-1,2)+min(kfqold-1,2)
            else
               ipos=3*3+3*max(0,min(kfqpop-2,1))+min(kfqold-1,2)
            endif
         else
            if(mstu(121).eq.0)then
               ipos=3*5+5*min(kfqpop-1,3)+min(kfqold-1,4)
            else
               ipos=3*5+5*4+min(kfqold-1,4)
            endif
         endif
         ipos=200+30*ipos+1
 
         imes=-1
         rmes=pjr(0)*parf(194)
  120    imes=imes+1
         rmes=rmes-parf(ipos+imes)
         if(imes.eq.30) then
            mstu(121)=-1
            kf=-111
            return
         endif
         if(rmes.gt.0d0) goto 120
         kmul=imes/5
         kfj=2*kmul+1
         if(kmul.eq.2) kfj=10003
         if(kmul.eq.3) kfj=10001
         if(kmul.eq.4) kfj=20003
         if(kmul.eq.5) kfj=5
         idiag=0
         kfqver=mod(imes,5)+1
         if(kfqver.ge.kfqold) kfqver=kfqver+1
         if(kfqver.gt.3)then
            idiag=kfqver-3
            kfqver=kfqold
         endif
      else
         if(mbary.eq.-1) idw=170
         sqwt=parf(idw+2)
         if(kfqpop.eq.3) sqwt=parf(idw+3)
         if(kfqpop.gt.3) sqwt=parf(idw+3)*(1d0/parf(idw+5)+1d0)/2d0
         kfqver=min(3,1+int((2d0+sqwt)*pjr(0)))
         if(kfqpop.lt.3.and.kfqver.lt.3)then
            kfqver=kfqpop
            if(pjr(0).gt.parf(idw+4)) kfqver=3-kfqpop
         endif
      endif
 
C..x->H+qq: form outgoing diquark with KFQPOP flag at 10000-pos
      kflds=3
      if(kfqpop.ne.kfqver)then
         swt=parf(idw+7)
         if(kfqver.eq.3) swt=parf(idw+6)
         if(kfqpop.ge.3) swt=parf(idw+5)
         if((1d0+swt)*pjr(0).lt.1d0) kflds=1
      endif
      kfdiq=900*max(kfqver,kfqpop)+100*(kfqver+kfqpop)+kflds
     &      +10000*kfqpop
      kfl3=isign(kfdiq,kfin)
 
C..x->M+y: flavour for meson.
  130 if(mbary.le.0)then
        kfla=max(kfqold,kfqver)
        kflb=min(kfqold,kfqver)
        kfs=isign(1,kfl1)
        if(kfla.ne.kfqold) kfs=-kfs
C... Form meson, with spin and flavour mixing for diagonal states.
        if(mbary.eq.-1.and.mstj(12).ge.5)then
           if(idiag.gt.0) kf=110*idiag+kfj
           if(idiag.eq.0) kf=(100*kfla+10*kflb+kfj)*kfs*(-1)**kfla
           return
        endif
        if(kfla.le.2) kmul=int(parj(11)+pjr(0))
        if(kfla.eq.3) kmul=int(parj(12)+pjr(0))
        if(kfla.ge.4) kmul=int(parj(13)+pjr(0))
        if(kmul.eq.0.and.parj(14).gt.0d0)then
          if(pjr(0).lt.parj(14)) kmul=2
        elseif(kmul.eq.1.and.parj(15)+parj(16)+parj(17).gt.0d0)then
          rmul=pjr(0)
          if(rmul.lt.parj(15)) kmul=3
          if(kmul.eq.1.and.rmul.lt.parj(15)+parj(16)) kmul=4
          if(kmul.eq.1.and.rmul.lt.parj(15)+parj(16)+parj(17)) kmul=5
        endif
        kfls=3
        if(kmul.eq.0.or.kmul.eq.3) kfls=1
        if(kmul.eq.5) kfls=5
        if(kfla.ne.kflb)then
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
 
C..Optional extra suppression of eta and eta'.
C..Allow shift to qq->B+q in old version (set IRANK to 0)
        if(kf.eq.221.or.kf.eq.331)then
           if(pjr(0).gt.parj(25+kf/300))then
              if(kf2a.gt.0) goto 130
              if(mstj(12).lt.4) irank=0
              goto 110
           endif
        endif
        mstu(121)=0
 
C.. x->B+y: Flavour for baryon
      else
        kfla=kfqver
        if(kf1a.le.10) kfla=kfqold
        kflb=mod(kfdiq/1000,10)
        kflc=mod(kfdiq/100,10)
        kflds=mod(kfdiq,10)
        kfld=max(kfla,kflb,kflc)
        kflf=min(kfla,kflb,kflc)
        kfle=kfla+kflb+kflc-kfld-kflf
 
C...  SU(6) factors for formation of baryon.
        kbary=3
        kdmax=5
        kflg=kflb
        if(kflb.ne.kflc)then
           kbary=2*kflds-1
           kdmax=1+kflds/2
           if(kflb.gt.2) kdmax=kdmax+2
        endif
        if(kfla.ne.kflb.and.kfla.ne.kflc)then
           kbary=kbary+1
           kflg=kfla
        endif
 
        su6max=parf(140+kdmax)
        su6dec=parj(18)
        su6s  =parf(146)
        if(mstj(12).ge.5.and.irank.eq.0) then
           su6max=1d0
           su6dec=1d0
           su6s  =1d0
        endif
        su6oct=parf(60+kbary)
        if(kflg.gt.max(kfla+kflb-kflg,2))then
           su6oct=su6oct*4*su6s/(3*su6s+1)
           if(kbary.eq.2) su6oct=parf(60+kbary)*4/(3*su6s+1)
        else
           if(kbary.eq.6) su6oct=su6oct*(3+su6s)/(3*su6s+1)
        endif
        su6wt=su6oct+su6dec*parf(70+kbary)
 
C..   SU(6) test. Old options enforce new baryon if q->B+qq is rejected.
        if(su6wt.lt.pjr(0)*su6max.and.kf2a.eq.0)then
           mstu(121)=0
           if(mstj(12).le.2.and.mbary.eq.1) mstu(121)=-1
           goto 110
        endif
 
C.. Form baryon. Distinguish Lambda- and Sigmalike baryons.
        ksig=1
        kfls=2
        if(su6wt*pjr(0).gt.su6oct) kfls=4
        if(kfls.eq.2.and.kfld.gt.kfle.and.kfle.gt.kflf)then
          ksig=kflds/3
          if(kfla.ne.kfld) ksig=int(3*su6s/(3*su6s+kflds**2)+pjr(0))
        endif
        kf=isign(1000*kfld+100*kfle+10*kflf+kfls,kfl1)
        if(ksig.eq.0) kf=isign(1000*kfld+100*kflf+10*kfle+kfls,kfl1)

cJAM++ baryon resonance
        if(pjr(0).le.parj(27)) then
          call jamkfres(kf,kf1)
          kf=kf1
        endif
cJAM--
      endif

      return
 
C...Use tabulated probabilities to select new flavour and hadron.
  140 if(ktab2.eq.0.and.mstj(12).le.0) then
        kt3l=1
        kt3u=6
      elseif(ktab2.eq.0.and.ktab1.ge.7.and.mstj(12).le.1) then
        kt3l=1
        kt3u=6
      elseif(ktab2.eq.0) then
        kt3l=1
        kt3u=22
      else
        kt3l=ktab2
        kt3u=ktab2
      endif
      rfl=0d0
      do 160 kts=0,2
        do 150 kt3=kt3l,kt3u
          rfl=rfl+parf(120+80*ktab1+25*kts+kt3)
  150   continue
  160 continue
      rfl=pjr(0)*rfl
      do 180 kts=0,2
        ktabs=kts
        do 170 kt3=kt3l,kt3u
          ktab3=kt3
          rfl=rfl-parf(120+80*ktab1+25*kts+kt3)
          if(rfl.le.0d0) goto 190
  170   continue
  180 continue
  190 continue
 
C...Reconstruct flavour of produced quark/diquark.
      if(ktab3.le.6) then
        kfl3a=ktab3
        kfl3b=0
        kfl3=isign(kfl3a,kfl1*(2*ktab1-13))
      else
        kfl3a=1
        if(ktab3.ge.8) kfl3a=2
        if(ktab3.ge.11) kfl3a=3
        if(ktab3.ge.16) kfl3a=4
        kfl3b=(ktab3-6-kfl3a*(kfl3a-2))/2
        kfl3=1000*kfl3a+100*kfl3b+1
        if(kfl3a.eq.kfl3b.or.ktab3.ne.6+kfl3a*(kfl3a-2)+2*kfl3b) kfl3=
     &  kfl3+2
        kfl3=isign(kfl3,kfl1*(13-2*ktab1))
      endif
 
C...Reconstruct meson code.
      if(kfl3a.eq.kfl1a.and.kfl3b.eq.kfl1b.and.(kfl3a.le.3.or.
     &kfl3b.ne.0)) then
        rfl=pjr(0)*(parf(143+80*ktab1+25*ktabs)+parf(144+80*ktab1+
     &  25*ktabs)+parf(145+80*ktab1+25*ktabs))
        kf=110+2*ktabs+1
        if(rfl.gt.parf(143+80*ktab1+25*ktabs)) kf=220+2*ktabs+1
        if(rfl.gt.parf(143+80*ktab1+25*ktabs)+parf(144+80*ktab1+
     &  25*ktabs)) kf=330+2*ktabs+1
      elseif(ktab1.le.6.and.ktab3.le.6) then
        kfla=max(ktab1,ktab3)
        kflb=min(ktab1,ktab3)
        kfs=isign(1,kfl1)
        if(kfla.ne.kf1a) kfs=-kfs
        kf=(100*kfla+10*kflb+2*ktabs+1)*kfs*(-1)**kfla
      elseif(ktab1.ge.7.and.ktab3.ge.7) then
        kfs=isign(1,kfl1)
        if(kfl1a.eq.kfl3a) then
          kfla=max(kfl1b,kfl3b)
          kflb=min(kfl1b,kfl3b)
          if(kfla.ne.kfl1b) kfs=-kfs
        elseif(kfl1a.eq.kfl3b) then
          kfla=kfl3a
          kflb=kfl1b
          kfs=-kfs
        elseif(kfl1b.eq.kfl3a) then
          kfla=kfl1a
          kflb=kfl3b
        elseif(kfl1b.eq.kfl3b) then
          kfla=max(kfl1a,kfl3a)
          kflb=min(kfl1a,kfl3a)
          if(kfla.ne.kfl1a) kfs=-kfs
        else
          call pjerrm(2,'(PYKFDI:) no matching flavours for qq -> qq')
          goto 100
        endif
        kf=(100*kfla+10*kflb+2*ktabs+1)*kfs*(-1)**kfla
 
C...Reconstruct baryon code.
      else
        if(ktab1.ge.7) then
          kfla=kfl3a
          kflb=kfl1a
          kflc=kfl1b
        else
          kfla=kfl1a
          kflb=kfl3a
          kflc=kfl3b
        endif
        kfld=max(kfla,kflb,kflc)
        kflf=min(kfla,kflb,kflc)
        kfle=kfla+kflb+kflc-kfld-kflf
        if(ktabs.eq.0) kf=isign(1000*kfld+100*kflf+10*kfle+2,kfl1)
        if(ktabs.ge.1) kf=isign(1000*kfld+100*kfle+10*kflf+2*ktabs,kfl1)

cJAM++ baryon resonance
        if(pjr(0).le.parj(27)) then
          call jamkfres(kf,kf1)
          kf=kf1
        endif
cJAM--
      endif
 
C...Check that constructed flavour code is an allowed one.
      if(kfl2.ne.0) kfl3=0
      kc=jamcomp(kf)
      if(kc.eq.0) then
        call pjerrm(2,'(PYKFDI:) user-defined flavour probabilities '//
     &  'failed')
        goto 100
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYNMES
C...Generates number of popcorn mesons and stores some relevant
C...parameters.
 
      subroutine pjnmes(kfdiq)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jydat1/,/jydat2/
 
      mstu(121)=0
      if(mstj(12).lt.2) return
 
C..Old version: Get 1 or 0 popcorn mesons
      if(mstj(12).lt.5)then
         popwt=parf(131)
         if(kfdiq.ne.0) then
            kfdiqa=iabs(kfdiq)
            kfa=mod(kfdiqa/1000,10)
            kfb=mod(kfdiqa/100,10)
            kfs=mod(kfdiqa,10)
            popwt=parf(132)
            if(kfa.eq.3) popwt=parf(133)
            if(kfb.eq.3) popwt=parf(134)
            if(kfs.eq.1) popwt=popwt*sqrt(parj(4))
         endif
         mstu(121)=int(popwt/(1d0+popwt)+pjr(0))
         return
      endif
 
C..New version: Store popcorn- or rank 0 diquark parameters
      mstu(122)=170
      parf(193)=parj(8)
      parf(194)=parf(139)
      if(kfdiq.ne.0) then
         mstu(122)=180
         parf(193)=parj(10)
         parf(194)=parf(140)
      endif
      if(parf(194).lt.1d-5.or.parf(194).gt.1d0-1d-5) then
         if(parf(194).gt.1d0-1d-5) call pjerrm(9,
     &        '(PYNMES:) Neglecting too large popcorn possibility')
         return
      endif
 
C..New version: Get number of popcorn mesons
  100 rtst=pjr(0)
      mstu(121)=-1
  110 mstu(121)=mstu(121)+1
      rtst=rtst/parf(194)
      if(rtst.lt.1d0) goto 110
      if(kfdiq.eq.0.and.pjr(0)*(2d0+parf(135)).gt.
     &     (2d0+parf(135)*parf(138)**mstu(121))) goto 100
      return
      end
 
C*********************************************************************
 
C...PYKFIN
C...Precalculates a set of diquark and popcorn weights.
C.. (Results stored in order SU0,US0,SS1,UU1,SU1,US1,UD1)
 
      subroutine pjkfin
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jydat1/,/jydat2/
 
      dimension su6(12),su6m(7)
 
      mstu(123)=1
C..Curtain tunneling factor T(D,q)/T(ud0,u).
      if(mstj(12).ge.5) then
         pmud0=pjmass(2101)
         pmud1=pjmass(2103)-pmud0
         pmus0=pjmass(3201)-pmud0
         pmus1=pjmass(3203)-pmus0-pmud0
         pmss1=pjmass(3303)-pmus0-pmud0
         parf(151)=exp(-(parj(9)+parj(8))*pmus0-parj(9)*parf(191))
         parf(152)=exp(-parj(8)*pmus0)
         parf(153)=exp(-(parj(9)+parj(8))*pmss1)*parf(151)
         parf(154)=exp(-parj(8)*pmud1)
         parf(155)=exp(-(parj(9)+parj(8))*pmus1)*parf(151)
         parf(156)=exp(-parj(8)*pmus1)*parf(152)
         parf(157)=parf(154)
      else
         par2m=sqrt(parj(2))
         par3m=sqrt(parj(3))
         par4m=sqrt(parj(4))
         parf(151)=par2m*par3m
         parf(152)=par3m
         parf(153)=par2m*parj(3)*par4m
         parf(154)=par4m
         parf(155)=par4m*parf(151)
         parf(156)=par4m*parf(152)
         parf(157)=par4m
      endif
 
C.. Total tunneling factor tau(D,q)=T*vertex*spin.
      parf(161)=parf(151)
      parf(162)=parj(2)*parf(152)
      parf(163)=parj(2)*6d0*parf(153)
      parf(164)=6d0*parf(154)
      parf(165)=3d0*parf(155)
      parf(166)=parj(2)*3d0*parf(156)
      parf(167)=3d0*parf(157)
 
      do 100 i=1,7
         parf(150+i)=parf(150+i)*parf(160+i)
  100 continue
 
C..Modified SU(6) factors.
      parf(146)=1d0
      if(mstj(12).ge.5) parf(146)=3d0*parj(18)/(2d0*parj(18)+1d0)
      if(parj(18).lt.1d0-1d-5.and.mstj(12).lt.5) call pjerrm(9,
     &     '(PYKFIN:) PARJ(18)<1 combined with 0<MSTJ(12)<5 option')
      do 110 i=1,6
         su6(i)=parf(60+i)
         su6(6+i)=su6(i)*4*parf(146)/(3*parf(146)+1)
  110 continue
      su6(8)=su6(2)*4/(3*parf(146)+1)
      su6(6)=su6(6)*(3+parf(146))/(3*parf(146)+1)
      do 120 i=1,6
         su6(i)=su6(i)+parj(18)*parf(70+i)
         su6(6+i)=su6(6+i)+parj(18)*parf(70+i)
  120 continue
 
C..Total diquark quark*SU(6).
      pud0=(2d0*su6(1)+parj(2)*su6(8))
      parf(171)=(su6(7)+su6(2)+parj(2)*su6(1))/pud0
      parf(172)=parf(171)
      parf(173)=(2d0*su6(4)+parj(2)*su6(3))/pud0
      parf(174)=(su6(3)+su6(4)+parj(2)*su6(10))/pud0
      parf(175)=(su6(11)+su6(6)+parj(2)*su6(5))/pud0
      parf(176)=parf(175)
      parf(177)=(2d0*su6(5)+parj(2)*su6(12))/pud0
 
C..SU(6)max         q       q'     s,c,b
      su6mud =max(su6(1) ,       su6(8) )
      su6m(7)=max(su6(5) ,       su6(12))
      su6m(1)=max(su6(7) ,su6(2),su6mud )
      su6m(4)=max(su6(3) ,su6(4),su6(10))
      su6m(5)=max(su6(11),su6(6),su6m(7))
      su6m(2)=su6m(1)
      su6m(3)=su6m(4)
      su6m(6)=su6m(5)
 
      if(mstj(12).ge.5)then
C..New version: tau for rank 0 diquark.
         parf(181)=exp(-parj(10)*pmus0)
         parf(182)=parj(2)*parf(181)
         parf(183)=6d0*parj(2)*exp(-parj(10)*pmss1)*parf(181)
         parf(184)=3d0*exp(-parj(10)*pmud1)
         parf(185)=3d0*exp(-parj(10)*pmus1)*parf(181)
         parf(186)=parj(2)*parf(185)
         parf(187)=2d0*parf(184)
 
C..New version: s/u curtain ratios.
         wu=1d0+parf(167)+parf(162)+parf(166)+parf(164)
         parf(135)=(2d0*(parf(161)+parf(165))+parf(163))/wu
         wu=1d0+parf(187)+parf(182)+parf(186)+parf(184)
         parf(136)=(2d0*(parf(181)+parf(185))+parf(183))/wu
         parf(137)=(parf(181)+parf(185))*
     &        (2d0+parf(183)/(2d0*parf(185)))/wu
      else
C..Old version: Shuffle PARJ(7) into tau
         parf(162)=parf(162)*parj(7)
         parf(163)=parf(163)*parj(7)
         parf(166)=parf(166)*parj(7)
 
C..Old version: s/u curtain ratios.
         wu=1d0+parf(167)+parf(162)+parf(166)+parf(164)
         parf(135)=(2d0*(parf(161)+parf(165))+parf(163))/wu
         parf(136)=parf(135)*parj(6)*parf(161)/parf(162)
         parf(137)=(1d0+parf(167))*(2d0+parf(162))/wu
      endif
 
C..Combine SU(6), SU(6)max, tau and T into proper products
      do 140 i=1,7
         parf(180+i)=parf(180+i)*parf(170+i)
         parf(170+i)=parf(170+i)*parf(160+i)
         parf(160+i)=parf(160+i)*su6m(i)/su6mud
         parf(150+i)=parf(150+i)*su6m(i)/su6mud
  140 continue
 
C..Store SU(6)max, in order UD0,UD1,US0,US1,QQ1
      parf(141)=su6mud
      parf(142)=su6m(7)
      parf(143)=su6m(1)
      parf(144)=su6m(5)
      parf(145)=su6m(3)
 
      if(mstj(12).lt.5)then
C.. Old version: Resulting popcorn weights.
         parf(138)=parj(6)
         ws=parf(135)*parf(138)
         wq=wu*parj(5)/3d0
         parf(132)=wq*parf(167)/parf(157)
         parf(133)=wq*(parf(166)/parf(156)+ws*parf(165)/parf(155))/2d0
         parf(134)=wq*ws*parf(163)/parf(153)
         parf(131)=wq*((1d0+parf(167))*(1d0+parf(162)+ws*parf(161))+
     &     parf(164)+ws*parf(163)/2d0)/
     &    ((1d0+parf(157))*(1d0+2d0*parf(152))+parf(154)+parf(153)/2d0)
      else
C..New version: Store weights for popcorn mesons,
C..get prel. popcorn weights.
         do 150 ipos=201,1400
            parf(ipos)=0d0
  150    continue
         do 160 i=138,140
            parf(i)=0d0
  160    continue
         ipos=200
         parf(193)=parj(8)
         do 240 mr=170,180,10
           if(mr.eq.180) parf(193)=parj(10)
           sqwt=2d0*(parf(mr+2)+parf(mr+6))/(1d0+parf(mr+7)+parf(mr+4))
           qqwt=parf(mr+4)/(1d0+parf(mr+7)+parf(mr+4))
           do 230 nmes=0,1
             if(nmes.eq.1) sqwt=parj(2)
             do 220 kfqpop=1,4
               if(mr.eq.170.and.kfqpop.gt.3) goto 220
               if(nmes.eq.0.and.kfqpop.ge.3)then
                  sqwt=parf(mr+3)/(parf(mr+1)+parf(mr+5))
                  qqwt=0.5d0
                  if(mr.eq.170) parf(193)=parj(8)+parj(9)
                  if(kfqpop.eq.4) sqwt=sqwt*(1d0/parf(185)+1d0)/2d0
               endif
               do 210 kfqold =1,5
                  if(mr.eq.170.and.kfqold.gt.3) goto 210
                  if(mr*nmes.eq.170.and.kfqpop.eq.1) goto 210
                  if(mr*nmes.eq.180.and.kfqpop.ne.1) goto 210
                  wttot=0d0
                  wtfail=0d0
      do 190 kmul=0,5
         pjwt=parj(12+kmul)
         if(kmul.eq.0) pjwt=1d0-parj(14)
         if(kmul.eq.1) pjwt=1d0-parj(15)-parj(16)-parj(17)
         if(pjwt.le.0d0) goto 190
         if(pjwt.gt.1d0) pjwt=1d0
         imes=5*kmul
         imix=2*kfqold+10*kmul
         kfj=2*kmul+1
         if(kmul.eq.2) kfj=10003
         if(kmul.eq.3) kfj=10001
         if(kmul.eq.4) kfj=20003
         if(kmul.eq.5) kfj=5
         do 180 kfqver =1,3
            kfla=max(kfqold,kfqver)
            kflb=min(kfqold,kfqver)
            swt=parj(11+kfla/3+kfla/4)
            if(kmul.eq.0.or.kmul.eq.2) swt=1d0-swt
            swt=swt*pjwt
            qwt=sqwt/(2d0+sqwt)
            if(kfqver.lt.3)then
               if(kfqver.eq.kfqpop) qwt=(1d0-qwt)*qqwt
               if(kfqver.ne.kfqpop) qwt=(1d0-qwt)*(1d0-qqwt)
            endif
            if(kfqver.ne.kfqold)then
               imes=imes+1
               kfm=100*kfla+10*kflb+kfj
               pmm=pmas(jamcomp(kfm),1)-pmas(jamcomp(kfm),3)
               parf(ipos+imes)=qwt*swt*exp(-parf(193)*pmm)
               wttot=wttot+parf(ipos+imes)
            else
               do 170 id=3,5
                  if(id.eq.3) dwt=1d0-parf(imix-1)
                  if(id.eq.4) dwt=parf(imix-1)-parf(imix)
                  if(id.eq.5) dwt=parf(imix)
                  kfm=110*(id-2)+kfj
                  pmm=pmas(jamcomp(kfm),1)-pmas(jamcomp(kfm),3)
                  parf(ipos+5*kmul+id)=qwt*swt*dwt*exp(-parf(193)*pmm)
                  if(kmul.eq.0.and.id.gt.3) then
                     wtfail=wtfail+qwt*swt*dwt*(1d0-parj(21+id))
                     parf(ipos+5*kmul+id)=
     &                    parf(ipos+5*kmul+id)*parj(21+id)
                  endif
                  wttot=wttot+parf(ipos+5*kmul+id)
  170          continue
            endif
  180    continue
  190 continue
                  do 200 imes=1,30
                     parf(ipos+imes)=parf(ipos+imes)/(1d0-wtfail)
  200             continue
                  if(mr.eq.180) parf(140)=
     &                 max(parf(140),wttot/(1d0-wtfail))
                  if(mr.eq.170) parf(139-kfqpop/3)=
     &                 max(parf(139-kfqpop/3),wttot/(1d0-wtfail))
                  ipos=ipos+30
  210           continue
  220         continue
  230       continue
  240    continue
         if(parf(139).gt.1d-10) parf(138)=parf(138)/parf(139)
         mstu(121)=0
 
         parf(186)=parf(186)/parf(182)
         parf(185)=parf(185)/parf(181)
      endif
 
C..Recombine diquark weights to flavour and spin ratios
      do 250 i=150,170,10
         wswq=(2d0*(parf(i+1)+parf(i+5))+parf(i+3))/
     &        (1d0+parf(i+7)+parf(i+4)+parf(i+2)+parf(i+6))
         wsswsq=parf(i+3)/(parf(i+1)+parf(i+5))
         wqswqq=2d0*(parf(i+2)+parf(i+6))/(1d0+parf(i+7)+parf(i+4))
         wuuwqq=parf(i+4)/(1d0+parf(i+7)+parf(i+4))
         parf(i+5)=parf(i+5)/parf(i+1)
         parf(i+6)=parf(i+6)/parf(i+2)
         parf(i+1)=wswq
         parf(i+2)=wqswqq
         parf(i+3)=wsswsq
         parf(i+4)=wuuwqq
  250 continue

      return
      end
 
C*********************************************************************
 
C...PYPTDI
C...Generates transverse momentum according to a Gaussian.
 
      subroutine pjptdi(kfl,px,py)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      save /jydat1/
 
C...Generate p_T and azimuthal angle, gives p_x and p_y.
      kfla=iabs(kfl)
      pt=parj(21)*sqrt(-log(max(1d-10,pjr(0))))
      if(parj(23).gt.pjr(0)) pt=parj(24)*pt
      if(mstj(91).eq.1) pt=parj(22)*pt
      if(kfla.eq.0.and.mstj(13).le.0) pt=0d0
      phi=paru(2)*pjr(0)
      px=pt*cos(phi)
      py=pt*sin(phi)
 
      return
      end
 
C*********************************************************************
 
C...PYZDIS
C...Generates the longitudinal splitting variable z.
 
      subroutine pjzdis(kfl1,kfl2,pr,z)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jydat1/,/jydat2/
 
C...Check if heavy flavour fragmentation.
      kfla=iabs(kfl1)
      kflb=iabs(kfl2)
      kflh=kfla
      if(kfla.ge.10) kflh=mod(kfla/1000,10)
 
C...Lund symmetric scaling function: determine parameters of shape.
      if(mstj(11).eq.1.or.(mstj(11).eq.3.and.kflh.le.3).or.
     &mstj(11).ge.4) then
        fa=parj(41)
        if(mstj(91).eq.1) fa=parj(43)
        if(kflb.ge.10) fa=fa+parj(45)
        fbb=parj(42)
        if(mstj(91).eq.1) fbb=parj(44)
        fb=fbb*pr
        fc=1d0
        if(kfla.ge.10) fc=fc-parj(45)
        if(kflb.ge.10) fc=fc+parj(45)
        if(mstj(11).ge.4.and.kflh.ge.4.and.kflh.le.5) then
          fred=parj(46)
          if(mstj(11).eq.5.and.kflh.eq.5) fred=parj(47)
          fc=fc+fred*fbb*parf(100+kflh)**2
        elseif(mstj(11).ge.4.and.kflh.ge.6.and.kflh.le.8) then
          fred=parj(46)
          if(mstj(11).eq.5) fred=parj(48)
          fc=fc+fred*fbb*pmas(kflh,1)**2
        endif
        mc=1
        if(abs(fc-1d0).gt.0.01d0) mc=2
 
C...Determine position of maximum. Special cases for a = 0 or a = c.
        if(fa.lt.0.02d0) then
          ma=1
          zmax=1d0
          if(fc.gt.fb) zmax=fb/fc
        elseif(abs(fc-fa).lt.0.01d0) then
          ma=2
          zmax=fb/(fb+fc)
        else
          ma=3
          zmax=0.5d0*(fb+fc-sqrt((fb-fc)**2+4d0*fa*fb))/(fc-fa)
          if(zmax.gt.0.9999d0.and.fb.gt.100d0) zmax=min(zmax,1d0-fa/fb)
        endif
 
C...Subdivide z range if distribution very peaked near endpoint.
        mmax=2
!///////////////////
!        write(0,*) ' fc =', fc, ' mc=',mc, ' zmax=',zmax
!c/////////////
        if(zmax.lt.0.1d0) then
!///////////////////
!        write(0,*) ' A' 
!/////////////
          mmax=1
          zdiv=2.75d0*zmax
          if(mc.eq.1) then
            fint=1d0-log(zdiv)
          else

!///////////////////
!        write(0,*) ' B' 
!/////////////

            zdivc=zdiv**(1d0-fc)

!////////////
!            write(0,*) ' zdiv=',zdiv, ' zdivc =',zdivc
!////////////////
!!!!!       fint=1d0+(1d0-1d0/zdivc)/(fc-1d0) !! relaced by next/ KK
            if( fc == 1.d0 ) then
               fint = 1.
            else
               fint=1d0+(1d0-1d0/zdivc)/(fc-1d0)
            endif
!!!
          endif
        elseif(zmax.gt.0.85d0.and.fb.gt.1d0) then

!///////////////////
!        write(0,*) ' C' 
!        write(0,*) ' fb=',fb, 'zmax=', zmax
!////////////////
          mmax=3
          fscb=sqrt(4d0+(fc/fb)**2)
          zdiv=fscb-1d0/zmax-(fc/fb)*log(zmax*0.5d0*(fscb+fc/fb))
          if(ma.ge.2) zdiv=zdiv+(fa/fb)*log(1d0-zmax)
          zdiv=min(zmax,max(0d0,zdiv))
          fint=1d0+fb*(1d0-zdiv)
        endif

!///////////////////
!        write(0,*) ' D' 
!////////////////
 
C...Choice of z, preweighted for peaks at low or high z.
 100    continue
!///////////////////
!        write(0,*) ' E' 
!////////////////
        z=pjr(0)
!/////////////
!        write(0,*) '2 fc=',fc, 'zdivc=',zdivc, ' z=',z
!        write(0,*) 'zmax =',zmax, ' nmax=',nmax
!///////////////
        fpre=1d0
        if(mmax.eq.1) then
          if(fint*pjr(0).le.1d0) then
            z=zdiv*z
          elseif(mc.eq.1) then
             if( zdiv == 0.) then  !! KK 
                write(0,*) ' zdiv=',zdiv, ' z=', z
                fpre = 1.
             else
                z=zdiv**z
                fpre=zdiv/z
             endif
          else
!/////////////
!             write(0,*) 'X  fc =', fc, ' zdivc=', zdivc
!////////////
             if( fc == 1.0 ) then
                z = zdivc
             else
                z=(zdivc+z*(1d0-zdivc))**(1d0/(1d0-fc))
                fpre=(zdiv/z)**fc
             endif
          endif
        elseif(mmax.eq.3) then
!////////////
!           write(0,*) ' posE , zdiv', zdiv, ' z=',z
!/////////////////
          if(fint*pjr(0).le.1d0) then
            z=zdiv+log(z)/fb
            fpre=exp(fb*(z-zdiv))
          else
            z=zdiv+z*(1d0-zdiv)
          endif
        endif
!////////////////
!        write(0,*) ' escape ? zmax=', zmax
!////////////
C...Weighting according to correct formula.
        if(z.le.0d0.or.z.ge.1d0) goto 100
        fexp=fc*log(zmax/z)+fb*(1d0/zmax-1d0/z)
        if(ma.ge.2) fexp=fexp+fa*log((1d0-z)/(1d0-zmax))
        fval=exp(max(-50d0,min(50d0,fexp)))
        if(fval.lt.pjr(0)*fpre) goto 100
 
C...Generate z according to Field-Feynman, SLAC, (1-z)**c OR z**c.
      else
        fc=parj(50+max(1,kflh))
        if(mstj(91).eq.1) fc=parj(59)
  110   z=pjr(0)
        if(fc.ge.0d0.and.fc.le.1d0) then
          if(fc.gt.pjr(0)) z=1d0-z**(1d0/3d0)
        elseif(fc.gt.-1.and.fc.lt.0d0) then
          if(-4d0*fc*z*(1d0-z)**2.lt.pjr(0)*((1d0-z)**2-fc*z)**2)
     &    goto 110
        else
          if(fc.gt.0d0) z=1d0-z**(1d0/fc)
          if(fc.lt.0d0) z=z**(-1d0/fc)
        endif
      endif
!//////////
!      write(0,*) ' return'
!//////////////

      return
      end
 
C*********************************************************************
 
C...PYSHOW
C...Generates timelike parton showers from given partons.
 
      subroutine pjshow(ip1,ip2,qmax)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jyjets/,/jydat1/,/jydat2/
C...Local arrays.
      dimension pmth(5,50),ps(5),pma(4),pmsd(4),iep(4),ipa(4),
     &kfla(4),kfld(4),kfl(4),itry(4),isi(4),isl(4),dp(4),dpt(5,4),
     &ksh(0:40),kcii(2),niis(2),iiis(2,2),theiis(2,2),phiiis(2,2),
     &isii(2)
 
C...Initialization of cutoff masses etc.
      if(mstj(41).le.0.or.(mstj(41).eq.1.and.qmax.le.parj(82)).or.
     &qmax.le.min(parj(82),parj(83))) return
      do 100 ifl=0,40
        ksh(ifl)=0
  100 continue
      ksh(21)=1
      pmth(1,21)=pjmass(21)
      pmth(2,21)=sqrt(pmth(1,21)**2+0.25d0*parj(82)**2)
      pmth(3,21)=2d0*pmth(2,21)
      pmth(4,21)=pmth(3,21)
      pmth(5,21)=pmth(3,21)
      pmth(1,22)=pjmass(22)
      pmth(2,22)=sqrt(pmth(1,22)**2+0.25d0*parj(83)**2)
      pmth(3,22)=2d0*pmth(2,22)
      pmth(4,22)=pmth(3,22)
      pmth(5,22)=pmth(3,22)
      pmqth1=parj(82)
      if(mstj(41).ge.2) pmqth1=min(parj(82),parj(83))
      pmqth2=pmth(2,21)
      if(mstj(41).ge.2) pmqth2=min(pmth(2,21),pmth(2,22))
      do 110 ifl=1,8
        ksh(ifl)=1
        pmth(1,ifl)=pjmass(ifl)
        pmth(2,ifl)=sqrt(pmth(1,ifl)**2+0.25d0*pmqth1**2)
        pmth(3,ifl)=pmth(2,ifl)+pmqth2
        pmth(4,ifl)=sqrt(pmth(1,ifl)**2+0.25d0*parj(82)**2)+pmth(2,21)
        pmth(5,ifl)=sqrt(pmth(1,ifl)**2+0.25d0*parj(83)**2)+pmth(2,22)
  110 continue
      do 120 ifl=11,17,2
        if(mstj(41).ge.2) ksh(ifl)=1
        pmth(1,ifl)=pjmass(ifl)
        pmth(2,ifl)=sqrt(pmth(1,ifl)**2+0.25d0*parj(83)**2)
        pmth(3,ifl)=pmth(2,ifl)+pmth(2,22)
        pmth(4,ifl)=pmth(3,ifl)
        pmth(5,ifl)=pmth(3,ifl)
  120 continue
      pt2min=max(0.5d0*parj(82),1.1d0*parj(81))**2
      alams=parj(81)**2
      alfm=log(pt2min/alams)
 
C...Store positions of shower initiating partons.
      if(ip1.gt.0.and.ip1.le.min(n,mstu(4)-mstu(32)).and.ip2.eq.0) then
        npa=1
        ipa(1)=ip1
      elseif(min(ip1,ip2).gt.0.and.max(ip1,ip2).le.min(n,mstu(4)-
     &  mstu(32))) then
        npa=2
        ipa(1)=ip1
        ipa(2)=ip2
      elseif(ip1.gt.0.and.ip1.le.min(n,mstu(4)-mstu(32)).and.ip2.lt.0
     &  .and.ip2.ge.-3) then
        npa=iabs(ip2)
        do 130 i=1,npa
          ipa(i)=ip1+i-1
  130   continue
      else
        call pjerrm(12,
     &  '(PYSHOW:) failed to reconstruct showering system')
        if(mstu(21).ge.1) return
      endif
 
C...Check on phase space available for emission.
      irej=0
      do 140 j=1,5
        ps(j)=0d0
  140 continue
      pm=0d0
      do 160 i=1,npa
        kfla(i)=iabs(k(ipa(i),2))
        pma(i)=p(ipa(i),5)
C...Special cutoff masses for t, l, h with variable masses.
        ifla=kfla(i)
        if(kfla(i).ge.6.and.kfla(i).le.8) then
          ifla=37+kfla(i)+isign(2,k(ipa(i),2))
          pmth(1,ifla)=pma(i)
          pmth(2,ifla)=sqrt(pmth(1,ifla)**2+0.25d0*pmqth1**2)
          pmth(3,ifla)=pmth(2,ifla)+pmqth2
          pmth(4,ifla)=sqrt(pmth(1,ifla)**2+0.25d0*parj(82)**2)+
     &    pmth(2,21)
          pmth(5,ifla)=sqrt(pmth(1,ifla)**2+0.25d0*parj(83)**2)+
     &    pmth(2,22)
        endif
        if(kfla(i).le.40) then
          if(ksh(kfla(i)).eq.1) pma(i)=pmth(3,ifla)
        endif
        pm=pm+pma(i)
        if(kfla(i).gt.40) then
          irej=irej+1
        else
          if(ksh(kfla(i)).eq.0.or.pma(i).gt.qmax) irej=irej+1
        endif
        do 150 j=1,4
          ps(j)=ps(j)+p(ipa(i),j)
  150   continue
  160 continue
      if(irej.eq.npa) return
      ps(5)=sqrt(max(0d0,ps(4)**2-ps(1)**2-ps(2)**2-ps(3)**2))
      if(npa.eq.1) ps(5)=ps(4)
      if(ps(5).le.pm+pmqth1) return
 
C...Check if 3-jet matrix elements to be used.
      m3jc=0
      if(npa.eq.2.and.mstj(47).ge.1) then
        if(kfla(1).ge.1.and.kfla(1).le.8.and.kfla(2).ge.1.and.
     &  kfla(2).le.8) m3jc=1
        if((kfla(1).eq.11.or.kfla(1).eq.13.or.kfla(1).eq.15.or.
     &  kfla(1).eq.17).and.kfla(2).eq.kfla(1)) m3jc=1
        if((kfla(1).eq.11.or.kfla(1).eq.13.or.kfla(1).eq.15.or.
     &  kfla(1).eq.17).and.kfla(2).eq.kfla(1)+1) m3jc=1
        if((kfla(1).eq.12.or.kfla(1).eq.14.or.kfla(1).eq.16.or.
     &  kfla(1).eq.18).and.kfla(2).eq.kfla(1)-1) m3jc=1
        if(mstj(47).eq.2.or.mstj(47).eq.4) m3jc=1
        m3jcm=0
        if(m3jc.eq.1.and.mstj(47).ge.3.and.kfla(1).eq.kfla(2)) then
          m3jcm=1
          qme=(2d0*pmth(1,kfla(1))/ps(5))**2
        endif
      endif
 
C...Find if interference with initial state partons.
      miis=0
      if(mstj(50).ge.1.and.mstj(50).le.3.and.npa.eq.2) miis=mstj(50)
      if(miis.ne.0) then
        do 180 i=1,2
          kcii(i)=0
          kca=jamcomp(kfla(i))
          if(kca.ne.0) kcii(i)=kchg(kca,2)*isign(1,k(ipa(i),2))
          niis(i)=0
          if(kcii(i).ne.0) then
            do 170 j=1,2
              icsi=mod(k(ipa(i),3+j)/mstu(5),mstu(5))
              if(icsi.gt.0.and.icsi.ne.ipa(1).and.icsi.ne.ipa(2).and.
     &        (kcii(i).eq.(-1)**(j+1).or.kcii(i).eq.2)) then
                niis(i)=niis(i)+1
                iiis(i,niis(i))=icsi
              endif
  170       continue
          endif
  180   continue
        if(niis(1)+niis(2).eq.0) miis=0
      endif
 
C...Boost interfering initial partons to rest frame
C...and reconstruct their polar and azimuthal angles.
      if(miis.ne.0) then
        do 200 i=1,2
          do 190 j=1,5
            k(n+i,j)=k(ipa(i),j)
            p(n+i,j)=p(ipa(i),j)
            v(n+i,j)=0d0
  190     continue
  200   continue
        do 220 i=3,2+niis(1)
          do 210 j=1,5
            k(n+i,j)=k(iiis(1,i-2),j)
            p(n+i,j)=p(iiis(1,i-2),j)
            v(n+i,j)=0d0
  210     continue
  220   continue
        do 240 i=3+niis(1),2+niis(1)+niis(2)
          do 230 j=1,5
            k(n+i,j)=k(iiis(2,i-2-niis(1)),j)
            p(n+i,j)=p(iiis(2,i-2-niis(1)),j)
            v(n+i,j)=0d0
  230     continue
  240   continue
        call pjrobo(n+1,n+2+niis(1)+niis(2),0d0,0d0,-ps(1)/ps(4),
     &  -ps(2)/ps(4),-ps(3)/ps(4))
        phi=pjangl(p(n+1,1),p(n+1,2))
        call pjrobo(n+1,n+2+niis(1)+niis(2),0d0,-phi,0d0,0d0,0d0)
        the=pjangl(p(n+1,3),p(n+1,1))
        call pjrobo(n+1,n+2+niis(1)+niis(2),-the,0d0,0d0,0d0,0d0)
        do 250 i=3,2+niis(1)
          theiis(1,i-2)=pjangl(p(n+i,3),sqrt(p(n+i,1)**2+p(n+i,2)**2))
          phiiis(1,i-2)=pjangl(p(n+i,1),p(n+i,2))
  250   continue
        do 260 i=3+niis(1),2+niis(1)+niis(2)
          theiis(2,i-2-niis(1))=paru(1)-pjangl(p(n+i,3),
     &    sqrt(p(n+i,1)**2+p(n+i,2)**2))
          phiiis(2,i-2-niis(1))=pjangl(p(n+i,1),p(n+i,2))
  260   continue
      endif
 
C...Define imagined single initiator of shower for parton system.
      ns=n
      if(n.gt.mstu(4)-mstu(32)-5) then
        call pjerrm(11,'(PYSHOW:) no more memory left in PYJETS')
        if(mstu(21).ge.1) return
      endif
      if(npa.ge.2) then
        k(n+1,1)=11
        k(n+1,2)=21
        k(n+1,3)=0
        k(n+1,4)=0
        k(n+1,5)=0
        p(n+1,1)=0d0
        p(n+1,2)=0d0
        p(n+1,3)=0d0
        p(n+1,4)=ps(5)
        p(n+1,5)=ps(5)
        v(n+1,5)=ps(5)**2
        n=n+1
      endif
 
C...Loop over partons that may branch.
      nep=npa
      im=ns
      if(npa.eq.1) im=ns-1
  270 im=im+1
      if(n.gt.ns) then
        if(im.gt.n) goto 510
        kflm=iabs(k(im,2))
        if(kflm.gt.40) goto 270
        if(ksh(kflm).eq.0) goto 270
        iflm=kflm
        if(kflm.ge.6.and.kflm.le.8) iflm=37+kflm+isign(2,k(im,2))
        if(p(im,5).lt.pmth(2,iflm)) goto 270
        igm=k(im,3)
      else
        igm=-1
      endif
      if(n+nep.gt.mstu(4)-mstu(32)-5) then
        call pjerrm(11,'(PYSHOW:) no more memory left in PYJETS')
        if(mstu(21).ge.1) return
      endif
 
C...Position of aunt (sister to branching parton).
C...Origin and flavour of daughters.
      iau=0
      if(igm.gt.0) then
        if(k(im-1,3).eq.igm) iau=im-1
        if(n.ge.im+1.and.k(im+1,3).eq.igm) iau=im+1
      endif
      if(igm.ge.0) then
        k(im,4)=n+1
        do 280 i=1,nep
          k(n+i,3)=im
  280   continue
      else
        k(n+1,3)=ipa(1)
      endif
      if(igm.le.0) then
        do 290 i=1,nep
          k(n+i,2)=k(ipa(i),2)
  290   continue
      elseif(kflm.ne.21) then
        k(n+1,2)=k(im,2)
        k(n+2,2)=k(im,5)
      elseif(k(im,5).eq.21) then
        k(n+1,2)=21
        k(n+2,2)=21
      else
        k(n+1,2)=k(im,5)
        k(n+2,2)=-k(im,5)
      endif
 
C...Reset flags on daughers and tries made.
      do 300 ip=1,nep
        k(n+ip,1)=3
        k(n+ip,4)=0
        k(n+ip,5)=0
        kfld(ip)=iabs(k(n+ip,2))
        if(kchg(jamcomp(kfld(ip)),2).eq.0) k(n+ip,1)=1
        itry(ip)=0
        isl(ip)=0
        isi(ip)=0
        if(kfld(ip).le.40) then
          if(ksh(kfld(ip)).eq.1) isi(ip)=1
        endif
  300 continue
      islm=0
 
C...Maximum virtuality of daughters.
      if(igm.le.0) then
        do 310 i=1,npa
          if(npa.ge.3) p(n+i,4)=(ps(4)*p(ipa(i),4)-ps(1)*p(ipa(i),1)-
     &    ps(2)*p(ipa(i),2)-ps(3)*p(ipa(i),3))/ps(5)
          p(n+i,5)=min(qmax,ps(5))
          if(npa.ge.3) p(n+i,5)=min(p(n+i,5),p(n+i,4))
          if(isi(i).eq.0) p(n+i,5)=p(ipa(i),5)
  310   continue
      else
        if(mstj(43).le.2) pem=v(im,2)
        if(mstj(43).ge.3) pem=p(im,4)
        p(n+1,5)=min(p(im,5),v(im,1)*pem)
        p(n+2,5)=min(p(im,5),(1d0-v(im,1))*pem)
        if(k(n+2,2).eq.22) p(n+2,5)=pmth(1,22)
      endif
      do 320 i=1,nep
        pmsd(i)=p(n+i,5)
        if(isi(i).eq.1) then
          ifld=kfld(i)
          if(kfld(i).ge.6.and.kfld(i).le.8) ifld=37+kfld(i)+
     &    isign(2,k(n+i,2))
          if(p(n+i,5).le.pmth(3,ifld)) p(n+i,5)=pmth(1,ifld)
        endif
        v(n+i,5)=p(n+i,5)**2
  320 continue
 
C...Choose one of the daughters for evolution.
  330 inum=0
      if(nep.eq.1) inum=1
      do 340 i=1,nep
        if(inum.eq.0.and.isl(i).eq.1) inum=i
  340 continue
      do 350 i=1,nep
        if(inum.eq.0.and.itry(i).eq.0.and.isi(i).eq.1) then
          ifld=kfld(i)
          if(kfld(i).ge.6.and.kfld(i).le.8) ifld=37+kfld(i)+
     &    isign(2,k(n+i,2))
          if(p(n+i,5).ge.pmth(2,ifld)) inum=i
        endif
  350 continue
      if(inum.eq.0) then
        rmax=0d0
        do 360 i=1,nep
          if(isi(i).eq.1.and.pmsd(i).ge.pmqth2) then
            rpm=p(n+i,5)/pmsd(i)
            ifld=kfld(i)
            if(kfld(i).ge.6.and.kfld(i).le.8) ifld=37+kfld(i)+
     &      isign(2,k(n+i,2))
            if(rpm.gt.rmax.and.p(n+i,5).ge.pmth(2,ifld)) then
              rmax=rpm
              inum=i
            endif
          endif
  360   continue
      endif
 
C...Store information on choice of evolving daughter.
      inum=max(1,inum)
      iep(1)=n+inum
      do 370 i=2,nep
        iep(i)=iep(i-1)+1
        if(iep(i).gt.n+nep) iep(i)=n+1
  370 continue
      do 380 i=1,nep
        kfl(i)=iabs(k(iep(i),2))
  380 continue
      itry(inum)=itry(inum)+1
      if(itry(inum).gt.200) then
        call pjerrm(14,'(PYSHOW:) caught in infinite loop')
        if(mstu(21).ge.1) return
      endif
      z=0.5d0
      if(kfl(1).gt.40) goto 430
      if(ksh(kfl(1)).eq.0) goto 430
      ifl=kfl(1)
      if(kfl(1).ge.6.and.kfl(1).le.8) ifl=37+kfl(1)+
     &isign(2,k(iep(1),2))
      if(p(iep(1),5).lt.pmth(2,ifl)) goto 430
 
C...Select side for interference with initial state partons.
      if(miis.ge.1.and.iep(1).le.ns+3) then
        iii=iep(1)-ns-1
        isii(iii)=0
        if(iabs(kcii(iii)).eq.1.and.niis(iii).eq.1) then
          isii(iii)=1
        elseif(kcii(iii).eq.2.and.niis(iii).eq.1) then
          if(pjr(0).gt.0.5d0) isii(iii)=1
        elseif(kcii(iii).eq.2.and.niis(iii).eq.2) then
          isii(iii)=1
          if(pjr(0).gt.0.5d0) isii(iii)=2
        endif
      endif
 
C...Calculate allowed z range.
      if(nep.eq.1) then
        pmed=ps(4)
      elseif(igm.eq.0.or.mstj(43).le.2) then
        pmed=p(im,5)
      else
        if(inum.eq.1) pmed=v(im,1)*pem
        if(inum.eq.2) pmed=(1d0-v(im,1))*pem
      endif
      if(mod(mstj(43),2).eq.1) then
        zc=pmth(2,21)/pmed
        zce=pmth(2,22)/pmed
      else
        zc=0.5d0*(1d0-sqrt(max(0d0,1d0-(2d0*pmth(2,21)/pmed)**2)))
        if(zc.lt.1d-4) zc=(pmth(2,21)/pmed)**2
        zce=0.5d0*(1d0-sqrt(max(0d0,1d0-(2d0*pmth(2,22)/pmed)**2)))
        if(zce.lt.1d-4) zce=(pmth(2,22)/pmed)**2
      endif
      zc=min(zc,0.491d0)
      zce=min(zce,0.491d0)
      if((mstj(41).eq.1.and.zc.gt.0.49d0).or.(mstj(41).ge.2.and.
     &min(zc,zce).gt.0.49d0)) then
        p(iep(1),5)=pmth(1,ifl)
        v(iep(1),5)=p(iep(1),5)**2
        goto 430
      endif
 
C...Integral of Altarelli-Parisi z kernel for QCD.
      if(mstj(49).eq.0.and.kfl(1).eq.21) then
        fbr=6d0*log((1d0-zc)/zc)+mstj(45)*(0.5d0-zc)
      elseif(mstj(49).eq.0) then
        fbr=(8d0/3d0)*log((1d0-zc)/zc)
 
C...Integral of Altarelli-Parisi z kernel for scalar gluon.
      elseif(mstj(49).eq.1.and.kfl(1).eq.21) then
        fbr=(parj(87)+mstj(45)*parj(88))*(1d0-2d0*zc)
      elseif(mstj(49).eq.1) then
        fbr=(1d0-2d0*zc)/3d0
        if(igm.eq.0.and.m3jc.eq.1) fbr=4d0*fbr
 
C...Integral of Altarelli-Parisi z kernel for Abelian vector gluon.
      elseif(kfl(1).eq.21) then
        fbr=6d0*mstj(45)*(0.5d0-zc)
      else
        fbr=2d0*log((1d0-zc)/zc)
      endif
 
C...Reset QCD probability for lepton.
      if(kfl(1).ge.11.and.kfl(1).le.18) fbr=0d0
 
C...Integral of Altarelli-Parisi kernel for photon emission.
      if(mstj(41).ge.2.and.kfl(1).ge.1.and.kfl(1).le.18) then
        fbre=(kchg(kfl(1),1)/3d0)**2*2d0*log((1d0-zce)/zce)
        if(mstj(41).eq.10) fbre=parj(84)*fbre
      endif
 
C...Inner veto algorithm starts. Find maximum mass for evolution.
  390 pms=v(iep(1),5)
      if(igm.ge.0) then
        pm2=0d0
        do 400 i=2,nep
          pm=p(iep(i),5)
          if(kfl(i).le.40) then
            ifli=kfl(i)
            if(kfl(i).ge.6.and.kfl(i).le.8) ifli=37+kfl(i)+
     &      isign(2,k(iep(i),2))
            if(ksh(kfl(i)).eq.1) pm=pmth(2,ifli)
          endif
          pm2=pm2+pm
  400   continue
        pms=min(pms,(p(im,5)-pm2)**2)
      endif
 
C...Select mass for daughter in QCD evolution.
      b0=27d0/6d0
      do 410 iff=4,mstj(45)
        if(pms.gt.4d0*pmth(2,iff)**2) b0=(33d0-2d0*iff)/6d0
  410 continue
      if(fbr.lt.1d-3) then
        pmsqcd=0d0
      elseif(mstj(44).le.0) then
        pmsqcd=pms*exp(max(-50d0,log(pjr(0))*paru(2)/(paru(111)*fbr)))
      elseif(mstj(44).eq.1) then
        pmsqcd=4d0*alams*(0.25d0*pms/alams)**(pjr(0)**(b0/fbr))
      else
        pmsqcd=pms*exp(max(-50d0,alfm*b0*log(pjr(0))/fbr))
      endif
      if(zc.gt.0.49d0.or.pmsqcd.le.pmth(4,ifl)**2) pmsqcd=pmth(2,ifl)**2
      v(iep(1),5)=pmsqcd
      mce=1
 
C...Select mass for daughter in QED evolution.
      if(mstj(41).ge.2.and.kfl(1).ge.1.and.kfl(1).le.18) then
        pmsqed=pms*exp(max(-50d0,log(pjr(0))*paru(2)/(paru(101)*fbre)))
        if(zce.gt.0.49d0.or.pmsqed.le.pmth(5,ifl)**2) pmsqed=
     &  pmth(2,ifl)**2
        if(pmsqed.gt.pmsqcd) then
          v(iep(1),5)=pmsqed
          mce=2
        endif
      endif
 
C...Check whether daughter mass below cutoff.
      p(iep(1),5)=sqrt(v(iep(1),5))
      if(p(iep(1),5).le.pmth(3,ifl)) then
        p(iep(1),5)=pmth(1,ifl)
        v(iep(1),5)=p(iep(1),5)**2
        goto 430
      endif
 
C...Select z value of branching: q -> qgamma.
      if(mce.eq.2) then
        z=1d0-(1d0-zce)*(zce/(1d0-zce))**pjr(0)
        if(1d0+z**2.lt.2d0*pjr(0)) goto 390
        k(iep(1),5)=22
 
C...Select z value of branching: q -> qg, g -> gg, g -> qqbar.
      elseif(mstj(49).ne.1.and.kfl(1).ne.21) then
        z=1d0-(1d0-zc)*(zc/(1d0-zc))**pjr(0)
        if(1d0+z**2.lt.2d0*pjr(0)) goto 390
        k(iep(1),5)=21
      elseif(mstj(49).eq.0.and.mstj(45)*(0.5d0-zc).lt.pjr(0)*fbr) then
        z=(1d0-zc)*(zc/(1d0-zc))**pjr(0)
        if(pjr(0).gt.0.5d0) z=1d0-z
        if((1d0-z*(1d0-z))**2.lt.pjr(0)) goto 390
        k(iep(1),5)=21
      elseif(mstj(49).ne.1) then
        z=zc+(1d0-2d0*zc)*pjr(0)
        if(z**2+(1d0-z)**2.lt.pjr(0)) goto 390
        kflb=1+int(mstj(45)*pjr(0))
        pmq=4d0*pmth(2,kflb)**2/v(iep(1),5)
        if(pmq.ge.1d0) goto 390
        pmq0=4d0*pmth(2,21)**2/v(iep(1),5)
        if(mod(mstj(43),2).eq.0.and.(1d0+0.5d0*pmq)*sqrt(1d0-pmq).lt.
     &  pjr(0)*(1d0+0.5d0*pmq0)*sqrt(1d0-pmq0)) goto 390
        k(iep(1),5)=kflb
 
C...Ditto for scalar gluon model.
      elseif(kfl(1).ne.21) then
        z=1d0-sqrt(zc**2+pjr(0)*(1d0-2d0*zc))
        k(iep(1),5)=21
      elseif(pjr(0)*(parj(87)+mstj(45)*parj(88)).le.parj(87)) then
        z=zc+(1d0-2d0*zc)*pjr(0)
        k(iep(1),5)=21
      else
        z=zc+(1d0-2d0*zc)*pjr(0)
        kflb=1+int(mstj(45)*pjr(0))
        pmq=4d0*pmth(2,kflb)**2/v(iep(1),5)
        if(pmq.ge.1d0) goto 390
        k(iep(1),5)=kflb
      endif
      if(mce.eq.1.and.mstj(44).ge.2) then
        if(z*(1d0-z)*v(iep(1),5).lt.pt2min) goto 390
        if(alfm/log(v(iep(1),5)*z*(1d0-z)/alams).lt.pjr(0)) goto 390
      endif
 
C...Check if z consistent with chosen m.
      if(kfl(1).eq.21) then
        kflgd1=iabs(k(iep(1),5))
        kflgd2=kflgd1
      else
        kflgd1=kfl(1)
        kflgd2=iabs(k(iep(1),5))
      endif
      if(nep.eq.1) then
        ped=ps(4)
      elseif(nep.ge.3) then
        ped=p(iep(1),4)
      elseif(igm.eq.0.or.mstj(43).le.2) then
        ped=0.5d0*(v(im,5)+v(iep(1),5)-pm2**2)/p(im,5)
      else
        if(iep(1).eq.n+1) ped=v(im,1)*pem
        if(iep(1).eq.n+2) ped=(1d0-v(im,1))*pem
      endif
      if(mod(mstj(43),2).eq.1) then
        iflgd1=kflgd1
        if(kflgd1.ge.6.and.kflgd1.le.8) iflgd1=ifl
        pmqth3=0.5d0*parj(82)
        if(kflgd2.eq.22) pmqth3=0.5d0*parj(83)
        pmq1=(pmth(1,iflgd1)**2+pmqth3**2)/v(iep(1),5)
        pmq2=(pmth(1,kflgd2)**2+pmqth3**2)/v(iep(1),5)
        zd=sqrt(max(0d0,(1d0-v(iep(1),5)/ped**2)*((1d0-pmq1-pmq2)**2-
     &  4d0*pmq1*pmq2)))
        zh=1d0+pmq1-pmq2
      else
        zd=sqrt(max(0d0,1d0-v(iep(1),5)/ped**2))
        zh=1d0
      endif
      zl=0.5d0*(zh-zd)
      zu=0.5d0*(zh+zd)
      if(z.lt.zl.or.z.gt.zu) goto 390
      if(kfl(1).eq.21) v(iep(1),3)=log(zu*(1d0-zl)/max(1d-20,zl*
     &(1d0-zu)))
      if(kfl(1).ne.21) v(iep(1),3)=log((1d0-zl)/max(1d-10,1d0-zu))
 
C...Width suppression for q -> q + g.
      if(mstj(40).ne.0.and.kfl(1).ne.21) then
        if(igm.eq.0) then
          eglu=0.5d0*ps(5)*(1d0-z)*(1d0+v(iep(1),5)/v(ns+1,5))
        else
          eglu=pmed*(1d0-z)
        endif
        chi=parj(89)**2/(parj(89)**2+eglu**2)
        if(mstj(40).eq.1) then
          if(chi.lt.pjr(0)) goto 390
        elseif(mstj(40).eq.2) then
          if(1d0-chi.lt.pjr(0)) goto 390
        endif
      endif
 
C...Three-jet matrix element correction.
      if(igm.eq.0.and.m3jc.eq.1) then
        x1=z*(1d0+v(iep(1),5)/v(ns+1,5))
        x2=1d0-v(iep(1),5)/v(ns+1,5)
        x3=(1d0-x1)+(1d0-x2)
        if(mce.eq.2) then
          ki1=k(ipa(inum),2)
          ki2=k(ipa(3-inum),2)
          qf1=kchg(iabs(ki1),1)*isign(1,ki1)/3d0
          qf2=kchg(iabs(ki2),1)*isign(1,ki2)/3d0
          wshow=qf1**2*(1d0-x1)/x3*(1d0+(x1/(2d0-x2))**2)+
     &    qf2**2*(1d0-x2)/x3*(1d0+(x2/(2d0-x1))**2)
          wme=(qf1*(1d0-x1)/x3-qf2*(1d0-x2)/x3)**2*(x1**2+x2**2)
        elseif(mstj(49).ne.1) then
          wshow=1d0+(1d0-x1)/x3*(x1/(2d0-x2))**2+
     &    (1d0-x2)/x3*(x2/(2d0-x1))**2
          wme=x1**2+x2**2
          if(m3jcm.eq.1) wme=wme-qme*x3-0.5d0*qme**2-
     &    (0.5d0*qme+0.25d0*qme**2)*((1d0-x2)/max(1d-7,1d0-x1)+
     &    (1d0-x1)/max(1d-7,1d0-x2))
        else
          wshow=4d0*x3*((1d0-x1)/(2d0-x2)**2+(1d0-x2)/(2d0-x1)**2)
          wme=x3**2
          if(mstj(102).ge.2) wme=x3**2-2d0*(1d0+x3)*(1d0-x1)*(1d0-x2)*
     &    parj(171)
        endif
        if(wme.lt.pjr(0)*wshow) goto 390
 
C...Impose angular ordering by rejection of nonordered emission.
      elseif(mce.eq.1.and.igm.gt.0.and.mstj(42).ge.2) then
        maom=1
        zm=v(im,1)
        if(iep(1).eq.n+2) zm=1d0-v(im,1)
        the2id=z*(1d0-z)*(zm*p(im,4))**2/v(iep(1),5)
        iaom=im
  420   if(k(iaom,5).eq.22) then
          iaom=k(iaom,3)
          if(k(iaom,3).le.ns) maom=0
          if(maom.eq.1) goto 420
        endif
        if(maom.eq.1) then
          the2im=v(iaom,1)*(1d0-v(iaom,1))*p(iaom,4)**2/v(iaom,5)
          if(the2id.lt.the2im) goto 390
        endif
      endif
 
C...Impose user-defined maximum angle at first branching.
      if(mstj(48).eq.1) then
        if(nep.eq.1.and.im.eq.ns) then
          the2id=z*(1d0-z)*ps(4)**2/v(iep(1),5)
          if(the2id.lt.1d0/parj(85)**2) goto 390
        elseif(nep.eq.2.and.iep(1).eq.ns+2) then
          the2id=z*(1d0-z)*(0.5d0*p(im,4))**2/v(iep(1),5)
          if(the2id.lt.1d0/parj(85)**2) goto 390
        elseif(nep.eq.2.and.iep(1).eq.ns+3) then
          the2id=z*(1d0-z)*(0.5d0*p(im,4))**2/v(iep(1),5)
          if(the2id.lt.1d0/parj(86)**2) goto 390
        endif
      endif
 
C...Impose angular constraint in first branching from interference
C...with initial state partons.
      if(miis.ge.2.and.iep(1).le.ns+3) then
        the2d=max((1d0-z)/z,z/(1d0-z))*v(iep(1),5)/(0.5d0*p(im,4))**2
        if(iep(1).eq.ns+2.and.isii(1).ge.1) then
          if(the2d.gt.theiis(1,isii(1))**2) goto 390
        elseif(iep(1).eq.ns+3.and.isii(2).ge.1) then
          if(the2d.gt.theiis(2,isii(2))**2) goto 390
        endif
      endif
 
C...End of inner veto algorithm. Check if only one leg evolved so far.
  430 v(iep(1),1)=z
      isl(1)=0
      isl(2)=0
      if(nep.eq.1) goto 460
      if(nep.eq.2.and.p(iep(1),5)+p(iep(2),5).ge.p(im,5)) goto 330
      do 440 i=1,nep
        if(itry(i).eq.0.and.kfld(i).le.40) then
          if(ksh(kfld(i)).eq.1) then
            ifld=kfld(i)
            if(kfld(i).ge.6.and.kfld(i).le.8) ifld=37+kfld(i)+
     &      isign(2,k(n+i,2))
            if(p(n+i,5).ge.pmth(2,ifld)) goto 330
          endif
        endif
  440 continue
 
C...Check if chosen multiplet m1,m2,z1,z2 is physical.
      if(nep.eq.3) then
        pa1s=(p(n+1,4)+p(n+1,5))*(p(n+1,4)-p(n+1,5))
        pa2s=(p(n+2,4)+p(n+2,5))*(p(n+2,4)-p(n+2,5))
        pa3s=(p(n+3,4)+p(n+3,5))*(p(n+3,4)-p(n+3,5))
        pts=0.25d0*(2d0*pa1s*pa2s+2d0*pa1s*pa3s+2d0*pa2s*pa3s-
     &  pa1s**2-pa2s**2-pa3s**2)/pa1s
        if(pts.le.0d0) goto 330
      elseif(igm.eq.0.or.mstj(43).le.2.or.mod(mstj(43),2).eq.0) then
        do 450 i1=n+1,n+2
          kflda=iabs(k(i1,2))
          if(kflda.gt.40) goto 450
          if(ksh(kflda).eq.0) goto 450
          iflda=kflda
          if(kflda.ge.6.and.kflda.le.8) iflda=37+kflda+
     &    isign(2,k(i1,2))
          if(p(i1,5).lt.pmth(2,iflda)) goto 450
          if(kflda.eq.21) then
            kflgd1=iabs(k(i1,5))
            kflgd2=kflgd1
          else
            kflgd1=kflda
            kflgd2=iabs(k(i1,5))
          endif
          i2=2*n+3-i1
          if(igm.eq.0.or.mstj(43).le.2) then
            ped=0.5d0*(v(im,5)+v(i1,5)-v(i2,5))/p(im,5)
          else
            if(i1.eq.n+1) zm=v(im,1)
            if(i1.eq.n+2) zm=1d0-v(im,1)
            pml=sqrt((v(im,5)-v(n+1,5)-v(n+2,5))**2-
     &      4d0*v(n+1,5)*v(n+2,5))
            ped=pem*(0.5d0*(v(im,5)-pml+v(i1,5)-v(i2,5))+pml*zm)/v(im,5)
          endif
          if(mod(mstj(43),2).eq.1) then
            pmqth3=0.5d0*parj(82)
            if(kflgd2.eq.22) pmqth3=0.5d0*parj(83)
            iflgd1=kflgd1
            if(kflgd1.ge.6.and.kflgd1.le.8) iflgd1=iflda
            pmq1=(pmth(1,iflgd1)**2+pmqth3**2)/v(i1,5)
            pmq2=(pmth(1,kflgd2)**2+pmqth3**2)/v(i1,5)
            zd=sqrt(max(0d0,(1d0-v(i1,5)/ped**2)*((1d0-pmq1-pmq2)**2-
     &      4d0*pmq1*pmq2)))
            zh=1d0+pmq1-pmq2
          else
            zd=sqrt(max(0d0,1d0-v(i1,5)/ped**2))
            zh=1d0
          endif
          zl=0.5d0*(zh-zd)
          zu=0.5d0*(zh+zd)
          if(i1.eq.n+1.and.(v(i1,1).lt.zl.or.v(i1,1).gt.zu)) isl(1)=1
          if(i1.eq.n+2.and.(v(i1,1).lt.zl.or.v(i1,1).gt.zu)) isl(2)=1
          if(kflda.eq.21) v(i1,4)=log(zu*(1d0-zl)/max(1d-20,
     &    zl*(1d0-zu)))
          if(kflda.ne.21) v(i1,4)=log((1d0-zl)/max(1d-10,1d0-zu))
  450   continue
        if(isl(1).eq.1.and.isl(2).eq.1.and.islm.ne.0) then
          isl(3-islm)=0
          islm=3-islm
        elseif(isl(1).eq.1.and.isl(2).eq.1) then
          zdr1=max(0d0,v(n+1,3)/max(1d-6,v(n+1,4))-1d0)
          zdr2=max(0d0,v(n+2,3)/max(1d-6,v(n+2,4))-1d0)
          if(zdr2.gt.pjr(0)*(zdr1+zdr2)) isl(1)=0
          if(isl(1).eq.1) isl(2)=0
          if(isl(1).eq.0) islm=1
          if(isl(2).eq.0) islm=2
        endif
        if(isl(1).eq.1.or.isl(2).eq.1) goto 330
      endif
      ifld1=kfld(1)
      if(kfld(1).ge.6.and.kfld(1).le.8) ifld1=37+kfld(1)+
     &isign(2,k(n+1,2))
      ifld2=kfld(2)
      if(kfld(2).ge.6.and.kfld(2).le.8) ifld2=37+kfld(2)+
     &isign(2,k(n+2,2))
      if(igm.gt.0.and.mod(mstj(43),2).eq.1.and.(p(n+1,5).ge.
     &pmth(2,ifld1).or.p(n+2,5).ge.pmth(2,ifld2))) then
        pmq1=v(n+1,5)/v(im,5)
        pmq2=v(n+2,5)/v(im,5)
        zd=sqrt(max(0d0,(1d0-v(im,5)/pem**2)*((1d0-pmq1-pmq2)**2-
     &  4d0*pmq1*pmq2)))
        zh=1d0+pmq1-pmq2
        zl=0.5d0*(zh-zd)
        zu=0.5d0*(zh+zd)
        if(v(im,1).lt.zl.or.v(im,1).gt.zu) goto 330
      endif
 
C...Accepted branch. Construct four-momentum for initial partons.
  460 mazip=0
      mazic=0
      if(nep.eq.1) then
        p(n+1,1)=0d0
        p(n+1,2)=0d0
        p(n+1,3)=sqrt(max(0d0,(p(ipa(1),4)+p(n+1,5))*(p(ipa(1),4)-
     &  p(n+1,5))))
        p(n+1,4)=p(ipa(1),4)
        v(n+1,2)=p(n+1,4)
      elseif(igm.eq.0.and.nep.eq.2) then
        ped1=0.5d0*(v(im,5)+v(n+1,5)-v(n+2,5))/p(im,5)
        p(n+1,1)=0d0
        p(n+1,2)=0d0
        p(n+1,3)=sqrt(max(0d0,(ped1+p(n+1,5))*(ped1-p(n+1,5))))
        p(n+1,4)=ped1
        p(n+2,1)=0d0
        p(n+2,2)=0d0
        p(n+2,3)=-p(n+1,3)
        p(n+2,4)=p(im,5)-ped1
        v(n+1,2)=p(n+1,4)
        v(n+2,2)=p(n+2,4)
      elseif(nep.eq.3) then
        p(n+1,1)=0d0
        p(n+1,2)=0d0
        p(n+1,3)=sqrt(max(0d0,pa1s))
        p(n+2,1)=sqrt(pts)
        p(n+2,2)=0d0
        p(n+2,3)=0.5d0*(pa3s-pa2s-pa1s)/p(n+1,3)
        p(n+3,1)=-p(n+2,1)
        p(n+3,2)=0d0
        p(n+3,3)=-(p(n+1,3)+p(n+2,3))
        v(n+1,2)=p(n+1,4)
        v(n+2,2)=p(n+2,4)
        v(n+3,2)=p(n+3,4)
 
C...Construct transverse momentum for ordinary branching in shower.
      else
        zm=v(im,1)
        pzm=sqrt(max(0d0,(pem+p(im,5))*(pem-p(im,5))))
        pmls=(v(im,5)-v(n+1,5)-v(n+2,5))**2-4d0*v(n+1,5)*v(n+2,5)
        if(pzm.le.0d0) then
          pts=0d0
        elseif(mod(mstj(43),2).eq.1) then
          pts=(pem**2*(zm*(1d0-zm)*v(im,5)-(1d0-zm)*v(n+1,5)-
     &    zm*v(n+2,5))-0.25d0*pmls)/pzm**2
        else
          pts=pmls*(zm*(1d0-zm)*pem**2/v(im,5)-0.25d0)/pzm**2
        endif
        pt=sqrt(max(0d0,pts))
 
C...Find coefficient of azimuthal asymmetry due to gluon polarization.
        hazip=0d0
        if(mstj(49).ne.1.and.mod(mstj(46),2).eq.1.and.k(im,2).eq.21
     &  .and.iau.ne.0) then
          if(k(igm,3).ne.0) mazip=1
          zau=v(igm,1)
          if(iau.eq.im+1) zau=1d0-v(igm,1)
          if(mazip.eq.0) zau=0d0
          if(k(igm,2).ne.21) then
            hazip=2d0*zau/(1d0+zau**2)
          else
            hazip=(zau/(1d0-zau*(1d0-zau)))**2
          endif
          if(k(n+1,2).ne.21) then
            hazip=hazip*(-2d0*zm*(1d0-zm))/(1d0-2d0*zm*(1d0-zm))
          else
            hazip=hazip*(zm*(1d0-zm)/(1d0-zm*(1d0-zm)))**2
          endif
        endif
 
C...Find coefficient of azimuthal asymmetry due to soft gluon
C...interference.
        hazic=0d0
        if(mstj(49).ne.2.and.mstj(46).ge.2.and.(k(n+1,2).eq.21.or.
     &  k(n+2,2).eq.21).and.iau.ne.0) then
          if(k(igm,3).ne.0) mazic=n+1
          if(k(igm,3).ne.0.and.k(n+1,2).ne.21) mazic=n+2
          if(k(igm,3).ne.0.and.k(n+1,2).eq.21.and.k(n+2,2).eq.21.and.
     &    zm.gt.0.5d0) mazic=n+2
          if(k(iau,2).eq.22) mazic=0
          zs=zm
          if(mazic.eq.n+2) zs=1d0-zm
          zgm=v(igm,1)
          if(iau.eq.im-1) zgm=1d0-v(igm,1)
          if(mazic.eq.0) zgm=1d0
          if(mazic.ne.0) hazic=(p(im,5)/p(igm,5))*
     &    sqrt((1d0-zs)*(1d0-zgm)/(zs*zgm))
          hazic=min(0.95d0,hazic)
        endif
      endif
 
C...Construct kinematics for ordinary branching in shower.
  470 if(nep.eq.2.and.igm.gt.0) then
        if(mod(mstj(43),2).eq.1) then
          p(n+1,4)=pem*v(im,1)
        else
          p(n+1,4)=pem*(0.5d0*(v(im,5)-sqrt(pmls)+v(n+1,5)-v(n+2,5))+
     &    sqrt(pmls)*zm)/v(im,5)
        endif
        phi=paru(2)*pjr(0)
        p(n+1,1)=pt*cos(phi)
        p(n+1,2)=pt*sin(phi)
        if(pzm.gt.0d0) then
          p(n+1,3)=0.5d0*(v(n+2,5)-v(n+1,5)-v(im,5)+
     &    2d0*pem*p(n+1,4))/pzm
        else
          p(n+1,3)=0d0
        endif
        p(n+2,1)=-p(n+1,1)
        p(n+2,2)=-p(n+1,2)
        p(n+2,3)=pzm-p(n+1,3)
        p(n+2,4)=pem-p(n+1,4)
        if(mstj(43).le.2) then
          v(n+1,2)=(pem*p(n+1,4)-pzm*p(n+1,3))/p(im,5)
          v(n+2,2)=(pem*p(n+2,4)-pzm*p(n+2,3))/p(im,5)
        endif
      endif
 
C...Rotate and boost daughters.
      if(igm.gt.0) then
        if(mstj(43).le.2) then
          bex=p(igm,1)/p(igm,4)
          bey=p(igm,2)/p(igm,4)
          bez=p(igm,3)/p(igm,4)
          ga=p(igm,4)/p(igm,5)
          gabep=ga*(ga*(bex*p(im,1)+bey*p(im,2)+bez*p(im,3))/(1d0+ga)-
     &    p(im,4))
        else
          bex=0d0
          bey=0d0
          bez=0d0
          ga=1d0
          gabep=0d0
        endif
        the=pjangl(p(im,3)+gabep*bez,sqrt((p(im,1)+gabep*bex)**2+
     &  (p(im,2)+gabep*bey)**2))
        phi=pjangl(p(im,1)+gabep*bex,p(im,2)+gabep*bey)
        do 480 i=n+1,n+2
          dp(1)=cos(the)*cos(phi)*p(i,1)-sin(phi)*p(i,2)+
     &    sin(the)*cos(phi)*p(i,3)
          dp(2)=cos(the)*sin(phi)*p(i,1)+cos(phi)*p(i,2)+
     &    sin(the)*sin(phi)*p(i,3)
          dp(3)=-sin(the)*p(i,1)+cos(the)*p(i,3)
          dp(4)=p(i,4)
          dbp=bex*dp(1)+bey*dp(2)+bez*dp(3)
          dgabp=ga*(ga*dbp/(1d0+ga)+dp(4))
          p(i,1)=dp(1)+dgabp*bex
          p(i,2)=dp(2)+dgabp*bey
          p(i,3)=dp(3)+dgabp*bez
          p(i,4)=ga*(dp(4)+dbp)
  480   continue
      endif
 
C...Weight with azimuthal distribution, if required.
      if(mazip.ne.0.or.mazic.ne.0) then
        do 490 j=1,3
          dpt(1,j)=p(im,j)
          dpt(2,j)=p(iau,j)
          dpt(3,j)=p(n+1,j)
  490   continue
        dpma=dpt(1,1)*dpt(2,1)+dpt(1,2)*dpt(2,2)+dpt(1,3)*dpt(2,3)
        dpmd=dpt(1,1)*dpt(3,1)+dpt(1,2)*dpt(3,2)+dpt(1,3)*dpt(3,3)
        dpmm=dpt(1,1)**2+dpt(1,2)**2+dpt(1,3)**2
        do 500 j=1,3
          dpt(4,j)=dpt(2,j)-dpma*dpt(1,j)/dpmm
          dpt(5,j)=dpt(3,j)-dpmd*dpt(1,j)/dpmm
  500   continue
        dpt(4,4)=sqrt(dpt(4,1)**2+dpt(4,2)**2+dpt(4,3)**2)
        dpt(5,4)=sqrt(dpt(5,1)**2+dpt(5,2)**2+dpt(5,3)**2)
        if(min(dpt(4,4),dpt(5,4)).gt.0.1d0*parj(82)) then
          cad=(dpt(4,1)*dpt(5,1)+dpt(4,2)*dpt(5,2)+
     &    dpt(4,3)*dpt(5,3))/(dpt(4,4)*dpt(5,4))
          if(mazip.ne.0) then
            if(1d0+hazip*(2d0*cad**2-1d0).lt.pjr(0)*(1d0+abs(hazip)))
     &      goto 470
          endif
          if(mazic.ne.0) then
            if(mazic.eq.n+2) cad=-cad
            if((1d0-hazic)*(1d0-hazic*cad)/(1d0+hazic**2-2d0*hazic*cad)
     &      .lt.pjr(0)) goto 470
          endif
        endif
      endif
 
C...Azimuthal anisotropy due to interference with initial state partons.
      if(mod(miis,2).eq.1.and.igm.eq.ns+1.and.(k(n+1,2).eq.21.or.
     &k(n+2,2).eq.21)) then
        iii=im-ns-1
        if(isii(iii).ge.1) then
          iaziid=n+1
          if(k(n+1,2).ne.21) iaziid=n+2
          if(k(n+1,2).eq.21.and.k(n+2,2).eq.21.and.
     &    p(n+1,4).gt.p(n+2,4)) iaziid=n+2
          theiid=pjangl(p(iaziid,3),sqrt(p(iaziid,1)**2+p(iaziid,2)**2))
          if(iii.eq.2) theiid=paru(1)-theiid
          phiiid=pjangl(p(iaziid,1),p(iaziid,2))
          hazii=min(0.95d0,theiid/theiis(iii,isii(iii)))
          cad=cos(phiiid-phiiis(iii,isii(iii)))
          phirel=abs(phiiid-phiiis(iii,isii(iii)))
          if(phirel.gt.paru(1)) phirel=paru(2)-phirel
          if((1d0-hazii)*(1d0-hazii*cad)/(1d0+hazii**2-2d0*hazii*cad)
     &    .lt.pjr(0)) goto 470
        endif
      endif
 
C...Continue loop over partons that may branch, until none left.
      if(igm.ge.0) k(im,1)=14
      n=n+nep
      nep=2
      if(n.gt.mstu(4)-mstu(32)-5) then
        call pjerrm(11,'(PYSHOW:) no more memory left in PYJETS')
        if(mstu(21).ge.1) n=ns
        if(mstu(21).ge.1) return
      endif
      goto 270
 
C...Set information on imagined shower initiator.
  510 if(npa.ge.2) then
        k(ns+1,1)=11
        k(ns+1,2)=94
        k(ns+1,3)=ip1
        if(ip2.gt.0.and.ip2.lt.ip1) k(ns+1,3)=ip2
        k(ns+1,4)=ns+2
        k(ns+1,5)=ns+1+npa
        iim=1
      else
        iim=0
      endif
 
C...Reconstruct string drawing information.
      do 520 i=ns+1+iim,n
        if(k(i,1).le.10.and.k(i,2).eq.22) then
          k(i,1)=1
        elseif(k(i,1).le.10.and.iabs(k(i,2)).ge.11.and.
     &    iabs(k(i,2)).le.18) then
          k(i,1)=1
        elseif(k(i,1).le.10) then
          k(i,4)=mstu(5)*(k(i,4)/mstu(5))
          k(i,5)=mstu(5)*(k(i,5)/mstu(5))
        elseif(k(mod(k(i,4),mstu(5))+1,2).ne.22) then
          id1=mod(k(i,4),mstu(5))
          if(k(i,2).ge.1.and.k(i,2).le.8) id1=mod(k(i,4),mstu(5))+1
          id2=2*mod(k(i,4),mstu(5))+1-id1
          k(i,4)=mstu(5)*(k(i,4)/mstu(5))+id1
          k(i,5)=mstu(5)*(k(i,5)/mstu(5))+id2
          k(id1,4)=k(id1,4)+mstu(5)*i
          k(id1,5)=k(id1,5)+mstu(5)*id2
          k(id2,4)=k(id2,4)+mstu(5)*id1
          k(id2,5)=k(id2,5)+mstu(5)*i
        else
          id1=mod(k(i,4),mstu(5))
          id2=id1+1
          k(i,4)=mstu(5)*(k(i,4)/mstu(5))+id1
          k(i,5)=mstu(5)*(k(i,5)/mstu(5))+id1
          if(iabs(k(i,2)).le.10.or.k(id1,1).ge.11) then
            k(id1,4)=k(id1,4)+mstu(5)*i
            k(id1,5)=k(id1,5)+mstu(5)*i
          else
            k(id1,4)=0
            k(id1,5)=0
          endif
          k(id2,4)=0
          k(id2,5)=0
        endif
  520 continue
 
C...Transformation from CM frame.
      if(npa.ge.2) then
        bex=ps(1)/ps(4)
        bey=ps(2)/ps(4)
        bez=ps(3)/ps(4)
        ga=ps(4)/ps(5)
        gabep=ga*(ga*(bex*p(ipa(1),1)+bey*p(ipa(1),2)+bez*p(ipa(1),3))
     &  /(1d0+ga)-p(ipa(1),4))
      else
        bex=0d0
        bey=0d0
        bez=0d0
        gabep=0d0
      endif
      the=pjangl(p(ipa(1),3)+gabep*bez,sqrt((p(ipa(1),1)
     &+gabep*bex)**2+(p(ipa(1),2)+gabep*bey)**2))
      phi=pjangl(p(ipa(1),1)+gabep*bex,p(ipa(1),2)+gabep*bey)
      if(npa.eq.3) then
        chi=pjangl(cos(the)*cos(phi)*(p(ipa(2),1)+gabep*bex)+cos(the)*
     &  sin(phi)*(p(ipa(2),2)+gabep*bey)-sin(the)*(p(ipa(2),3)+gabep*
     &  bez),-sin(phi)*(p(ipa(2),1)+gabep*bex)+cos(phi)*(p(ipa(2),2)+
     &  gabep*bey))
        mstu(33)=1
        call pjrobo(ns+1,n,0d0,chi,0d0,0d0,0d0)
      endif
      mstu(33)=1
      call pjrobo(ns+1,n,the,phi,bex,bey,bez)
 
C...Decay vertex of shower.
      do 540 i=ns+1,n
        do 530 j=1,5
          v(i,j)=v(ip1,j)
  530   continue
  540 continue
 
C...Delete trivial shower, else connect initiators.
      if(n.eq.ns+npa+iim) then
        n=ns
      else
        do 550 ip=1,npa
          k(ipa(ip),1)=14
          k(ipa(ip),4)=k(ipa(ip),4)+ns+iim+ip
          k(ipa(ip),5)=k(ipa(ip),5)+ns+iim+ip
          k(ns+iim+ip,3)=ipa(ip)
          if(iim.eq.1.and.mstu(16).ne.2) k(ns+iim+ip,3)=ns+1
          if(k(ns+iim+ip,1).ne.1) then
            k(ns+iim+ip,4)=mstu(5)*ipa(ip)+k(ns+iim+ip,4)
            k(ns+iim+ip,5)=mstu(5)*ipa(ip)+k(ns+iim+ip,5)
          endif
  550   continue
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYBOEI
C...Modifies an event so as to approximately take into account
C...Bose-Einstein effects according to a simple phenomenological
C...parametrization.
 
      subroutine pjboei(nsav)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      save /jyjets/,/jydat1/
C...Local arrays and data.
      dimension dps(4),kfbe(9),nbe(0:9),bei(100)
      data kfbe/211,-211,111,321,-321,130,310,221,331/
 
C...Boost event to overall CM frame. Calculate CM energy.
      if((mstj(51).ne.1.and.mstj(51).ne.2).or.n-nsav.le.1) return
      do 100 j=1,4
        dps(j)=0d0
  100 continue
      do 120 i=1,n
        kfa=iabs(k(i,2))
        if(k(i,1).le.10.and.((kfa.gt.10.and.kfa.le.20).or.kfa.eq.22)
     &  .and.k(i,3).gt.0) then
          kfma=iabs(k(k(i,3),2))
          if(kfma.gt.10.and.kfma.le.80) k(i,1)=-k(i,1)
        endif
        if(k(i,1).le.0.or.k(i,1).gt.10) goto 120
        do 110 j=1,4
          dps(j)=dps(j)+p(i,j)
  110   continue
  120 continue
      call pjrobo(0,0,0d0,0d0,-dps(1)/dps(4),-dps(2)/dps(4),
     &-dps(3)/dps(4))
      pecm=0d0
      do 130 i=1,n
        if(k(i,1).ge.1.and.k(i,1).le.10) pecm=pecm+p(i,4)
  130 continue
 
C...Reserve copy of particles by species at end of record.
      nbe(0)=n+mstu(3)
      do 160 ibe=1,min(9,mstj(52))
        nbe(ibe)=nbe(ibe-1)
        do 150 i=nsav+1,n
          if(k(i,2).ne.kfbe(ibe)) goto 150
          if(k(i,1).le.0.or.k(i,1).gt.10) goto 150
          if(nbe(ibe).ge.mstu(4)-mstu(32)-5) then
            call pjerrm(11,'(PYBOEI:) no more memory left in PYJETS')
            return
          endif
          nbe(ibe)=nbe(ibe)+1
          k(nbe(ibe),1)=i
          do 140 j=1,3
            p(nbe(ibe),j)=0d0
  140     continue
  150   continue
  160 continue
      if(nbe(min(9,mstj(52)))-nbe(0).le.1) goto 280
 
C...Tabulate integral for subsequent momentum shift.
      do 220 ibe=1,min(9,mstj(52))
        if(ibe.ne.1.and.ibe.ne.4.and.ibe.le.7) goto 180
        if(ibe.eq.1.and.max(nbe(1)-nbe(0),nbe(2)-nbe(1),nbe(3)-nbe(2))
     &  .le.1) goto 180
        if(ibe.eq.4.and.max(nbe(4)-nbe(3),nbe(5)-nbe(4),nbe(6)-nbe(5),
     &  nbe(7)-nbe(6)).le.1) goto 180
        if(ibe.ge.8.and.nbe(ibe)-nbe(ibe-1).le.1) goto 180
        if(ibe.eq.1) pmhq=2d0*pjmass(211)
        if(ibe.eq.4) pmhq=2d0*pjmass(321)
        if(ibe.eq.8) pmhq=2d0*pjmass(221)
        if(ibe.eq.9) pmhq=2d0*pjmass(331)
        qdel=0.1d0*min(pmhq,parj(93))
        if(mstj(51).eq.1) then
          nbin=min(100,nint(9d0*parj(93)/qdel))
          beex=exp(0.5d0*qdel/parj(93))
          bert=exp(-qdel/parj(93))
        else
          nbin=min(100,nint(3d0*parj(93)/qdel))
        endif
        do 170 ibin=1,nbin
          qbin=qdel*(ibin-0.5d0)
          bei(ibin)=qdel*(qbin**2+qdel**2/12d0)/sqrt(qbin**2+pmhq**2)
          if(mstj(51).eq.1) then
            beex=beex*bert
            bei(ibin)=bei(ibin)*beex
          else
            bei(ibin)=bei(ibin)*exp(-(qbin/parj(93))**2)
          endif
          if(ibin.ge.2) bei(ibin)=bei(ibin)+bei(ibin-1)
  170   continue
 
C...Loop through particle pairs and find old relative momentum.
  180   do 210 i1m=nbe(ibe-1)+1,nbe(ibe)-1
          i1=k(i1m,1)
          do 200 i2m=i1m+1,nbe(ibe)
            i2=k(i2m,1)
            q2old=max(0d0,(p(i1,4)+p(i2,4))**2-(p(i1,1)+p(i2,1))**2-
     &      (p(i1,2)+ p(i2,2))**2-(p(i1,3)+p(i2,3))**2-
     &      (p(i1,5)+p(i2,5))**2)
            qold=sqrt(q2old)
 
C...Calculate new relative momentum.
            if(qold.lt.1d-3*qdel) then
              goto 200
            elseif(qold.le.qdel) then
              qmov=qold/3d0
            elseif(qold.lt.(nbin-0.1d0)*qdel) then
              rbin=qold/qdel
              ibin=rbin
              rinp=(rbin**3-ibin**3)/(3*ibin*(ibin+1)+1)
              qmov=(bei(ibin)+rinp*(bei(ibin+1)-bei(ibin)))*
     &        sqrt(q2old+pmhq**2)/q2old
            else
              qmov=bei(nbin)*sqrt(q2old+pmhq**2)/q2old
            endif
            q2new=q2old*(qold/(qold+3d0*parj(92)*qmov))**(2d0/3d0)
 
C...Calculate and save shift to be performed on three-momenta.
            hc1=(p(i1,4)+p(i2,4))**2-(q2old-q2new)
            hc2=(q2old-q2new)*(p(i1,4)-p(i2,4))**2
            ha=0.5d0*(1d0-sqrt(hc1*q2new/(hc1*q2old-hc2)))
            do 190 j=1,3
              pd=ha*(p(i2,j)-p(i1,j))
              p(i1m,j)=p(i1m,j)+pd
              p(i2m,j)=p(i2m,j)-pd
  190       continue
  200     continue
  210   continue
  220 continue
 
C...Shift momenta and recalculate energies.
      do 240 im=nbe(0)+1,nbe(min(9,mstj(52)))
        i=k(im,1)
        do 230 j=1,3
          p(i,j)=p(i,j)+p(im,j)
  230   continue
        p(i,4)=sqrt(p(i,5)**2+p(i,1)**2+p(i,2)**2+p(i,3)**2)
  240 continue
 
C...Rescale all momenta for energy conservation.
      pes=0d0
      pqs=0d0
      do 250 i=1,n
        if(k(i,1).le.0.or.k(i,1).gt.10) goto 250
        pes=pes+p(i,4)
        pqs=pqs+p(i,5)**2/p(i,4)
  250 continue
      fac=(pecm-pqs)/(pes-pqs)
      do 270 i=1,n
        if(k(i,1).le.0.or.k(i,1).gt.10) goto 270
        do 260 j=1,3
          p(i,j)=fac*p(i,j)
  260   continue
        p(i,4)=sqrt(p(i,5)**2+p(i,1)**2+p(i,2)**2+p(i,3)**2)
  270 continue
 
C...Boost back to correct reference frame.
  280 call pjrobo(0,0,0d0,0d0,dps(1)/dps(4),dps(2)/dps(4),dps(3)/dps(4))
      do 290 i=1,n
        if(k(i,1).lt.0) k(i,1)=-k(i,1)
  290 continue
 
      return
      end

C***********************************************************************
 
C...PYMASS
C...Gives the mass of a particle/parton.
 
      function pjmass(kf)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jydat1/,/jydat2/
 
C...Reset variables. Compressed code. Special case for popcorn diquarks.
      pjmass=0d0
      kfa=iabs(kf)
      kc=jamcomp(kf)
      if(kc.eq.0) then
        mstj(93)=0
        return
      endif
 
C...Guarantee use of constituent masses for internal checks.
      if((mstj(93).eq.1.or.mstj(93).eq.2).and.
     &(kfa.le.10.or.mod(kfa/10,10).eq.0)) then
        parf(106)=pmas(6,1)
        parf(107)=pmas(7,1)
        parf(108)=pmas(8,1)
        if(kfa.le.10) then
          pjmass=parf(100+kfa)
          if(mstj(93).eq.2) pjmass=max(0d0,pjmass-parf(121))
        elseif(mstj(93).eq.1) then
          pjmass=parf(100+mod(kfa/1000,10))+parf(100+mod(kfa/100,10))
        else
          pjmass=max(0d0,pmas(kc,1)-parf(122)-2d0*parf(112)/3d0)
        endif
 
C...Other masses can be read directly off table.
      else
        pjmass=pmas(kc,1)
      endif
 
C...Optional mass broadening according to truncated Breit-Wigner
C...(either in m or in m^2).
      if(mstj(24).ge.1.and.pmas(kc,2).gt.1d-4) then
        if(mstj(24).eq.1.or.(mstj(24).eq.2.and.kfa.gt.100)) then
          pjmass=pjmass+0.5d0*pmas(kc,2)*tan((2d0*pjr(0)-1d0)*
     &    atan(2d0*pmas(kc,3)/pmas(kc,2)))
        else
          pm0=pjmass
          pmlow=atan((max(0d0,pm0-pmas(kc,3))**2-pm0**2)/
     &    (pm0*pmas(kc,2)))
          pmupp=atan(((pm0+pmas(kc,3))**2-pm0**2)/(pm0*pmas(kc,2)))
          pjmass=sqrt(max(0d0,pm0**2+pm0*pmas(kc,2)*tan(pmlow+
     &    (pmupp-pmlow)*pjr(0))))
        endif
      endif
      mstj(93)=0
 
      return
      end
 
C*********************************************************************
 
C...PYNAME
C...Gives the particle/parton name as a character string.
 
      subroutine pjname(kf,chau)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/jydat4/chaf(500,2)
      character chaf*16
      save /jydat1/,/jydat2/,/jydat4/
C...Local character variable.
      character chau*16
 
C...Read out code with distinction particle/antiparticle.
      chau=' '
      kc=jamcomp(kf)
      if(kc.ne.0) chau=chaf(kc,(3-isign(1,kf))/2)
 
      return
      end
 
C*********************************************************************
 
C...PYERRM
C...Informs user of errors in program execution.
 
      subroutine pjerrm(merr,chmess)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      save /jyjets/,/jydat1/
C...Local character variable.
      character chmess*(*)
 
C...Write first few warnings, then be silent.
      if(merr.le.10) then
        mstu(27)=mstu(27)+1
        mstu(28)=merr
        if(mstu(25).eq.1.and.mstu(27).le.mstu(26)) write(mstu(11),5000)
     &  merr,mstu(31),chmess
 
C...Write first few errors, then be silent or stop program.
      elseif(merr.le.20) then
        mstu(23)=mstu(23)+1
        mstu(24)=merr-10
        if(mstu(21).ge.1.and.mstu(23).le.mstu(22)) write(mstu(11),5100)
     &  merr-10,mstu(31),chmess
        if(mstu(21).ge.2.and.mstu(23).gt.mstu(22)) then
          write(mstu(11),5100) merr-10,mstu(31),chmess
          write(mstu(11),5200)
          if(merr.ne.17) call pjlist(2)
          stop
        endif
 
C...Stop program in case of irreparable error.
      else
        write(mstu(11),5300) merr-20,mstu(31),chmess
        stop
      endif
 
C...Formats for output.
 5000 format(/5x,'Advisory warning type',i2,' given after',i9,
     &' PYEXEC calls:'/5x,a)
 5100 format(/5x,'Error type',i2,' has occured after',i9,
     &' PYEXEC calls:'/5x,a)
 5200 format(5x,'Execution will be stopped after listing of last ',
     &'event!')
 5300 format(/5x,'Fatal error type',i2,' has occured after',i9,
     &' PYEXEC calls:'/5x,a/5x,'Execution will now be stopped!')
 
      return
      end
 
C*********************************************************************
 
C...PYALEM
C...Calculates the running alpha_electromagnetic.
 
      function pjalem(q2)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      save /jydat1/
 
C...Calculate real part of photon vacuum polarization.
C...For leptons simplify by using asymptotic (Q^2 >> m^2) expressions.
C...For hadrons use parametrization of H. Burkhardt et al.
C...See R. Kleiss et al, CERN 89-08, vol. 3, pp. 129-131.
      aempi=paru(101)/(3d0*paru(1))
      if(mstu(101).le.0.or.q2.lt.2d-6) then
        rpigg=0d0
      elseif(mstu(101).eq.2.and.q2.lt.paru(104)) then
        rpigg=0d0
      elseif(mstu(101).eq.2) then
        rpigg=1d0-paru(101)/paru(103)
      elseif(q2.lt.0.09d0) then
        rpigg=aempi*(13.4916d0+log(q2))+0.00835d0*log(1d0+q2)
      elseif(q2.lt.9d0) then
        rpigg=aempi*(16.3200d0+2d0*log(q2))+
     &  0.00238d0*log(1d0+3.927d0*q2)
      elseif(q2.lt.1d4) then
        rpigg=aempi*(13.4955d0+3d0*log(q2))+0.00165d0+
     &  0.00299d0*log(1d0+q2)
      else
        rpigg=aempi*(13.4955d0+3d0*log(q2))+0.00221d0+
     &  0.00293d0*log(1d0+q2)
      endif
 
C...Calculate running alpha_em.
      pjalem=paru(101)/(1d0-rpigg)
      paru(108)=pjalem
 
      return
      end
 
C*********************************************************************
 
C...PYALPS
C...Gives the value of alpha_strong.
 
      function pjalps(q2)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jydat1/,/jydat2/
 
C...Constant alpha_strong trivial. Pick artificial Lambda.
      if(mstu(111).le.0) then
        pjalps=paru(111)
        mstu(118)=mstu(112)
        paru(117)=0.2d0
        if(q2.gt.0.04d0) paru(117)=sqrt(q2)*exp(-6d0*paru(1)/
     &  ((33d0-2d0*mstu(112))*paru(111)))
        paru(118)=paru(111)
        return
      endif
 
C...Find effective Q2, number of flavours and Lambda.
      q2eff=q2
      if(mstu(115).ge.2) q2eff=max(q2,paru(114))
      nf=mstu(112)
      alam2=paru(112)**2
  100 if(nf.gt.max(2,mstu(113))) then
        q2thr=paru(113)*pmas(nf,1)**2
        if(q2eff.lt.q2thr) then
          nf=nf-1
          alam2=alam2*(q2thr/alam2)**(2d0/(33d0-2d0*nf))
          goto 100
        endif
      endif
  110 if(nf.lt.min(8,mstu(114))) then
        q2thr=paru(113)*pmas(nf+1,1)**2
        if(q2eff.gt.q2thr) then
          nf=nf+1
          alam2=alam2*(alam2/q2thr)**(2d0/(33d0-2d0*nf))
          goto 110
        endif
      endif
      if(mstu(115).eq.1) q2eff=q2eff+alam2
      paru(117)=sqrt(alam2)
 
C...Evaluate first or second order alpha_strong.
      b0=(33d0-2d0*nf)/6d0
      algq=log(max(1.0001d0,q2eff/alam2))
      if(mstu(111).eq.1) then
        pjalps=min(paru(115),paru(2)/(b0*algq))
      else
        b1=(153d0-19d0*nf)/6d0
        pjalps=min(paru(115),paru(2)/(b0*algq)*(1d0-b1*log(algq)/
     &  (b0**2*algq)))
      endif
      mstu(118)=nf
      paru(118)=pjalps
 
      return
      end
 
C*********************************************************************
 
C...PYANGL
C...Reconstructs an angle from given x and y coordinates.
 
      function pjangl(x,y)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      save /jydat1/
 
      pjangl=0d0
      r=sqrt(x**2+y**2)
      if(r.lt.1d-20) return
      if(abs(x)/r.lt.0.8d0) then
        pjangl=sign(acos(x/r),y)
      else
        pjangl=asin(y/r)
        if(x.lt.0d0.and.pjangl.ge.0d0) then
          pjangl=paru(1)-pjangl
        elseif(x.lt.0d0) then
          pjangl=-paru(1)-pjangl
        endif
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYRND
C...Generates random numbers uniformly distributed between
C...0 and 1, excluding the endpoints.
 
      function pjrnd(idummy)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydatr/mrpy(6),rrpy(100)
      save /jydatr/
C...Equivalence between commonblock and local variables.
      equivalence (mrpy1,mrpy(1)),(mrpy2,mrpy(2)),(mrpy3,mrpy(3)),
     &(mrpy4,mrpy(4)),(mrpy5,mrpy(5)),(mrpy6,mrpy(6)),
     &(rrpy98,rrpy(98)),(rrpy99,rrpy(99)),(rrpy00,rrpy(100))
 
C...Initialize generation from given seed.
      if(mrpy2.eq.0) then
        ij=mod(mrpy1/30082,31329)
        kl=mod(mrpy1,30082)
        i=mod(ij/177,177)+2
        j=mod(ij,177)+2
        k=mod(kl/169,178)+1
        l=mod(kl,169)
        do 110 ii=1,97
          s=0d0
          t=0.5d0
          do 100 jj=1,48
            m=mod(mod(i*j,179)*k,179)
            i=j
            j=k
            k=m
            l=mod(53*l+1,169)
            if(mod(l*m,64).ge.32) s=s+t
            t=0.5d0*t
  100     continue
          rrpy(ii)=s
  110   continue
        twom24=1d0
        do 120 i24=1,24
          twom24=0.5d0*twom24
  120   continue
        rrpy98=362436d0*twom24
        rrpy99=7654321d0*twom24
        rrpy00=16777213d0*twom24
        mrpy2=1
        mrpy3=0
        mrpy4=97
        mrpy5=33
      endif
 
C...Generate next random number.
  130 runi=rrpy(mrpy4)-rrpy(mrpy5)
      if(runi.lt.0d0) runi=runi+1d0
      rrpy(mrpy4)=runi
      mrpy4=mrpy4-1
      if(mrpy4.eq.0) mrpy4=97
      mrpy5=mrpy5-1
      if(mrpy5.eq.0) mrpy5=97
      rrpy98=rrpy98-rrpy99
      if(rrpy98.lt.0d0) rrpy98=rrpy98+rrpy00
      runi=runi-rrpy98
      if(runi.lt.0d0) runi=runi+1d0
      if(runi.le.0d0.or.runi.ge.1d0) goto 130
 
C...Update counters. Random number to output.
      mrpy3=mrpy3+1
      if(mrpy3.eq.1000000000) then
        mrpy2=mrpy2+1
        mrpy3=0
      endif
      pjrnd=runi
 
      return
      end
 
C*********************************************************************
 
C...PYRGET
C...Dumps the state of the random number generator on a file
C...for subsequent startup from this state onwards.
 
      subroutine pjrget(lfn,move)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydatr/mrpy(6),rrpy(100)
      save /jydatr/
C...Local character variable.
      character cherr*8
 
C...Backspace required number of records (or as many as there are).
      if(move.lt.0) then
        nbck=min(mrpy(6),-move)
        do 100 ibck=1,nbck
          backspace(lfn,err=110,iostat=ierr)
  100   continue
        mrpy(6)=mrpy(6)-nbck
      endif
 
C...Unformatted write on unit LFN.
      write(lfn,err=110,iostat=ierr) (mrpy(i1),i1=1,5),
     &(rrpy(i2),i2=1,100)
      mrpy(6)=mrpy(6)+1
      return
 
C...Write error.
  110 write(cherr,'(I8)') ierr
      call pjerrm(18,'(PYRGET:) error when accessing file, IOSTAT ='//
     &cherr)
 
      return
      end
 
C*********************************************************************
 
C...PYRSET
C...Reads a state of the random number generator from a file
C...for subsequent generation from this state onwards.
 
      subroutine pjrset(lfn,move)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydatr/mrpy(6),rrpy(100)
      save /jydatr/
C...Local character variable.
      character cherr*8
 
C...Backspace required number of records (or as many as there are).
      if(move.lt.0) then
        nbck=min(mrpy(6),-move)
        do 100 ibck=1,nbck
          backspace(lfn,err=120,iostat=ierr)
  100   continue
        mrpy(6)=mrpy(6)-nbck
      endif
 
C...Unformatted read from unit LFN.
      nfor=1+max(0,move)
      do 110 ifor=1,nfor
        read(lfn,err=120,iostat=ierr) (mrpy(i1),i1=1,5),
     &  (rrpy(i2),i2=1,100)
  110 continue
      mrpy(6)=mrpy(6)+nfor
      return
 
C...Write error.
  120 write(cherr,'(I8)') ierr
      call pjerrm(18,'(PYRSET:) error when accessing file, IOSTAT ='//
     &cherr)
 
      return
      end
 
C*********************************************************************
 
C...PYROBO
C...Performs rotations and boosts.
 
      subroutine pjrobo(imi,ima,the,phi,bex,bey,bez)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      save /jyjets/,/jydat1/
C...Local arrays.
      dimension rot(3,3),pr(3),vr(3),dp(4),dv(4)
 
C...Find and check range of rotation/boost.
      imin=imi
      if(imin.le.0) imin=1
      if(mstu(1).gt.0) imin=mstu(1)
      imax=ima
      if(imax.le.0) imax=n
      if(mstu(2).gt.0) imax=mstu(2)
      if(imin.gt.mstu(4).or.imax.gt.mstu(4)) then
        call pjerrm(11,'(PYROBO:) range outside PYJETS memory')
        return
      endif
 
C...Optional resetting of V (when not set before.)
      if(mstu(33).ne.0) then
        do 110 i=min(imin,mstu(4)),min(imax,mstu(4))
          do 100 j=1,5
            v(i,j)=0d0
  100     continue
  110   continue
        mstu(33)=0
      endif
 
C...Rotate, typically from z axis to direction (theta,phi).
      if(the**2+phi**2.gt.1d-20) then
        rot(1,1)=cos(the)*cos(phi)
        rot(1,2)=-sin(phi)
        rot(1,3)=sin(the)*cos(phi)
        rot(2,1)=cos(the)*sin(phi)
        rot(2,2)=cos(phi)
        rot(2,3)=sin(the)*sin(phi)
        rot(3,1)=-sin(the)
        rot(3,2)=0d0
        rot(3,3)=cos(the)
        do 140 i=imin,imax
          if(k(i,1).le.0) goto 140
          do 120 j=1,3
            pr(j)=p(i,j)
            vr(j)=v(i,j)
  120     continue
          do 130 j=1,3
            p(i,j)=rot(j,1)*pr(1)+rot(j,2)*pr(2)+rot(j,3)*pr(3)
            v(i,j)=rot(j,1)*vr(1)+rot(j,2)*vr(2)+rot(j,3)*vr(3)
  130     continue
  140   continue
      endif
 
C...Boost, typically from rest to momentum/energy=beta.
      if(bex**2+bey**2+bez**2.gt.1d-20) then
        dbx=bex
        dby=bey
        dbz=bez
        db=sqrt(dbx**2+dby**2+dbz**2)
        eps1=1d0-1d-12
        if(db.gt.eps1) then
C...Rescale boost vector if too close to unity.
          call pjerrm(3,'(PYROBO:) boost vector too large')
          dbx=dbx*(eps1/db)
          dby=dby*(eps1/db)
          dbz=dbz*(eps1/db)
          db=eps1
        endif
        dga=1d0/sqrt(1d0-db**2)
        do 160 i=imin,imax
          if(k(i,1).le.0) goto 160
          do 150 j=1,4
            dp(j)=p(i,j)
            dv(j)=v(i,j)
  150     continue
          dbp=dbx*dp(1)+dby*dp(2)+dbz*dp(3)
          dgabp=dga*(dga*dbp/(1d0+dga)+dp(4))
          p(i,1)=dp(1)+dgabp*dbx
          p(i,2)=dp(2)+dgabp*dby
          p(i,3)=dp(3)+dgabp*dbz
          p(i,4)=dga*(dp(4)+dbp)
          dbv=dbx*dv(1)+dby*dv(2)+dbz*dv(3)
          dgabv=dga*(dga*dbv/(1d0+dga)+dv(4))
          v(i,1)=dv(1)+dgabv*dbx
          v(i,2)=dv(2)+dgabv*dby
          v(i,3)=dv(3)+dgabv*dbz
          v(i,4)=dga*(dv(4)+dbv)
  160   continue
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYEDIT
C...Performs global manipulations on the event record, in particular
C...to exclude unstable or undetectable partons/particles.
 
      subroutine pjedit(medit)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jyjets/,/jydat1/,/jydat2/
C...Local arrays.
      dimension ns(2),pts(2),pls(2)
 
C...Remove unwanted partons/particles.
      if((medit.ge.0.and.medit.le.3).or.medit.eq.5) then
        imax=n
        if(mstu(2).gt.0) imax=mstu(2)
        i1=max(1,mstu(1))-1
        do 110 i=max(1,mstu(1)),imax
          if(k(i,1).eq.0.or.k(i,1).gt.20) goto 110
          if(medit.eq.1) then
            if(k(i,1).gt.10) goto 110
          elseif(medit.eq.2) then
            if(k(i,1).gt.10) goto 110
            kc=jamcomp(k(i,2))
            if(kc.eq.0.or.kc.eq.12.or.kc.eq.14.or.kc.eq.16.or.kc.eq.18)
     &      goto 110
          elseif(medit.eq.3) then
            if(k(i,1).gt.10) goto 110
            kc=jamcomp(k(i,2))
            if(kc.eq.0) goto 110
            if(kchg(kc,2).eq.0.and.jamchge(k(i,2)).eq.0) goto 110
          elseif(medit.eq.5) then
            if(k(i,1).eq.13.or.k(i,1).eq.14) goto 110
            kc=jamcomp(k(i,2))
            if(kc.eq.0) goto 110
            if(k(i,1).ge.11.and.kchg(kc,2).eq.0) goto 110
          endif
 
C...Pack remaining partons/particles. Origin no longer known.
          i1=i1+1
          do 100 j=1,5
            k(i1,j)=k(i,j)
            p(i1,j)=p(i,j)
            v(i1,j)=v(i,j)
  100     continue
          k(i1,3)=0
  110   continue
        if(i1.lt.n) mstu(3)=0
        if(i1.lt.n) mstu(70)=0
        n=i1
 
C...Selective removal of class of entries. New position of retained.
      elseif(medit.ge.11.and.medit.le.15) then
        i1=0
        do 120 i=1,n
          k(i,3)=mod(k(i,3),mstu(5))
          if(medit.eq.11.and.k(i,1).lt.0) goto 120
          if(medit.eq.12.and.k(i,1).eq.0) goto 120
          if(medit.eq.13.and.(k(i,1).eq.11.or.k(i,1).eq.12.or.
     &    k(i,1).eq.15).and.k(i,2).ne.94) goto 120
          if(medit.eq.14.and.(k(i,1).eq.13.or.k(i,1).eq.14.or.
     &    k(i,2).eq.94)) goto 120
          if(medit.eq.15.and.k(i,1).ge.21) goto 120
          i1=i1+1
          k(i,3)=k(i,3)+mstu(5)*i1
  120   continue
 
C...Find new event history information and replace old.
        do 140 i=1,n
          if(k(i,1).le.0.or.k(i,1).gt.20.or.k(i,3)/mstu(5).eq.0)
     &    goto 140
          id=i
  130     im=mod(k(id,3),mstu(5))
          if(medit.eq.13.and.im.gt.0.and.im.le.n) then
            if((k(im,1).eq.11.or.k(im,1).eq.12.or.k(im,1).eq.15).and.
     &      k(im,2).ne.94) then
              id=im
              goto 130
            endif
          elseif(medit.eq.14.and.im.gt.0.and.im.le.n) then
            if(k(im,1).eq.13.or.k(im,1).eq.14.or.k(im,2).eq.94) then
              id=im
              goto 130
            endif
          endif
          k(i,3)=mstu(5)*(k(i,3)/mstu(5))
          if(im.ne.0) k(i,3)=k(i,3)+k(im,3)/mstu(5)
          if(k(i,1).ne.3.and.k(i,1).ne.13.and.k(i,1).ne.14) then
            if(k(i,4).gt.0.and.k(i,4).le.mstu(4)) k(i,4)=
     &      k(k(i,4),3)/mstu(5)
            if(k(i,5).gt.0.and.k(i,5).le.mstu(4)) k(i,5)=
     &      k(k(i,5),3)/mstu(5)
          else
            kcm=mod(k(i,4)/mstu(5),mstu(5))
            if(kcm.gt.0.and.kcm.le.mstu(4)) kcm=k(kcm,3)/mstu(5)
            kcd=mod(k(i,4),mstu(5))
            if(kcd.gt.0.and.kcd.le.mstu(4)) kcd=k(kcd,3)/mstu(5)
            k(i,4)=mstu(5)**2*(k(i,4)/mstu(5)**2)+mstu(5)*kcm+kcd
            kcm=mod(k(i,5)/mstu(5),mstu(5))
            if(kcm.gt.0.and.kcm.le.mstu(4)) kcm=k(kcm,3)/mstu(5)
            kcd=mod(k(i,5),mstu(5))
            if(kcd.gt.0.and.kcd.le.mstu(4)) kcd=k(kcd,3)/mstu(5)
            k(i,5)=mstu(5)**2*(k(i,5)/mstu(5)**2)+mstu(5)*kcm+kcd
          endif
  140   continue
 
C...Pack remaining entries.
        i1=0
        mstu90=mstu(90)
        mstu(90)=0
        do 170 i=1,n
          if(k(i,3)/mstu(5).eq.0) goto 170
          i1=i1+1
          do 150 j=1,5
            k(i1,j)=k(i,j)
            p(i1,j)=p(i,j)
            v(i1,j)=v(i,j)
  150     continue
          k(i1,3)=mod(k(i1,3),mstu(5))
          do 160 iz=1,mstu90
            if(i.eq.mstu(90+iz)) then
              mstu(90)=mstu(90)+1
              mstu(90+mstu(90))=i1
              paru(90+mstu(90))=paru(90+iz)
            endif
  160     continue
  170   continue
        if(i1.lt.n) mstu(3)=0
        if(i1.lt.n) mstu(70)=0
        n=i1
 
C...Fill in some missing daughter pointers (lost in colour flow).
      elseif(medit.eq.16) then
        do 220 i=1,n
          if(k(i,1).le.10.or.k(i,1).gt.20) goto 220
          if(k(i,4).ne.0.or.k(i,5).ne.0) goto 220
C...Find daughters who point to mother.
          do 180 i1=i+1,n
            if(k(i1,3).ne.i) then
            elseif(k(i,4).eq.0) then
              k(i,4)=i1
            else
              k(i,5)=i1
            endif
  180     continue
          if(k(i,5).eq.0) k(i,5)=k(i,4)
          if(k(i,4).ne.0) goto 220
C...Find daughters who point to documentation version of mother.
          im=k(i,3)
          if(im.le.0.or.im.ge.i) goto 220
          if(k(im,1).le.20.or.k(im,1).gt.30) goto 220
          if(k(im,2).ne.k(i,2).or.abs(p(im,5)-p(i,5)).gt.1d-2) goto 220
          do 190 i1=i+1,n
            if(k(i1,3).ne.im) then
            elseif(k(i,4).eq.0) then
              k(i,4)=i1
            else
              k(i,5)=i1
            endif
  190     continue
          if(k(i,5).eq.0) k(i,5)=k(i,4)
          if(k(i,4).ne.0) goto 220
C...Find daughters who point to documentation daughters who,
C...in their turn, point to documentation mother.
          id1=im
          id2=im
          do 200 i1=im+1,i-1
            if(k(i1,3).eq.im.and.k(i1,1).gt.20.and.k(i1,1).le.30) then
              id2=i1
              if(id1.eq.im) id1=i1
            endif
  200     continue
          do 210 i1=i+1,n
            if(k(i1,3).ne.id1.and.k(i1,3).ne.id2) then
            elseif(k(i,4).eq.0) then
              k(i,4)=i1
            else
              k(i,5)=i1
            endif
  210     continue
          if(k(i,5).eq.0) k(i,5)=k(i,4)
  220   continue
 
C...Save top entries at bottom of PYJETS commonblock.
      elseif(medit.eq.21) then
        if(2*n.ge.mstu(4)) then
          call pjerrm(11,'(PYEDIT:) no more memory left in PYJETS')
          return
        endif
        do 240 i=1,n
          do 230 j=1,5
            k(mstu(4)-i,j)=k(i,j)
            p(mstu(4)-i,j)=p(i,j)
            v(mstu(4)-i,j)=v(i,j)
  230     continue
  240   continue
        mstu(32)=n
 
C...Restore bottom entries of commonblock PYJETS to top.
      elseif(medit.eq.22) then
        do 260 i=1,mstu(32)
          do 250 j=1,5
            k(i,j)=k(mstu(4)-i,j)
            p(i,j)=p(mstu(4)-i,j)
            v(i,j)=v(mstu(4)-i,j)
  250     continue
  260   continue
        n=mstu(32)
 
C...Mark primary entries at top of commonblock PYJETS as untreated.
      elseif(medit.eq.23) then
        i1=0
        do 270 i=1,n
          kh=k(i,3)
          if(kh.ge.1) then
            if(k(kh,1).gt.20) kh=0
          endif
          if(kh.ne.0) goto 280
          i1=i1+1
          if(k(i,1).gt.10.and.k(i,1).le.20) k(i,1)=k(i,1)-10
  270   continue
  280   n=i1
 
C...Place largest axis along z axis and second largest in xy plane.
      elseif(medit.eq.31.or.medit.eq.32) then
        call pjrobo(1,n+mstu(3),0d0,-pjangl(p(mstu(61),1),
     &  p(mstu(61),2)),0d0,0d0,0d0)
        call pjrobo(1,n+mstu(3),-pjangl(p(mstu(61),3),
     &  p(mstu(61),1)),0d0,0d0,0d0,0d0)
        call pjrobo(1,n+mstu(3),0d0,-pjangl(p(mstu(61)+1,1),
     &  p(mstu(61)+1,2)),0d0,0d0,0d0)
        if(medit.eq.31) return
 
C...Rotate to put slim jet along +z axis.
        do 290 is=1,2
          ns(is)=0
          pts(is)=0d0
          pls(is)=0d0
  290   continue
        do 300 i=1,n
          if(k(i,1).le.0.or.k(i,1).gt.10) goto 300
          if(mstu(41).ge.2) then
            kc=jamcomp(k(i,2))
            if(kc.eq.0.or.kc.eq.12.or.kc.eq.14.or.kc.eq.16.or.
     &      kc.eq.18) goto 300
            if(mstu(41).ge.3.and.kchg(kc,2).eq.0.and.jamchge(k(i,2))
     &      .eq.0) goto 300
          endif
          is=2d0-sign(0.5d0,p(i,3))
          ns(is)=ns(is)+1
          pts(is)=pts(is)+sqrt(p(i,1)**2+p(i,2)**2)
  300   continue
        if(ns(1)*pts(2)**2.lt.ns(2)*pts(1)**2)
     &  call pjrobo(1,n+mstu(3),paru(1),0d0,0d0,0d0,0d0)
 
C...Rotate to put second largest jet into -z,+x quadrant.
        do 310 i=1,n
          if(p(i,3).ge.0d0) goto 310
          if(k(i,1).le.0.or.k(i,1).gt.10) goto 310
          if(mstu(41).ge.2) then
            kc=jamcomp(k(i,2))
            if(kc.eq.0.or.kc.eq.12.or.kc.eq.14.or.kc.eq.16.or.
     &      kc.eq.18) goto 310
            if(mstu(41).ge.3.and.kchg(kc,2).eq.0.and.jamchge(k(i,2))
     &      .eq.0) goto 310
          endif
          is=2d0-sign(0.5d0,p(i,1))
          pls(is)=pls(is)-p(i,3)
  310   continue
        if(pls(2).gt.pls(1)) call pjrobo(1,n+mstu(3),0d0,paru(1),
     &  0d0,0d0,0d0)
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYLIST
C...Gives program heading, or lists an event, or particle
C...data, or current parameter values.
 
      subroutine pjlist(mlist)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Parameter statement to help give large particle numbers.
      parameter (ksusy1=1000000,ksusy2=2000000,kexcit=4000000)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/jydat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
      save /jyjets/,/jydat1/,/jydat2/,/jydat3/
C...Local arrays, character variables and data.
      character chap*16,chac*16,chan*16,chad(5)*16,chdl(7)*4
      dimension ps(6)
      data chdl/'(())',' ','()','!!','<>','==','(==)'/
 
C...Initialization printout: version number and date of last change.
      if(mlist.eq.0.or.mstu(12).eq.1) then
        call pjlogo
        mstu(12)=0
        if(mlist.eq.0) return
      endif
 
C...List event data, including additional lines after N.
      if(mlist.ge.1.and.mlist.le.3) then
        if(mlist.eq.1) write(mstu(11),5100)
        if(mlist.eq.2) write(mstu(11),5200)
        if(mlist.eq.3) write(mstu(11),5300)
        lmx=12
        if(mlist.ge.2) lmx=16
        istr=0
        imax=n
        if(mstu(2).gt.0) imax=mstu(2)
        do 120 i=max(1,mstu(1)),max(imax,n+max(0,mstu(3)))
          if((i.gt.imax.and.i.le.n).or.k(i,1).lt.0) goto 120
 
C...Get particle name, pad it and check it is not too long.
          call pjname(k(i,2),chap)
          len=0
          do 100 lem=1,16
            if(chap(lem:lem).ne.' ') len=lem
  100     continue
          mdl=(k(i,1)+19)/10
          ldl=0
          if(mdl.eq.2.or.mdl.ge.8) then
            chac=chap
            if(len.gt.lmx) chac(lmx:lmx)='?'
          else
            ldl=1
            if(mdl.eq.1.or.mdl.eq.7) ldl=2
            if(len.eq.0) then
              chac=chdl(mdl)(1:2*ldl)//' '
            else
              chac=chdl(mdl)(1:ldl)//chap(1:min(len,lmx-2*ldl))//
     &        chdl(mdl)(ldl+1:2*ldl)//' '
              if(len+2*ldl.gt.lmx) chac(lmx:lmx)='?'
            endif
          endif
 
C...Add information on string connection.
          if(k(i,1).eq.1.or.k(i,1).eq.2.or.k(i,1).eq.11.or.k(i,1).eq.12)
     &    then
            kc=jamcomp(k(i,2))
            kcc=0
            if(kc.ne.0) kcc=kchg(kc,2)
            if(iabs(k(i,2)).eq.39) then
              if(len+2*ldl+3.le.lmx) chac(lmx-1:lmx-1)='X'
            elseif(kcc.ne.0.and.istr.eq.0) then
              istr=1
              if(len+2*ldl+3.le.lmx) chac(lmx-1:lmx-1)='A'
            elseif(kcc.ne.0.and.(k(i,1).eq.2.or.k(i,1).eq.12)) then
              if(len+2*ldl+3.le.lmx) chac(lmx-1:lmx-1)='I'
            elseif(kcc.ne.0) then
              istr=0
              if(len+2*ldl+3.le.lmx) chac(lmx-1:lmx-1)='V'
            endif
          endif
 
C...Write data for particle/jet.
          if(mlist.eq.1.and.abs(p(i,4)).lt.9999d0) then
            write(mstu(11),5400) i,chac(1:12),(k(i,j1),j1=1,3),
     &      (p(i,j2),j2=1,5)
          elseif(mlist.eq.1.and.abs(p(i,4)).lt.99999d0) then
            write(mstu(11),5500) i,chac(1:12),(k(i,j1),j1=1,3),
     &      (p(i,j2),j2=1,5)
          elseif(mlist.eq.1) then
            write(mstu(11),5600) i,chac(1:12),(k(i,j1),j1=1,3),
     &      (p(i,j2),j2=1,5)
          elseif(mstu(5).eq.10000.and.(k(i,1).eq.3.or.k(i,1).eq.13.or.
     &      k(i,1).eq.14)) then
            write(mstu(11),5700) i,chac,(k(i,j1),j1=1,3),
     &      k(i,4)/100000000,mod(k(i,4)/10000,10000),mod(k(i,4),10000),
     &      k(i,5)/100000000,mod(k(i,5)/10000,10000),mod(k(i,5),10000),
     &      (p(i,j2),j2=1,5)
          else
            write(mstu(11),5800) i,chac,(k(i,j1),j1=1,5),
     &      (p(i,j2),j2=1,5)
          endif
          if(mlist.eq.3) write(mstu(11),5900) (v(i,j),j=1,5)
 
C...Insert extra separator lines specified by user.
          if(mstu(70).ge.1) then
            isep=0
            do 110 j=1,min(10,mstu(70))
              if(i.eq.mstu(70+j)) isep=1
  110       continue
            if(isep.eq.1.and.mlist.eq.1) write(mstu(11),6000)
            if(isep.eq.1.and.mlist.ge.2) write(mstu(11),6100)
          endif
  120   continue
 
C...Sum of charges and momenta.
        do 130 j=1,6
          ps(j)=pjp(0,j)
  130   continue
        if(mlist.eq.1.and.abs(ps(4)).lt.9999d0) then
          write(mstu(11),6200) ps(6),(ps(j),j=1,5)
        elseif(mlist.eq.1.and.abs(ps(4)).lt.99999d0) then
          write(mstu(11),6300) ps(6),(ps(j),j=1,5)
        elseif(mlist.eq.1) then
          write(mstu(11),6400) ps(6),(ps(j),j=1,5)
        else
          write(mstu(11),6500) ps(6),(ps(j),j=1,5)
        endif
 
C...Give simple list of KF codes defined in program.
      elseif(mlist.eq.11) then
        write(mstu(11),6600)
        do 140 kf=1,80
          call pjname(kf,chap)
          call pjname(-kf,chan)
          if(chap.ne.' '.and.chan.eq.' ') write(mstu(11),6700) kf,chap
          if(chan.ne.' ') write(mstu(11),6700) kf,chap,-kf,chan
  140   continue
        do 170 kfls=1,3,2
          do 160 kfla=1,5
            do 150 kflb=1,kfla-(3-kfls)/2
              kf=1000*kfla+100*kflb+kfls
              call pjname(kf,chap)
              call pjname(-kf,chan)
              write(mstu(11),6700) kf,chap,-kf,chan
  150       continue
  160     continue
  170   continue
        kf=130
        call pjname(kf,chap)
        write(mstu(11),6700) kf,chap
        kf=310
        call pjname(kf,chap)
        write(mstu(11),6700) kf,chap
        do 200 kmul=0,5
          kfls=3
          if(kmul.eq.0.or.kmul.eq.3) kfls=1
          if(kmul.eq.5) kfls=5
          kflr=0
          if(kmul.eq.2.or.kmul.eq.3) kflr=1
          if(kmul.eq.4) kflr=2
          do 190 kflb=1,5
            do 180 kflc=1,kflb-1
              kf=10000*kflr+100*kflb+10*kflc+kfls
              call pjname(kf,chap)
              call pjname(-kf,chan)
              write(mstu(11),6700) kf,chap,-kf,chan
  180       continue
            kf=10000*kflr+110*kflb+kfls
            call pjname(kf,chap)
            write(mstu(11),6700) kf,chap
  190     continue
  200   continue
        kf=100443
        call pjname(kf,chap)
        write(mstu(11),6700) kf,chap
        kf=100553
        call pjname(kf,chap)
        write(mstu(11),6700) kf,chap
        do 240 kflsp=1,3
          kfls=2+2*(kflsp/3)
          do 230 kfla=1,5
            do 220 kflb=1,kfla
              do 210 kflc=1,kflb
                if(kflsp.eq.1.and.(kfla.eq.kflb.or.kflb.eq.kflc))
     &          goto 210
                if(kflsp.eq.2.and.kfla.eq.kflc) goto 210
                if(kflsp.eq.1) kf=1000*kfla+100*kflc+10*kflb+kfls
                if(kflsp.ge.2) kf=1000*kfla+100*kflb+10*kflc+kfls
                call pjname(kf,chap)
                call pjname(-kf,chan)
                write(mstu(11),6700) kf,chap,-kf,chan
  210         continue
  220       continue
  230     continue
  240   continue
        do 250 kf=ksusy1+1,ksusy1+40
          call pjname(kf,chap)
          call pjname(-kf,chan)
          if(chap.ne.' '.and.chan.eq.' ') write(mstu(11),6700) kf,chap
          if(chan.ne.' ') write(mstu(11),6700) kf,chap,-kf,chan
  250   continue
        do 260 kf=ksusy2+1,ksusy2+40
          call pjname(kf,chap)
          call pjname(-kf,chan)
          if(chap.ne.' '.and.chan.eq.' ') write(mstu(11),6700) kf,chap
          if(chan.ne.' ') write(mstu(11),6700) kf,chap,-kf,chan
  260   continue
        do 270 kf=kexcit+1,kexcit+40
          call pjname(kf,chap)
          call pjname(-kf,chan)
          if(chap.ne.' '.and.chan.eq.' ') write(mstu(11),6700) kf,chap
          if(chan.ne.' ') write(mstu(11),6700) kf,chap,-kf,chan
  270   continue
 
C...List parton/particle data table. Check whether to be listed.
      elseif(mlist.eq.12) then
        write(mstu(11),6800)
        do 300 kc=1,mstu(6)
          kf=kchg(kc,4)
          if(kf.eq.0) goto 300
          if(kf.lt.mstu(1).or.(mstu(2).gt.0.and.kf.gt.mstu(2)))
     &    goto 300
 
C...Find particle name and mass. Print information.
          call pjname(kf,chap)
          if(kf.le.100.and.chap.eq.' '.and.mdcy(kc,2).eq.0) goto 300
          call pjname(-kf,chan)
          write(mstu(11),6900) kf,kc,chap,chan,(kchg(kc,j1),j1=1,3),
     &    (pmas(kc,j2),j2=1,4),mdcy(kc,1)
 
C...Particle decay: channel number, branching ratios, matrix element,
C...decay products.
          do 290 idc=mdcy(kc,2),mdcy(kc,2)+mdcy(kc,3)-1
            do 280 j=1,5
              call pjname(kfdp(idc,j),chad(j))
  280       continue
            write(mstu(11),7000) idc,mdme(idc,1),mdme(idc,2),brat(idc),
     &      (chad(j),j=1,5)
  290     continue
  300   continue
 
C...List parameter value table.
      elseif(mlist.eq.13) then
        write(mstu(11),7100)
        do 310 i=1,200
          write(mstu(11),7200) i,mstu(i),paru(i),mstj(i),parj(i),parf(i)
  310   continue
      endif
 
C...Format statements for output on unit MSTU(11) (by default 6).
 5100 format(///28x,'Event listing (summary)'//4x,'I particle/jet KS',
     &5x,'KF  orig    p_x      p_y      p_z       E        m'/)
 5200 format(///28x,'Event listing (standard)'//4x,'I  particle/jet',
     &'  K(I,1)   K(I,2) K(I,3)     K(I,4)      K(I,5)       P(I,1)',
     &'       P(I,2)       P(I,3)       P(I,4)       P(I,5)'/)
 5300 format(///28x,'Event listing (with vertices)'//4x,'I  particle/j',
     &'et  K(I,1)   K(I,2) K(I,3)     K(I,4)      K(I,5)       P(I,1)',
     &'       P(I,2)       P(I,3)       P(I,4)       P(I,5)'/73x,
     &'V(I,1)       V(I,2)       V(I,3)       V(I,4)       V(I,5)'/)
 5400 format(1x,i4,1x,a12,1x,i2,i8,1x,i4,5f9.3)
 5500 format(1x,i4,1x,a12,1x,i2,i8,1x,i4,5f9.2)
 5600 format(1x,i4,1x,a12,1x,i2,i8,1x,i4,5f9.1)
 5700 format(1x,i4,2x,a16,1x,i3,1x,i9,1x,i4,2(3x,i1,2i4),5f13.5)
 5800 format(1x,i4,2x,a16,1x,i3,1x,i9,1x,i4,2(3x,i9),5f13.5)
 5900 format(66x,5(1x,f12.3))
 6000 format(1x,78('='))
 6100 format(1x,130('='))
 6200 format(19x,'sum:',f6.2,5x,5f9.3)
 6300 format(19x,'sum:',f6.2,5x,5f9.2)
 6400 format(19x,'sum:',f6.2,5x,5f9.1)
 6500 format(19x,'sum charge:',f6.2,3x,'sum momentum and inv. mass:',
     &5f13.5)
 6600 format(///20x,'List of KF codes in program'/)
 6700 format(4x,i9,4x,a16,6x,i9,4x,a16)
 6800 format(///30x,'Particle/parton data table'//8x,'KF',5x,'KC',4x,
     &'particle',8x,'antiparticle',6x,'chg  col  anti',8x,'mass',7x,
     &'width',7x,'w-cut',5x,'lifetime',1x,'decay'/11x,'IDC',1x,'on/off',
     &1x,'ME',3x,'Br.rat.',4x,'decay products')
 6900 format(/1x,i9,3x,i4,4x,a16,a16,3i5,1x,f12.5,2(1x,f11.5),
     &1x,1p,e13.5,3x,i2)
 7000 format(10x,i4,2x,i3,2x,i3,2x,f10.6,4x,5a16)
 7100 format(///20x,'Parameter value table'//4x,'I',3x,'MSTU(I)',
     &8x,'PARU(I)',3x,'MSTJ(I)',8x,'PARJ(I)',8x,'PARF(I)')
 7200 format(1x,i4,1x,i9,1x,f14.5,1x,i9,1x,f14.5,1x,f14.5)
 
      return
      end
 
C*********************************************************************
 
C...PYLOGO
C...Writes a logo for the program.
 
      subroutine pjlogo
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Parameter for length of information block.
      parameter (irefer=17)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      save /jydat1/,/pjpars/
C...Local arrays and character variables.
      integer idati(6)
      character month(12)*3, logo(48)*32, refer(2*irefer)*36, line*79,
     &vers*1, subv*3, date*2, year*4, hour*2, minu*2, seco*2
 
C...Data on months, logo, titles, and references.
      data month/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',
     &'Oct','Nov','Dec'/
      data (logo(j),j=1,19)/
     &'            *......*            ',
     &'       *:::!!:::::::::::*       ',
     &'    *::::::!!::::::::::::::*    ',
     &'  *::::::::!!::::::::::::::::*  ',
     &' *:::::::::!!:::::::::::::::::* ',
     &' *:::::::::!!:::::::::::::::::* ',
     &'  *::::::::!!::::::::::::::::*! ',
     &'    *::::::!!::::::::::::::* !! ',
     &'    !! *:::!!:::::::::::*    !! ',
     &'    !!     !* -><- *         !! ',
     &'    !!     !!                !! ',
     &'    !!     !!                !! ',
     &'    !!                       !! ',
     &'    !!        ep             !! ',
     &'    !!                       !! ',
     &'    !!                 pp    !! ',
     &'    !!   e+e-                !! ',
     &'    !!                       !! ',
     &'    !!                          '/
      data (logo(j),j=20,38)/
     &'Welcome to the Lund Monte Carlo!',
     &'                                ',
     &'PPP  Y   Y TTTTT H   H III   A  ',
     &'P  P  Y Y    T   H   H  I   A A ',
     &'PPP    Y     T   HHHHH  I  AAAAA',
     &'P      Y     T   H   H  I  A   A',
     &'P      Y     T   H   H III A   A',
     &'                                ',
     &'This is PYTHIA version x.xxx    ',
     &'Last date of change: xx xxx 199x',
     &'                                ',
     &'Now is xx xxx 199x at xx:xx:xx  ',
     &'                                ',
     &'Disclaimer: this program comes  ',
     &'without any guarantees. Beware  ',
     &'of errors and use common sense  ',
     &'when interpreting results.      ',
     &'                                ',
     &'Copyright T. Sjostrand (1997)   '/
      data (refer(j),j=1,18)/
     &'An archive of program versions and d',
     &'ocumentation is found on the web:   ',
     &'http://www.thep.lu.se/tf2/staff/torb',
     &'jorn/Pythia.html                    ',
     &'                                    ',
     &'                                    ',
     &'When you cite this program, currentl',
     &'y the official reference is         ',
     &'T. Sjostrand, Computer Physics Commu',
     &'n. 82 (1994) 74.                    ',
     &'The supersymmetry extensions are des',
     &'cribed in                           ',
     &'S. Mrenna, Computer Physics Commun. ',
     &'101 (1997) 232                      ',
     &'Also remember that the program, to a',
     &' large extent, represents original  ',
     &'physics research. Other publications',
     &' of special relevance to your       '/
      data (refer(j),j=19,2*irefer)/
     &'studies may therefore deserve separa',
     &'te mention.                         ',
     &'                                    ',
     &'                                    ',
     &'Main author: Torbjorn Sjostrand; Dep',
     &'artment of Theoretical Physics 2,   ',
     &'  Lund University, Solvegatan 14A, S',
     &'-223 62 Lund, Sweden;               ',
     &'  phone: + 46 - 46 - 222 48 16; e-ma',
     &'il: torbjorn@thep.lu.se             ',
     &'SUSY author: Stephen Mrenna, Argonne',
     &' National Laboratory,               ',
     &'  9700 South Cass Avenue, Argonne, I',
     &'L 60439, USA;                       ',
     &'  phone: + 1 - 630 - 252 - 7615; e-m',
     &'ail: mrenna@hep.anl.gov             '/
 
C...Check that PYDATA linked.
      if(mstp(183)/10.ne.199) then
        write(mstu(11),'(1X,A)')
     &  'Error: PYDATA has not been linked.'
        write(mstu(11),'(1X,A)') 'Execution stopped!'
        stop
 
C...Write current version number and current date+time.
      else
        write(vers,'(I1)') mstp(181)
        logo(28)(24:24)=vers
        write(subv,'(I3)') mstp(182)
        logo(28)(26:28)=subv
        if(mstp(182).lt.100) logo(28)(26:26)='0'
        write(date,'(I2)') mstp(185)
        logo(29)(22:23)=date
        logo(29)(25:27)=month(mstp(184))
        write(year,'(I4)') mstp(183)
        logo(29)(29:32)=year
        call pjtime(idati)
        if(idati(1).le.0) then
          logo(31)='                                '
        else
          write(date,'(I2)') idati(3)
          logo(31)(8:9)=date
          logo(31)(11:13)=month(max(1,min(12,idati(2))))
          write(year,'(I4)') idati(1)
          logo(31)(15:18)=year
          write(hour,'(I2)') idati(4)
          logo(31)(23:24)=hour
          write(minu,'(I2)') idati(5)
          logo(31)(26:27)=minu
          if(idati(5).lt.10) logo(31)(26:26)='0'
          write(seco,'(I2)') idati(6)
          logo(31)(29:30)=seco
          if(idati(6).lt.10) logo(31)(29:29)='0'
        endif
      endif
 
C...Loop over lines in header. Define page feed and side borders.
      do 100 ilin=1,29+irefer
        line=' '
        if(ilin.eq.1) then
          line(1:1)='1'
        else
          line(2:3)='**'
          line(78:79)='**'
        endif
 
C...Separator lines and logos.
        if(ilin.eq.2.or.ilin.eq.3.or.ilin.ge.28+irefer) then
          line(4:77)='***********************************************'//
     &    '***************************'
        elseif(ilin.ge.6.and.ilin.le.24) then
          line(6:37)=logo(ilin-5)
          line(44:75)=logo(ilin+14)
        elseif(ilin.ge.26.and.ilin.le.25+irefer) then
          line(5:40)=refer(2*ilin-51)
          line(41:76)=refer(2*ilin-50)
        endif
 
C...Write lines to appropriate unit.
        write(mstu(11),'(A79)') line
  100 continue
 
      return
      end
 
C*********************************************************************
 
C...PYUPDA
C...Facilitates the updating of particle and decay data
C...by allowing it to be done in an external file.
 
      subroutine pjupda(mupda,lfn)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/jydat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
      common/jydat4/chaf(500,2)
      character chaf*16
      common/pjint4/mwid(500),wids(500,5)
      save /jydat1/,/jydat2/,/jydat3/,/jydat4/,/pjint4/
C...Local arrays, character variables and data.
      character chinl*120,chkf*9,chvar(22)*9,chlin*72,
     &chblk(20)*72,chold*16,chtmp*16,chnew*16,chcom*24
      data chvar/ 'KCHG(I,1)','KCHG(I,2)','KCHG(I,3)','KCHG(I,4)',
     &'PMAS(I,1)','PMAS(I,2)','PMAS(I,3)','PMAS(I,4)','MDCY(I,1)',
     &'MDCY(I,2)','MDCY(I,3)','MDME(I,1)','MDME(I,2)','BRAT(I)  ',
     &'KFDP(I,1)','KFDP(I,2)','KFDP(I,3)','KFDP(I,4)','KFDP(I,5)',
     &'CHAF(I,1)','CHAF(I,2)','MWID(I)  '/
 
C...Write header if not yet done.
      if(mstu(12).ge.1) call pjlist(0)
 
C...Write information on file for editing.
      if(mupda.eq.1) then
        do 110 kc=1,500
          if(kchg(kc,4).eq.0) return
          write(lfn,5000) kchg(kc,4),(chaf(kc,j1),j1=1,2),
     &    (kchg(kc,j2),j2=1,3),(pmas(kc,j3),j3=1,4),
     &    mwid(kc),mdcy(kc,1)
          do 100 idc=mdcy(kc,2),mdcy(kc,2)+mdcy(kc,3)-1
            write(lfn,5100) mdme(idc,1),mdme(idc,2),brat(idc),
     &      (kfdp(idc,j),j=1,5)
  100     continue
  110   continue
 
C...Read complete set of information from edited file or
C...read partial set of new or updated information from edited file.
      elseif(mupda.eq.2.or.mupda.eq.3) then
 
C...Reset counters.
        kcc=100
        ndc=0
        chkf='         '
        if(mupda.eq.2) then
          do 120 i=1,mstu(6)
c           kchg(i,4)=0
          do j=1,4
           kchg(i,j)=0
           pmas(i,j)=0.0d0
          end do
          chaf(i,1)=' '
          chaf(i,2)=' '
          mwid(i)=0
          mdcy(i,1)=0
          mdcy(i,2)=0
          mdcy(i,3)=0
  120     continue
        else
          do 130 kc=1,mstu(6)
            if(kc.gt.100.and.kchg(kc,4).gt.100) kcc=kc
            ndc=max(ndc,mdcy(kc,2)+mdcy(kc,3)-1)
  130     continue
        endif
 
C...Begin of loop: read new line; unknown whether particle or
C...decay data.
  140   read(lfn,5200,end=190) chinl
 
C...Identify particle code and whether already defined  (for MUPDA=3).
        if(chinl(2:10).ne.'         ') then
          chkf=chinl(2:10)
          read(chkf,5300) kf
          if(mupda.eq.2) then
            if(kf.le.100) then
              kc=kf
            else
              kcc=kcc+1
              kc=kcc
            endif
          else
            kcrep=0
            if(kf.le.100) then
              kcrep=kf
            else
              do 150 kcr=101,kcc
                if(kchg(kcr,4).eq.kf) kcrep=kcr
  150         continue
            endif
C...Remove duplicate old decay data.
            if(kcrep.ne.0) then
              idcrep=mdcy(kcrep,2)
              ndcrep=mdcy(kcrep,3)
              do 160 i=1,kcc
                if(mdcy(i,2).gt.idcrep) mdcy(i,2)=mdcy(i,2)-ndcrep
  160         continue
              do 180 i=idcrep,ndc-ndcrep
                mdme(i,1)=mdme(i+ndcrep,1)
                mdme(i,2)=mdme(i+ndcrep,2)
                brat(i)=brat(i+ndcrep)
                do 170 j=1,5
                  kfdp(i,j)=kfdp(i+ndcrep,j)
  170           continue
  180         continue
              ndc=ndc-ndcrep
              kc=kcrep
            else
              kcc=kcc+1
              kc=kcc
            endif
          endif
 
C...Study line with particle data.
          if(kc.gt.mstu(6)) call pjerrm(27,
     &    '(PYUPDA:) Particle arrays full by KF ='//chkf)
          read(chinl,5000) kchg(kc,4),(chaf(kc,j1),j1=1,2),
     &    (kchg(kc,j2),j2=1,3),(pmas(kc,j3),j3=1,4),
     &    mwid(kc),mdcy(kc,1)
          mdcy(kc,2)=0
          mdcy(kc,3)=0
 
C...Study line with decay data.
        else
          ndc=ndc+1
          if(ndc.gt.mstu(7)) call pjerrm(27,
     &    '(PYUPDA:) Decay data arrays full by KF ='//chkf)
          if(mdcy(kc,2).eq.0) mdcy(kc,2)=ndc
          mdcy(kc,3)=mdcy(kc,3)+1
          read(chinl,5100) mdme(ndc,1),mdme(ndc,2),brat(ndc),
     &    (kfdp(ndc,j),j=1,5)
        endif
 
C...End of loop; ensure that PYCOMP tables are updated.
        goto 140
  190   continue
        mstu(20)=0
 
C...Perform possible tests that new information is consistent.
        mstj24=mstj(24)
        mstj(24)=0
        do 220 kc=1,mstu(6)
          kf=kchg(kc,4)
          if(kf.eq.0) goto 220
          write(chkf,5300) kf
          if(min(pmas(kc,1),pmas(kc,2),pmas(kc,3),pmas(kc,1)-pmas(kc,3),
     &    pmas(kc,4)).lt.0d0.or.mdcy(kc,3).lt.0) call pjerrm(17,
     &    '(PYUPDA:) Mass/width/life/(# channels) wrong for KF ='//chkf)
          brsum=0d0
          do 210 idc=mdcy(kc,2),mdcy(kc,2)+mdcy(kc,3)-1
            if(mdme(idc,2).gt.80) goto 210
            kq=kchg(kc,1)
            pms=pmas(kc,1)-pmas(kc,3)-parj(64)
            merr=0
            do 200 j=1,5
              kp=kfdp(idc,j)
              if(kp.eq.0.or.kp.eq.81.or.iabs(kp).eq.82) then
                if(kp.eq.81) kq=0
              elseif(jamcomp(kp).eq.0) then
                merr=3
              else
                kq=kq-jamchge(kp)
                pms=pms-pjmass(kp)
c               kpc=jamcomp(kp)
c               pms=pms-pmas(kpc,1)
c               if(mstj(24).gt.0) pms=pms+0.5d0*min(pmas(kpc,2),
c    &          pmas(kpc,3))
              endif
  200       continue
            if(kq.ne.0) merr=max(2,merr)
            if(mwid(kc).eq.0.and.kf.ne.311.and.pms.lt.0d0)
     &      merr=max(1,merr)
            if(merr.eq.3) call pjerrm(17,
     &      '(PYUPDA:) Unknown particle code in decay of KF ='//chkf)
            if(merr.eq.2) call pjerrm(17,
     &      '(PYUPDA:) Charge not conserved in decay of KF ='//chkf)
            if(merr.eq.1) call pjerrm(7,
     &      '(PYUPDA:) Kinematically unallowed decay of KF ='//chkf)
            brsum=brsum+brat(idc)
  210     continue
          write(chtmp,5500) brsum
          if(abs(brsum).gt.0.0005d0.and.abs(brsum-1d0).gt.0.0005d0)
     &    call pjerrm(7,'(PYUPDA:) Sum of branching ratios is '//
     &    chtmp(9:16)//' for KF ='//chkf)
  220   continue
        mstj(24)=mstj24
 
C...Write DATA statements for inclusion in program.
      elseif(mupda.eq.4) then
 
C...Find out how many codes and decay channels are actually used.
        kcc=0
        ndc=0
        do 230 i=1,mstu(6)
          if(kchg(i,4).ne.0) then
            kcc=i
            ndc=max(ndc,mdcy(i,2)+mdcy(i,3)-1)
          endif
  230   continue
 
C...Initialize writing of DATA statements for inclusion in program.
        do 300 ivar=1,22
          ndim=mstu(6)
          if(ivar.ge.12.and.ivar.le.19) ndim=mstu(7)
          nlin=1
          chlin=' '
          chlin(7:35)='DATA ('//chvar(ivar)//',I=   1,    )/'
          llin=35
          chold='START'
 
C...Loop through variables for conversion to characters.
          do 280 idim=1,ndim
            if(ivar.eq.1) write(chtmp,5400) kchg(idim,1)
            if(ivar.eq.2) write(chtmp,5400) kchg(idim,2)
            if(ivar.eq.3) write(chtmp,5400) kchg(idim,3)
            if(ivar.eq.4) write(chtmp,5400) kchg(idim,4)
            if(ivar.eq.5) write(chtmp,5500) pmas(idim,1)
            if(ivar.eq.6) write(chtmp,5500) pmas(idim,2)
            if(ivar.eq.7) write(chtmp,5500) pmas(idim,3)
            if(ivar.eq.8) write(chtmp,5500) pmas(idim,4)
            if(ivar.eq.9) write(chtmp,5400) mdcy(idim,1)
            if(ivar.eq.10) write(chtmp,5400) mdcy(idim,2)
            if(ivar.eq.11) write(chtmp,5400) mdcy(idim,3)
            if(ivar.eq.12) write(chtmp,5400) mdme(idim,1)
            if(ivar.eq.13) write(chtmp,5400) mdme(idim,2)
            if(ivar.eq.14) write(chtmp,5600) brat(idim)
            if(ivar.eq.15) write(chtmp,5400) kfdp(idim,1)
            if(ivar.eq.16) write(chtmp,5400) kfdp(idim,2)
            if(ivar.eq.17) write(chtmp,5400) kfdp(idim,3)
            if(ivar.eq.18) write(chtmp,5400) kfdp(idim,4)
            if(ivar.eq.19) write(chtmp,5400) kfdp(idim,5)
            if(ivar.eq.20) chtmp=chaf(idim,1)
            if(ivar.eq.21) chtmp=chaf(idim,2)
            if(ivar.eq.22) write(chtmp,5400) mwid(idim)
 
C...Replace variables beyond what is properly defined.
            if(ivar.le.4) then
              if(idim.gt.kcc) chtmp='               0'
            elseif(ivar.le.8) then
              if(idim.gt.kcc) chtmp='             0.0'
            elseif(ivar.le.11) then
              if(idim.gt.kcc) chtmp='               0'
            elseif(ivar.le.13) then
              if(idim.gt.ndc) chtmp='               0'
            elseif(ivar.le.14) then
              if(idim.gt.ndc) chtmp='             0.0'
            elseif(ivar.le.19) then
              if(idim.gt.ndc) chtmp='               0'
            elseif(ivar.le.21) then
              if(idim.gt.kcc) chtmp='                '
            else
              if(idim.gt.kcc) chtmp='               0'
            endif
 
C...Length of variable, trailing decimal zeros, quotation marks.
            llow=1
            lhig=1
            do 240 ll=1,16
              if(chtmp(17-ll:17-ll).ne.' ') llow=17-ll
              if(chtmp(ll:ll).ne.' ') lhig=ll
  240       continue
            chnew=chtmp(llow:lhig)//' '
            lnew=1+lhig-llow
            if((ivar.ge.5.and.ivar.le.8).or.ivar.eq.14) then
              lnew=lnew+1
  250         lnew=lnew-1
              if(lnew.ge.2.and.chnew(lnew:lnew).eq.'0') goto 250
              if(chnew(lnew:lnew).eq.'.') lnew=lnew-1
              if(lnew.eq.0) then
                chnew(1:3)='0D0'
                lnew=3
              else
                chnew(lnew+1:lnew+2)='D0'
                lnew=lnew+2
              endif
            elseif(ivar.eq.20.or.ivar.eq.21) then
              do 260 ll=lnew,1,-1
                if(chnew(ll:ll).eq.'''') then
                  chtmp=chnew
                  chnew=chtmp(1:ll)//''''//chtmp(ll+1:11)
                  lnew=lnew+1
                endif
  260         continue
              lnew=min(14,lnew)
              chtmp=chnew
              chnew(1:lnew+2)=''''//chtmp(1:lnew)//''''
              lnew=lnew+2
            endif
 
C...Form composite character string, often including repetition counter.
            if(chnew.ne.chold) then
              nrpt=1
              chold=chnew
              chcom=chnew
              lcom=lnew
            else
              lrpt=lnew+1
              if(nrpt.ge.2) lrpt=lnew+3
              if(nrpt.ge.10) lrpt=lnew+4
              if(nrpt.ge.100) lrpt=lnew+5
              if(nrpt.ge.1000) lrpt=lnew+6
              llin=llin-lrpt
              nrpt=nrpt+1
              write(chtmp,5400) nrpt
              lrpt=1
              if(nrpt.ge.10) lrpt=2
              if(nrpt.ge.100) lrpt=3
              if(nrpt.ge.1000) lrpt=4
              chcom(1:lrpt+1+lnew)=chtmp(17-lrpt:16)//'*'//chnew(1:lnew)
              lcom=lrpt+1+lnew
            endif
 
C...Add characters to end of line, to new line (after storing old line),
C...or to new block of lines (after writing old block).
            if(llin+lcom.le.70) then
              chlin(llin+1:llin+lcom+1)=chcom(1:lcom)//','
              llin=llin+lcom+1
            elseif(nlin.le.19) then
              chlin(llin+1:72)=' '
              chblk(nlin)=chlin
              nlin=nlin+1
              chlin(6:6+lcom+1)='&'//chcom(1:lcom)//','
              llin=6+lcom+1
            else
              chlin(llin:72)='/'//' '
              chblk(nlin)=chlin
              write(chtmp,5400) idim-nrpt
              chblk(1)(30:33)=chtmp(13:16)
              do 270 ilin=1,nlin
                write(lfn,5700) chblk(ilin)
  270         continue
              nlin=1
              chlin=' '
              chlin(7:35+lcom+1)='DATA ('//chvar(ivar)//
     &        ',I=    ,    )/'//chcom(1:lcom)//','
              write(chtmp,5400) idim-nrpt+1
              chlin(25:28)=chtmp(13:16)
              llin=35+lcom+1
            endif
  280     continue
 
C...Write final block of lines.
          chlin(llin:72)='/'//' '
          chblk(nlin)=chlin
          write(chtmp,5400) ndim
          chblk(1)(30:33)=chtmp(13:16)
          do 290 ilin=1,nlin
            write(lfn,5700) chblk(ilin)
  290     continue
  300   continue
      endif
 
C...Formats for reading and writing particle data.
 5000 format(1x,i9,2x,a16,2x,a16,3i3,3f12.5,1p,e13.5,2i3)
 5100 format(10x,2i5,f12.6,5i10)
 5200 format(a120)
 5300 format(i9)
 5400 format(i16)
 5500 format(f16.5)
 5600 format(f16.6)
 5700 format(a72)
 
      return
      end
 
C*********************************************************************
 
C...PYK
C...Provides various integer-valued event related data.
 
      function pjk(i,j)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jyjets/,/jydat1/,/jydat2/
 
C...Default value. For I=0 number of entries, number of stable entries
C...or 3 times total charge.
      pjk=0
      if(i.lt.0.or.i.gt.mstu(4).or.j.le.0) then
      elseif(i.eq.0.and.j.eq.1) then
        pjk=n
      elseif(i.eq.0.and.(j.eq.2.or.j.eq.6)) then
        do 100 i1=1,n
          if(j.eq.2.and.k(i1,1).ge.1.and.k(i1,1).le.10) pjk=pjk+1
          if(j.eq.6.and.k(i1,1).ge.1.and.k(i1,1).le.10) pjk=pjk+
     &    jamchge(k(i1,2))
  100   continue
      elseif(i.eq.0) then
 
C...For I > 0 direct readout of K matrix or charge.
      elseif(j.le.5) then
        pjk=k(i,j)
      elseif(j.eq.6) then
        pjk=jamchge(k(i,2))
 
C...Status (existing/fragmented/decayed), parton/hadron separation.
      elseif(j.le.8) then
        if(k(i,1).ge.1.and.k(i,1).le.10) pjk=1
        if(j.eq.8) pjk=pjk*k(i,2)
      elseif(j.le.12) then
        kfa=iabs(k(i,2))
        kc=jamcomp(kfa)
        kq=0
        if(kc.ne.0) kq=kchg(kc,2)
        if(j.eq.9.and.kc.ne.0.and.kq.ne.0) pjk=k(i,2)
        if(j.eq.10.and.kc.ne.0.and.kq.eq.0) pjk=k(i,2)
        if(j.eq.11) pjk=kc
        if(j.eq.12) pjk=kq*isign(1,k(i,2))
 
C...Heaviest flavour in hadron/diquark.
      elseif(j.eq.13) then
        kfa=iabs(k(i,2))
        pjk=mod(kfa/100,10)*(-1)**mod(kfa/100,10)
        if(kfa.lt.10) pjk=kfa
        if(mod(kfa/1000,10).ne.0) pjk=mod(kfa/1000,10)
        pjk=pjk*isign(1,k(i,2))
 
C...Particle history: generation, ancestor, rank.
      elseif(j.le.15) then
        i2=i
        i1=i
  110   pjk=pjk+1
        i2=i1
        i1=k(i1,3)
        if(i1.gt.0.and.k(i1,1).gt.0.and.k(i1,1).le.20) goto 110
        if(j.eq.15) pjk=i2
      elseif(j.eq.16) then
        kfa=iabs(k(i,2))
        if(k(i,1).le.20.and.((kfa.ge.11.and.kfa.le.20).or.kfa.eq.22.or.
     &  (kfa.gt.100.and.mod(kfa/10,10).ne.0))) then
          i1=i
  120     i2=i1
          i1=k(i1,3)
          if(i1.gt.0) then
            kfam=iabs(k(i1,2))
            ilp=1
            if(kfam.ne.0.and.kfam.le.10) ilp=0
            if(kfam.eq.21.or.kfam.eq.91.or.kfam.eq.92.or.kfam.eq.93)
     &      ilp=0
            if(kfam.gt.100.and.mod(kfam/10,10).eq.0) ilp=0
            if(ilp.eq.1) goto 120
          endif
          if(k(i1,1).eq.12) then
            do 130 i3=i1+1,i2
              if(k(i3,3).eq.k(i2,3).and.k(i3,2).ne.91.and.k(i3,2).ne.92
     &        .and.k(i3,2).ne.93) pjk=pjk+1
  130       continue
          else
            i3=i2
  140       pjk=pjk+1
            i3=i3+1
            if(i3.lt.n.and.k(i3,3).eq.k(i2,3)) goto 140
          endif
        endif
 
C...Particle coming from collapsing jet system or not.
      elseif(j.eq.17) then
        i1=i
  150   pjk=pjk+1
        i3=i1
        i1=k(i1,3)
        i0=max(1,i1)
        kc=jamcomp(k(i0,2))
        if(i1.eq.0.or.k(i0,1).le.0.or.k(i0,1).gt.20.or.kc.eq.0) then
          if(pjk.eq.1) pjk=-1
          if(pjk.gt.1) pjk=0
          return
        endif
        if(kchg(kc,2).eq.0) goto 150
        if(k(i1,1).ne.12) pjk=0
        if(k(i1,1).ne.12) return
        i2=i1
  160   i2=i2+1
        if(i2.lt.n.and.k(i2,1).ne.11) goto 160
        k3m=k(i3-1,3)
        if(k3m.ge.i1.and.k3m.le.i2) pjk=0
        k3p=k(i3+1,3)
        if(i3.lt.n.and.k3p.ge.i1.and.k3p.le.i2) pjk=0
 
C...Number of decay products. Colour flow.
      elseif(j.eq.18) then
        if(k(i,1).eq.11.or.k(i,1).eq.12) pjk=max(0,k(i,5)-k(i,4)+1)
        if(k(i,4).eq.0.or.k(i,5).eq.0) pjk=0
      elseif(j.le.22) then
        if(k(i,1).ne.3.and.k(i,1).ne.13.and.k(i,1).ne.14) return
        if(j.eq.19) pjk=mod(k(i,4)/mstu(5),mstu(5))
        if(j.eq.20) pjk=mod(k(i,5)/mstu(5),mstu(5))
        if(j.eq.21) pjk=mod(k(i,4),mstu(5))
        if(j.eq.22) pjk=mod(k(i,5),mstu(5))
      else
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYP
C...Provides various real-valued event related data.
 
      function pjp(i,j)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jyjets/,/jydat1/,/jydat2/
C...Local array.
      dimension psum(4)
 
C...Set default value. For I = 0 sum of momenta or charges,
C...or invariant mass of system.
      pjp=0d0
      if(i.lt.0.or.i.gt.mstu(4).or.j.le.0) then
      elseif(i.eq.0.and.j.le.4) then
        do 100 i1=1,n
          if(k(i1,1).gt.0.and.k(i1,1).le.10) pjp=pjp+p(i1,j)
  100   continue
      elseif(i.eq.0.and.j.eq.5) then
        do 120 j1=1,4
          psum(j1)=0d0
          do 110 i1=1,n
            if(k(i1,1).gt.0.and.k(i1,1).le.10) psum(j1)=psum(j1)+
     &      p(i1,j1)
  110     continue
  120   continue
        pjp=sqrt(max(0d0,psum(4)**2-psum(1)**2-psum(2)**2-psum(3)**2))
      elseif(i.eq.0.and.j.eq.6) then
        do 130 i1=1,n
          if(k(i1,1).gt.0.and.k(i1,1).le.10)
     $                         pjp=pjp+jamchge(k(i1,2))/3d0
  130   continue
      elseif(i.eq.0) then
 
C...Direct readout of P matrix.
      elseif(j.le.5) then
        pjp=p(i,j)
 
C...Charge, total momentum, transverse momentum, transverse mass.
      elseif(j.le.12) then
        if(j.eq.6) pjp=jamchge(k(i,2))/3d0
        if(j.eq.7.or.j.eq.8) pjp=p(i,1)**2+p(i,2)**2+p(i,3)**2
        if(j.eq.9.or.j.eq.10) pjp=p(i,1)**2+p(i,2)**2
        if(j.eq.11.or.j.eq.12) pjp=p(i,5)**2+p(i,1)**2+p(i,2)**2
        if(j.eq.8.or.j.eq.10.or.j.eq.12) pjp=sqrt(pjp)
 
C...Theta and phi angle in radians or degrees.
      elseif(j.le.16) then
        if(j.le.14) pjp=pjangl(p(i,3),sqrt(p(i,1)**2+p(i,2)**2))
        if(j.ge.15) pjp=pjangl(p(i,1),p(i,2))
        if(j.eq.14.or.j.eq.16) pjp=pjp*180d0/paru(1)
 
C...True rapidity, rapidity with pion mass, pseudorapidity.
      elseif(j.le.19) then
        pmr=0d0
        if(j.eq.17) pmr=p(i,5)
        if(j.eq.18) pmr=pjmass(211)
        pr=max(1d-20,pmr**2+p(i,1)**2+p(i,2)**2)
        pjp=sign(log(min((sqrt(pr+p(i,3)**2)+abs(p(i,3)))/sqrt(pr),
     &  1d20)),p(i,3))
 
C...Energy and momentum fractions (only to be used in CM frame).
      elseif(j.le.25) then
        if(j.eq.20) pjp=2d0*sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)/paru(21)
        if(j.eq.21) pjp=2d0*p(i,3)/paru(21)
        if(j.eq.22) pjp=2d0*sqrt(p(i,1)**2+p(i,2)**2)/paru(21)
        if(j.eq.23) pjp=2d0*p(i,4)/paru(21)
        if(j.eq.24) pjp=(p(i,4)+p(i,3))/paru(21)
        if(j.eq.25) pjp=(p(i,4)-p(i,3))/paru(21)
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYEEVT
C...Handles the generation of an e+e- annihilation jet event.
 
      subroutine pjeevt(kfl,ecm)
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jyjets/,/jydat1/,/jydat2/
 
C...Check input parameters.
      if(mstu(12).ge.1) call pjlist(0)
      if(kfl.lt.0.or.kfl.gt.8) then
        call pjerrm(16,'(PYEEVT:) called with unknown flavour code')
        if(mstu(21).ge.1) return
      endif
      if(kfl.le.5) ecmmin=parj(127)+2.02d0*parf(100+max(1,kfl))
      if(kfl.ge.6) ecmmin=parj(127)+2.02d0*pmas(kfl,1)
      if(ecm.lt.ecmmin) then
        call pjerrm(16,'(PYEEVT:) called with too small CM energy')
        if(mstu(21).ge.1) return
      endif
 
C...Check consistency of MSTJ options set.
      if(mstj(109).eq.2.and.mstj(110).ne.1) then
        call pjerrm(6,
     &  '(PYEEVT:) MSTJ(109) value requires MSTJ(110) = 1')
        mstj(110)=1
      endif
      if(mstj(109).eq.2.and.mstj(111).ne.0) then
        call pjerrm(6,
     &  '(PYEEVT:) MSTJ(109) value requires MSTJ(111) = 0')
        mstj(111)=0
      endif
 
C...Initialize alpha_strong and total cross-section.
      mstu(111)=mstj(108)
      if(mstj(108).eq.2.and.(mstj(101).eq.0.or.mstj(101).eq.1))
     &mstu(111)=1
      paru(112)=parj(121)
      if(mstu(111).eq.2) paru(112)=parj(122)
      if(mstj(116).gt.0.and.(mstj(116).ge.2.or.abs(ecm-parj(151)).ge.
     &parj(139).or.10*mstj(102)+kfl.ne.mstj(119))) call pjxtee(kfl,ecm,
     &xtot)
      if(mstj(116).ge.3) mstj(116)=1
      parj(171)=0d0
 
C...Add initial e+e- to event record (documentation only).
      ntry=0
  100 ntry=ntry+1
      if(ntry.gt.100) then
        call pjerrm(14,'(PYEEVT:) caught in an infinite loop')
        return
      endif
      mstu(24)=0
      nc=0
      if(mstj(115).ge.2) then
        nc=nc+2
        call pj1ent(nc-1,11,0.5d0*ecm,0d0,0d0)
        k(nc-1,1)=21
        call pj1ent(nc,-11,0.5d0*ecm,paru(1),0d0)
        k(nc,1)=21
      endif
 
C...Radiative photon (in initial state).
      mk=0
      ecmc=ecm
      if(mstj(107).ge.1.and.mstj(116).ge.1) call pjradk(ecm,mk,pak,
     &thek,phik,alpk)
      if(mk.eq.1) ecmc=sqrt(ecm*(ecm-2d0*pak))
      if(mstj(115).ge.1.and.mk.eq.1) then
        nc=nc+1
        call pj1ent(nc,22,pak,thek,phik)
        k(nc,3)=min(mstj(115)/2,1)
      endif
 
C...Virtual exchange boson (gamma or Z0).
      if(mstj(115).ge.3) then
        nc=nc+1
        kf=22
        if(mstj(102).eq.2) kf=23
        mstu10=mstu(10)
        mstu(10)=1
        p(nc,5)=ecmc
        call pj1ent(nc,kf,ecmc,0d0,0d0)
        k(nc,1)=21
        k(nc,3)=1
        mstu(10)=mstu10
      endif
 
C...Choice of flavour and jet configuration.
      call pjxkfl(kfl,ecm,ecmc,kflc)
      if(kflc.eq.0) goto 100
      call pjxjet(ecmc,njet,cut)
      kfln=21
      if(njet.eq.4) call pjx4jt(njet,cut,kflc,ecmc,kfln,x1,x2,x4,
     &x12,x14)
      if(njet.eq.3) call pjx3jt(njet,cut,kflc,ecmc,x1,x3)
      if(njet.eq.2) mstj(120)=1
 
C...Fill jet configuration and origin.
      if(njet.eq.2.and.mstj(101).ne.5) call pj2ent(nc+1,kflc,-kflc,ecmc)
      if(njet.eq.2.and.mstj(101).eq.5) call pj2ent(-(nc+1),kflc,-kflc,
     &ecmc)
      if(njet.eq.3) call pj3ent(nc+1,kflc,21,-kflc,ecmc,x1,x3)
      if(njet.eq.4.and.kfln.eq.21) call pj4ent(nc+1,kflc,kfln,kfln,
     &-kflc,ecmc,x1,x2,x4,x12,x14)
      if(njet.eq.4.and.kfln.ne.21) call pj4ent(nc+1,kflc,-kfln,kfln,
     &-kflc,ecmc,x1,x2,x4,x12,x14)
      if(mstu(24).ne.0) goto 100
      do 110 ip=nc+1,n
        k(ip,3)=k(ip,3)+min(mstj(115)/2,1)+(mstj(115)/3)*(nc-1)
  110 continue
 
C...Angular orientation according to matrix element.
      if(mstj(106).eq.1) then
        call pjxdif(nc,njet,kflc,ecmc,chi,the,phi)
        call pjrobo(nc+1,n,0d0,chi,0d0,0d0,0d0)
        call pjrobo(nc+1,n,the,phi,0d0,0d0,0d0)
      endif
 
C...Rotation and boost from radiative photon.
      if(mk.eq.1) then
        dbek=-pak/(ecm-pak)
        nmin=nc+1-mstj(115)/3
        call pjrobo(nmin,n,0d0,-phik,0d0,0d0,0d0)
        call pjrobo(nmin,n,alpk,0d0,dbek*sin(thek),0d0,dbek*cos(thek))
        call pjrobo(nmin,n,0d0,phik,0d0,0d0,0d0)
      endif
 
C...Generate parton shower. Rearrange along strings and check.
      if(mstj(101).eq.5) then
        call pjshow(n-1,n,ecmc)
        mstj14=mstj(14)
        if(mstj(105).eq.-1) mstj(14)=-1
        if(mstj(105).ge.0) mstu(28)=0
        call pjprep(0)
        mstj(14)=mstj14
        if(mstj(105).ge.0.and.mstu(28).ne.0) goto 100
      endif
 
C...Fragmentation/decay generation. Information for PYTABU.
      if(mstj(105).eq.1) call pjexec
      mstu(161)=kflc
      mstu(162)=-kflc
 
      return
      end
 
C*********************************************************************
 
C...PYXTEE
C...Calculates total cross-section, including initial state
C...radiation effects.
 
      subroutine pjxtee(kfl,ecm,xtot)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jydat1/,/jydat2/
 
C...Status, (optimized) Q^2 scale, alpha_strong.
      parj(151)=ecm
      mstj(119)=10*mstj(102)+kfl
      if(mstj(111).eq.0) then
        q2r=ecm**2
      elseif(mstu(111).eq.0) then
        parj(168)=min(1d0,max(parj(128),exp(-12d0*paru(1)/
     &  ((33d0-2d0*mstu(112))*paru(111)))))
        q2r=parj(168)*ecm**2
      else
        parj(168)=min(1d0,max(parj(128),paru(112)/ecm,
     &  (2d0*paru(112)/ecm)**2))
        q2r=parj(168)*ecm**2
      endif
      alspi=pjalps(q2r)/paru(1)
 
C...QCD corrections factor in R.
      if(mstj(101).eq.0.or.mstj(109).eq.1) then
        rqcd=1d0
      elseif(iabs(mstj(101)).eq.1.and.mstj(109).eq.0) then
        rqcd=1d0+alspi
      elseif(mstj(109).eq.0) then
        rqcd=1d0+alspi+(1.986d0-0.115d0*mstu(118))*alspi**2
        if(mstj(111).eq.1) rqcd=max(1d0,rqcd+(33d0-2d0*mstu(112))/12d0*
     &  log(parj(168))*alspi**2)
      elseif(iabs(mstj(101)).eq.1) then
        rqcd=1d0+(3d0/4d0)*alspi
      else
        rqcd=1d0+(3d0/4d0)*alspi-(3d0/32d0+0.519d0*mstu(118))*alspi**2
      endif
 
C...Calculate Z0 width if default value not acceptable.
      if(mstj(102).ge.3) then
        rva=3d0*(3d0+(4d0*paru(102)-1d0)**2)+6d0*rqcd*(2d0+
     &  (1d0-8d0*paru(102)/3d0)**2+(4d0*paru(102)/3d0-1d0)**2)
        do 100 kflc=5,6
          vq=1d0
          if(mod(mstj(103),2).eq.1) vq=sqrt(max(0d0,1d0-
     &    (2d0*pjmass(kflc)/ ecm)**2))
          if(kflc.eq.5) vf=4d0*paru(102)/3d0-1d0
          if(kflc.eq.6) vf=1d0-8d0*paru(102)/3d0
          rva=rva+3d0*rqcd*(0.5d0*vq*(3d0-vq**2)*vf**2+vq**3)
  100   continue
        parj(124)=paru(101)*parj(123)*rva/(48d0*paru(102)*
     &  (1d0-paru(102)))
      endif
 
C...Calculate propagator and related constants for QFD case.
      poll=1d0-parj(131)*parj(132)
      if(mstj(102).ge.2) then
        sff=1d0/(16d0*paru(102)*(1d0-paru(102)))
        sfw=ecm**4/((ecm**2-parj(123)**2)**2+(parj(123)*parj(124))**2)
        sfi=sfw*(1d0-(parj(123)/ecm)**2)
        ve=4d0*paru(102)-1d0
        sf1i=sff*(ve*poll+parj(132)-parj(131))
        sf1w=sff**2*((ve**2+1d0)*poll+2d0*ve*(parj(132)-parj(131)))
        hf1i=sfi*sf1i
        hf1w=sfw*sf1w
      endif
 
C...Loop over different flavours: charge, velocity.
      rtot=0d0
      rqq=0d0
      rqv=0d0
      rva=0d0
      do 110 kflc=1,max(mstj(104),kfl)
        if(kfl.gt.0.and.kflc.ne.kfl) goto 110
        mstj(93)=1
        pmq=pjmass(kflc)
        if(ecm.lt.2d0*pmq+parj(127)) goto 110
        qf=kchg(kflc,1)/3d0
        vq=1d0
        if(mod(mstj(103),2).eq.1) vq=sqrt(1d0-(2d0*pmq/ecm)**2)
 
C...Calculate R and sum of charges for QED or QFD case.
        rqq=rqq+3d0*qf**2*poll
        if(mstj(102).le.1) then
          rtot=rtot+3d0*0.5d0*vq*(3d0-vq**2)*qf**2*poll
        else
          vf=sign(1d0,qf)-4d0*qf*paru(102)
          rqv=rqv-6d0*qf*vf*sf1i
          rva=rva+3d0*(vf**2+1d0)*sf1w
          rtot=rtot+3d0*(0.5d0*vq*(3d0-vq**2)*(qf**2*poll-
     &    2d0*qf*vf*hf1i+vf**2*hf1w)+vq**3*hf1w)
        endif
  110 continue
      rsum=rqq
      if(mstj(102).ge.2) rsum=rqq+sfi*rqv+sfw*rva
 
C...Calculate cross-section, including QCD corrections.
      parj(141)=rqq
      parj(142)=rtot
      parj(143)=rtot*rqcd
      parj(144)=parj(143)
      parj(145)=parj(141)*86.8d0/ecm**2
      parj(146)=parj(142)*86.8d0/ecm**2
      parj(147)=parj(143)*86.8d0/ecm**2
      parj(148)=parj(147)
      parj(157)=rsum*rqcd
      parj(158)=0d0
      parj(159)=0d0
      xtot=parj(147)
      if(mstj(107).le.0) return
 
C...Virtual cross-section.
      xkl=parj(135)
      xku=min(parj(136),1d0-(2d0*parj(127)/ecm)**2)
      ale=2d0*log(ecm/pjmass(11))-1d0
      sigv=ale/3d0+2d0*log(ecm**2/(pjmass(13)*pjmass(15)))/3d0-4d0/3d0+
     &1.526d0*log(ecm**2/0.932d0)
 
C...Soft and hard radiative cross-section in QED case.
      if(mstj(102).le.1) then
        sigv=1.5d0*ale-0.5d0+paru(1)**2/3d0+2d0*sigv
        sigs=ale*(2d0*log(xkl)-log(1d0-xkl)-xkl)
        sigh=ale*(2d0*log(xku/xkl)-log((1d0-xku)/(1d0-xkl))-(xku-xkl))
 
C...Soft and hard radiative cross-section in QFD case.
      else
        szm=1d0-(parj(123)/ecm)**2
        szw=parj(123)*parj(124)/ecm**2
        parj(161)=-rqq/rsum
        parj(162)=-(rqq+rqv+rva)/rsum
        parj(163)=(rqv*(1d0-0.5d0*szm-sfi)+rva*(1.5d0-szm-sfw))/rsum
        parj(164)=(rqv*szw**2*(1d0-2d0*sfw)+rva*(2d0*sfi+szw**2-
     &  4d0+3d0*szm-szm**2))/(szw*rsum)
        sigv=1.5d0*ale-0.5d0+paru(1)**2/3d0+((2d0*rqq+sfi*rqv)/
     &  rsum)*sigv+(szw*sfw*rqv/rsum)*paru(1)*20d0/9d0
        sigs=ale*(2d0*log(xkl)+parj(161)*log(1d0-xkl)+parj(162)*xkl+
     &  parj(163)*log(((xkl-szm)**2+szw**2)/(szm**2+szw**2))+
     &  parj(164)*(atan((xkl-szm)/szw)-atan(-szm/szw)))
        sigh=ale*(2d0*log(xku/xkl)+parj(161)*log((1d0-xku)/
     &  (1d0-xkl))+parj(162)*(xku-xkl)+parj(163)*
     &  log(((xku-szm)**2+szw**2)/((xkl-szm)**2+szw**2))+
     &  parj(164)*(atan((xku-szm)/szw)-atan((xkl-szm)/szw)))
      endif
 
C...Total cross-section and fraction of hard photon events.
      parj(160)=sigh/(paru(1)/paru(101)+sigv+sigs+sigh)
      parj(157)=rsum*(1d0+(paru(101)/paru(1))*(sigv+sigs+sigh))*rqcd
      parj(144)=parj(157)
      parj(148)=parj(144)*86.8d0/ecm**2
      xtot=parj(148)
 
      return
      end
 
C*********************************************************************
 
C...PYRADK
C...Generates initial state photon radiation.
 
      subroutine pjradk(ecm,mk,pak,thek,phik,alpk)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      save /jydat1/
 
C...Function: cumulative hard photon spectrum in QFD case.
      fxk(xx)=2d0*log(xx)+parj(161)*log(1d0-xx)+parj(162)*xx+
     &parj(163)*log((xx-szm)**2+szw**2)+parj(164)*atan((xx-szm)/szw)
 
C...Determine whether radiative photon or not.
      mk=0
      pak=0d0
      if(parj(160).lt.pjr(0)) return
      mk=1
 
C...Photon energy range. Find photon momentum in QED case.
      xkl=parj(135)
      xku=min(parj(136),1d0-(2d0*parj(127)/ecm)**2)
      if(mstj(102).le.1) then
  100   xk=1d0/(1d0+(1d0/xkl-1d0)*((1d0/xku-1d0)/(1d0/xkl-1d0))**pjr(0))
        if(1d0+(1d0-xk)**2.lt.2d0*pjr(0)) goto 100
 
C...Ditto in QFD case, by numerical inversion of integrated spectrum.
      else
        szm=1d0-(parj(123)/ecm)**2
        szw=parj(123)*parj(124)/ecm**2
        fxkl=fxk(xkl)
        fxku=fxk(xku)
        fxkd=1d-4*(fxku-fxkl)
        fxkr=fxkl+pjr(0)*(fxku-fxkl)
        nxk=0
  110   nxk=nxk+1
        xk=0.5d0*(xkl+xku)
        fxkv=fxk(xk)
        if(fxkv.gt.fxkr) then
          xku=xk
          fxku=fxkv
        else
          xkl=xk
          fxkl=fxkv
        endif
        if(nxk.lt.15.and.fxku-fxkl.gt.fxkd) goto 110
        xk=xkl+(xku-xkl)*(fxkr-fxkl)/(fxku-fxkl)
      endif
      pak=0.5d0*ecm*xk
 
C...Photon polar and azimuthal angle.
      pme=2d0*(pjmass(11)/ecm)**2
  120 cthm=pme*(2d0/pme)**pjr(0)
      if(1d0-(xk**2*cthm*(1d0-0.5d0*cthm)+2d0*(1d0-xk)*pme/max(pme,
     &cthm*(1d0-0.5d0*cthm)))/(1d0+(1d0-xk)**2).lt.pjr(0)) goto 120
      cthe=1d0-cthm
      if(pjr(0).gt.0.5d0) cthe=-cthe
      sthe=sqrt(max(0d0,(cthm-pme)*(2d0-cthm)))
      thek=pjangl(cthe,sthe)
      phik=paru(2)*pjr(0)
 
C...Rotation angle for hadronic system.
      sgn=1d0
      if(0.5d0*(2d0-xk*(1d0-cthe))**2/((2d0-xk)**2+(xk*cthe)**2).gt.
     &pjr(0)) sgn=-1d0
      alpk=asin(sgn*sthe*(xk-sgn*(2d0*sqrt(1d0-xk)-2d0+xk)*cthe)/
     &(2d0-xk*(1d0-sgn*cthe)))
 
      return
      end
 
C*********************************************************************
 
C...PYXKFL
C...Selects flavour for produced qqbar pair.
 
      subroutine pjxkfl(kfl,ecm,ecmc,kflc)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jydat1/,/jydat2/
 
C...Calculate maximum weight in QED or QFD case.
      if(mstj(102).le.1) then
        rfmax=4d0/9d0
      else
        poll=1d0-parj(131)*parj(132)
        sff=1d0/(16d0*paru(102)*(1d0-paru(102)))
        sfw=ecmc**4/((ecmc**2-parj(123)**2)**2+(parj(123)*parj(124))**2)
        sfi=sfw*(1d0-(parj(123)/ecmc)**2)
        ve=4d0*paru(102)-1d0
        hf1i=sfi*sff*(ve*poll+parj(132)-parj(131))
        hf1w=sfw*sff**2*((ve**2+1d0)*poll+2d0*ve*(parj(132)-parj(131)))
        rfmax=max(4d0/9d0*poll-4d0/3d0*(1d0-8d0*paru(102)/3d0)*hf1i+
     &  ((1d0-8d0*paru(102)/3d0)**2+1d0)*hf1w,1d0/9d0*poll+2d0/3d0*
     &  (-1d0+4d0*paru(102)/3d0)*hf1i+((-1d0+4d0*paru(102)/3d0)**2+
     &  1d0)*hf1w)
      endif
 
C...Choose flavour. Gives charge and velocity.
      ntry=0
  100 ntry=ntry+1
      if(ntry.gt.100) then
        call pjerrm(14,'(PYXKFL:) caught in an infinite loop')
        kflc=0
        return
      endif
      kflc=kfl
      if(kfl.le.0) kflc=1+int(mstj(104)*pjr(0))
      mstj(93)=1
      pmq=pjmass(kflc)
      if(ecm.lt.2d0*pmq+parj(127)) goto 100
      qf=kchg(kflc,1)/3d0
      vq=1d0
      if(mod(mstj(103),2).eq.1) vq=sqrt(max(0d0,1d0-(2d0*pmq/ecmc)**2))
 
C...Calculate weight in QED or QFD case.
      if(mstj(102).le.1) then
        rf=qf**2
        rfv=0.5d0*vq*(3d0-vq**2)*qf**2
      else
        vf=sign(1d0,qf)-4d0*qf*paru(102)
        rf=qf**2*poll-2d0*qf*vf*hf1i+(vf**2+1d0)*hf1w
        rfv=0.5d0*vq*(3d0-vq**2)*(qf**2*poll-2d0*qf*vf*hf1i+vf**2*hf1w)+
     &  vq**3*hf1w
        if(rfv.gt.0d0) parj(171)=min(1d0,vq**3*hf1w/rfv)
      endif
 
C...Weighting or new event (radiative photon). Cross-section update.
      if(kfl.le.0.and.rf.lt.pjr(0)*rfmax) goto 100
      parj(158)=parj(158)+1d0
      if(ecmc.lt.2d0*pmq+parj(127).or.rfv.lt.pjr(0)*rf) kflc=0
      if(mstj(107).le.0.and.kflc.eq.0) goto 100
      if(kflc.ne.0) parj(159)=parj(159)+1d0
      parj(144)=parj(157)*parj(159)/parj(158)
      parj(148)=parj(144)*86.8d0/ecm**2
 
      return
      end
 
C*********************************************************************
 
C...PYXJET
C...Selects number of jets in matrix element approach.
 
      subroutine pjxjet(ecm,njet,cut)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      save /jydat1/
C...Local array and data.
      dimension zhut(5)
      data zhut/3.0922d0, 6.2291d0, 7.4782d0, 7.8440d0, 8.2560d0/
 
C...Trivial result for two-jets only, including parton shower.
      if(mstj(101).eq.0.or.mstj(101).eq.5) then
        cut=0d0
 
C...QCD and Abelian vector gluon theory: Q^2 for jet rate and R.
      elseif(mstj(109).eq.0.or.mstj(109).eq.2) then
        cf=4d0/3d0
        if(mstj(109).eq.2) cf=1d0
        if(mstj(111).eq.0) then
          q2=ecm**2
          q2r=ecm**2
        elseif(mstu(111).eq.0) then
          parj(169)=min(1d0,parj(129))
          q2=parj(169)*ecm**2
          parj(168)=min(1d0,max(parj(128),exp(-12d0*paru(1)/
     &    ((33d0-2d0*mstu(112))*paru(111)))))
          q2r=parj(168)*ecm**2
        else
          parj(169)=min(1d0,max(parj(129),(2d0*paru(112)/ecm)**2))
          q2=parj(169)*ecm**2
          parj(168)=min(1d0,max(parj(128),paru(112)/ecm,
     &    (2d0*paru(112)/ecm)**2))
          q2r=parj(168)*ecm**2
        endif
 
C...alpha_strong for R and R itself.
        alspi=(3d0/4d0)*cf*pjalps(q2r)/paru(1)
        if(iabs(mstj(101)).eq.1) then
          rqcd=1d0+alspi
        elseif(mstj(109).eq.0) then
          rqcd=1d0+alspi+(1.986d0-0.115d0*mstu(118))*alspi**2
          if(mstj(111).eq.1) rqcd=max(1d0,rqcd+
     &    (33d0-2d0*mstu(112))/12d0*log(parj(168))*alspi**2)
        else
          rqcd=1d0+alspi-(3d0/32d0+0.519d0*mstu(118))*(4d0*alspi/3d0)**2
        endif
 
C...alpha_strong for jet rate. Initial value for y cut.
        alspi=(3d0/4d0)*cf*pjalps(q2)/paru(1)
        cut=max(0.001d0,parj(125),(parj(126)/ecm)**2)
        if(iabs(mstj(101)).le.1.or.(mstj(109).eq.0.and.mstj(111).eq.0))
     &  cut=max(cut,exp(-sqrt(0.75d0/alspi))/2d0)
        if(mstj(110).eq.2) cut=max(0.01d0,min(0.05d0,cut))
 
C...Parametrization of first order three-jet cross-section.
  100   if(mstj(101).eq.0.or.cut.ge.0.25d0) then
          parj(152)=0d0
        else
          parj(152)=(2d0*alspi/3d0)*((3d0-6d0*cut+2d0*log(cut))*
     &    log(cut/(1d0-2d0*cut))+(2.5d0+1.5d0*cut-6.571d0)*
     &    (1d0-3d0*cut)+5.833d0*(1d0-3d0*cut)**2-3.894d0*
     &    (1d0-3d0*cut)**3+1.342d0*(1d0-3d0*cut)**4)/rqcd
          if(mstj(109).eq.2.and.(mstj(101).eq.2.or.mstj(101).le.-2))
     &    parj(152)=0d0
        endif
 
C...Parametrization of second order three-jet cross-section.
        if(iabs(mstj(101)).le.1.or.mstj(101).eq.3.or.mstj(109).eq.2.or.
     &  cut.ge.0.25d0) then
          parj(153)=0d0
        elseif(mstj(110).le.1) then
          ct=log(1d0/cut-2d0)
          parj(153)=alspi**2*ct**2*(2.419d0+0.5989d0*ct+0.6782d0*ct**2-
     &    0.2661d0*ct**3+0.01159d0*ct**4)/rqcd
 
C...Interpolation in second/first order ratio for Zhu parametrization.
        elseif(mstj(110).eq.2) then
          iza=0
          do 110 iy=1,5
            if(abs(cut-0.01d0*iy).lt.0.0001d0) iza=iy
  110     continue
          if(iza.ne.0) then
            zhurat=zhut(iza)
          else
            iz=100d0*cut
            zhurat=zhut(iz)+(100d0*cut-iz)*(zhut(iz+1)-zhut(iz))
          endif
          parj(153)=alspi*parj(152)*zhurat
        endif
 
C...Shift in second order three-jet cross-section with optimized Q^2.
        if(mstj(111).eq.1.and.iabs(mstj(101)).ge.2.and.mstj(101).ne.3
     &  .and.cut.lt.0.25d0) parj(153)=parj(153)+
     &  (33d0-2d0*mstu(112))/12d0*log(parj(169))*alspi*parj(152)
 
C...Parametrization of second order four-jet cross-section.
        if(iabs(mstj(101)).le.1.or.cut.ge.0.125d0) then
          parj(154)=0d0
        else
          ct=log(1d0/cut-5d0)
          if(cut.le.0.018d0) then
            xqqgg=6.349d0-4.330d0*ct+0.8304d0*ct**2
            if(mstj(109).eq.2) xqqgg=(4d0/3d0)**2*(3.035d0-2.091d0*ct+
     &      0.4059d0*ct**2)
            xqqqq=1.25d0*(-0.1080d0+0.01486d0*ct+0.009364d0*ct**2)
            if(mstj(109).eq.2) xqqqq=8d0*xqqqq
          else
            xqqgg=-0.09773d0+0.2959d0*ct-0.2764d0*ct**2+0.08832d0*ct**3
            if(mstj(109).eq.2) xqqgg=(4d0/3d0)**2*(-0.04079d0+
     &      0.1340d0*ct-0.1326d0*ct**2+0.04365d0*ct**3)
            xqqqq=1.25d0*(0.003661d0-0.004888d0*ct-0.001081d0*ct**2+
     &      0.002093d0*ct**3)
            if(mstj(109).eq.2) xqqqq=8d0*xqqqq
          endif
          parj(154)=alspi**2*ct**2*(xqqgg+xqqqq)/rqcd
          parj(155)=xqqqq/(xqqgg+xqqqq)
        endif
 
C...If negative three-jet rate, change y' optimization parameter.
        if(mstj(111).eq.1.and.parj(152)+parj(153).lt.0d0.and.
     &  parj(169).lt.0.99d0) then
          parj(169)=min(1d0,1.2d0*parj(169))
          q2=parj(169)*ecm**2
          alspi=(3d0/4d0)*cf*pjalps(q2)/paru(1)
          goto 100
        endif
 
C...If too high cross-section, use harder cuts, or fail.
        if(parj(152)+parj(153)+parj(154).ge.1) then
          if(mstj(110).eq.2.and.cut.gt.0.0499d0.and.mstj(111).eq.1.and.
     &    parj(169).lt.0.99d0) then
            parj(169)=min(1d0,1.2d0*parj(169))
            q2=parj(169)*ecm**2
            alspi=(3d0/4d0)*cf*pjalps(q2)/paru(1)
            goto 100
          elseif(mstj(110).eq.2.and.cut.gt.0.0499d0) then
            call pjerrm(26,
     &      '(PYXJET:) no allowed y cut value for Zhu parametrization')
          endif
          cut=0.26d0*(4d0*cut)**(parj(152)+parj(153)+
     &    parj(154))**(-1d0/3d0)
          if(mstj(110).eq.2) cut=max(0.01d0,min(0.05d0,cut))
          goto 100
        endif
 
C...Scalar gluon (first order only).
      else
        alspi=pjalps(ecm**2)/paru(1)
        cut=max(0.001d0,parj(125),(parj(126)/ecm)**2,exp(-3d0/alspi))
        parj(152)=0d0
        if(cut.lt.0.25d0) parj(152)=(alspi/3d0)*((1d0-2d0*cut)*
     &  log((1d0-2d0*cut)/cut)+0.5d0*(9d0*cut**2-1d0))
        parj(153)=0d0
        parj(154)=0d0
      endif
 
C...Select number of jets.
      parj(150)=cut
      if(mstj(101).eq.0.or.mstj(101).eq.5) then
        njet=2
      elseif(mstj(101).le.0) then
        njet=min(4,2-mstj(101))
      else
        rnj=pjr(0)
        njet=2
        if(parj(152)+parj(153)+parj(154).gt.rnj) njet=3
        if(parj(154).gt.rnj) njet=4
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYX3JT
C...Selects the kinematical variables of three-jet events.
 
      subroutine pjx3jt(njet,cut,kfl,ecm,x1,x2)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      save /jydat1/
C...Local array.
      dimension zhup(5,12)
 
C...Coefficients of Zhu second order parametrization.
      data ((zhup(ic1,ic2),ic2=1,12),ic1=1,5)/
     &18.29d0,  89.56d0,  4.541d0,  -52.09d0, -109.8d0,  24.90d0,
     &11.63d0,  3.683d0,  17.50d0,0.002440d0, -1.362d0,-0.3537d0,
     &11.42d0,  6.299d0, -22.55d0,  -8.915d0,  59.25d0, -5.855d0,
     &-32.85d0, -1.054d0, -16.90d0,0.006489d0,-0.8156d0,0.01095d0,
     &7.847d0, -3.964d0, -35.83d0,   1.178d0,  29.39d0, 0.2806d0,
     &47.82d0, -12.36d0, -56.72d0, 0.04054d0,-0.4365d0, 0.6062d0,
     &5.441d0, -56.89d0, -50.27d0,   15.13d0,  114.3d0, -18.19d0,
     &97.05d0, -1.890d0, -139.9d0, 0.08153d0,-0.4984d0, 0.9439d0,
     &-17.65d0,  51.44d0, -58.32d0,   70.95d0, -255.7d0, -78.99d0,
     &476.9d0,  29.65d0, -239.3d0,  0.4745d0, -1.174d0,  6.081d0/
 
C...Dilogarithm of x for x<0.5 (x>0.5 obtained by analytic trick).
      dilog(x)=x+x**2/4d0+x**3/9d0+x**4/16d0+x**5/25d0+x**6/36d0+
     &x**7/49d0
 
C...Event type. Mass effect factors and other common constants.
      mstj(120)=2
      mstj(121)=0
      pmq=pjmass(kfl)
      qme=(2d0*pmq/ecm)**2
      if(mstj(109).ne.1) then
        cutl=log(cut)
        cutd=log(1d0/cut-2d0)
        if(mstj(109).eq.0) then
          cf=4d0/3d0
          cn=3d0
          tr=2d0
          wtmx=min(20d0,37d0-6d0*cutd)
          if(mstj(110).eq.2) wtmx=2d0*(7.5d0+80d0*cut)
        else
          cf=1d0
          cn=0d0
          tr=12d0
          wtmx=0d0
        endif
 
C...Alpha_strong and effects of optimized Q^2 scale. Maximum weight.
        als2pi=paru(118)/paru(2)
        wtopt=0d0
        if(mstj(111).eq.1) wtopt=(33d0-2d0*mstu(112))/6d0*
     &  log(parj(169))*als2pi
        wtmax=max(0d0,1d0+wtopt+als2pi*wtmx)
 
C...Choose three-jet events in allowed region.
  100   njet=3
  110   y13l=cutl+cutd*pjr(0)
        y23l=cutl+cutd*pjr(0)
        y13=exp(y13l)
        y23=exp(y23l)
        y12=1d0-y13-y23
        if(y12.le.cut) goto 110
        if(y13**2+y23**2+2d0*y12.le.2d0*pjr(0)) goto 110
 
C...Second order corrections.
        if(mstj(101).eq.2.and.mstj(110).le.1) then
          y12l=log(y12)
          y13m=log(1d0-y13)
          y23m=log(1d0-y23)
          y12m=log(1d0-y12)
          if(y13.le.0.5d0) y13i=dilog(y13)
          if(y13.ge.0.5d0) y13i=1.644934d0-y13l*y13m-dilog(1d0-y13)
          if(y23.le.0.5d0) y23i=dilog(y23)
          if(y23.ge.0.5d0) y23i=1.644934d0-y23l*y23m-dilog(1d0-y23)
          if(y12.le.0.5d0) y12i=dilog(y12)
          if(y12.ge.0.5d0) y12i=1.644934d0-y12l*y12m-dilog(1d0-y12)
          wt1=(y13**2+y23**2+2d0*y12)/(y13*y23)
          wt2=cf*(-2d0*(cutl-y12l)**2-3d0*cutl-1d0+3.289868d0+
     &    2d0*(2d0*cutl-y12l)*cut/y12)+
     &    cn*((cutl-y12l)**2-(cutl-y13l)**2-(cutl-y23l)**2-
     &    11d0*cutl/6d0+67d0/18d0+1.644934d0-(2d0*cutl-y12l)*cut/y12+
     &    (2d0*cutl-y13l)*cut/y13+(2d0*cutl-y23l)*cut/y23)+
     &    tr*(2d0*cutl/3d0-10d0/9d0)+
     &    cf*(y12/(y12+y13)+y12/(y12+y23)+(y12+y23)/y13+(y12+y13)/y23+
     &    y13l*(4d0*y12**2+2d0*y12*y13+4d0*y12*y23+y13*y23)/
     &    (y12+y23)**2+y23l*(4d0*y12**2+2d0*y12*y23+4d0*y12*y13+
     &    y13*y23)/(y12+y13)**2)/wt1+
     &    cn*(y13l*y13/(y12+y23)+y23l*y23/(y12+y13))/wt1+(cn-2d0*cf)*
     &    ((y12**2+(y12+y13)**2)*(y12l*y23l-y12l*y12m-y23l*
     &    y23m+1.644934d0-y12i-y23i)/(y13*y23)+(y12**2+(y12+y23)**2)*
     &    (y12l*y13l-y12l*y12m-y13l*y13m+1.644934d0-y12i-y13i)/
     &    (y13*y23)+(y13**2+y23**2)/(y13*y23*(y13+y23))-
     &    2d0*y12l*y12**2/(y13+y23)**2-4d0*y12l*y12/(y13+y23))/wt1-
     &    cn*(y13l*y23l-y13l*y13m-y23l*y23m+1.644934d0-y13i-y23i)
          if(1d0+wtopt+als2pi*wt2.le.0d0) mstj(121)=1
          if(1d0+wtopt+als2pi*wt2.le.wtmax*pjr(0)) goto 110
          parj(156)=(wtopt+als2pi*wt2)/(1d0+wtopt+als2pi*wt2)
 
        elseif(mstj(101).eq.2.and.mstj(110).eq.2) then
C...Second order corrections; Zhu parametrization of ERT.
          zx=(y23-y13)**2
          zy=1d0-y12
          iza=0
          do 120 iy=1,5
            if(abs(cut-0.01d0*iy).lt.0.0001d0) iza=iy
  120     continue
          if(iza.ne.0) then
            iz=iza
            wt2=zhup(iz,1)+zhup(iz,2)*zx+zhup(iz,3)*zx**2+(zhup(iz,4)+
     &      zhup(iz,5)*zx)*zy+(zhup(iz,6)+zhup(iz,7)*zx)*zy**2+
     &      (zhup(iz,8)+zhup(iz,9)*zx)*zy**3+zhup(iz,10)/(zx-zy**2)+
     &      zhup(iz,11)/(1d0-zy)+zhup(iz,12)/zy
          else
            iz=100d0*cut
            wtl=zhup(iz,1)+zhup(iz,2)*zx+zhup(iz,3)*zx**2+(zhup(iz,4)+
     &      zhup(iz,5)*zx)*zy+(zhup(iz,6)+zhup(iz,7)*zx)*zy**2+
     &      (zhup(iz,8)+zhup(iz,9)*zx)*zy**3+zhup(iz,10)/(zx-zy**2)+
     &      zhup(iz,11)/(1d0-zy)+zhup(iz,12)/zy
            iz=iz+1
            wtu=zhup(iz,1)+zhup(iz,2)*zx+zhup(iz,3)*zx**2+(zhup(iz,4)+
     &      zhup(iz,5)*zx)*zy+(zhup(iz,6)+zhup(iz,7)*zx)*zy**2+
     &      (zhup(iz,8)+zhup(iz,9)*zx)*zy**3+zhup(iz,10)/(zx-zy**2)+
     &      zhup(iz,11)/(1d0-zy)+zhup(iz,12)/zy
            wt2=wtl+(wtu-wtl)*(100d0*cut+1d0-iz)
          endif
          if(1d0+wtopt+2d0*als2pi*wt2.le.0d0) mstj(121)=1
          if(1d0+wtopt+2d0*als2pi*wt2.le.wtmax*pjr(0)) goto 110
          parj(156)=(wtopt+2d0*als2pi*wt2)/(1d0+wtopt+2d0*als2pi*wt2)
        endif
 
C...Impose mass cuts (gives two jets). For fixed jet number new try.
        x1=1d0-y23
        x2=1d0-y13
        x3=1d0-y12
        if(4d0*y23*y13*y12/x3**2.le.qme) njet=2
        if(mod(mstj(103),4).ge.2.and.iabs(mstj(101)).le.1.and.qme*x3+
     &  0.5d0*qme**2+(0.5d0*qme+0.25d0*qme**2)*((1d0-x2)/(1d0-x1)+
     &  (1d0-x1)/(1d0-x2)).gt.(x1**2+x2**2)*pjr(0)) njet=2
        if(mstj(101).eq.-1.and.njet.eq.2) goto 100
 
C...Scalar gluon model (first order only, no mass effects).
      else
  130   njet=3
  140   x3=sqrt(4d0*cut**2+pjr(0)*((1d0-cut)**2-4d0*cut**2))
        if(log((x3-cut)/cut).le.pjr(0)*log((1d0-2d0*cut)/cut)) goto 140
        yd=sign(2d0*cut*((x3-cut)/cut)**pjr(0)-x3,pjr(0)-0.5d0)
        x1=1d0-0.5d0*(x3+yd)
        x2=1d0-0.5d0*(x3-yd)
        if(4d0*(1d0-x1)*(1d0-x2)*(1d0-x3)/x3**2.le.qme) njet=2
        if(mstj(102).ge.2) then
          if(x3**2-2d0*(1d0+x3)*(1d0-x1)*(1d0-x2)*parj(171).lt.
     &    x3**2*pjr(0)) njet=2
        endif
        if(mstj(101).eq.-1.and.njet.eq.2) goto 130
      endif
 
      return
      end
 
C*********************************************************************
 
C...PYX4JT
C...Selects the kinematical variables of four-jet events.
 
      subroutine pjx4jt(njet,cut,kfl,ecm,kfln,x1,x2,x4,x12,x14)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      save /jydat1/
C...Local arrays.
      dimension wta(4),wtb(4),wtc(4),wtd(4),wte(4)
 
C...Common constants. Colour factors for QCD and Abelian gluon theory.
      pmq=pjmass(kfl)
      qme=(2d0*pmq/ecm)**2
      ct=log(1d0/cut-5d0)
      if(mstj(109).eq.0) then
        cf=4d0/3d0
        cn=3d0
        tr=2.5d0
      else
        cf=1d0
        cn=0d0
        tr=15d0
      endif
 
C...Choice of process (qqbargg or qqbarqqbar).
  100 njet=4
      it=1
      if(parj(155).gt.pjr(0)) it=2
      if(mstj(101).le.-3) it=-mstj(101)-2
      if(it.eq.1) wtmx=0.7d0/cut**2
      if(it.eq.1.and.mstj(109).eq.2) wtmx=0.6d0/cut**2
      if(it.eq.2) wtmx=0.1125d0*cf*tr/cut**2
      id=1
 
C...Sample the five kinematical variables (for qqgg preweighted in y34).
  110 y134=3d0*cut+(1d0-6d0*cut)*pjr(0)
      y234=3d0*cut+(1d0-6d0*cut)*pjr(0)
      if(it.eq.1) y34=(1d0-5d0*cut)*exp(-ct*pjr(0))
      if(it.eq.2) y34=cut+(1d0-6d0*cut)*pjr(0)
      if(y34.le.y134+y234-1d0.or.y34.ge.y134*y234) goto 110
      vt=pjr(0)
      cp=cos(paru(1)*pjr(0))
      y14=(y134-y34)*vt
      y13=y134-y14-y34
      vb=y34*(1d0-y134-y234+y34)/((y134-y34)*(y234-y34))
      y24=0.5d0*(y234-y34)*(1d0-4d0*sqrt(max(0d0,vt*(1d0-vt)*
     &vb*(1d0-vb)))*cp-(1d0-2d0*vt)*(1d0-2d0*vb))
      y23=y234-y34-y24
      y12=1d0-y134-y23-y24
      if(min(y12,y13,y14,y23,y24).le.cut) goto 110
      y123=y12+y13+y23
      y124=y12+y14+y24
 
C...Calculate matrix elements for qqgg or qqqq process.
      ic=0
      wttot=0d0
  120 ic=ic+1
      if(it.eq.1) then
        wta(ic)=(y12*y34**2-y13*y24*y34+y14*y23*y34+3d0*y12*y23*y34+
     &  3d0*y12*y14*y34+4d0*y12**2*y34-y13*y23*y24+2d0*y12*y23*y24-
     &  y13*y14*y24-2d0*y12*y13*y24+2d0*y12**2*y24+y14*y23**2+2d0*y12*
     &  y23**2+y14**2*y23+4d0*y12*y14*y23+4d0*y12**2*y23+2d0*y12*y14**2+
     &  2d0*y12*y13*y14+4d0*y12**2*y14+2d0*y12**2*y13+2d0*y12**3)/
     &  (2d0*y13*y134*y234*y24)+(y24*y34+y12*y34+y13*y24-
     &  y14*y23+y12*y13)/(y13*y134**2)+2d0*y23*(1d0-y13)/
     &  (y13*y134*y24)+y34/(2d0*y13*y24)
        wtb(ic)=(y12*y24*y34+y12*y14*y34-y13*y24**2+y13*y14*y24+2d0*y12*
     &  y14*y24)/(y13*y134*y23*y14)+y12*(1d0+y34)*y124/(y134*y234*y14*
     &  y24)-(2d0*y13*y24+y14**2+y13*y23+2d0*y12*y13)/(y13*y134*y14)+
     &  y12*y123*y124/(2d0*y13*y14*y23*y24)
        wtc(ic)=-(5d0*y12*y34**2+2d0*y12*y24*y34+2d0*y12*y23*y34+
     &  2d0*y12*y14*y34+2d0*y12*y13*y34+4d0*y12**2*y34-y13*y24**2+
     &  y14*y23*y24+y13*y23*y24+y13*y14*y24-y12*y14*y24-y13**2*y24-
     &  3d0*y12*y13*y24-y14*y23**2-y14**2*y23+y13*y14*y23-
     &  3d0*y12*y14*y23-y12*y13*y23)/(4d0*y134*y234*y34**2)+
     &  (3d0*y12*y34**2-3d0*y13*y24*y34+3d0*y12*y24*y34+
     &  3d0*y14*y23*y34-y13*y24**2-y12*y23*y34+6d0*y12*y14*y34+
     &  2d0*y12*y13*y34-2d0*y12**2*y34+y14*y23*y24-3d0*y13*y23*y24-
     &  2d0*y13*y14*y24+4d0*y12*y14*y24+2d0*y12*y13*y24+
     &  3d0*y14*y23**2+2d0*y14**2*y23+2d0*y14**2*y12+
     &  2d0*y12**2*y14+6d0*y12*y14*y23-2d0*y12*y13**2-
     &  2d0*y12**2*y13)/(4d0*y13*y134*y234*y34)
        wtc(ic)=wtc(ic)+(2d0*y12*y34**2-2d0*y13*y24*y34+y12*y24*y34+
     &  4d0*y13*y23*y34+4d0*y12*y14*y34+2d0*y12*y13*y34+2d0*y12**2*y34-
     &  y13*y24**2+3d0*y14*y23*y24+4d0*y13*y23*y24-2d0*y13*y14*y24+
     &  4d0*y12*y14*y24+2d0*y12*y13*y24+2d0*y14*y23**2+4d0*y13*y23**2+
     &  2d0*y13*y14*y23+2d0*y12*y14*y23+4d0*y12*y13*y23+2d0*y12*y14**2+
     &  4d0*y12**2*y13+4d0*y12*y13*y14+2d0*y12**2*y14)/
     &  (4d0*y13*y134*y24*y34)-(y12*y34**2-2d0*y14*y24*y34-
     &  2d0*y13*y24*y34-y14*y23*y34+y13*y23*y34+y12*y14*y34+
     &  2d0*y12*y13*y34-2d0*y14**2*y24-4d0*y13*y14*y24-
     &  4d0*y13**2*y24-y14**2*y23-y13**2*y23+y12*y13*y14-
     &  y12*y13**2)/(2d0*y13*y34*y134**2)+(y12*y34**2-
     &  4d0*y14*y24*y34-2d0*y13*y24*y34-2d0*y14*y23*y34-
     &  4d0*y13*y23*y34-4d0*y12*y14*y34-4d0*y12*y13*y34-
     &  2d0*y13*y14*y24+2d0*y13**2*y24+2d0*y14**2*y23-
     &  2d0*y13*y14*y23-y12*y14**2-6d0*y12*y13*y14-
     &  y12*y13**2)/(4d0*y34**2*y134**2)
        wttot=wttot+y34*cf*(cf*wta(ic)+(cf-0.5d0*cn)*wtb(ic)+
     &  cn*wtc(ic))/8d0
      else
        wtd(ic)=(y13*y23*y34+y12*y23*y34-y12**2*y34+y13*y23*y24+2d0*y12*
     &  y23*y24-y14*y23**2+y12*y13*y24+y12*y14*y23+y12*y13*y14)/(y13**2*
     &  y123**2)-(y12*y34**2-y13*y24*y34+y12*y24*y34-y14*y23*y34-y12*
     &  y23*y34-y13*y24**2+y14*y23*y24-y13*y23*y24-y13**2*y24+y14*
     &  y23**2)/(y13**2*y123*y134)+(y13*y14*y12+y34*y14*y12-y34**2*y12+
     &  y13*y14*y24+2d0*y34*y14*y24-y23*y14**2+y34*y13*y24+y34*y23*y14+
     &  y34*y13*y23)/(y13**2*y134**2)-(y34*y12**2-y13*y24*y12+y34*y24*
     &  y12-y23*y14*y12-y34*y14*y12-y13*y24**2+y23*y14*y24-y13*y14*y24-
     &  y13**2*y24+y23*y14**2)/(y13**2*y134*y123)
        wte(ic)=(y12*y34*(y23-y24+y14+y13)+y13*y24**2-y14*y23*y24+y13*
     &  y23*y24+y13*y14*y24+y13**2*y24-y14*y23*(y14+y23+y13))/(y13*y23*
     &  y123*y134)-y12*(y12*y34-y23*y24-y13*y24-y14*y23-y14*y13)/(y13*
     &  y23*y123**2)-(y14+y13)*(y24+y23)*y34/(y13*y23*y134*y234)+
     &  (y12*y34*(y14-y24+y23+y13)+y13*y24**2-y23*y14*y24+y13*y14*y24+
     &  y13*y23*y24+y13**2*y24-y23*y14*(y14+y23+y13))/(y13*y14*y134*
     &  y123)-y34*(y34*y12-y14*y24-y13*y24-y23*y14-y23*y13)/(y13*y14*
     &  y134**2)-(y23+y13)*(y24+y14)*y12/(y13*y14*y123*y124)
        wttot=wttot+cf*(tr*wtd(ic)+(cf-0.5d0*cn)*wte(ic))/16d0
      endif
 
C...Permutations of momenta in matrix element. Weighting.
  130 if(ic.eq.1.or.ic.eq.3.or.id.eq.2.or.id.eq.3) then
        ysav=y13
        y13=y14
        y14=ysav
        ysav=y23
        y23=y24
        y24=ysav
        ysav=y123
        y123=y124
        y124=ysav
      endif
      if(ic.eq.2.or.ic.eq.4.or.id.eq.3.or.id.eq.4) then
        ysav=y13
        y13=y23
        y23=ysav
        ysav=y14
        y14=y24
        y24=ysav
        ysav=y134
        y134=y234
        y234=ysav
      endif
      if(ic.le.3) goto 120
      if(id.eq.1.and.wttot.lt.pjr(0)*wtmx) goto 110
      ic=5
 
C...qqgg events: string configuration and event type.
      if(it.eq.1) then
        if(mstj(109).eq.0.and.id.eq.1) then
          parj(156)=y34*(2d0*(wta(1)+wta(2)+wta(3)+wta(4))+4d0*(wtc(1)+
     &    wtc(2)+wtc(3)+wtc(4)))/(9d0*wttot)
          if(wta(2)+wta(4)+2d0*(wtc(2)+wtc(4)).gt.pjr(0)*(wta(1)+wta(2)+
     &    wta(3)+wta(4)+2d0*(wtc(1)+wtc(2)+wtc(3)+wtc(4)))) id=2
          if(id.eq.2) goto 130
        elseif(mstj(109).eq.2.and.id.eq.1) then
          parj(156)=y34*(wta(1)+wta(2)+wta(3)+wta(4))/(8d0*wttot)
          if(wta(2)+wta(4).gt.pjr(0)*(wta(1)+wta(2)+wta(3)+wta(4))) id=2
          if(id.eq.2) goto 130
        endif
        mstj(120)=3
        if(mstj(109).eq.0.and.0.5d0*y34*(wtc(1)+wtc(2)+wtc(3)+
     &  wtc(4)).gt.pjr(0)*wttot) mstj(120)=4
        kfln=21
 
C...Mass cuts. Kinematical variables out.
        if(y12.le.cut+qme) njet=2
        if(njet.eq.2) goto 150
        q12=0.5d0*(1d0-sqrt(1d0-qme/y12))
        x1=1d0-(1d0-q12)*y234-q12*y134
        x4=1d0-(1d0-q12)*y134-q12*y234
        x2=1d0-y124
        x12=(1d0-q12)*y13+q12*y23
        x14=y12-0.5d0*qme
        if(y134*y234/((1d0-x1)*(1d0-x4)).le.pjr(0)) njet=2
 
C...qqbarqqbar events: string configuration, choose new flavour.
      else
        if(id.eq.1) then
          wtr=pjr(0)*(wtd(1)+wtd(2)+wtd(3)+wtd(4))
          if(wtr.lt.wtd(2)+wtd(3)+wtd(4)) id=2
          if(wtr.lt.wtd(3)+wtd(4)) id=3
          if(wtr.lt.wtd(4)) id=4
          if(id.ge.2) goto 130
        endif
        mstj(120)=5
        parj(156)=cf*tr*(wtd(1)+wtd(2)+wtd(3)+wtd(4))/(16d0*wttot)
  140   kfln=1+int(5d0*pjr(0))
        if(kfln.ne.kfl.and.0.2d0*parj(156).le.pjr(0)) goto 140
        if(kfln.eq.kfl.and.1d0-0.8d0*parj(156).le.pjr(0)) goto 140
        if(kfln.gt.mstj(104)) njet=2
        pmqn=pjmass(kfln)
        qmen=(2d0*pmqn/ecm)**2
 
C...Mass cuts. Kinematical variables out.
        if(y24.le.cut+qme.or.y13.le.1.1d0*qmen) njet=2
        if(njet.eq.2) goto 150
        q24=0.5d0*(1d0-sqrt(1d0-qme/y24))
        q13=0.5d0*(1d0-sqrt(1d0-qmen/y13))
        x1=1d0-(1d0-q24)*y123-q24*y134
        x4=1d0-(1d0-q24)*y134-q24*y123
        x2=1d0-(1d0-q13)*y234-q13*y124
        x12=(1d0-q24)*((1d0-q13)*y14+q13*y34)+q24*((1d0-q13)*y12+
     &  q13*y23)
        x14=y24-0.5d0*qme
        x34=(1d0-q24)*((1d0-q13)*y23+q13*y12)+q24*((1d0-q13)*y34+
     &  q13*y14)
        if(pmq**2+pmqn**2+min(x12,x34)*ecm**2.le.
     &  (parj(127)+pmq+pmqn)**2) njet=2
        if(y123*y134/((1d0-x1)*(1d0-x4)).le.pjr(0)) njet=2
      endif
  150 if(mstj(101).le.-2.and.njet.eq.2) goto 100
 
      return
      end
 
C*********************************************************************
 
C...PYXDIF
C...Gives the angular orientation of events.
 
      subroutine pjxdif(nc,njet,kfl,ecm,chi,the,phi)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jyjets/,/jydat1/,/jydat2/
 
C...Charge. Factors depending on polarization for QED case.
      qf=kchg(kfl,1)/3d0
      poll=1d0-parj(131)*parj(132)
      pold=parj(132)-parj(131)
      if(mstj(102).le.1.or.mstj(109).eq.1) then
        hf1=poll
        hf2=0d0
        hf3=parj(133)**2
        hf4=0d0
 
C...Factors depending on flavour, energy and polarization for QFD case.
      else
        sff=1d0/(16d0*paru(102)*(1d0-paru(102)))
        sfw=ecm**4/((ecm**2-parj(123)**2)**2+(parj(123)*parj(124))**2)
        sfi=sfw*(1d0-(parj(123)/ecm)**2)
        ae=-1d0
        ve=4d0*paru(102)-1d0
        af=sign(1d0,qf)
        vf=af-4d0*qf*paru(102)
        hf1=qf**2*poll-2d0*qf*vf*sfi*sff*(ve*poll-ae*pold)+
     &  (vf**2+af**2)*sfw*sff**2*((ve**2+ae**2)*poll-2d0*ve*ae*pold)
        hf2=-2d0*qf*af*sfi*sff*(ae*poll-ve*pold)+2d0*vf*af*sfw*sff**2*
     &  (2d0*ve*ae*poll-(ve**2+ae**2)*pold)
        hf3=parj(133)**2*(qf**2-2d0*qf*vf*sfi*sff*ve+(vf**2+af**2)*
     &  sfw*sff**2*(ve**2-ae**2))
        hf4=-parj(133)**2*2d0*qf*vf*sfw*(parj(123)*parj(124)/ecm**2)*
     &  sff*ae
      endif
 
C...Mass factor. Differential cross-sections for two-jet events.
      sq2=sqrt(2d0)
      qme=0d0
      if(mstj(103).ge.4.and.iabs(mstj(101)).le.1.and.mstj(102).le.1.and.
     &mstj(109).ne.1) qme=(2d0*pjmass(kfl)/ecm)**2
      if(njet.eq.2) then
        sigu=4d0*sqrt(1d0-qme)
        sigl=2d0*qme*sqrt(1d0-qme)
        sigt=0d0
        sigi=0d0
        siga=0d0
        sigp=4d0
 
C...Kinematical variables. Reduce four-jet event to three-jet one.
      else
        if(njet.eq.3) then
          x1=2d0*p(nc+1,4)/ecm
          x2=2d0*p(nc+3,4)/ecm
        else
          ecmr=p(nc+1,4)+p(nc+4,4)+sqrt((p(nc+2,1)+p(nc+3,1))**2+
     &    (p(nc+2,2)+p(nc+3,2))**2+(p(nc+2,3)+p(nc+3,3))**2)
          x1=2d0*p(nc+1,4)/ecmr
          x2=2d0*p(nc+4,4)/ecmr
        endif
 
C...Differential cross-sections for three-jet (or reduced four-jet).
        xq=(1d0-x1)/(1d0-x2)
        ct12=(x1*x2-2d0*x1-2d0*x2+2d0+qme)/sqrt((x1**2-qme)*(x2**2-qme))
        st12=sqrt(1d0-ct12**2)
        if(mstj(109).ne.1) then
          sigu=2d0*x1**2+x2**2*(1d0+ct12**2)-qme*(3d0+ct12**2-x1-x2)-
     &    qme*x1/xq+0.5d0*qme*((x2**2-qme)*st12**2-2d0*x2)*xq
          sigl=(x2*st12)**2-qme*(3d0-ct12**2-2.5d0*(x1+x2)+x1*x2+qme)+
     &    0.5d0*qme*(x1**2-x1-qme)/xq+0.5d0*qme*((x2**2-qme)*ct12**2-
     &    x2)*xq
          sigt=0.5d0*(x2**2-qme-0.5d0*qme*(x2**2-qme)/xq)*st12**2
          sigi=((1d0-0.5d0*qme*xq)*(x2**2-qme)*st12*ct12+
     &    qme*(1d0-x1-x2+0.5d0*x1*x2+0.5d0*qme)*st12/ct12)/sq2
          siga=x2**2*st12/sq2
          sigp=2d0*(x1**2-x2**2*ct12)
 
C...Differential cross-sect for scalar gluons (no mass effects).
        else
          x3=2d0-x1-x2
          xt=x2*st12
          ct13=sqrt(max(0d0,1d0-(xt/x3)**2))
          sigu=(1d0-parj(171))*(x3**2-0.5d0*xt**2)+
     &    parj(171)*(x3**2-0.5d0*xt**2-4d0*(1d0-x1)*(1d0-x2)**2/x1)
          sigl=(1d0-parj(171))*0.5d0*xt**2+
     &    parj(171)*0.5d0*(1d0-x1)**2*xt**2
          sigt=(1d0-parj(171))*0.25d0*xt**2+
     &    parj(171)*0.25d0*xt**2*(1d0-2d0*x1)
          sigi=-(0.5d0/sq2)*((1d0-parj(171))*xt*x3*ct13+
     &    parj(171)*xt*((1d0-2d0*x1)*x3*ct13-x1*(x1-x2)))
          siga=(0.25d0/sq2)*xt*(2d0*(1d0-x1)-x1*x3)
          sigp=x3**2-2d0*(1d0-x1)*(1d0-x2)/x1
        endif
      endif
 
C...Upper bounds for differential cross-section.
      hf1a=abs(hf1)
      hf2a=abs(hf2)
      hf3a=abs(hf3)
      hf4a=abs(hf4)
      sigmax=(2d0*hf1a+hf3a+hf4a)*abs(sigu)+2d0*(hf1a+hf3a+hf4a)*
     &abs(sigl)+2d0*(hf1a+2d0*hf3a+2d0*hf4a)*abs(sigt)+2d0*sq2*
     &(hf1a+2d0*hf3a+2d0*hf4a)*abs(sigi)+4d0*sq2*hf2a*abs(siga)+
     &2d0*hf2a*abs(sigp)
 
C...Generate angular orientation according to differential cross-sect.
  100 chi=paru(2)*pjr(0)
      cthe=2d0*pjr(0)-1d0
      phi=paru(2)*pjr(0)
      cchi=cos(chi)
      schi=sin(chi)
      c2chi=cos(2d0*chi)
      s2chi=sin(2d0*chi)
      the=acos(cthe)
      sthe=sin(the)
      c2phi=cos(2d0*(phi-parj(134)))
      s2phi=sin(2d0*(phi-parj(134)))
      sig=((1d0+cthe**2)*hf1+sthe**2*(c2phi*hf3-s2phi*hf4))*sigu+
     &2d0*(sthe**2*hf1-sthe**2*(c2phi*hf3-s2phi*hf4))*sigl+
     &2d0*(sthe**2*c2chi*hf1+((1d0+cthe**2)*c2chi*c2phi-2d0*cthe*s2chi*
     &s2phi)*hf3-((1d0+cthe**2)*c2chi*s2phi+2d0*cthe*s2chi*c2phi)*hf4)*
     &sigt-2d0*sq2*(2d0*sthe*cthe*cchi*hf1-2d0*sthe*(cthe*cchi*c2phi-
     &schi*s2phi)*hf3+2d0*sthe*(cthe*cchi*s2phi+schi*c2phi)*hf4)*sigi+
     &4d0*sq2*sthe*cchi*hf2*siga+2d0*cthe*hf2*sigp
      if(sig.lt.sigmax*pjr(0)) goto 100
 
      return
      end
 
C*********************************************************************
 
C...PYONIA
C...Generates Upsilon and toponium decays into three gluons
C...or two gluons and a photon.
 
      subroutine pjonia(kfl,ecm)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /jyjets/,/jydat1/,/jydat2/
 
C...Printout. Check input parameters.
      if(mstu(12).ge.1) call pjlist(0)
      if(kfl.lt.0.or.kfl.gt.8) then
        call pjerrm(16,'(PYONIA:) called with unknown flavour code')
        if(mstu(21).ge.1) return
      endif
      if(ecm.lt.parj(127)+2.02d0*parf(101)) then
        call pjerrm(16,'(PYONIA:) called with too small CM energy')
        if(mstu(21).ge.1) return
      endif
 
C...Initial e+e- and onium state (optional).
      nc=0
      if(mstj(115).ge.2) then
        nc=nc+2
        call pj1ent(nc-1,11,0.5d0*ecm,0d0,0d0)
        k(nc-1,1)=21
        call pj1ent(nc,-11,0.5d0*ecm,paru(1),0d0)
        k(nc,1)=21
      endif
      kflc=iabs(kfl)
      if(mstj(115).ge.3.and.kflc.ge.5) then
        nc=nc+1
        kf=110*kflc+3
        mstu10=mstu(10)
        mstu(10)=1
        p(nc,5)=ecm
        call pj1ent(nc,kf,ecm,0d0,0d0)
        k(nc,1)=21
        k(nc,3)=1
        mstu(10)=mstu10
      endif
 
C...Choose x1 and x2 according to matrix element.
      ntry=0
  100 x1=pjr(0)
      x2=pjr(0)
      x3=2d0-x1-x2
      if(x3.ge.1d0.or.((1d0-x1)/(x2*x3))**2+((1d0-x2)/(x1*x3))**2+
     &((1d0-x3)/(x1*x2))**2.le.2d0*pjr(0)) goto 100
      ntry=ntry+1
      njet=3
      if(mstj(101).le.4) call pj3ent(nc+1,21,21,21,ecm,x1,x3)
      if(mstj(101).ge.5) call pj3ent(-(nc+1),21,21,21,ecm,x1,x3)
 
C...Photon-gluon-gluon events. Small system modifications. Jet origin.
      mstu(111)=mstj(108)
      if(mstj(108).eq.2.and.(mstj(101).eq.0.or.mstj(101).eq.1))
     &mstu(111)=1
      paru(112)=parj(121)
      if(mstu(111).eq.2) paru(112)=parj(122)
      qf=0d0
      if(kflc.ne.0) qf=kchg(kflc,1)/3d0
      rgam=7.2d0*qf**2*paru(101)/pjalps(ecm**2)
      mk=0
      ecmc=ecm
      if(pjr(0).gt.rgam/(1d0+rgam)) then
        if(1d0-max(x1,x2,x3).le.max((parj(126)/ecm)**2,parj(125)))
     &  njet=2
        if(njet.eq.2.and.mstj(101).le.4) call pj2ent(nc+1,21,21,ecm)
        if(njet.eq.2.and.mstj(101).ge.5) call pj2ent(-(nc+1),21,21,ecm)
      else
        mk=1
        ecmc=sqrt(1d0-x1)*ecm
        if(ecmc.lt.2d0*parj(127)) goto 100
        k(nc+1,1)=1
        k(nc+1,2)=22
        k(nc+1,4)=0
        k(nc+1,5)=0
        if(mstj(101).ge.5) k(nc+2,4)=mstu(5)*(nc+3)
        if(mstj(101).ge.5) k(nc+2,5)=mstu(5)*(nc+3)
        if(mstj(101).ge.5) k(nc+3,4)=mstu(5)*(nc+2)
        if(mstj(101).ge.5) k(nc+3,5)=mstu(5)*(nc+2)
        njet=2
        if(ecmc.lt.4d0*parj(127)) then
          mstu10=mstu(10)
          mstu(10)=1
          p(nc+2,5)=ecmc
          call pj1ent(nc+2,83,0.5d0*(x2+x3)*ecm,paru(1),0d0)
          mstu(10)=mstu10
          njet=0
        endif
      endif
      do 110 ip=nc+1,n
        k(ip,3)=k(ip,3)+(mstj(115)/2)+(kflc/5)*(mstj(115)/3)*(nc-1)
  110 continue
 
C...Differential cross-sections. Upper limit for cross-section.
      if(mstj(106).eq.1) then
        sq2=sqrt(2d0)
        hf1=1d0-parj(131)*parj(132)
        hf3=parj(133)**2
        ct13=(x1*x3-2d0*x1-2d0*x3+2d0)/(x1*x3)
        st13=sqrt(1d0-ct13**2)
        sigl=0.5d0*x3**2*((1d0-x2)**2+(1d0-x3)**2)*st13**2
        sigu=(x1*(1d0-x1))**2+(x2*(1d0-x2))**2+(x3*(1d0-x3))**2-sigl
        sigt=0.5d0*sigl
        sigi=(sigl*ct13/st13+0.5d0*x1*x3*(1d0-x2)**2*st13)/sq2
        sigmax=(2d0*hf1+hf3)*abs(sigu)+2d0*(hf1+hf3)*abs(sigl)+2d0*(hf1+
     &  2d0*hf3)*abs(sigt)+2d0*sq2*(hf1+2d0*hf3)*abs(sigi)
 
C...Angular orientation of event.
  120   chi=paru(2)*pjr(0)
        cthe=2d0*pjr(0)-1d0
        phi=paru(2)*pjr(0)
        cchi=cos(chi)
        schi=sin(chi)
        c2chi=cos(2d0*chi)
        s2chi=sin(2d0*chi)
        the=acos(cthe)
        sthe=sin(the)
        c2phi=cos(2d0*(phi-parj(134)))
        s2phi=sin(2d0*(phi-parj(134)))
        sig=((1d0+cthe**2)*hf1+sthe**2*c2phi*hf3)*sigu+2d0*(sthe**2*hf1-
     &  sthe**2*c2phi*hf3)*sigl+2d0*(sthe**2*c2chi*hf1+((1d0+cthe**2)*
     &  c2chi*c2phi-2d0*cthe*s2chi*s2phi)*hf3)*sigt-
     &  2d0*sq2*(2d0*sthe*cthe*cchi*hf1-2d0*sthe*
     &  (cthe*cchi*c2phi-schi*s2phi)*hf3)*sigi
        if(sig.lt.sigmax*pjr(0)) goto 120
        call pjrobo(nc+1,n,0d0,chi,0d0,0d0,0d0)
        call pjrobo(nc+1,n,the,phi,0d0,0d0,0d0)
      endif
 
C...Generate parton shower. Rearrange along strings and check.
      if(mstj(101).ge.5.and.njet.ge.2) then
        call pjshow(nc+mk+1,-njet,ecmc)
        mstj14=mstj(14)
        if(mstj(105).eq.-1) mstj(14)=-1
        if(mstj(105).ge.0) mstu(28)=0
        call pjprep(0)
        mstj(14)=mstj14
        if(mstj(105).ge.0.and.mstu(28).ne.0) goto 100
      endif
 
C...Generate fragmentation. Information for PYTABU:
      if(mstj(105).eq.1) call pjexec
      mstu(161)=110*kflc+3
      mstu(162)=0
 
      return
      end
 
C********************************************************************* 

      function pjr(idum)
      real*8 pjr,rn
      pjr=rn(idum)
      end
