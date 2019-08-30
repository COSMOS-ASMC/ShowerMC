C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(1000,5),P(1000,5),V(1000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,7),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,3),BRAT(4000),KFDP(4000,5)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
c     SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYDAT3/,/PYSUBS/,/PYPARS/
      common/pypos1/jqconst(2),kfcq(3),icq(3),icms
      save /pypos1/
C...Local arrays.
      DIMENSION PSUM(5),PINI(6),PFIN(6)

      icms=0
      jqconst(1)=1
      jqconst(2)=1

c     call init
c...Choice of baryon production model.
      mstj(12)=2
c...Form of particle decay.
      mstj(21)=0
c....Hadron formation point.
      mstj(10)=3

c...Baryon resonance production parameter.
      parj(27)=0.0d0


      ymin=-2.0D0
      ymax=2.0D0
      wy=0.1D0
      nymx=(ymax-ymin)/wy
      call vbook1(11,'dN/dy proton',nymx,ymin,ymax)
      call vbook1(12,'dN/dy pi+',nymx,ymin,ymax)
      em=3.15d0

      nev=10000
      do iev=1,nev
        call pj2ent(0,1,2101,em)
        call pjexec
c       call pjlist(1)

        do  i=1,n
c......This is already removed particle.
        if(k(i,1).ge.10.or.k(i,1).eq.0) go to 390
         kf=k(i,2)
         rap=0.5D0*log( max(p(i,4)+p(i,3),1.D-8)/max(p(i,4)-p(i,3), 
     & 1.D-8) )
          if(kf.eq.2212) call vfill1(11,rap,1.D0/wy)
          if(kf.eq.221) call vfill1(12,rap,1.D0/wy)
 390    end do

c....Loop over event
      end do

c...Event weight
      fac=1.d0/dble(nev)

c...Mass distributions.
      do i=1,2
       call vscale(10+i,fac)
       call vprint(10+i,0,0)
      end do

      end

c*******************************************************************
      subroutine init

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)

        mdcy(jamcomp(211),1)=0  ! pi+
        mdcy(jamcomp(-211),1)=0 ! pi-
        mdcy(jamcomp(111),1)=0  ! pi0
        mdcy(jamcomp(311),1)=0  ! k0
        mdcy(jamcomp(-311),1)=0 ! ak0
        mdcy(jamcomp(321),1)=0  ! k+
        mdcy(jamcomp(-321),1)=0 ! k-
        mdcy(jamcomp(411),1)=0    ! D+
        mdcy(jamcomp(-411),1)=0   ! D-
        mdcy(jamcomp(421),1)=0    ! D0
        mdcy(jamcomp(-421),1)=0   ! aD0
        mdcy(jamcomp(221),1)=0    ! eta
        mdcy(jamcomp(331),1)=0    ! eta'
        mdcy(jamcomp(441),1)=0    ! eta_c

        mdcy(jamcomp(310),1)=0
        mdcy(jamcomp(431),1)=0
        mdcy(jamcomp(-431),1)=0
        mdcy(jamcomp(511),1)=0
        mdcy(jamcomp(-511),1)=0
        mdcy(jamcomp(521),1)=0
        mdcy(jamcomp(-521),1)=0
        mdcy(jamcomp(531),1)=0
        mdcy(jamcomp(-531),1)=0
        mdcy(jamcomp(3122),1)=0
        mdcy(jamcomp(-3122),1)=0
        mdcy(jamcomp(3112),1)=0
        mdcy(jamcomp(-3112),1)=0
        mdcy(jamcomp(3212),1)=0
        mdcy(jamcomp(-3212),1)=0
        mdcy(jamcomp(3222),1)=0
        mdcy(jamcomp(-3222),1)=0
        mdcy(jamcomp(3312),1)=0
        mdcy(jamcomp(-3312),1)=0
        mdcy(jamcomp(3322),1)=0
        mdcy(jamcomp(-3322),1)=0
        mdcy(jamcomp(3334),1)=0
        mdcy(jamcomp(-3334),1)=0

      end
