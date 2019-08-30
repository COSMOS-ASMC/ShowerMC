
      subroutine gheini
coff  implicit none
      integer *4 maxmix
      parameter (maxmix=5)
      integer *4 idum
      character *20 name
c
c *** initialization of relevant gheisha variables ***
c *** interface with gheisha8 ***
c *** nve 20-may-1988 cern geneva ***
c
c called by : gpghei, gheish
c origin : f.carminati
c
c
c
      common/gsecti/ aiel(20),aiin(20),aifi(20),aica(20),alam,k0flag
      integer k0flag
      real aiel,aiin,aifi,aica,alam
c
c
c --- gheisha commons ---
c --- initialization flags for various gheisha routines ---
      common /kginit/ kginit(50)
c
      common/consts/ pi,twpi,pibtw,mp,mpi,mmu,mel,mkch,mk0,smp,smpi,
     $               smu,ct,ctkch,ctk0,
     $               ml0,msp,ms0,msm,mx0,mxm,ctl0,ctsp,ctsm,ctx0,ctxm,
     $               rmass(35),rcharg(35)
c
                     real mp,mpi,mmu,mel,mkch,mk0,
     *                    ml0,msp,ms0,msm,mx0,mxm
c
      common/event / nsize,ncur,next,ntot,eve(1200)
c
      common/prntfl/inbcd,newbcd,inbin,newbin,npevt,nevtp,lprt,nprt(10)
                    logical lprt,nprt
c
c --- boundary limits for arguments of intrinsic functions ---
c --- xl denotes lower bound whereas xu denotes upper bound ---
      common /limits/ expxl,expxu
c
c
c --- "nevent" changed to "kevent" in common /curpar/ due to clash ---
c --- with variable "nevent" in geant common ---
c
      common /curpar/ weight(10),ddeltn,ifile,irun,nevt,kevent,shflag,
     $                ithst,ittot,itlst,ifrnd,tofcut,cmom(5),ceng(5),
     $                rs,s,enp(10),np,nm,nn,nr,no,nz,ipa(200),
     $                atno2,zno2
c
      data clight /2.99792458e10/
      newbcd=6 !standard output
      inbcd=5 !terminal
c
      lprt=.false.
c
c --- initialise kginit array ---
      do 20 j=1,50
      kginit(j)=0
20    continue
      call vzero(nprt,10)
c
c --- initialize some cut-off parameters with geant values ---
      tofcut=1.0e+20
      nsize=1200
      k0flag=0
      ceng(3)=0.0 !let gismo apply cuts
      ceng(4)=0.0 !let gismo apply cuts
coff  ceng(3)=cuthad
coff  ceng(4)=cutneu
c
c --- initialize pi, 2*pi, pi/2 and particle parameters ---
      pi=acos(-1.0)
      twpi=2.0*pi
      pibtw=pi/2.0
c *** gamma ***
      call gfpart(1,name,idum,rmass(1),rcharg(1),tlife)
c *** neutrino ***
      call gfpart(4,name,idum,rmass(2),rcharg(2),tlife)
c *** e+ ***
      call gfpart(2,name,idum,rmass(3),rcharg(3),tlife)
c *** e- ***
      call gfpart(3,name,idum,rmass(4),rcharg(4),tlife)
c *** mu+ ***
      call gfpart(5,name,idum,rmass(5),rcharg(5),tlife)
c *** mu- ***
      call gfpart(6,name,idum,rmass(6),rcharg(6),tlife)
c *** pi+ ***
      call gfpart(8,name,idum,rmass(7),rcharg(7),tlife)
      ct=tlife*clight
c *** pi0 ***
      call gfpart(7,name,idum,rmass(8),rcharg(8),tlife)
c *** pi- ***
      call gfpart(9,name,idum,rmass(9),rcharg(9),tlife)
c *** k+ ***
      call gfpart(11,name,idum,rmass(10),rcharg(10),tlife)
     
      ctkch=clight*tlife
c *** k0 short (==> k0) ***
      rmass(11)=0.49772
      rcharg(11)=0.0
      ctk0=2.675
      call gfpart(16,name,idum,rmass(11),rcharg(11),tlife)
      ctk0=clight*tlife
c *** k0 long (==> k0 bar) ***
      rmass(12)=-0.49772
      rcharg(12)=0.0
      call gfpart(10,name,idum,rmass(12),rcharg(12),tlife)
      rmass(12)=-rmass(12)
c *** k- ***
      call gfpart(12,name,idum,rmass(13),rcharg(13),tlife)
c *** p ***
      call gfpart(14,name,idum,rmass(14),rcharg(14),tlife)
c *** p bar ***
      call gfpart(15,name,idum,rmass(15),rcharg(15),tlife)
      rmass(15)=-rmass(15)
c *** n ***
      call gfpart(13,name,idum,rmass(16),rcharg(16),tlife)
c *** n bar ***
      call gfpart(25,name,idum,rmass(17),rcharg(17),tlife)
      rmass(17)=-rmass(17)
c *** l0 ***
      call gfpart(18,name,idum,rmass(18),rcharg(18),tlife)
      ctlo=clight*tlife
c *** l0 bar ***
      call gfpart(26,name,idum,rmass(19),rcharg(19),tlife)
      rmass(19)=-rmass(19)
c *** s+ ***
      call gfpart(19,name,idum,rmass(20),rcharg(20),tlife)
      ctsp=clight*tlife
c *** s0 ***
      call gfpart(20,name,idum,rmass(21),rcharg(21),tlife)
c ctau ?
c *** s- ***
      call gfpart(21,name,idum,rmass(22),rcharg(22),tlife)
      ctsm=clight*tlife
c *** s+ bar ***
      call gfpart(27,name,idum,rmass(23),rcharg(23),tlife)
      rmass(23)=-rmass(23)
c *** s0 bar ***
      call gfpart(28,name,idum,rmass(24),rcharg(24),tlife)
      rmass(24)=-rmass(24)
c *** s- bar ***
      call gfpart(29,name,idum,rmass(25),rcharg(25),tlife)
      rmass(25)=-rmass(25)
c *** xi0 ***
      call gfpart(22,name,idum,rmass(26),rcharg(26),tlife)
      ctx0=clight*tlife
c *** xi- ***
      call gfpart(23,name,idum,rmass(27),rcharg(27),tlife)
      ctxm=clight*tlife
c *** xi0 bar ***
      call gfpart(30,name,idum,rmass(28),rcharg(28),tlife)
      rmass(28)=-rmass(28)
c *** xi- bar ***
      call gfpart(31,name,idum,rmass(29),rcharg(29),tlife)
      rmass(29)=-rmass(29)
c *** deuteron ***
      call gfpart(45,name,idum,rmass(30),rcharg(30),tlife)
c *** triton ***
      call gfpart(46,name,idum,rmass(31),rcharg(31),tlife)
c *** alpha ***
      call gfpart(47,name,idum,rmass(32),rcharg(32),tlife)
c *** omega- ***
      call gfpart(24,name,idum,rmass(33),rcharg(33),tlife)
c *** omega- bar ***
      call gfpart(32,name,idum,rmass(34),rcharg(34),tlife)
      rmass(34)=-rmass(34)
c *** new particle (geantino) ***
      rmass(35)=0.0
      rcharg(35)=0.0
c
c       write(6,1000) (i,rmass(i),rcharg(i),i=1,33),
c     $            ct,ctkch,ctk0,ctl0,ctsp,ctsm,ctx0,ctxm
c 1000 format(' *gheini* === gheisha particle properties ==='/
c     $ '0index',5x,'mass (gev)',5x,'charge'/1h /
c     $ 33(1h ,1x,i3,5x,f11.6,6x,f5.2/),
c     $ '0pi +-  ct = ',g12.5,' k  +-  ct = ',g12.5/
c     $ ' k0     ct = ',g12.5,' l0     ct = ',g12.5/
c     $ ' s+     ct = ',g12.5,' s-     ct = ',g12.5/
c     $ ' x0     ct = ',g12.5,' x-     ct = ',g12.5)
c
      mp=rmass(14)
      mpi=rmass(7)
      mmu=rmass(5)
      mel=rmass(3)
      mkch=rmass(10)
      mk0=rmass(11)
      smp=mp**2
      smpi=mpi**2
      smu=mmu**2
      ml0=rmass(18)
      msp=rmass(20)
      ms0=rmass(21)
      msm=rmass(22)
      mx0=rmass(26)
      mxm=rmass(27)
c
c --- load limits for intrinsic function arguments ---
      expxl = - 82.0
      expxu =   82.0
      return
      end
