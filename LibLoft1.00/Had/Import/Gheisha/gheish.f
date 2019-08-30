c
      subroutine gheish
     x(imat,ipart,pvec,mprod,nprod,iprod,tprod,pprod,code,stop)
c    integer     imat:input. matter index
c    integer     ipart:input. particle index; giant code
c    real*4      pvec(3): input.  incident particle momentum
c    integer     mprod:input.  max particle number treatable
c    integer     nprod:output.  produced particle number
c    integer     iprod(mprod):output. each particle type produced,giant code
c    real*4      tprod(mprod): output. ? (time?)
c    real*4      pprod(mprod): output. 3 momentum of each produced ptcl.
c    integer     code:  output. interaction type occurred.
c                      13: elestic scattering
c                      20: inelastic scattering
c                      15: inelascit + nuclear fission
c                      18: neutron capture
c 
c    integer     stop:  output. if non 0, the particle stopped.


c
coff  implicit none
      integer *4 maxmix
      parameter (maxmix=5)
      integer *4 nprod,mprod,code,stop
      integer *4 iprod(mprod)
      real *4 pprod(3,mprod)
      real *4 tprod(mprod)
      real *4 destep
      integer *4 imat,ipart
      real *4 pvec(3)
      integer *4 ndum
      real *4 dum(10)
      character *(*) todo
      real *4 arg
      real *4 value
      integer *4 ivalue
      equivalence (value,ivalue)
c
      integer *4 ihadr
      integer *4 ipfis
      real *4 cuthad
      real *4 cutneu
c
c *** main steering for hadron shower development ***
c *** nve 15-jun-1988 cern geneva ***
c
c called by : guhadr (user routine)
c origin : f.carminati, h.fesefeldt
c                       routines : calim  16-sep-1987
c                                  setres 19-aug-1985
c                                  intact 06-oct-1987
      real *4 a,z,dens,radl,absl
      character *20 name
      real *4 wmix(maxmix)
      integer *4 nmix,imix(maxmix)
      integer *4 ind
c
c
c
c
c
c
c
      common/gsecti/ aiel(20),aiin(20),aifi(20),aica(20),alam,k0flag
      integer k0flag
      real aiel,aiin,aifi,aica,alam
c
c --- gheisha commons ---
      common /vecuty/ pv(10,200)
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
c
c
      common /curpar /weight(10),ddeltn,ifile,irun,nevt,kevent,shflag,
     $                ithst,ittot,itlst,ifrnd,tofcut,cmom(5),ceng(5),
     $                rs,s,enp(10),np,nm,nn,nr,no,nz,ipa(200),
     $                atno2,zno2
c
c --- "ipart" changed to "kpart" in common /result/ due to clash ---
c --- with variable "ipart" in geant common ---
c
      common /result/ xend,yend,zend,rca,rce,amas,nch,tof,px,py,pz,
     $                userw,intct,p,en,ek,amasq,deltn,itk,ntk,kpart,ind,
     $                lcalo,icel,sinl,cosl,sinp,cosp,
     $                xold,yold,zold,pold,pxold,pyold,pzold,
     $                xscat,yscat,zscat,pscat,pxscat,pyscat,pzscat
                      real nch,intct
c
c --- "absl(21)" changed to "abslth(21)" in common /mat/ due to clash ---
c --- with variable "absl" in geant common ---
c
      common /mat/ lmat,
     $             den(21),radlth(21),atno(21),zno(21),abslth(21),
     $             cden(21),mden(21),x0den(21),x1den(21),rion(21),
     $             matid(21),matid1(21,24),parmat(21,10),
     $             ifrat,ifrac(21),frac1(21,10),den1(21,10),
     $             atno1(21,10),zno1(21,10)
c
      dimension ipelos(35)
      save ideol
c
c --- random number array                          --
      dimension rndm(1)
c
c --- dimension stmts. for geant/gheisha particle code conversions ---
c --- kipart(i)=gheisha code corresponding to geant   code i ---
c --- ikpart(i)=geant   code corresponding to gheisha code i ---
c
      integer kipart(48),ikpart(35)
c
c --- data stmts. for geant/gheisha particle code conversions ---
c --- kipart(i)=gheisha code corresponding to geant   code i ---
c --- ikpart(i)=geant   code corresponding to gheisha code i ---
c
      data kipart/
     $               1,   3,   4,   2,   5,   6,   8,   7,
     $               9,  12,  10,  13,  16,  14,  15,  11,
     $              35,  18,  20,  21,  22,  26,  27,  33,
     $              17,  19,  23,  24,  25,  28,  29,  34,
     $              35,  35,  35,  35,  35,  35,  35,  35,
     $              35,  35,  35,  35,  30,  31,  32,  35/
c
      data ikpart/
     $               1,   4,   2,   3,   5,   6,   8,   7,
     $               9,  11,  16,  10,  12,  14,  15,  13,
     $              25,  18,  26,  19,  20,  21,  27,  28,
     $              29,  22,  23,  30,  31,  45,  46,  47,
     $              24,  32,  48/
c
c
c --- denote stable particles according to gheisha code ---
c --- stable : gamma, neutrino, electron, proton and heavy fragments ---
c --- when stopping these particles only loose their kinetic energy ---
      data ipelos/
     $             1,   1,   0,   1,   0,   0,   0,   0,
     $             0,   0,   0,   0,   0,   1,   0,   0,
     $             0,   0,   0,   0,   0,   0,   0,   0,
     $             0,   0,   0,   0,   0,   1,   1,   1,
     $             0,   0,   1/
c
c --- lowerbound of kinetic energy bin in n cross-section tables ---
      data teklow /0.0001/
c
c --- kinetic energy to switch from "casn" to "gnslwd" for n cascade ---
      data swtekn /0.05/
c
      data ideol/0/
c
c !low energy hadronic interaction
      data ihadr/1/
      data ipfis/0/
c !let gismo do it
      data cuthad/0.0/
      data cutneu/0.0/
c
c --- set the interaction mechanism to "hadr" ---
c
c --- init output ---
      stop=0
      code=0
      nprod=0
      destep=0.0
c
c --- update x-sections ---
      p=sqrt(pvec(1)**2+pvec(2)**2+pvec(3)**2)
      if(p.ne.0.0) call gpghei(imat,ipart,p,dum)
c
c
 9004 continue
      kpart=kipart(ipart)
      kkpart=kpart
c
c --- transport the track number to gheisha and initialise some numbers
      ntk=0
      intct=0.0
      next=1
      ntot=0
      intval=0
      tof=0.0
c
c --- fill result common for this track with geant values ---
c --- calim code ---
      xend=0.0
      yend=0.0
      zend=0.0
      amas=rmass(kpart)
      nch=rcharg(kpart)
      charge=rcharg(kpart)
      if(p.lt.1.0e-10) then
       px=0.0
       py=0.0
       pz=0.0
       stop=1
      else
      px=pvec(1)/p
      py=pvec(2)/p
      pz=pvec(3)/p
      endif !p.lt.1.0e-10)
c --- setres code ---
      amasq=amas*amas
      en=sqrt(amasq+p*p)
      ek=abs(en-abs(amas))
      enold=en
c
      if(stop.ne.0) then
       call ghstop(ipart,code,stop)
       if(code.eq.5) go  to 9999
       if(ihadr.ne.2) go to 40
       if(ipelos(kpart).eq.0) then
        destep=destep+en !unstable deposit all energy
       else
        destep=destep+ek !stable deposit kin energy only
       endif !ipelos(kpart).eq.0
       go to 9999
      endif !stop.ne.0
      stop=0
c
      sinl=0.0
      cosl=1.0
      sinp=0.0
      cosp=1.0
c
      if (abs(p) .le. 1.0e-10) go to 1
      sinl=pz
      cosl=sqrt(abs(1.0-sinl**2))
c
 1    continue
      call grndm(rndm,1)
      phi=rndm(1)*twpi
      if ((px .eq. 0.0) .and. (py .eq. 0.0)) goto 3
      if (abs(px) .lt. 1.e-10) goto 2
      phi=atan2(py,px)
      goto 3
c
 2    continue
      if (py .gt. 0.0) phi=pi/2.0
      if (py .le. 0.0) phi=3.0*pi/2.0
c
 3    continue
      sinp=sin(phi)
      cosp=cos(phi)
c
c --- set gheisha index for the current medium always to 1 ---
      ind=1
c
c --- transfer global material constants for current medium ---
c --- detailed data for compounds is obtained via routine compo ---
      nmix=maxmix
      call gfmate(imat,name,a,z,dens,radl,absl,nmix,imix,wmix)
      atno(ind+1)=a
      zno(ind+1)=z
      den(ind+1)=dens
      radlth(ind+1)=radl
      abslth(ind+1)=absl
      call compoi(imat,name,a,z,nmix,imix,wmix)
c
c --- setup parmat for physics steering ---
      parmat(ind+1,5)=0.0
      parmat(ind+1,8)=ipfis
      parmat(ind+1,9)=0.0
      parmat(ind+1,10)=0.0
      parmat(ind+1,5)=0.0 !qcor 0 for time being
c
c --- indicate light (<= pi) and heavy particles (historically) ---
c --- calim code ---
      j=2
      test=rmass(7)-0.001
      if (abs(amas) .lt. test) j=1
c
c *** division into various interaction channels denoted by "intval" ***
c the convention for "intval" is the following
c
c intval  = -1 reaction cross sections not yet tabulated/programmed
c      =  0 no interaction
c      =  1 eleastic scattering
c      =  2 inelastic scattering
c      =  3 nuclear fission with ineleastic scattering
c      =  4 neutron capture
c
c --- intact code ---
      kk=nmix
      alam1=0.0
      call grndm(rndm,1)
      rat=rndm(1)*alam
      atno2=a
      zno2 =z
c
      do 6 k=1,kk
      if (kk .le. 0) go to 6
c
      if (kk .eq. 1) go to 7
      ndum=10
      call gfmate(imix(k),name,atno2,zno2,dum,dum,dum,ndum,ndum,dum)
c
 7    continue
c
c --- try for elastic scattering ---
      intval=1
      code=13
      alam1=alam1+aiel(k)
      if (rat .lt. alam1) go to 8
c
c --- try for inelastic scattering ---
      intval=2
      code=20
      alam1=alam1+aiin(k)
      if (rat .lt. alam1) go to 8
c
c --- try for nuclear fission with inelastic scattering ---
      intval=3
      code=15
      alam1=alam1+aifi(k)
      if (rat .lt. alam1) go to 8
c
c --- try for neutron capture ---
      intval=4
      code=18
      alam1=alam1+aica(k)
      if (rat .lt. alam1) go to 8
c
 6    continue
c --- no reaction selected ==> elastic scattering ---
      intval=1
      code=13
c
c *** take action according to selected reaction channel ***
c --- following code is a translation of "calim" into geant jargon ---
c
 8    continue
c
c --- in case of no interaction or unknown cross sections ==> done ---
      if (intval .le. 0) go to 40
c
c --- in case of non-elastic scattering and no generation of sec. ---
c --- particles deposit total particle energy and return ---
      if ((intval .eq. 1) .or. (ihadr .ne. 2)) go to 9
      stop=2
      destep=destep+en
      go to 9999
c
 9    continue
      if (intval .ne. 4) go to 10
c
c --- neutron capture ---
      stop=1
      call captur(nopt)
      go to 40
c
 10   continue
      if (intval .ne. 3) go to 11
c --- nuclear fission ---
      stop=1
      tkin=fissio(ek)
      intval=0
      go to 40
c
 11   continue
c
c --- inelastic scattering at energies below 5 gev are treated ---
c --- by the nucrin/hadrin package in case ihadr=3 is selected ---
c --- reference : infn/be-88/3 by paolo pedroni (pavia italia) ---
c
      iok=0
      if ((ihadr .eq. 3) .and. (intval .eq. 2))
     $ call hadnuc
     $ (pvec,amas,atno2,zno2,ipart,iok,stop,mprod,nprod,iprod,pprod)
      if (iok .ne. 0) go to 9999
c
c --- elastic and inelastic scattering ---
      pv(1,200)=p*px
      pv(2,200)=p*py
      pv(3,200)=p*pz
      pv(4,200)=en
      pv(5,200)=amas
      pv(6,200)=nch
      pv(7,200)=tof
      pv(8,200)=kpart
      pv(9,200)=0.
      pv(10,200)=userw
c
c --- additional parameters to simulate fermi motion and evaporation ---
      enp(5)=ek
      enp(6)=en
      enp(7)=p
c
      if (intval .ne. 1) go to 12
c
c *** elastic scattering processes ***
c
c --- only nuclear interactions for heavy fragments ---
      if ((kpart .ge. 30) .and. (kpart .le. 32)) go to 35
c
c --- normal elastic scattering for light media ---
      if (atno2 .lt. 1.5) go to 35
c
c --- coherent elastic scattering for heavy media ---
      call coscat
      go to 40
c
c *** non-elastic scattering processes ***
 12   continue
c
c --- only nuclear interactions for heavy fragments ---
      if ((kpart .ge. 30) .and. (kpart .le. 32)) go to 35
c
c *** use sometimes nuclear reaction routine "nucrec" for low energy ***
c *** proton and neutron scattering ***
      call grndm(rndm,1)
      test1=rndm(1)
      test2=4.5*(ek-0.01)
      if ((kpart .eq. 14) .and. (test1 .gt. test2)) go to 85
      if ((kpart .eq. 16) .and. (test1 .gt. test2)) go to 86
c
c *** fermi motion and evaporation ***
      tkin=cinema(ek)
      enp(5)=ek+tkin
c --- check for lowerbound of ekin in cross-section tables ---
      if (enp(5) .le. teklow) enp(5)=teklow
      enp(6)=enp(5)+abs(amas)
      enp(7)=(enp(6)-amas)*(enp(6)+amas)
      enp(7)=sqrt(abs(enp(7)))
      tkin=fermi(enp(5))
      enp(5)=enp(5)+tkin
c --- check for lowerbound of ekin in cross-section tables ---
      if (enp(5) .le. teklow) enp(5)=teklow
      enp(6)=enp(5)+abs(amas)
      enp(7)=(enp(6)-amas)*(enp(6)+amas)
      enp(7)=sqrt(abs(enp(7)))
      tkin=exnu(enp(5))
      enp(5)=enp(5)-tkin
c --- check for lowerbound of ekin in cross-section tables ---
      if (enp(5) .le. teklow) enp(5)=teklow
      enp(6)=enp(5)+abs(amas)
      enp(7)=(enp(6)-amas)*(enp(6)+amas)
      enp(7)=sqrt(abs(enp(7)))
c
c *** in case of energy above cut-off let the particle cascade ***
      test=abs(charge)
      if ((test .gt. 1.0e-10) .and. (enp(5) .gt. cuthad)) go to 35
      if ((test .le. 1.0e-10) .and. (enp(5) .gt. cutneu)) go to 35
c
c --- second chance for anti-baryons due to possible annihilation ---
      if ((amas .ge. 0.0) .or. (kpart .le. 14)) go to 13
      anni=1.3*p
      if (anni .gt. 0.4) anni=0.4
      call grndm(rndm,1)
      test=rndm(1)
      if (test .gt. anni) go to 35
c
c *** particle with energy below cut-off ***
c --- ==> only nuclear evaporation and quasi-elastic scattering ---
 13   continue
c
      stop=3
c
c
      if ((kpart .ne. 14) .and. (kpart .ne. 16)) go to 14
      if (kpart .eq. 16) go to 86
c
c --- slow proton ---
 85   continue
      call nucrec(nopt,2)
c
      if (nopt .ne. 0) go to 50
c
      call coscat
      go to 40
c
c --- slow neutron ---
 86   continue
      nucflg=0
      call gnslwd(nucflg,intval,nfl,teklow,stop)
      if (nucflg .ne. 0) go to 50
      go to 40
c
c --- other slow particles ---
 14   continue
      ipa(1)=kpart
c --- decide for proton or neutron target ---
      ipa(2)=16
      call grndm(rndm,1)
      test1=rndm(1)
      test2=zno2/atno2
      if (test1 .lt. test2) ipa(2)=14
      avern=0.0
      nfl=1
      if (ipa(2) .eq. 16) nfl=2
      ippp=kpart
c      if (nprt(9)) print 2005
 2005 format(' *gheish* routine twob will be called')
      call twob(ippp,nfl,avern)
      goto 40
c
c --- initialisation of cascade quantities ---
 35   continue
c
c *** cascade generation ***
c --- calculate final state multiplicity and longitudinal and ---
c --- transverse momentum distributions ---
c
c --- fixed particle type to steer the cascade ---
      kkpart=kpart
c
c --- no cascade for leptons ---
      if (kkpart .le. 6) go to 9999
c
c *** what to do with "new particles" for gheisha ?????? ***
c --- return for the time being ---
      if (kkpart .ge. 35) go to 9999
c
c --- cascade of heavy fragments
      if ((kkpart .ge. 30) .and. (kkpart .le. 32)) go to 390
c
c --- initialize the ipa array ---
      call vzero(ipa(1),100)
c
c --- cascade of omega - and omega - bar ---
      if (kkpart .eq. 33) go to 330
      if (kkpart .eq. 34) go to 331
c
      nvepar=kkpart-17
      if (nvepar .le. 0) go to 15
      go to (318,319,320,321,322,323,324,325,326,327,328,329),nvepar
c
 15   continue
      nvepar=kkpart-6
      go to (307,308,309,310,311,312,313,314,315,316,317,318),nvepar
c
c --- pi+ cascade ---
 307  continue
      call caspip(j,intval,nfl)
      go to 40
c
c --- pi0 ==> no cascade ---
 308  continue
      go to 40
c
c --- pi- cascade ---
 309  continue
      call caspim(j,intval,nfl)
      go to 40
c
c --- k+ cascade ---
 310  continue
      call caskp(j,intval,nfl)
      go to 40
c
c --- k0 cascade ---
 311  continue
      call cask0(j,intval,nfl)
      go to 40
c
c --- k0 bar cascade ---
 312  continue
      call cask0b(j,intval,nfl)
      go to 40
c
c --- k- cascade ---
 313  continue
      call caskm(j,intval,nfl)
      go to 40
c
c --- proton cascade ---
 314  continue
      call casp(j,intval,nfl)
      go to 40
c
c --- proton bar cascade ---
 315  continue
c      if (nprt(9)) print 2013
 2013 format(' *gheish* routine caspb will be called')
      call caspb(j,intval,nfl)
      go to 40
c
c --- neutron cascade ---
 316  continue
      nucflg=0
      if (ek .gt. swtekn) call casn(j,intval,nfl)
      if (ek .le. swtekn) call gnslwd(nucflg,intval,nfl,teklow,stop)
      if (nucflg .ne. 0) go to 50
      go to 40
c
c --- neutron bar cascade ---
 317  continue
      call casnb(j,intval,nfl)
      go to 40
c
c --- lambda cascade ---
 318  continue
      call casl0(j,intval,nfl)
      go to 40
c
c --- lambda bar cascade ---
 319  continue
c      if (nprt(9)) print 2018
 2018 format(' *gheish* routine casal0 will be called')
      call casal0(j,intval,nfl)
      go to 40
c
c --- sigma + cascade ---
 320  continue
      call cassp(j,intval,nfl)
      go to 40
c
c --- sigma 0 ==> no cascade ---
 321  continue
      go to 40
c
c --- sigma - cascade ---
 322  continue
      call cassm(j,intval,nfl)
      go to 40
c
c --- sigma + bar cascade ---
 323  continue
      call casasp(j,intval,nfl)
      go to 40
c
c --- sigma 0 bar ==> no cascade ---
 324  continue
      go to 40
c
c --- sigma - bar cascade ---
 325  continue
      call casasm(j,intval,nfl)
      go to 40
c
c --- xi 0 cascade ---
 326  continue
c      if (nprt(9)) print 2023
 2023 format(' *gheish* routine casx0 will be called')
      call casx0(j,intval,nfl)
      go to 40
c
c --- xi - cascade ---
 327  continue
c      if (nprt(9)) print 2024
 2024 format(' *gheish* routine casxm will be called')
      call casxm(j,intval,nfl)
      go to 40
c
c --- xi 0 bar cascade ---
 328  continue
      call casax0(j,intval,nfl)
      go to 40
c
c --- xi - bar cascade ---
 329  continue
c      if (nprt(9)) print 2026
 2026 format(' *gheish* routine casaxm will be called')
      call casaxm(j,intval,nfl)
      go to 40
c
c --- omega - cascade ---
 330  continue
c      if (nprt(9)) print 2027
 2027 format(' *gheish* routine casom will be called')
      call casom(j,intval,nfl)
      go to 40
c
c --- omega - bar cascade ---
 331  continue
c      if (nprt(9)) print 2028
 2028 format(' *gheish* routine casaom will be called')
      call casaom(j,intval,nfl)
      go to 40
c
c --- heavy fragment cascade ---
 390  continue
      nucflg=0
      call casfrg(nucflg,intval,nfl)
      if (nucflg .ne. 0) go to 50
c
c *** check whether there are new particles generated ***
 40   continue
      if ((ntot .ne. 0) .or. (kkpart .ne. kpart)) go to 50
       nprod=1
       iprod(1)=ikpart(kpart)
       pprod(1,1)=px*p
       pprod(2,1)=py*p
       pprod(3,1)=pz*p
       edep=enold-en
       if(edep.gt.0.0.and.stop.eq.0) destep=destep+edep
      go to 9999
c
c *** current particle is not the same as in the beginning or/and ***
c *** one or more secondaries have been generated ***
 50   continue
c
      nvedum=kipart(ipart)
c
c --- initial particle type has been changed ==> put new type on ---
c --- the geant temporary stack ---
c
c
c
c --- put particle on the stack ---
       nprod=1
       iprod(1)=ikpart(kpart)
       pprod(1,1)=px*p
       pprod(2,1)=py*p
       pprod(3,1)=pz*p
       tprod(1)=tprod(1)+tof*0.5e-10
c
c
c *** check whether secondaries have been generated and copy them ***
c *** also on the geant stack ***
 60   continue
c
c --- all quantities are taken from the gheisha stack where the ---
c --- convention is the following ---
c
c eve(index+ 1)= x
c eve(index+ 2)= y
c eve(index+ 3)= z
c eve(index+ 4)= ncal
c eve(index+ 5)= ncell
c eve(index+ 6)= mass
c eve(index+ 7)= charge
c eve(index+ 8)= tof
c eve(index+ 9)= px
c eve(index+10)= py
c eve(index+11)= pz
c eve(index+12)= type
c
      if (ntot .le. 0) go to 9999
c
c --- one or more secondaries have been generated ---
      do 61 l=1,ntot
      index=(l-1)*12
      jnd=eve(index+12)
c
c --- make choice between k0 long / k0 short ---
      if ((jnd .ne. 11) .and. (jnd .ne. 12)) go to 63
      call grndm(rndm,1)
      jnd=11.5+rndm(1)
c
c --- forget about neutrinos ---
 63   continue
      if (jnd .eq. 2) go to 61
c
c --- swith to geant quantities ---
      ity=ikpart(jnd)
      plx=eve(index+9)
      ply=eve(index+10)
      plz=eve(index+11)
c
c
c --- add particle to the stack if stack not yet full ---
      fail=1301
      if (nprod.ge.mprod) go to 1313
      nprod=nprod+1
      iprod(nprod)=ity
      pprod(1,nprod)=plx
      pprod(2,nprod)=ply
      pprod(3,nprod)=plz
      tprod(nprod)=tprod(nprod)+eve(index+8)*0.5e-10
c
c
 61   continue
c
 9999 continue
      if(destep.gt.0.0) then !add pseudo particle
       nprod=nprod+1
       iprod(nprod)=200
       pprod(1,nprod)=0.0
       pprod(2,nprod)=0.0
       pprod(3,nprod)=destep
      endif !destep.gt.0.0
      return
c
1313  continue
      write(6,*) 'fail=',fail,' in gheish'
      if(fail.eq.1301) then
       write(6,*) 'produced particle array overflowed'
       write(6,*) 'results will be truncated'
      endif !fail.eq.1301
      return
c
      entry gheiset(todo,arg)
      value=arg
      if(todo.eq.'ihadr') then
       ihadr=ivalue
      else if (todo.eq.'ipfis') then
       ipfis=ivalue
      endif !todo.eq.'ihadr'
      end
