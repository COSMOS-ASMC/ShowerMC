c     subroutine nucrinC(itttt,elab,cx,cy,cz,anuc,znuc,rhoo)
      subroutine nucrinC(itttt,elab,cx,cy,cz,anuc,znuc)
c$$$$$$$$$$$$$$$$$$$
c          gauss---->aagausC
c$$$$$$$$$$$$$$$$

      common/hfragC/toeff,atar,sqatar,tpnuc,tanuc
      common /tpnkcC/tpnk
      common /fermtC/ elabke,a,ifert,tpnvk,elate
c
c--------------------------------------------------
c
c***
c***
c***  sampling of a hadron-nucleus collision event
c***
c   itttt - type of incoming hadron
c
c--------------------------------------------------
c***  possible: n, an,pi,k,ak,y,ay by the index itttt=it=1,2,8-25
c   plab - momentum of incoming hadron
c
c***  range up to 5. gev/c (from about 0. ... 0.1 gev/c)
c   elab - total energy of incoming hadron in gev
c   cx,cy,cz - direction cosines
c***  anuc,znuc = numbers of nucleons and protons
c*** rhoo = material density (g/cm**3)
c
c*** final state particle characteristics table in /finuc/:
c*** irn   - number of final state particles
c*** itrn(i) - final state particle type index
c*** cxrn,cyrn,czrn (i) - direction cosines of f.s.p. (lab syst.)
c*** elr,plr (i) - lab.energy and momentum of f.s.p. (gev, gev/c)
c*** tv    - excitation energy (gev)
c
c--------------------------------------------------
      dimension thresr(30)
      dimension ixpi(30)

      common /regimC/ xregimC(4),signa
      common /percoC/ iperco
      common /encoteC/ etest,tnkte
     *      /nucparC/cddt
      common /probaC/ jnuc,is1,nq,iq,bit,w(32)

c
c--------------------------------------------------
c*** finuc: nucrinC final state particle list with kinem. variables
c*** finlspC hadrin final state particle list with kinem. variables
c*** (number of particles,particle type index,direction cosines,energy
c*** absolute momentum and in finuc in addit. excitation energy
c--------------------------------------------------
     *   /finucC/irn,itrn(60),cxrn(60),cyrn(60),czrn(60),
     *elr(60),plr(60),tv
     *      /finlspC/ir,itr(20),cxr(20),cyr(20),czr(20),
     *el(20),pl(20)
c      integer*2ich,ibar,k1,k2,nzk,nrk
      integer*2 ich,ibar,k1,k2
c
c--------------------------------------------------
c***  particle characteristics: masses, decay width, life time,elect.and
c*** baryonic charge,decay channel indicees
c--------------------------------------------------
      common/abltisC/am(110),ga(110),tau(110),ich(110)
     *,ibar(110),k1(110),k2(110)
     *   /thrC/ ethr,pthr
      dimension ins(30)

      common /eventC/ ievent
      integer*2jnuc,is1,nq,iq
      dimension amhh(15),iamc(15)
      dimension itprf(110)
c############################
      real*8  v, vz, elran, vt
c########################

      real amsave(15)

      data cdddd/0./
      data ceet/0.5/

      

      data iamc/1,2,8,9,12,15,16,17,18,19,20,21,22,24,25/
      data ins/13,13,5*32,14,14,3*32,10,10
     *,12,11,7*32,15,15,5*32/


      data itprf/-1,-1,5*1,-1,-1,1,1,1,-1,-1,-1,-1,6*1,-1,-1,-1,85*1/
      data
     * thresr/1.9,0.,5*9.,1.9,0.,3*9.,1.08,1.08,1.44,1.08,6*9.,1.08,
     *1.44,1.08,5*9./

      data ixpi/1,-1,5*2,1,-1,3*2,1,0,1,0,4,4,2,3*4,1,1,1,5*2/
      data ieveno/0/
      save ieveno, amsave



      itnucr=itttt
      if (itprf( itttt).lt.0) go to 99999
      write (6,99998)itttt
      stop
99998 format (3(5h ****/),
     *45h false use of the particle type index, value ,i4,3(5h ****/))
99999 continue
      sico=1.
      elph=0.
      plabco=1.
c  &&&&&&&&&&&&&6666666 since am is changed inside
c           and might not be restored, we every time
c           restore it from the second call  KK v.6.16
c
      if(ieveno .eq. 0 ) then
c              save mass
         do i = 1, 15
            amsave(i)=  am(iamc(i))
         enddo
      else
c               restore  
         do i=1,15
            am(iamc(i)) = amsave(i)
         enddo
      endif
c    &&&&&&&&&&&&&&
c
c
c--------------------------------------------------
c*** sico,plabco - effective cross section- and effect. labmomentum
c***      correction factors for y-a- and ay-a-collisions
c*** itnucr - hyperonC nucleus coll. particle type for use in nizl,
c***      it - used in hadrin
c--------------------------------------------------
      ibat=0
      ieveno=ieveno+1
      ievent=ieveno
c
c--------------------------------------------------
c*** inxpi:
c*** ordering indicees for stopping particles (0), other mesons (1),
c*** antinucleons (-1), hyperoCns (4),leptons,ko,ako(2),other part. (1)
c
c--------------------------------------------------
c*** inclusion of hyperoCn-nucleus-collisions
c*** preliminary for sigma+ no strangenes-conservation
c
c--------------------------------------------------
         it=itttt
         nxpi=0
         inxpi=ixpi(it)
         if (inxpi.eq.4) call hyperoC(it,itnucr,sico,plabco)
c
c--------------------------------------------------
c*** cut off energy constants:
c--------------------------------------------------
      tnkte=0.001
      ethr=0.001
c
c--------------------------------------------------
c*** cddt,ceet parameters for gaussian width'
c--------------------------------------------------
      n=itttt
      if (itttt.eq.18) inxpi=-1
      do 1513 iam=1,15
1513  amhh(iam)=am(iamc(iam))
c
c--------------------------------------------------
c*** etest = energy conservation test variable (should be 0 at return)
c*** tlab = kinetic energy
c--------------------------------------------------
      jj=jnuc
      cddt=0.5
      ihaca=0
      in=ins(n)
c
c--------------------------------------------------
c*** jj=1:
c*** option 1: go to nizl, cross section calculation, replace "ccc"
c*** jj=2:
c*** option 2: event sampling
c--------------------------------------------------
      etest=elab
      tlab=max(elab-am(it), 0.)
c
c---------------------------------------------------
c*** decide, whether p or neu in h-n-coll. inside a is hitten
c*** itta will be the target nucleon index:
c------------------------------------
      itta=8
      call rndc(v)
      if (v.le.znuc/anuc) itta=1
      plab=sqrt(tlab*(elab+am(it)))
      if (abs(plab-5.).lt.4.99999) go to 99997
      if(it  .ne. 2 .and. it  .ne. 9) then
         write (6,99996) plab
         write(6,  *) ' it, elab, mass it',
     *    itttt, elab, am(it), it
c///////////////
         plab =sin(0.)
         plab = (1.+plab)/plab
         write(*,*) ' plab=', plab
c////////////
         stop
      endif
99996 format (3(5h ****/),64
     *h projectile nucleus momentum outside of the allowed region, plab,
     *1e15.5/,3(5h ****/))
99997 continue
      if (jj.eq.1) go to 1000
c
c--------------------------------------------------
c*** if low kinetic energy (tlab<ethr): only excitation
c*** stopping particles
c*** and annihilation , inxpi=ixpi(it) le 0
c--------------------------------------------------
         if (inxpi.le.0) go to 402
c
c--------------------------------------------------
c*** ethrr=threshold value, for tlab<ethrr is no h-n collision possible,
c*** only excitation + cascade and/or annihilation corresp to the hadrn
c
c--------------------------------------------------
c*** for ap, pi-, k- go to 402
c--------------------------------------------------
         if (tlab.gt.ethr) go to 401
         tv=tlab
         irn=0
         return
 402     continue
         ethrr=ethr
         if (inxpi.lt.0) ethrr=0.08
c
c--------------------------------------------------
c*** up to ethrr=.08 gev no an, ap in final state
c--------------------------------------------------
         aio=10.
         nxpi=1
          xpi=1.-aio*(tlab-ethrr)
c
c--------------------------------------------------
c*** nxpi=1 pi-,k- absorption possible up to xpi=0
c--------------------------------------------------
         call rndc(vz)
         if (vz.gt.xpi) nxpi=0
c
c--------------------------------------------------
c*** nxpi=0 means hadrin-call is possible, else impossible
c--------------------------------------------------
401      continue
      if (nxpi.le.0.or.ibar(it).ne.0) go to 4401
c
c--------------------------------------------------
c*** give the total lab energy of stopping pi-, k- into cascade and exit
c*** energy calculation for conservation tests in fermi momentum version
c--------------------------------------------------
      elpp=elab
      tpk=0.
      tnk=0.
      tv =0.
      etest=elpp
      tpnuc=znuc
      tanuc=anuc
      toeff=elpp
      atar=anuc
      sqatar=sqrt(atar)
      go to 405
4401  continue
      cdddd=cddt
      inuc=anuc-znuc
      if (itta.le.1) inuc=znuc
c
c--------------------------------------------------
c***  elabke = invariant kinetic energy of the h-a-system+proj.h.mass
c*** a = total energy of the primary particle in labsystem
c*** iforbi,ifert loop counting indicees for choice of fermion momentum
c*** akmas - nucleus mass
c*** umojan - h-a-c.m.s.-energy
c*** elabke - kinetic h-a-energy in c.m.s
c--------------------------------------------------
      akmas=0.938*znuc+.940*(anuc-znuc)
      umojan=sqrt(am(it)**2+akmas**2+2.*elab*akmas)
c
c--------------------------------------------------
c***  start with the event simulation
c--------------------------------------------------
      elabke=umojan-akmas
      a=elab
      iforbi=0
      ifert=0
1651  continue
      irn=0
      tanuc=anuc
      tpnuc=znuc
      tnnuc=anuc-znuc
      ddt=cddt*tlab
      ammit=am(it)
      ammta=am(itta)
      ait2=ammit**2
c
c--------------------------------------------------
c*** effective collision energy calculation for use in cascade and
c*** and excitation energy calculation for annihilation
c--------------------------------------------------
      ata2=ammta**2
c
c--------------------------------------------------
c*** i653 loop counter integer for bad energy correction in annihilation
c--------------------------------------------------
      i653=0
      ecmo=ait2+ata2+2.*ammta*elab
      emd=0.
653   continue
c
c--------------------------------------------------
c*** toeff = effective kinetic collision energy
c--------------------------------------------------
      i653=i653+1
         toeff=tlab
c
c--------------------------------------------------
c*** toeff only for annihilation ne tlab possible, else just selection
c*** of energy fractions
c--------------------------------------------------
         if (inxpi.gt.0) goto 6513
      ibl=ibar(itta)*ibar(it)
      ibl=iabs(ibl)
      icl=ich(it)
      emd=am(it)*(((ibar(itta)-ibar(it))*ibl
     * + (1-ibl)*iabs(icl)))
c
c--------------------------------------------------
c*** emd=2*am(it) for annihilation, =am(it) for mesons, =0else
c--------------------------------------------------

      thres=thresr(it)
      ecms=sqrt(ecmo)
      if (ecms.lt.thres) emd=0.
         if (nxpi.lt.1) emd=0.
c
c--------------------------------------------------
c*** emd=0 for mesons and baryons,if no stopping of particles if h-n-cms
c*** -energy lt tabulated threshold
c--------------------------------------------------
      if (ibar(it).ge.0) emd=0.
6512   continue
c
c--------------------------------------------------
c*** eci = total effective energy in h-n-cms-system
c--------------------------------------------------
         eci=((sqrt(ecmo)+emd)**2-ait2-ata2)/(2.*ammta)
         toeff=eci-ammit
      if (toeff.lt.tlab) toeff=tlab
6513     continue
      t1=anuc
      t2=sqrt(t1)
      atar=t1
c
c--------------------------------------------------
c*** tpk0,tnk0 - proton/neutron average cascade energy
c*** eexo average evaporation energy, eex average nuclear fragment energ
c*** +average cascade energy
c
c--------------------------------------------------
      tpk0=ekekaC( 2,toeff,t1,t2)
      tnk0=ekekaC( 3,toeff,t1,t2)
      eexo=ekekaC( 1,toeff,t1,t2)
      eex=eexiC(1,tlab,t1)
54    continue
      dex=ceet*eexo
c
c--------------------------------------------------
c*** excitation energy calculation
c*** gaussian distribution around eex
c--------------------------------------------------
      if (eexo.gt.eex) dex=ceet*eexo
      call aagausC(tvv,eex,dex)
      tv=tvv-tpk0-tnk0
      tvh=tvv*0.35-tpk0-tnk0
59    if (tvh.gt.0.) go to 58
      call aagausC(tv,eexo,dex)
c
c--------------------------------------------------
c*** gaussian distribution around eexo
c--------------------------------------------------
      if (tv) 59,59,599
58    continue
      if(tvv.ge.tlab) go to 54
599   continue
      t1=tpk0+tnk0
      t2=toeff-tv
      itpn=0
      if(t2.lt.t1) go to 54
50    continue
c
c--------------------------------------------------
c*** if sampled excitation energy tv is not inside the kinematic region
c*** start again tv-random choice
c--------------------------------------------------
      itpn=itpn+1
c
c--------------------------------------------------
c*** gaussian random choice of the cascade energy tpnk
c--------------------------------------------------
      ddt=0.5*t1
      call aagausC(tpnk,t1,ddt)
      if(itpn.gt.1) tpnk=t1+abs(tpnk-t1)*(stpn/atpn)
c
c--------------------------------------------------
c*** equal weight of each branch of the gaussian in the allowed region
c--------------------------------------------------
      stpn=tpnk-t1+1.e-18
      atpn=abs(stpn)
      if (tpnk.gt.t2) go to 54
      if(tpnk.lt.0.) go to 50
c
c--------------------------------------------------
c*** now tpnk within allowed region found
c*** calculation of proton and neutron cascade energy (tpk,tnk)
c--------------------------------------------------
      tpp=tlab-t1
      tpk=tpnk*tpk0/t1
      tnk=tpnk*tnk0/t1
c
c--------------------------------------------------
c
c   event simulation
c   hadron-nucleon interaction
c
c--------------------------------------------------
c*** elpp=total projectile energy,plpp=absolut momentum of the projectil
c*** in the labsystem
c
c--------------------------------------------------
      t3=am(it)
      elpp=elab-tpk-tnk-tv
         if (inxpi.ne.0.or.nxpi.eq.0) go to 403
c
c--------------------------------------------------
c*** if no primary pi-,k- or if hadrin for pi-,k- will be called
c*** (no stopping particle case) go to 403
c--------------------------------------------------
      if (ibar(itttt).lt.0) go to 403
c
c--------------------------------------------------
c*** if itttt is a stopped negative meson,(but also for annihilation)
c*** distribute the available total energy of the particle to the cascad
c*** and excitation
c--------------------------------------------------
405      continue
      elpo=elpp
      if (ibar(it).gt.0) elpo=elpp-am(itttt)
      tpk=tpk+.3*elpo
      tnk=tnk+.3*elpo
      tv= tv +.4*elpo
         elpp=0.000001
         go to 407
c
c--------------------------------------------------
c*** go to the cascade simulation
c--------------------------------------------------
403      continue
c
c--------------------------------------------------
c*** case of energy correction if loop in annihilation
c--------------------------------------------------
      if (i653.lt.10) go to 1653
c
c--------------------------------------------------
c*** for less then 10 iterations (i653), start again with effective coll
c*** energy calculation, else correction of tpk,tnk,tv if elpp le 0.
c--------------------------------------------------
      if (elpp.gt.0.) go to 1653
      a65=1.+elpp/(tpk+tnk+tv+1.e-18       )
      tpk=tpk*a65
      tnk=tnk*a65
      tv =tv *a65
      elpp=0.
1653  continue
      if (elpp.lt.0.) go to 653
      if(elpp.gt.t3)go to 510
c
c--------------------------------------------------
c*** if elpp greaterthen the mass of it, go to momentum limitation for h
c*** else if no effective  energy was calculated do the same
c*** else if it is no antinucleon and elpp lt mass of it, start again wi
c***   effective collision energy calculation
c*** else new energy definition for it and new mass definition for baryo
c*** and antibaryons
c--------------------------------------------------
      if (emd .le.0.) go to 1510
      if (ibar(it).ge.0.and.elpp.lt.t3) go to 653
      elph=elpp+am(itta)-1.e-6
2511  continue
      t3=elph/2.
c
c--------------------------------------------------
c*** case of annihilation,pseudo mass definition
c--------------------------------------------------
      elpp=t3
      elpp=elph
      if (t3.lt.0.14.and.inxpi.gt.0.and.nxpi.gt.0) go to 405
      if (ibar(itttt).ge.0) go to 5010
c
c--------------------------------------------------
c*** if in annihilation for a stopping antinucleon not more then 2 pion
c*** masses available, give this energy to cascade ane eccitation
c--------------------------------------------------
         etdd=0.
         etkor=tv+tpk+tnk+ammit
         if (nxpi) 15010,15010,15011
15010 continue
         et=elab
c
c--------------------------------------------------
c*** the particle it (=an) is not stopped in the nucleus
c--------------------------------------------------
         etest=et
         etrest=et-etkor
         if (etrest.lt.0.) etdd=etrest
         tkor=1.+etdd/(tv+tnk+tpk+1.e-10)
         if (tkor.lt.1.e-6) go to 15011
         tv=tv*tkor
         tpk=tpk*tkor
         tnk=tnk*tkor
         elpp=et-tv-tpk-tnk
         go to 510
15011 continue
         etdd=0.
c
c--------------------------------------------------
c*** the particle it (=an) is     stopped in the nucleus
c--------------------------------------------------
         et=elab    +ammta
         etest=et
         etrest=et-etkor
         if (etrest.lt.0.) etdd=etrest
         tkor=1.+etdd/(tv+tnk+tpk+1.e-10)
         if (tkor.lt.1.e-6) go to 653
         tv=tv*tkor
         tpk=tpk*tkor
         tnk=tnk*tkor
         elpp=et-tv-tpk-tnk
         t3=elph/2.
5010  continue
c
c--------------------------------------------------
c*** new mass definition of baryons and antibaryons
c*** (only for remaining effective cm-energies,lower than the total
c*** 1.88gev-threshold cm-energy in annihilation)
c--------------------------------------------------
      do 1511 iam=1,15
      jam=iamc(iam)
1511  am(jam)=t3

1510  continue
      elpp=t3+1.e-6
510   continue
      plpp=sqrt(elpp**2-t3**2)
c
c--------------------------------------------------
c*** momentum limitation for hadrin in fermi-mom.-version
c--------------------------------------------------
      plabou=15.
      if (plpp.lt.plabou) go to 5101
      apl=plpp-plabou
      elpp=sqrt(plabou**2+am(it)**2)
      tpk=tpk+apl*.3
      tnk=tnk+apl*.3
      tv =tv +apl*.4
5101  continue
      iforbi=iforbi+1
      if (iforbi.gt.10) go to 405
      tpnvk=tpk+tnk+tv
      pplab=plab
c
c--------------------------------------------------
c***   peripheral collision
c--------------------------------------------------
      eelab=elab
c
c--------------------------------------------------
c*** i.e.,h-collision with n in the most outside nucleus shell
c*** with the thickness of one nucleon radius
c*** ration of the shell volume to nucleus volume is elanu3
c*** the h-n-coll. will be sampled firstly with the total primary energy
c*** elab, the remaining energy is used in the above sampeld ratio for
c*** cascade and excitation
c*** switch variable is iperco=1
c*** random choice by elran
c--------------------------------------------------
      iperco= 1
      call rndc(elran)
      elanu3=anuc**(-0.7)
c
c--------------------------------------------------
c*** if the pripheral case should be switched off, set elanu3=0.
c--------------------------------------------------
      if (elran.gt.elanu3) go to 882
      if (abs(elpp-t3).gt.0.0001) go to 881
882   continue
      pplab=plpp
c
c--------------------------------------------------
c*** no peripheral collission
c--------------------------------------------------
      eelab=elpp
c
c--------------------------------------------------
c*** i.e., th e remaining energy(elab-cascade-excitation) is used for
c*** h-n-collision
c*** switch variable is iperco =-1
c--------------------------------------------------
      iperco=-1
881   continue
      ihaca=ihaca+1
      call ferhadC(anuc,znuc,it,pplab,eelab,cx,cy,cz,itta)
c
c--------------------------------------------------
c***
c***
c*** hadron nucleon collision simulation
c*** eelab=labenergy of h,pplab=absolut momentum of h, cx,cy,cz directio
c*** cosines of h,in the lab system, itta target nucleon index,
c*** anuc,znuc nucleon and proton numbers
c***
c*** regarding fermi momentum of the nucleon
c***
c--------------------------------------------------
      if (ihaca.gt.30) go to 88
c
c--------------------------------------------------
c***
c***
c*** ihaca times called h-n-coll., if in final state no parti cle,call a
c*** again,
c*** if in annihilation of a stopping nucleon,   in the final
c*** state is a baryon, aall again h-n-coll.
c*** (not more than 30 times, else take this state)
c***
c***
c--------------------------------------------------
      ihaca=ihaca+1
      if (ir.eq.0) go to 510
      if (ihaca.gt.1000) go to 88
      if (inxpi.lt.0.and.ir.eq.2.and.ibar(itr(1)).ne.0
     *.and.nxpi.gt.0) go to 510
c
c--------------------------------------------------
c*** if no anihilation into mesons for an+n
c--------------------------------------------------
      if (ir.lt.2.and.inxpi.lt.0) go to 510
88    continue

      do 1512 iam=1,15
c
c--------------------------------------------------
c*** redefinition of masses, if thes were changed for annihilation
c*** beyond 1.88gev-threshold
c--------------------------------------------------
1512  am(iamc(iam))=amhh(iam)
c
c--------------------------------------------------
c*** correction of cascade and excitation energies for fermi-momentum
c*** if energy conservation allows no h-n-collision, go to the start
c*** point of event simul.
c*** else energy correction by variable ekikor
c--------------------------------------------------
c--------------------------------------------------
      if (ifert.lt.0) go to 1651
      ekikor=elpp-elate
      if (-ekikor.gt.tpnvk) go to 1651
      aa65=1.+ekikor/tpnvk
      tpk=tpk*aa65
      tnk=tnk*aa65
      tv =tv *aa65
11651 continue
c
c--------------------------------------------------
c
c*** redefine the energies of sampled particles in case of changed masse
c
c--------------------------------------------------
      meson=0
      do 2510 irz=1,ir
      itrz=itr(irz)
      ibarot=ibar(itrz)
      if (ibarot.eq.0) meson=meson+1
      if (iabs(ibarot).lt.1.or.inxpi.ge.0) go to 2510
      el(irz)=sqrt(pl(irz)**2+am(itrz)**2)
2510  continue
c
c--------------------------------------------------
c
c** give random angle like for casc. part. to single sec-ry from hadrin
c--------------------------------------------------
      if (ir.ge.2) go to 1221
      tki=el(1)-am(itr(1))
     *+1.e-6
      if(tki.le.0.) goto 1221
      ade=0.12*(1.+0.003*atar)/tki
      dex=exp(-1.57*1.57/ade)
      an1=(1.-dex)*ade/2.
      an2=dex
     **1.57
      an=an1+an2
      an1=an1/an
      call rndc(v)
      if(v.gt.an1) goto 1222
 1223 continue
         call rndc(v)
c
c--------------------------------------------------
c*** teta angle determination (single part.case), = de
c--------------------------------------------------
      de=sqrt(-ade*log(1.-v*(1.-dex)))
      if(de.gt.1.57) goto 1223
      goto 1224
 1222 continue
      call rndc(v)
      de=-v
      de=atan2(sqrt(1.-de**2),de)
 1224 continue
      call cosiC(sfe,cfe)
c
c--------------------------------------------------
c*** cos phi, sin phi determination(s.p.case)
c--------------------------------------------------
      sid=sin(de)
      cod=cos(de)
      call transC(cxr(1),cyr(1),czr(1),cod,sid,sfe,cfe,ccxr,ccyr,cczr)
c
c--------------------------------------------------
c*** turning of angle s into the labsystem (s.p.case)
c--------------------------------------------------
      cxr(1)=ccxr
      cyr(1)=ccyr
      czr(1)=cczr
 1221 continue
      etest=elab
c
c--------------------------------------------------
c*** handling of energy conservation test (etest) and
c*** nucleon number correction because of itta choice
c--------------------------------------------------
      if (ir.lt.2) go to 12
      tanuc=tanuc-1.
      i1=itta/7
      a1=i1
      tpnuc=tpnuc-1.+a1
      tnnuc=tnnuc-a1
      etest=etest+am(itta)
   12 continue
      irn=0
      irtvi=0
      irtvj=0
      ithinv=0
      ithinr=0
      itan=0
      do 10 i=1,ir
c
c--------------------------------------------------
c*** storage of h-n-coll. products into nucleus final state c./finucC/
c*** cut off for kinetic energies less then ethr, push them into ex.en.
c*** tv
c--------------------------------------------------
      t1=el(i)
      i1=itr(i)
      if (t1.gt.ethr+am(i1)) go to 11
      tv=tv+t1
      itvt=0
      if (irtvi.gt.0) goto 1201
      if (i   .gt.ir) go to 1205
      istv=0
      do 1203 irt=i,ir
c
c--------------------------------------------------
c*** tv,etest-calculation with respect to  particles (stoppingp.,
c*** no stopping p., produce4 or only scattered particles)
c--------------------------------------------------
      istv=istv+1
      itrirt=itr(irt)
      if (ithinr .eq. itrirt .or. ithinv .gt. 0)goto 1203
      if (ibar(itrirt).ne.ibar(it)  .or .abs(am(itrirt)-am(it  )).gt.0.0
     *1) go to 1203
      itvt=1
      irtvi=irtvi+1
      ithinv=itrirt
      go to 1207
1203  continue
1207  continue
      if (irtvi.gt.0) tv=tv-am(it)
      if (irtvi.gt.0) etest=etest-am(it)
1205  continue
1201  continue
      if (irtvj.gt.0) goto 1202
      if (i   .gt.ir) go to 1206
      istv=0
      do 1204 irt=i,ir
      istv=istv+1
      itrirt=itr(irt)
      if (ithinv .eq. itrirt .or. ithinr.gt.0) goto 1204
      if (ibar(itrirt).ne.ibar(itta) .or. abs(am(itrirt)-am(itta))
     *  .gt. 0.01) go to 1204
      irtvj=irtvj+1
      ithinr=itrirt
      go to 1208
1204  continue
1208  continue
      if (itvt.gt.0) go to 1206
      if (irtvj.gt.0) tv=tv-am(itta)
1206  continue
1202  continue
      if(abs(abs(i1-4.5)-3.5).gt.0.1)go to 10
      i12=i1/7
c
c--------------------------------------------------
c*** nucleon number correction for nucleons, given into tv
c*** then go to the  next particle
c--------------------------------------------------
      a1=i12
      itan=itan+1
      if(itan-1.gt.0)go to 10
      tpnuc=tpnuc+1.-a1
      tnnuc=tnnuc+a1
      tanuc=tanuc+1.
       go to 10
   11 continue
      irn=irn+1
c
c--------------------------------------------------
c*** store the particle characteristics into c./finucC/, then take next
c*** particle
c--------------------------------------------------
      itrn(irn)=i1
      cxrn(irn)=cxr(i )
      cyrn(irn)=cyr(i )
      czrn(irn)=czr(i )
      elr(irn)=el(i )
      plr(irn)=pl(i )
      ibat=ibat+1
      etest=etest-el(i)
   10 continue
      go to 406
407   continue
c
c--------------------------------------------------
c***  zero secondary part. case (h-n-coll.), stopping mesons
c--------------------------------------------------
      irn=0
406   continue
      iu=irn
c
c--------------------------------------------------
c
c   simulation of cascade neutrons
c
c--------------------------------------------------
      ibat=ibat-ibar(it)-ibar(itta)
      if (ibat+1) 4406,4407,4408
4406  if (ibar(it).eq.0) etest=etest-0.939
4407  continue
      if (ibar(it).eq.0) etest=etest-am(itta)
4408  continue
c
c--------------------------------------------------
c*** nucleon number and cascade energy(versus cut off)-test
c*** give the energy to  the excitation , if test demands this
c--------------------------------------------------
      if ((tpk.le.tnkte).or.(tpnuc.lt.0.5)) tv=tv+tpk
      if ((tnk.le.tnkte).or.(tnnuc.lt.0.5)) tv=tv+tnk
   20 continue
      t2=am(8)
c
c--------------------------------------------------
c*** neutron number and neutron cascade energy test
c--------------------------------------------------
      if(irn.gt.58)go to 30
      if (tnnuc.lt.0.5) go to 30
      if (tnk.le.tnkte) go to 30
c
c--------------------------------------------------
c*** sampling of kinetic energy tn and teta angle dn of the emitted neu
c--------------------------------------------------
      call rakekaC(2,toeff,atar,sqatar,tn,pn,dn)
      tnk=tnk-tn
c
c--------------------------------------------------
c*** test: if the remaining energy and corresp.nucleon number allows fur
c*** ther corresp. nucleon emission, else give the  the energy or takesrnc
c*** the missing form the excitation
c*** if for last chosen nuc leon tn>remaining rest, take in 90% the ener
c*** difference from excit.energy, in 10% or if tv would be less then 0.
c***  add tn to tv
c--------------------------------------------------
      if(irn.gt.58)go to 211
      if (tnnuc.lt.1.5) go to 211
      if (tnk.gt.tnkte) go to 21
  211 continue
      call rndc(vt)
      thk=tv+tnk
      if (vt.gt.0.1.or.thk.lt.0.) tv=tv+tnk+tn
      if (vt.gt.0.1) go to 30
      if (thk.lt.0.) go to 30
      tnk=0.
      tv=thk
   21 continue
c
c--------------------------------------------------
c*** corresp. nucleon number counting
c--------------------------------------------------
      tnnuc=tnnuc-1.
      tanuc=tanuc-1.
      if(tn.ge.ethr)go to 22
      tv=tv+tn
      tnnuc=tnnuc+1.
      tanuc=tanuc+1.
      go to 20
   22 continue
c
c--------------------------------------------------
c*** take the emitted corresp. nucleon into the final particle table
c*** c./finucC/ after choice of phi and turning into the lab system
c--------------------------------------------------
      irn=irn+1
      itrn(irn)=8
      elr(irn)=tn+t2
      etest=etest-tn
      pn=sqrt((tn+t2)**2-t2**2)
      plr(irn)=pn
      call cosiC(sfe,cfe)
      sid=sin(dn)
      cod=cos(dn)
      call transC(cx,cy,cz,cod,sid,sfe,cfe,
     *cxrn(irn),cyrn(irn),czrn(irn))
      go to 20
   30 continue
      if(irn.gt.59)go to 40
      if (tpnuc.lt.0.5) go to 40
      if (tpk.le.tnkte) go to 40
c
c--------------------------------------------------
c
c   simulation of cascade protons
c
c--------------------------------------------------
c*** proton number  and proton cascade energy test
c
c
c--------------------------------------------------
c*** sampling of kinetic energy tn and teta angle dn of the emitted p
c--------------------------------------------------
c--------------------------------------------------
      t2=am(1)
      call rakekaC(1,toeff,atar,sqatar,tn,pn,dn)
c
c--------------------------------------------------
c*** test: if the remaining energy and corresp.nucleon number allows fur
c*** ther corresp. nucleon emission, else give the  the energy or take
c*** the missing form the excitation
c*** if for last chosen nucleon tn>remaining rest, take in 90% the energ
c*** difference from excit.energy, in 10% or if tv would be less then 0.
c***  add tn to tv
c--------------------------------------------------
      tpk=tpk-tn
      if(irn.gt.59)go to 311
      if (tpnuc.lt.1.5) go to 311
      if (tpk.gt.tnkte) go to 31
  311 continue
      call rndc(vt)
      thk=tv+tpk
         if (vt.gt.0.1.or.thk.lt.0.) tv=tv+tpk+tn
      if (vt.gt.0.1) go to 40
      if (thk.lt.0.) go to 40
         tpk=0.
      tv=thk
c
c--------------------------------------------------
c*** corresp. nucleon number counting
c--------------------------------------------------
   31 continue
      tanuc=tanuc-1.
      tpnuc=tpnuc-1.
      if(tn.ge.ethr)go to 32
      tv=tv+tn
      tanuc=tanuc+1.
      tpnuc=tpnuc+1.
      go to 30
   32 continue
c
c--------------------------------------------------
c*** take the emitted corresp. nucleon into the final particle table
c*** c./finucC/ after chi ce of phi and turning into the lab system
c--------------------------------------------------
      irn=irn+1
      itrn(irn)=1
      elr(irn)=tn+t2
      etest=etest-tn
      pn=sqrt((tn+t2)**2-t2**2)
      plr(irn)=pn
      call cosiC(sfe,cfe)
      sid=sin(dn)
      cod=cos(dn)
      call transC(cx,cy,cz,cod,sid,sfe,cfe,
     *cxrn(irn),cyrn(irn),czrn(irn))
      go to 30
c
c--------------------------------------------------
c*** end of cascade nucleon emission
c*** handling of the cases: no or one particle in h-a-final state
c*** final calculation of tv and etest
c--------------------------------------------------
   40 htv=0.
      esuev=0.
      ibari=ibar(itttt)
      if (ibari.gt.0) htv=am(itttt)
      if (ibari.lt.0) htv=-am(itttt)
c/////////////KK
c      if (irn.lt.1) go to 4014
c      htv=0.
c      ibai=0
c////////////// ibai may be used without
c               value assinged. so order changed.
c        
      ibai=0
      if (irn.lt.1) go to 4014
      htv=0.
c//////////////
      ianti=0
      do 4011 ies=1,irn
      iti=itrn(ies)
      ibi=ibar(iti)
      if (ibi.gt.0) ibai=ibai+1
      if (ibi.lt.0) ianti=ianti+1
      amim=am(iti)*(1-iabs(ibi))
      esuev=esuev+elr(ies)-am(iti)+amim
4011  continue
      htv=esuev
      if (ibari.lt.0.and.ianti.le.0) htv=htv-am(itta)
      inuoc=0
      inuo=1+ibari
      if (ianti.gt.0) inuo=2
      do 4013 ies=1,irn
      iti=itrn(ies)
      ibi=ibar(iti)
      if (ibi.eq.0) go to 4013
      inuoc=inuoc+1
      if (inuoc.eq.1) htv=htv-am(itta)
      if (inuoc.gt.inuo) go to 4014
      if (ianti.gt.0.and.ibai.le.0) htv=htv+am(itta)
      htv=htv+am(iti)
4013  continue
4014  if (ibari.gt.0.and.ibai.eq.1) htv=htv+am(itta)
      tv=elab-htv
      etest=elab-tv-htv
      if (tv.lt.0.) go to 99999
      return
 1000 continue
c
c--------------------------------------------------
c   sampling random absorption length
c--------------------------------------------------
      zl=1.
      it=itnucr
      plabcr=plab*plabco
ccc          call nizl(it,anuc,plabcr,sinuc,zl)
c for nizl-call replace "ccc" by " "
c--------------------------------------------------
c*** calculated is the crosssection sinuc, the absorption length zl,
c*** fo  hyperoCn inclusion: use effective labmomentum and correct the cr
c*** section and absorption length
c--------------------------------------------------
ccc          in=ins(n)
ccc          sinuc=sinuc/(sico+1.e-14)
ccc          zl=zl* sico
ccc          w(in)=rhoo/zl
      return
c
      end
      function akekaC(it,to,amss)
c*** nuclear cascade nucleon emission parameters
      dimension a(2),b(2)
      data a/0.11,0.1/
      data b/0.21,0.20/
      if (to-10.) 1,2,2
    1 continue
      akekaC=(1.-0.001*amss)*(a(it)+0.01*to)
      return
    2 continue
      akekaC=b(it)*(1.-0.001*amss)
      return
      end
c$$$  subroutine altrafC(ga,bga,cx,cy,cz,cod,cof,sif,pc,ec,p,px,py,pz,e)
      subroutine altrafC(ga,bga,cx,cy,cz,dx, dy, dz, pc,ec,p,px,py,pz,e)
c
c--------------------------------------------------
c*** (s1 in s2 / particle in s1 / result: particle in s2)
c*** arbitrary lorentz transform
c*** si are the different lorentz-systems
c--------------------------------------------------
c$$$$$$$$$$
c         write(*,*) ' cx,cy,cz in altrafC=',cx,cy,cz
c         write(*,*) ' dx,dy,dz in altrafC=',dx, dy, dz
c         write(*,*) ' pc=',pc,' ec=',ec
c$$$$$$$$$
      bgx=bga*cx
      bgy=bga*cy
      bgz=bga*cz
c$$$$$$$$$$$$$
c       write(*,*) ' bgx,y,z=',bgx, bgy, bgz
c$$$$$$$$$$$
c$$   cod2=cod*cod
c$$   if (cod2.gt.0.999999) cod2=0.999999
c$$   sid=sqrt(1.-cod2)*pc
c$$   pcx=sid*cof
c$$   pcy=sid*sif
c$$   pcz=cod*pc
      pcx=pc*dx
      pcy=pc*dy
      pcz=pc*dz
      ep=pcx*bgx+pcy*bgy+pcz*bgz
      pe=ep/(ga+1.)+ec
      px=pcx+bgx*pe
      py=pcy+bgy*pe
      pz=pcz+bgz*pe
      p=sqrt(px*px+py*py+pz*pz)
      if(p .eq. 0.) then
         px=0.
         py=0.
         pz=1.
      else
         pm=1./p
         px=px*pm
         py=py*pm
         pz=pz*pm
      endif
      e=ga*ec+ep
c$$$$$$$$$$$$
c        write(*,*) ' p=',p, ' e=',e
c$$$$$$$$$$
      return
      end
      function ankekaC(it,to,amss,sqamss )
c*** cascade parameter calculation
      dimension a(2),b(2),c(2),d(2)
      data a/1.,1.3/
      data b/.15,.225/
      data c /1.35,2.28/
      data d/.3,.4/
      if(to-0.1) 1,1,2
    1 continue
      ankekaC=sqamss*b(it)*0.1
      return
    2 continue
      if (to-10.) 3,4,4
    3 continue
      ankekaC=0.1*sqamss*(0.5+a(it)*(1.+log10(to))**2)
     **d(it)
      return
    4 continue
      ankekaC=sqamss*0.1*c(it)
      return
      end
      subroutine chanwnC
      integer * 2 ich,ibar,k1,k2,nzk,nrk
     *,ieii,ikii,nure
      common/abltisC/am(110),ga(110),tau(110),ich(110)
     *,ibar(110),k1(110),k2(110)
      common/splitC/ nzk(460,3),wt(460)
      dimension hwt(460)
      common /redverC/ irii(17),ikii(17),ieii(17)
     *,thresh( 268)
      dimension hwk(40)
      common /reacC/umo( 296),plabf( 296),siin( 296),wk( 5184),
     *nrk(2, 268),nure(30,2)
      dimension si(5184)
      equivalence (wk(1),si(1))
c*** weights for the sampling procedure (added one to each other in
c*** corresp. channels) specific for nucrinC only
c*** calculation of threshold energy of the reaction channels
c
      ireg=16
      do 222 ire=1,ireg
      iwko=irii(ire)
      iee=ieii(ire+1)-ieii(ire)
      ike=ikii(ire+1)-ikii(ire)
      ieo=ieii(ire)+1
      iika=ikii(ire)
      do 221 ie=1,iee
      sis=1.e-14
      sinorc=0.1
      do 223 ik=1,ike
      iwk=iwko+iee*(ik-1)+ie
      if (nrk(2,iika+ik).eq.0) sinorc=1.
      if (sinorc.lt.1..and.ie.gt.3) si(iwk)=si(iwk)*0.1**(ie-3)
      if (sinorc.lt.1..and.ie.gt.8) si(iwk)=        1.e-5
  223 sis=sis+si(iwk)
     **sinorc
      siin(ieo+ie-1)=sis
      sio=0.
      if (sis.ge.1.e-12) go to 2231
      sis=1.
      sio=1.
 2231 continue
      sinorc=0.1
      do 224 ik=1,ike
      if (nrk(2,iika+ik).eq.0) sinorc=1.
      iwk=iwko+iee*(ik-1)+ie
      sio=sio+si(iwk)/sis
     **sinorc
  224 hwk(ik)=sio
      do 225 ik=1,ike
      iwk=iwko+iee*(ik-1)+ie
  225 wk(iwk)=hwk(ik)
      iiki=ikii(ire)
      do 226 ik=1,ike
      am111=0.
      inrk1=nrk(1,iiki+ik)
      if (inrk1.gt.0) am111=am(inrk1)
      am222=0.
      inrk2=nrk(2,iiki+ik)
      if (inrk2.gt.0) am222=am(inrk2)
      thresh(iiki+ik)=am111 +am222
      if (inrk2-1.ge.0) go to 227
      inrkk=k1(inrk1)
      amss=5.
      inrko=k2(inrk1)
      do 228 inrk1=inrkk,inrko
      inzk1=nzk(inrk1,1)
      inzk2=nzk(inrk1,2)
      inzk3=nzk(inrk1,3)
      ams=am(inzk1)+am(inzk2)
      if (inzk3-1.ge.0) ams=ams+am(inzk3)
      if (amss.gt.ams) amss=ams
  228 continue
      ams=amss
      if (ams.lt.umo(ieo)) ams=umo(ieo)
      thresh(iiki+ik)=ams
  227 continue
  226 continue
  221 continue
  222 continue
      do 3 j=1,460
    3 hwt(j)=0.
      do 1 i=1,110
      ik1=k1(i)
      ik2=k2(i)
      hv=0.
      do 2 j=ik1,ik2
      hv=hv+wt(j)
      hwt(j)=hv
      ji=j
    2 continue
      if (abs(hv-1.).gt.1.e-4)write(6,101)
  101 format(45h error in hwt because of false use of chanwnC)
    1 continue
      do 4 j=1,460
    4 wt(j)=hwt(j)
      return
      end
      function eexiC(ij,e,a)
c-----------------------------------------------
c     nuclear excitation energy including all low energy secondaries in
c     collisions of particle ij of energy e gev on nuclei a
c-----------------------------------------------
      if (e.le.0.125) go to 20
      b=sqrt(a)/9.-.2
      if (b.gt.a*0.01) b=a*0.01
      if (e.ge.3.) go to 10
      if (b.lt.0.125) b=0.125
      eexiC=0.125+(e-0.125)*(b-0.125)/2.875
      return
   10 continue
      eexiC=b
      return
   20 continue
      eexiC=e
      return
      end
      function ekekaC(ix,to,amss,sqamss)
c*** calculation of average cascade and excitation energy
c*****ix=1 eev    2epk    3 enk     4 eex=epk+eev    5 eext=eex+enk
      go to (1,2,3,1,1),ix
    1 continue
      if(to-0.1) 11,11,12
   11 continue
      aa=0.001*sqamss
      go to 19
   12 continue
         apar=0.035
      bpar=3.
      cpar=0.1
      aa=cpar*sqamss*(0.01+apar*(bpar+log10(to))**2)
   19 continue
      if(ix.gt.3) go to 2
      ekekaC=aa
      return
    2 continue
      an=ankekaC(1,to,amss,sqamss)
      a=akekaC(1,to,amss )
      extoa=0.
      if(to.lt.5.*a) extoa=exp(-to/a)
      tpkav=a*(1.-(to/a+1.)*extoa)/(1.-extoa)
      bb=tpkav*an
      if(ix.eq.4) go to 4
      if(ix.eq.5) go to 3
      ekekaC=bb
      return
    4 continue
      ekekaC=aa+bb
      return
    3 continue
      an=ankekaC(2,to,amss,sqamss)
      a=akekaC(2,to,amss)
      extoa=0.
      if(to.lt.5.*a) extoa=exp(-to/a)
      tnkav=a*(1.-(to/a+1.)*extoa)/(1.-extoa)
      cc=tnkav*an
      if(ix.eq.5) go to 5
      ekekaC=cc
      return
    5 continue
      ekekaC=aa+bb+cc
      return
      end
      subroutine ferhadC(anuc,znuc,it,plab,elab,cx,cy,cz,itta)
c
c--------------------------------------------------
c***  collision of hadron it with nucleon itta, it has labenergy elab,
c*** momentum plab, directions cx,cy,cz, the nucleus has anuc nucleons
c*** and znuc protons
c*** itta gets a fermimomentum
c*** conserved is the energy, the momentum, electric and baryon. charge
c*** and strangeness
c--------------------------------------------------
      common /percoC/ iperco

      common /kinpaC/amta,amit,pxta,pyta,pzta,psx,psy,psz,pges2,pges,
     *pgesm,pta2,ptam,umo2,ccx,ccy,ccz,cxta,cyta,czta,ga,bga
c      dimension pxyz(10),e(2),ferf(240)
      dimension pxyz(10),e(2)
      integer * 2 ich,ibar,k1,k2
      common /finlspC/ ir,itr(20),cxr(20),cyr(20),czr(20),elr(20),
     * plr(20)
      common/abltisC/am(110),gaa(110),tau(110),ich(110),ibar(110),
     * k1(110) ,k2(110)

      common /fermtC/elabke,a,ifert,tpnvk
     *,elate

      data itesv/0/
      data pxyz/10*0./

c$$$$$$$$$$$$$$
c       write(*,*) ' elab, plab in ferhadC=',elab, plab
c$$$$$$$$$$
      elate=elab
      ia=anuc
      ir=0
      iz=znuc
      ifert=-1
c$$$  az=1.
c
c--------------------------------------------------
c*** calculation of fermimomentum
c*** ferp,fern = upper limits of proton or neutron momenta
c--------------------------------------------------
      if (ia.le.iz) ia=iz+1
      an=0.4
      if (ia.le.1) an=0.
      ferp=an*(znuc/anuc)**0.3333333
      fern=an*((anuc-znuc)/anuc)**0.3333333
      ihadc=0
  21  continue
      ihadc=ihadc+1
      ir=1
      itr(1)=it
      cxr(1)=cx
      cyr(1)=cy
      czr(1)=cz
      elr(1)=elab
      plr(1)=plab
c$$$  if (ihadc.gt.20) return
      if (ihadc.gt.10) return
      nfer=0
      ferm=ferp
c
c--------------------------------------------------
c***  absolut momentum choice, distribution propor. momentum **2
c--------------------------------------------------
      if (itta.eq.8) ferm=fern
c
c--------------------------------------------------
c*** greatest of 3 random numbers
c--------------------------------------------------
      call frmicC(p2)
c$$   if (anuc.le.1.) ferm=0.000001
c$$   if (iperco.gt.0) p2=0.0000001
      if (anuc.le.1.) ferm=0.0
      if (iperco.gt.0) p2=0.0
      p2=ferm*p2
c$$$$$$$$$$$$$$$$
c      p2=0.
c       write(*,*) ' anuc =',anuc,' iperco=',iperco, ' ferm=',ferm
c       if(iperco .gt. 0.) then
c            p2=0.
c       else
c            p2=1.
c       endif
c       p2=p2*ferm
c$$$$$$$$$$$$
c--------------------------------------------------
c***  angle choice
c--------------------------------------------------
    1 continue
      call poliC(polc,pols)
      call cosiC(sfe,cfe)
c$$$$$$$$$$$$$$$$$
c        polc= 1.
c        pols=0.
c$$$$$$$$$$$$$$$$
c
c--------------------------------------------------
c***  fermimomentum components pxyz(2,6,10),fermi-energy e(2) of itta
c--------------------------------------------------
c$$$$$$$$$$$ added
        cxta=pols*cfe
        cyta=pols*sfe
        czta=polc
c$$$$$$$$$$$
c$$$  aaa=p2*pols
c$$$  if (nfer.lt.10.and.ihadc.lt.10) go to 12
      if (nfer.lt.5 .and.ihadc.lt. 5) go to 12
c*** passage for a higher run speed at low energies
      if(plab .lt. ferm) then
         iplab=plab*100.+10
         ifpla=ferm*100.+10
         ferml=ferm
         if(itta .eq. 8) then
              ferm=fern*iplab/ifpla
         else
              ferm=ferp*iplab/ifpla
         endif
c$$      if (plab.lt.ferml) ferm=ferp*iplab/ifpla
c$$      if (plab.lt.ferml.and.itta.eq.8) ferm=fern*iplab/ifpla
      endif
12    continue
c$$   pxyz(2)=aaa*cfe
      pxyz(6)=p2*cxta
c$$   pxyz(2)=aaa*cfe
      pxyz(6)=p2*cyta
c$$$  pxyz(10)=p2*polc
      pxyz(10)=p2*polc
c
c--------------------------------------------------
c*** collision kinematics
c--------------------------------------------------
      e(2)=sqrt(p2*p2+am(itta)**2)
c$$$$ rlke=e(2)-am(itta)
      amta=am(itta)
      amit=am(it)
      tlab=elab-amit
c$$$$ pot=anuc*0.007
c$$$$ ipot=ipot+1
      e(1)=elab
      pxyz(1)=cx*plab
      pxyz(5)=cy*plab
      pxyz(9)=cz*plab
      px=pxyz(1)+pxyz(2)
      py=pxyz(5)+pxyz(6)
      pz=pxyz(9)+pxyz(10)
c
c--------------------------------------------------
c*** calcul. of kinem. parameters
c--------------------------------------------------
c$$$  call kinparC(pxyz,e,0.001)
      call kinparC(pxyz,e)
c
c--------------------------------------------------
c*** lorentz-transform into target-nucleon-rest-system
c--------------------------------------------------
      ga=e(2)/amta
      bga=-sqrt(ga*ga-0.999999)
c$$$$$$$$$$$
c      write(*,*) ' amta=',amta,' ga=',ga,' bga=',bga
c$$$$$$$$$$$
c$$$  cod=cz
c$$$  sid=sqrt(1.-cod*cod+1.e-12)
c$$$  if (sid.le.0.) sid=1.e-10
c$$$  sidm=1./sid
c$$$  sif=cy*sidm
c$$$  cof=cx*sidm
c
c--------------------------------------------------
c***  momentum limitation of the projectile hadron in target nucleon res
c*** system, plabs = momentum of the hadron, elabs = energy of the hadro
c*** in target nucleon rest system
c--------------------------------------------------
      call altrafC(ga,bga,cxta,cyta,czta,cx, cy, cz, plab,elab,
     *plabs,cxct,cyct,czct,elabs)
      plabou=15.
c$$$  if (elabs.gt.plabou) az=0.
      if (elabs.gt.plabou) elabs=plabou
c      ******
      temp = elabs*elabs-amit*amit
      if(temp .le. 0.) then
         plabs = 0.0001
      else
         plabs=sqrt(temp)
      endif
c
c--------------------------------------------------
c*** h-n-interaction in target rest system (t.r.s.)
c***  of the target nucleon
c*** (elabt = invariant kinetic h-n-energy+mass(h))
c*** if (elabt is grater then the available invariant projectile energy
c*** in the h-a-system, go back to fermimomentum choice
c*** counting index = nfer
c*** if (nfer.gt.10 or 20 , start newprojectile labmomentum choice in
c*** nucrinC, else interaction simulation
c--------------------------------------------------
      elabt=sqrt(umo2)-amta
      nfer=nfer+1
c     if (nfer.gt.20) go to 21
      if (nfer.gt. 5) go to 21
      if (elabt.gt.elabke) go to 1
      ifert= 1
c
c--------------------------------------------------
c***  for particles of the h-n-collision, store the kinem.variables in
c*** common /finlspC/
c--------------------------------------------------
      call hadrinC(it,plabs,elabs,cxct,cyct,czct,itta)
c
c--------------------------------------------------
c*** lorentz-transform from trs into ls
c--------------------------------------------------
      elate=0.
c$$$  ga=e(2)/amta
c$$$  bga=sqrt(ga*ga-1.+1.e-6)
      bga=-bga
      do 11 iii=1,ir
          itp=itr(iii)
          ami=am(itp)
          el=elr(iii)
          pl=plr(iii)
c$$$      cod=czr(iii)
c$$$      sid=sqrt(1.-cod*cod+1.e-12)
c$$$      sidm=1./sid
c$$$      sif=cyr(iii)*sidm
c$$$      cof=cxr(iii)*sidm
          call altrafC(ga,bga,cxta,cyta,czta, cxr(iii),
     *    cyr(iii), czr(iii),
     *    pl,el,
     *    plr(iii),cxr(iii),cyr(iii),czr(iii),elr(iii))
          elate=elate+elr(iii)
          px=px-plr(iii)*cxr(iii)
          py=py-plr(iii)*cyr(iii)
          pz=pz-plr(iii)*czr(iii)
11    continue
      nfe=nfer+1
      e2=sqrt(p2**2+amta**2)
c
c--------------------------------------------------
c***  if in the lab system the sampled event has more energy then the
c*** primary hadron in h-a-collision, start again with fermi mom. choice
c--------------------------------------------------
      if (elate-e(2).gt.a) go to 1
      if (ir.gt.1) elate=elate-am(itta)
      return
      end
      subroutine frmicC(gpart)
      dimension g(3)
      real*8 u
      do 5 i=1,3
      call rndc(u)
      g(i) = u
5     continue
c
c--------------------------------------------------
c*** find largest of 3 random numbers
c--------------------------------------------------
      if (g(3).lt.g(2)) go to 25
      if (g(3 ).lt.g(1)) go to 20
      gpart=g(3)
10    return
20    gpart=g(1)
      go to 10
25    if (g(2).lt.g(1)) go to 20
      gpart=g(2)
      go to 10
      end
      subroutine hyperoC(it,itnucr,sico,plabco)
      common /abltisC/ am(330),ich(110),ibar(110),k(220)
      integer*2 ich,ibar,k,ihatrs,inutrs
c      dimension ihatrs(30),inutrs(30),sicor(30)
      dimension ihatrs(30),inutrs(30)
      data ihatrs/16*30,25,24,30,16,15,25,8*30/
c
c--------------------------------------------------
c*** attension, in case of it=sig+=21 no strangenes-cons.,differ.1unit
c--------------------------------------------------
      data inutrs/16*30,8,2,30,8,1,8,8*30/
      itnucr=inutrs(it)
      ithacr=ihatrs(it)
      plabco=1.
      ithatr=ithacr
      if (itnucr.ge.30) ithatr=it
      if (itnucr.ge.30) itnucr=it
      sico=1.
      if (ithatr.lt.30) sico=42/39.
      if (it.eq.21) sico=42/40.
      if (ithatr.lt.30) plabco=am(it)/(am(itnucr)+1.e-10)
      it=ithatr
      return
      end
c$$$$ subroutine kinparC(pxyz,e,rcpmv)
      subroutine kinparC(pxyz,e)
      common /kinpaC/amta,amit,pxta,pyta,pzta,psx,psy,psz,pges2,pges,
     *pgesm,pta2,ptam,umo2,ccx,ccy,ccz,cxta,cyta,czta,ga,bga
c*** calculation of kinematic parameters
      dimension pxyz(10),e(2)
      pxta=pxyz(  2)
      pyta=pxyz(  6)
      pzta=pxyz( 10)
c$$$$$$$$$$$$
c$$      write(*,*) ' pxta,..z=',pxta, pyta, pzta
c$$$$$$$$
      psx=pxyz( 1)+pxta
      psy=pxyz( 5)+pyta
      psz=pxyz( 9)+pzta
      pges2=psx*psx+psy*psy+psz*psz
      pges=sqrt(pges2)
      pgesm=1./pges
c$$$  umo2=((e(1)+e(2))**2-pges2)/(rcpmv*rcpmv)
      umo2=((e(1)+e(2))**2-pges2)
c     umo2=umo2*0.000001
      ccx=psx*pgesm
      ccy=psy*pgesm
      ccz=psz*pgesm
c$$$  pta2=pxta*pxta+pyta*pyta+pzta*pzta
c$$   if(pta2 .eq. 0.) then
c$$       cxta=0.
c$$       cyta=0.
c$$       czta=1.
c$$   else
c$$       ptam=1./sqrt(pta2)
c$$       cxta=pxta*ptam
c$$       cyta=pyta*ptam
c$$       czta=pzta*ptam
c$$$  endif
c$$$  ga=(e(1)+e(2))*0.001/(rcpmv*sqrt(umo2))
c$$$  bga=sqrt(ga*ga-1.)
      return
      end
      function pofekC(ekin,it)
c*** nucleon/pion mass calculation
      dimension am(3)
      data am/0.9382,0.9382,0.13959/
      pofekC=am(it)
      if(ekin.gt.3.*pofekC)go to 100
      pofekC=sqrt((ekin+pofekC)**2-pofekC**2)
      return
  100 continue
      pofekC=ekin+pofekC
      return
      end
      subroutine poliC(cso,sio)
c
c--------------------------------------------------
c*** random choice of angle teta (cs = cos(teta),si = sin(teta)
c--------------------------------------------------
c###############
      real*8  u, cs
c#################
      call rndc(u)
      call rndc(cs)
      if (u.lt.0.5) cs=-cs
      sio = sqrt(1.-cs*cs+1.e-10)
      cso = cs 
      end
      subroutine ranchaC(u,ir)
      real * 8 ir
      ir=u*1.e+9
      return
      end
      subroutine rakekaC(it,to,amss,sqamss,t,p,de)
c###############
      real*8  v, u
c############
c*** sampling of cascade nucleon angle and energy
      ap=akekaC(it,to,amss)
      tmin=0.020
      ap=ap-(ap+tmin)*0.15
      call rndc(v)
      if (v.lt.0.6) ap=ap/6.
    1 continue
c          below is to produce a random variable t with
c          exp(-app t)dt. But   this may result in
c          semi infinite loop when a very small
c          random number is produced in rexpC.
c          So the method is replaced by usual one.
c          K.K 98/03/05
c      app=1./ap
c      t=rexpC(app)
c
      call rndc(u)
      t = -ap * log(u)
c
      t=t + tmin
      if(t.gt.to) go to 1
      p=pofekC(t,it)
      ade=0.090*(1.+0.081*amss**.3333)
     */t
      dex=exp(-1.57*1.57/ade)
      an1=(1.-dex)*ade/2.
      an2=dex*1.57
      an=an1+an2
      an1=an1/an
      call rndc(v)
      if(v.gt.an1)  go to 3
    2 continue
      call rndc(v)
      de=sqrt(-ade*log(1.-v*(1.-dex)))
      if(de.gt.1.57) go to 2
      return
    3 continue
      call rndc(v)
      de=-v
      de=atan2(sqrt(1.-de**2),de)
      return
      end
