c       program testkdecay
c      implicit real*8 (a-h,o-z)
cc      parameter (aMk0=0.497672d0, aMkc=0.493677d0, aMpi0=0.1349764d0, 
cc     &  aMpic=0.13956995d0,aMmu=0.105658389d0, aMe=0.00051099907d0,
cc     &  aMnu=0.d0)
c      parameter (kep=1,kem=2,kgamma=3,
c     &     kkp=4,kkm=5,kk0l=6,kk0s=7,
c     &     kneut=8,kneutbar=9,kpro=10,kprobar=11,
c     &     kpip=12,kpim=13,kpi0=14,
c     &     kmup=15, kmum=16, knue=17, knuebar=18, knumu=19,
c     &     knumubar=20)
c
c      character*4 cptl(20)
c      integer*4 kchild(3)
c      real*8 pchild(5,3), spinmu(5), pmother(5), zeta(3), sum(4),
c     &     am(20)
c      data cptl/'E+  ', 'E-  ','GAMM','K+  ','K-  ','K0L ','K0S ',
c     &'N+  ','N-  ','P+  ','P-  ','PI+ ','PI- ','PI0 ','MU+ ',
c     &'MU- ','NUE ','NUEB','NUM ','NUMB'/
c      data am/2*0.51099907d-3,0.0,2*0.493677d0,2*0.497672d0,
c     &   2*0.93956563d0,2*0.93827231d0,2*0.13956995d0,0.1349764d0,
c     &   2*0.105658389d0,4*0.0/
c
c      data zeta/0,0.8660254,0.5/
c      read(*,*) kin0, e0
c
c      do 100 i=1,4
c         sum(i)= 0.
c 100  continue
c      nsum = 0
c      sumenu=0.
c      nsumenu=0.
c      summunu=0.
c      nsummunu=0.
c
c      do 1 i=1,1000
c         kin = kin0
c         pmother(1) = 0.
c         pmother(2) = 0.
c         pmother(3) = sqrt(e0**2 - pmother(5)**2)
c         pmother(5) = am(kin0)
c         pmother(4) = e0
c         call zeta2spin(zeta,pmother,spinmu)
c         call decaycontrol
c     &        (kin, pmother, nchild, kchild, pchild, spinmu, kdcy)
c         write(*,*) i, nchild, kdcy
c
c         do 10 k=1,nchild
c            write(*,'(a,10F12.7)') cptl(kchild(k)),(pchild(j,k),j=1,5)
c
c            if(kchild(k).eq.kmup.or.kchild(k).eq.kmum) then
c               write(*,'(a,10F12.7)') 'mu spin',(spinmu(j),j=1,4)
c     &        ,spinmu(1)**2+spinmu(2)**2+spinmu(3)**2-spinmu(4)**2
c
c               call spin2zeta(spinmu, pchild(1,k), zeta)
c               write(*,'(a,10F12.7)') 'at rest',zeta 
c     &              ,zeta(1)**2+zeta(2)**2+zeta(3)**2
cc     if mu, succeed to mudecay
c               do 21 j=1,5
c                  pmother(j)=pchild(j,nchild)
c 21            continue
c               kin = kchild(nchild)
c               call decaycontrol
c     &       (kin, pmother, nchild, kchild, pchild, spinmu, kdcy)
c               do 20 l=1,nchild
c                  write(*,'(a,10F12.7)') 
c     &                 cptl(kchild(l)),(pchild(j,l),j=1,5)
c                  do 101 j=1,4
c                     sum(j) = sum(j) + pchild(j,l)
c 101              continue
c 20            continue
c            else
c               if(kchild(k).eq.knue.or.kchild(k).eq.knuebar) then
c                  sumenu=sumenu+pchild(4,k)
c                  nsumenu=nsumenu+1
c               else if(kchild(k).eq.knumu.or.kchild(k).eq.knumubar) then
c                  summunu=summunu+pchild(4,k)
c                  nsummunu=nsummunu+1
c               end if
c               do 102 j=1,4
c                  sum(j) = sum(j) + pchild(j,k)
c 102           continue
c            end if
c 10      continue
c         nsum = nsum + 1
c 1       continue
c      write(*,*) nsum, (sum(j)/nsum, j=1,4)
c      write(*,*)  nsumenu,sumenu/nsumenu,nsummunu,summunu/nsummunu
c      end
cccc 1000  0.000611414343  0.00118362141  1.00214135  1.11746537
cc
      subroutine decaycontrol(kin, pmother, nchild, kchild, pchild, 
     &     spinmu, kinddecay)
      implicit real*8 (a-h,o-z)
c Control particl decay at rest and boost to lab frame:
c inputs :
c      kin: kind of particle
c      pmother (1-5): momentum energy, mass
c outputs :
c      nchild : number of decay particle (=3, fixed)
c     kchild(3) : kind
c     pchild(3) : momentum
c     spinmu: (1-5) spacial spin vector of muon (Note, this is input for muon)
c       for mu-, ture direction, but for mu+ inverse direction.
c       I.E. in both case, nue tends to run toward 'smu' direction
c     kinddecay : (See the comments in following code. Too many)

      parameter (aMk0=0.497672d0, aMkc=0.493677d0, aMpi0=0.1349764d0, 
     &  aMpic=0.13956995d0,aMmu=0.105658389d0, aMe=0.00051099907d0,
     &  aMnu=0.d0, aMg=0.d0)
      parameter (kep=1,kem=2,kgamma=3,
     &     kkp=4,kkm=5,kk0l=6,kk0s=7,
     &     kneut=8,kneutbar=9,kpro=10,kprobar=11,
     &     kpip=12,kpim=13,kpi0=14,
     &     kmup=15, kmum=16, knue=17, knuebar=18, knumu=19,
     &     knumubar=20)
      integer positive, negative
      parameter(positive=+1, negative=-1)

      integer*4 kkin, nchild, kchild(3)
      real*8 pchild(5,3), spinmu(5)
      real*8 p0(5), p1(5), spinmu0(5), telement(5)

c                           get decay particle at rest frame of mother

      if(kin.eq.kmup .or. kin.eq.kmum) then
c                                                 mu -> numu nue e
         call mudecay(kin, pmother, spinmu, nchild, kchild, pchild,
     &        kinddecay) 
      else if(kin.eq.kpip .or. kin.eq.kpim .or. kin.eq.kpi0) then
c                                                 pic --> mu + neu
c                                                 pi0 --> g + g
        call pidecay(kin, pmother, nchild, kchild, pchild, spinmu0,
     &     kinddecay)
      else if(kin.eq.kk0s) then
c                                           K0s -> 2pic, 2pi0
         call k0sdecay(nchild, kchild, pchild, kinddecay)
      else if(kin.eq.kk0l) then
c                                           K0l -> 3l nu, 3pi
         call k0ldecay(nchild, kchild, pchild, spinmu0, kinddecay)
      else if(kin.eq.kkp) then
c                                     K+ -> mu+nu, 3lnup, 3pi
         call kcdecay(positive, nchild, kchild, pchild, spinmu0, 
     &        kinddecay)
      else if(kin.eq.kkm) then
c                                     K- -> mu-nu, 3lnup, 3pi
         call kcdecay(negative, nchild, kchild, pchild, spinmu0,
     &        kinddecay)
      else
         write(0,*)' This particl will not decay : ',kin
         return
      endif
      if(nchild.gt.3) then
         stop ' too many child in kdecay'
      end if
c
c                           lorenz boost to mother particle moving frame
      call mom2telement(pmother, telement)
      do 1 i=1, nchild
         do 10 j=1, 5
            p0(j) = pchild(j,i)
 10      continue

         call lboost(telement, p0, p1)
         do 11 j=1, 5
            pchild(j,i) = p1(j)
 11      continue
         if(p1(4).lt.0) then
            write(*,*) telement
            write(*,*) p0
            write(*,*) p1
         end if
         if(kchild(i).eq.kmup .or. kchild(i).eq.kmum) then
            call lboost(telement, spinmu0, spinmu)
         end if
 1    continue
      end

      subroutine k0sdecay(nchild, kchild, pchild,
     &     kinddecay)
      implicit real*8 (a-h,o-z)
      parameter (aMk0=0.497672d0, aMkc=0.493677d0, aMpi0=0.1349764d0, 
     &  aMpic=0.13956995d0,aMmu=0.105658389d0, aMe=0.00051099907d0,
     &  aMnu=0.d0)
      parameter (kep=1,kem=2,kgamma=3,
     &     kkp=4,kkm=5,kk0l=6,kk0s=7,
     &     kneut=8,kneutbar=9,kpro=10,kprobar=11,
     &     kpip=12,kpim=13,kpi0=14,
     &     kmup=15, kmum=16, knue=17, knuebar=18, knumu=19,
     &     knumubar=20)

c control k0s-decay:
c no input,
c output  n : number of child particles
c    kchild(n) : kind of child particles
c    pchild(5,n) : 1-3 momentum, 4 energy, 5 mass of child particle
c
c            k0s decay
c  -- process --
c         1)  ---->pi+ + pi-       68.61% (kinddecay=1)
c         2)  ---->pi0 + pi0       31.39  (kinddecay=2)
      integer*4 kchild(2)
      real*8 pchild(5,2), p1(5), p2(5)

      call rndc(u1)
      if(u1.lt. 0.6861) then
c                                            1)  ---->pi+ + pi-       68.61%
         kinddecay = 1
         nchild=2
         call sampdecay2(aMk0, aMpic, aMpic, p1, p2)
         do 101, i=1, 5
            pchild(i,1)=p1(i)
            pchild(i,2)=p2(i)
 101     continue
         
         kchild(1)=kpip
         kchild(2)=kpim
      else
c                                            2)  ---->pi0 + pi0       31.39
         kinddecay = 2
         nchild=2
         call sampdecay2(aMk0, aMpi0, aMpi0, p1, p2)
         do 201, i=1, 5
            pchild(i,1)=p1(i)
            pchild(i,2)=p2(i)
 201     continue
         kchild(1)=kpi0
         kchild(2)=kpi0
      end if
      end

      subroutine k0ldecay(nchild, kchild, pchild, spinmu,
     &     kinddecay)
      implicit real*8 (a-h,o-z)
      parameter (aMk0=0.497672d0, aMkc=0.493677d0, aMpi0=0.1349764d0, 
     &  aMpic=0.13956995d0,aMmu=0.105658389d0, aMe=0.00051099907d0,
     &  aMnu=0.d0)
      parameter (kep=1,kem=2,kgamma=3,
     &     kkp=4,kkm=5,kk0l=6,kk0s=7,
     &     kneut=8,kneutbar=9,kpro=10,kprobar=11,
     &     kpip=12,kpim=13,kpi0=14,
     &     kmup=15, kmum=16, knue=17, knuebar=18, knumu=19,
     &     knumubar=20)

c control k0l-decay:
c no input,
c output  n : number of child particles
c    kchild(n) : kind of child particles (mu is the last(3) ptl if exist.)
c    pchild(5,n) : 1-3 momentum, 4 energy, 5 mass of child particle
c    spinmu(5) : if mu exist in child, give 4-dim spin vector(5=0)
c            k0l decay
c  -- process --
c         1)  ---->e  pi neue   38.7 %(kinddecay=1)
c         2)  ---->mu pi neum   27.1  (kinddecay=2, k0==>mu+, k0bar==>mu-)
c         3)  ---->3 pi0        21.5  (kinddecay=3)
c         4)  ---->pi+ pi- pi0  12.4  (kinddecay=4)
      integer*4  kchild(3)
      real*8 pchild(5,3), spinmu(5), p1(5), p2(5), p3(5)

      call rndc(u)
      if(u .lt. .387) then
c                                                  1) k0l -> e  + neue + pi
         kinddecay = 1
         nchild=3
         call k03e(p3, p2, p1, pol)
c                                                   3:e,   2:nue,   3:pi
         do 101, i=1, 5
            pchild(i,1)=p1(i)
            pchild(i,2)=p2(i)
            pchild(i,3)=p3(i)
 101     continue
         call rndc(u)
         if(u.gt.0.5) then
            kchild(1)=kpim
            kchild(2)=knue
            kchild(3)=kep
         else
            kchild(1)=kpip
            kchild(2)=knuebar
            kchild(3)=kem
         end if
      elseif(u.lt. .658) then
c                                                  2)  mu + neumu + pi
         kinddecay = 2
         nchild=3
         call k03mu(p3, p2, p1, pol)
c                                                  3: mu, 2: munu, 1: pi 
         do 201, i=1, 5
            pchild(i,1)=p1(i)
            pchild(i,2)=p2(i)
            pchild(i,3)=p3(i)
 201     continue
         call rndc(u)
         if(u.gt.0.5) then
c                                                 2a) treat as k0
            kchild(1)=kpim
            kchild(2)=knumu
            kchild(3)=kmup
         else
c                                                 2b) treat as k0bar
            kchild(1)=kpip
            kchild(2)=knumubar
            kchild(3)=kmum
         end if
         call samplespin(pol, p3, spinmu)
      elseif(u .lt. .873) then
c                                                  3)  3 pi0
         kinddecay = 3
         nchild = 3
         pchild(5,1)=aMpi0
         pchild(5,2)=aMpi0
         pchild(5,3)=aMpi0
         call cnbdcy(3, aMk0, pchild, 0, w, icon)
         if(icon.ne.0) then
            stop 'at k0ldecay, case 3. Icon>0'
         end if
         kchild(1)=kpi0
         kchild(2)=kpi0
         kchild(3)=kpi0
      else
c                                                  4)  pi+ pi- pi0
         kinddecay = 4
         nchild = 3
         pchild(5,1)=aMpic
         pchild(5,2)=aMpic
         pchild(5,3)=aMpi0
         call cnbdcy(3, aMk0, pchild, 0, w, icon)
         if(icon.ne.0) then
            stop 'at k0ldecay, case 4. Icon>0'
         end if
         kchild(1)=kpip
         kchild(2)=kpim
         kchild(3)=kpi0
      endif
      end

      subroutine kcdecay(motherc, nchild, kchild, pchild, spinmu,
     &     kinddecay)
      implicit real*8 (a-h,o-z)
      parameter (aMk0=0.497672d0, aMkc=0.493677d0, aMpi0=0.1349764d0, 
     &  aMpic=0.13956995d0,aMmu=0.105658389d0, aMe=0.00051099907d0,
     &  aMnu=0.d0)
      parameter (kep=1,kem=2,kgamma=3,
     &     kkp=4,kkm=5,kk0l=6,kk0s=7,
     &     kneut=8,kneutbar=9,kpro=10,kprobar=11,
     &     kpip=12,kpim=13,kpi0=14,
     &     kmup=15, kmum=16, knue=17, knuebar=18, knumu=19,
     &     knumubar=20)
c     control kc-decay:
c   input,
c     motherc : kaon charge
c   outputs;
c     nchild : number of child particles
c     kchild(nchild) : kind of child particles
c     pchild(5,nchild) : 1-3 momentum, 4 energy, 5 mass of child particle
c     spinmu(5) : muon spin if there is muon in children
c      spin direction for mu+ and antispin direction for mu-
c      I.e. zeta represent the direction for which numu and e are 
c      emitted less energetically, when mu decays.
c            k+- decay
c  -- process --
c         1) k---->mu+neu          63.5 % (kinddecay = 1)
c         2)  ---->pic + pi0       21     (kinddecay = 2)
c         3)  ---->pi0+e+neu       4.8    (kinddecay = 3)
c         4)  ---->pi0+mu+neu      3.2    (kinddecay = 4)
c         5)  ---->pic+pic+pic     5.6    (kinddecay = 5)
c         6)  ---->pic+pi0+pi0     1.7    (kinddecay = 6)
      integer*4 n, kchild(3) 
      real*8 pchild(5, 3), spinmu(5), p1(5),p2(5),p3(5)

      call rndc(u)
      if(u  .lt. 0.635) then
c                                            1) k---->mu+neu          63.5 %
         kinddecay = 1
         nchild = 2
         call sampdecay2(aMkc, aMmu, aMnu, p2, p1)
c                                            1: numu, 2:mu
         do 101, i=1, 5
            pchild(i,1)=p1(i)
            pchild(i,2)=p2(i)
 101     continue
         if(motherc.gt.0) then
            kchild(1)=knumu
            kchild(2)=kmup
         else
            kchild(1)=knumubar
            kchild(2)=kmum
         end if
         pol = 1.0d0
         call samplespin(pol, p2, spinmu)
      elseif(u .lt. .845) then
c                                            2)  ---->pic + pi0       21
         kinddecay = 2
         nchild = 2
         call sampdecay2(aMkc, aMpic, aMpi0, p1, p2)
         do 201, i=1, 5
            pchild(i,1)=p1(i)
            pchild(i,2)=p2(i)
 201     continue
         if(motherc.gt.0) then
            kchild(1)=kpip
            kchild(2)=kpi0
         else
            kchild(1)=kpim
            kchild(2)=kpi0
         end if
      elseif(u .lt. .893) then
c                                             3)  ---->pi0+e+neu       4.8
         kinddecay = 3
         nchild = 3
         call kc3e(p3, p2, p1, pol)
         do 301, i=1, 5
            pchild(i,1)=p1(i)
            pchild(i,2)=p2(i)
            pchild(i,3)=p3(i)
 301     continue
         kchild(1)=kpi0
         if(motherc.gt.0) then
            kchild(2)=knue
            kchild(3)=kep
         else
            kchild(2)=knuebar
            kchild(3)=kem
         end if
      elseif(u .lt. .925) then
c                                             4)  ---->pi0+mu+neu      3.2
         kinddecay = 4
         nchild = 3
         call kc3mu(p3, p2, p1, pol)
         do 401, i=1, 5
            pchild(i,1)=p1(i)
            pchild(i,2)=p2(i)
            pchild(i,3)=p3(i)
 401     continue
         kchild(1)=kpi0
         if(motherc.gt.0) then
            kchild(2)=knumu
            kchild(3)=kmup
         else
            kchild(2)=knumubar
            kchild(3)=kmum
         end if
         call samplespin(pol, p3, spinmu)
      elseif(u .lt. .981) then
c                                            5)  ---->pic+pic+pic     5.6
         kinddecay = 5
         nchild = 3
         pchild(5,1)=aMpic
         pchild(5,2)=aMpic
         pchild(5,3)=aMpic
         call cnbdcy(3, aMkc, pchild, 0, w, icon)
         if(icon.ne.0) then
            stop 'at kcdecay, case 5. Icon>0'
         end if
         if(motherc.gt.0) then
            kchild(1)=kpip
         else
            kchild(1)=kpim
         end if            
         kchild(2)=kpip
         kchild(3)=kpim
      else
c                                             6)  ---->pic+pi0+pi0     1.7
         kinddecay = 6
         nchild = 3
         pchild(5,1)=aMpic
         pchild(5,2)=aMpi0
         pchild(5,3)=aMpi0
         call cnbdcy(3, aMkc, pchild, 0, w, icon)
         if(icon.ne.0) then
            stop 'at kcdecay, case 6. Icon>0'
         end if
         if(motherc.gt.0) then
            kchild(1)=kpip
         else
            kchild(1)=kpim
         end if
         kchild(2)=kpi0
         kchild(3)=kpi0
      end if
      end

      subroutine mudecay(kmuin, pmu, smu, nchild, kchild, pchild,
     &     kinddecay) 
      implicit real*8 (a-h,o-z)
c control muon decay:
c inputs :
c     kmuin (should be kmup(15) or kmum(16))
c     pmu (1-5) momentum energy, mass
c     smu (1-5) spacial spin vector of muon 
c     for mu-, ture direction, but for mu+ inverse direction.
c     I.E. in both case, nue tends to run toward 'smu' direction
c outputs :
c     nchild : number of decay particle (=3, fixed)
c     kchild(3) : kind 
c     pchild(3) : momentum
c     kinddecay = 1
      parameter (aM=0.0d0)
      parameter (kep=1, kem=2,
     &     kmup=15, kmum=16, knue=17, knuebar=18, knumu=19,
     &     knumubar=20)

      integer*4 nchild, kchild(3)
      real*8 pmu(5), smu(5), pchild(5,3)
      real*8 telement(5), zeta(3), p0(5), p1(5),p2(5),p3(5)

      if(kmuin.eq.kmum) then
c                              mu- -> numu + e- + nuebar
         kchild(1)=knumu
         kchild(2)=kem
         kchild(3)=knuebar
      else if(kmuin.eq.kmup) then
c                              mu+ -> numubar + e+ + nue
         kchild(1)=knumubar
         kchild(2)=kep
         kchild(3)=knue
      else
         write(0,*) kmuin
         stop 'irregal kmuin iput to mudecay'
      end if

      nchild = 3
      kinddecay = 1
c                              muon spin direction at muon rest frame
      call spin2zeta(smu, pmu, zeta)
c                              mu-decay in muon rest frame 
      call sampmudecay(Enumu,tetnumu, Ee,tete, Enue,tetnue)
      if(Enumu.lt.0.or.Ee.lt.0.or.Enue.lt.0) then
         write(*,*) 'before sampmo', Enumu,tetnumu, Ee,tete, Enue,tetnue
      end if
c                        sample momentum of decay particle in muon rest frame,
      call sampmoment(Enumu, aM, zeta, tetnumu, p1)
      call sampmoment(Ee,    aM, zeta, tete,    p2)
      call sampmoment(Enue,  aM, zeta, tetnue,  p3)
      if(p1(4).lt.0.or.p1(4).lt.0.or.p1(4).lt.0) then
         write(*,*) 'after sampmo',Enumu,tetnumu, Ee,tete, Enue,tetnue
      end if
      do 1 i=1, 5
         pchild(i,1) = p1(i)
         pchild(i,2) = p2(i)
         pchild(i,3) = p3(i)
 1    continue
      end

      subroutine sampmoment(E, aM, zeta, costet, p)
      implicit real*8 (a-h, o-z)
c **Sample muon spin direction for given p and polarization
      dimension p(5), zeta(3), u(3), spin(5)

      p(4) = E
      p(5) = aM
      ap = sqrt(E**2 - aM**2) 

      if(costet .gt. 0.999d0) then
c                                         0-1: if costeta = 1
         do 14 i = 1, 3
            p(i) = zeta(i)*ap
 14      continue
      else
c                                      if costeta != 1, then 
c                                      1: Sample a vector perpendicular to dirp
         ru =0.0
         do 1 while(ru.lt.1e-3)
c                                         1-1: sample a vector inside of 
c                                               a unit sphere
            rr=2.0
            do 3 while (rr.gt.1.0d0 .or. rr.lt.1.0d-2)
               rr=0.0
               do 13 i=1,3
                  call rndc(u1)
                  u(i) = 2*u1 - 1.0d0
                  rr = rr + u(i)**2
 13            continue
 3          continue
c                                     1-2: substruct the component 
c                                                    parallel to zeta
            pinner = zeta(1)*u(1) + zeta(2)*u(2) + zeta(3)*u(3)
            ru = 0.0
            do 2 i=1,3
               u(i) = u(i) - zeta(i)*pinner
               ru   = ru + u(i)**2
 2          continue
c                                         1: Resulted vector shouldn't be zero 
c                                             (see condition in 'do 1 while')
 1       continue
c                                         2: sin(teta)* normalization
         tfact = sqrt((1.d0 - costet**2)/ru)
         do 4 i = 1, 3
            p(i) = ap*(zeta(i)*costet + u(i)*tfact)
 4       continue
      end if
      end

      subroutine sampmudecay(Enumu,tetnumu, Ee, tete, Enue, tetnue)
      implicit real*8 (a-h,o-z)
c spin direction is assumed z+
      parameter (aMmu=0.105658389d0,aMe=0.00051099907d0)
      real*8 zeta(3), pnumu(5), pnue(5), pe(5)

c Function form for Enumu and Ee are the same (assuming Me = 0).
c kk=1 calculate Enumu, and kk=2  Ee
      do 2000 kk=1, 2
c        sample decay angle : costet
         u2=2.1
         u1=1.0
         do 1 while (2*u2. gt. 2.d0 - u1)
            call rndc(u1)
            call rndc(u2)
 1       continue
c  Costet direction is different from Nue
         costet = 2.d0*u1 - 1.0d0
c  Then to coefs of energy spectrum of NUmu and E
c      g     = (3.d0 + costet)*x**2   - 2*(1.d0 + costet)*x**3
         coe1 =   (3.d0 + costet)
         coe2 = 2*(1.d0 + costet)

c select energy integral by a random number
         call rndc(u3)
c sample a spectrum in inegral [ g_int(1)*u3 (random number)]
         Gi0 = (coe1/3 - coe2/4)*u3
c Newtonian method to find solution for g_int(x) = gi0
         dx = 1.
         x = 0.5
         do 2 while (abs(dx) .gt. 1e-10)
            Gi1 = coe1*x**3/3 - coe2*x**4/4
            g0  = coe1*x**2   - coe2*x**3
            diff = Gi0 - Gi1
c For stability
            if(g0.lt.1d-2) then
               g0=1.d-2
            end if

            x1 = x + diff/g0
C limit x region
            if(x.lt.0.0) then
               x1 = 0.0
            else if(x.gt.1.0) then
               x1 = 1.0d0
            end if

            dx = x1 - x
            x = x1
 2       continue
         if(kk.eq.1) then
c            Enumu = ((aMmu**2 - aMe**2)/2/aMmu)*x
            Enumu = (aMmu/2)*x
            tetnumu = costet
         else
            Ee = (AMmu/2)*x
            tete = costet
         end if
 2000 continue

c For Nue
c        sample decay angle : costet
         u2=2
         u1=1
         do 11 while (2*u2. gt. 1.d0 - u1)
            call rndc(u2)
            call rndc(u1)
 11      continue
c  Costet: direction is different from Numu or E
         costet = -2.d0*u1 + 1.0d0
c  Coeficients for Nue
c         coe1 = 1.d0
c         coe2 = 1.d0

c select energy integral by a random number
         call rndc(u3)

c energy-integral(of allowed region)*u3 (random number)
c         Gi0 = (coe1/3 - coe2/4)*u3
         Gi0 = (1.d0/3 - 1.d0/4)*u3
         dx = 1.
         x = 0.5
c Newtonian method to find solution for g_int(x) = gi0
         do 12 while (abs(dx) .gt. 1e-10)
c            Gi1 = coe1*x**3/3 - coe2*x**4/4
c            g0  = coe1*x**2   - coe2*x**3
            Gi1 = x**3/3 - x**4/4
            g0  = x**2   - x**3
            diff = Gi0 - Gi1
c For stability
            if(g0.lt.1d-2) then
               g0=1d-2
            end if

            x1 = x + diff/g0
C limit x region
            if(x.lt.0.0) then
               x1 = 0.0
            else if(x.gt.1.0) then
               x1 = 1.0d0
            end if
            dx = x1 - x
            x = x1
 12      continue
c            Enue = ((aMmu**2 - aMe**2)/2/aMmu)*x
         Enue = (AMmu/2)*x
         tetnue = costet
         end

      subroutine pidecay(kpiin, ppi, nchild, kchild, pchild, spinmu,
     &     kinddecay)
      implicit real*8 (a-h,o-z)
c     control pi-decay:
c   input,
c     kpiin : kind of pion
c   outputs;
c     nchild : number of child particles
c     kchild(nchild) : kind of child particles
c     pchild(5,nchild) : 1-3 momentum, 4 energy, 5 mass of child particle
c     spinmu(5) : muon spin if there is muon in children
c      spin direction for mu+ and antispin direction for mu-
c      I.e. zeta represent the direction for which numu and e are 
c      emitted less energetically, when mu decays.
c     kinddecay: 1 mu + nu (pic), 2:  g + g (pi0)
      parameter (aMk0=0.497672d0, aMkc=0.493677d0, aMpi0=0.1349764d0, 
     &  aMpic=0.13956995d0,aMmu=0.105658389d0, aMe=0.00051099907d0,
     &  aMnu=0.d0, aMg=0.d0)
      parameter (kep=1, kem=2, kgamma=3,
     &     kkp=4, kkm=5, kk0l=6, kk0s=7,
     &     kneut=8, kneutbar=9, kpro=10, kprobar=11,
     &     kpip=12, kpim=13, kpi0=14,
     &     kmup=15, kmum=16, knue=17, knuebar=18, knumu=19,
     &     knumubar=20)

      integer*4 n, kchild(3) 
      real*8 pchild(5, 3), spinmu(5), p1(5), p2(5)
c                                                       decay at rest of pi
      if(kpiin.eq.kpip .or. kpiin.eq.kpim) then
c                                                       pic----> mu + neu
         kinddecay = 1
         nchild = 2
         call sampdecay2(aMpic, aMmu, aMnu, p2, p1)
         if(kpiin.eq.kpip) then
            kchild(1)=knumu
            kchild(2)=kmup
         else
            kchild(1)=knumubar
            kchild(2)=kmum
         end if
         pol = 1.0d0
         call samplespin(pol, p2, spinmu)
      else if(kpiin.eq.kpi0) then
c                                                       pi0----> g + g
         kinddecay = 2
         nchild = 2
         call sampdecay2(aMpi0, aMg, aMg, p1, p2)
         kchild(1)=kgamma
         kchild(2)=kgamma
      else
         write(0,*)'irregal kpiin in pidecay:',kpiin
         stop
      end if
c                                             save in child array
      do 1 i=1, nchild
         do 11 j=1, 5
            pchild(j,1) = p1(j)
            pchild(j,2) = p2(j)
 11      continue
 1    continue
      end

       subroutine cnbdcy(n, ecm, p, jw,  w, icon)
       implicit real*8 (a-h,o-z)
c      ***********************************************************
c
c   Originally coded by K.Kasahara for Cosmos (That's huge code!)
C
c        ref:  CPC.  40(1986)p359.  Kleiss, Stirling and Ellis
c
c       n: input.  number of ptcls >=2 (see however, for n=2,
c                  c2bdcy and for n=3, c3bdcy)
c     ecm: input.  cms energy.
c     p(5,n): input.  gives the mass of the i-th ptcl
c     p(1-4,n) : output undefined at input
c      jw: input.  0--->unweighted event (w=1) obtained. the event
c                  generated need not be discarded.
c                  1--> weighted event( w changes event to event )
c                  the event must be discarded according to the
c                  acceptance probability of w=weight/wax weight).
c      w: output.  see jw
c   icon: output.  0-->event generated successfully
c                  1-->ecm < sum of mass
c                  2-->could not generate (weight problem)
c
c
       integer n, jw, icon
       real*8 p(5,n)
       real*8 ecm, w
c      ------------------
       real*8 mu(1, 1)/0.d0/, wx, w0, wmax, gzai, u
c
       logical ok
       integer nc
c
	nc = 0    ! counter to break inf. loop
c       *** until loop*** 
       do while (.true.)
c            generate massless ptcls isotropically without conservation
           call cnbdc1(n, p)
c             conformal transformation to conserve 4-momentum
           call cnbdc2(n, ecm, p)
c$$$$$$$$$$$
c          call cnbdct(n, p)
c$$$$$$$$$$
c             get gzai to transform massive case
           call cnbdc3(n, ecm, p, mu, 0,  gzai, icon)
c          **********************
           if(icon .ne. 0) return
c          **********************
c             tranform to massive case
           call cnbdc4(n, p,  mu, 0, gzai)
c$$$$$$$$$$$
c          call cnbdct(n, p)
c$$$$$$$$$$
c             compute weight for massive  case
           call cnbdc5(n, ecm,p, wx)
c             compute weight for massless case
           call cnbdc6(n, ecm, w0)
c$$$$$$$$$$$$$$
c          write(*,*) ' wx=',wx,' w0=',w0
c$$$$$$$$$$$
           w=wx*w0
c                compute max possible weight
           call cnbdc7(n, ecm, p,  wmax)
           wmax=wmax*w0
           if(jw .eq. 0) then
c$$$$$$$$$$$$$$
c          write(*,*) ' wmax=',wmax
c$$$$$$$$$$$
c                judge if the event is to be accepted
              call rndc(u)
              if(wmax .eq. 0.d0) then
                 ok=.true.
              else
                 ok = u .lt. w/wmax
              endif
              w=1.
           else
              if(wmax .eq. 0.d0) then
                 w=1.d0
              else
                 w=w/wmax
              endif
              ok=.true.
           endif
           if(ok) goto 100
	   nc = nc +1
	   if(nc .gt. 20) then
	      icon = 2
	      goto 100
           endif
       enddo
  100  continue
       end

       subroutine cnbdct(n, p)
       implicit real*8 (a-h,o-z)
       integer n
c       record /ptcl/ p(n)
       real*8 p(5,n)
c
       real*8 sumx, sumy, sumz, sume
       integer i
c
           sumx=0.d0
           sumy=0.d0
           sumz=0.d0
           sume=0.d0
           do   i=1, n
              sumx=sumx+p(1,i)
              sumy=sumy+p(2,i)
              sumz=sumz+p(3,i)
              sume=sume+p(4,i)
           enddo
           write(*,*) ' sumx,y,z=',sumx, sumy, sumz, ' sume=',sume
       end

       subroutine cnbdc1(n, p)
       implicit real*8 (a-h, o-z)
c            generate massless ptcls isotropically without conservation
       integer n
c       record /ptcl/ p(n)
       real*8 p(5,n)
c
       integer i
       real*8 u1, u2, u3, x1,y1,z1
       do   i=1, n
c             *** until loop*** 
             do while (.true.)
                 call rndc(u1)
                 call rndc(u2)
                 u=u1*u2
                if(u .gt. 0.) goto 10
             enddo
   10        continue
             p(4,i) = -log(u)

             rr = 2.0
             do 20 while(rr.gt.1.d0.or.rr.lt.1d-2)
                call rndc(u1)
                call rndc(u2)
                call rndc(u3)

                x1 = 2*u1-1.d0
                y1 = 2*u2-1.d0
                z1 = 2*u3-1.d0
                rr = x1**2 + y1**2 + z1**2
 20          continue
             rr=sqrt(rr)
             p(1,i) = p(4,i)*x1/rr
             p(2,i) = p(4,i)*y1/rr
             p(3,i) = p(4,i)*z1/rr
           enddo
       end

       subroutine cnbdc2(n, ecm, p)
       implicit real*8 (a-h, o-z)
c             conformal transformation to conserve 4-momentum
       integer n
c       record /ptcl/ p(n)
       real*8 p(5,n)
       real*8 ecm
c
       real*8 sumx, sumy, sumz, sume, em, g
       real*8 a, x, bx, by, bz, bq, pe, tmp, px, py, pz
       integer i
c
       sumx=0.d0
       sumy=0.d0
       sumz=0.d0
       sume=0.d0
       do   i=1, n
          sumx=sumx+p(1,i)
          sumy=sumy+p(2,i)
          sumz=sumz+p(3,i)
          sume=sume+p(4,i)
       enddo
       em=sqrt( sume**2 - (sumx**2+sumy**2+sumz**2) )
       g=sume/em

       a=1.d0/(1.d0+g)
       x=ecm/em
       bx=-sumx/em
       by=-sumy/em
       bz=-sumz/em
c
       do   i=1, n
          bq=bx*p(1,i) + by*p(2,i) + bz*p(3,i)
          pe=x*(g*p(4,i) +bq)
          tmp=p(4,i)+a*bq
          px=x*(p(1,i) +   tmp*bx)
          py=x*(p(2,i) +   tmp*by)
          pz=x*(p(3,i) +   tmp*bz)
          p(1,i)=px
          p(2,i)=py
          p(3,i)=pz
          p(4,i)=pe
        enddo
       end
c      ***********************************************
       subroutine  cnbdc3(n, ecm, p, mu, inm, gzai, icon)
c             get gzai to transform massive case
c             put inm=0 if all mu are the same.
c      ***********************************************
       implicit real*8 (a-h,o-z)
       integer n, inm, icon
c       record /ptcl/ p(n)
       real*8 p(5,n)
       real*8 ecm, mu(inm, n), gzai
c
       real*8 eps/1.d-3/, f, fp, fow
       integer nr
c           initial guess of gzai
       gzai=.85d0
       nr=0
c          *** until loop*** 
       do while (.true.)
               call cnbdcf(n, ecm, p,  mu, inm,  gzai, f, fp)
               gzai=  gzai - f/fp
               fow=f
               nr=nr+1
c$$$$$$$$$$$$
c              write(*,*) ' fow=',fow
c$$$$$$$$$$$$
              if(abs(fow) .lt. eps .or. nr .gt. 15) goto 100
       enddo
  100  continue
       if(nr .gt. 15) then
           icon=1
       else
           icon=0
       endif
       end
c      *********************************************
       subroutine cnbdcf(n, ecm, p, 
     *            mu, inm, gzai, f, fp)
c      *********************************************
       implicit real*8 (a-h,o-z)
       integer n, inm
c       record /ptcl/ p(n)
       real*8 p(5,n)
       real*8 gzai, f, fp, mu(inm, n)
c
       real*8 mux, fx, tmp, ecm
       integer i
c
       fx=0.d0
       fp=0.d0
       do   i=1, n
c                  if compiler is good, we can use mu(1,i)
c                  even for inm=0; next is for safty.
          if(inm .eq. 0) then
c                 mux=mu(1,1)
             mux=0.d0
          else
             mux=mu(1,i)
          endif

          tmp=  sqrt(p(5,i)**2+
     *         gzai**2 *( p(4,i)**2 -mux**2 ) )
          fx=fx + tmp
          fp=fp + ( p(4,i)**2- mux**2)/ tmp
       enddo
       f=log(fx/ecm)
       fp=fp*gzai/fx
       end
c      *********************************************
       subroutine  cnbdc4(n, p, mu, inm,  gzai)
c             tranform to massive case
c      *********************************************
       implicit real*8 (a-h,o-z)
       integer n, inm
c       record /ptcl/ p(n)
       real*8 p(5,n)
       real*8 mu(inm, n), gzai

       real*8 mux
       integer i
       do   i=1,n
             p(1,i) = gzai*p(1,i)
             p(2,i) = gzai*p(2,i)
             p(3,i) = gzai*p(3,i)
c                next treatment is for safty
             if(inm .eq. 0) then
c                mux=mu(1,1)
                 mux=0.
             else
                 mux=mu(1,i)
             endif
             p(4,i) = sqrt(p(5,i)**2 +
     *       gzai**2*( p(4,i)**2-mux**2 ) )
           enddo
       end
       subroutine  cnbdc5(n, ecm, p, wx)
c            compute weig.p(4) for massive case
       implicit real*8 (a-h,o-z)
       integer n
c       record /ptcl/ p(n)
       real*8 p(5,n)
       real*8 ecm, wx
c
       real*8 sum1, pro2, sum3, pab
       integer i
c
          sum1=0.
          pro2=1.
          sum3=0.
          do   i=1, n
             pab =sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2)
             sum1=sum1+pab
             pro2=pro2* pab/p(4,i)
             sum3=sum3+pab**2/p(4,i)
          enddo
          wx=(sum1/ecm)**(2*n-3)*pro2 /sum3
       end
       subroutine cnbdc6(n, ecm, w0)
c             compute weig.p(4) for massless case
       implicit real*8 (a-h,o-z)
       integer n
       real*8 ecm, w0
c
       real*8 pi, hpi, gn1
       integer i
       parameter (pi=3.14159265d0, hpi=pi/2)
c
          gn1=1.
          do   i=1, n-2
             gn1=gn1*i
          enddo
          w0=  hpi**(n-1) * ecm**(2*n-4)/(n-1)/gn1/gn1
       end
       subroutine cnbdc7(n, ecm, p,  wmax)
c                compute max possible weig.p(4)
       implicit real*8 (a-h,o-z)
       integer n
       real*8 p(5,n)
       real*8 ecm, wmax

       integer idx(2), nm, i
       real*8 summ, en, beta
c          count massive ptcls
         nm=0
         do   i=1,n
            if(p(5,i) .gt. 0.d0) then
               nm=nm+1
               if(nm .le. 2) then
                  idx(nm)=i
               endif
            endif
         enddo
         if(nm .eq. 1) then
             wmax=(1. - p(5,idx(1))/ecm)**(2*n-3)
         elseif(nm .eq. 2)  then
             wmax=
     *        (1. + (p(5,idx(1))/ecm)**2 -
     *        (p(5,idx(2))/ecm)**2 )**2
     *       -4*(p(5,idx(1))/ecm)**2
             if(wmax .le. 0.d0)then
                wmax=1.d-30
             else
                wmax= wmax**(n-1.5d0)
             endif
         else
             summ=0.d0
             do   i=1, n
                 if(p(5,i) .gt. 0.d0) then
                    en=p(5,i)/ecm
                    summ=summ+en
                 endif
             enddo
             beta=1. - summ**2
             if(beta .le. 0.d0) then
                wmax=1.d-30
             else
                wmax=sqrt(beta)**(2*n+nm-5)
             endif
         endif
       end

      subroutine sampdecay2(aM0, am1, am2, p1, p2)
      implicit real*8 (a-h,o-z)
c   Sample momemtums for 2 body decay in mother particle system
c   am0,am1,am2: inputs
c   p1,p2 : 4+1 vector
c
      dimension p1(5), p2(5), u(3)

      e1 = 0.5d0*(Am0**2 + am1**2 - am2**2)/am0
      e2 = 0.5d0*(Am0**2 + am2**2 - am1**2)/am0

      call  rnduvec(u)

      ap1 = sqrt(e1**2 - am1**2)
      do 3 i=1, 3
         p1(i) =   ap1*u(i)
         p2(i) =  -ap1*u(i)
 3    continue
      p1(4) = e1
      p2(4) = e2

      p1(5) = am1
      P2(5) = am2
      end

c Energy valence
ckc3mu:
c  0.4936258  0.166662326  0.114956492  0.212006981
ck03mu:
c  0.497595354  0.166944719  0.115720306  0.214930328
ck03e:
c  0.497408122  0.131486635  0.131406588  0.2345149
ckc3e:
c  0.493411368  0.131190333  0.131110321  0.231110715
c
      subroutine rnduvec(vec)
      implicit real*8 (a-h,o-z)
      real*8 rvec, vec(3), u1, rr
c generage
      rr=2.0
      do 11 while(rr .gt. 1.d0 .or. rr .lt. 1.e-2)
         rr=0.
         do 12 i=1,3
            call rndc(u1)
            vec(i) = 2*u1 - 1.d0
            rr=rr+vec(i)**2
 12      continue
 11   continue

      rr=sqrt(rr)
      do 13 i=1,3
         vec(i) = vec(i)/rr
 13   continue
      end

      subroutine k03e(pl, pnu, ppi, pol)
      implicit real*8 (a-h,o-z)
      parameter (aMk0=0.497672d0, aMkc=0.493677d0, aMpi0=0.1349764d0,
     &  aMpic=0.13956995d0,aMmu = 0.105658389d0, aMe=0.00051099907d0,
     &  aMnu=0.d0)
      real*8 pl(5), pnu(5), ppi(5), u(3)

      El  = Samplek03e()
      apl = sqrt(el**2 - aMe**2)
      call rnduvec(u)
      do 1 i=1, 3
         pl(i)=apl*u(i)
 1    continue
      pl(4)= El
      pl(5)= aMe

      Enu = Samplek03enu()
      apnu = enu
      call rnduvec(u)
      do 2 i=1, 3
         pnu(i)=apnu*u(i)
 2    continue
      pnu(4)= Enu
      pnu(5)= aMnu

      Epi = Samplek03epi()
      appi = sqrt(epi**2 - aMpic**2)
      call rnduvec(u)
      do 3 i=1, 3
         ppi(i)=appi*u(i)
 3    continue
      ppi(4)= Epi
      ppi(5)= aMpic

      pol=0
      end

      subroutine kc3e(pl, pnu, ppi, pol)
      implicit real*8 (a-h,o-z)
      parameter (aMk0=0.497672d0, aMkc=0.493677d0, aMpi0=0.1349764d0,
     &  aMpic=0.13956995d0,aMmu = 0.105658389d0, aMe=0.00051099907d0,
     &  aMnu=0.d0)
      real*8 pl(5), pnu(5), ppi(5), u(3)

      El  = Samplekc3e()
      apl = sqrt(el**2 - aMe**2)
      call rnduvec(u)
      do 1 i=1, 3
         pl(i)=apl*u(i)
 1    continue
      pl(4)= El
      pl(5)= aMe

      Enu = Samplekc3enu()
      apnu = sqrt(enu**2 - aMnu**2)
      call rnduvec(u)
      do 2 i=1, 3
         pnu(i)=apnu*u(i)
 2    continue
      pnu(4)= Enu
      pnu(5)= aMnu

      Epi = Samplekc3epi()
      appi = sqrt(epi**2 - aMpi0**2)
      call rnduvec(u)
      do 3 i=1, 3
         ppi(i)=appi*u(i)
 3    continue
      ppi(4)= Epi
      ppi(5)= aMpi0

      pol=0
      end

      subroutine k03mu(pl, pnu, ppi, pol)
      implicit real*8 (a-h,o-z)
      parameter (aMk0=0.497672d0, aMkc=0.493677d0, aMpi0=0.1349764d0,
     &  aMpic=0.13956995d0,aMmu = 0.105658389d0, aMe=0.00051099907d0,
     &  aMnu=0.d0)
      real*8 pl(5), pnu(5), ppi(5), u(3)
      real*8 k03mupol

      El  = Samplek03mu()
      apl = sqrt(el**2 - aMmu**2)
      call rnduvec(u)
      do 1 i=1, 3
         pl(i)=apl*u(i)
 1    continue
      pl(4)= El
      pl(5)= aMmu

      Enu = Samplek03munu()
      apnu = sqrt(enu**2 - aMnu**2)
      call rnduvec(u)
      do 2 i=1, 3
         pnu(i)=apnu*u(i)
 2    continue
      pnu(4)= Enu
      pnu(5)= aMnu

      Epi = Samplek03mupi()
      appi = sqrt(epi**2 - aMpic**2)
      call rnduvec(u)
      do 3 i=1, 3
         ppi(i)=appi*u(i)
 3    continue
      ppi(4)= Epi
      ppi(5)= aMpic
      pol = k03mupol(El)
      end

      subroutine kc3mu(pl, pnu, ppi, pol)
      implicit real*8 (a-h,o-z)
      parameter (aMk0=0.497672d0, aMkc=0.493677d0, aMpi0=0.1349764d0,
     &  aMpic=0.13956995d0,aMmu = 0.105658389d0, aMe=0.00051099907d0,
     &  aMnu=0.d0)
      real*8 pl(5), pnu(5), ppi(5), u(3)
      real*8 kc3mupol

      El  = Samplekc3mu()
      apl = sqrt(el**2 - aMmu**2)
      call rnduvec(u)
      do 1 i=1, 3
         pl(i)=apl*u(i)
 1    continue
      pl(4)= El
      pl(5)= aMmu

      Enu = Samplekc3munu()
      apnu = sqrt(enu**2 - aMnu**2)
      call rnduvec(u)
      do 2 i=1, 3
         pnu(i)=apnu*u(i)
 2    continue
      pnu(4)= Enu
      pnu(5)= aMnu

      Epi = Samplekc3mupi()
      appi = sqrt(epi**2 - aMpi0**2)
      call rnduvec(u)
      do 3 i=1, 3
         ppi(i)=appi*u(i)
 3    continue
      ppi(4)= Epi
      ppi(5)= aMpi0
      pol = kc3mupol(El)
      end

c K0l -> e+  +  nu_e + pi-  &  K0l -> e- +  nu_e~ + pi+
      function Samplek03e()
      real*8 Samplek03e
      real*8 k03ee(101),u, x
      integer i
       save
c k03e-e
c    set Minimum etot = Me + eps (allways increment last digit)
      data (k03ee(i),i=1, 50)/
     &.000511,.030122,.039491,.044693,.050061,.053949,.057807,.061065,
     &.064192,.067025,.069725,.072252,.074667,.076965,.079173,.081295,
     &.083344,.085326,.087248,.089117,.090936,.092711,.094444,.096140,
     &.097801,.099430,.101029,.102601,.104147,.105669,.107169,.108649,
     &.110109,.111550,.112975,.114385,.115779,.117160,.118527,.119883,
     &.121227,.122561,.123886,.125201,.126508,.127807,.129100,.130385,
     &.131665,.132940/
c k03e-e
      data (k03ee(i),i=101, 51, -1)/
     &.229265,.214314,.209088,.205773,.202519,.199937,.197429,.195204,
     &.193068,.191076,.189162,.187335,.185571,.183868,.182213,.180602,
     &.179031,.177494,.175988,.174510,.173056,.171626,.170215,.168823,
     &.167448,.166088,.164742,.163408,.162086,.160774,.159471,.158177,
     &.156890,.155609,.154334,.153065,.151800,.150538,.149280,.148024,
     &.146770,.145517,.144265,.143014,.141761,.140508,.139254,.137997,
     &.136738,.135476,.134210/
      call rndc(u)
      x = 100*u
      i = int(x)
      x = x - real(i)
      i = i + 1
      Samplek03e = k03ee(i)*(1.0d0 -x) + k03ee(i+1)*x
      end

      function Samplek03enu()
      real*8 Samplek03enu
      real*8 k03enu(101),u, x
      integer i
       save
c k03e-nu
c    set Minimum etot = Mnu + eps (always increment last digit)
      data (k03enu(i),i=1, 50)/
     &.000001,.030091,.039512,.044674,.050067,.053941,.057808,.061061,
     &.064190,.067023,.069723,.072250,.074665,.076963,.079171,.081293,
     &.083342,.085324,.087247,.089115,.090934,.092709,.094443,.096138,
     &.097800,.099429,.101028,.102600,.104146,.105668,.107168,.108647,
     &.110107,.111549,.112974,.114383,.115778,.117158,.118526,.119882,
     &.121226,.122560,.123885,.125200,.126507,.127806,.129098,.130384,
     &.131664,.132939/
c k03e-nu
      data (k03enu(i),i=101, 51, -1)/
     &.229265,.214314,.209087,.205772,.202518,.199936,.197429,.195203,
     &.193067,.191075,.189161,.187335,.185571,.183867,.182212,.180602,
     &.179030,.177493,.175987,.174509,.173056,.171625,.170214,.168823,
     &.167447,.166087,.164741,.163408,.162085,.160773,.159470,.158176,
     &.156889,.155608,.154333,.153064,.151799,.150537,.149279,.148023,
     &.146769,.145516,.144264,.143013,.141760,.140507,.139253,.137996,
     &.136737,.135475,.134209/
      call rndc(u)
      x = 100*u
      i = int(x)
      x = x - real(i)
      i = i + 1
      Samplek03enu = k03enu(i)*(1.0d0 -x) + k03enu(i+1)*x
      end

      function Samplek03epi()
      real*8 Samplek03epi
      real*8 k03epi(101),u,x
      integer i
       save
c k03e-pi
c    set Minimum etot = Mpic + eps (always increment last digit)
      data (k03epi(i),i=1, 50)/
     &.139570,.162673,.170460,.175057,.179502,.182821,.186006,.188710,
     &.191262,.193567,.195737,.197755,.199665,.201469,.203187,.204825,
     &.206393,.207898,.209346,.210741,.212089,.213393,.214656,.215882,
     &.217073,.218231,.219359,.220459,.221532,.222579,.223603,.224604,
     &.225584,.226544,.227484,.228406,.229311,.230198,.231070,.231927,
     &.232769,.233597,.234411,.235213,.236002,.236779,.237544,.238299,
     &.239042,.239776/
c k03e-pi
      data (k03epi(i),i=101, 51, -1)/
     &.268407,.267929,.267476,.267020,.266561,.266099,.265633,.265164,
     &.264691,.264215,.263735,.263252,.262765,.262274,.261779,.261280,
     &.260777,.260270,.259759,.259244,.258724,.258200,.257671,.257137,
     &.256599,.256056,.255508,.254955,.254397,.253833,.253264,.252690,
     &.252109,.251523,.250931,.250333,.249728,.249118,.248500,.247876,
     &.247244,.246606,.245960,.245306,.244645,.243976,.243298,.242612,
     &.241917,.241213,.240499/
      call rndc(u)
      x = 100*u
      i = int(x)
      x = x - real(i)
      i = i + 1
      Samplek03epi = k03epi(i)*(1.0d0 -x) + k03epi(i+1)*x
      end

c K0l -> mu+  +  nu_mu + pi-  &  K0l -> mu- +  nu_mu~ + pi+
      function Samplek03mu()
      real*8 Samplek03mu
      real*8 k03mumu(101),u,x
      integer i
       save
c k03mu-mu
c    set Minimum etot = Mmu + eps (always increment last digit)
      data (k03mumu(i),i=1, 50)/
     &.105659,.110121,.112741,.114895,.116837,.118602,.120257,.121815,
     &.123301,.124725,.126096,.127422,.128709,.129960,.131181,.132374,
     &.133542,.134687,.135812,.136917,.138005,.139077,.140135,.141179,
     &.142210,.143229,.144238,.145236,.146225,.147206,.148178,.149143,
     &.150100,.151052,.151997,.152937,.153871,.154801,.155727,.156649,
     &.157567,.158482,.159394,.160304,.161212,.162117,.163021,.163924,
     &.164826,.165728/
c k03mu-mu
      data (k03mumu(i),i=101, 51, -1)/
     &.240481,.228132,.223870,.221210,.218585,.216526,.214523,.212756,
     &.211062,.209488,.207979,.206544,.205160,.203828,.202536,.201283,
     &.200062,.198871,.197707,.196567,.195448,.194349,.193269,.192205,
     &.191156,.190121,.189099,.188089,.187090,.186101,.185121,.184150,
     &.183186,.182230,.181281,.180337,.179400,.178467,.177539,.176616,
     &.175696,.174780,.173867,.172956,.172048,.171142,.170237,.169334,
     &.168432,.167530,.166629/
      call rndc(u)
      x = 100*u
      i = int(x)
      x = x - real(i)
      i = i + 1
      Samplek03mu = k03mumu(i)*(1.0d0 -x) + k03mumu(i+1)*x
      end

      function Samplek03munu()
      real*8 Samplek03munu
      real*8 k03munu(101),u,x
      integer i
       save
c k03mu-nu
c    set Minimum etot = Mnu + eps (allways increment last digit)
      data (k03munu(i),i=1, 50)/
     &.000001,.027984,.036752,.041562,.046585,.050197,.053801,.056835,
     &.059754,.062397,.064916,.067276,.069530,.071677,.073739,.075722,
     &.077637,.079489,.081287,.083034,.084736,.086396,.088018,.089606,
     &.091161,.092686,.094184,.095656,.097105,.098532,.099938,.101325,
     &.102695,.104047,.105384,.106707,.108016,.109313,.110597,.111871,
     &.113134,.114388,.115634,.116871,.118100,.119323,.120540,.121750,
     &.122956,.124157/
c k03mu-nu
      data (k03munu(i),i=101, 51, -1)/
     &.218049,.202530,.197245,.193985,.190748,.188222,.185758,.183588,
     &.181505,.179570,.177713,.175945,.174240,.172596,.171001,.169451,
     &.167940,.166464,.165019,.163603,.162211,.160842,.159494,.158165,
     &.156852,.155555,.154272,.153002,.151743,.150495,.149256,.148026,
     &.146804,.145588,.144379,.143175,.141976,.140781,.139589,.138400,
     &.137214,.136029,.134845,.133662,.132479,.131296,.130111,.128926,
     &.127738,.126547,.125354/
      call rndc(u)
      x = 100*u
      i = int(x)
      x = x - real(i)
      i = i + 1
c      Samplek03munu = k03munu(i)*(1.0d0 -x) + k03munu(i+1)*x
c This factor 0.941 is to retain the energy conservation 
      Samplek03munu = (k03munu(i)*(1.0d0 -x) + k03munu(i+1)*x)*0.941d0
      end

      function Samplek03mupi()
      real*8 Samplek03mupi
      real*8 k03mupi(101),u,x
      integer i
       save
c k03mu-pi
c    set Minimum etot = Mpic + eps (allways increment last digit)
      data (k03mupi(i),i=1, 50)/
     &.139570,.150435,.158561,.160463,.165180,.167321,.170341,.172439,
     &.174722,.176639,.178541,.180271,.181940,.183509,.185015,.186451,
     &.187833,.189162,.190444,.191683,.192882,.194046,.195176,.196275,
     &.197346,.198390,.199409,.200404,.201377,.202330,.203263,.204178,
     &.205076,.205957,.206822,.207673,.208509,.209333,.210143,.210941,
     &.211728,.212504,.213269,.214024,.214770,.215506,.216234,.216953,
     &.217665,.218368/
c k03mu-pi
      data (k03mupi(i),i=101, 51, -1)/
     &.257191,.252920,.251334,.250247,.249192,.248307,.247446,.246656,
     &.245890,.245159,.244448,.243757,.243080,.242416,.241764,.241120,
     &.240484,.239854,.239230,.238610,.237993,.237380,.236769,.236160,
     &.235552,.234945,.234339,.233732,.233125,.232517,.231908,.231298,
     &.230686,.230072,.229455,.228836,.228215,.227590,.226962,.226330,
     &.225694,.225054,.224410,.223761,.223107,.222448,.221783,.221113,
     &.220437,.219754,.219065/
      call rndc(u)
      x = 100*u
      i = int(x)
      x = x - real(i)
      i = i + 1
      Samplek03mupi = k03mupi(i)*(1.0d0 -x) + k03mupi(i+1)*x
      end

      function k03mupol(Emu)
      real*8 Emu, k03mupol, emumin, emumax, k0mupol(101), x
      integer*4 i
      save
      data emumin, emumax /.105658389,.240481023/
      data (k0mupol(i),i=1,101)/
     &.000000,.158923,.222525,.269875,.308624,.341774,.370886,.396898,
     &.420429,.441914,.461674,.479954,.496946,.512804,.527655,.541603,
     &.554738,.567134,.578857,.589963,.600502,.610516,.620045,.629123,
     &.637781,.646046,.653943,.661495,.668722,.675644,.682277,.688637,
     &.694739,.700595,.706218,.711620,.716810,.721799,.726595,.731207,
     &.735643,.739909,.744013,.747960,.751757,.755409,.758920,.762296,
     &.765540,.768658,.771651,.774525,.777281,.779924,.782455,.784876,
     &.787191,.789400,.791506,.793510,.795413,.797216,.798920,.800526,
     &.802033,.803441,.804751,.805961,.807072,.808081,.808988,.809790,
     &.810485,.811071,.811545,.811902,.812140,.812252,.812234,.812080,
     &.811781,.811331,.810720,.809938,.808972,.807809,.806433,.804826,
     &.802967,.800831,.798392,.795616,.792465,.788893,.784846,.780261,
     &.775058,.769145,.762405,.754694,.746145/
      if(Emu.lt.emumin .or. Emu.gt.emumax) then
         k03mupol=0.0
         return
      end if
      x = 100.d0*(Emu - emumin)/(emumax - emumin)
      i = int(x)
      x = x - real(i)
      i = i + 1
      k03mupol=k0mupol(i)*(1.0d0-x) + k0mupol(i+1)*x
      end

c K+ -> e+ +  nu_e + pi0 & K- -> e- +  nu_e~ + pi0
      function Samplekc3e()
      real*8 Samplekc3e
      real*8 kc3ee(101), u, x
      integer i
       save
c kc3e-e
c    set Minimum etot = Me + eps (always increment last digit)
      data (kc3ee(i),i=1, 50)/
     &.000511,.030071,.039423,.044615,.049973,.053854,.057705,.060956,
     &.064077,.066905,.069599,.072121,.074531,.076825,.079028,.081145,
     &.083190,.085167,.087086,.088950,.090765,.092536,.094266,.095958,
     &.097615,.099240,.100836,.102404,.103946,.105465,.106962,.108437,
     &.109894,.111332,.112753,.114159,.115550,.116927,.118291,.119643,
     &.120984,.122315,.123636,.124948,.126251,.127547,.128835,.130118,
     &.131394,.132665/
c kc3e-e
      data (kc3ee(i),i=101, 51, -1)/
     &.228387,.213635,.208462,.205170,.201944,.199378,.196888,.194677,
     &.192554,.190573,.188670,.186854,.185099,.183404,.181757,.180155,
     &.178591,.177061,.175562,.174090,.172643,.171218,.169814,.168428,
     &.167058,.165704,.164363,.163034,.161717,.160409,.159111,.157822,
     &.156539,.155263,.153993,.152728,.151467,.150210,.148955,.147704,
     &.146454,.145205,.143957,.142709,.141461,.140211,.138960,.137708,
     &.136452,.135194,.133931/
      call rndc(u)
      x = 100*u
      i = int(x)
      x = x - real(i)
      i = i + 1
      Samplekc3e = kc3ee(i)*(1.0d0 -x) + kc3ee(i+1)*x
      end

      function Samplekc3enu()
      real*8 Samplekc3enu
      real*8 kc3enu(101),u, x
      integer i
       save
c kc3e-nu
c    set Minimum etot = Mpi + eps (always increment last digit)
      data (kc3enu(i),i=1, 50)/
     &.000001,.030039,.039444,.044596,.049979,.053846,.057705,.060952,
     &.064075,.066902,.069597,.072119,.074529,.076823,.079026,.081143,
     &.083188,.085166,.087084,.088948,.090764,.092534,.094264,.095956,
     &.097614,.099239,.100834,.102403,.103945,.105464,.106960,.108436,
     &.109892,.111331,.112752,.114158,.115549,.116926,.118290,.119642,
     &.120983,.122314,.123634,.124946,.126250,.127545,.128834,.130116,
     &.131393,.132664/
c kc3e-nu
      data (kc3enu(i),i=101, 51, -1)/
     &.228386,.213634,.208461,.205169,.201943,.199377,.196888,.194676,
     &.192553,.190572,.188669,.186853,.185098,.183403,.181756,.180154,
     &.178590,.177060,.175561,.174089,.172642,.171218,.169813,.168427,
     &.167057,.165703,.164362,.163033,.161716,.160409,.159110,.157821,
     &.156538,.155262,.153992,.152727,.151466,.150208,.148954,.147702,
     &.146453,.145204,.143956,.142708,.141459,.140210,.138959,.137706,
     &.136451,.135192,.133930/
      call rndc(u)
      x = 100*u
      i = int(x)
      x = x - real(i)
      i = i + 1
      Samplekc3enu = kc3enu(i)*(1.0d0 -x) + kc3enu(i+1)*x
      end

      function Samplekc3epi()
      real*8 Samplekc3epi
      real*8 kc3epi(101),u, x
      integer i
       save
c kc3e-pi
c     Minimum etot = Mpi0 + eps (always increment last digit)
      data (kc3epi(i),i=1, 50)/
     &.134977,.158463,.166367,.171024,.175531,.178891,.182117,.184855,
     &.187438,.189770,.191965,.194007,.195939,.197763,.199500,.201155,
     &.202740,.204261,.205724,.207134,.208495,.209813,.211089,.212327,
     &.213530,.214700,.215839,.216949,.218032,.219090,.220123,.221134,
     &.222123,.223091,.224040,.224971,.225883,.226779,.227659,.228523,
     &.229372,.230208,.231029,.231838,.232634,.233418,.234190,.234950,
     &.235700,.236440/
c kc3e-pi
      data (kc3epi(i),i=101, 51, -1)/
     &.265290,.264809,.264353,.263894,.263432,.262966,.262497,.262024,
     &.261548,.261069,.260586,.260099,.259608,.259114,.258615,.258113,
     &.257606,.257095,.256581,.256061,.255538,.255010,.254477,.253939,
     &.253397,.252850,.252298,.251741,.251178,.250610,.250037,.249458,
     &.248873,.248283,.247686,.247083,.246474,.245858,.245235,.244606,
     &.243970,.243326,.242675,.242016,.241349,.240675,.239991,.239299,
     &.238599,.237889,.237169/
      call rndc(u)
      x = 100*u
      i = int(x)
      x = x - real(i)
      i = i + 1
      Samplekc3epi = kc3epi(i)*(1.0d0 -x) + kc3epi(i+1)*x
      end

c------------ K+ -> mu+ +  nu_mu + pi0 & K- -> mu- +  nu_mu~ + pi0
      function Samplekc3mu()
      real*8 Samplekc3mu
      real*8 kc3mumu(101),u, x
      integer i
       save
c kc3mu-mu
c    set Minimum etot = Mmu + eps (always increment last digit)
      data (kc3mumu(i),i=1, 50)/
     &.105659,.110135,.112759,.114915,.116856,.118620,.120274,.121829,
     &.123313,.124733,.126100,.127422,.128704,.129952,.131168,.132356,
     &.133518,.134658,.135777,.136877,.137960,.139026,.140077,.141115,
     &.142140,.143153,.144156,.145148,.146130,.147104,.148070,.149028,
     &.149979,.150924,.151862,.152795,.153723,.154646,.155564,.156479,
     &.157390,.158298,.159203,.160105,.161005,.161904,.162800,.163696,
     &.164590,.165484/
c kc3mu-mu
      data (kc3mumu(i),i=101, 51, -1)/
     &.239693,.227371,.223128,.220488,.217879,.215836,.213847,.212094,
     &.210414,.208853,.207357,.205933,.204562,.203241,.201961,.200718,
     &.199508,.198328,.197174,.196044,.194936,.193848,.192777,.191723,
     &.190684,.189658,.188646,.187645,.186655,.185675,.184704,.183742,
     &.182787,.181840,.180899,.179964,.179035,.178111,.177192,.176276,
     &.175365,.174457,.173552,.172650,.171749,.170851,.169955,.169059,
     &.168165,.167271,.166378/
      call rndc(u)
      x = 100*u
      i = int(x)
      x = x - real(i)
      i = i + 1
      Samplekc3mu = kc3mumu(i)*(1.0d0 -x) + kc3mumu(i+1)*x
      end

      function Samplekc3munu()
      real*8 Samplekc3munu
      real*8 kc3munu(101),u, x
      integer i
       save
c kc3mu-nu
c    set Minimum etot = Mnu + eps (always increment last digit)
      data (kc3munu(i),i=1, 50)/
     &.000001,.027901,.036642,.041437,.046444,.050045,.053637,.056662,
     &.059571,.062205,.064717,.067068,.069316,.071455,.073511,.075487,
     &.077395,.079242,.081033,.082775,.084471,.086125,.087742,.089324,
     &.090873,.092393,.093886,.095353,.096797,.098219,.099620,.101002,
     &.102367,.103714,.105047,.106364,.107669,.108960,.110240,.111509,
     &.112768,.114018,.115258,.116491,.117715,.118934,.120145,.121351,
     &.122552,.123749/
c kc3mu-nu
      data (kc3munu(i),i=101, 51, -1)/
     &.217080,.201708,.196466,.193227,.190013,.187503,.185055,.182899,
     &.180828,.178904,.177058,.175300,.173604,.171969,.170382,.168841,
     &.167338,.165869,.164432,.163022,.161638,.160276,.158934,.157611,
     &.156305,.155014,.153736,.152472,.151219,.149977,.148743,.147519,
     &.146302,.145092,.143888,.142689,.141495,.140305,.139118,.137935,
     &.136753,.135573,.134394,.133216,.132038,.130859,.129680,.128498,
     &.127315,.126130,.124941/
      call rndc(u)
      x = 100*u
      i = int(x)
      x = x - real(i)
      i = i + 1
c      Samplekc3munu = kc3munu(i)*(1.0d0 -x) + kc3munu(i+1)*x
c This factor 0.938 is to retain the energy conservation 
      Samplekc3munu = (kc3munu(i)*(1.0d0 -x) + kc3munu(i+1)*x)*0.938d0
      end

      function Samplekc3mupi()
      real*8 Samplekc3mupi
      real*8 kc3mupi(101),u, x
      integer i
       save
c kc3mu-pi
c    set Minimum etot = Mpi + eps (always increment last digit)
      data (kc3mupi(i),i=1, 50)/
     &.134977,.148401,.154255,.158346,.161862,.164775,.167415,.169764,
     &.171934,.173932,.175803,.177559,.179220,.180796,.182299,.183737,
     &.185116,.186443,.187722,.188957,.190154,.191313,.192439,.193534,
     &.194600,.195639,.196653,.197643,.198611,.199558,.200486,.201395,
     &.202287,.203163,.204022,.204867,.205698,.206516,.207320,.208113,
     &.208894,.209664,.210423,.211172,.211912,.212643,.213365,.214078,
     &.214784,.215481/
c kc3mu-pi
      data (kc3mupi(i),i=101, 51, -1)/
     &.253984,.249734,.248156,.247079,.246030,.245153,.244299,.243515,
     &.242756,.242031,.241325,.240640,.239970,.239312,.238665,.238027,
     &.237397,.236773,.236154,.235540,.234929,.234322,.233717,.233113,
     &.232511,.231910,.231309,.230708,.230106,.229504,.228901,.228296,
     &.227690,.227081,.226471,.225857,.225241,.224622,.224000,.223374,
     &.222744,.222109,.221471,.220828,.220180,.219526,.218868,.218203,
     &.217532,.216856,.216172/
      call rndc(u)
      x = 100*u
      i = int(x)
      x = x - real(i)
      i = i + 1
      Samplekc3mupi = kc3mupi(i)*(1.0d0 -x) + kc3mupi(i+1)*x
      end

      function kc3mupol(Emu)
      real*8 Emu, kc3mupol, emumin, emumax, kcmupol(101), x
      integer*4 i
      save
      data emumin, emumax /.105658389,.239693207/
      data (kcmupol(i),i=1,101)/
     &.000000,.172241,.241026,.292143,.333901,.369568,.400839,.428739,
     &.453941,.476921,.498028,.517531,.535638,.552520,.568313,.583134,
     &.597078,.610229,.622658,.634426,.645588,.656191,.666278,.675885,
     &.685048,.693796,.702157,.710155,.717813,.725152,.732190,.738945,
     &.745433,.751669,.757665,.763435,.768990,.774341,.779497,.784470,
     &.789266,.793894,.798362,.802676,.806844,.810873,.814766,.818531,
     &.822172,.825695,.829103,.832401,.835594,.838684,.841676,.844573,
     &.847378,.850095,.852725,.855271,.857737,.860123,.862434,.864669,
     &.866832,.868924,.870946,.872900,.874787,.876608,.878364,.880056,
     &.881685,.883251,.884754,.886194,.887572,.888887,.890138,.891326,
     &.892447,.893502,.894488,.895404,.896246,.897011,.897695,.898294,
     &.898802,.899213,.899518,.899707,.899771,.899693,.899458,.899045,
     &.898428,.897574,.896443,.894983,.893213/
      if(Emu.lt.emumin .or. Emu.gt.emumax) then
         kc3mupol=0.0
         return
      end if
      x = 100.d0*(Emu - emumin)/(emumax - emumin)
      i = int(x)
      x = x - real(i)
      i = i + 1
      kc3mupol=kcmupol(i)*(1.0d0-x) + kcmupol(i+1)*x
      end

