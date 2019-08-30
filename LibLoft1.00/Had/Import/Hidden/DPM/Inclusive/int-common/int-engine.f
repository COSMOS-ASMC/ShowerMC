c      parameter (kep=1,kem=2,kgamma=3,
c     &     kkp=4,kkm=5,kk0l=6,kk0s=7,
c     &     kneut=8,kneutbar=9,kpro=10,kprobar=11,
c     &     kpip=12,kpim=13,kpi0=14,
c     &     kmup=15, kmum=16, knue=17, knuebar=18, knumu=19,
c     &     knumubar=20)
c
c      dimension cptl(20)
c      data cptl/'E+  ', 'E-  ','GAMM','K+  ','K-  ','K0L ','K0S ',
c     &'N+  ','N-  ','P+  ','P-  ','PI+ ','PI- ','PI0 ','MU+ ',
c     &'MU- ','NUE ','NUEB','NUM ','NUMB'/
c
      subroutine hadronint(kin, pin, nchild, kchild, pchild, nextc)
      implicit doubleprecision (a-h,o-z)
      parameter (maxchild=1000)

c  Ratio of production to total-inelastic corss-section 
c  For quasi-elastic correction
      common /prod_ratio/ prodratio

      dimension pin(5), pchild(5,maxchild), nextc(maxchild)
      dimension  dirp(3), eko(10)
      integer kin, kchild(maxchild), kout
      integer toin(20), toout(10)
      logical quasi_elastic

      include '../include/sim-control.inc' ! next action, boundary
      include '../include/atmnc-particle-code.inc' 
      include '../include/atmnc-particle-mass2.inc' 
c
c  iin, iout 内部での粒子コード
c  kin, kout 外部での粒子コード
c  toin : 外部から内部へのコード変更(ハドロン相互作用しない粒子は0)
c                 1  2  3  4  5  6  7  8  9 10 11 12 13 14
      data toin / 0, 0, 0, 3, 4, 5, 5, 8, 9, 6, 7, 1, 2, 0,
c                15 16 17 18 19 20
     &            0, 0, 0, 0, 0, 0/
c  toout : 内部から外部へのコード変更(γ、レプトンは出力されない)
c         修正 内部での K0l と K0s の区別を出力でもやめる。 (1999.3.25)
c                 1  2  3  4  5  6  7  8  9 10
      data toout/14,12,13, 4, 5, 0,10,11, 8, 9/

      iin = toin(kin)
      ekin = pin(4) - pin(5)
      nchild=0

      if(iin.eq.0 .or. ekin.lt.0.2d0) then
C ハドロン相互作用しない粒子、又はエネルギーが低すぎる
c         write(0,*)
c     &   ' not supported ptl :', kin,' Child = mother ptl'
         nchild=nchild + 1
         do 231 k=1,5
            pchild(k, nchild) = pin(k)
 231     continue
         kchild(nchild) = kin
         return
      end if

c 4+1元運動量から、運動の方向ベクトル、力学エネルギーを求める。
      if(Ekin.lt.0.2d0) then
         e00 = 0.2d0
      elseif(Ekin.gt.1.d6) then
         e00 = 1.d6
      else
         e00 = ekin
      end if
      el = log10(e00)

      rr=sqrt(pin(1)**2 + pin(2)**2 + pin(3)**2)
      do 1 i=1, 3
         dirp(i) = pin(i)/rr
 1    continue

c Spallation and Quasi-elastic energyloss prosess,
      call rndc(u)
      quasi_elastic = (u.gt.prodratio)
c      write(*,*) ekin, u, quasi_elastic

      if( kin.gt.7 .and. kin.lt.12.and.e00.gt.3.1622) then ! nucleon
c         write(*,*) kin, ekin
c         r_spa = 0.6 ! step @3

         r_spa = 2*log10(el + 1.d0) - 0.35218d0 ! a1
         r_spa = 2*r_spa ! a2
c         write(*,*) ekin, el, r_spa
         nspallation = int(r_spa)
         d_spa = r_spa -real(nspallation)

         call rndc(u)
         if(u.lt.d_spa) then
            nspallation = nspallation + 1
         end if

c         write(*,*) nchild0
         do nc = 1, nspallation
            nchild = nchild + 1

            call rndc(u)
            if(u.gt.0.5d0) then
               kout = 10 ! proton
            else
               kout = 8 ! neutron
            end if

            Ekspa = 1e10
            do while (Ekspa.ge.e00)

               call rndc(u)
               Ekspal = 1.2*u - 1.
c               Ekspal = 0.7*u - .5

c     Ekspal = (2.2*u - 1.84056)/0.622148
               Ekspa = 10.D0**Ekspal
c               if(e00 .lt. 3.78216D0) then
c                  Ekspa = e00*(Ekspa/3.78216D0)
c               end if
            enddo
c            e00 = e00 - Ekspa

            echild = Ekspa + am(kout)
c            write(*,*)ekspal, ekspa

            gamma = echild/am(kout)
            beta = sqrt(1.d0 - 1.d0/gamma**2)
            p_abs = gamma*beta*am(kout)

            amass = am(kout)
            call sample_intcospt(kin, amass, kout, Ekspa, p_abs, cospt)
            call sample_angvec(cospt, pin, dirp)
c            write(*,'(3I3,5F15.5)') kin,kout,nchild,ekspa, cospt

            do 111 j=1, 3
               pchild(j,nchild) = p_abs*dirp(j)
 111        continue
            pchild(4,nchild) = echild
            pchild(5,nchild) = am(kout)
            kchild(nchild) = kout
            nextc(nchild) = iact_yet ! next action unknown
         enddo
      end if

! quasi-elastic interaction
      if(quasi_elastic) then
         call rndc(u)
         eko(1) = e00*u
         if(e00.gt.2.) then
            eko(2) = e00*(1.d0 -u)
         else
            eko(2) = 0.
         end if

c         write(*,*) 'quasi-elastic', ekin
         do 1002 nc=1, 2
            ekout = eko(nc)
c            if(ekout.lt.0.5d0) go to 1002

            nchild = nchild + 1
c  charge exchange and punched out particle
            call rndc(u)
            if(nc.eq.1) then    ! charge exchange
c               if(3*u.gt.1.) then
                  kout = kin
c               else 
c                  if(kin.eq.kpro) then
c                     kout = kneut
c                  else if(kin.eq.kneut) then
c                     kout = kpro
c                  else if(kin.eq.kprobar)
c                     kout = kneutbar
c                  else
c                     kout = kprobar
c                  end if
c               end if
            else                ! punched out
               if(2*u.gt.1) then
                  kout = kpro
               else
                  kout = kneut
               end if
            end if
c elastic (box) energyloss and compensating particle
            echild = ekout + am(kout)

            gamma = echild/am(kout)
            beta = sqrt(1.d0 - 1.d0/gamma**2)
            p_abs = gamma*beta*am(kout)

            amass = am(kout)
            call sample_intcospt(kin, amass, kout, ekout, p_abs, cospt)
            call sample_angvec(cospt, pin, dirp)

            do 112 j=1, 3
               pchild(j,nchild) = p_abs*dirp(j)
 112        continue
            pchild(4,nchild) = echild
            pchild(5,nchild) = am(kout)
            kchild(nchild) = kout
            nextc(nchild) = iact_yet  ! next action unknown
c     write(*,*) ekin, e00, nchild
 1002    enddo
         return
      end if

      do 2 iout=1,10
         call hadronmulti(iin, iout, el, amult, eldist)

c     The simplest distribution of multiplicity (takes only 2 values)
         nrp = int(amult)
         u1 = real(nrp + 1) - amult
         call rndc(u)
         if(u.gt.u1) then
            nrp = nrp + 1
         end if

         do 100 nev = 1, nrp
c           used in newer version : call ie2use(eldist, ie)
            call rndc(u)
            call hadronxdist(iin, iout, eldist, u, x)
c            write(*,*) iin,iout,x
            nchild = nchild + 1
            if(nchild.gt.maxchild) then
               stop 'Too Many child in the hadroninc interaction'
            end if

c           selection of K0l and/or K0s by random number
            if(iout.eq.6) then
               call rndc(u)
               if(u.gt.0.5d0) then
                  kout = 6
               else
                  kout = 7
               end if
            else
               kout = toout(iout)
            end if

            nextc(nchild) = iact_yet
            kchild(nchild) = kout

c           Create child particle energy from Feynman-x
            echild = e00*x + am(kout)

            if(am(kout).gt.0) then
               gamma = echild/am(kout)
               beta = sqrt(1.d0 - 1.d0/gamma**2)
               p_abs = gamma*beta*am(kout) ! momentum absolute value.
            else
c           Should not be here !!! (Gamma is discarded)
               p_abs = echild
            end if
c
c  routines for transverse momentums
c
c            write(*,*) kin,kout
            amass = am(kout)
            call sample_intcospt(kin, amass, kout, e00, p_abs, cospt)
c            write(*,*) kout, avrpt, cospt
c            write(*,'(5F10.5)') pin 
            call sample_angvec(cospt, pin, dirp)
c            write(*,*) dirp
c
c  Construction of 4+1 vector
c
            do 101 j=1, 3
               pchild(j,nchild) = p_abs*dirp(j)
 101        continue
            pchild(4,nchild) = echild
            pchild(5,nchild) = am(kout)
c            write(*,'(i3,5F10.5)') kout, (pchild(j,nchild),j=1,5)
 100     continue
 2    continue
c      write(*,*) 'hadronint ',iin,el, eksum, etsum
      end

c-------------------- 
c そのエネルギーに近い二つのエネルギー分点を探し、
c 直線補間に必要な量の計算の爲のサブルーチン
      block data beldiv
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (ndiv = 16)
      common /celdiv/ el0, eldiv(ndiv)
      data eldiv/-0.7,-0.6,-0.4, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 
     &     3.5, 4.0, 4.5, 5.0, 5.5, 6.0/
      data el0/-10./
      end

      subroutine findEdiv(el, iea, ieb, dle, elu)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (ndiv = 16)
      common /celdiv/ el0, eldiv(ndiv)
      save iea0, ieb0, dle0, elu0

      if(el .eq. el0)then
         iea = iea0
         ieb = ieb0
         dle = dle0
         elu = elu0
         return
      end if

      if(el.lt.eldiv(1)) then
         iea = 1
         ieb = 2
         dle = 0.d0
         elu = eldiv(1)
      else if(el.gt.eldiv(ndiv)) then
         iea = 15
         ieb = 16
         dle = 1.d0
         elu = eldiv(ndiv)
      else
         iea = 1
         ieb = ndiv
         do 1 while(ieb-iea.gt.1)
            iec = (iea+ieb)/2
            if(eldiv(iec).gt.el) then
               ieb = iec
            else
               iea = iec
            end if
 1       continue
         dle = (el - eldiv(iea))/(eldiv(ieb)-eldiv(iea)) 
         elu = el
      end if

      iea0 = iea
      ieb0 = ieb
      dle0 = dle
      elu0 = elu
      end
c----------------
      subroutine hadronmulti(iin,iout,El,amult,eldist)
      implicit real*8 (a-h,o-z)
      parameter (ndiv = 16)
      common /celdiv/ el0, eldiv(ndiv)

c どのエネルギー固定点のx分布を用いるか決定する。
c まずそのエネルギーに近い二つのエネルギー点を探す。
      call findEdiv(el, ie1, ie2, dle, elu)

C まず両端のエネルギー固定点での多重度のチェック

      almlt1 = hlmulti(iin,iout,ie1)
      almlt2 = hlmulti(iin,iout,ie2)

c      write(*,'(10(E10.3,I5))') el, ie1,almlt1,ie2,almlt2

c 一方のみの多重度が0の時、x分布のためのエネルギーを高い方の
c エネルギー固定点をとる。
      if(almlt1.le.-10 .and. almlt2.gt.-10)then
         eldist = eldiv(ie2)
c         eldist = el2
      else
         eldist = el
      end if

C 直線内挿で与えられた、エネルギーでの多重度を求める

      amult = 10.d0**((1.d0 - dle)*almlt1 + dle*almlt2)

      varim = (1.d0 + (1.d0 - dle)*varimlt(iin,iout,ie1)
     &     + dle*varimlt(iin,iout,ie2))

      amult = varim*amult
      end
c---- Next hadronxdist is used for no variation, such as vari-parameter
c     construction
c      subroutine hadronxdist(Iin, Iout, Einl, u0, x)
c      IMPLICIT REAL*8 (A-H,O-Z)
c      logical yet
c
cc どのエネルギー固定点のx分布を用いるか決定する。
c  まずそのエネルギーに近い二つのエネルギー点を探す。
c      call findEdiv(einl, ie1, ie2, dle, elu)
cc 乱数によって確率的にどちらを使うか決定。
c      call rndc(u)
c      if(u.gt.dle) then
c         ie = ie1
c      else
c         ie = ie2
c      end if
c
c      call distbndry(iin,iout,ie,rlmin,rlmax)
c      if(rlmax.le.rlmin) then
c         x = 1.d-6
c         return
c      end if
c
c      call dist2v(iin,iout,ie,rlmin, vmaxl)
c      vmax = 10.d0**vmaxl
c      vmax = (vmax - 1.d0/vmax)/2
c      vin = vmax*u0
c      vinl = log10(vin + sqrt(vin**2+1.d0))
c      rl = (rlmin - rlmax)*u0 + rlmax
c
c      call newtonxsolver(iin,iout,ie, vinl, rl, icon)
c      if(icon.lt.0 .or. rl.lt.rlmin .or. rl.gt.rlmax) then
cc               ニュートン法で振動状態に陥ったと判断される時はサンドイッチ法で
c         call sandwichxsolver(iin,iout,ie,rlmin,rlmax,vinl, rl)
c      end if
c
c      call dist2x(iin,iout,ie,rl, xl)
c      x = 10.d0**xl
c      icon0 = ie
c      end
c
      subroutine hadronxdist(Iin, Iout, Einl, u0, x)
      IMPLICIT REAL*8 (A-H,O-Z)
      logical yet, variation

c     どのエネルギー固定点のx分布を用いるか決定する。
c     まずそのエネルギーに近い二つのエネルギー点を探す。
      call findEdiv(einl, ie1, ie2, dle, elu)
c     乱数によって確率的にどちらを使うか決定。
      call rndc(u)
      if(u.gt.dle) then
         ie = ie1
      else
         ie = ie2
      end if

      call distbndry(iin,iout,ie,rlmin,rlmax)
      if(rlmax.le.rlmin) then
         x = 1.d-6
         return
      end if

      varia = variamp(iin,iout,ie)

      if(abs(varia).ge.1e-2) then
         variation = .true.

         call getvarip(iin,iout,ie,xl0, xint,vxint)
         famax = abs(varia) + 1.d0
      else
         variation = .false.
      end if
c      write(*,*) iin,iout,ie,varia, variation

      call dist2v(iin,iout,ie, rlmin, vmaxl)
      vmax = 10.d0**vmaxl
      vmax = (vmax - 1.d0/vmax)/2

      ntry = 0
      yet = .true.
      do 1 while(yet)
         call rndc(u0)
         vin = vmax*u0
         vinl = log10(vin + sqrt(vin**2+1.d0))
         rl = (rlmin - rlmax)*u0 + rlmax

         call newtonxsolver(iin,iout,ie, vinl, rl, icon)
         if(icon.lt.0) then
c     ニュートン法で振動状態に陥ったと判断される時はサンドイッチ法で
            call sandwichxsolver(iin,iout,ie,rlmin,rlmax,vinl, rl)
         end if
         call dist2x(iin,iout,ie,rl, xl)
c         write(*,*) u0, rl
         if(variation) then
            fa = varia*varifunc(xl0, xl) + 1.d0

            call rndc(u1)
            yet = (u1*famax .gt. fa)
c            write(*,*) iin,iout,ie,xl,varia,fa,famax,ntry,yet
            ntry = ntry + 1
         else
            yet = .false.
         end if
 1    continue
      x = 10.d0**xl
      icon0 = ie
      end

      subroutine newtonxsolver(iin,iout,ie, vin, rl, icon)
      implicit real*8 (a-h,o-z)
c     ニュートン法により解を求める。速いが、微係数の振舞いによっては
c     不安定、 また微係数を同時に求める必要あり。

      call distbndry(iin,iout,ie,rlmin,rlmax)
      call dist2v(iin,iout,ie, rl, val)
      call dist2dv(iin,iout,ie, rl, dvdx)

c      call dist2vdv(iin,iout,ie, rl, val, dvdx)
      dval = vin - val
      nturn = 0
      do 1000 while(abs(dval).gt.1e-8)

         if(abs(dvdx).gt.1e-6) then
            rl2 = rl + dval/dvdx
         else if(dval.gt.0.) then
            rl2 = (rlmin + rl)/2
         else
            rl2 = (rlmax + rl)/2
         end if

         if(rl2.lt.rlmin) then
            rl2 = (rlmin + rl)/2
         end if

         if(rl2.gt.rlmax) then
            rl2 = (rlmax + rl)/2
         end if

         rl = rl2
         call dist2v(iin,iout,ie, rl, val)
         call dist2dv(iin,iout,ie, rl, dvdx)

c         call dist2vdv(iin,iout,ie, rl, val, dvdx)
         dval = vin - val

C 10回以上繰り返すのは、振動状態に入ったと判断して、エラーフラグを立て戻る
         nturn = nturn + 1
         if(nturn.gt.20) then
            icon = -1
            return
         end if
 1000 continue
      icon = nturn
      end

      subroutine sandwichxsolver(iin,iout,ie,rlmin,rlmax,vinl,rl)
      implicit real*8 (a-h,o-z)
C 両側からはさんで解を求める。遅いが、ニュートン法より安定。
c 解の範囲を限定できる。微係数がいらない。

      rl2  = rlmin
      call dist2v(iin,iout,ie, rl2, val2)

      rl1 =  rlmax
      call dist2v(iin,iout,ie, rl1, val1)

      if(vinl.gt.val2.or.vinl.lt.val1) then
         write(0,*) 'funny in sanwitch',iin,iout,val2,vinl,val1
         stop
      end if

      do 1000 while(abs(rl1-rl2) .gt. 1e-12)
         rl3  = (rl1+rl2)/2
         call dist2v(iin,iout,ie, rl3, val3)

         if(val3.gt.vinl) then
            rl2  = rl3
            val2 = val3
         else
            rl1  = rl3
            val1 = val3
         end if
 1000 continue
      rl = rl3
      end
