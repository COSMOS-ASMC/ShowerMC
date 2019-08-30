      subroutine energyloss(kin, colden, p)
      implicit real*8 (a-h,o-z)
c definition of particles
      include '../include/atmnc-particle-code.inc'
      real*8 p(5), p1(5)
      integer icharge(21)
c                   1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 
      data icharge/ 1,-1, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 1,-1,
c                  17 18 19 20 21
     &              0, 0, 0, 0, 2/

      if(icharge(kin).eq.0) then
         return
      end if

      ek = p(4) - p(5)
      if(abs(colden) .lt. 100.d0) then
         de = dedx(kin, ek)*colden
         ek1 = ek - de

         if(ek1 .lt. 1.d-3) then
            ek1 = 1.d-3
         end if
      else
         call elmagrange(kin, ek, range)

         range1 = range - colden
         if(range1.gt.0) then
            call invemrange(kin, range1, ek1)
         else
            ek1 = 1.d-3
         end if
      end if
c                                 When particl Ek < 1.e-3, it's stopped
      if(ek1 .gt. 1.d-3) then
         gamma0 = p(4)/p(5)
         gamma1 = (ek1 + p(5))/p(5)
         fact = sqrt((gamma1**2 - 1.d0)/(gamma0**2 - 1.d0))
      else
         gamma1 = 1.d0
         fact = 0.
      end if

      do 1 i=1, 3
         p(i) = p(i)*fact
 1    continue
      p(4) = p(5)*gamma1
      end

      function dedx(kin0, ek0)
      implicit real*8 (a-h,o-z)
      include '../include/atmnc-particle-code.inc'
      logical alpha

      kin = kin0
      ek = ek0
      if(kin.eq.kalpha) then
         alpha = .true.
         kin = kpro
         ek = ek/4
      else
         alpha = .false.
      end if

      if(ek.lt.1.d-3) then
         ek = 1.d-3
      end if

      if(kin.eq.kep) then
         DEDX=dex_ep (ek)
      else if(kin.eq.kem) then
         DEDX=dex_em (ek)
      else if(kin.eq.kkp.or.kin.eq.kkm) then
         DEDX=dex_k (ek)
      else if(kin.eq.kpip.or.kin.eq.kpim) then
         DEDX=dex_pi(ek)
      else if(kin.eq.kpro.or.kin.eq.kprobar) then
         DEDX=dex_p (ek)
      else if(kin.eq.kmup.or.kin.eq.kmum) then
         DEDX=dex_mu(ek)
      else
         dedx = 0.
      end if

      if(alpha) then
         dedx = 4*dedx
      end if
      end

      subroutine elmagrange(kin0, ek0, cold) ! range in column density
      implicit real*8 (a-h,o-z)
      include '../include/atmnc-particle-code.inc' ! particl code for atmnc
      logical alpha

      ek = ek0
      kin = kin0
      if(kin.eq.kalpha) then
         alpha = .true.
         kin = 10
         ek = ek/4
      else
         alpha = .false.
      end if

      if(ek.le.1e-3) then
         cold = 0.0
         return
      else if(ek.gt.1e6) then
         cold = 1.e10
         return
      end if

      if(kin.eq.kep) then
         COLD=rng_ep (ek)
      else if(kin.eq.kem) then
         COLD=rng_em (ek)
      else if(kin.eq.kkp.or.kin.eq.kkm) then
         COLD=rng_k (ek)
      else if(kin.eq.kpip.or.kin.eq.kpim) then
         COLD=rng_pi(ek)
      else if(kin.eq.kpro.or.kin.eq.kprobar) then
         COLD=rng_p (ek)
      else if(kin.eq.kmup.or.kin.eq.kmum) then
         COLD=rng_mu(ek)
      else
         cold = 1.e10
      end if

      end

      subroutine invemrange(kin0, cold, Ek)
      implicit real*8 (a-h,i,o-z)
      include '../include/atmnc-particle-code.inc' ! particl code for atmnc

      external rng_ep, rng_em, rng_k, rng_p, rng_pi, rng_mu
      external dex_ep, dex_em, dex_k, dex_p, dex_pi, dex_mu

      logical alpha

      kin = kin0
      if(kin.eq.kalpha) then
         alpha = .true.
         kin = kpro
      else
         alpha = .false.
      end if

      if(cold.le.1e-5) then
         ek = 1.d-3
         return
      end if

      if(kin.eq.kep) then
         ek0 = irg_ep (cold)
         call newtonrngsolver(rng_ep, dex_ep, cold, ek0, ek)
      else if(kin.eq.kem) then
         ek0 = irg_em (cold)
         call newtonrngsolver(rng_em, dex_em, cold, ek0, ek)
      else if(kin.eq.kkp.or.kin.eq.kkm) then
         ek0 = irg_k (cold)
         call newtonrngsolver(rng_k, dex_k, cold, ek0, ek)
      else if(kin.eq.kpip.or.kin.eq.kpim) then
         ek0 = irg_pi(cold)
         call newtonrngsolver(rng_pi, dex_pi, cold, ek0, ek)
      else if(kin.eq.kpro.or.kin.eq.kprobar) then
         ek0 = irg_p (cold)
         call newtonrngsolver(rng_p, dex_p, cold, ek0, ek)
      else if(kin.eq.kmup.or.kin.eq.kmum) then
         ek0 = irg_mu(cold)
         call newtonrngsolver(rng_mu, dex_mu, cold, ek0, ek)
      else
         ek = 1.D-3
      end if

      if(alpha) then
         ek = 4*ek
      end if
C               For the stability of the simulation give nonzero ek
      if(ek.lt.1.d-3) then
         ek = 1.d-3
      end if
c      write(*,*)'invemrange', kin,cold,ek0,ek
      end

      subroutine newtonrngsolver(func, difunc, cold, ek0, ek)
      implicit real*8 (a-h,o-z)
c ニュートン法により解を求める。速いが、微係数の振舞いによっては不安定?
c また微係数を同時に求める必要あり。
      parameter (ekmin = 1.d-3, ekmax=1.d6)

      ek = ek0
      if(ek0.lt.ekmin .or. ek.gt.ekmax) then
         return
      end if

      cold1 = func(ek)
      dval = cold - cold1

      nturn = 0
      do 1000 while(abs(dval)/cold.gt.1e-10)
c         write(*,*) nturn,cold1,cold, ek0, ek
         ek = ek + dval*difunc(ek)

         cold1 = func(ek)
         dval = cold - cold1
         nturn = nturn + 1
         if(nturn.gt.10) go to 2000
 1000 continue
      return
c                  for non conversing case
 2000 continue
      ek1 = ekmin
      ek2 = ekmax
      do 2100 while(ek2.gt.ek1+1.d-3)
         ek3 = (ek1+ek2)/2
         cld3 = func(ek3)
         if(cld3.gt.cold) then
            ek2 = ek3
         else
            ek1 = ek3
         end if
 2100 continue
      ek = (ek1+ek2)/2
      end

      FUNCTION rng_ep(X)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 CF(200),X,X0,DX,rng_ep
      INTEGER*4 NC, K
      LOGICAL XL,YL
      Data (CF(K),K=1, 43)/
     &-.984350E+00,-.257273E+00,-.139656E+00,-.642706E-01,-.846426E-02,
     &0.388565E-01,0.777208E-01,0.121683E+00,0.160241E+00,0.198969E+00,
     &0.236160E+00,0.273391E+00,0.309615E+00,0.344687E+00,0.377618E+00,
     &0.407699E+00,0.434254E+00,0.457133E+00,0.476552E+00,0.493000E+00,
     &0.507007E+00,0.519063E+00,0.529564E+00,0.538821E+00,0.547073E+00,
     &0.554506E+00,0.561259E+00,0.567442E+00,0.573143E+00,0.578430E+00,
     &0.583358E+00,0.587973E+00,0.592311E+00,0.596404E+00,0.600278E+00,
     &0.603955E+00,0.607455E+00,0.610793E+00,0.613983E+00,0.617039E+00,
     &0.619971E+00,0.622789E+00,0.625501E+00/
      Data NC/  43/
      Data X0, DX/-.375000E+01,0.250000E+00/
      Data XL,YL/.TRUE.,.TRUE./
      rng_ep=BSPFIT(NC,CF,X0,DX,XL,YL,X)
      END
      FUNCTION rng_em(X)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 CF(200),X,X0,DX,rng_em
      INTEGER*4 NC, K
      LOGICAL XL,YL
      Data (CF(K),K=1, 43)/
     &-.988729E+00,-.261656E+00,-.144307E+00,-.690476E-01,-.133159E-01,
     &0.341475E-01,0.733791E-01,0.117368E+00,0.156160E+00,0.195079E+00,
     &0.232508E+00,0.269954E+00,0.306417E+00,0.341747E+00,0.374975E+00,
     &0.405381E+00,0.432270E+00,0.455460E+00,0.475147E+00,0.491814E+00,
     &0.505997E+00,0.518191E+00,0.528800E+00,0.538145E+00,0.546468E+00,
     &0.553959E+00,0.560760E+00,0.566984E+00,0.572720E+00,0.578036E+00,
     &0.582990E+00,0.587627E+00,0.591986E+00,0.596097E+00,0.599987E+00,
     &0.603678E+00,0.607191E+00,0.610541E+00,0.613742E+00,0.616808E+00,
     &0.619749E+00,0.622576E+00,0.625296E+00/
      Data NC/  43/
      Data X0, DX/-.375000E+01,0.250000E+00/
      Data XL,YL/.TRUE.,.TRUE./
      rng_em=BSPFIT(NC,CF,X0,DX,XL,YL,X)
      END
      FUNCTION rng_mu(X)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 CF(200),X,X0,DX,rng_mu
      INTEGER*4 NC, K
      LOGICAL XL,YL
      Data (CF(K),K=1, 53)/
     &-.557632E+00,-.324041E+00,-.231293E+00,-.158540E+00,-.940131E-01,
     &-.315557E-01,0.295517E-01,0.902187E-01,0.150203E+00,0.209337E+00,
     &0.267048E+00,0.322691E+00,0.375478E+00,0.424723E+00,0.470027E+00,
     &0.511446E+00,0.549435E+00,0.584685E+00,0.617901E+00,0.649711E+00,
     &0.680585E+00,0.710912E+00,0.740877E+00,0.770882E+00,0.801078E+00,
     &0.831448E+00,0.861976E+00,0.892565E+00,0.923073E+00,0.953277E+00,
     &0.982878E+00,0.101148E+01,0.103863E+01,0.106384E+01,0.108674E+01,
     &0.110711E+01,0.112498E+01,0.114052E+01,0.115402E+01,0.116578E+01,
     &0.117610E+01,0.118522E+01,0.119335E+01,0.120066E+01,0.120727E+01,
     &0.121331E+01,0.121884E+01,0.122395E+01,0.122869E+01,0.123310E+01,
     &0.123723E+01,0.124110E+01,0.124475E+01/
      Data NC/  53/
      Data X0, DX/-.360000E+01,0.200000E+00/
      Data XL,YL/.TRUE.,.TRUE./
      rng_mu=BSPFIT(NC,CF,X0,DX,XL,YL,X)
      END
      FUNCTION rng_pi(X)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 CF(200),X,X0,DX,rng_pi
      INTEGER*4 NC, K
      LOGICAL XL,YL
      Data (CF(K),K=1, 53)/
     &-.572719E+00,-.339387E+00,-.246912E+00,-.174423E+00,-.110120E+00,
     &-.478211E-01,0.132272E-01,0.739740E-01,0.134232E+00,0.193907E+00,
     &0.252514E+00,0.309495E+00,0.364092E+00,0.415549E+00,0.463256E+00,
     &0.506989E+00,0.546975E+00,0.583789E+00,0.618146E+00,0.650731E+00,
     &0.682108E+00,0.712713E+00,0.742854E+00,0.772755E+00,0.802822E+00,
     &0.833060E+00,0.863500E+00,0.894091E+00,0.924746E+00,0.955312E+00,
     &0.985573E+00,0.101522E+01,0.104387E+01,0.107106E+01,0.109632E+01,
     &0.111925E+01,0.113966E+01,0.115756E+01,0.117312E+01,0.118663E+01,
     &0.119841E+01,0.120873E+01,0.121786E+01,0.122598E+01,0.123329E+01,
     &0.123990E+01,0.124592E+01,0.125145E+01,0.125655E+01,0.126127E+01,
     &0.126568E+01,0.126979E+01,0.127365E+01/
      Data NC/  53/
      Data X0, DX/-.360000E+01,0.200000E+00/
      Data XL,YL/.TRUE.,.TRUE./
      rng_pi=BSPFIT(NC,CF,X0,DX,XL,YL,X)
      END
      FUNCTION rng_k (X)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 CF(200),X,X0,DX,rng_k 
      INTEGER*4 NC, K
      LOGICAL XL,YL
      Data (CF(K),K=1, 53)/
     &-.635542E+00,-.404121E+00,-.313627E+00,-.243068E+00,-.180524E+00,
     &-.119708E+00,-.597791E-01,0.301807E-03,0.604703E-01,0.120814E+00,
     &0.181112E+00,0.241145E+00,0.300528E+00,0.358764E+00,0.415193E+00,
     &0.469051E+00,0.519585E+00,0.566270E+00,0.608997E+00,0.648096E+00,
     &0.684200E+00,0.718025E+00,0.750233E+00,0.781343E+00,0.811774E+00,
     &0.841766E+00,0.871621E+00,0.901663E+00,0.931882E+00,0.962335E+00,
     &0.993001E+00,0.102386E+01,0.105483E+01,0.108581E+01,0.111663E+01,
     &0.114708E+01,0.117682E+01,0.120547E+01,0.123256E+01,0.125764E+01,
     &0.128034E+01,0.130050E+01,0.131816E+01,0.133351E+01,0.134684E+01,
     &0.135845E+01,0.136864E+01,0.137766E+01,0.138569E+01,0.139291E+01,
     &0.139946E+01,0.140542E+01,0.141090E+01/
      Data NC/  53/
      Data X0, DX/-.360000E+01,0.200000E+00/
      Data XL,YL/.TRUE.,.TRUE./
      rng_k =BSPFIT(NC,CF,X0,DX,XL,YL,X)
      END
      FUNCTION rng_p (X)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 CF(200),X,X0,DX,rng_p 
      INTEGER*4 NC, K
      LOGICAL XL,YL
      Data (CF(K),K=1, 53)/
     &-.661678E+00,-.432141E+00,-.343483E+00,-.274652E+00,-.213667E+00,
     &-.154198E+00,-.953893E-01,-.361958E-01,0.233211E-01,0.832664E-01,
     &0.143462E+00,0.203769E+00,0.263928E+00,0.323615E+00,0.382374E+00,
     &0.439607E+00,0.494581E+00,0.546516E+00,0.594770E+00,0.639052E+00,
     &0.679528E+00,0.716735E+00,0.751381E+00,0.784156E+00,0.815646E+00,
     &0.846295E+00,0.876439E+00,0.906275E+00,0.936226E+00,0.966348E+00,
     &0.996686E+00,0.102725E+01,0.105804E+01,0.108902E+01,0.112014E+01,
     &0.115133E+01,0.118244E+01,0.121331E+01,0.124367E+01,0.127319E+01,
     &0.130145E+01,0.132799E+01,0.135239E+01,0.137435E+01,0.139377E+01,
     &0.141073E+01,0.142547E+01,0.143830E+01,0.144950E+01,0.145937E+01,
     &0.146812E+01,0.147595E+01,0.148301E+01/
      Data NC/  53/
      Data X0, DX/-.360000E+01,0.200000E+00/
      Data XL,YL/.TRUE.,.TRUE./
      rng_p =BSPFIT(NC,CF,X0,DX,XL,YL,X)
      END
      FUNCTION irg_ep(X)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 CF(200),X,X0,DX,irg_ep
      INTEGER*4 NC, K
      LOGICAL XL,YL
      Data (CF(K),K=1, 31)/
     &-.927946E+01,0.000000E+00,-.577878E+00,-.566527E+00,-.560215E+00,
     &-.549695E+00,-.536484E+00,-.519730E+00,-.500048E+00,-.476725E+00,
     &-.450879E+00,-.421795E+00,-.386243E+00,-.354532E+00,-.320434E+00,
     &-.285070E+00,-.248627E+00,-.211589E+00,-.174189E+00,-.136150E+00,
     &-.970060E-01,-.559108E-01,-.112993E-01,0.393160E-01,0.100890E+00,
     &0.181755E+00,0.298582E+00,0.477552E+00,0.761910E+00,0.120858E+01,
     &0.197850E+01/
      Data NC/  31/
      Data X0, DX/-.200000E+01,0.200000E+00/
      Data XL,YL/.TRUE.,.TRUE./
      irg_ep=BSPFIT(NC,CF,X0,DX,XL,YL,X)
      END
      FUNCTION irg_em(X)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 CF(200),X,X0,DX,irg_em
      INTEGER*4 NC, K
      LOGICAL XL,YL
      Data (CF(K),K=1, 31)/
     &-.593650E+01,0.000000E+00,-.573753E+00,-.566256E+00,-.558914E+00,
     &-.548058E+00,-.534365E+00,-.517094E+00,-.496978E+00,-.473032E+00,
     &-.447092E+00,-.416931E+00,-.381834E+00,-.350344E+00,-.315952E+00,
     &-.280744E+00,-.244342E+00,-.207524E+00,-.170274E+00,-.132404E+00,
     &-.933914E-01,-.524158E-01,-.789805E-02,0.426404E-01,0.104159E+00,
     &0.184987E+00,0.301796E+00,0.480761E+00,0.765125E+00,0.121176E+01,
     &0.198376E+01/
      Data NC/  31/
      Data X0, DX/-.200000E+01,0.200000E+00/
      Data XL,YL/.TRUE.,.TRUE./
      irg_em=BSPFIT(NC,CF,X0,DX,XL,YL,X)
      END
      FUNCTION irg_mu(X)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 CF(200),X,X0,DX,irg_mu
      INTEGER*4 NC, K
      LOGICAL XL,YL
      Data (CF(K),K=1, 54)/
     &-.735444E+03,0.000000E+00,0.956559E-01,-.591126E+00,-.563802E+00,
     &-.555937E+00,-.545559E+00,-.533133E+00,-.518952E+00,-.503410E+00,
     &-.486891E+00,-.469708E+00,-.452090E+00,-.434193E+00,-.416115E+00,
     &-.397912E+00,-.379608E+00,-.361209E+00,-.342700E+00,-.324052E+00,
     &-.305222E+00,-.286151E+00,-.266759E+00,-.246945E+00,-.226579E+00,
     &-.205495E+00,-.183491E+00,-.160319E+00,-.135701E+00,-.109352E+00,
     &-.810284E-01,-.506151E-01,-.181779E-01,0.159973E-01,0.515258E-01,
     &0.879316E-01,0.124921E+00,0.162010E+00,0.198817E+00,0.235404E+00,
     &0.271780E+00,0.308104E+00,0.344590E+00,0.381601E+00,0.419694E+00,
     &0.459789E+00,0.503435E+00,0.553369E+00,0.614638E+00,0.696860E+00,
     &0.818580E+00,0.101353E+01,0.134034E+01,0.203985E+01/
      Data NC/  54/
      Data X0, DX/-.280000E+01,0.200000E+00/
      Data XL,YL/.TRUE.,.TRUE./
      irg_mu=BSPFIT(NC,CF,X0,DX,XL,YL,X)
      END
      FUNCTION irg_pi(X)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 CF(200),X,X0,DX,irg_pi
      INTEGER*4 NC, K
      LOGICAL XL,YL
      Data (CF(K),K=1, 55)/
     &-.836398E+01,0.000000E+00,-.572907E+00,-.567645E+00,-.560603E+00,
     &-.551384E+00,-.540001E+00,-.526663E+00,-.511737E+00,-.495628E+00,
     &-.478698E+00,-.461226E+00,-.443411E+00,-.425384E+00,-.407223E+00,
     &-.388974E+00,-.370656E+00,-.352268E+00,-.333798E+00,-.315220E+00,
     &-.296494E+00,-.277570E+00,-.258377E+00,-.238829E+00,-.218811E+00,
     &-.198181E+00,-.176759E+00,-.154323E+00,-.130612E+00,-.105340E+00,
     &-.782310E-01,-.490803E-01,-.178410E-01,0.153375E-01,0.501087E-01,
     &0.860607E-01,0.122743E+00,0.159888E+00,0.196892E+00,0.233638E+00,
     &0.270130E+00,0.306427E+00,0.342676E+00,0.379104E+00,0.416073E+00,
     &0.454146E+00,0.494253E+00,0.537949E+00,0.588016E+00,0.649525E+00,
     &0.732277E+00,0.855037E+00,0.105244E+01,0.138428E+01,0.398583E+01/
      Data NC/  55/
      Data X0, DX/-.280000E+01,0.200000E+00/
      Data XL,YL/.TRUE.,.TRUE./
      irg_pi=BSPFIT(NC,CF,X0,DX,XL,YL,X)
      END
      FUNCTION irg_k (X)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 CF(200),X,X0,DX,irg_k 
      INTEGER*4 NC, K
      LOGICAL XL,YL
      Data (CF(K),K=1, 61)/
     &0.000000E+00,-.608410E+00,-.573422E+00,-.567979E+00,-.561002E+00,
     &-.551719E+00,-.540151E+00,-.526488E+00,-.511125E+00,-.494510E+00,
     &-.477053E+00,-.459070E+00,-.440786E+00,-.422348E+00,-.403846E+00,
     &-.385332E+00,-.366835E+00,-.348364E+00,-.329921E+00,-.311496E+00,
     &-.293078E+00,-.274647E+00,-.256181E+00,-.237648E+00,-.219012E+00,
     &-.200226E+00,-.181230E+00,-.161949E+00,-.142288E+00,-.122125E+00,
     &-.101308E+00,-.796452E-01,-.569050E-01,-.328173E-01,-.709336E-02,
     &0.205332E-01,0.502330E-01,0.820086E-01,0.115653E+00,0.150802E+00,
     &0.187017E+00,0.223913E+00,0.261155E+00,0.298168E+00,0.334936E+00,
     &0.371396E+00,0.407582E+00,0.443546E+00,0.479399E+00,0.515306E+00,
     &0.551528E+00,0.588460E+00,0.626737E+00,0.667391E+00,0.712208E+00,
     &0.764315E+00,0.829644E+00,0.919487E+00,0.105587E+01,0.127755E+01,
     &0.187210E+01/
      Data NC/  61/
      Data X0, DX/-.320000E+01,0.200000E+00/
      Data XL,YL/.TRUE.,.TRUE./
      irg_k =BSPFIT(NC,CF,X0,DX,XL,YL,X)
      END
      FUNCTION irg_p (X)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 CF(200),X,X0,DX,irg_p 
      INTEGER*4 NC, K
      LOGICAL XL,YL
      Data (CF(K),K=1, 64)/
     &-.254163E+02,0.000000E+00,-.580805E+00,-.568344E+00,-.562261E+00,
     &-.553137E+00,-.541649E+00,-.527895E+00,-.512283E+00,-.495306E+00,
     &-.477428E+00,-.459010E+00,-.440308E+00,-.421485E+00,-.402640E+00,
     &-.383828E+00,-.365075E+00,-.346393E+00,-.327780E+00,-.309231E+00,
     &-.290735E+00,-.272278E+00,-.253843E+00,-.235409E+00,-.216953E+00,
     &-.198447E+00,-.179855E+00,-.161133E+00,-.142227E+00,-.123069E+00,
     &-.103571E+00,-.836241E-01,-.630872E-01,-.417856E-01,-.195034E-01,
     &0.401505E-02,0.290554E-01,0.558988E-01,0.847660E-01,0.115730E+00,
     &0.148674E+00,0.183274E+00,0.219134E+00,0.255790E+00,0.292986E+00,
     &0.330136E+00,0.367023E+00,0.403636E+00,0.439941E+00,0.475973E+00,
     &0.511775E+00,0.547434E+00,0.583078E+00,0.618904E+00,0.655211E+00,
     &0.692470E+00,0.731432E+00,0.773348E+00,0.820368E+00,0.876430E+00,
     &0.948878E+00,0.105203E+01,0.121200E+01,0.150981E+01/
      Data NC/  64/
      Data X0, DX/-.340000E+01,0.200000E+00/
      Data XL,YL/.TRUE.,.TRUE./
      irg_p =BSPFIT(NC,CF,X0,DX,XL,YL,X)
      END
      FUNCTION dex_ep(X)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 CF(200),X,X0,DX,dex_ep
      INTEGER*4 NC, K
      LOGICAL XL,YL
      Data (CF(K),K=1, 43)/
     &0.439985E+00,0.445814E+00,0.462752E+00,0.462457E+00,0.470111E+00,
     &0.450836E+00,0.455395E+00,0.453685E+00,0.446324E+00,0.441973E+00,
     &0.436391E+00,0.431049E+00,0.423577E+00,0.413641E+00,0.399451E+00,
     &0.380106E+00,0.355039E+00,0.324813E+00,0.290506E+00,0.253370E+00,
     &0.214405E+00,0.174323E+00,0.133573E+00,0.924351E-01,0.510700E-01,
     &0.957603E-02,-.319940E-01,-.736057E-01,-.115243E+00,-.156894E+00,
     &-.198553E+00,-.240216E+00,-.281883E+00,-.323551E+00,-.365220E+00,
     &-.406890E+00,-.448560E+00,-.490230E+00,-.531901E+00,-.573572E+00,
     &-.615243E+00,-.656914E+00,-.698594E+00/
      Data NC/  43/
      Data X0, DX/-.375000E+01,0.250000E+00/
      Data XL,YL/.TRUE.,.TRUE./
      dex_ep=BSPFIT(NC,CF,X0,DX,XL,YL,X)
      END
      FUNCTION dex_em(X)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 CF(200),X,X0,DX,dex_em
      INTEGER*4 NC, K
      LOGICAL XL,YL
      Data (CF(K),K=1, 43)/
     &0.435461E+00,0.441545E+00,0.457820E+00,0.457708E+00,0.465022E+00,
     &0.446916E+00,0.451220E+00,0.449620E+00,0.442646E+00,0.438498E+00,
     &0.433174E+00,0.428054E+00,0.420875E+00,0.411283E+00,0.397517E+00,
     &0.378630E+00,0.354003E+00,0.324134E+00,0.290086E+00,0.253119E+00,
     &0.214259E+00,0.174239E+00,0.133526E+00,0.924080E-01,0.510552E-01,
     &0.956710E-02,-.319988E-01,-.736083E-01,-.115244E+00,-.156895E+00,
     &-.198553E+00,-.240217E+00,-.281883E+00,-.323551E+00,-.365221E+00,
     &-.406890E+00,-.448560E+00,-.490230E+00,-.531901E+00,-.573572E+00,
     &-.615243E+00,-.656914E+00,-.698594E+00/
      Data NC/  43/
      Data X0, DX/-.375000E+01,0.250000E+00/
      Data XL,YL/.TRUE.,.TRUE./
      dex_em=BSPFIT(NC,CF,X0,DX,XL,YL,X)
      END
      FUNCTION dex_mu(X)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 CF(200),X,X0,DX,dex_mu
      INTEGER*4 NC, K
      LOGICAL XL,YL
      Data (CF(K),K=1, 53)/
     &-.313915E+00,-.339331E+00,-.365431E+00,-.392039E+00,-.418985E+00,
     &-.446084E+00,-.473106E+00,-.499744E+00,-.525573E+00,-.550002E+00,
     &-.572252E+00,-.591380E+00,-.606441E+00,-.616788E+00,-.622409E+00,
     &-.624008E+00,-.622743E+00,-.619758E+00,-.615928E+00,-.611778E+00,
     &-.607627E+00,-.603573E+00,-.599698E+00,-.596886E+00,-.594109E+00,
     &-.591613E+00,-.589064E+00,-.586335E+00,-.583079E+00,-.578909E+00,
     &-.573265E+00,-.565473E+00,-.554784E+00,-.540546E+00,-.522399E+00,
     &-.500427E+00,-.475102E+00,-.447108E+00,-.417146E+00,-.385811E+00,
     &-.353557E+00,-.320703E+00,-.287462E+00,-.253974E+00,-.220328E+00,
     &-.186580E+00,-.152765E+00,-.118905E+00,-.850131E-01,-.510966E-01,
     &-.171611E-01,0.167902E-01,0.507549E-01/
      Data NC/  53/
      Data X0, DX/-.360000E+01,0.200000E+00/
      Data XL,YL/.TRUE.,.TRUE./
      dex_mu=BSPFIT(NC,CF,X0,DX,XL,YL,X)
      END
      FUNCTION dex_pi(X)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 CF(200),X,X0,DX,dex_pi
      INTEGER*4 NC, K
      LOGICAL XL,YL
      Data (CF(K),K=1, 53)/
     &-.298948E+00,-.323855E+00,-.349566E+00,-.375888E+00,-.402652E+00,
     &-.429683E+00,-.456783E+00,-.483696E+00,-.510074E+00,-.535430E+00,
     &-.559102E+00,-.580234E+00,-.597857E+00,-.611100E+00,-.619534E+00,
     &-.623442E+00,-.623752E+00,-.621680E+00,-.618253E+00,-.614255E+00,
     &-.610016E+00,-.605946E+00,-.601782E+00,-.598355E+00,-.595586E+00,
     &-.592914E+00,-.590483E+00,-.588002E+00,-.585327E+00,-.582123E+00,
     &-.577999E+00,-.572399E+00,-.564650E+00,-.554004E+00,-.539806E+00,
     &-.521699E+00,-.499762E+00,-.474460E+00,-.446476E+00,-.416511E+00,
     &-.385163E+00,-.352887E+00,-.320006E+00,-.286734E+00,-.253213E+00,
     &-.219533E+00,-.185750E+00,-.151899E+00,-.118002E+00,-.840731E-01,
     &-.501199E-01,-.161478E-01,0.178400E-01/
      Data NC/  53/
      Data X0, DX/-.360000E+01,0.200000E+00/
      Data XL,YL/.TRUE.,.TRUE./
      dex_pi=BSPFIT(NC,CF,X0,DX,XL,YL,X)
      END
      FUNCTION dex_k (X)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 CF(200),X,X0,DX,dex_k 
      INTEGER*4 NC, K
      LOGICAL XL,YL
      Data (CF(K),K=1, 53)/
     &-.237000E+00,-.258163E+00,-.281043E+00,-.305205E+00,-.330344E+00,
     &-.356233E+00,-.382688E+00,-.409540E+00,-.436612E+00,-.463693E+00,
     &-.490508E+00,-.516677E+00,-.541673E+00,-.564780E+00,-.585104E+00,
     &-.601672E+00,-.613692E+00,-.620899E+00,-.623767E+00,-.623335E+00,
     &-.620812E+00,-.617148E+00,-.613035E+00,-.608743E+00,-.604635E+00,
     &-.600418E+00,-.597207E+00,-.594336E+00,-.591736E+00,-.589393E+00,
     &-.587213E+00,-.585104E+00,-.582909E+00,-.580413E+00,-.577309E+00,
     &-.573175E+00,-.567452E+00,-.559452E+00,-.548449E+00,-.533838E+00,
     &-.515316E+00,-.493004E+00,-.467396E+00,-.439180E+00,-.409047E+00,
     &-.377582E+00,-.345224E+00,-.312283E+00,-.278967E+00,-.245412E+00,
     &-.211704E+00,-.177898E+00,-.144028E+00/
      Data NC/  53/
      Data X0, DX/-.360000E+01,0.200000E+00/
      Data XL,YL/.TRUE.,.TRUE./
      dex_k =BSPFIT(NC,CF,X0,DX,XL,YL,X)
      END
      FUNCTION dex_p (X)
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 CF(200),X,X0,DX,dex_p 
      INTEGER*4 NC, K
      LOGICAL XL,YL
      Data (CF(K),K=1, 53)/
     &-.211700E+00,-.229202E+00,-.249503E+00,-.271765E+00,-.295461E+00,
     &-.320243E+00,-.345859E+00,-.372110E+00,-.398825E+00,-.425833E+00,
     &-.452940E+00,-.479899E+00,-.506379E+00,-.531916E+00,-.555872E+00,
     &-.577418E+00,-.595590E+00,-.609486E+00,-.618590E+00,-.623080E+00,
     &-.623811E+00,-.621986E+00,-.618670E+00,-.614683E+00,-.610407E+00,
     &-.606233E+00,-.601993E+00,-.598275E+00,-.595387E+00,-.592582E+00,
     &-.590164E+00,-.587915E+00,-.585876E+00,-.583919E+00,-.581960E+00,
     &-.579817E+00,-.577258E+00,-.573942E+00,-.569409E+00,-.563061E+00,
     &-.554213E+00,-.542168E+00,-.526401E+00,-.506738E+00,-.483428E+00,
     &-.457041E+00,-.428277E+00,-.397796E+00,-.366136E+00,-.333692E+00,
     &-.300739E+00,-.267460E+00,-.233973E+00/
      Data NC/  53/
      Data X0, DX/-.360000E+01,0.200000E+00/
      Data XL,YL/.TRUE.,.TRUE./
      dex_p =BSPFIT(NC,CF,X0,DX,XL,YL,X)
      END
