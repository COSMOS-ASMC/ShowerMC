       subroutine cpimuPolari(pion, muon,  polari)
!        pion: /ptcl/. input.  one pion data. charge and energy are used.
!        muon: /ptcl/. input.  muon from the pion. energy is used.
!      polari: real*8. output. polarization of muon in lab frame along
!              its momentum.
!
!        polarization of mu at lab is obtained for pi--->mu decay.
!    *** note ***
!       For mu-, polari is mostly positive hence decay eletron
!       goes opposit side of mu and neutrino goes same
!       side (use csampNeueEMu and 
!                 (1+P cos)dcos for neue;
!                 csNeumuEMu and
!                 (1+ XPcos)dcos for neumu;
!                X=(1-2f)/(3-2f)
!       for energy and angle sampling, use porali as it is)
!       For mu+,  polari is mostly negative  but positron goes
!       the opposit side of mu and hence
!       neutrinos goes like mu- case so that use csampNeumuCos etc
!       by reversing the sign of polarization).

       implicit none
!----       include '../../../Zptcl.h'
#include  "Zptcl.h"
!----       include '../../../Zmass.h'
#include  "Zmass.h"
       type(ptcl):: pion
       type(ptcl):: muon
       real*8  polari
!
       real*8  masmu2, est, pst
       parameter (masmu2 = masmu**2, est=(maspic**2+ masmu2)/2/maspic,
     *     pst=(maspic**2-masmu2)/2/maspic)
!
       real*8 g, pmu

           g=pion%fm%p(4)/maspic
           pmu = muon%fm%p(4)**2- masmu2
           if(pmu .lt. 0.) then           
!//////////////////
!             some of stopping pi result in pmu< 0 (Ppi 
!              write(0,*) ' pmu =',pmu, muon.fm.p(4), muon.mass
!              call checkstat("in cpimuPolari")
!///////////////
              pmu =0.
           endif

           pmu=sqrt(pmu)
           polari=(muon%fm%p(4) * est - g * masmu2)/pmu/pst
           if(muon%charge .gt. 0) then
               polari = -polari
           endif
           if(abs(polari) .gt. 1.) then
              polari = sign(1.d0, polari)
           endif
        end
       subroutine ckmuPolari(kaon,  muon,  polari)
!         k----->mu+nuew decay.  polarization of mu at lab.
!        
!         kaon: /ptcl/.  input. charge and energy is used.
!         muon: /ptcl/.  input. enegy is used.
!       polari: real*8.  output. muon polarizaiton.
       implicit none
!----       include '../../../Zptcl.h'
#include  "Zptcl.h"
!----       include '../../../Zmass.h'
#include  "Zmass.h"
       type(ptcl):: kaon, muon
       real*8 polari

       real*8  masmu2,  est, pst

       parameter(masmu2=masmu**2, est=(maskc**2+ masmu2)/2/maskc,
     *     pst=(maskc**2-masmu2)/2/maskc)
!
       real*8  g, pmu
!
           g = kaon%fm%p(4)/maskc
           pmu=sqrt(muon%fm%p(4)**2- masmu2)
           polari = (muon%fm%p(4)*est - g * masmu2)/pmu/pst
           if(kaon%charge .gt. 0) then
               polari = -polari
           endif
           if(abs(polari) .gt. 1.) then
              polari = sign(1.d0, polari)
           endif
        end
!       ****************************************************************
!       *
!       * csampNeuEKl3:  sample energy of neutrino from kl3 decay.
!       * csampMuEKl3:  sample energy of muon from kl3 decay.
!       *          approx by  k mass is k+-, electron mass=0
!       * cmuPolAtK:  longitudinal polarization of mu at k-rest
!       * cmuPolAtLabK:  //                                 lab.
!       *
!       ************************** tested 88.07.27 ***********k.k*******
        subroutine csampNeuEKl3(f)
        implicit none
!----        include '../../../Zmass.h'
#include  "Zmass.h"

        integer i
        real*8 f
          real*8 mpmk2, snorm, f1, u
          integer l
          parameter ( snorm=2.43e-2,
     *        mpmk2=(maspic/maskc)**2,
     *        f1=1.7678*snorm/(1.-mpmk2) )

!            f=e/mk  (=0 to (1-(mass_pi/mass_k)**2))/2)
!        neutrino energy sampling table
      real*8 fn(101)
      data (fn    (i),i=   1,  72)/
     1 0.0000, 0.0618, 0.0789, 0.0911, 0.1010, 0.1094, 0.1168, 0.1236,
     2 0.1300, 0.1356, 0.1412, 0.1462, 0.1512, 0.1558, 0.1603, 0.1645,
     3 0.1688, 0.1727, 0.1766, 0.1805, 0.1841, 0.1877, 0.1912, 0.1946,
     4 0.1980, 0.2013, 0.2045, 0.2077, 0.2109, 0.2139, 0.2170, 0.2200,
     5 0.2229, 0.2258, 0.2288, 0.2316, 0.2344, 0.2372, 0.2400, 0.2427,
     6 0.2454, 0.2482, 0.2508, 0.2535, 0.2561, 0.2588, 0.2614, 0.2640,
     7 0.2666, 0.2692, 0.2717, 0.2743, 0.2768, 0.2794, 0.2819, 0.2845,
     8 0.2870, 0.2895, 0.2921, 0.2946, 0.2971, 0.2997, 0.3022, 0.3048,
     9 0.3073, 0.3099, 0.3124, 0.3150, 0.3176, 0.3202, 0.3228, 0.3255/
      data (fn    (i),i=  73, 101)/
     1 0.3281, 0.3308, 0.3335, 0.3362, 0.3389, 0.3417, 0.3446, 0.3474,
     2 0.3503, 0.3533, 0.3562, 0.3592, 0.3624, 0.3656, 0.3688, 0.3721,
     3 0.3756, 0.3791, 0.3829, 0.3867, 0.3907, 0.3951, 0.3995, 0.4046,
     4 0.4098, 0.4162, 0.4236, 0.4335, 0.4632/
        call rndc(u)
        if(u .lt. .007) then
           f= (u*f1)**.4
        else
           l=u*100.+1
           f=(fn(l+1)-fn(l))*100.*(u-(l-1)/100.) + fn(l)
        endif
       end
       subroutine csampMuEKl3(f)
       implicit none
!----       include '../../../Zmass.h'
#include  "Zmass.h"
       real *8  f, p
!            muon energy sampling table
       integer jpa, i
          real*8 alfa, a2, gz, gzs, gz2, u, ff
          integer l
          real*8  tmp, pp
          save ff

          parameter (alfa=masmu/maskc,
     *    a2=alfa**2,
     *    gz=-.35, gzs=gz**2, gz2=2.*gz)

!
      real*8 fb(101)
      data (fb    (i),i=   1,  72)/
     1 0.2140, 0.2232, 0.2285, 0.2329, 0.2368, 0.2404, 0.2438, 0.2470,
     2 0.2500, 0.2529, 0.2557, 0.2584, 0.2610, 0.2635, 0.2660, 0.2684,
     3 0.2708, 0.2731, 0.2754, 0.2777, 0.2799, 0.2820, 0.2842, 0.2863,
     4 0.2884, 0.2905, 0.2925, 0.2945, 0.2965, 0.2985, 0.3005, 0.3024,
     5 0.3044, 0.3063, 0.3082, 0.3101, 0.3120, 0.3139, 0.3157, 0.3176,
     6 0.3195, 0.3213, 0.3232, 0.3250, 0.3268, 0.3287, 0.3305, 0.3323,
     7 0.3341, 0.3359, 0.3378, 0.3396, 0.3414, 0.3432, 0.3451, 0.3469,
     8 0.3487, 0.3505, 0.3524, 0.3542, 0.3561, 0.3579, 0.3598, 0.3617,
     9 0.3635, 0.3654, 0.3673, 0.3692, 0.3712, 0.3731, 0.3751, 0.3770/
      data (fb    (i),i=  73, 101)/
     1 0.3790, 0.3810, 0.3831, 0.3851, 0.3872, 0.3893, 0.3915, 0.3936,
     2 0.3959, 0.3981, 0.4004, 0.4027, 0.4051, 0.4076, 0.4101, 0.4127,
     3 0.4154, 0.4182, 0.4211, 0.4241, 0.4273, 0.4306, 0.4342, 0.4381,
     4 0.4425, 0.4474, 0.4533, 0.4611, 0.4861/

           call rndc(u)
           l=u*100.+1
           f = (fb(l+1)-fb(l))*100.*(u-(l-1)/100.) + fb(l)
           ff = f
          return
!      **************************
       entry cmuPolAtK(jpa, p)
!      **************************
!            this must be called after csampMuEKl3
!         jpa: -1 for k- and k0
!              +1 for k+ and k0bar
!          p: real*8. output
!            this is for k at rest; and for k- and k0
!
           tmp=ff**2-a2
           if(tmp .le. 0.) then
              pp=0.
           else
              pp=sqrt(tmp)* ( -4*(1.-2*ff)+ (gzs+gz2-3)*a2) /
     *       (4*ff*(1.-2*ff) + a2* (5*ff-a2+ gz*(4-6.*ff+2*a2) +
     *        gzs*(ff-a2) ) )
           endif
           if(jpa .lt. 0) then
              p=pp
           else
              p=-pp
           endif
       end
!          muon polarization in lab for k-->pi+mu+neu;
      subroutine cmuPolAtLabK(jpa, muon, kaon, p)
          implicit none
!----          include '../../../Zmass.h'
#include  "Zmass.h"
!----          include '../../../Zptcl.h'
#include  "Zptcl.h"
          type(ptcl):: muon, kaon
          integer jpa
          real*8 p
!             jpa: intege. input. -1 for k- k0 bar
!                          1  for k+ k0
!            muon: /ptcl/  input.  muon energy in lab is used
!            kaon: /ptcl/ input.  kaon energy in lab is used.
!              p:  real*8. output. longitudianl polarization of muon
!             k0==> mu+ neu + pi-
!             k+==> mu+ neu + pi0
!             k0bar==>mu-+neu+pi+
!             k-==>mu- + neu+pi0
!                   gzai=-.35 is assumed.
              real*8 pa( 7), pb(18)
              real*8 x
!                  table is for k0 or k+      see n.p vol22 p553
!               pa: p of lab for 40*muon.e/ek*mk=0.8 to 2.0 step .2
              data pa/-.74,-.74, -.735, -.72, -.70, -.64, -.61/
!               pb: same for  2.0 to 18 step 1
              data pb/-.61, -.40, -.22, -.08, 0.03, .12, .19, .263,.31,
     *               .365,.435, .48,.54, .58, .63, .675,.72,.78/
!


              x=muon%fm%p(4)/kaon%fm%p(4)*40.*maskc
              if(x .gt. 19.) then
!                     very high energy same as cms
                  call cmuPolAtK(jpa, p)
              elseif(x .lt. .8) then
!                     very low energy same as cms with oppsit sign
                  call cmuPolAtK(jpa, p)
                  p=-p
              elseif(x .lt. 2.) then
                  call kintp3(pa, 1,  7, 0.8d0, .2d0, x, p)
              else
                  call kintp3(pb, 1, 18, 2.d0, 1.d0, x, p)
              endif
              if(jpa .lt. 0) then
                   p=-p
              endif
       end
