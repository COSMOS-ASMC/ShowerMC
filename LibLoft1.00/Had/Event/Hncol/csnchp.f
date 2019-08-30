!      ******************************************************
       subroutine csnchp(icon)
!            sample # of charged ptcls
!      ******************************************************
!        Nch: integer. Output. sampled charged
!                              particle number (excludeing
!                              leadings)
!       icon: integer. Output. 0 --> o.k
!                              1 --> n.g (missing mass
!                                       is too small)
       implicit none

#include  "Zptcl.h"
#include  "Zmass.h"
#include  "Zevhnv.h"
!
        integer  icon
        real*8 redf
        real*8 missgm     ! missing mass
        real*8 roots, parmk

!
        roots =  Cmsp%mass
        missgm = Missingp%mass
        if(missgm  .gt. maspic*1.1) then
!             Efrs = missgm* 2.5 *(roots/200.)**0.05
!             Efrs = missgm* 2.5 *(roots/200.)**0.06 
              Efrs =max(
     *            missgm* 2.5 *(roots/200.)**0.06,
     *            missgm*2. + Pjtatr%mass + masp)

!                <nch> as a funcion of roots
            call ccpmul(Efrs, Avncharged)
!            call ccpmul(roots, Avncharged)
!                see p.l 116b p195; correction : NOW  OBSOLUTE. in this case
!                                           we must use Avncharged by roots.  
!                for available energy (reduction factor)
!            redf=1.5d0*sqrt(missgm/roots)   good at 53 gev
!     *        /(1.d0+ 95.d0/Pjlab.fm.p(4))**0.30d0 *.977d0  
!            redf=1.560d0*sqrt(missgm/roots)      ! good at 900 gev
!     *        /(1.d0+ 95.d0/Pjlab.fm.p(4))**0.30d0 
!               compromise  two above.
!            redf = 1.4385 * (Pjlab.fm.p(4)/1490)**0.0143 *
!     *             sqrt(missgm/roots)

!              <n> from ccpmul is for NSD.
!              For inclusion of SD events,
!              aven must be corrected at low energies
!            Avncharged = Avncharged* redf    ! effective <N>

!
!                fix n_b (=negative binormial) parameter k
            call cnbk(Efrs,  parmk)
!                sample n_charge
            if(parmk .le. 0.d0) then
              call cknoNarrow(Avncharged, Nch)
!               call kpoisn(Avncharged, Nch)
            elseif(Avncharged .lt. 8.) then
!                 normal distribution is too wide
               call cknoNarrow(Avncharged, Nch)
!               call kpoisn(Avncharged, Nch)
            else
!                    negative binomial
               call knbino(parmk, Avncharged, Nch)
            endif
            icon=0
        else
            Nch=0
            icon=1
        endif

        end
!       *************** simplest kno
      subroutine ckno(ave, sampled)
!
!          use z*exp(-pi/4 *x**2) dz
!  
      implicit none
      real*8  ave  ! input. average number
      integer sampled
!  
      real*8 u
      real*8  sqrtpi/1.772453851/  ! sqrt(pi)

      call rndc(u)
!           not add 0.5
      sampled = max(0.d0, sqrt( -log(u) )* 2.0/sqrtpi * ave )
      end
!       *************** narrow kno good at low energies.
      subroutine cknoNarrow(ave, sampled)
!
!          use z**2 exp(-0.73 *x**3) dz
!  
      implicit none
      real*8  ave  ! input. average number
      integer sampled
!  
      real*8 u
      call rndc(u)
      sampled = ( -log(u)/0.73)**0.333 * ave
      end
!     *****************************************************
      subroutine cnbk(roots,  ak)
        implicit none
        real*8 roots, ak
!             Negative binormial parameter k.(UA5 parameterization)
!             slog= log(s/gev**2); effective s
        real*8 slog
!
          slog = log(roots)*2
          if(slog .gt. 5.3d0)then
             ak= 1.d0/ (slog * 0.029d0 - 0.104d0)
          else
             ak=-1.d0
          endif
        end
!       ******************
        subroutine cfnptc(a,  ntot)
!       ******************
!          fix # of ptcls of each type, give mass, code
!
!  nch: integer.  Input. # of charged ptcls to be sampled.

!    a: /ptcl/     Output. to get ptcls. (mass, code, 
!                  charge) are assigned. some of subcode  is
!                  also assigned. (nn~, dd~ mass should be
!                  refixed later)
! ntot: integer.   Output. to get the toal # of ptcls to be
!                  produced.
!                 
!  *** Note ***
!       After this call, the # of particle of pi+-0, K+-0,
!       etc can be obtained as Npi0 etc which are in
!       ../Zevhnv.h
!       ----------------------------
        implicit none
!----        include '../../Zptcl.h'
#include  "Zptcl.h"
!----        include '../../Zcode.h'
#include  "Zcode.h"
!----        include '../../Zmass.h'
#include  "Zmass.h"
!----        include '../Zevhnp.h'
#include  "Zevhnp.h"
!----        include '../Zevhnv.h' 
#include  "Zevhnv.h"
        type(ptcl)::  a(*)
        integer  ntot
!
!
        real*8 missml, rnnb, rkc, p, exe, ddb
        integer nchc,  ntp, i
!

        missml = log( Missingp%mass ) * 2
!             get average fraction of n-n~ pair to Nch (non leading)
!             lamda decay product. (exclude lamda c)
        call cfrnnb(missml, rnnb)
!                get average fraction; (k+ + k-)/(pi+ + pi-)
!c        call cfrkc(missml, rkc)
        call cfrkc(log(Efrs)*2, rkc)
!                get average # of dd~ pair
!                this is to account for the prompt muon production
        call cnddb(Efrs, ddb)
        ddb=ddb * Mudirp      ! mudirp;  default is 1.0
!                 fix the # of particles of each type -------------
        if(Nch .eq. 0) then
                Nnnb = 0
                Nkaon = 0
                Nddb = 0
                Npic = 0
                Nk0 = 0
                Nkch = 0
                call kpoisn(Avncharged/2, Npi0)
!                call ckno(Avncharged/2, Npi0)
         else
                nchc=Nch
                p =rnnb*nchc
                call kpoisn(p, Nnnb)
!                call ckno(p, Nnnb)
                call kpoisn(ddb, Nddb)
!                call ckno(ddb, Nddb)
!                  the number of remaining charge(statistically)
                nchc = nchc-Nnnb - Nddb
!                   k+,k-,k0,k0~ (eqaul number in each type)
                p =rkc/(1.+rkc)*nchc
                call kpoisn(p , Nkch)
!                call ckno(p , Nkch)
                Nkaon = Nkch*2
                Nk0 = Nkaon- Nkch
                Npic = max(nchc- Nkch, 0)
                p = Nch*.51
                call kpoisn(p, Npi0)
!                call ckno(p, Npi0)
         endif
         if(Npi0 .gt. 10) then
!            assume some of them are eta.   the pi0/eta ratio is
!            a parameter. normally 0.2. which means Neta= 0.16*Npi0 
!           this is only to see the effect at >>10^19 ev region where
!           the decay of pi0 is inhibited and only eta can be the
!           source of h.e gamma.
            Neta =  Eta2Pi0 / (1 . + Eta2Pi0)  * Npi0
            Npi0 = Npi0 - Neta
         else
            Neta = 0
         endif
!             
         ntp=0      ! counter for storing ptlcs in a.
         do   i=1, Nnnb
!                     sample additional excitation mass of nn~
!                      <>=400 MeV
                  call ksgmim(1, 400.d-3, exe)
                  ntp=ntp+1
                  call cmkptc(knnb, 0, 0,  a(ntp)) 
                  a(ntp)%mass=exe + a(ntp)%mass
         enddo
!                (d,d~)
         do   i=1, Nddb
!                     sample additional excitation mass of dd~=
                  call ksgmim(1, 400.d-3, exe)
                  ntp=ntp+1
                  call cmkptc(kddb, 0, 0,  a(ntp))
                  a(ntp)%mass=exe + a(ntp)%mass
         enddo
!                 pi+/-
         do   i=1, Npic
                  ntp=ntp+1
!                     + or - is determined later (set tentative +)
                  call cmkptc(kpion, 0, 1,  a(ntp))
         enddo
!                 pi0
         do   i=1, Npi0
                  ntp=ntp+1
                  call cmkptc(kpion, 0, 0, a(ntp))
         enddo
         do  i = 1, Neta
            ntp = ntp + 1
            call cmkptc(keta, 0, 0, a(ntp))
         enddo
!                 kaon +/ -
         do   i=1, Nkch
                  ntp=ntp+1
!                     + or - is determined later (set tentative +)
                  call cmkptc(kkaon,  0, 1,  a(ntp))
         enddo
!                 k0
         do   i=1, Nk0
                  ntp=ntp+1
!                     k0,k0~ (long,short) is determined later
                  call cmkptc(kkaon, 0, 0, a(ntp))
         enddo
         ntot=ntp
        end
!     *****************************************************
        subroutine cfrnnb(efsl, rn)
!         fraction of (nn~) pairs; including lamda decay products
!         (but not lamda_c) to the total charged ptcls
!          efsl: log(s/gev**2). s is effective s. based on UA5
        implicit none
        real*8 efsl, rn
!           rn= 0.0115*efsl - 0.015  ; this is # of n + n~
           rn= 0.0057d0*efsl - 0.0075d0
        end
!     *****************************************************
        subroutine cfrkc(efsl, rk)
        implicit none
!----        include '../../Zptcl.h' 
#include  "Zptcl.h"
!----        include '../Zevhnp.h' 
#include  "Zevhnp.h"
!----        include '../Zevhnv.h'
#include  "Zevhnv.h"
!          fraction of k_charge to the pi_charge
!          efsl=log(s/gev**2). Cmsp.mass=root(s)
!          rk=0.07, 12.3, 14, 21 at 10, 100, 10000 GeV, 10**18 eV.

        real*8 efsl, rk, tmp
          tmp = Cmsp%mass**2-4.63
          if(tmp .le. 0.) then
             rk = 0.
          else
             rk=(Kpilog*(efsl+0.069) + Kpicns)*
     *        exp(-8.0/ tmp)
          endif
        end
!     ************************************************
        subroutine cnddb(efrs,  ddb)
!            average # of ddb pairs for p-p collisions
!          efrs: effective roots in GeV
!         *** for p-p or pi-p collision. the number of ddb
!          ---------old----------------
!            is assumed to be 1.e-3*log(roots)* exp(-78/roots)
!            at 1 TeV lab., this is about 6.2e-4
!            at 10000 TeV lab.,     8e-3
!         ---------after v3.0 ----------
!             is assumed to be 3.e-3*roots**0.25* exp(-78/roots)
!             at 1 TeV lab., this is about 1.2e-3
!             at 10000 TeV lab.         2.e-2
        implicit none

        real*8 efrs,  ddb
!
           ddb=3.d-3 * efrs ** 0.25 * exp(-78.d0/efrs)
        end










