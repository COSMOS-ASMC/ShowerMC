!             hadron A collision by Lund
        subroutine chALund(pj, ia, iz,  a, np)
        implicit none

#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zevhnv.h"
!
        type (ptcl):: pj, a(*)
        integer ia,iz, np
!             to assure the inclusion of block data
         external cblkfritio
!
!         lund or nucrin branch point: total  energy

         real*8 eintmn/0.2/
!
         real*8 elab, ek, u, alpha

         elab=pj%fm%p(4)
         ek = elab - pj%mass
!         if(elab .gt. 5. ) then
         if(ActiveMdl .eq. 'fritiof1.6'
     *    .or. ( ActiveMdl .eq. 'byenergy' .and.
     *           ek .gt. 4.0 ) ) then
!                lund monte carlo ; fritiof
            call cfritiof(pj, ia, iz,  a, np)
!              in  above, momentum  is measured from pj.fm so 
!              we need rotate it.
!       ////////
!cc            call ccheckMom('after fritio: before rot', a, np)
!      /////////            
            call crot3mom(pj, a, np)
         elseif(ek  .gt. eintmn .and.
     *         ( ActiveMdl .eq. 'nucrin' .or.
     *          ActiveMdl .eq. 'byenergy' ) ) then
            call cnucrin(pj, ia, iz,  a, np)
            call crot3mom(pj, a, np)
         elseif(pj%code .eq. knuc .and. pj%subcode .eq. antip
     *      .and. (ActiveMdl .eq. 'nucrin'  .or.
     *             ActiveMdl .eq. 'byenergy' ) ) then
!              for anti-barion, no threshold.  completely 0
!              K.E must be avoided for nucrin
            if( pj%fm%p(4) - pj%mass .le. 0.0 ) then
               pj%fm%p(4) = pj%mass*1.0001
            endif
            call cnucrin(pj, ia, iz,  a, np)
            call crot3mom(pj, a, np)
         else
!              move incident to a:
            np=1
            a(1) = pj
!              reduce energy (evaporation etc)
            alpha =( (ia-1.0)/(ia+1.0) )**2
!              Elastic scattering: E1--> E2
!             E2/E1 = (A**2 + 2Acos + 1)/(A+1)**2 = 1/2 ((alpha+1) + (1-alpha)*cos)
!            assume isotroic; then E2  is uniform in (Alpha, 1)E1
!              
            call rndc(u)
            a(1)%fm%p(4) = a(1)%mass + ek*((1.0-alpha)*u+alpha)
            call cadjm(a(1), a(1))
         endif
       end
!      ////////////
       subroutine ccheckMom(com, a, np)
        implicit none

#include  "Zptcl.h"
#include  "Zcode.h"
!
        character*(*) com
        integer np
        type (ptcl):: a(np)
        real*8 temp
        integer i
        do i =1, np
           temp= sqrt(a(i)%fm%p(4)**2  - a(i)%mass**2)
           temp = (a(i)%fm%p(1)/temp)**2
     *         +  (a(i)%fm%p(2)/temp)**2
     *         +  (a(i)%fm%p(3)/temp)**2
           if(abs(temp -1.) .gt. 0.1 .or. 
     *             a(i)%code .eq. 2) then
              write(*, *) com, temp, a(i)%code, a(i)%fm%p(4)
           endif
        enddo
        end
!      ///////////////////
!         *****************************************************
!         *
!         * cfritiof: lund monte carlo at Ek>4 GeV  :fritiof
!         *
!         *****************************************************
!
       subroutine cfritiof(pj, ia, iz,  a, np)
       implicit none
#include  "ZstrangeCode.h"
#include  "Zptcl.h"
#include  "Zcode.h"
       integer ia, iz, np
       type (ptcl):: pj, a(*)
!
!
      real*4 elab, rots, r0p, r0t, bmin, bmax
      integer nap,  nat, nzt, iflspv, neve, iproty, ifermi,
     *          iflout, nzp
      common /indataC/ elab,rots,nap,nzp,r0p,nat,nzt,r0t,iflspv,bmin,
     *                bmax,neve,iproty,ifermi,iflout

      integer*4 nevent, isppp, isppn, isptp, isptn, idi, ipr
      real*4 bimp
      common /evevecC/ nevent,isppp,isppn,isptp,isptn,bimp,
     *                idi(2000),ipr(2000)

      real*4 p
      integer n, k
      common /lujetsC/ n,k(2000,2),p(2000,5)

      real*4 dpar, cbr
      integer idb, kdp
      common/ludat3C/dpar(20),idb(120),cbr(400),kdp(1600)

      integer mst
      real*4 par
      character*150   msg
! ---------------------------------------------------------------------
!              Debug-def    discard-def  notify-def use-strg-flg  halt
!     debug       o             x           o           o          x
!     discard     x             o           x           o          x
!     other       x             x           o           x          o
! ---------------------------------------------------------------------
!
#ifdef USE_STRANGEFLAG
      logical strangeCode
#endif

#ifdef NOTIFY_STRANGECODE
      integer irnow(2)  ! to store random num. seed for this event
#endif

      common/ludat1C/mst(40),par(80)
!  -----------------------------------------------
      integer kidn, icg, l, kk, ktmp, ctmp, pap
      real*8 x, et
!   ---------------------------------------------

!
!                                               r0 when used as target
!                                                   together with
!                    na    nz      r0             meson projectile
!     neutron         1     0    1.0    n               1.0
!     proton          1     1    1.0    p               1.0
!     deuterum        2     1     .69   d
!     helium          4     2     .81   he
!     bor            11     5     .62   b
!     carbon         12     6     .64   c                .54
!     oxygen         16     8     .72   o
!     aluminum       27    13     .86   al               .77
!     silicon        28    14     .87   si
!     argon          40    18     .94   ar
!     calcium        40    20     .94   ca
!     copper         64    29    1.01   cu               .92
!     silver        108    47    1.06   ag              1.00
!     xenon         131    54    1.08   xe
!     gold          197    79    1.12   au
!     lead          207    82    1.12   pb
!     uranium       238    92    1.13   u               1.07
!-----------------------------------------------------------------------
!$$$$$$$$ at a>=11 these can be approximated by
!         r0=-0.3664+1.196x-0.2396x**2 (x=log10(a))
!         r0=-0.4515+1.158x-0.2184x**2 (meson projectile)
!$$$$$$$$$
!-----------------------------------------------------------------------
! the following parameters should be used for projectiles only
!-----------------------------------------------------------------------
!     antiproton      1    -1    1.0    p-bar
!     pion (positive) 1     1    1.0    pi+
!     pion (negative) 1    -1    1.0    pi-
!     kaon (positive) 1     1    1.0    k+
!     kaon (negative) 1    -1    1.0    k-
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! ingebo creates an event and fills the commonblocks evevec and lujetsaa
!-----------------------------------------------------------------------
! evevec : nevent => event number (consecutive numbering)
!          isppp  => number of non-interacting protons associated with
!                    the projectile
!          isppn  => number of non-interacting neutrons associated with
!                    the projectile
!          isptp  => number of non-interacting protons associated with
!                    the target
!          isptn  => number of non-interacting neutrons associated with
!                    the target
!          bimp   => impact parameter given in fermi
!          idi    => =1 if a diffractively produced particle
!                    =0 otherwise
!          ipr    => =1 if the produced particle comes from a
!                       projectile string.
!                    =0 if it comes from a target string.
!          note!     the validity of idi and ipr are not kept if calls
!                    to luedit are made
!
!-----------------------------------------------------------------------
! lujets : n      => number of entries in the lujets arrays
!          p      => array containing the three momentum components
!                    (x, y and z, where z is the component along the
!                    beam-axis), total relativistic energy, and the
!                    invariant mass of the particle
!          k      => k(j,1) contains the history of the j:th particle
!                    k(j,2) contains the kf flavour code of the j:th
!                    particle
!-----------------------------------------------------------------------
!

!  
!          make pi0 and k0s stable
         idb(23)=0
         idb(37)=0
         mst(19)=0
         kidn=pj%code
         icg=pj%charge
!              
         elab=pj%fm%p(4)
!       ---  fix target: ------------------
         nat=ia
         nzt=iz
         if(ia .eq. 14) then
             if(kidn .ne. knuc) then
                 r0t=.588
             else
                 r0t=.689
             endif
         elseif(ia .eq. 16) then
!                oxigen
             if(kidn .ne. knuc) then
                 r0t=.658
             else
                 r0t=.72
             endif
         elseif(ia .eq. 1)then
!               nucleon (comes here if used in our model)
             r0t=1.
         elseif(ia .eq. 40) then
!               argon
             if(kidn .ne. knuc) then
                 r0t=.843
             else
                 r0t=.94
             endif
         elseif(ia .ge. 11) then
             x=log10(float(ia))
             if(kidn .eq. knuc) then
                r0t=-0.3664+1.196*x-0.2396*x**2
             else
                r0t=-0.4515+1.158*x-0.2184*x**2
             endif
         elseif(ia .eq. 4) then
             if(kidn .eq. knuc) then
                r0t=.81
             else
                r0t=.65
             endif
         elseif(ia .ge. 5) then
            r0t = (5-ia)*0.04   + .8
         else
             write(msg,*) ' target mass=',ia, ' for cfritiof invalid'
             call cerrorMsg(msg, 0)
         endif
!        --- fix projectile --------------
         if(kidn .eq. knuc) then
             if(icg .eq. 1) then
                 iproty=41
             elseif(icg .eq. 0) then
                 if(pj%subcode .eq. regptcl) then
                     iproty=42
                 else
                     iproty=-42
                 endif
             else
                 iproty=-41
             endif
          elseif(kidn .eq. kpion) then
             if(icg .eq. 1) then
                 iproty=17
             elseif(icg .eq. -1) then
                 iproty=-17
             else
                 iproty=23
             endif
          elseif(kidn .eq. kkaon) then
             if(icg .eq. 1) then
                 iproty=18
             elseif(icg .eq. -1)then
                 iproty=-18
             elseif(abs(pj%subcode) .eq. k0l) then
!                  k0l
                 iproty=sign(38, pj%subcode*1)
             elseif(abs(pj%subcode) .eq. k0s) then
!                  k0s
                 iproty=sign(37, pj%subcode*1)
             endif
          else
             write(msg,*) ' ptcl code=',kidn,' not acceptable in ',
     *                  ' cfritiof; event neglectd '
             call cerrorMsg(msg, 1)
             np = 0
             return  !  exit
          endif
#ifdef USE_STRANGEFLAG
!           skip particles of undefined code; 
          do while (.true.)
             strangeCode = .false.
#endif
#ifdef NOTIFY_STRANGECODE
          call rnd1s(irnow)  ! save current seed for strange code
#endif
!                generate ptcls
          call ingeboC
!                make decay of unstable ptcls
          call lueditC(1)
!              move patcl into a
          np=0
          do   l=1, n
              kk=k(l,2)
              et=p(l, 4)
!                  ------------ pi --------------
              if(kk .eq. 17) then
                  ktmp=kpion
                  ctmp=1
                  pap=regptcl
              elseif(kk .eq. -17) then
                  ktmp=kpion
                  ctmp=-1
                  pap=antip
              elseif(kk .eq.  23) then
                  ktmp=kpion
                  ctmp=0
                  pap=0
              elseif(kk .eq. 18) then
!              -------   k charge ----------
                  ktmp=kkaon
                  ctmp=1
                  pap=regptcl
              elseif(kk .eq. -18) then
                  ktmp=kkaon
                  ctmp=-1
                  pap=antip
              elseif(kk .eq. 1) then
!              ----------------- gamma
!                  gamma; say from eta
                  ktmp=kphoton
                  ctmp=0
                  pap=kdirectg
!                  --------------  nucleon---
              elseif(kk .eq. 41) then
!                   proton
                  ktmp=knuc
                  ctmp=1
                  pap=regptcl
              elseif(kk .eq. -41) then
!                    anti proton
                  ktmp=knuc
                  ctmp=-1
                  pap=antip
              elseif(kk .eq. 42) then
!                   nutron
                  ktmp=knuc
                  ctmp=0
                  pap=regptcl
              elseif(kk .eq. -42) then
!                    n_b
                  ktmp=knuc
                  ctmp=0
                  pap=antip
!                -------------------  kaon  0
              elseif(abs(kk) .eq. 37) then
!                     kos
                  ktmp=kkaon
                  ctmp=0
                  pap=sign(k0s, kk)
              elseif(abs(kk) .eq. 38) then
!                    k0l
                  ktmp=kkaon
                  ctmp=0
                  pap=sign(k0l, kk)
              elseif(kk .eq.  7) then
!                    e-
                  ktmp=kelec
                  ctmp=-1
                  pap=regptcl
              elseif(kk .eq. -7) then
!                    e+
                  ktmp=kelec
                  ctmp= 1
                  pap= antip
              elseif(kk .eq.  9) then
!                   mu -
                  ktmp=kmuon
                  ctmp=-1
                  pap=regptcl
              elseif(kk .eq. -9) then
!                   mu +
                  ktmp=kmuon
                  ctmp= 1
                  pap= antip
              else
#ifdef NOTIFY_STRANGECODE
                  write(msg,*) ' ptcl code=',kk, ' from ingebo',
     *            ' not acceptable in cfritiof'
                  call cerrorMsg(msg, 1)
                  write(msg, *) ' input cosmos ptcl code=',
     *             pj%code, ' subcode=', pj%subcode, ' chg=',
     *             pj%charge, ' target A,Z=',ia, iz
                  call cerrorMsg(msg, 1)
                  write(msg, *)  ' E=', pj%fm%p(4), ' mass=',pj%mass,
     *            ' seed of the event= ', irnow
                  call cerrorMsg(msg, 1)
                  write(msg, *) ' its E=', et,' it is ', l, '-th ',
     *            ' among ',n, ' ptlcs'
                  call cerrorMsg(msg, 1)
#ifdef DEBUG_STRANGECODE
                  msg = ' The event is discarded'
                  call cerrorMsg(msg, 1)
#endif
#endif
#ifdef USE_STRANGEFLAG
                  strangeCode = .true.
#else
                  call cerrorMsg('execution halted', 0)
#endif
              endif
!              
              np=np+1
              call cmkptc(ktmp, pap, ctmp, a(np))
              a(np)%fm%p(4) = et
              a(np)%fm%p(1) = p(l,1)
              a(np)%fm%p(2) = p(l,2)
              a(np)%fm%p(3) = p(l,3)
           enddo
#ifdef  USE_STRANGEFLAG
           if(.not. strangeCode) goto 100
           enddo
 100       continue
#endif

        end
        block data cblkfritio
         common /indataC/ elab,rots,nap,nzp,r0p,nat,nzt,r0t,iflspv,bmin,
     *                bmax,neve,iproty,ifermi,iflout
        data r0p/1./
        data iflspv,bmin,bmax/0,0.,0./
        data rots,neve/0.,1/
        data ifermi/0/
        data iflout/-1/
        data nap/1/, nzp/1/
       end
!      **********************************************************
!      *
!      *  cnucrin:  nucrin for energy less than 4 GeV (k.e)
!      *
!      **********************************************************
       subroutine cnucrin(pj, ia, iz, a, np)
       implicit none
!----       include '../../Zptcl.h'
#include  "Zptcl.h"
!----       include '../../Zcode.h'
#include  "Zcode.h"
       integer ia, iz, np
       type (ptcl):: pj, a(*)

!       real*4 am, ga, tau
!       integer*2 ich,ibar,k1,k2
!       common/abltisC/am(110),ga(110),tau(110),ich(110),
!     * ibar(110),k1(110),k2(110)
       
       real*4 cxr, cyr, czr, elr, plr, tv
       integer ir, itr
       common/finucC/ir,itr(60),cxr(60),cyr(60),czr(60),elr(60),plr(60),
     *  tv
   
!
       real*4 elab, plab, anuc, znuc
       integer kidn, icg, iproty, itta, l, kk, ktmp, ctmp, pap
       real*8 et, u
       character*100  msg
!
!
         elab=pj%fm%p(4)
         plab=sqrt(pj%fm%p(4)**2 - pj%mass**2)
         kidn=pj%code
         icg=pj%charge
!       ---  fix target: ------------------
         anuc=ia
         znuc=iz
!        --- fix projectile --------------
         if(kidn .eq. knuc) then
             if(icg .eq. 1) then
                 iproty=1
             elseif(icg .eq. 0) then
                 if(pj%subcode .eq. regptcl) then
                     iproty=8
                 else
                     iproty=9
                 endif
             else
                 iproty=2
             endif
          elseif(kidn .eq. kpion) then
             if(icg .eq. 1) then
                 iproty=13
             elseif(icg .eq. -1) then
                 iproty=14
             else
                 iproty=23
             endif
          elseif(kidn .eq. kkaon) then
             if(icg .eq. 1) then
                 iproty=15
             elseif(icg .eq. -1)then
                 iproty=16
             else
                 iproty=24
             endif
          else
             write(msg,*) ' ptcl code=',kidn,' not acceptable in ',
     *        ' cnucrin; envet neglected '
             call cerrorMsg(msg, 1)
             np = 0
             return  ! ***********
          endif
!                generate ptcls
          if(ia .gt. 1) then
              call nucrinC(iproty, elab, 0., 0., 1., anuc, znuc)
          else
!                if target is nucleon, use hadrin below(high speed)
              if(iz .eq. 1) then
                  itta=1
              else
                  itta=8
              endif
              call chadrin(iproty, plab, elab,  itta)
          endif
!              move patcl into a
          np=0
          do   l=1, ir
              kk=itr(l)
              et=elr(l)
!                  ------------ pi --------------
              if(kk .eq. 13) then
                  ktmp=kpion
                  ctmp=1
                  pap=regptcl
              elseif(kk .eq. 14) then
                  ktmp=kpion
                  ctmp=-1
                  pap=antip
              elseif(kk .eq. 23) then
                  ktmp=kpion
                  ctmp=0
                  pap=0
              elseif(kk .eq. 15) then
                  ktmp=kkaon
                  ctmp=1
                  pap=regptcl
              elseif(kk .eq.  16) then
                  ktmp=kkaon
                  ctmp=-1
                  pap=antip
              elseif(kk .eq. 7) then
!              ----------------- gamma
!                  gamma; say from eta
                  ktmp=kphoton
                  ctmp=0
                  pap=kdirectg
!                  --------------  nucleon---
              elseif(kk .eq. 1 ) then
!                   proton
                  ktmp=knuc
                  ctmp=1
                  pap=regptcl
              elseif(kk .eq.  2 ) then
!                    anti proton
                  ktmp=knuc
                  ctmp=-1
                  pap=antip
              elseif(kk .eq. 8 ) then
!                   nutron
                  ktmp=knuc
                  ctmp=0
                  pap=regptcl
              elseif(kk .eq.  9 ) then
!                    n_b
                  ktmp=knuc
                  ctmp=0.
                  pap=antip
!                -------------------  k0
              elseif(kk .eq. 24) then
                  ktmp=kkaon
                  ctmp =0
                  call rndc(u)
                  if(u .lt. 0.5) then
                     pap = k0s
                  else
                     pap = k0l
                  endif
               elseif(kk .eq. 25) then
!                     k0b
                  ktmp=kkaon
                  ctmp =0
                  call rndc(u)
                  if(u .lt. 0.5) then
                     pap = -k0s
                  else
                     pap = -k0l
                  endif
              elseif(kk .eq. 3)  then
!                    e-
                  ktmp = kelec
                  ctmp = -1
                  pap = regptcl
              elseif(kk .eq.  4)  then
!                    e+
                  ktmp=kelec
                  ctmp=1
                  pap=antip
              elseif(kk .eq. 10)  then
!                    mu+
                  ktmp=kmuon
                  ctmp=1
                  pap=antip
              elseif(kk .eq. 11)  then
!                    mu-
                  ktmp=kmuon
                  ctmp=-1
                  pap=regptcl
              else
                 write(msg,*) ' ptcl code=',kk,' from nucrin',
     *            ' not acceptable in cnucrin'
                  call cerrorMsg(msg, 1)
                  call cerrorMsg(' it is neglected' ,1)
                  goto 10    ! *************
              endif
              np=np+1
              call cmkptc(ktmp, pap, ctmp, a(np))
              a(np)%fm%p(4) = et
              a(np)%fm%p(1) = plr(l)*cxr(l)
              a(np)%fm%p(2) = plr(l)*cyr(l)
              a(np)%fm%p(3) = plr(l)*czr(l)
!                  sometimes nucrin give unbalanced energy mom. at low E
              if( abs( 
     *        a(np)%fm%p(1)**2 + a(np)%fm%p(2)**2 + a(np)%fm%p(3)**2
     *        + a(np)%mass**2
     *        - a(np)%fm%p(4)**2)/a(np)%fm%p(4)  .gt. .005) then
                 call cadjm(a(np), a(np))
              endif
 10           continue
           enddo
        end
        subroutine chadrin(iproty, plab, elab, itta)
        real*4 plab, elab
!              bypass
!*** finuc: nucrin final state particle list with kinem. variables
!*** finlsp hadrin final state particle list with kinem. variables
!*** (number of particles,particle type index,direction cosines,energy
!*** absolute momentum and in finuc in addit. excitation energy
!--------------------------------------------------
        common    /finucC/irn,itrn(60),cxrn(60),cyrn(60),czrn(60),
     *       elr(60),plr(60),tv
        common    /finlspC/ir,itr(20),cxr(20),cyr(20),czr(20),
     *      el(20),pl(20)

          call hadrinC(iproty, plab, elab, 0., 0., 1., itta)

!             move /finlsp/ to /finuc/
           do   i=1, ir
             itrn(i)=itr(i)
             cxrn(i)=cxr(i)
             cyrn(i)=cyr(i)
             czrn(i)=czr(i)
             elr(i)=el(i)
             plr(i)=pl(i)
           enddo
          irn=ir
        end
