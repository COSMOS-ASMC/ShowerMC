!    ******************************************************************
!    *                                                                *
!    *   ckaonDecay
!    *                                                                *
!    ******************************************************************
!
       subroutine ckaonDecay(pj, mupol, a,  np, polari)
       implicit none
!----       include '../../Zptcl.h'
#include  "Zptcl.h"
!----       include '../../Zcode.h'
#include  "Zcode.h"

       integer np               !output. no. of ptcls produced
       type(ptcl):: pj         ! input. kaon
       logical mupol            ! input. if T, muon polarization is considered
       type(ptcl):: a(*)      ! output. produced ptcls
       real*8  polari         ! output. polarization of the muon.  if a containes
!                               muon. muon is put in a(np), if any.
!
!
      if(pj%charge .ne. 0) then
!           k+-
         call ckChgDcy(pj, mupol, a, np, polari)
      elseif(abs(pj%subcode) .eq. k0s) then
!           k0 short
         call ckShortDecay(pj,mupol,  a, np, polari)
      else
!           k0 long
         call ckLongDecay(pj, mupol, a, np, polari)
      endif
      end
      subroutine ckChgDcy(pj, mupol,  a, np, polari)
!
!            k+- decay
!  -- process --
!         1) k---->mu+neu         (63.5 %)
!         2)  ---->pic + pi0       21
!         3)  ---->pi0+e+neu       4.8
!         4)  ---->pi0+mu+neu      3.2
!         5)  ---->pic+pic+pic     5.6
!         6)  ---->pic+pi0+pi0     1.7
      implicit none
!
!----      include '../../Zptcl.h'
#include  "Zptcl.h"
      
      integer np
      type(ptcl):: pj
      type(ptcl):: a(*)
      logical mupol
      real*8 polari
!
      real*8 u
      integer icp(2)/1, 0/, icm(2)/-1, 0/
      integer icp3(3)/1, 1, -1/
      integer icm3(3)/-1, 1, -1/
      integer ic3p0(3)/1, 0, 0/, ic3m0(3)/-1,0,0/
!
      call rndc(u)

      if(u  .lt. 0.635) then
!//////////
!         write(0,'(a, 1p,e14.4)') 'kch 1 ', pj.fm.p(4)
!////////////
!           k-->mu + neu
          call ckMuDecay(pj, mupol, a, np, polari)
      elseif(u .lt. .845) then
!           k+ --> pi+ + pi0  or  c.c
         if(pj%charge .gt. 0) then
            call ck2piDecay(pj, icp,  a, np)
         else
            call ck2piDecay(pj, icm,  a, np)
         endif
      elseif(u .lt. .893) then
!           k+ ---> pi0 +  e+ + neue
!//////////
!         write(0,'(a, 1p,e14.4)') 'kch 2 ', pj.fm.p(4)
!////////////
          call ckPiENeuDecay(pj, a, np)
      elseif(u .lt. .925) then
!           k+ --->mu+ + Neumu +  pi0 
!//////////
!         write(0,'(a, 1p,e14.4)') 'kch 3 ', pj.fm.p(4)
!////////////
         if(mupol) then
            call ckMuNeuPiDcy(pj,  a, np, polari)
         else
            call ckMuNeuPiDcy2(pj, a, np)
            polari = 0.
         endif
      elseif(u .lt. .981) then
!           k+ ---> 3*pic
         if(pj%charge .eq. 1) then
            call ck3PiDecay(pj, icp3, a, np)
         else
            call ck3PiDecay(pj, icm3, a, np)
         endif
      else
         if(pj%charge .eq. 1.) then
            call ck3PiDecay(pj,  ic3p0, a, np)
         else
            call ck3PiDecay(pj,  ic3m0, a, np)
         endif
      endif
      end
!     *******************************************
!            k0s decay
!  -- process --                             2014.sep
!         1)  ---->pi+ + pi-       68.61%-->69.20
!         2)  ---->pi0 + pi0       31.39 -->30.69
!         3)  ---->pic + e + neu_e          7.04e-4
!         4)  ---->pic + mu + neu_mu        4.69e-4
      subroutine ckShortDecay(pj, mupol, a, np, polari)
      implicit none
!----      include '../../Zptcl.h'
#include  "Zptcl.h"
#include  "Zevhnp.h"
      integer np
      type(ptcl):: pj, a(*)
      logical,intent(in):: mupol
      real(8),intent(out):: polari
      
      integer ic1(2)/1, -1/, ic2(2)/0, 0/
      real*8 u

      real(8),parameter::br(4)=(/69.20d0,30.69d0,7.04d-2,4.69d-2/)
      logical,save::first=.true.
      real(8),save:: cbr(4)
      integer::i
      if( first ) then
         cbr(:)=br(:)
         do i = 2, 4
            cbr(i) = cbr(i-1)+ cbr(i)
         enddo
         if( K0sSemiLD == 1) then
!           3)   pi+e +nue is considered
            cbr(:) = cbr(:)/cbr(3)
         elseif(K0sSemiLD == 2) then
!            4)  pi+mu +numu is considered
            cbr(3) = cbr(2)
            cbr(:) = cbr(:)/cbr(4)
         elseif( K0sSemiLD == 12) then
!            3),4) are considerd
            cbr(:) = cbr(:)/cbr(4) 
         elseif( K0sSemiLD == 0) then
!            3),4) are neglected
            cbr(:) = cbr(:)/cbr(2)
         else
            write(0,*) ' K0sSemiLD =', K0sSemiLD, ' invalid'
            stop
         endif
!///////////////
!         write(0,*) ' cbr =', cbr(:)
!/////////////
         first = .false.
      endif

      call rndc(u)
      do i = 1, 4
         if(u .lt. cbr(i)) exit
      enddo
!/////
!      write(0,'(a, i3, 1p, e14.4)') 'k0s  ', i, pj.fm.p(4)
!///////////
      polari = 0.
      if(i == 1 ) then
         call ck2piDecay(pj, ic1, a, np)
      elseif( i == 2 ) then
         call ck2piDecay(pj, ic2, a, np)
      elseif( i == 3 ) then
         call ckPiENeuDecay(pj, a, np)  ! use K0l case
      elseif( i == 4)  then  ! use k0l case
         if(mupol) then
            call ckMuNeuPiDcy(pj,  a, np, polari)
         else
            call ckMuNeuPiDcy2(pj, a, np)
            polari = 0.
         endif
      endif
      end
!            k0l decay
!  -- process --
!         1)  ---->e  pi neue   38.7
!         2)  ---->mu pi neum   27.1 %.    (k0==>mu+, k0bar==>mu-)
!         3)  ---->3 pi0        21.5
!         4)  ---->pi+ pi- pi0  12.4
!
      subroutine ckLongDecay(pj, mupol, a, np, polari)
      implicit none
!----      include '../../Zptcl.h'

#include  "Zptcl.h"
      logical mupol
      integer np
      type(ptcl):: pj       ! kaon
      type(ptcl):: a(*)      ! outputn
      real*8 polari            ! output

      real*8 u
      integer ic(3)/0, 0, 0/, ic2(3)/1, -1, 0/
!
      call rndc(u)

      if(u .lt. .387) then
!           e  + neue + pi
!//////////
!         write(0,'(a, 1p, e14.4)') 'k0l 1 ', pj.fm.p(4)
!///////////
          call ckPiENeuDecay(pj, a, np)
      elseif(u .lt. .658) then
!//////////
!         write(0,'(a, 1p, e14.4)') 'k0l 2 ',  pj.fm.p(4)
!//////////
!           mu + neumu + pi
          if(mupol) then
               call ckMuNeuPiDcy(pj,  a, np, polari)
          else
               call ckMuNeuPiDcy2(pj, a, np)
               polari = 0.
          endif
      elseif(u .lt. .873) then
!           3 pi0
         call ck3PiDecay(pj, ic, a, np)
      else
!           pi+ pi- pi0
         call ck3PiDecay(pj, ic2, a, np)
      endif
      end
!     ****************************************************
      subroutine ckMuNeuPiDcy(pj, a, np, polari)
!     ****************************************************
!       k->    mu + neum + pi (parent may be k charge,or k0)
!             
      implicit none
!----      include '../../Zptcl.h'
#include  "Zptcl.h"
!----      include '../../Zcode.h'
#include  "Zcode.h"

      integer np
      type(ptcl):: pj       ! kaon
      type(ptcl):: a(*)      ! outputn
      real*8 polari            ! output

      real*8 u, ecm, cosa, f, pcm
      integer jpa, i
      logical ok
!            make ptcl
!             k+==>pi0  mu+  nue(m)    k- ==> pi0 mu- nue_b(m)
!             k0l==> pi- mu+ nue(m) or  => pi+ mu- nue_b(m)


      if(pj%charge .eq. -1.) then
         call cmkptc(kneumu, antip, 0, a(1))
         call cmkptc(kpion, 0, 0, a(2))
         call cmkptc(kmuon, 0, -1, a(3))
         jpa = -1
      elseif(pj%charge .eq. 1) then
         call cmkptc(kneumu, regptcl, 0, a(1))
         call cmkptc(kpion, 0, 0, a(2))
         call cmkptc(kmuon, 0, 1, a(3))
         jpa = 1
      else
!           k0l; anti or regptcl
         call rndc(u)
         if(u .lt. 0.5) then
            call cmkptc(kneumu, antip, 0, a(1))
            call cmkptc(kpion, 0, 1, a(2))
            call cmkptc(kmuon, 0, -1, a(3))
            jpa = 1
         else
            call cmkptc(kneumu, regptcl, 0, a(1))
            call cmkptc(kpion, 0, -1, a(2))
            call cmkptc(kmuon, 0, 1, a(3))
            jpa = -1
         endif
      endif

      ok = .false.
      do while (.not. ok) 
!              sample energy of neum at rest of k
         call csampNeuEKl3(f)

         ecm=f*pj%mass
!           angle
         call rndc(u)
         cosa=2*u-1.
!          set px,py,pz
         call cpCos2pxyz(cosa, ecm, a(1)%fm)
         a(1)%fm%p(4) = ecm

!           muon ; should be put in the last place in a
         np =3
         call csampMuEKl3(f)
         ecm=max(f*pj%mass, a(np)%mass*1.0001d0)
         pcm=sqrt(ecm**2- a(np)%mass**2)
!           angle
         call rndc(u)
         cosa=2*u-1.
         call cpCos2pxyz(cosa, pcm, a(np)%fm)
         a(np)%fm%p(4) = ecm
         a(2)%fm%p(4) = pj%mass- ecm - a(1)%fm%p(4)
         if( a(2)%fm%p(4) .gt. a(2)%mass) then
!              simply to satisfy the conservation
            a(2)%fm%p(1) = -a(1)%fm%p(1) - a(np)%fm%p(1)
            a(2)%fm%p(2) = -a(1)%fm%p(2) - a(np)%fm%p(2)
            a(2)%fm%p(3) = -a(1)%fm%p(3) - a(np)%fm%p(3)
            ok = .true.
         endif
      enddo

      do i = 1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
!            set muon polarization
      call  cmuPolAtLabK(jpa, a(np), pj, polari)


      end
!     ******************************
      subroutine ckMuNeuPiDcy2(pj,  a, np)
!     ******************************
!        k->   mu + neummu + pi (parent may be k charge,or k0)
!           all is considered but not polarization
      implicit none
!----      include '../../Zptcl.h'
#include  "Zptcl.h"
!----      include '../../Zcode.h'
#include  "Zcode.h"

      integer np
      type(ptcl):: pj       ! kaon
      type(ptcl):: a(*)      ! outputn

      real*8  w, u
      integer i, icon

!           make 3 ptcls
!             k+==>pi0  mu+  nue(m)    k- ==> pi0 mu- nue_b(m)
!             k0l==> pi- mu+ nue(m) or => pi+ mu- nue_b(m)
         if(pj%charge .eq. 1) then
            call cmkptc(kneumu, regptcl, 0, a(1))
            call cmkptc(kpion, 0, 0, a(2))
            call cmkptc(kmuon, 0, 1, a(3))
         elseif(pj%charge .eq. -1) then
            call cmkptc(kneumu, antip, 0, a(1))
            call cmkptc(kpion, 0, 0, a(2))
            call cmkptc(kmuon, 0, -1, a(3))
         else
            call rndc(u)
            if(u .lt. 0.5) then
               call cmkptc(kneumu, antip, 0, a(1))
               call cmkptc(kpion, 0, 1, a(2))
               call cmkptc(kmuon, 0, -1, a(3))
            else
               call cmkptc(kneumu, regptcl, 0, a(1))
               call cmkptc(kpion, 0, -1, a(2))
               call cmkptc(kmuon, 0, 1, a(3))
            endif
         endif
         np=3
!           3  body pure phase space
         call cnbdcy(3, pj%mass, a, 0, w, icon)
!              boost to lab
         do i = 1, np
            call cibst1(i, pj, a(i), a(i))
         enddo
      end
      subroutine ckPiENeuDecay(pj,  a, np )
!           e  + neue + pi  ; 
!              
!          k+==>pi0 e+ neue;   k-==>pi0 e- neue_b
!          k0l==>pi+ e- +neue_b or ==>pi- e+ neue
      implicit none
!----      include '../../Zptcl.h'
#include  "Zptcl.h"
!----      include '../../Zcode.h'
#include  "Zcode.h"
      
      integer np
      type(ptcl):: pj, a(*)
      type(ptcl)::piesys
      real*8 u, f, ecm, cosa
      integer i,  echg
      logical ok


      ok = .false.
      echg = pj%charge
      if(pj%charge .eq. -1.) then
         call cmkptc(kneue, antip, 0, a(1))
         call cmkptc(kelec, 0, echg, a(2))
         call cmkptc(kpion, 0,  0, a(3))
      elseif(pj%charge .eq. 1) then
         call cmkptc(kneue, regptcl, 0, a(1))
         call cmkptc(kelec, 0 , echg, a(2))
         call cmkptc(kpion, 0,  0, a(3))
      else
         call rndc(u)
         if(u .lt. 0.5) then
            call cmkptc(kneue, regptcl, 0, a(1))
            call cmkptc(kelec,  0,  1,  a(2))
            call cmkptc(kpion,  0,  -1, a(3))
         else
            call cmkptc(kneue, antip, 0, a(1))
            call cmkptc(kelec,  0,  -1,  a(2))
            call cmkptc(kpion,  0,  1, a(3))
         endif
      endif         

      do while ( .not. ok)
!              sample energy of neue at rest of k
         call csampNeuEKl3(f)
         ecm=f*pj%mass
!           angle
         call rndc(u)
         cosa=2*u-1.
!          set px,py,pz
         a(1)%fm%p(4) = ecm
         call cpCos2pxyz(cosa, ecm, a(1)%fm)

!         make pi + e  system 
         piesys%fm%p(4) = pj%mass - ecm
         if( piesys%fm%p(4) .gt. a(2)%mass+ a(3)%mass) then
            piesys%fm%p(1) = -a(1)%fm%p(1)
            piesys%fm%p(2) = -a(1)%fm%p(2)
            piesys%fm%p(3) = -a(1)%fm%p(3)
            piesys%mass = piesys%fm%p(4)**2 -(
     *        piesys%fm%p(1)**2+ piesys%fm%p(2)**2+piesys%fm%p(3)**2 )
            if(piesys%mass .gt. 0.) then
               piesys%mass = sqrt(piesys%mass)
               call c2bdcy(piesys, a(2), a(3))
               ok = .true.
            endif
         endif
      enddo
      np=3
!              boost to lab
      do i = 1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      end
!    ******************************************************************
!    *   ckMuDecay:   k -> mu + neumu
!    ******************************************************************
!
!
!     decay of k ---> mu + neu ( b.r=.635 ) is treated.
!
      subroutine ckMuDecay(pj, mupol,  a, np, polari)
      implicit none
!----      include '../../Zptcl.h'
#include  "Zptcl.h"
!----      include '../../Zcode.h'
#include  "Zcode.h"
      
      integer np
      type(ptcl):: pj  ! input. kaon
      type(ptcl):: a(*)  ! output. ptcls produced
      logical mupol     ! input.  T==>  muon polarization taken into acc.
      real*8 polari     ! output. muon polarizaton.
!                         muon must be put a(np).
      integer charge, subcode
!
!
      
!             make muon neutrino : muon set last
      subcode =  -pj%charge
      call cmkptc(kneumu, subcode, 0, a(1))
!             make muon
      charge = pj%charge
      call cmkptc(kmuon, 0, charge, a(2))
!
!           k-->mu + neu
      call c2bdcy(pj, a(1), a(2))
!               set polarization of muon
      if(mupol) then
         call ckmuPolari(pj, a(2),  polari)
      else
         polari=0.
      endif
      np=2
      end
!     **************************************
      subroutine ck2piDecay(pj, ic,  a, np )
!     **************************************
!            k--> 2 pi (k+-->pi+ pi0 or k- --> pi-  pi0)
!                      (k0-->pi+ pi- or k0 --> pi0 +pi0)
      implicit none
!----      include '../../Zptcl.h'
#include  "Zptcl.h"
!----      include '../../Zcode.h'
#include  "Zcode.h"
      
      integer np, ic(2)
      type(ptcl):: pj, a(*)
      call cmkptc(kpion, 0, ic(1), a(1))
      call cmkptc(kpion, 0, ic(2), a(2))

      call c2bdcy(pj, a(1), a(2))
      np=2
      end
!     **************************************
      subroutine ck3PiDecay(pj, ic, a, np )
!     **************************************
!          k--> 3 pi;

      implicit none
!----      include '../../Zptcl.h'
#include  "Zptcl.h"
!----      include '../../Zcode.h'
#include  "Zcode.h"
      
      integer np, ic(3)
      type(ptcl):: pj, a(*)
      
      real*8 w
      integer icon, i
!
      call cmkptc(kpion, 0, ic(1), a(1))
      call cmkptc(kpion, 0, ic(2), a(2))
      call cmkptc(kpion, 0, ic(3), a(3))
      call cnbdcy(3, pj%mass, a, 0, w, icon)
      np = 3
      do   i=1, np
         call cibst1(i, pj, a(i), a(i))
      enddo
      end
