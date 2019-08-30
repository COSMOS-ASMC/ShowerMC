      subroutine cAAXsec2(pj, tgA, tgZ,  xs)
!     compute AA' inelastice cross -section 
!           At low energies, Shen's formulation is used  as is used
!         in PHITS. Over 5 GeV/n, Frier et al's formula is used
!         but it is normalized to the Shen's value at 5 GeV/n.
!
!
!     formula: xs = pi(r1 + r2 - d)**2
!        r = r0 A**(1/3)
!        d = 1.189 r0 exp(-0.055min(A1,A2))
!       r0 = (1.29f to 1.41f)
!  Therefore
!      xs = pi r0^2 (A1^0.333 + A2^0.333 - 1.189exp(-0.055min(A1,A2)))**2
!  This is by Frier et al in ICRC Paris conf. (from Uchu Hosha Sen Edited
!      by Nishimura, p.170)
!
!  At high energies (sqrt(s) > 80 GeV for E/A ?; we use Ek=50 GeV), we include energy dependence
!  of cross-section as follows.
!    Let the pp cross-section increases as E**delta, then
!    AB crosssection is well fitted by the dependence of E**alfa with
!    
!     alfa = 2.5* delta/(p + t + p*t/2)
!
!   where  p = A**(1/3) and t = B**(1/3)
!
!  
!       

      implicit none
#include "Zcode.h"
#include "Zmass.h"
#include "Zptcl.h"
#include "Ztrackp.h"
      type(ptcl):: pj    ! input projectile
      real(8),intent(in):: tgA  ! target mass #
      real(8),intent(in):: tgZ  ! target charge #
      real(8),intent(out):: xs  ! inela xs in mb

      real(8)::  ek
      real(8):: pjA, pjZ
!
      real(8):: p, t

      real(8):: Ekt, Ekpn, Etot, EktMeV

      real(8),parameter:: Enorm=5. ! KE/n. GeV.  normalize to shen at this E
      
      real(8):: elxs, elxs5  ! dummy. elastic cross section is 0
      real(8):: impactp, impactp5 ! impact parameter fm
      real(8):: xs5s, xs5c
      integer:: tcode, tsubcode, tcharge
      type(ptcl):: target
      real(8)::g, beta 
!
      if(pj%code /= kgnuc ) then
         write(0, *) ' cAAxsec2 is for heavy ion while '
         write(0, *) ' projectile =',pj%code, pj%subcode, pj%charge
         stop
      endif
      Etot = pj%fm%p(4)  ! GeV
      Ekt =  Etot - pj%mass
      EktMeV=Ekt*1000.  ! MeV
      Ekpn = Ekt/pj%subcode
      pjA = pj%subcode
      pjZ = pj%charge
      if( tgA == 1.0 ) then
!            use pA (nA) system ; will be done in cinelx
         call cinelx(pj, tgA, tgZ, xs)
!         tcode = knuc
!         tsubcode = -1
!         tcharge = tgZ
!         call cmkptc(tcode, tsubcode, tcharge, target)
!         g = Etot/pj.mass
!         target.fm.p(4) = g* masp
!         call cinelx(target, pjA, pjZ, xs)
      else
         if( Ekt <=  0. ) then
            xs = 0.
         elseif( AAXsec == 1 ) then   ! normalize to Shen's Xs
            if( Ekpn <  Enorm ) then
!                              ! this energy total kinetic E in MeV !!
               call shen(pjA, pjZ, EktMeV, tgA, tgZ, xs, elxs, impactp)
               xs = xs*1000.    ! output is in b
            else
               call cAAxsec0(pjA, Ekpn, tgA, xs)
               EktMeV= Enorm*pj%subcode*1000.
               call shen(pjA, pjZ, EktMeV,
     *           tgA, tgZ, xs5s, elxs5, impactp5)
               xs5s = xs5s*1000. ! output is in b
               call cAAxsec0(pjA, Enorm, tgA, xs5c)
               xs = xs * xs5s/xs5c
            endif
         elseif( AAXsec == 0 ) then  ! normalize to Cosmos xs
            if( Ekpn >  Enorm ) then
               call cAAxsec0(pjA, Ekpn, tgA, xs)
            else
!                              ! this energy total kinetic E in MeV !!
               call shen(pjA, pjZ, EktMeV, tgA, tgZ, xs, elxs, impactp)
               xs = xs*1000.    ! output is in b
               EktMeV= Enorm*pj%subcode*1000.
               call shen(pjA, pjZ, EktMeV, 
     *              tgA, tgZ, xs5s, elxs5, impactp5)
               xs5s = xs5s*1000. ! output is in b
               call cAAxsec0(pjA, Enorm, tgA, xs5c)
               xs = xs * xs5c/xs5s
            endif
         else
            write(0,*) 'AAXsec =',AAXsec, ' not usable'
            stop
         endif
      endif
      end
      subroutine cAAxsec0(pjA, Ekpn, tgA,  xs)
      implicit none
      real(8),intent(in):: pjA ! proj. A
      real(8),intent(in):: Ekpn ! proj.'s kinetic E /n GeV
      real(8),intent(in):: tgA ! target A
      real(8),intent(out):: xs ! in mb

      real(8):: p, t
      real(8),parameter:: Ekb(5) = (/5., 10., 100., 700., 3000./)  
!                            sigma( pp)  Ek power dependence above
!                            energy of Ekb(i) : Ek**d
      real(8),parameter:: delta(5)=(/-0.05,-0.01, 0.02, 0.05, 0.083/)
      real(8)::d
      real(8):: pw
      integer:: i, j

      p = pjA**0.3333
      t = tgA**0.3333 
      xs = 52.2 *( p + t - 
     *     1.189 * exp(- 0.055*min( pjA, tgA)))**2
      d = 0.
      do i = 5, 1, -1 
         if(Ekpn > Ekb(i)) then
            pw  = 1.
            do j = 1, i-1
               d = delta(j)
               pw = pw * (Ekb(j+1)/Ekb(j))**(2.5* d/(p+t+p*t/2.))
            enddo
            d = delta(i)
            pw = pw *(Ekpn/Ekb(i))**(2.5* d/(p+t+p*t/2.))
            xs = xs *pw
            exit
         endif
      enddo

      end
