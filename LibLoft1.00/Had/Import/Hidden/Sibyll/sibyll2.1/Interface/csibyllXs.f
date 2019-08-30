      module modsibyllXs
      logical,save:: sibyllXsUsed
      end      module modsibyllXs
      subroutine csibyllXs(pj, tg, xs)
!  sibyll can compute  xs for
!      p-p, pi-p, K-p
!      p-Air, pi-Air, K-Air, A-Air (A=2~56)
!      for other  targets, we employ cosmsos standard
!  However, it seems possilbe to generate hadronic 
!   collision events for targets A'=2~??(56?)
!    projectile of n, nba, pbar, Kch (Kl, Ks) are also allowed
!
      use modsibyllXs
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmass.h"
      type(ptcl):: pj  !input projectile particle
      type(ptcl):: tg  ! //   target particle.
                        ! if subcode  is 0, Air is assumed 
                        ! in sibyll
      real(8),intent(out)::xs  ! obtained x-section in mb
      
      integer::TA, PA, L
      real(4):: roots, sxs
      real(8)::At, Zt 
      real(4):: SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO
      logical ok
      logical waround
      integer,save:: waroundn=0
      integer:: IA

      L = 0   ! p,pi,K code; not given yet 
      PA = 0  ! projectile is A; not given yet
      ok = .false.  ! sibyll could  give xs; not yet 
      waround =.false.
      
      if( pj%code == knuc) then 
         L = 1
      elseif( pj%code == kpion  ) then
         L = 2
      elseif( pj%code == kkaon ) then
         L = 3
      elseif( pj%code == kgnuc .and. pj%subcode <=56 ) then
         PA = pj%subcode
      else
         if( pj%code == keta  .or.  pj%code == ksigma .or.
     *        pj%code == klambda ) then
         else
            if( waroundn < 10 ) then
               write(0,*)
     *       '********* projectile is N.G for XS calc.in Sibyll2.3c:'
               write(0,*) '********* pj%code, subcode, charge=',
     *              pj%code, pj%subcode, pj%charge
               write(0,*) '********* E=',pj%fm%p(4)
               write(0,*) '********* some workaround is taken'
               waroundn = waroundn + 1
            endif
         endif
         waround = .true.
      endif
      if(.not. waround) then
         if( tg%code ==9 .and. tg%subcode == 0 ) then ! air target               
            if( PA > 0 ) then
               call  csibyllXsAAir(PA, pj%fm%p(4), xs)
               ok = .true.
            else                !  p,pi,K- air xsec.
!           s = ( Epj + massT )**2 - Ppj**2
               roots =
     *              sqrt(2*masp*pj%fm%p(4) + masp**2 + pj%mass**2)
               call SIB_SIGMA_HAIR (L, roots, sxs)
               xs = sxs
               ok = .true.
            endif
         elseif( tg%code == knuc .and. tg%charge == 1 ) then
                   ! target p
            if( L > 0 ) then    !  p,pi,K-p xsec. 
               SQS =
     *              sqrt(2*masp*pj%fm%p(4) + masp**2 + pj%mass**2)
               call SIB_SIGMA_HP(
     *              L,SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)
               xs = SIGINEL
               ok = .true.
            else
!             not ok. see later
            endif
         elseif( tg%code == kgnuc ) then
            if( L > 0 ) then
!                no practical routine for hadron -A
!                cross-section.  see later
            else
!                 not ok. see later
            endif
         else
            write(0,*) ' target error in csibyllXs'
            write(0,*) ' tg code, sub, chge=',tg%code,
     *              tg%subcode, tg%charge            
         endif
      endif
      
      if( .not. ok ) then
         if( tg%subcode == 0 ) then
            At = 14.45
         elseif(tg%code /= kgnuc ) then
            At = 1.
         else
            At = tg%subcode 
         endif
         Zt = tg%charge ! not important
         call cinelx(pj, At, Zt, xs)
      endif

      sibyllXsUsed = ok
      end  subroutine csibyllXs

