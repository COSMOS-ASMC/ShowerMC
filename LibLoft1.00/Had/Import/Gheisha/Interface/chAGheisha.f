!             hadron A collision by Gheisha
        subroutine chAGheisha(pj, ia, iz,  a, np)
        implicit none

#include  "Zptcl.h"
#include  "Zcode.h"
!
        type (ptcl):: pj  ! input. projectile particle
        type (ptcl):: a(*) ! output. produced particles
        integer ia    !  input target nucleus mass number.
        integer iz    ! input. target nucleus charge
        integer np    !  output. number of particles produced.

!

!
        integer maxp
        parameter (maxp = 100)  ! maximum particle number allowed.

        type (ptcl)::aPtcl, bPtcl(maxp), cPtcl(maxp)

        integer ghecode,  retcode, stopped, j, nc
        integer  subcode, charge, k, i

        real*4 rho/1./, radl/1./, absl/1./  ! dummy setting
        real*4 aa4, zz4

        real*4 pin(3), tprod(maxp), pprod(3, maxp)
        integer iprod(maxp)

         
         call ccos2gheCode(pj, ghecode)
         aa4 = ia
         zz4 = iz
!           inform target matter
         call gsmate(1, aa4, zz4, rho, radl, absl)
         pin(3) = sqrt(abs(pj%fm%p(4)**2 - pj%mass**2))
         pin(1) =  0.
         pin(2) = 0.
!         pin(2) = pj.fm.p(2)
!         pin(3) = pj.fm.p(3)


!           generate an event
         call gheish(1, ghecode, pin, maxp, np, iprod, tprod,
     *   pprod, retcode, stopped)
!         retcode 13:  elastic scattering
!                 20:  inelastic scattering
!                 15:  nuclear fission + inelastic scat 
!                 18:  neutron capture
!             
        do j = 1, np
           call cghe2cos(iprod(j), k, subcode, charge)
           call cmkptc(k, subcode, charge, aPtcl)
           aPtcl%fm%p(1) = pprod(1,j)
           aPtcl%fm%p(2) = pprod(2,j)
           aPtcl%fm%p(3) = pprod(3,j)
           call cpm2e(aPtcl, aPtcl)
           bPtcl(j) = aPtcl
        enddo
!             make  simga etc decay. Since Gheisha is
!        used at low energies, there is no problem of
!        doing this here.

        do while(.true.)

           call cgheDecay(bPtcl, np, cPtcl, nc)
           if(nc .eq. np) goto 10
           call cgheDecay(cPtcl, nc, bPtcl, np)
           if(nc .eq. np) goto 10

        enddo
 10     continue

!           move bPtcl to a(*)
        np = 0
        do i =1, nc
!              omit very rare particles such as tau.
           if(bPtcl(i)%code .ne. krare) then
              np = np + 1
              a(np) = bPtcl(i)
           endif
        enddo
!        rotate
        call crot3mom(pj, a, np)
      end
