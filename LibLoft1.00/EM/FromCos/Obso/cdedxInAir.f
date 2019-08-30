!#include  "BlockData/cblkdedx.h"
!#include  "BlockData/cblkElemag.h"
!           test cdedxInAir
!      dedt in  (GeV)/(kg/m2)
!       so it is 1000 MeV/ ((1000g)/ 10000 cm2) = 1.e4 MeV/(g/cm2)
!      implicit none
!
!#include "Zptcl.h"
!#include "Zmass.h"
!#include  "Zcode.h"
!#include "Zelemagp.h"

!      type(ptcl):: aPtcl

!      integer j
!      real*8  rho, dedt,  ek
!
!      
!      Knockon =.false.
!      call cmkptc(kelec, regptcl, -1, aPtcl)
!      do  j=1, 10
!             rho=2.e-3*10.**(-(j-1)/2.)
!             ek = aPtcl.mass*1.00001 - aPtcl.mass
!             do  while (ek .lt. 100.)
!                 aPtcl.fm.p(4) = ek + aPtcl.mass 
!                 call cdedxInAir(aPtcl, rho, dedt)
!                 write(*, *) sngl(ek), sngl(dedt*1.e4)
!                 ek=ek*10.**0.1
!             enddo
!             write(*, *)
!      enddo
!      end
!     ****************************************************************
!     *                                                              *
!     * cdedxInAir: gives -de/dx (GeV/(kg/m**2)) of e+/e- mu, pi, k,
!     *              p etc in air.
!     ****************************************************************
!
! /usage/  call cdedxInAir(ptcl, rhoin, dedt)
! -- input--
!     ptcl:  type ptcl
!    rhoin: density of air in kg/m**3
! -- output --
!     dedt; energy loss   GeV/(kg/m**2) 
!
!
!
         subroutine cdedxInAir(aPtcl, rhoin, dedt, dedtF)
         implicit none

#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zelemagp.h"

         type(ptcl):: aPtcl ! input paticle
         real*8  rhoin   !  input.  density 
         real*8  dedt    !  output. restriced energy  loss rate
         real*8  dedtF   !  output. full energy loss rate
!
         real*8 erg
         integer charge

!            next is moved to cbeginRun
!         if(first) then
!            call cdedxEleci(RecoilKineMinE, Knockon)
!            first = .false.
!         endif

         erg = aPtcl%fm%p(4)
!
         if(aPtcl%code .eq. kelec) then
            charge = aPtcl%charge
            call cdedxe(aPtcl, rhoin, dedt, dedtF)  ! in GeV/(g/cm2)
            dedt = dedt / 10.                ! to GeV/(kg/m2)
            dedtF = dedtF/10. 
         else
            call cdedxNone(aPtcl, rhoin,  dedt, dedtF) ! in GeV/(g/cm2)
            dedt = dedt / 10.                ! in GeV/(kg/m2)
            dedtF = dedtF/10.
         endif

         end

            

