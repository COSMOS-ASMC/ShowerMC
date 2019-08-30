      subroutine cgheDecay(a,  na, b,  nb)
!         to make particles in 'a' decay, if some of them 
!        are one of  sigma, gzai, lambda or  eta, 
!        No multistep decay is treated.  The 
!        decay product is placed in b.
!        In normal usage, if you want to get final particles
!
!        do while (.true.)
!           call cgheDecay(A, Na, B, Nb)
!           if(Na .eq. Nb) goto 10
!           call cgheDecay(B, Nb, A, Na)
!           if(Na .eq. Nb) goto 10
!        enddo
! 10     continue
!        you may use A and Na as the final decay product.
!        or B and Nb: They are the same.
!      
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
      integer na   ! input. number of particles in a. 
      type (ptcl):: a(500)   ! input particle list in Cosmos format
      integer nb   ! output. nubmer of particle in b.
      type (ptcl):: b(500)  ! output.  Decay product in Comsos format.

      type (ptcl):: aPtcl
      integer i, ny

      nb = 0
      do i = 1, na
         aPtcl = a(i)
         if(aPtcl%code .eq. ksigma) then
            call csigmaDecay(aPtcl, b(nb+1), ny)
         elseif(aPtcl%code .eq. klambda) then
            call clambdaDcy(aPtcl, b(nb+1), ny)
         elseif(aPtcl%code .eq. kgzai) then
            call cgzaiDecay(aPtcl, b(nb+1), ny)
         elseif(aPtcl%code .eq. keta) then
            call cetaDecay(aPtcl, b(nb+1), ny)
         else
            b(nb+1) = aPtcl
            ny = 1
         endif
         nb = nb + ny
      enddo
      end
