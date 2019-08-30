      subroutine cnelas(pj, ia, iz, a, ntp)
!      implicit none
#include "param00.inc"
#include "Zptcl.h"
#include "Zmass.h"
      type(ptcl):: pj
      integer,intent(in):: ia, iz ! target A,Z  
      type(ptcl):: a(*)  ! normally size should  be  2
      integer,intent(out):: ntp

      integer ityp, ktyp
      real(8):: ata, atz, epin

      integer ncount
      
      common /clustf/ nclst, iclust(nnn)


      integer:: code, subcode, charge
      code = pj%code
      subcode = pj%subcode
      charge = pj%charge

      call ccos2phits(code, subcode, charge, ityp, ktyp)
      ata=ia
      atz =iz
      epin =( pj%fm%p(4) - pj%mass)*1000.   !  K%E in MeV
      ncount = 0
      do while (ncount < 100)
         call nelst(ityp,ktyp,ata,atz,epin)
         if( nclst > 0 ) exit
         ncount = ncount + 1
      enddo
!///////////
!      write(0,*) 'nelstLoop=', ncount
!////////////
!      call nevap(0)    ! Is it ok to call evaporation routine
!                       !  after elastic scattering
!                       ! otherwise we must move array
      call cphitsOut(1, pj, ia, iz, ntp, a)
      call crot3mom(pj, a, ntp)  ! to the  mother system.
      end subroutine


