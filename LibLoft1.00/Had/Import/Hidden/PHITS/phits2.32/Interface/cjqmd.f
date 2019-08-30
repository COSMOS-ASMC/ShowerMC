	subroutine cjqmd(pj, ia, iz, sig, a, ntp)

	use modnevap
	implicit none
#include "Zcode.h"
#include "Zptcl.h"
	type(ptcl)::pj   ! input. projectile.
	integer,intent(in)::ia ! target A (nucleon #)
	integer,intent(in)::iz ! target Z
	real(8),intent(in)::sig ! cross-setion (mb) on this target
                         ! current  qmd00  dose not use but may  be
                         ! better to give it
	type(ptcl):: a(*)  ! output. produced particles
	      ! supoose Fe+Bi collsion.  If all Fe 
              ! nucleon of energy 4 GeV collide and produce
              !( 2pions + n  )*56 + 200 n; 3*56+200 = 380; max size of a.

	integer,intent(out):: ntp   ! total number of produced ptcls

!	external qmddata    ! shoul be put in main.

	integer n, i, j

	integer code, subcode, charge
	integer ityp,ktyp,mmas,mchg
	real(8):: eein, bmax0
	integer nout

	code = pj%code
	subcode  = pj%subcode
	charge = pj%charge
	call ccos2phits(code, subcode, charge, ityp, ktyp)
!	eein = 500*56    total kinetic energy in MeV
	eein = pj%fm%p(4)-pj%mass
	eein = eein *1000.d0    ! MeV

!	mmas = 14
!	mchg =7
!	mmas = 1
!	mchg = 1
!           sig = pi r**2*10 mb r in fm.
	bmax0= sqrt(sig/3.141592/10.)



	call cjqmdin(ityp,ktyp,eein, ia, iz, bmax0)
        call cphitsADJnp(pj, ia, iz) ! may try to adjust # of n,p
	call nevap(0)
        call cphitsOut(2, pj, ia, iz,  n, a)
	call crot3mom( pj, a, n ) ! rotate to the current  cooord
!	   do j=1, n
!	      write(*,'(4i4, 1p, 5g12.4)')
!     *        j, a(j).code, a(j).subcode, a(j).charge, a(j).fm.p(:),
!     *        a(j).mass
!	   enddo
!	   write(*,*) ' '
	ntp = n
	end subroutine 

	subroutine cjqmdin(ityp,ktyp,eein, ia, iz, bmax0)
	use jqmd
!ccc	implicit none
#include "param00.inc"
	! interface with jqmdin
	integer,intent(in):: ityp, ktyp  
	real(8),intent(in)::eein
	integer,intent(in):: ia, iz
	real(8),intent(in)::bmax0


	integer ncount
        common /clustf/ nclst, iclust(nnn)
	ncount = 0
	do while (ncount < 100)
	   call jqmdin(ityp,ktyp,eein, ia, iz, bmax0)
	   if( nclst > 0 ) exit 
	   ncount = ncount +1
	enddo
!////////////
!	write(0,*) 'jqmdinLoop=',ncount
!///////////
	end subroutine

	

	
	
