	use jqmd
	implicit none
#include "Zptcl.h"

	external qmddata

	record /ptcl/ a(1000)
	integer n, i, j
	integer ityp,ktyp,mmas,mchg
	real(8):: eein, bmax0
	call cprePhits
	call ccos2phits(9, 56, 26, ityp, ktyp)
c/////////
	write(0,*) ' itype, ktyp=', ityp, ktyp
c////////
	eein = 500*56
c	mmas = 14
c	mchg =7
	mmas = 1
	mchg = 1
	bmax0=7.4
	do i = 1, 1000
	   call jqmdin(ityp,ktyp,eein,mmas,mchg,bmax0)
	   call cphitsOut(n, a)
	   do j=1, n
	      write(*,'(4i4, 1p, 5g12.4)')
     *        j, a(j).code, a(j).subcode, a(j).charge, a(j).fm.p(:),
     *        a(j).mass
	   enddo
	   write(*,*) ' '
	enddo
       end program
#include "qmddflt.f"
	
