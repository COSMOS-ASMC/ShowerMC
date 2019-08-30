#include "Zptcl.h"
	record /ptcl/ a(100)
	integer n
	integer ityp,ktypmmas,mchg
	real(8):: eein, bmax0
	ityp =19
	call ccos2kf(9, 56, 26, ktyp)
	eein =500
	mmas = 14
	mchg =7
	bmax0=8.
	call jqmdin(ityp,ktyp,eein,mmas,mchg,bmax0)
        call cjqmdout(n, a)
	end
c	integer function  ichgf(ityp,ktyp)
c	ichgf = 0
c	end
c	integer function ibryf(ityp,ktyp)
c	ibryf = 1
c	end
c	real*8 function pcmsr(a,b,c)
c	real*8 a, b, c
c	pcmsr = 0.0
c	end
c	real*8 function bindeg(nz,nn)
c	bindeg = 0.
c	end
	real*8 function rn(idumm)
	integer idumm
	real*8 u
	call rndc(u)
	rn =  u
	end
	subroutine parastop(idm)
	stop 'parastop'
	end
	real*8 function  unirn(idumm)
	integer idumm
	real*8 u
	call rndc(u)
	unirn =  u
	end





