c	real*8 function rn(idumm)
c	integer idumm
c	real*8 u
c	call rndc(u)
c	rn =  u
c       end  function
       subroutine parastop(idm)
	write(0,*) ' parastop called with id=',idm
	stop 'parastop'
       end subroutine
	real*8 function  unirn(idumm)
	integer idumm
	real*8 u
	call rndc(u)
	unirn =  u
       end function
	real*8 function rang()
	real*8 u
	call rndc(u)
	rang = u
	end function

