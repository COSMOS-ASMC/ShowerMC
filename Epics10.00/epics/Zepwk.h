!                for 
!                   stacking area
        common /Zepwk/ xa(nstkmx), ya(nstkmx), za(nstkmx), ta(nstkmx),
     1     ea(nstkmx), wxa(nstkmx), wya(nstkmx), wza(nstkmx), dl,
     2     dt, xp, yp, zp, xpp, ypp, zpp, w1, w2, w3, ee1,
     3     bx, by, bz, ex, ey, ez, aary(nstkmx), firsti,
     4     firstc, Abort, papa(nstkmx)
!
        real*8   xa, ya, za, ea, wxa, wya, wza, dl, dt,
     1           xp, yp, zp, xpp, ypp, zpp, w1, w2, w3,
     2           ee1, ta, bx, by, bz, ex, ey, ez, aary

	integer papa, Abort
        logical firsti, firstc
!
        common /Zepwkc/ mat, matp, proc, matpp, matns
!
        character*8 mat, matp, matpp, matns
        character*4  proc
!
        common /Zepwk1/ fpth, xbmv, ybmv, zbmv,
     1   xinc,yinc,zinc,einc,wxinc,wyinc,wzinc,
     2   tinc, ekmin, ebmv, maxPathL(ncmax), maxPath,
     3   ica(nstkmx), ka(nstkmx), cna(nstkmx), nstk,
     4   jbnd, jde,  cnp, trunc, userb, usere, cnpp,
     5   jbnda(ncmax), jdea(ncmax), 
     6   nevrun, icinc, kinc,
     7   icninc

	real*8  fpth, xbmv, ybmv, zbmv, xinc, yinc, zinc, einc
	real*8 wxinc, wyinc, wzinc, tinc, ekmin, maxPathL
	real*8 maxPath, ebmv
        integer ica, ka, cna, cnp, cnpp, nevrun, nstk, jbnd, jde
	integer jbnda, jdea, icinc, kinc, icninc
        logical trunc, userb, usere

! nevrun: # of events created in the current run
!  nstk:  stack area counter
! matns:  to keep the previous non space non void matter.
!  trunc: becomes t if path is truncated eles f.
!  fpth:  mean free path of the particle for the process given in proc.
!         (cm)
!  userb: becomes t if user wants to count particles at the component
!         boundary (call epuif1(xxx, yyy) must be issued to make
!         userb=t, where xxx is the subroutine name to count forward
!         going ptcls, and yyy the one for backward going ptcls.
!  usere: becomes t if user wants to measure energy deposit in some
!         component.  call epuif2(xxxx) must be issued to make
!         usere=t, where xxxx is the subroutine name to count energy
!         deposit.
!
