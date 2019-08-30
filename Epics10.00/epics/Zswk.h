!             sepics common
!
         real*8   rpipo, rpipi, rcyl


         logical cont
	 integer  Ir1st(2),  nwld, nevc
         character*4 formx
         character*8 wstruct 
       type(epPos)::  orgw
       type(ep3Vec)::  abcw
         common /Zswk/ orgw, abcw,
     *    rcyl, 
     *    rpipo, rpipi,  
     *    cont, Ir1st, nwld, nevc
         common /Zswkc/ wstruct, formx


!
!
! cont:  to be f for the first run.  t for the cont run.
! Ir1st: to keep the initial random # for each event.
!
