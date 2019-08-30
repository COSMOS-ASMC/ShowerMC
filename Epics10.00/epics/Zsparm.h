       type(epDirec)::  DCInp
       type(epPos)::   PosInp

        complex(8):: Xrange, Yrange, Zrange
	real*8  E0Max
	real*8 Hwhm, ProfR
        complex(8)::  CosNormal
        logical LogIr
        character*12 InputA
        character*8 InputP

!             sepics common
         common /Zspc/ 
     3   InputA,  InputP

         common /Zsp/ Xrange, Yrange, Zrange,
     1   DCInp, PosInp, CosNormal, 
     2   E0Max, Hwhm,  ProfR,  LogIr

!
! IoCont: fortran dev # for cont dataset.
!


