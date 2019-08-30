!   This version is faster than cl2tT and accurracy is comparable.
!
!     *********************************************
      real*8 function clen2thickTA(z,  leng)
!
!     z: real*8. input.  vertical height in m.
!   leng: real*8. input.  length along cosz direction.
!     function value.      thickness of air in kg/m2. for leng
      use modAtmosDef
      implicit none
!  #include "Zatmos.h"
      real*8 z,  leng
      real*8  s1, t1, s2, t2

!          get slant length along cosz for z:   s1
!          get slant length ; s2 = s1 + leng
!          get thickness <-- by T vs slant length
      call cl2tIntp(HeightTbl, LenTbl, z, s1)
      call cl2tIntp(LenTbl, ThickTbl, s1, t1)
      s2 = s1 + leng
      call cl2tIntp(LenTbl, ThickTbl, s2, t2)
      clen2thickTA = t2 - t1   		
      if(clen2thickTA  .le. 0.) then
!          actually case of < 0 dose not happen.
!          when particle comes to the last segment of 
!          the table,  this happens.  and the
!          particle does not move at all. To avoid this
!          give some thikness. (actually this is not
!          used but this is for comfort only)
         clen2thickTA  = 10.
      endif
      end
!          
      subroutine cl2tIntp(x, y, xx, ans)
      use modAtmosDef
      implicit none
! #include "Zatmos.h"

	real*8 x(*)  !  input. 
        real*8 y(*)  !  input.
        real*8  xx   !  input. some value inside x(*).
        real*8  ans  !  output.  interpolated value of y at xx.

	real*8 error
	integer  m
	integer loca, k
	parameter (m = 3)  ! old value was 5
!            where is xx in x
	call kdwhereis(xx, NumStep, x, 1, loca)
        if(loca .ge. NumStep) then
!            give arbitray large value.
  	     ans = y(NumStep)
        else
!              use max of m points around loca for interpolation
             k = min(max(loca - (m-1)/2,1), NumStep+1-m) 
                            ! max of m points from k
             call kpolintp(x(k), 1, y(k), 1, m, xx, ans, error)
	endif
        end
!         get slant length  corresponding to slant thickness
!
      real*8 function ct2lTA(z,  t)
      use  modAtmosDef
      implicit none
! #include "Zatmos.h"

      real*8 z  ! input vertical height in m
      real*8 t  !  input.  thickness of air in kg/m^2 along cosz
!
      real*8 t1, t2, s1, s2
!          z ---> thick 
      call cl2tIntp(HeightTbl, ThickTbl, z, t1)
      call cl2tIntp(HeightTbl, LenTbl,  z, s1)
      t2 = t1 + t
      call cl2tIntp(ThickTbl,  LenTbl, t2, s2)
      ct2lTA = s2- s1	
      if(ct2lTA .le. 0.) then
!           actually case of < 0 dose not happen.
!           when a particle comes to an end segment of
!           a table, this happnes and the particle 
!           cannot move. To avoid this, give some
!           finite value to move it.
!           clen2thickTA  a thickness corresponding to this length
!           but is not used since the particle reaches
!           the bottom observation depth.
         ct2lTA = 10.
      endif
      end
