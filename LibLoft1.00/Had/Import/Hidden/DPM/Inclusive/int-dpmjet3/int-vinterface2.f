      subroutine dist2v(iin,iout,ie, xx, v)
      IMPLICIT REAL*8 (A-H,O-Z)
      if(iin.eq.1) then
         call dist2v1(iout,ie,xx, v)
      else if(iin.eq.2) then
         call dist2v2(iout,ie,xx, v)
      else if(iin.eq.3) then
         call dist2v3(iout,ie,xx, v)
      else if(iin.eq.4) then
         call dist2v4(iout,ie,xx, v)
      else if(iin.eq.5) then
         call dist2v5(iout,ie,xx, v)
      else if(iin.eq.6) then
         call dist2v6(iout,ie,xx, v)
      else if(iin.eq.7) then
         call dist2v7(iout,ie,xx, v)
      else if(iin.eq.8) then
         call dist2v8(iout,ie,xx, v)
      else if(iin.eq.9) then
         call dist2v9(iout,ie,xx, v)
      else
         stop 'input ptl error in hadron int: dist2v'
      end if
      END

      subroutine dist2dv(iin,iout,ie, xx, dv)
      IMPLICIT REAL*8 (A-H,O-Z)
      if(iin.eq.1) then
        call dist2dv1(iout,ie,xx, dv)
      else if(iin.eq.2) then
        call dist2dv2(iout,ie,xx, dv)
      else if(iin.eq.3) then
        call dist2dv3(iout,ie,xx, dv)
      else if(iin.eq.4) then
        call dist2dv4(iout,ie,xx, dv)
      else if(iin.eq.5) then
        call dist2dv5(iout,ie,xx, dv)
      else if(iin.eq.6) then
        call dist2dv6(iout,ie,xx, dv)
      else if(iin.eq.7) then
        call dist2dv7(iout,ie,xx, dv)
      else if(iin.eq.8) then
        call dist2dv8(iout,ie,xx, dv)
      else if(iin.eq.9) then
        call dist2dv9(iout,ie,xx, dv)
      else
         stop ' input ptl error in hadron int: dist2dv'
      end if
      END
