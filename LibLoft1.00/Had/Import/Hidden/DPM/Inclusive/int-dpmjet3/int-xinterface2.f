      subroutine dist2x(iin,iout,ie,xx, x)
      IMPLICIT REAL*8 (A-H,O-Z)
      if(iin.eq.1) then
         call dist2x1(iout,ie,xx, x)
      else if(iin.eq.2) then
         call dist2x2(iout,ie,xx, x)
      else if(iin.eq.3) then
         call dist2x3(iout,ie,xx, x)
      else if(iin.eq.4) then
         call dist2x4(iout,ie,xx, x)
      else if(iin.eq.5) then
         call dist2x5(iout,ie,xx, x)
      else if(iin.eq.6) then
         call dist2x6(iout,ie,xx, x)
      else if(iin.eq.7) then
         call dist2x7(iout,ie,xx, x)
      else if(iin.eq.8) then
         call dist2x8(iout,ie,xx, x)
      else if(iin.eq.9) then
         call dist2x9(iout,ie,xx, x)
      else
         stop 'input ptl error in hadron int'
      end if
      END

      subroutine dist2dx(iin,iout,ie, xx, dx)
      IMPLICIT REAL*8 (A-H,O-Z)
      if(iin.eq.1) then
         call dist2dx1(iout,ie,xx, dx)
      else if(iin.eq.2) then
         call dist2dx2(iout,ie,xx, dx)
      else if(iin.eq.3) then
         call dist2dx3(iout,ie,xx, dx)
      else if(iin.eq.4) then
         call dist2dx4(iout,ie,xx, dx)
      else if(iin.eq.5) then
         call dist2dx5(iout,ie,xx, dx)
      else if(iin.eq.6) then
         call dist2dx6(iout,ie,xx, dx)
      else if(iin.eq.7) then
         call dist2dx7(iout,ie,xx, dx)
      else if(iin.eq.8) then
         call dist2dx8(iout,ie,xx, dx)
      else if(iin.eq.9) then
         call dist2dx9(iout,ie,xx, dx)
      else
         stop ' input ptl error in hadron int'
      end if
      END

