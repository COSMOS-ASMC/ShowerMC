!   Manipulate input data (Lightmn.dat or Sencormn.dat)
!   after it is read.  (but refraction  index, if not, given,
!   is manipulated inside epLightCoutDE since it needs Media info
!   (interplay between old and new fortran.)
!
subroutine epLightManipInp(cmpno, pty, each )
  use modepLight
  use modepLightPty
  implicit none
  integer,intent(in):: cmpno ! component number
  type(property),intent(in)::pty  ! property in comInfo
  type(eachInfo),target,intent(out)::each  ! eachInfo in Lcomp



  integer:: na    ! volume attribute (box->3, cyl->2, ecyl->3
  real(8),target:: vol(maxVolAttr)  ! at least with na dim.
  character(len=20)::struct  !  actual max is 12 now
  real(8), pointer::boxa, boxb, boxc, r1, r2, height

  integer,pointer :: nx, ny, nz
  real(8),pointer ::dx, dy, dz
  real(8) ::safefac=(1.0d0 + 1.d-12)

  real(8)::wl0   ! dummy usse


  if(each%mn >= LightNoMin .and. each%mn < SensorNoMin) then
     call epqvolatr(cmpno, na, vol)  ! vol(i), i = 1, na  
     call epqstruc(cmpno, struct)    ! vol structure ; box etc

     nx => each%nx
     ny => each%ny
     nz => each%nz

     dx => each%dx
     dy => each%dy
     dz => each%dz

     if( struct(1:3) == "box" )  then
        boxa => vol(1)
        boxb => vol(2)
        boxc => vol(3)

        nx = 100./pty%CellSize(1)
        dx = boxa/nx 
        ny = 100./pty%CellSize(2)
        dy = boxb/ny

        nz = 100./pty%CellSize(3)
        dz = boxc/nz
     elseif( struct(1:7) == "octagon" )  then
        boxa => vol(1)
        boxb => vol(2)
        boxc => vol(3)

        nx = 100./pty%CellSize(1)
        dx = boxa/nx 
        ny = 100./pty%CellSize(2)
        dy = boxb/ny

        nz = 100./pty%CellSize(3)
        dz = boxc/nz
     elseif(struct(1:3) == "cyl" .or. struct(1:4) == "ecyl"  ) then
        r1 => vol(1)
        if(struct(1:3) == "cyl" ) then
           r2 => r1
           height => vol(2)
        else
           r2 = vol(2)
           height => vol(3)
        endif
        nx = 100./pty%CellSize(1)
        dx = r1/nx
        
        ny = 100./pty%CellSize(2)
        dy = r2/ny

        nz = 100./pty%CellSize(3)
        dz = height/nz
     elseif(struct(1:4) == "pipe"  ) then
        r1 => vol(1)    ! inner radius
        r2 => vol(2)     ! outer /
        height => vol(3)
        nx = 100./pty%CellSize(1)
     !    dx = (r2-r1)/nx
        dx = r2/nx
        ny = 100./pty%CellSize(2)
        dy = (r2-r1)/ny
        nz = 100./pty%CellSize(3)
        dz = height/nz
     else
        write(0,*) " compoent #=",cmpno, " structure=",struct
        write(0,*) " is not yet supported for ray tracing"
        stop
     endif
     dx = dx * safefac
     dy = dy * safefac
     dz = dz * safefac

     each%dmin=min(dx,dy,dz)

         !         ElightH, ElightL for light w.l
     call epLightwl2E(pty%minmaxWL(1),  1.d0, wl0, pty%ElightH)
     call epLightwl2E(pty%minmaxWL(2),  1.d0, wl0, pty%ElightL)

     if( pty%WLS(1) == 1. ) then
        if(pty%WLS(2) > pty%WLS(3)) then
           write(0,*) ' 2 and 3rd of  WLS values=',pty%WLS(1:4),  ' wrong '
           stop
        endif
        if( pty%WLS(4) < pty%WLS(3) ) then
           write(0,*) ' 4th of WLS values=',pty%WLS(1:4),  ' wrong '
           stop
        endif
     endif
!     /////////
!     write(0,*) 'pty%minmaxWL(1),  1.d0, wl0, pty%ElightH'
!     write(0,*) pty%minmaxWL(1),  1.d0, wl0, pty%ElightH
!     write(0,*) 'pty%minmaxWL(2),  1.d0, wl0, pty%ElightL'
!     write(0,*) pty%minmaxWL(2),  1.d0, wl0, pty%ElightL
!      /////////
  else
!     sensor. nothing to do at present
  endif


end subroutine epLightManipInp
