  subroutine kxplPipe(p, dir, r1, r2, h, el, icon1, icon2)
    implicit none
    real(8),intent(in):: p(3)
    real(8),intent(in):: dir(3)
    real(8),intent(in):: r1, r2, h
    real(8),intent(out):: el
    integer,intent(out):: icon1
    integer,intent(out):: icon2

    integer:: icon, nx, icon2a(2)
    real(8):: ela(2), eltemp
    
    call kspipe(p, r1, r2, h, icon1) 
    
    nx = 0
    call kxplCylWall(p, dir, r2, h, eltemp, icon)
    if( icon == 0 ) then
       nx = nx + 1
       ela(nx) = eltemp
       icon2a(nx) = 6
    endif
    call kxplCylWall(p, dir, r1, h, eltemp, icon)
    if( icon == 0 ) then
       nx = nx + 1
       ela(nx) = eltemp
       icon2a(nx) = 4
    endif
    if( nx == 2 ) goto 20
    call kxplPipeCeiling(p, dir, r1, r2, h, eltemp, icon)
    if(icon == 0 ) then
       nx = nx + 1
       ela(nx) = eltemp
       icon2a(nx) = 1
    endif
    if( nx == 2 ) goto 20

    call kxplPipeFloor(p, dir, r1, r2,  eltemp, icon)
    if(icon == 0 ) then
       nx = nx + 1
       ela(nx) = eltemp
       icon2a(nx) = 2
    endif
20  continue
    if(nx == 0 ) then
       icon2 = -1
    elseif(nx == 1) then
       el = ela(1)
       icon2 = icon2a(1)
    elseif( ela(1) < ela(2) ) then
       icon2 = icon2a(1)
       el = ela(1)
    else
       icon2 = icon2a(2)
       el = ela(2)
    endif

  end subroutine kxplPipe

  subroutine kspipe(pos, r1, r2, h, icon)
    implicit none
    real(8),intent(in):: pos(3)
    real(8),intent(in):: r1, r2, h
    integer,intent(out):: icon

    real(8)::x2y2

    x2y2 = pos(1)**2 + pos(2)**2
    if( x2y2 < r1**2 ) then
       icon = 1
    elseif( x2y2 > r2**2 ) then
       icon = 1
    elseif( pos(3) < 0 ) then
       icon = 1
    elseif( pos(3) > h ) then
       icon = 1
    else
       icon = 0
    endif
  end subroutine kspipe

  program testPipe
    implicit none
    real(8):: u, p(3), dir(3), cs, sn, s
    real(8):: r1=3, r2=5, h=10

    integer::i, icon1, icon2, ir(2)
    
    real(8):: eldummy=10, el
!      ir(1)=1688068097
!      ir(2)=1694128957 
      ir(1)= 2021831529 
      ir(2) =    561239215
      call rnd1r(ir)
    do i = 1, 100000
!       call rnd1s(ir)
!       write(0,*) 'ir=',ir(:)
       call randompdir(r2+r1, r2+r1, h+3, p,  dir)
       call kxplPipe(p, dir, r1, r2, h, el, icon1,icon2)
       if(icon2 /= -1) then
          write(*,'(0p,3g15.5,2i3)') p(:)+ el*dir(:), icon2, icon1
!       if(icon1 /= -1) then
!          write(*,'(1p,3g14.4, 0p,i3)')  p(:), icon2
!          if(icon2 /= -1) then
!             write(*,'(1p,3g14.4,0p,i3)') p(:) + el*dir(:), icon2
!          else   
!             write(*,'(1p,3g14.4,0p,i3)') p(:) + eldummy*dir(:), icon2
!          endif
!          write(*,*) 
!          write(*,*) 
       endif
    enddo
  end program testPipe
       
    

  subroutine randompdir(a, b, c, p, dir)
    implicit none
    real(8),intent(in):: a, b, c
    real(8),intent(out):: p(3), dir(3)
    real(8):: u, cs, sn, s

    call rndc(u)
    p(1) = (2*u-1)*a
    call rndc(u)
    p(2) = (2*u-1)*b
    call rndc(u)
    p(3) = u*(c+c/4) -c/4
    call kcossn(cs, sn)
    call rndc(u)
    s=sqrt(1.-u**2)
    dir(1) = s*cs
    dir(2) = s*sn
    dir(3) = 2*u-1
  end subroutine randompdir
       
