  subroutine kxplCylWall(p, dir, r, h, el, icon)
    implicit none
    real(8),intent(in):: p(3)
    real(8),intent(in):: dir(3)
    real(8),intent(in):: r, h
    real(8),intent(out):: el
    integer,intent(out):: icon
    
    real(8):: a, b, c, D, sqrtD 
    real(8):: el1, el2, z1, z2

    a =dir(1)**2 + dir(2)**2
    if(a == 0.) then
       icon = -1
    else
       b = p(1)*dir(1) + p(2)*dir(2)
       c = p(1)**2 + p(2)**2 - r**2
!       solve a*el**2 + 2b*el + c = 0
       D = b**2 - a*c
       if(D >= 0.) then
          sqrtD = sqrt(D)
          el1 = (-b + sqrtD) / a
          z1 = el1*dir(3) + p(3)
          el2 = (-b - sqrtD) / a
          z2 = el2*dir(3) + p(3)

          if( el1 > 0. .and. z1 > 0. .and. z1 < h) then 
             el=el1
             icon = 0
             if( el2 > 0. .and. z2 > 0. .and. z2 < h) then
                el = min(el, el2)
             endif
          elseif( el2 > 0. .and. z2 > 0. .and. z2 < h) then
             el = el2
             icon = 0
          else
             icon = -1
          endif
       else
          icon= -1
       endif
    end if

  end subroutine kxplCylWall

  subroutine kxplPipeCeiling(p, dir, r1, r2, h, el, icon)
    implicit none
    real(8),intent(in):: p(3)
    real(8),intent(in):: dir(3)
    real(8),intent(in):: r1, r2, h
    real(8),intent(out):: el
    integer,intent(out):: icon
    
    real(8):: x, y, z, temp
    if(dir(3) == 0.) then
       icon = -1
    else
       !   el*dir(3) + p(3) = h
       el =(h-p(3))/dir(3)
       x = el*dir(1) + p(1)
       y = el*dir(2) + p(2)
       temp = x**2 + y**2
       if( el <= 0. ) then
          icon = -1
       elseif( temp < r1**2 .or. temp > r2**2) then
          icon = -1
       else
          icon = 0
       endif
    endif

  end subroutine kxplPipeCeiling

  subroutine kxplPipeFloor(p, dir, r1, r2,  el, icon)
    implicit none
    real(8),intent(in):: p(3)
    real(8),intent(in):: dir(3)
    real(8),intent(in):: r1, r2
    real(8),intent(out):: el
    integer,intent(out):: icon
    
    real(8):: x, y, temp
    if(dir(3) == 0.) then
       icon = -1
    else
       !   el*dir(3) + p(3) = 0
       el =-p(3)/dir(3)
       x = el*dir(1) + p(1)
       y = el*dir(2) + p(2)
       temp = x**2 + y**2
       if(el <= 0.) then
          icon = -1
       elseif( temp <  r1**2 .or. temp > r2**2) then
          icon = -1
       else
          icon =  0
       endif
    endif


  end subroutine kxplPipeFloor

   
!  program  test
!    implicit none
!    real(8):: p(3), dir(3)
!    real(8):: u,  cs, sn, s, el
!    real(8):: r1=3, r2=5, h=10
!    integer:: icon, i
!
!    do i = 1, 100000
!!       write(0,*) 'Wall test r=', r1, ' h=',h
!!       write(0,*) ' Enter point p(:)'
!!       read(*,*) p(:)
!!       write(0,*) ' p=', p(:)
!!       write(0,*) ' Enter dir. dir(3) is only -1 or 1 sign'
!!       read(*,*) dir(:)
!!       dir(3) = sign(sqrt(1.-dir(1)**2 - dir(2)**2), dir(3))
!!       write(0,*) ' dir=', dir(:)
!       call randompdir(r1+2, r2+2, h+5, p, dir)
!!       call kxplCylWall(p, dir, r1, h, el, icon)
!       call kxplPipeCeiling(p, dir, r1, r2, h, el, icon)
!       if(icon == 0 ) then
!          write(*,'(1p,3g15.5)') p(:)+el*dir(:)
!       endif
!       call kxplPipeFloor(p, dir, r1, r2,  el, icon)
!       if(icon == 0 ) then
!          write(*,'(1p,3g15.5)') p(:)+el*dir(:)
!       endif
!!       write(0,*) 'Wall icon=',icon, ' el=',el
!    enddo
!    
!
!  end program test
!
!  subroutine randompdir(a, b, c, p, dir)
!    implicit none
!    real(8),intent(in):: a, b, c
!    real(8),intent(out):: p(3), dir(3)
!    real(8):: u, cs, sn, s
!
!    call rndc(u)
!    p(1) = (2*u-1)*a
!    call rndc(u)
!    p(2) = (2*u-1)*b
!    call rndc(u)
!    p(3) = u*(c+c/4) -c/4
!    call kcossn(cs, sn)
!    call rndc(u)
!    s=sqrt(1.-u**2)
!    dir(1) = s*cs
!    dir(2) = s*sn
!    dir(3) = 2*u-1
!  end subroutine randompdir
!       

    
    



    
    
    


