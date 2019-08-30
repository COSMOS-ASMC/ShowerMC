  subroutine kxpLand3vpB(p1, p2, p3,  &
       pos, dir, el,  towhich, icon)
!  Given  3 points p1,p2,p3 make a plane. 
!  Examin a given line passing pos with dir direction cosine
!  cross the triangle.   If so, the distance to the x-poit 
!  el is obtained.
    real(8),intent(in):: p1(3), p2(3), p3(3)
    real(8),intent(in):: pos(3), dir(3)

    real(8),intent(out):: el
    real(8),intent(out):: towhich
    integer,intent(out):: icon

    real(8)::xyz(3,3), nv(3), k, sumteta, xp(3)
    integer:: jcon
    character(4):: cond
    
    xyz(:,1) = p1(:)
    xyz(:,2) = p2(:)
    xyz(:,3) = p3(:)

    icon = -1

    call kgetNormalVec(xyz, 3, 1, nv, k, jcon)
    if(jcon /= 0 ) then
       write(0,*) ' error of p1,p2,p3 for kxpLand3vpB'
       stop
    endif
    call kxplp(pos(1), pos(2), pos(3), dir(1), dir(2), dir(3), &
      nv(1), nv(2), nv(3), k,  el, jcon)
    if( jcon == 0 .and. el > 0.) then
       ! xp found
       xp(:) = pos(:) + el*dir(:)
       call kinout3(xyz, 3, xp, nv, sumteta, cond)
       if( abs(sumteta) > 3.0d0 )  then
          ! inside (if outside ~0) )             
          icon = 0
       endif
       towhich = dot_product(dir(:), nv(:)) 
    else

       towhich = 0.
    endif
  end subroutine kxpLand3vpB
