module modShiftInciPos
  real(8),save:: eps=0.d0
  real(8),save:: minv1(3), maxv1(3), minv2(3), maxv2(3)
  character(len=12),save::ShiftInciPos=" " 
  character(len=2),save:: plane="xy"
end  module modShiftInciPos
subroutine epShiftInciPos( Trackpos )
  use modShiftInciPos
  implicit none
!! #include "ZepTrack.h"  
!!      record /epTrack/ aTrack  ! incident track with world coord.
!!         next is Track.pos. so (1) is pos.x
  real(8),intent(inout)::Trackpos(3)
  
  logical,save:: first=.true.
  integer::n
  character(12),save:: wstruc
  real(8):: org(3), abc(3)
  real(8):: u1, u2, s1, s2
  integer,save:: i, j, nf
  character(len=9):: field(2)
  
  if( first ) then
!     get eps and plane  from data like 
!           1.e-5           ===> eps=1e-5 plane="xy" 
!       or
!           1.e-5 yz        ===> eps=1e-5 plane="yz" 
     call ksplit(ShiftInciPos, 9, 2, field, nf)
     read(field(1), *) eps
     if( nf > 1 ) then
        plane =trim(field(2))
     endif
     
     call epqworld(n, wstruc)
     if( n == 0 ) then
        ! box
        call epenvlpAll
        ! be carefule org and ,abc are defined as epPos
        call epqcnf(org, abc) 
        wstruc = 'box_w'
     else
        call epqwcoord(org, abc)
        if( wstruc /= "box_w") then
           if( wstruc /= "sphere_w" ) then
              write(0,*) 'For epShiftInciPos, world must be'
              write(0,*) 'box_w(better) or sphere_w  but world is', &
                 wstruc
              stop
           endif
        endif
     endif
     if( plane == "xy") then
        i = 1
        j = 2
     elseif( plane == "yz" ) then
        i = 2
        j = 3
     elseif( plane == "zx" ) then
        i = 3
        j = 1
     else
        write(0,*) ' in epShiftInciPos'
        write(0,*) ' plane=',plane, ' invalid'
        stop
     endif
     minv1(i) = org(i) + eps*2
     maxv1(i) = org(i) + abc(i) - eps*2
     minv2(j) = org(j) + eps*2
     maxv2(j) = org(j) + abc(j) - eps*2
     first = .false.
  endif
  if( wstruc == 'box_w' ) then
     call rndc(u1)
     call rndc(u2)
     if(u1 < 0.5d0) then
        s1 = 1.0d0
     else
        s1 = -1.0d0
     endif
     if( u2 < 0.5d0) then
        s2 = 1.0d0
     else
        s2 = -1.0d0
     endif
     if(Trackpos(i) > minv1(i) .and. Trackpos(i) < maxv1(i)) then
        Trackpos(i) = Trackpos(i) + s1*eps
     elseif( Trackpos(i) <= minv1(1) ) then
        Trackpos(i) = Trackpos(i) + eps
     else
        Trackpos(i) = Trackpos(i) - eps
     endif

     if(Trackpos(j) > minv2(j) .and. Trackpos(j) < maxv2(j)) then
        Trackpos(j) = Trackpos(j) + s2*eps
     elseif( Trackpos(j) <= minv2(j) ) then
        Trackpos(j) = Trackpos(j) + eps
     else
        Trackpos(j) = Trackpos(j) - eps
     endif
  else
     if( Trackpos(i) < org(i) ) then
        Trackpos(i) = Trackpos(i) + eps
     else
        Trackpos(i) = Trackpos(i) - eps
     endif

     if(Trackpos(j) < org(j) ) then
        Trackpos(j) = Trackpos(j) + eps
     else
        Trackpos(j) = Trackpos(j) - eps
     endif
  endif
end   subroutine epShiftInciPos

