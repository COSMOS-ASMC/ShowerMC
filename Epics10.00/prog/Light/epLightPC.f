!       photon counter at the exit of a "light" component
!            or at the entrance of a  "light" component

      subroutine epLightPC( info )
      use modepLightPty
      use modepLightCounter
      implicit none
#include  "ZepTrackv.h"
      integer,intent(in):: info  !  0 : exiting to void
                                 ! <0 : exiting to |info| comp.
                                 ! >0 : entering from comp # info  to Cn
                                 !  

      integer:: sub             ! scinti or Ceren or ..

      integer::surfN

      if(info <=0 ) then
         call epLightGetSurfN(Cn, Move%boundary,  surfN)
      else
         call epLightGetSurfN(Cn, Move%Track%pos,  surfN)
      endif

      sub = Move%Track%p%subcode
      if(info > 0 ) then
         call epLightAllocPC( entering )
         entering( cLcompNo )%pc(sub, surfN)  =
     *        entering( cLcompNo )%pc(sub, surfN)   + 
     *              Move%Track%wgt
      else
         call epLightAllocPC( exiting )
         exiting( cLcompNo )%pc(sub, surfN)  =
     *   exiting( cLcompNo )%pc(sub, surfN)   + 
     *        Move%Track%wgt
      endif
      end

      subroutine epLightAllocPC( pcc )
!          we cannot use exiting and  entering as argument
!        so  we make EntPC and ExtPC separately
      use modepLightCounter
      implicit none
#include  "ZepMaxdef.h"
#include  "ZepTrackv.h"

      type(photonc),intent(inout)::pcc(*)
      character(len=MAX_STRUCCHR)::struc, basename  ! v9.164


!///////////////
      if( Move%Track%cn == 0 ) then
         write(0,*)  ' MOve cn=0 iLightCounter '
      endif
!/////////////
      call epqstruc(Move%Track%cn, struc)
      call epGetBaseStrucName(struc, basename)
      if( .not. allocated( pcc(cLcompNo)%pc )) then
         if( basename == 'box') then
            allocate(  pcc(cLcompNo)%pc(1:2,1:6)  )
            allocate(  pcc(cLcompNo)%pct(1:6) )
            pcc(cLcompNo)%pc(:,:) = 0.
         elseif( basename == 'octagon') then
            allocate(  pcc(cLcompNo)%pc(1:2,1:10)  )
            allocate(  pcc(cLcompNo)%pct(1:10) )
            pcc(cLcompNo)%pc(:,:) = 0.
         elseif( basename == 'cyl') then 
            allocate( pcc(cLcompNo)%pc(1:2,1:3)  )
            allocate( pcc(cLcompNo)%pct(1:3) )
            pcc(cLcompNo)%pc(:,:) = 0.
         elseif( basename == 'ecyl') then 
            allocate( pcc(cLcompNo)%pc(1:2,1:3)  )
            allocate( pcc(cLcompNo)%pct(1:3) )
            pcc(cLcompNo)%pc(:,:) =0.
         elseif( basename == 'pipe') then 
            allocate( pcc(cLcompNo)%pc(1:2,1:4)  )
            allocate( pcc(cLcompNo)%pct(1:4) )
            pcc(cLcompNo)%pc(:,:) = 0.
         elseif( basename == 'prism') then
            allocate(  pcc(cLcompNo)%pc(1:2,1:5)  )
            allocate(  pcc(cLcompNo)%pct(1:5) )
            pcc(cLcompNo)%pc(:,:) = 0.
         else
            write(0,*) ' structure =', struc, ' is ',
     *           ' not yet supported for light transport '
         endif
      endif
      end
