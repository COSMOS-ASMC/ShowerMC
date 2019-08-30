      subroutine epMinMaxQuenchComp(compno1, compno2)
!        This is to be inside of epmodUI.f90 but due to
!        include  files, separated.
      implicit none
#include "ZmediaLoft.h"
#include "ZepTrackv.h"      
#include "Zcnfig.h"
      integer,intent(out):: compno1, compno2 ! min/max component #
         !  for which
         ! quenching effect must be considered.
         ! if there is no, 0 is returned. 
         ! Such scintillators are assumed to be put at the
         ! top part of the detetor, so their component # 
         ! would be comparatively small.

      integer::i, mn
      integer:: modif

      compno1 = 0
      compno2 = 0
      do i = Det%nct, 1, -1
         modif = Det%cmp(i)%modifier
         mn = Det%Cn2media(i)
         if(modif > 0 .or. 
     *      Media(mn)%Birks /= " ") then
            compno2 =  i
            exit
         endif
      enddo
      do i = 1, Det%nct
         modif = Det%cmp(i)%modifier
         mn = Det%Cn2media(i)
         if(modif > 0 .or. 
     *      Media(mn)%Birks /= " ") then
            compno1 =  i
            exit
         endif
      enddo

      end  subroutine epMinMaxQuenchComp

         
