      subroutine cinteElec(pj)
      use modSetIntInf
      implicit none
#include  "Zptcl.h"
! #include  "Ztrack.h"
!     #include  "Ztrackv.h"

      type(ptcl),intent(inout):: pj  ! in some rare case, out 
      
         character*70 msg
!
         if(IntInfArray(ProcessNo)%process .eq. 'brems') then
            call  epbrems(pj)
         elseif(IntInfArray(ProcessNo)%process .eq. 'mscat') then
            call epknoc( pj )
         elseif(IntInfArray(ProcessNo)%process .eq. 'bscat') then
            call epknoc( pj )
         elseif(IntInfArray(ProcessNo)%process .eq. 'anihi') then
            call epanih( pj )
         elseif(IntInfArray(ProcessNo)%process .eq. 'mbrem') then
            call epsync( pj ) 
         else

              write(msg, *) ' process#=', ProcessNo,
     *       'and process for electron=',
     *        IntInfArray(ProcessNo)%process, ' undef'
              call cerrorMsg(msg, 0)
         endif
         end
!     cbrems, cknockon etc are now not used and removed ; now stored in Obso/cbremsETC.f
