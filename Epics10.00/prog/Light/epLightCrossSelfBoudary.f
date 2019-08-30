      subroutine epLightCrossSelfBoudary(cnx, icon)
      use modepLightPty
      implicit none
#include "ZepTrackv.h"
#include "ZepTrackp.h"
#include "Zcnfig.h"
      integer,intent(inout):: cnx ! outside is this comp. #. 
                         ! if reflection happens, will become Cn
      integer,intent(out)::icon  !  0 if not absorbed
                                 !  1 if absorbed  

      
       ! --------------------------
      integer::surfN 
      integer::jcon
      ! get surface #

      call epLightGetSurfN(Cn, Move%boundary,  surfN)
!//////////////////
!      call Lcompchk('in self', cPtyNo)
!////////////////
       ! see wrapper of the current comp. is reflector ?      
      if( comInfo( cPtyNo )%refMode(surfN) == 1 ) then
        !  current comp. is surrounded by a reflector 
         call epLightReflector(surfN, jcon)  ! jcon=0. reflected,
                ! jcon=1 abosrbed
         if(jcon == 0 ) then
              ! angle has already been set.
            cnx = Cn
            cTrack%pos = Move%Track%pos
            cTrack%cn = cnx
            Move%Cross = .false.
            icon = 0
         else
            icon = 1
            Move%Cross = .false.
         endif
      elseif( comInfo( cPtyNo )%refMode(surfN) == 0 ) then
        !  current comp. is surrounded by a grease or something
        !  very thin material not in config.
         call epLightReflecOrRefrac(cnx, surfN, jcon) ! jcon=0, reflected else
                               !      jcon = 1, refracted
          ! angle should have already been set
         if( jcon == 0 ) then
            cnx = Cn
            cTrack%pos = Move%Track%pos
            cTrack%cn = cnx
            Move%Cross =.false.
            icon= 0   ! /////////// 
         elseif( jcon ==1 ) then
            ! refraction. if at void
            if(cnx <= Det%nct) then
               icon = 0 
            else
               icon = 1
!               Move.Cross = .false.   ! must keep Cross
            endif
         endif
      elseif( comInfo( cPtyNo )%refMode(surfN) == -1 ) then
        ! no wrapper;  see if next media is "light" comp. or sensor or
        ! other.  If "light" component (&& refractive index > 0),
        ! do similar to epLightReflectionOrRefrac
         if(cnx > Det%nct) then
            icon = 1
         else
            call epLightSeeNextComp(cnx, surfN, jcon)
            if(jcon == 0 ) then
            !  reflection
               cTrack%pos = Move%Track%pos
               cTrack%cn = Cn
               cnx = Cn
               Move%Cross = .false.
               icon = 0
            elseif( jcon == 1 ) then
             ! refraction
               icon = 0
            else
               icon = 1
               Move%Cross = .false.
            endif
         endif
      else
         write(0,*) ' comInfo( cPtyNo )%refMode =',
     *      comInfo( cPtyNo )%refMode, ' strange'
         write(0,*) ' detected in epLightAtBndry'
         stop
      endif
      end
