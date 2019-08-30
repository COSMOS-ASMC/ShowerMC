      subroutine epinteLight( pj )
!     electron interaction.   
      use modSetIntInf
      implicit none
#include "Zptcl.h"
      type(pj),intent(in):: pj

      if(IntInfArray(ProcessNo)%process == "pe" ) then
           ! photo electron generation at sensor
         call epLightAtSensor(pj)
      elseif( IntInfArray(ProcessNo)%process == "rayl" ) then
                  ! Rayleigh scattering; use Xray region fomulat since
                  ! (1+cos^2)dcos
         call epcoher(pj)
      elseif( IntInfArray(ProcessNo)%process == "absorb" ) then
                  ! absorbed. nothing to do;  not push any thing
      elseif( IntInfArray(ProcessNo)%process == "wls" ) then
                   ! wave length shift
         call epLightPreWLS(pj)
      else
         write(0,*) ' light interacion=', IntInfArray(ProcessNo)%process
         write(0,*) ' not defined '
         stop
      endif
      end
