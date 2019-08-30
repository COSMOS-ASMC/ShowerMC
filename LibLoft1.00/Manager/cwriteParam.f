!             To write parameters on the error output
!
      subroutine cwriteParam(io, force)
      use modEfield
      use modMCSparam
      use modEMcontrol
      use modMuNucontrol
      implicit none
#include "ZincForNameL.h"
      integer io      !  input.  output logical dev. #. ErrorOut --> stderr
      integer force   !  input.  if non zero, Hidden parameters are written.
                      !          hidden ones are also written when Hidden=T

      write(io, Param)
      if(Hidden .or. force .ne. 0 ) then
         write(io, HParam)
      endif
      end
