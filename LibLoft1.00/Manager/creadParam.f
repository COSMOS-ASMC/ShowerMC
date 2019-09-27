!             read Parameters for the cosmos simulation.
!
!     *********************************
      subroutine creadParam(io)
      use modEfield
      use modMCSparam
      use modEMcontrol
      use modMuNucontrol
      implicit none
#include  "ZincForNameL.h"

      integer io    !  logical file #.  5 is standard input.
!                      even if it is not stdin on some system.
!                      If this is not 5, the corresponding 
!     disk file must have been opened.
      write(0,*) "**********  NameList reading  ******************"
      write(0,*) "If error happens examine next:  "
      write(0,*) '  If &Param contains "PhotoProd", please omit it'
      write(0,*) '  If &HParam contains "ALateCor", it is now '  
      write(0,*) "  integer (old one is logical);"
      write(0,*) "  Ecrit, X0, MagChgDist must be removed"      
      write(0,*) "else examin the system message "
      write(0,*) " "

      if(io .eq. 5) then
         read(*, Param, end=100)
      else
         read(io, Param, end=100)
      endif
!      if(hidden) then
      if(io .eq. 5) then
         read(*, HParam, end=80)
      else
         read(io, HParam, end=80)
      endif
 80   continue
      return
 100  continue
      write(0,*) 'creadParam. Namelist assumed  on I/O unit number=',
     *   io
      call 
     * cerrorMsg('but dose not exist', 1)
      end
