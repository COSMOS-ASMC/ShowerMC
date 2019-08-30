      subroutine epMuBrEcheck(Emu, media)
      implicit none
#include "Zmedia.h"
      real(8),intent(in):: Emu
      type(epmedia),intent(in):: media

      if( Emu < media%cnst%muBrEmin*0.75 ) then
         write(0,*)' Emu=',Emu, ' for muon brems is too small '
         write(0,*) ' must be > ~0.8*',media%cnst%muBrEmin
         stop
      endif
      end   subroutine epMuBrEcheck
      subroutine epMuPrEcheck(Emu, media)
      implicit none
#include "Zmedia.h"
      real(8),intent(in):: Emu
      type(epmedia),intent(in):: media

      if( Emu < media%cnst%muPrEmin*0.75 ) then
         write(0,*)' Emu=',Emu, ' for muon pair cre. is too small '
         write(0,*) ' must be > ~0.8*',media%cnst%muPrEmin
         stop
      endif
      end   subroutine epMuPrEcheck
      subroutine epMuNiEcheck(Emu, media)
      implicit none
#include "Zmedia.h"
      real(8),intent(in):: Emu
      type(epmedia),intent(in):: media

      if( Emu < media%cnst%muNEmin ) then
         write(0,*)' Emu=',Emu, ' for muon N.I is too small '
         write(0,*) ' must be > ',media%cnst%muNEmin
         stop
      endif
      end   subroutine epMuNiEcheck
