!  This is an unrealistic  example of cxsSpecial.
!     next is for before Cosmos version 9
#define BEFCOSV9
      subroutine cxsSpecial(pj, media, model)
!     Normally, in default, the hadronic intraction x-section is provided by the
!     currently active model of hadronic interaction (event generator).
!     In some case, e.g, if the user wants to use model "A" above 80 GeV/n
!     and "dpmjet3" below 80GeV/n,  for the event generation, while for hadronic
!     interaction x-section,  the user wants to use model "B" above 80 GeV/n,
!     Then, the user may give
!        IntModel='"dpmjet3", 80, "A"',
!        XsecModel='dpmjet3", 80, "B"'
!     In this case, model "A" and "B" must be pre-registered interaction models
!     (If XsecModel is not given,   XsecModel=IntModel is assumed.)
!     These functionality hasp been available in recent versions.
!     However, in some case, instead of model "B", the user may want to use peculiar
!     x-section values not providable by existing models (for some test etc).
!     In this case, the userp must give
!        XsecModel='dpmjet3", 80, "special"'
!     and must give the x-sction values through the cxsSpecial subroutine.
!
!     This file is an exampale of such a subroutine, but is non physical
!     and should not be used except fo test purpose. To use this, the user
!     may copy this to own application area and should compile together
!     with its programs (e.g  add a line (from the 1st columun)
!          #include "cxsSpeical.f"
!      to  cphook.f or  ephook.f) (
!
!     
!     This program is used only when "special" is specified in XsecModel.

!
#if defined (BEFCOSV9)
      use modXsecMedia, xmedia=>media
      implicit none
#else
      implicit none
#include "Zmedia.h"
#endif      


#include "Zcode.h"
#include "Zptcl.h"
#include "Zevhnp.h"
      
      character(*),intent(out):: model
      type(ptcl)::pj            ! input.   projectile
#if defined (BEFCOSV9)
      type(xsmedia),intent(inout):: media ! media
#else
      type(epmedia),intent(inout):: media
#endif

      real(8):: xs, sumxs
      integer:: subc
      logical needown
      integer:: i

      needown = .false.
      subc = pj%subcode
      needown = subc > 8 .and. pj%fm%p(4)/subc > 100.

      if( .not. needown) then
          ! use current active model for x-section
         call cqActiveMdl(model)
      else
         !  my own fake x-section
         do i = 1, media%noOfElem
!            xs = 0.3* media%elem(i)%A**0.5 ! very small xs
            xs = 300. * media%elem(i)%A**0.85 ! very large xs
            media%elem(i)%nsigma = xs * media%elem(i)%No
            sumxs = sumxs + media%elem(i)%nsigma 
         enddo
         media%xs = sumxs
      endif
      end
