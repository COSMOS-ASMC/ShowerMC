!     next is for  Cosmos bofore version 9
!     This is a dummy routine to avoid "unresolved external name"
!     at link time when the user does not supply own cxsSpecial.f
!      
      subroutine cxsSpecial(pj, media, model)
!     Normally, in default, the hadronic intraction x-section is provided by the
!     currently active model of hadronic interaction (event generator).
!     In some case, e.g, if the user wants to use model "A" above 80 GeV/n
!     and "dpmjet3" below 80GeV/n,  for the event generation, while for hadronic
!     interaction x-section,  the user want to use model "B" above 80 GeV/n,
!     Then, the user may give
!        IntModel='"dpmjet3", 80, "A"',
!        XsecModel='dpmjet3", 80, "B"'
!     In this case, model "A" and "B" must be pre-registered interaction models
!     (If XsecModel is not given,   XsecModel=IntModel is assumed.)
!     These functionality has been available in recent versions.
!     However, in some case, instead of model "B", the user may want to use peculiar
!     x-section values not providable by existing models (for some test etc).
!     In this case, the user must give
!        XsecModel='dpmjet3", 80, "special"'
!     and must give the x-sction values through the cxsSpecial subroutine.
!
!     This file contains a kind of dummy prototype of such a subroutine;
!     although it is dummy,
!     it must exist and included in the library,  even if the user
!     does  not give "special" to "B".  This is because there is a calling
!     sequence to cxsSpecial in cGetXsec.f; without this, error will happen
!     when the user mekes an application.
!
!     Meanwhile if the user gives "special", then the user must  provide
!     cxsSpecial.f  to replace this dummy,  otherwise, program 
!     exectution will stop. 
!     In summary,
!     A) If the user does not use  "special" in XsecModel, the user need not
!        prepaare own cxsSpecial.f.  Even if the user prepares own routine,
!        it is neglected.
!     B) If the user specifies "special" in XsecModel, the user must prepare
!        own  cxsSpecial.f and must give x-section values in that subroutine
!        Otherwise, rogram execution will stop.   
!       
!     Inside the user's  own cxsSpecial,  the user may skip calculatiion of
!     x-section depending on projectile, target and projectile energy.
!     In such a case, the user may employ values by some registered model. 
!
      implicit none

#include "Zmedia.h"
#include "Zcode.h"
#include "Zptcl.h"
#include "Zevhnp.h"
      

      type(ptcl)::pj            ! input.   projectile
      type(epmedia),intent(inout):: media

      character(*),intent(out):: model
!         next check is only for prototype.  don't use in your own cxsSpecial
      if( model == "special" ) then
!        "special" case should not be called for this program;
!        the user must prepare own cxsSpecial
         write(0,*)
     *   ' The dummy cxsSpecial is called for XsecModel="special"'
         write(0,*)
     *   ' must prepare own cxsSpecial and  compile with your apps.'
         write(9,*) ' see  Cosmos/UserHook/cxsSpecial.f' 
         stop
      endif

!     In your own cxsSpecial,  you must give media%xs, media%elem(:)%nsigmas.
!     Let xs for pj and the i-th target element (with mass number
!     media%elem(i)%A)  be XSi,  then
!       media%elem(i)%nsigma= XSi*media%elem(i)%No
!     and
!        media%xs = sum( media%elem(1:n)%nsigma )
!     where  n = media%noOfElem      

!      *************************    
!     If pj, media, pj%fm%p(4)/n are not suited for "special", next may be
!     activated,      
!!!!      call cqActiveMdl(model)
!     to use ActiveMdl.
!     (You may also set "model" with another known model without calling cqA...

!        Then you may simply exit the subrouine.
!      ******************a      
      end
