!        cxai:  compute xai
!
!            compute critical value Xai
!           = Eg/m *  H/Hcr
      real*8  function cxai(gam, mag)
      implicit none
#include "Zglobalc.h"
#include "Zptcl.h"
#include "Zmass.h"
#include "Zmagfield.h"

      type(ptcl):: gam ! gamma ray
      type(magfield):: mag  ! magnetic field.
      real*8 bsin, cgetBsin

      bsin = cgetBsin(gam, mag)
!            Eg/m  * Bsin/Bcr      
      cxai = gam%fm%p(4) /masele * bsin/Bcr/ 2.
      end
