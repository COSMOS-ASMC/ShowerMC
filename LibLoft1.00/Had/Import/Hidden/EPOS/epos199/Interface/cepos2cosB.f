!           cos code <--> epos code via /ptcl/
      subroutine cepos2cosB(eposid, p)
      implicit none
#include "Zptcl.h"
      integer,intent(in):: eposid  ! ptcl id of epos
      type(ptcl)::p   ! output  p%code, p%subcode, p%charge will be set

      integer code, subcode, charge

      call cepos2cos(eposid, code, subcode, charge)
      p%code = code
      p%subcode = subcode
      p%charge = charge
      end
      subroutine ccos2eposB(p, eposid)
      implicit none
#include "Zptcl.h"
      type(ptcl)::p   ! input  p%code, p%subcode, p%charge is used
      integer,intent(out):: eposid  ! ptcl id of epos

      integer code, subcode, charge
      
      code = p%code
      subcode = p%subcode
      charge = p%charge
      call ccos2epos(code, subcode, charge, eposid)
      end
      

