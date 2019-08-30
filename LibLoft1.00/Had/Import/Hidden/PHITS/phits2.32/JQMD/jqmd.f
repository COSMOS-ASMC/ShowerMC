      module  jqmd
      private
      public::  jqmdin,  ground, chname,   
     *         qmdjudge, qmdint,  caldisa,
     *         rboost, sect_a, cput_a
c     *         ichgf,  ibryf,  rboost
      contains
      include "qmd00.f"
      include "qmdcoll.f"
C      include "qmddflt.f"  ! this must be outside of module
      include "qmdgrnd.f"
      include "qmdinit.f"
      include "qmdmfld.f"
      include "mdp-uni.f" 
      include "utlnmtc.f"
      end module jqmd
