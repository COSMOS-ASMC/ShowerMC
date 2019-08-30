!              print primary information
      subroutine cprintPrim(out)
      implicit none

#include  "Zptcl.h"
#include  "Zprimary.h"
#include  "Zprimaryv.h"
      integer out  ! input.  output logical device #

      call cprintPrim0(Prim, out)
      end
!
!              print primary information
      subroutine cprintPrim0(prm, out)
      implicit none
#include  "Zptcl.h"
#include  "Zprimary.h"
#include  "Zprimaryv.h"
      type(primaries):: prm     ! input. primary data
      integer out                ! input. output logical dev. #
!
!
      integer i
!
      if(out .eq. 5) then
         write(*,'("# -----------------# of Component defined=",
     *           i2,"----------------------")')  prm%no_of_comps
         write(*,
     *     '("# Seq.#  Code  Comp.    E_ Unit  E_type Diff/Integ",
     *     "  Emin    Emax   # of seg.  inte%flux")')
      else
         write(out,'("# -----------------# of Component defined=",
     *           i2,"----------------------")')  prm%no_of_comps
         write(out,
     *     '("# Seq.#  Code  Comp.    E_ Unit  E_type Diff/Integ",
     *        "  Emin    Emax   # of seg. inte%flux ")')
      endif
      do i = 1, prm%no_of_comps
         call cprintEachCmp(prm%each(i), out)
      enddo
      end
      subroutine cprintEachCmp(each, out)
      implicit none

#include  "Zptcl.h"
#include  "Zprimary.h"
      type(component):: each
      integer out
!
      if(out .eq. 5) then
         write(*, '("#",i4,i6, 2x, a10, a7, a9, a8, 5x, 2g9.3,i4,
     *   1x, g9.3)')
     *   each%label, each%code, each%symb, each%eunit, each%etype,
     *   each%diff_or_inte,
     *   each%emin, each%emax, each%no_of_seg, each%inte_value
      else
         write(out, '("#",i4,i6, 2x, a10, a7, a9, a8, 5x, 2g9.3,i4,
     *    1x,  g9.3)')
     *   each%label, each%code, each%symb, each%eunit, each%etype,
     *   each%diff_or_inte,
     *   each%emin, each%emax, each%no_of_seg, each%inte_value
      endif
      end
