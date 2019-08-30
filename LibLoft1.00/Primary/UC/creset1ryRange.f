      subroutine creset1ryRange(cmpn, E1, E2)
!     this subroutine may be called when the user wants to
!     reset the 1ry energy range during the execution (normally
!     at start time)      
!     If called, those given in the table is overriden.
#include "Zprimary.h"
#include "Zprimaryv.h"
      implicit none

      integer,intent(in):: cmpn ! 1 ry component #
      ! 1 to Prim%no_of_comps
      real(8),intent(in):: E1, E2 ! min and max of
              !   1ry energy to be smapled:  in the same unit as
              !   gvien in the input table. must be E1<E2.
              ! If E1<min in the table, the  min is used.
              ! If E2>max in the table, the max is used.


      integer::nseg
      
      if( E2 <= E1 ) then
         write(0,*) ' E1=',E1, ' is not <  E2=',E2
         write(0,*) ' error in creset1ryRange '
         stop
      endif
      if( cmpn < 1 .or. cmpn > Prim%no_of_comps) then
         write(0,*) ' cmpn given to creset1ryRange=',cmpn,
     *        ' must be in 1 to', Prim%no_of_comps
         stop
      endif
      if( E1< Prim%each(cmpn)%energy(1) ) then
         Prim%each(cmpn)%cut = Prim%each(cmpn)%energy(1)
      else
         Prim%each(cmpn)%cut = E1
      endif

      nseg =  Prim%each(cmpn)%no_of_seg
      if( E2> Prim%each(cmpn)%energy(nseg) ) then
         Prim%each(cmpn)%cut2 = Prim%each(cmpn)%energy(nseg)
      else
         Prim%each(cmpn)%cut2 = E2
      endif
      call cprocPrimDt(Prim)
      end   subroutine creset1ryRange
      
      
