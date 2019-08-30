      module  modcMCS
      use modTPXS
      use modDCS
      use modMCS
      implicit none

      type(MCSconst),target,allocatable:: MCSnega(:)
      type(MCSconst),target,allocatable:: MCSposi(:)
!   currently assigned TPXS,DCS,MCS const (may be  for e-
!   e+ for a given media)
!   Must be fixed thru cfixMixedConst after starting
!   simulation and when a media and/or e+/e- is changed.
      type(TPXSconst),pointer::TPXSnow  
      type(DCSconst),pointer::DCSnow
      type(MCSconst),pointer::MCSnow
      end module modcMCS


      subroutine cmkMCSconst
!  use modXsecMedia
      use modcMCS
      implicit none
#include "ZmediaLoft.h"  
      integer::i
      call ciniSmpTab  ! set uarray 
      allocate( MCSnega(NoOfMedia) )  
      allocate( MCSposi(NoOfMedia) )
      do i = 1, NoOfMedia 
         call cMCSconstForMedia(media(i), -1,  MCSnega(i))
         call cMCSconstForMedia(media(i),  1,  MCSposi(i))
      enddo
      end subroutine cmkMCSconst


!!!!! In normal situation, next are to be placed
!  in cMCSconstForMedia.f but it is not f90 so
!  these are placed here.
      subroutine creadMCSTab(io, aMCS)
!   read one set of MCS table  for e- or e+
!   MCS data file is organized 
!   E  muc  lh  ls1 ls2 ---------  e-
!    .... data nEneg lines
!   blank
!   header for sampleV
!   usize x nEneg lines
!   blank
!   head for e+ (E muc ...)
!   ...
      use  modTPXS
      use  modcMCS
      implicit none
      integer,intent(in):: io ! opened file logical #
                      ! current head should be at the first
                      ! data position for e- or e+
      type(MCSconst),intent(out):: aMCS

      integer:: i, j
      real(8):: temp

      read(io, *) aMCS%minNon0mucEindex, aMCS%minNon0mucE
      read(io, *)  ! skep 1 line
      do i = 1, nEneg
         read(io, *)  
     *         temp, aMCS%muc(i), aMCS%lambdah(i),  
     *         aMCS%lambdas1(i),  aMCS%lambdas2(i)
         if( abs(temp-KEele(i))/ temp > 1.d-4) then
            write(0,*) ' table energy and KEele differ'
            write(0,*) ' in creadMCSTab'
            write(0,*) ' table : =', temp
            write(0,*) ' KEele : =', KEele(i)
            stop
         endif
      enddo
      read(io,*)
      read(io,*)
      do j = 1, nEneg
         read(io, *) (aMCS%sampleV(i, j), i=1, usize) 
      enddo
      end subroutine creadMCSTab

      subroutine cfixMixedConst(mn, pm)
      use  modcTPXS
      use  modcDCS
      use  modcMCS

      implicit none
      integer,intent(in):: mn  ! media #
      integer,intent(in):: pm  ! -1 or +1 for e-/e=
      
      if(pm == -1) then
         TPXSnow => TPXSnega(mn)
         DCSnow => DCSnega(mn)
         MCSnow => MCSnega(mn)
      elseif( pm == 1) then
         TPXSnow => TPXSposi(mn)
         DCSnow => DCSposi(mn)
         MCSnow => MCSposi(mn)
      else
         write(0,*) 'In cfixMixedCost:  pm=',pm,' invalid'
         write(0,*) 'mn=',mn
         stop
      endif
      end subroutine cfixMixedConst



