!  c4to3mom:    extract 3 momentum from 4 momentum record
!  c3to4momentum:    move    3 momentum to   4 momentum record
      subroutine c4to3mom(four, three)
      implicit none
!----      include '../Particle/Zptcl.h'
#include  "Zptcl.h"
!----      include 'Zcoord.h'
#include  "Zcoord.h"
!
      type(fmom)::four  ! input.
      type(coord)::three   ! output.
!
      three%r(1) = four%p(1)
      three%r(2) = four%p(2)
      three%r(3) = four%p(3)
      end
!     ***************************
      subroutine c3to4momentum(three, four)
      implicit none
!----      include '../Particle/Zptcl.h'
#include  "Zptcl.h"
!----      include 'Zcoord.h'
#include  "Zcoord.h"
!
      type(coord)::three   ! input
      type(fmom)::four  ! output
      
      four%p(1) = three%r(1)
      four%p(2) = three%r(2)
      four%p(3) = three%r(3)
      end
