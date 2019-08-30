      subroutine setpcode( p,  code)
      implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
       type(ptcl):: p
       integer code
!                set parent particle code by using last bit
!                      bit 2    bit 1
!                 pi    0         0
!                 K     0         1
!                 p     1         0
!             pbar/nbar 1         1 
       if( p.code .eq. 4 ) then
          call krsetbit(code, 1)
          call krsetbit(code, 2)
       elseif(p.code .eq. 5) then
          call ksetbit(code, 1)
          call krsetbit(code, 2)
       elseif( (p.code .eq. 6 .and. p.charge .eq. -1 ) .or.
     *         (p.code .eq. 6 .and. p.charge .eq. 0 .and.
     *          p.subcode .eq. antip ) ) then
          call ksetbit(code, 1)
          call ksetbit(code, 2)
       else
          call krsetbit(code, 1)
          call ksetbit(code, 2)
       endif
       end
      subroutine getpcode(code, codex)
      implicit none
      integer code  ! input.
      integer codex ! output.  4: pi, 5: K,  7: anti p/n 6: p/n and others
      logical kbitest
      if(.not.  kbitest(code, 2) ) then
         if( .not. kbitest(code, 1) ) then
            codex = 4
         else
            codex = 5
         endif
      else
         if(.not. kbitest(code, 1) ) then
            codex = 6
         else
            codex = 7
         endif
      endif
      end
