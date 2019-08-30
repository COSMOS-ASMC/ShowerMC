      subroutine eppos2B(aTrack, B)
#include "ZepTrack.h"

      record /epTrack/ aTrack  ! input. charged partcle.
      record /epDirec/ B       ! ouput. magnetic strength in Tesla.
                               !        (B.x, B.y, B.z) must be given in
                               !   the local coordinate of the 
                               !   current component.
c     If you know B in the world coordinate( BinWorld), you can translate it
c    into the local coordinate value by 
c              call epw2ld(aTrack.cn, BinWorld,  B)
      end
c     ***************************
      subroutine eppos2E(aTrack, E)
#include "ZepTrack.h"

      record /epTrack/ aTrack  ! input. charged partcle.
      record /epDirec/ E       ! ouput. Electric field in  V/m
                               !        (E.x, E.y, E.z) must be given in
                               !   the local coordinate of the 
                               !   current component.
c     If you know E in the world coordinate( EinWorld), you can translate it
c    into the local coordinate value by 
c              call epw2ld(aTrack.cn, EinWorld,  E)
      end
