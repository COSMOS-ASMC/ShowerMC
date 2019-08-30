      subroutine eppos2B(aTrack, B)
#include "ZepTrack.h"

       type(epTrack)::  aTrack  ! input. charged partcle.
       type(epDirec)::  B       ! ouput. magnetic strength in Tesla.
                               !        (B.x, B.y, B.z) must be given in
                               !   the local coordinate of the 
                               !   current component.
!     If you know B in the world coordinate( BinWorld), you can translate it
!    into the local coordinate value by 
!              call epw2ld(aTrack.cn, BinWorld,  B)
      end
!     ***************************
      subroutine eppos2E(aTrack, E)
#include "ZepTrack.h"

       type(epTrack)::  aTrack  ! input. charged partcle.
       type(epDirec)::  E       ! ouput. Electric field in  V/m
                               !        (E.x, E.y, E.z) must be given in
                               !   the local coordinate of the 
                               !   current component.
!     If you know E in the world coordinate( EinWorld), you can translate it
!    into the local coordinate value by 
!              call epw2ld(aTrack.cn, EinWorld,  E)
      end
