#include  "ZepMaxdef.h"
          integer MaxStackSize
          parameter (MaxStackSize = EPMAX_STACK )

          integer Stack_pos,  DiskPos, DdiskNo
       type(epTrack)::  Stack(MaxStackSize)
          common /StackC/ Stack, Stack_pos, DiskPos, DdiskNo 

    !     DdiskNo   ! direct access disk #.  0 means, it is not
                    ! yet opened. 

