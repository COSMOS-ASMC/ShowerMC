!   ###include  "Zmaxdef.h"
!   ###include  "Ztrack.h"
          integer Max_stack_size
          parameter (
#ifdef MAX_STACK
     *    Max_stack_size = MAX_STACK
#else
     *    Max_stack_size =15000
#endif
     *    )
          integer Stack_pos
      type(track):: Stack(Max_stack_size)
      common /Zstack/ Stack, Stack_pos
