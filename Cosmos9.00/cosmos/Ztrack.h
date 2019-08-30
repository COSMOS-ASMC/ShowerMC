!   #ifndef Ztrack_
!   #define Ztrack_
!    structure used when tracking a particle
!    *************************
#include  "Zcondc.h"
#include  "Zptcl.h"
#include  "Zcoord.h"
#include  "Zpos.h"
#include  "Zdirec.h"
#include  "Zmagfield.h"
!     ---------------------
      type track      ! full particle attributes in Cosmos
        sequence

          type(ptcl):: p    ! basic ptcl attributes.

!               position and time
          type(position):: pos
          real*8 t           ! time in length/beta (m)
          type(direc):: vec
          real*4 wgt         ! weight for thin sampling
           integer*2 where    ! current obsSite no. (0 is initial value)
          integer*2 asflag   ! non 0, if As has been generated from this
!                             ptcl (only for electrons)
          real*8   user      ! user use
#if LABELING > 0
          integer  label     ! put a label (1,2,...) on each particle.
                             ! There is a global label_counter which
              ! is cleared at the start of 1 event generation.
              ! it is counted up when a particle is poped up from the
              ! stack. The label_counter is given to the label of 
              ! the poped up particle. This may be needed to judge
              ! if the same particle crosses a given observation place
              ! more than once.
           integer  info      ! for each particle,  when a particle is born
              !     this is initialized to 0.  If the ptcl goes higher than
              !     380km, 1 is added. This is for AMS observation.
#endif
      end type track
!   #endif  /* Ztrack.h */


