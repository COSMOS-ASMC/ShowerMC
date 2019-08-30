/*
c    structure used when tracking a particle
c    *************************
*/
#include  "Zcondc.h"
#include  "Zptcl.h"
#include  "Zcoord.h"
#include  "Zpos.h"
#include  "Zdirec.h"
#include  "Zmagfield.h"

/*
c     ---------------------
      structure /track/      ! full particle attributes in Cosmos

          record /ptcl/ p    ! basic ptcl attributes

c               position and time
          record /position/ pos
          real*8 t           ! time in length/beta (m)
          record /direc/ vec
          real*4 wgt         ! weight for thin sampling
	  integer*2 where    ! current obsSite no. (0 is initial value)
          integer*2 asflag   ! non 0, if As has been generated from this
c                             ptcl (only for electrons)
          real*8  user       ! user can use this 
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
      end structure   

*/

struct track {
  struct ptcl p;
  struct position pos;
  double  t;
  struct direc vec;
  float  wgt;
  short int where;
  short int asflag;
#if LABELING > 0
  int   label;
  int   info;
#endif
};




