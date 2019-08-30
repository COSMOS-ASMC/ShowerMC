#include "Zmaxdef.h"
/*
  c                this is Primary class data type definition
  c         Note that primary angle information is not here.
  c
     maxSegments; //    ! max # of segments in each primary
     maxNoOfComps; //  ! max # of components usable
     at a time 
     maxErgUnit; //       ! max # of energy unit symbols
     parameter (NoOfSymbols = 53,
*/


const int NoOfSymbols = 53; //   ! max # of primary type symbols
#ifdef MAX_SEGMENTS
const int  maxsegments = MAX_SEGMENTS;
#else
const int  maxsegments = 40;
#endif
#ifdef MAX_NO_OF_COMPS
const int maxnoofcomps = MAX_NO_OF_COMPS;
#else
const int  maxnoofcomps = 8;
#endif
const int  maxergunit=7 ;

struct component { //     ! 1 component of 1ry 
  int label;//          ! composition label number
  char symb[12]; //        ! 'p', 'gamma' etc
  char eunit[3]; //        ! 'gev' etc
  char etype[4]; //        ! 'ke/n' etc
  char diff_or_inte[1];//   ! 'd' or 'i'
  double flatterer;   //            ! di/de*e**flatterer
  double cut; //                  !  lower cut off.
  double cut2; //                 ! upper cut off.
  double energy[maxsegments+1]; //  ! segment left energy
  double flux[maxsegments+1] ; //   ! input flux
  /*
    c                                            dI/dE * E**flatterer
    c                          above: from  input table directly
    c                          below: made by subroutines
  */
  int code;
  int subcode;
  int charge;  //  ! particle code
  double togev;
  /*
    c!!!!        real*8 diff(maxSegments+1)      ! diff. flux
  */
  double norm_inte[maxsegments+1];  // ! normalized 
  //c                      integral flux > E at segment left value
  double beta[maxsegments+1]; //     ! dI/dE=(true flux)=
  //c                                            const*E**(-beta)
  int no_of_seg; //         ! no of segments given
  double  inte_value; //       !  integral flux from min. E
  double   emin, emax; //     ! min and max energy defined
  /*
    c            integer histnbin
    c            parameter (histnbin = 30)
    c            record /histgl/ comphist
  */
};
//c         *******************************************
struct primaries { //   ! 1 set of primaries.
  struct component each[maxnoofcomps];
  double cumminteflux[maxnoofcomps];
  int  no_of_comps; //    ! how many diff. compositions
  int  noofsamplings; //  ! total number of samplings including ones discarded by cutoff
  int noofsampcomp[2][maxnoofcomps];// ! 1 is for number of sampling including
  //     ! discarded ones due to cutoff. 2 is only for
  //                                                  ! employed ones.
  //c                        after a sampling of a 1ry, the following 
  //c                        is fixed.
  int label; //            ! sampled primary label
  double sampled_e; //     ! sampled energy(or momentum) as
  /*
    c                                    defined in etype. if this is
    c                                    in total energy in gev, it is the
    c                                    same as particle.fm.e
  */
  struct ptcl particle;
};




