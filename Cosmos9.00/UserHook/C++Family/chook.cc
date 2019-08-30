#include "Zdef.h"
#include "Zcode.h"
#include "Zmanagerp.h"
#include "Ztrack.h"
#include "Ztrackv.h"
#include "Ztrackp.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"
#include <iostream.h>
#include <math.h> 
#include "Zprivate.h"
extern "C" {
  extern void cwriteparam_(int *, int *);
  extern void cprintprim_(int *);
  extern void cprintobs_(int *);
  extern void cqincident_(track *, coord *);
  extern void cqfirstid_( double *);
  extern void cecent2sph_(coord *, coord *);
  extern double cgetbsin_( ptcl *, magfield * );

  //  user defined routines follow:
  //  since these are defined in C++ only, we can use
  //  upper and lower cases mixed.

  int cdropBig( float * x, float * y, int n, float rmax) {
    int i, no;
    double r;
    no = 0;
    // we use the same counter (normally from 1 and say n)
    // (i.e, not 0 to n-1). and adjust the index adding -1
    for( i = 1; i <= n; i++ ) {
      r = sqrt( pow(x[i-1],2) + pow(y[i-1],2));
      if( r < rmax ) {
	no++;
	x[no-1] = x[i];
	y[no-1] = y[i];
      }
    }
    return (no);
  }
  float cgetave( float *x, int n) {
    int i;
    float ave = 0.;
    for(i = 1; i<= n; i++ ) {
      ave += x[i-1];
    }
    if( n >  1 ) ave /= n;
    return (ave);
  }
  void ccount( int nc[][nth], track *aTrack) {
    for(int i=1; i<=com.ntha; i++) {
      if( aTrack->p.fm.p[3] - aTrack->p.mass < com.eth[i-1] ) break;
      nc[aTrack->where-1][i-1]++;
    }
  }


  //  Cosmos interface routines

  void time_(int *xxx){
    *xxx = 1;
  }

/*
 *************************************** hook for Beginning of a Run
 * At this moment, all (system-level) initialization for this run
 * has been ended.  After this routine is executed, the system goes into the
 * event creation loop.
 *
*/

  void chookbgrun_(){
    //            namelist output
    int itemp=0;
    cwriteparam_( &ErrorOut, &itemp);
    //            primary information
    cprintprim_( &ErrorOut);
    //            observation level information
    cprintobs_( &ErrorOut);



    float temp[nth]={0.3e-3, 0.5e-3, 1.e-3, 2.e-3, 5.e-3, 1.e-2,
		     2.e-2, 5.e-2, 0.1, 0.2, 0.5, 1};
    for(int i=1; i<=nth; i++)  com.eth[i-1] = temp[i-1];
    com.ntha = 1; //use only 1 threshold for output
  }
  /*     
	 c*********************************** hook for Beginning of  1 event
	 c*  All system-level initialization for 1 event generation has been
	 c*  eneded at this moment.
	 c*  After this is executed, event generation starts.
	 c*
  */
  using namespace std;
  void  chookbgevent_(){
    struct track inci;
    struct coord angle;

    int i, j;
    for(i=1; i <= nl; i++) {
      for(j=1; j <= com.ntha; j++) {
	com.ng[i-1][j-1] = 0;
	com.ne[i-1][j-1] = 0;
	com.nmu[i-1][j-1] = 0;
      }
    }
    com.nnn = 0;

    cqincident_( &inci, &angle);
    cout << "p "
         <<  EventNo+1 <<
      " " <<  inci.p.code  <<
      " " <<  inci.p.subcode <<
      " " << inci.vec.coszenith <<
      " " <<  inci.p.fm.p[3]<< endl;


  }
  /*     ************************************ hook for observation
     *  One particle information is brought here by the system.
     *  All information of the particle is in aTrack
  */
  void chookobs_(track *atrack, int *id){
    /*
      c
      c     Note that every real variable is in double  precision so
      c     that you may output it in sigle precision to save the memory.
      c     In some cases it is essential to put it in sigle (say,
      c     for gnuplot).
      c 

      integer id  ! input.  1 ==> aTrack is going out from
      c                                 outer boundery.
      c                           2 ==> reached at an observation level
      c                           3 ==> reached at inner boundery.
      record /track/ aTrack
      c
      c     For id =2, you need not output the z value, because it is always
      c     0 (within the computational accuracy).
    */
    struct track inci;
    struct coord angle;

    cqincident_( &inci, &angle); // not used here
    int iij = atrack->p.code;
    switch (iij){
    case kelec:
      ccount( com.ne, atrack);
      com.nnn++;
      com.x[com.nnn-1] = atrack->pos.xyz.r[0];
      com.y[com.nnn-1] = atrack->pos.xyz.r[1];
      com.erg[com.nnn-1] = atrack->p.fm.p[3];
      break;
    case kphoton:
      ccount( com.ng, atrack);
      com.nnn++;
      com.x[com.nnn-1] = atrack->pos.xyz.r[0];
      com.y[com.nnn-1] = atrack->pos.xyz.r[1];
      com.erg[com.nnn-1] = atrack->p.fm.p[3];
      break;
    case kmuon:
      ccount( com.nmu, atrack);
      break;
    case kpion:
    case kkaon:
    case knuc:
      //      if( atrack->p.charge != 0) {
      ccount( com.nh, atrack);
	//      }
    default:;
    }

    //    /*

    if( *id == 2) {
      //      output typical quantities.
      cout << "o  " << atrack->where //observation level
	   << " " << atrack->p.code //  ptcl code.
	   << " " << atrack->p.charge //  charge, 
	// << " " << atrack->t   //  relateive arrival time in nsec (not sec).
                        	//  if timestructure is f, nonsense.
	   << " " << atrack->p.fm.p[3] //  ! total energy in gev.
	   << " " << atrack->pos.xyz.r[0]  //  x, y
	   << " " << atrack->pos.xyz.r[1] 
	   << " " << atrack->vec.w.r[0]   //direc. cos.x in the current detector system.
		
	   << " " << atrack->vec.w.r[1] // direc. cos.y
	   << " " << atrack->vec.w.r[2] // direc. cos.z
	   << " " << atrack->vec.coszenith;  //cos of zenith angle
      cout << endl;
    }
    //    */

    /*
        you may need in some case other information such as
      atrack->p.subcode   // sub code of the particle integer*2
      atrack->p.mass      // mass 
      atrack->wgt         // weight of the particle (may not be 1. if
        // thinsampling =t)
      atrack->p.fm.p[0]   // momentum x component.  note. momentum is
        // given in the  earth xyz system.
      atrack->p.fm.p[1]     //          y
      atrack->p.fm.p[2]     //          z
    */
  }

  /*
    c    *********************************** hook for end of 1 event
    c    * At this moment, 1 event generation has been ended.
    c    *
  */
  void chookenevent_(){
    track inci;
    coord angle;
    coord tetafai;
    int i, j, nnew;
    double  fdepth;
    float avex, avey, sume;
    double sumx;

    cqincident_( &inci, &angle);


    if( ObserveAS ) {
      cqfirstid_( &fdepth );   // first int. depth in kg/m^2
      fdepth *=0.1;    // to g/cm^2
      angle.r[0] =-  angle.r[0] ; // for our normal sense
      angle.r[1] =-  angle.r[1] ;
      angle.r[2] =-  angle.r[2] ;
      cecent2sph_( &angle, &tetafai);  //to polar angle
      double teta = tetafai.r[0];
      double fai =  tetafai.r[1];
      if( fai < 0.0) fai  = 360.0 - fai;
      // for GEOMAG effect on photon primary
      double bsin = cgetbsin_( &inci.p, &Mag_ ); 

      double sumsize = 0.;
      for( j = 1; j <= com.ntha; j++) {
	for( i = 1; i<= NoOfASSites; i++) {
	  sumsize += ASObsSites[i-1].esize;
	  cout << "a ";
	  cout << ASObsSites[i-1].pos.depth/10.0
	       << " " << ASObsSites[i-1].esize
	       << " " << ASObsSites[i-1].age
	       << " " << fdepth
	       << " " << bsin 
	       << " " << sumsize
	       << " " << teta
	       << " " << fai
	       << " " << com.ng[i-1][j-1]
	       << " " << com.ne[i-1][j-1]
	       << " " << com.nmu[i-1][j-1]
	       << " " << com.nh[i-1][j-1]
	       << " " << com.eth[j-1]<< endl;
	}
      }
    }
    // iterate 3 times to discard too distant partilces
    for(j = 1; j<= 3; j++) {
      avex = cgetave( com.x, com.nnn);
      avey = cgetave( com.y, com.nnn);
      sume = cgetave( com.erg, com.nnn);
      sume *= com.nnn;

      if (com.nnn > 4) {
	sumx = 0.0;
	for(i=1; i<=com.nnn; i++) {
	  com.x[i-1] -=  avex;
	  com.y[i-1] -=  avey;
	  sumx += sqrt(pow(com.x[i-1],2) + pow(com.y[i-1], 2));
	}
	nnew = cdropBig( com.x, com.y, com.nnn, (float)1.5e-1);
	if( nnew == com.nnn ) break;
	com.nnn = nnew;
      }
      else break;
    }
    //  why ?
    if( com.nnn < 0) {
      cqincident_( &inci, &angle);
      cqfirstid_( &fdepth );
      fdepth *= 0.1;

      cout << avex*100.0 << " " << avey*100.0 
	   << " " << sumx/com.nnn*100.0 
	   << " " << com.nnn  << " " << sume/1000. 
	   << " " << inci.p.code << " " << inci.p.fm.p[3]/1000.
	   << " " << -angle.r[2] << endl;
    }
  }
  /*
    c     ********************************* hook for end of a run
    c     *  all events have been created or time lacks
    c     *
  */
  void chookenrun_(){
    cout << "end of run" << endl;
  }
  /*
    c     ********************************* hook for trace
    c     *  This is called only when Trace > 60
    c     *  User should manage the trace information here.
    c     *  If you use this, you may need some output for trace
    c     *  at the beginning of 1 event generatio and at the end of  1 event
    c     *  generation so that you can identfy each event.
    c     *
    c     *
  */
  void  chooktrace_(){

      /*
	c
	c    Every time a particle is moved in the atmosphere, this routine is called,
	c    if trace > 60. 
	c         For a one track segment,
	c     TrackBefMove  has  track information at the beginning of the segment.
	c     MoveTrack    has   track information at the end of the segment.
	c   
	c     You can know the  information a track contains in the 
	c     chookObs routine. (Note however, no conversion of coordinate
	c     has been done.  The values are in the Earth xyz system.)
	c     Besides quantities explained there, you can use, for a  given 'track'
	c
	c     atrack.pos.xyz.x, atrack.pos.xyz.y, atrack.pos.xyz.z    (x,y.z)
	c     atrack.pos.radiallen   (distance from the center of the earth)
	c     atrack.pos.depth       (vertical depth)
	c     atrack.pos.height      (vertical heigth from sea level)  
	c
      */
    float h1,  h2;
    //   Fortran;  index 0, 1, 2,    NoOfSites, NoOfSintes+1
    //                   LB                      UB
    //    C;             0  1, 2,       ...         
    //    so, the index can be the same for the both
    h1 = TrackBefMove.pos.height -
      ObsSites[NoOfSites].pos.height;
    h2 = MovedTrack.pos.height -
      ObsSites[NoOfSites].pos.height;
  }
  /*
    c     ********************* this is the hook called when
    c       an electron made an interaction.
    c
  */
  void chookeint_(int * never) {
    /*
      integer never   ! input & output
      c         don't make never = 1, if you want to get
      c         information after an electron made interaction
      c         if this is made non zero, this routine will never be called.
      c
      c   MovedTrack is the electron that made interaction
      c   Pwork contains produced particles.
      c   Nproduced has the number of particles in Pwork
      c   IntInfArray(ProcessNo) contains the type of interaction
      c
      c        default setting
    */
    *never = 1;
    /*
	c
	c        IntInfArray(ProcessNo).process will have one of
	c       'brems', 'mscat', 'bscat',or  'anihi'
	c
    */
  }
  /*
    c     ********************* this is the hook called when
    c       a gamma ray made an interaction.
    c
  */
  void chookgint_(int * never) {
    /*
      integer never   ! input & output
      c         don't make never = 1, if you want to get
      c         information after a gamma ray made interaction
      c         if this is made non zero, this routine will never be called.
      c
      c   MovedTrack is the gamma that made interaction
      c   Pwork contains produced particles.
      c   Nproduced has the number of particles in Pwork
      c   IntInfArray(ProcessNo) contains the type of interaction
      c
      c        default setting
    */
    *never = 1;
    /*
      c         IntInfArray(ProcessNo).process will have one of
      c        'pair', 'comp', 'photoe' or 'photop'
      c       
    */
  }
  /*
    c     ********************* this is the hook called when
    c       non e-g particle made an interaction.
    c
  */
  void chooknepint_(int   * never){
    /*
      integer never   ! input & output
      c         don't make never = 1, if you want to get
      c         information after a non-e-g particle  made interaction
      c         if this is made non zero, this routine will never be called.
      c
      c   MovedTrack is the particle that made interaction
      c   Pwork contains produced particles.
      c   Nproduced has the number of particles in Pwork
      c   IntInfArray(ProcessNo) contains the type of interaction
      c
      c        default setting
    */
    *never = 1;
    /*
      c
      c        IntInfArray(ProcessNo).process  will have
      c             'col' or 'decay'
    */
  }
}
