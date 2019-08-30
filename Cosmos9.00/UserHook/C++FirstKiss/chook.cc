/*
  *************************************** hook for Beginning of a Run
  * At this moment, all (system-level) initialization for this run
  * has been ended.  After this routine is executed, the system goes into the
  * event creation loop.
  *
*/

#include "Zdef.h"
#include "Zmanagerp.h"
#include "Ztrack.h"
#include "Ztrackv.h"
#include "Ztrackp.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"
#include <iostream>
#include <string.h>
//         next is needed for C++
extern "C" {
  extern void cwriteparam_(int *, int *);
  extern void cprintprim_(int *);
  extern void cprintobs_(int *);

  void time_(int *xxx){
    *xxx = 1;
  }

  void chookbgrun_(){
    /*
     *         If you feel writing the parameters on cerr is
     *         a bother, comment out the next or
     *         use other device than ErrorOut.
     *         Also you may comment out all output routines below.
     */
    //            namelist output
    int temp=0;
    cwriteparam_( &ErrorOut, &temp);
    //            primary information
    cprintprim_( &ErrorOut );
    //            observation level information
    cprintobs_( &ErrorOut );
  }
  /*     
	 c*********************************** hook for Beginning of  1 event
	 c*  All system-level initialization for 1 event generation has been
	 c*  eneded at this moment.
	 c*  After this is executed, an event generation starts.
	 c*
  */
  void  chookbgevent_(){

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

    if( *id == 2) {
      //      output typical quantities.
      std::cout<< " " << atrack->where //observation level. 1 is highest.
       << " " << atrack->p.code //  ptcl code.  
       << " " << atrack->p.charge //  charge,  
       << " " << atrack->t   //  relateive arrival time in nsec (not sec).
                        	//  if timestructure is f, nonsense.
       << " " << atrack->p.fm.p[3] //  ! total energy in GeV.
       << " " << atrack->pos.xyz.r[0]  //  x
       << " " << atrack->pos.xyz.r[1]  //  y
       << " " << atrack->vec.w.r[0]   //direc. cos.x in the current detector system.
       << " " << atrack->vec.w.r[1] // direc. cos.y
       << " " << atrack->vec.w.r[2] // direc. cos.z
       << " " << atrack->vec.coszenith;  //cos of zenith angle; due to the
      // earth curverture, this is little bit diff. from w.r[2]

	//       << " " << atrack->pos.xyz.r[2] // z; this is almost 0
	//       << " " << atrack->pos.radiallen  
	//       << " " << atrack->pos.depth
	//       << " " << atrack->pos.height
	//       << " " << atrack->pos.colheight
      std::cout<< std::endl;
    }
    
  }
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

  /*
    c    *********************************** hook for end of 1 event
    c    * At this moment, 1 event generation has been ended.
    c    *
  */
  void chookenevent_(){
    int i;
    // hybrid AS observation or not
    if( ObserveAS ) {
      //                   electron size in B approx. and age
      for(i=1; i<=NoOfASSites ; i++) {
	std::cout<<  ASObsSites[i-1].esize << " " <<
	  ASObsSites[i-1].age << " ";
	std::cout<< std::endl;
      }
    }
  }

  /*
    c     ********************************* hook for end of a run
    c     *  all events have been created or time lacks
    c     *
  */
  using namespace std;
  extern int klena_(char *, int );
  void chookenrun_(){
     char  tracefile[74];
     cout<< 
      "****** Congratulations: Cosmos is now your friend *******"<< endl;
    if( Trace > 0 ) {
      strcpy(tracefile, TraceDir );
      tracefile[ klena_(tracefile, strlen(tracefile)) ]='\0';
      strcat(tracefile, "/trace?");
      cout<< "   particle trace data has been created" 
	  << " in " << tracefile << endl;
      cout<< " where ?=1,2..you can see it by gnuplot: for that, in gnuplot do"
	   << endl;
      cout<<"  set para" << endl;
      cout<<"  splot " << tracefile <<" w  l " << endl;
      cout<<"  to see charged particles only, use following: "<< endl;
      cout<<"  splot  \" < awk '$6 != 0 ; " << endl;
      cout<<"         NF == 0 ;' " <<  tracefile << " \" w l"<< endl;
      cout<< endl;
      cout<<"For a far better plot use Geomview in Util/Geomview"
	   << endl;
      cout<< "****REMEMBER You have to make Trace=0 in Param for NON-Demo run****";
      cout<< endl;
      cout<<"************      Have a nice day !!      **************"
      << endl;
   }
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

        c  atrack. here must be one either of  ztrackv_.trackbefmove.    or
        c                      ztrackv_.movedtrack.

	c     atrack.pos.radiallen   (distance from the center of the earth)
	c     atrack.pos.depth       (vertical depth)
	c     atrack.pos.height      (vertical heigth from sea level)  
	c
      */
    float h1,  h2;

    //h1 height of the  beginning of the segment   from the base ;
    //h2 the same for the end of the segment.
    //  Fortran index 0, 1, 2,   ... NoOfSites, NoOfSites+1
    //                LB                         UB
    //  C             0  1, 2,    ...
    //  They are the same in this case !!
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

      
