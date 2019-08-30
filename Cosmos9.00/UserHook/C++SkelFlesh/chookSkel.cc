/*
  #include "../cmain.cc"
  #include "chookHybAS.cc"
  #include "../ctemplCeren.cc"

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
#include "Zcode.h"
#include <iostream>
#include <string>
#include <strstream>
#include <fstream.h>
#include "Zprivate.h"



using namespace std;
extern "C" {
  extern void cwriteparam_(int *, int *);
  extern void cprintprim_(int *);
  extern void cprintstatus_();
  extern void cprintobs_(int *);
  extern void cquhooki_(int *, int *);
  extern void cquhookr_(int *, double *);
  extern void cquhookc_(int *, char *, int );
  extern int klena_(char *, int );
  extern void cgetfname_(char *, char *, int ,int );
  extern void cerrormsg_(char *, int *, int);
  extern void cqincident_(track *, coord *);
  extern void csetemin_( char *, double *,  double *,  double *,  int );
  extern void cwriteseed_();
  extern void cqeventno_(int  *, int *);
  extern void cqinirn_(int *); 

void cmemorize( fstream *, fstream *);
void cmemoNode( fstream *, int);
void cputHES(fstream *) ;
void cputNodInfo(fstream *,  fstream *);
void time_(int *xxx){
  *xxx = 1;
}
#include "cputNull.cc"
//  put '\0' at the end of a string (after last non blank character).


void chookbgrun_(){
  /*
   *         If you feel writing the parameters on cerr is
   *         a bother, comment out the next or
   *         use other device than ErrorOut.
   *         Also you may comment out all output routines below.
   */
  const int apenl=10;
  char append[apenl];
  const int msgl=100;
  char msg[msgl];
  logical ex, open;
  int dummy=100; 
  
  //            namelist output
  int itemp=0;
  cwriteparam_( &ErrorOut, &itemp);
  //            primary information
  cprintprim_( &ErrorOut);
  //            observation level information
  cprintobs_( &ErrorOut);
  
  int temp[5]={1,2,3,4,5};
  //      Mdev, Wdev are integer in Fortran; but in C++ it must be
  //      fstream ; defined in Zprivate.h
  //    cquhooki_(&temp[0], &Mdev); //      ! get skeleton memo dev #
  //    cquhooki_(&temp[1], &Wdev); //      ! get working disk dev #

  
  cquhooki_(&temp[2], &NgMin); //     ! get Nh min
  cquhooki_(&temp[3], &NhMin); //     ! get Ng min
  cquhooki_(&temp[4], &Where); //     ! where to check
  cquhookr_(&temp[0], &SumegMin);//    ! sum E min
  cquhookr_(&temp[1], &SumehMin); 
  cquhookc_( &temp[0], msg, msgl );// ! get file name for sekelton memo
  cgetfname_( msg, Mskel, msgl, dummy );//  ! add host name etc if needed
  cquhookc_(&temp[1], msg, msgl );//       ! get file name for working
  cgetfname_( msg, Wskel, msgl, dummy) ;//  ! add host name etc if needed
  cquhookc_(&temp[2], append, apenl );//! append data, if Mskel already exists

  ostrstream mout( msg, sizeof msg );  // putting something in mout will make
  // msg string.
  mout << "Skeleton is judged at obs.pos=" << Where << ends;
  itemp = 1;
  cerrormsg_( msg, &itemp, strlen(msg) );
  /*     probably it is  easier to use cerr   like:
	 cerr << "Skeleton is judged at obs.pos=" << Where << endl;
  */
  mout.seekp(0);
  mout << " Ngmin=" << NgMin << " SumEgmin=" << SumegMin/1000. 
       << " TeV" << ends;
  cerrormsg_( msg, &itemp, strlen(msg) );
  mout.seekp(0);
  mout << " Nhmin=" << NhMin << " SumEhmin=" << SumehMin/1000.
       << " TeV" << ends;
  cerrormsg_( msg, &itemp, strlen(msg));
  mout.seekp(0);

  //        binary output file definition


  cputNull(Mskel); // put \0 at the end of file name;
  if( strncmp(append, "append", 6) == 0) {
    // open append mode
    Mdev.open(Mskel, ios::out | ios::app | ios::binary );
    cerr << "skeleton node info. will be appended\n" << endl;
  }
  else {
    // open as new file
    Mdev.open(Mskel, ios::out |  ios::binary );
  }

  if( Mdev.fail() ) {
    cerr << Mskel << " cannot be opened " << endl;
    exit (1);
  }

  cputNull(Wskel);   // put \0 at the end of file name
  //  Wdev.open(Wskel, ios::out | ios::in | ios::binary); // open the file
  Wdev.open(Wskel, ios::out  | ios::binary); // open the file
  // if 'in' is not given, operation equivalent to rewind not possible
  Accepted = 0;  //   ! counter;  accepted as skeleton
}
/*     
       c*********************************** hook for Beginning of  1 event
       c*  All system-level initialization for 1 event generation has been
       c*  eneded at this moment.
       c*  After this is executed, event generation starts.
       c*
*/

void  chookbgevent_(){
  struct track incident;
  struct coord angle;

  static int EventNoLocal = 0;
  int i, j;
  int seed[2];
  double svEasWait, svEthin, kepn;

  Np = 0;
  cqincident_( &incident, &angle);
  kepn = incident.p.fm.p[3];
  if( incident.p.code == kgnuc ) kepn /= incident.p.subcode;
  Ethresh = kepn * WaitRatio;
  
  svEasWait = EasWait;//       ! for safety save
  svEthin = Ethin; //           ! //
  csetemin_( Generate2, &KEminObs2, 
	     &Cutneg, &Cuteg, (int )sizeof( Generate2)  );
  EasWait = svEasWait ;//      ! restore
  Ethin = svEthin ; //
  //  rewind Wdev.  
  //  Wdev.flush();
  //Wdev.seekg(0);
  //  next is ok
  Wdev.close();
  Wdev.open(Wskel,  ios::out | ios::binary);

  //c     ===================================================
  /*
  EventNoLocal++;
  for(int  i = 1; i <= NoOfSites; i++)  
    cout <<  ObsSites[i].pos.depth  //  For obssite, indexing is special
	 << " " <<  EventNoLocal 
	 << " " <<  incident.p.code
	 << " " <<  incident.p.subcode
	 << " " <<  incident.p.charge
	 << " " <<  incident.p.fm.p[3]
	 << " " <<  -angle.r[0]
	 << " " <<  -angle.r[1]	  
	 << " " <<  -angle.r[2] <<  endl;
  */
  //1000       format(f10.3,i9,3i4,e15.5,3(1x,f12.8))
  //c     ===================================================
}
/*     ************************************ hook for observation
 *  One particle information is brought here by the system.
 *  All information of the particle is in aTrack
 */

void chookobs_(track *aTrack, int *id){
  /*
    c
    c     Note that every real variable is in double  precision so
    c     that you may output it in single precision to save the memory.
    c     In some cases it is essential to put it in sigle (say,
    c     for gnuplot).
    c 
    
    integer id  ! input.  1 ==> aTrack is going out from
    c                                 outer boundery.
    c                     2 ==> reached at an observation level
    c                     3 ==> reached at inner boundery.
    

  cout <<  aTrack->where << " " <<  aTrack->p.code << " " << aTrack->p.fm.p[3] << endl;
  */
  //===================================================
  if( *id == 2  &&  aTrack->p.code !=  kneumu  &&
      aTrack->p.code != kneue) {
    // to avoid complexy of indexing, we increment Np after
    // storing ; diff. from Fortran.
    if( Np+1 >  NpMax) {
      cerr << "# of particles >NpMax in observation" << endl;
      exit (1);
    }
    o[Np].where = aTrack->where;
    o[Np].code = aTrack->p.code;
    o[Np].subcode = aTrack->p.subcode;
    o[Np].charge = aTrack->p.charge;
    o[Np].atime = aTrack->t;
    o[Np].erg = aTrack->p.fm.p[3];
    o[Np].mass = aTrack->p.mass;
    o[Np].x = aTrack->pos.xyz.r[0];
    o[Np].y = aTrack->pos.xyz.r[1];
    o[Np].wx =aTrack->vec.w.r[0];
    o[Np].wy =aTrack->vec.w.r[1];
    o[Np].wz =aTrack->vec.w.r[2];
    o[Np].zenith = aTrack->vec.coszenith;
    //     ===================================================

    //    if( o[Np].code <= 6  &&  o[Np].code != 3 ){
    //      cout << o[Np].where
    //   <<  " " <<  o[Np].code
    //   << " " <<   o[Np].charge
    //   << " " <<   o[Np].erg << endl;
    /*
	   << " " <<   o[Np].x 
	   << " " <<   o[Np].y
	   << " " <<   o[Np].wx 
	   << " " <<   o[Np].wy 
	   << " " <<   o[Np].wz 
	   << " " <<   o[Np].zenith << endl;
    }
    */
    Np ++;
    //     ===================================================
  }
  /*
    c
    c     For id =2, you need not output the z value, because it is always
    c     0 (within the computational accuracy).
  */
  
  
  /*  
    if( *id == 2) {
    //      output typical quantities.
      cout << " " << aTrack->where //observation level. int*2. 1 is highest.
	   << " " << aTrack->p.code //  ptcl code.  integer*2.
	   << " " << aTrack->p.charge //  charge,  integer*2 
	   << " " << aTrack->t   //  relateive arrival time in nsec (not sec).
    //  if timestructure is f, nonsense.
	   << " " << aTrack->p.fm.p[3]<< endl; //  ! total energy in gev.

	<< " " << aTrack->pos.xyz.r[0]  //  x, y
	<< " " << aTrack->pos.xyz.r[1] 
	<< " " << aTrack->vec.w.r[0]   //direc. cos.x in the current detector system.
		
	<< " " << aTrack->vec.w.r[1] // direc. cos.y
	<< " " << aTrack->vec.w.r[2] // direc. cos.z
	<< " " << aTrack->vec.coszenith;  //cos of zenith angle
	cout << endl;
    }
  */

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
  integer i;
  integer seed[2];
  real  sumeg, sumeh;
  logical memorize;
  
  integer ng, nh;
      
  ng = 0;
  nh = 0;

  sumeg = 0;
  sumeh =0;

  //c         count sum Eg etc
  for(i =1; i <= Np; i++){
    if( o[i-1].where == Where ) {
      if( o[i-1].code <=  kelec) {
	ng++; 
	sumeg +=  o[i-1].erg;
      }
      else if(  ( o[i-1].code >= kpion &&
		  o[i-1].code <= knuc )  ||
		o[i-1].code ==  kgnuc) {
	nh ++;
	sumeh += o[i-1].erg;
      }
    }
  }
  //     ===================================================
  memorize =(ng >=  NgMin &&  sumeg >= SumegMin) ||
    ( nh >= NhMin && sumeh >= SumehMin);
  //c     ===================================================
  if( memorize ) {
    Accepted++;
    cwriteseed_();//        !  SeedFile
    //c          flag for end of 1 event on working disk
    int minus=-1;
    //	write(Wdev)  -1, p    is now as below; is this the
    //      right way ?
    Wdev.write( (char *) (&minus), sizeof( minus ));
    Wdev.write( (char *) (&p),     sizeof( p ));
    cmemorize(&Wdev, &Mdev) ;//    !  reocord this event
  }
}

/*
  c     ********************************* hook for end of a run
  c     *  all events have been created or time lacks
  c     *
*/

void chookenrun_(){
  cerr << "+++++++++++++++" << endl;
  cerr << Accepted << " events are memorized as skeleton " << endl;
  cerr << "Their  seeds are also memorized" << endl;
  cerr << "--------------" << endl;
  cprintstatus_();
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
void chookeint_(int * never){
  /*
    c 
    c    If Job = 'newskel', input  "never" may be < 0,  and MovedTrack
    c    may not be an electron.
    c      never = -1 :  MovedTrack K.E becomes < KEminObs due to
    c                    energy loss during traveling.
    c      never = -2 :  The same as above, but the particle crosses an
    c                    observation depth, so the calling point to this
    c                    routine is different from the never = -1 case.
    c      never = -3 :  K.E >= KEminiObs.  The ptcl is observed at the
    c                    current deepest Obs. level. but at the flesh time
    c                    the deepest level will be more deep so that
    c                    this must be memorized.
    c
    c         For above cases, the product is only MovedTrack and 
    c         no particle is in PWork.
    c  Otherwise,
    c   MovedTrack is the electron that made interaction
    c   Pwork contains produced particles (normally gamma, but mayb  e).
    c   Nproduced has the number of particles in Pwork (normally 1)
    c   IntInfArray(ProcessNo) contains the type of interaction
    c        IntInfArray(ProcessNo).process will have one of
    c       'brems', 'mscat', 'bscat'  'a1nihi' or 'mbrem'
    c     
  */
  if( *never < 0 ) {
    Nproduced = 1;
    Pwork[0] = MovedTrack.p;
  }
  /*
    c         high energy parent at node might be used
    c        for hybrid AS generation if it is in some
    c        energy region.
  */
  if( MovedTrack.p.code == kelec ) {
    if( MovedTrack.asflag == 0 ) {
      if( MovedTrack.p.fm.p[3] <  Ethresh ) {
	MovedTrack.asflag = -1;
      }
    }
  }     
  cmemoNode( &Wdev, *never );
  if( MovedTrack.asflag == -1) {

    //            af flesh time, decendent should not be used to
    //            generae A.S

    MovedTrack.asflag  = -2;
  }

  *never = 0;
}
/*
  c     ********************* this is the hook called when
  c       a gamma ray made an interaction.
  c
*/
void chookgint_(int * never){
  /*      
	  c         don't make never = 1, if you want to get
	  c         information after a gamma ray made interaction
	  c         if this is made non zero, this routine will never be called.
	  c
	  c   MovedTrack is the gamma that made interaction
	  c   Pwork contains produced particles.
	  c   Nproduced has the number of particles in Pwork
	  c   IntInfArray(ProcessNo) contains the type of interaction
	  c         IntInfArray(ProcessNo).process will have one of
	  c        'pair', 'comp', 'photoe' 'photop' 'mpair'
  */
  cmemoNode( &Wdev, 1 );
  * never = 0;
}
/*
  c     ********************* this is the hook called when
  c       non e-g particle made an interaction.
  c
*/
void chooknepint_(int *never) {
  /*
    c         don't make never = 1, if you want to get
    c         information after a non-e-g particle  made interaction
    c         if this is made non zero, this routine will never be called.
    c
    c   MovedTrack is the particle that made interaction
    c   Pwork contains produced particles.
    c   Nproduced has the number of particles in Pwork
    c   IntInfArray(ProcessNo) contains the type of interaction
    c
  */

  cmemoNode( &Wdev, 1 );
  *never = 0;
}

//    memorize skelton dagta
void cmemoNode( fstream *dev, int flag){
  /*
      integer dev  ! input. fortran logical dev. number for data write
      integer flag ! input. If this is -3, high energy particles must
                   !       also be memorized.
		   c
		   c
		   c       memorize nodal information at interaction.
		   c     q
		   c        \
		   c         \  projectile    = MovedTrack
		   c         /|\  produced particles.: Pwork(i), i=1~Nproduced
		   c        / | \
		   c
		   c   From projectile, the track information is extracted
		   c   From produced particles, only those of K.E< KEminObs is
		   c   extracted and information needed for further tracking is 
		   c   obtaned and memorized. The position information is common
		   c   to all the children.
		   c
		   c   Track   information;   pos =
		   c                                xyz = r(1~3), sys
		   c                                radiallen, depth, height, colheight
		   c                            t
		   c                          vec  =
		   c                                 w = r(1~3), sys
		   c                                 coszenith
		   c                          wgt
		   c                         where
		   c                        asflag 
		   c
		   c      Among these, we don't memorize 
		   c         sys which is always 'xyz'
		   c       radiallen: can be computed from height
		   c         vec; children knows  their direction
		   c         wgt: normally not needed, it should be 1 for skeleton making
		   c              So thinsamping must not be used when making skeleton.
		   c      asflab: should be always F, (assume for skeleton making, hybrid
		   c              air shower is not requested)       
		   c
		   c   write  track info. of projectile
		   c
  */
  int i, nlow;
  double ke;

  nlow = 0 ;
  for(i = 1; i<= Nproduced; i++){ 
    ke = Pwork[i-1].fm.p[3] - Pwork[i-1].mass ;
    if( ( Pwork[i-1].code <= kelec &&  ke >=  Cuteg ) ||
	( Pwork[i-1].code > kelec  && ke >=  Cutneg ) ) {
      //              count low energy ptcls
      if( flag != -3 ) {
	if( ke < KEminObs ){
	  nlow++;
	}
      }
      else {
	  //         all ptcl must be memorized
	  nlow++;
      }
    }
  }
  if( nlow == 0 && MovedTrack.asflag != -1 )  return ; //  *************
  p.posx = MovedTrack.pos.xyz.r[0];
  p.posy = MovedTrack.pos.xyz.r[1];
  p.posz = MovedTrack.pos.xyz.r[2];
  p.depth = MovedTrack.pos.depth;
  p.height = MovedTrack.pos.height;

  if( MovedTrack.pos.colheight > 1.e36) p.colHeight = 1.e36; // no col. yet.
  else  p.colHeight =  MovedTrack.pos.colheight;

  p.atime = MovedTrack.t;
  p.where = MovedTrack.where;
  p.coszenith = MovedTrack.vec.coszenith;
  p.code = MovedTrack.p.code;
  p.erg  = MovedTrack.p.fm.p[3];
  p.asflag = MovedTrack.asflag;

  /*
    cc
    cc     *           p.posx, p.posy, p.posz, p.depth, p.height, 
    cc     *          p.colHeight, p.atime, p.where
    cc
    c       write particle info
    c           p(1~4)
    c           mass
    c           code
    c           subcode
    c           charge
  */

  // write(dev) nlow, p
  (*dev).write( (char *) (&nlow), sizeof( nlow ));
  (*dev).write( (char *) (&p), sizeof( p ));


  for(i = 1; i<=  Nproduced; i++) {
    ke = Pwork[i-1].fm.p[3] - Pwork[i-1].mass ;
    if( (Pwork[i-1].code <= kelec &&  ke >=  Cuteg ) || 
	( Pwork[i-1].code > kelec &&  ke >= Cutneg ) ){
      if(flag == -3  ||  ke < KEminObs ) {
	c.code = Pwork[i-1].code;
	c.subcode =  Pwork[i-1].subcode ;
	c.charge =  Pwork[i-1].charge;
	c.fm[0] = Pwork[i-1].fm.p[0];
	c.fm[1] = Pwork[i-1].fm.p[1];
	c.fm[2] = Pwork[i-1].fm.p[2];
	c.fm[3] = Pwork[i-1].fm.p[3];
	c.mass = Pwork[i-1].mass;
	(*dev).write( (char *) (&c), sizeof( c ));
	
      }
    }
  }
}


void cmemorize(fstream *from, fstream *to){
// from:   working disk file 
// to  :   permanent disk file where skeleton is sotered
  
  int num, cumnum, irevent[2];


  //      rewind from ; rewind member function is not
  // availble in c++ so we use following two; for this
  //  file must be opened with in/out mode
  //  (*from).flush();
  //  (*from).seekg(0);
  (*from).close();
  (*from).open(Wskel, ios::in  | ios::binary); // open the file
  
  //      need not memorize incident, can be generated at flesh
  //      call cqIncident(incident, angle)      

  cqeventno_( &num, &cumnum);
  cqinirn_( &irevent[0] ); // seed of the event
  //      write(to) cumnum, num, irevent,  Zfirst
  (*to).write( (char *)(&cumnum) , sizeof( cumnum ));
  (*to).write( (char *)(&num) , sizeof( num ));
  (*to).write( (char *)(&irevent) , sizeof( irevent ));
  (*to).write( (char *)(&Zfirst) , sizeof( Zfirst ));
  //      call cputHES(to)   ! put high energy showers.
  cputHES( to );
  //      call cputNodInfo(from, to)  ! put nordal info.
  cputNodInfo( from, to) ;
}

void cputHES( fstream *to) {
  //c
  //c        write high energy sekeleton data into 'to'      
  //c
  integer i;
  
  (*to).write( (char *)(&Np), sizeof(Np));

  for(i = 1;  i<=Np; i++) {
    (*to).write( (char *)( &o[i-1] ), sizeof( o[i-1] ));
  }
}

void cputNodInfo( fstream *from,  fstream *to){
  integer i, nlow;

  nlow = 1;
  while ( nlow >= 0 ) {
    //    read(from) nlow,  p
    (*from).read( (char *)(&nlow), sizeof(nlow));
    (*from).read( (char *)(&p), sizeof(p));
    //    write(to)  nlow,  p
    (*to).write( (char *)(&nlow), sizeof(nlow));
    (*to).write( (char *)(&p), sizeof(p));
    for(i = 1; i<=nlow; i++) {
      (*from).read( (char *)(&c), sizeof(c));
      (*to).write( (char *)(&c), sizeof(c));
    }
  }
}

}

      
