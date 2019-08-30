#include "Zmaxdef.h"
/*
c          common variables used in tracking ptcls.
       integer ToInteract, ToBeObserved, Truncated, Dead,
     *         BorderL,  BorderH, AngleLimit
       parameter(ToInteract = 1, ToBeObserved = 2, Truncated = 3,
     *  BorderL = 4, BorderH =5,  Dead = 6, AngleLimit = 7)

       integer MaxInte
       parameter(MaxInte = 6)  !  Max number of kinds of interactions a particle can
                               !  take. (such as brems, knockon, anihilation)
       structure /intinf/      ! Interaction information
           real*8  thickness   ! in kg/m2 set if decay is F
           real*8  length      ! in m, set if decay is T.  or eventually by cfixProc
           character*8 process ! process id string such as brems, pair
	   logical decay       ! if decay, T, else F
       end structure   
c          define array of intinf
       record /intinf/ IntInfArray(MaxInte)
c
       record /track/ TrackBefMove  ! track before moved       	
       record /track/ MovedTrack      ! to contain track moved
       record /coord/ Offset        ! the primary is directed to 
c                                     deepest detector origin + Offset
c                                     (in 'xyz')
       record /track/ Zfirst ! to keep first interaction info. V7.0
c       real*8 Zfirst       ! to keep first interaction slant depth
       integer MoveStat    ! status code for moving a particle
       real*8 TargetMassN  ! average target mass number. To be updated
                           ! after a ptcl is moved.
       real*8 TargetAtomicN   ! Average Z of the target.
       real*8 TargetZ2        ! <Z^2> of the target
       integer TargetNucleonNo  !  target nucleon number at a collision
       integer TargetProtonNo  !  target proton number  //
       integer NumberOfInte ! Number of different kind of interactions 
                            ! considered for the current particle.
       integer ProcessNo   ! The process really happend is the
                            ! ProcessNo-th process in  IntInfArray.
	logical ObserveAS     ! made to be T, if AS is to be generated
        logical Upgoing     ! if primary is going upward, made to be t
        logical UseTbl       ! becomes T, 
                            ! if length <--> thickness conv. is by table
	real*8  EminAS      ! minimum energy of e for AS generation.
	real*8  EasWait     ! for AS generation, must wait until e 
                            ! energy becomes < EasWait
        real*8 EnergyLoss   !  energy loss 
        real*8 Upsilon      ! Upsilon value 
        real*8 Xai          ! Xai value B x Eg/m /2
                   
	real*8  KEmin         ! min kinetic energy to be tracked
        real*8  KEminCas      ! //          (for em-cascade)
        real*8  KEmin2        ! min kinetic energy to be tracked.  for skeleton/flesh use.
        real*8  KEminCas2     ! for skeleton/flesh use.
	real*8  Ethin(2)      ! Thin sampling threshold.
	real*8 Beta  !  v/c for MovedTrack; given if TimeStrucrue=T.
	record /magfield/ Mag
	integer MaxPtcl

        logical FromEpics     !  to control muon iteraction (pair,brem,nuci)
                              !  must be made t, when Epics treats muon.
                              !  if Cosmos uses Epics, this must be made to
                              !  be t/f depending on Epics mode, or Cosmos mode
#if LABELING > 0
        integer Labelcounter   ! label counter to put a lalel on each patcl.
#endif

	parameter (
#ifdef MAX_PTCL
     *     MaxPtcl = MAX_PTCL
#else
     *     MaxPtcl = 8000
#endif
     *         )   ! max # of ptcls producable in coll.
	record /ptcl/ Pwork(MaxPtcl)  ! working array to store ptcls.
	integer Nproduced   ! no. of ptcls produced and stored in Pwork.
	real*8  MuonPolarization  ! muon polarization value.
c
*/







#ifdef MAX_PTCL
const int maxptcl =  MAX_PTCL;
#else
const int maxptcl = 8000;
#endif
const int maxinte = 6;

const int tointeract = 1;
const int tobeobserved = 2;
const int truncated = 3;
const int dead = 6;
const int borderl = 4;
const int borderh =5;
const int anglelimit = 7;

struct intinf {
  double  thickness; //   ! in kg/m2 set if decay is F
  double  length; //      ! in m, set if decay is T.  or eventually by cfixProc
  char    process[8]; // ! process id string such as brems, pair
  logical     decay;  //      ! if decay, T, else F
  DUMMYCHAR
};

extern struct ztrackv {
  struct ptcl pwork[maxptcl];
  struct intinf intinfarray[maxinte];
  struct track  trackbefmove;
  struct track  movedtrack;
  struct track  zfirst;
  struct coord  offset;
  struct magfield mag;
  double muonpolarization;
  double eminas;
  double easwait;
  double targetmassn;
  double targetatomicn;
  double targetz2;
  double energyloss;
  double kemin;
  double kemincas;
  double beta;
  double kemin2;
  double kemincas2;
  double ethin[2];
  double upsilon;
  double xai;
  logical    observeas;
  int    nproduced;
  int    movestat;
  int    numberofinte;
  int    processno;
  int    targetnucleonno;
  int    targetprotonno;
  logical    upgoing;
  logical    usetbl;
  logical    fromepics;
#if LABELING > 0
  int     labelcounter;
#endif
} ztrackv_;



#define Pwork ztrackv_.pwork
#define IntInfArray ztrackv_.intinfarray
#define TrackBefMove ztrackv_.trackbefmove
#define MovedTrack ztrackv_.movedtrack
#define Zfirst ztrackv_.zfirst
#define Offset ztrackv_.offset
//    Mag cannot be used so you have to use Mag_ instead of it
#define Mag_ ztrackv_.mag
#define MuonPolarization ztrackv_.muonpolarization
#define EminAS ztrackv_.eminas
#define EasWait ztrackv_.easwait
#define TargetMassN ztrackv_.targetmassn
#define TargetAtomicN ztrackv_.targetatomicn
#define TargetZ2 ztrackv_.targetz2
#define EnergyLoss ztrackv_.energyloss
#define KEminCas2 ztrackv_.kemincas2
#define KEminCas ztrackv_.kemincas
#define KEmin2 ztrackv_.kemin2
#define KEmin ztrackv_.kemin
#define Beta ztrackv_.beta
#define Ethin ztrackv_.ethin
#define Upsilon ztrackv_.upsilon
#define Xai ztrackv_.xai
#define ObserveAS ztrackv_.observeas
#define Nproduced ztrackv_.nproduced
#define MoveStat ztrackv_.movestat
#define NumberOfInte ztrackv_.numberofinte
#define ProcessNo ztrackv_.processno
#define TargetNucleonNo ztrackv_.targetnucleonno
#define TargetProtonNo ztrackv_.targetprotonno
#define Upgoing ztrackv_.upgoing
#define UseTbl ztrackv_.usetbl
#define FromEpics ztrackv_.fromepics
#define Labelcounter ztrackv_.labelcounter

