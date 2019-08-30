#include "Zmaxdef.h"
!          common variables used in tracking ptcls.
       integer ToInteract, ToBeObserved, Truncated, Dead,
     *         BorderL,  BorderH, AngleLimit
       parameter(ToInteract = 1, ToBeObserved = 2, Truncated = 3,
     *  BorderL = 4, BorderH =5,  Dead = 6, AngleLimit = 7)
       integer  BitPhotoElec, BitPhoton,
     *    BitElectron, BitPositron, BitProton,
     *    BitNeutron, BitAntiNuc, BitDecay,  BitOther, BitEconsv
       parameter( BitPhotoElec=1, BitPhoton=2, BitElectron=3, 
     *  BitPositron=3, BitProton=4, BitNeutron=5, BitAntiNuc=6,
     *  BitDecay=7, BitOther=8, BitEconsv=9  )

!!!       integer MaxInte
!!!       parameter(MaxInte = 6)  !  Max number of kinds of interactions a particle can
!!!                               !  take. (such as brems, knockon, anihilation)
  !!!             moved to csetIntInf
!!!       type intinf      ! Interaction information
!!!         sequence
!!!           real*8  thickness   ! in kg/m2 set if decay is F
!!!           real*8  length      ! in m, set if decay is T.  or eventually by cfixProc
!!!           character*8 process ! process id string such as brems, pair
!!!            logical decay       ! if decay, T, else F
!!!       end type intinf
!!!!          define array of intinf
!!!       type(intinf) IntInfArray(MaxInte)
!!!
!!!       integer NumberOfInte ! Number of different kind of interactions 
!!!                            ! considered for the current particle.
!!!       integer ProcessNo   ! The process really happend is the
!!!                            ! ProcessNo-th process in  IntInfArray.
!!!
!
       type(track):: TrackBefMove  ! track before moved       	
       type(track):: MovedTrack      ! to contain track moved
       type(coord):: Offset        ! the primary is directed to 
!                                     deepest detector origin + Offset
!                                     (in 'xyz')
       type(track):: Zfirst ! to keep first interaction info. V7.0
!       real*8 Zfirst       ! to keep first interaction slant depth
       integer MoveStat    ! status code for moving a particle
         logical ObserveAS     ! made to be T, if AS is to be generated
        logical Upgoing     ! if primary is going upward, made to be t
        logical UseTbl       ! becomes T, 
                            ! if length <--> thickness conv. is by table
         real*8  EminAS      ! minimum energy of e for AS generation.
         real*8  EasWait     ! for AS generation, must wait until e 
                            ! energy becomes < EasWait
        real*8 EnergyLoss   !  energy loss 
!        real*8 Upsilon      ! Upsilon value   --> modEMcontrol
!        real*8 Xai          ! Xai value B x Eg/m /2  --> modEmcontrol
        real*8 maxstep(0:MAX_NO_OF_SITES+1) ! used to cut the path
                ! 1/5 of the depth step. this is necessary cond.
         real*8  KEmin         ! min kinetic energy to be tracked
        real*8  KEminCas      ! //          (for em-cascade)
        real*8  KEmin2        ! min kinetic energy to be tracked.  for skeleton/flesh use.
        real*8  KEminCas2     ! for skeleton/flesh use.
         real*8  Ethin(4)      ! Thin sampling threshold and max weight.  for e/g and hadrons/muons
         real*8 Beta  !  v/c for MovedTrack; given if TimeStrucrue=T.
        type(magfield):: Mag
!         integer MaxPtcl   --> Zpwork.h
        logical FromEpics     !  to control muon iteraction (pair,brem,nuci)
                              !  must be made t, when Epics treats muon.
                              !  if Cosmos uses Epics, this must be made to
                              !  be t/f depending on Epics mode, or Cosmos mode
#if LABELING > 0
        integer Labelcounter   ! label counter to put a lalel on each patcl.
#endif


!->MuNucon         real*8  MuonPolarization  ! muon polarization value.
        integer FirstColA, FirstColZ
        real(8)::FirstColXs
!
!       common /Ztrackv/Pwork, IntInfArray, TrackBefMove,
!        common /Ztrackv/Pwork, TrackBefMove,  ! Pwork--> modColInfo
        common /Ztrackv/ TrackBefMove,        
     *  MovedTrack, Zfirst, Offset,
     *  Mag, EminAS, EasWait,
!->MuNucon	MuonPolarization,  
     *  EnergyLoss,  KEmin, KEminCas, Beta, 
     *  KEmin2, KEminCas2, Ethin, !!!! Upsilon, Xai, 
     *  maxstep,
!     *  ObserveAS, Nproduced, Nstacked,
     *  ObserveAS, 
!     *  NumberOfInte, ProcessNo,  !-->modColInfo
     *  MoveStat, 
     *  Upgoing, UseTbl, FirstColXs,
     *  FromEpics, FirstColA, FirstColZ
#if LABELING > 0
     * , Labelcounter
#endif

! #include "Zair.h"






